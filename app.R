# ============================================================================
# Interactive Empirical Variogram Explorer
# ----------------------------------------------------------------------------
# A Shiny app to explore how distance-binning strategies (cutoff, width)
# affect estimation of the empirical semi-variogram, and to select an
# "optimal" number of distance bins using a model-fit criterion.
#
# This app has NO built-in dataset — you must upload your own CSV.
#
# REQUIRED CSV FORMAT
# --------------------
#   - A plain comma-separated CSV with a header row (column names in row 1).
#   - Must contain at least three numeric columns:
#       one for the X coordinate, one for the Y coordinate, and one for the
#       attribute (Z) whose spatial correlation you want to study.
#   - Column names can be anything (e.g. "longitude", "station_x", "value")
#     — after upload, you map which column plays which role using the
#     X column / Y column / Z (attribute) column dropdowns in the sidebar.
#   - X, Y — numeric coordinates:
#       * If they are already projected/planar coordinates (e.g. UTM easting/
#         northing, State Plane, any coordinate system measured in meters),
#         leave "Coordinates are lon/lat" UNCHECKED.
#       * If they are geographic coordinates (decimal-degree longitude and
#         latitude, e.g. longitude in [-180, 180], latitude in [-90, 90]),
#         CHECK "Coordinates are lon/lat" and choose either:
#           - "Project to UTM"      -> distances computed in meters, or
#           - "Great-circle (haversine)" -> distances computed in kilometers.
#   - Z — the numeric attribute being analyzed (e.g. pollutant concentration,
#     soil moisture, elevation, disease rate, etc.). Must be numeric.
#   - Rows with a missing/non-numeric value in the X, Y, or Z column are
#     dropped automatically before computing anything.
#   - Duplicate coordinates are fine — they simply contribute zero-distance
#     pairs, which help estimate the nugget effect.
#   - Rule of thumb: aim for at least ~50-100 points. Variograms computed
#     from very few points are unstable and hard to interpret.
#
#   Example CSV (three columns named x, y, z; any names work once mapped):
#     x,y,z
#     178605,331035,4.6
#     178650,330995,4.9
#     178720,331090,5.3
#     179010,330810,4.1
#     ...
#
#   Example CSV with lon/lat and a differently-named attribute column:
#     station_longitude_deg,station_latitude_deg,NO2
#     4.8952,52.3702,34.2
#     4.9000,52.3650,31.8
#     ...
#
# To run:
#   install.packages(c("shiny","gstat","sf","ggplot2","dplyr","sp"))
#   shiny::runApp("variogram_explorer_app.R")
# ============================================================================

library(shiny)
library(gstat)
library(sf)
library(sp)
library(ggplot2)
library(dplyr)

# ---------------------------------------------------------------------------
# UI
# ---------------------------------------------------------------------------
ui <- fluidPage(
  
  titlePanel("Interactive Empirical Variogram Explorer"),
  
  sidebarLayout(
    
    sidebarPanel(
      width = 3,
      
      h4("Upload spatial data (CSV)"),
      fileInput("file", NULL, accept = c(".csv")),
      uiOutput("uploadStatus"),
      
      helpText("Needs an X column, a Y column, and a numeric attribute (Z)."),
      
      tags$details(
        style = "font-size: 85%; color:#333; margin-bottom:8px;",
        tags$summary(style = "cursor:pointer; color:#2b6cb0;",
                     "See CSV format details"),
        tags$div(
          style = "background:#f7f7f7; border:1px solid #ddd; border-radius:4px;
                   padding:8px 10px; margin-top:6px;",
          tags$ul(
            style = "margin:6px 0 6px 18px; padding:0;",
            tags$li("Header row + at least 3 numeric columns: X, Y, and an attribute Z."),
            tags$li("Your file can have any number of other columns/variables — ",
                    "just pick which one is X, Y, and Z below after upload."),
            tags$li("X/Y: projected meters, OR lon/lat degrees (tick the box below)."),
            tags$li("Z: numeric attribute to analyze (e.g. concentration, elevation)."),
            tags$li("Rows with missing/non-numeric X, Y, or Z are dropped automatically."),
            tags$li("Aim for \u2265 ~50 points for a stable variogram.")
          ),
          tags$b("Example:"),
          tags$pre(
            style = "font-size:80%; background:white; padding:4px 6px; margin:4px 0;",
            "x,y,z\n178605,331035,4.6\n178650,330995,4.9\n179010,330810,4.1"
          )
        )
      ),
      
      uiOutput("colSelectors"),
      
      checkboxInput("lonlat", "Coordinates are lon/lat", value = FALSE),
      
      conditionalPanel(
        condition = "input.lonlat == true",
        radioButtons(
          "projMethod", NULL,
          choices = c(
            "Project to UTM (accurate distances in meters)" = "utm",
            "Great-circle distance (haversine, km)"          = "haversine"
          ),
          selected = "utm"
        )
      ),
      
      actionButton("compute", "Compute", class = "btn-primary"),
      
      tags$hr(),
      verbatimTextOutput("dataSummary")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = "mainTabs",
        
        # ---------------------------------------------------------------
        tabPanel(
          "1. Variogram Cloud",
          br(),
          p("Every pair of points contributes one point to the cloud: ",
            "x = distance between the pair, y = semi-variance ",
            "\u00bd(z\u1d62 \u2212 z\u2c7c)\u00b2. This is the raw, unbinned material",
            " that the binned variogram in Tab 2 is built from."),
          plotOutput("cloudPlot", height = "480px")
        ),
        
        # ---------------------------------------------------------------
        tabPanel(
          "2. Binned Semi-variogram",
          br(),
          fluidRow(
            column(
              4,
              wellPanel(
                checkboxInput("compareMode", "Compare multiple binning settings",
                              value = FALSE),
                
                conditionalPanel(
                  condition = "!input.compareMode",
                  h4("Binning parameters"),
                  numericInput("cutoff", "Cutoff (max distance):", value = NA, min = 0),
                  numericInput("width", "Bin width:", value = NA, min = 0),
                  actionButton("resetDefault", "Reset to default (cutoff/3, width=cutoff/15)"),
                  tags$hr(),
                  h4("Find optimal number of bins"),
                  p(style = "font-size:90%;",
                    "Fixes the cutoff above, sweeps the number of bins from ",
                    "5 to 30, fits a spherical model to each binned variogram, ",
                    "and reports the bin count that minimizes the weighted ",
                    "model-fit error (SSErr), subject to a minimum of ~30 ",
                    "pairs per bin (Journel & Huijbregts rule of thumb)."),
                  actionButton("findOptimal", "Find optimal bins", class = "btn-primary"),
                  br(), br(),
                  actionButton("applyOptimal", "Apply optimal width to width box")
                ),
                
                conditionalPanel(
                  condition = "input.compareMode",
                  h4("Compare up to 4 settings"),
                  p(style = "font-size:85%; color:#555;",
                    "Tick the settings you want to overlay on one plot, and ",
                    "type a cutoff/width for each. Values outside the ",
                    "allowed range are clamped automatically."),
                  actionButton("fillDefaultsCompare", "Fill with sensible defaults"),
                  tags$hr(style = "margin:8px 0;"),
                  
                  checkboxInput("cmp_on1", "Setting 1", value = TRUE),
                  fluidRow(
                    column(6, numericInput("cmp_cutoff1", "Cutoff", value = NA, min = 0)),
                    column(6, numericInput("cmp_width1", "Width", value = NA, min = 0))
                  ),
                  tags$hr(style = "margin:8px 0;"),
                  
                  checkboxInput("cmp_on2", "Setting 2", value = TRUE),
                  fluidRow(
                    column(6, numericInput("cmp_cutoff2", "Cutoff", value = NA, min = 0)),
                    column(6, numericInput("cmp_width2", "Width", value = NA, min = 0))
                  ),
                  tags$hr(style = "margin:8px 0;"),
                  
                  checkboxInput("cmp_on3", "Setting 3", value = FALSE),
                  fluidRow(
                    column(6, numericInput("cmp_cutoff3", "Cutoff", value = NA, min = 0)),
                    column(6, numericInput("cmp_width3", "Width", value = NA, min = 0))
                  ),
                  tags$hr(style = "margin:8px 0;"),
                  
                  checkboxInput("cmp_on4", "Setting 4", value = FALSE),
                  fluidRow(
                    column(6, numericInput("cmp_cutoff4", "Cutoff", value = NA, min = 0)),
                    column(6, numericInput("cmp_width4", "Width", value = NA, min = 0))
                  )
                )
              )
            ),
            column(
              8,
              conditionalPanel(
                condition = "!input.compareMode",
                plotOutput("binnedPlot", height = "380px"),
                tags$div(
                  style = "background:#f0f4f8; border-left:4px solid #2b6cb0;
                           padding:8px 12px; margin-top:8px; font-size:90%;",
                  tags$b("What do the numbers next to each point mean? "),
                  "Each label is ", tags$code("np"), " — the number of ",
                  "point-pairs (observations) that fell into that distance ",
                  "bin and were averaged to get \u03b3(h) for that bin. More ",
                  "pairs \u2192 a more statistically reliable estimate. As a ",
                  "rule of thumb, bins with fewer than ~30 pairs (common at ",
                  "larger distances, where fewer pairs exist) should be read ",
                  "with caution — they can look noisy or unreliable on the plot."
                )
              ),
              conditionalPanel(
                condition = "input.compareMode",
                plotOutput("comparePlot", height = "440px"),
                tags$div(
                  style = "background:#f0f4f8; border-left:4px solid #2b6cb0;
                           padding:8px 12px; margin-top:8px; font-size:90%;",
                  tags$b("Comparing settings: "),
                  "Each colored line/curve is the binned variogram produced ",
                  "by one cutoff/width combination. Use this to see directly ",
                  "how a smaller bin width gives a noisier but more detailed ",
                  "curve, while a larger width gives a smoother but coarser one."
                )
              )
            )
          ),
          tags$hr(),
          conditionalPanel(
            condition = "!input.compareMode",
            fluidRow(
              column(6, plotOutput("optimalPlot", height = "300px")),
              column(6,
                     h4("Recommendation"),
                     textOutput("optimalText"),
                     tags$br(),
                     tableOutput("optimalTable"))
            )
          )
        ),
        
        # ---------------------------------------------------------------
        tabPanel(
          "3. How Binning Works",
          br(),
          p("Walks through how the binned variogram is actually built from ",
            "the cloud, step by step. Set your own cutoff/width here, or ",
            "they start out matching whatever is set in Tab 2."),
          fluidRow(
            column(4, numericInput("cutoffStep", "Cutoff:", value = NA, min = 0)),
            column(4, numericInput("widthStep", "Width:", value = NA, min = 0))
          ),
          sliderInput("step", "Step:", min = 1, max = 4, value = 1, step = 1,
                      width = "60%"),
          plotOutput("stepPlot", height = "480px"),
          conditionalPanel(
            condition = "input.step >= 3",
            tags$div(
              style = "background:#f0f4f8; border-left:4px solid #c0392b;
                       padding:8px 12px; margin-top:8px; font-size:90%;",
              tags$b("Note on the numbers shown at each bin average: "),
              "The number next to each red point is ", tags$code("np"),
              " — how many point-pairs (observations) were averaged to ",
              "produce that bin's \u03b3(h). Bins built from more pairs are ",
              "more trustworthy; a bin with very few pairs is a less ",
              "reliable estimate of the true spatial correlation at that distance."
            )
          )
        )
      )
    )
  )
)

# ---------------------------------------------------------------------------
# SERVER
# ---------------------------------------------------------------------------
server <- function(input, output, session) {
  
  # ---- Upload status banner ----------------------------------------------
  output$uploadStatus <- renderUI({
    req(input$file)
    tags$div(
      style = "background:#2b6cb0; color:white; padding:6px; margin:6px 0;
               text-align:center; border-radius:4px; font-size:90%;",
      "Upload complete"
    )
  })
  
  # ---- Column selectors, populated once a CSV is uploaded -----------------
  output$colSelectors <- renderUI({
    req(input$file)
    hdr <- tryCatch(read.csv(input$file$datapath, nrows = 1), error = function(e) NULL)
    req(hdr)
    cols <- names(hdr)
    
    guess_x <- cols[grepl("lon|longitude|^x$", cols, ignore.case = TRUE)]
    guess_y <- cols[grepl("lat|latitude|^y$",  cols, ignore.case = TRUE)]
    guess_x <- if (length(guess_x) > 0) guess_x[1] else cols[1]
    guess_y <- if (length(guess_y) > 0) guess_y[1] else cols[min(2, length(cols))]
    remaining <- setdiff(cols, c(guess_x, guess_y))
    guess_z <- if (length(remaining) > 0) remaining[1] else cols[min(3, length(cols))]
    
    tagList(
      tags$b("X column"),
      selectInput("xCol", NULL, choices = cols, selected = guess_x),
      tags$b("Y column"),
      selectInput("yCol", NULL, choices = cols, selected = guess_y),
      tags$b("Z (attribute) column"),
      selectInput("zCol", NULL, choices = cols, selected = guess_z)
    )
  })
  
  # ---- Load & validate data; fires only when "Compute" is clicked --------
  raw_df <- eventReactive(input$compute, {
    validate(need(input$file, "Please upload a CSV file, then click Compute."))
    validate(need(input$xCol, "Select the X column."))
    validate(need(input$yCol, "Select the Y column."))
    validate(need(input$zCol, "Select the Z (attribute) column."))
    
    full <- read.csv(input$file$datapath)
    validate(need(all(c(input$xCol, input$yCol, input$zCol) %in% names(full)),
                  "Selected columns were not found in the uploaded file."))
    
    df <- data.frame(x = suppressWarnings(as.numeric(full[[input$xCol]])),
                     y = suppressWarnings(as.numeric(full[[input$yCol]])),
                     z = suppressWarnings(as.numeric(full[[input$zCol]])))
    df <- df[stats::complete.cases(df), ]
    
    validate(need(nrow(df) >= 10,
                  "Need at least 10 valid numeric rows (after removing missing/non-numeric values) to compute a variogram."))
    
    list(df = df,
         lonlat = isTRUE(input$lonlat),
         proj = if (isTRUE(input$lonlat)) input$projMethod else "none")
  }, ignoreNULL = FALSE)
  
  # ---- Build the sf object, handling lon/lat projection choice -----------
  data_sf <- reactive({
    rd <- raw_df()
    df <- rd$df
    
    if (rd$lonlat) {
      sfobj <- st_as_sf(df, coords = c("x", "y"), crs = 4326, remove = FALSE)
      if (identical(rd$proj, "utm")) {
        lon_mean <- mean(df$x, na.rm = TRUE)
        lat_mean <- mean(df$y, na.rm = TRUE)
        zone <- floor((lon_mean + 180) / 6) + 1
        epsg <- if (lat_mean >= 0) 32600 + zone else 32700 + zone
        sfobj <- st_transform(sfobj, crs = epsg)
      }
      # else "haversine": keep geographic CRS 4326 — gstat's variogram()
      # automatically computes great-circle distances (in km) for
      # unprojected (lon/lat) coordinates.
    } else {
      sfobj <- st_as_sf(df, coords = c("x", "y"), remove = FALSE)
    }
    sfobj
  })
  
  output$dataSummary <- renderText({
    rd <- raw_df()
    df <- rd$df
    dist_unit <- if (rd$lonlat && identical(rd$proj, "haversine")) "km (great-circle)" else "m (projected)"
    paste0(
      "Points: ", nrow(df), "\n",
      "X range: ", round(min(df$x), 4), " to ", round(max(df$x), 4), "\n",
      "Y range: ", round(min(df$y), 4), " to ", round(max(df$y), 4), "\n",
      "Z range: ", round(min(df$z), 3), " to ", round(max(df$z), 3), "\n",
      "Distance unit used: ", dist_unit
    )
  })
  
  # ---- Default cutoff/width per gstat convention (from your slides) ------
  # cutoff  = 1/3 of the bounding box diagonal
  # width   = cutoff / 15
  default_params <- reactive({
    df <- raw_df()$df
    diagn <- sqrt((max(df$x) - min(df$x))^2 + (max(df$y) - min(df$y))^2)
    cutoff <- diagn / 3
    width  <- cutoff / 15
    list(cutoff = cutoff, width = width)
  })
  
  # Set numeric-box default values whenever new data is computed
  # (users are free to type any value they want afterward — no artificial limits)
  observeEvent(data_sf(), {
    dp <- default_params()
    updateNumericInput(session, "cutoff", value = round(dp$cutoff, 2))
    updateNumericInput(session, "width", value = round(dp$width, 2))
    updateNumericInput(session, "cutoffStep", value = round(dp$cutoff, 2))
    updateNumericInput(session, "widthStep", value = round(dp$width, 2))
  }, ignoreInit = TRUE)
  
  observeEvent(input$resetDefault, {
    dp <- default_params()
    updateNumericInput(session, "cutoff", value = round(dp$cutoff, 2))
    updateNumericInput(session, "width",  value = round(dp$width, 2))
  })
  
  # ---- Tab 1: Variogram cloud --------------------------------------------
  cloud_data <- reactive({
    variogram(z ~ 1, data_sf(), cloud = TRUE)
  })
  
  output$cloudPlot <- renderPlot({
    vc <- cloud_data()
    ggplot(vc, aes(x = dist, y = gamma)) +
      geom_point(color = "blue", alpha = 0.35, size = 1.3) +
      labs(title = "Variogram Cloud",
           x = "Distance, h", y = expression(gamma(h))) +
      theme_minimal(base_size = 14)
  })
  
  # ---- Tab 2: Binned variogram (single setting) ---------------------------
  binned_variogram <- reactive({
    req(input$cutoff, input$width)
    variogram(z ~ 1, data_sf(), cutoff = input$cutoff, width = input$width)
  })
  
  output$binnedPlot <- renderPlot({
    v <- binned_variogram()
    ggplot(v, aes(x = dist, y = gamma)) +
      geom_line(color = "darkred") +
      geom_point(shape = 1, size = 3, stroke = 1.2, color = "darkred") +
      geom_text(aes(label = np), hjust = -0.35, vjust = -0.35, size = 3.3,
                color = "black") +
      labs(title = paste0("Binned Semi-variogram  (cutoff = ", round(input$cutoff, 1),
                          ",  width = ", round(input$width, 1), ")"),
           subtitle = "Labels show number of point-pairs (np) per bin",
           x = "Distance, h", y = expression(gamma(h))) +
      theme_minimal(base_size = 14)
  })
  
  # ---- Compare mode: overlay several binned variograms ---------------------
  observeEvent(input$fillDefaultsCompare, {
    dp <- default_params()
    cutoff_val <- round(dp$cutoff, 2)
    widths <- round(c(dp$width / 2, dp$width, dp$width * 2, dp$width * 4), 2)
    
    for (i in 1:4) {
      updateNumericInput(session, paste0("cmp_cutoff", i), value = cutoff_val)
      updateNumericInput(session, paste0("cmp_width", i), value = widths[i])
    }
    updateCheckboxInput(session, "cmp_on1", value = TRUE)
    updateCheckboxInput(session, "cmp_on2", value = TRUE)
    updateCheckboxInput(session, "cmp_on3", value = TRUE)
    updateCheckboxInput(session, "cmp_on4", value = FALSE)
  })
  
  output$comparePlot <- renderPlot({
    sfobj <- data_sf()
    settings <- list()
    
    for (i in 1:4) {
      on <- isTRUE(input[[paste0("cmp_on", i)]])
      co <- input[[paste0("cmp_cutoff", i)]]
      wi <- input[[paste0("cmp_width", i)]]
      if (on && !is.null(co) && !is.null(wi) && !is.na(co) && !is.na(wi) &&
          co > 0 && wi > 0) {
        v <- tryCatch(variogram(z ~ 1, sfobj, cutoff = co, width = wi),
                      error = function(e) NULL)
        if (!is.null(v) && nrow(v) > 0) {
          v$setting <- paste0("Setting ", i, "  (cutoff=", round(co, 1),
                              ", width=", round(wi, 1), ")")
          settings[[length(settings) + 1]] <- v
        }
      }
    }
    
    validate(need(length(settings) > 0,
                  "Tick at least one setting and enter valid cutoff/width values (> 0)."))
    
    combined <- do.call(rbind, settings)
    ggplot(combined, aes(x = dist, y = gamma, color = setting)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.4) +
      labs(title = "Comparing Binned Variograms Across Settings",
           x = "Distance, h", y = expression(gamma(h)), color = "Binning setting") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "bottom")
  })
  
  # ---- Optimal bin-count search (single-setting mode only) -----------------
  optimal_results <- eventReactive(input$findOptimal, {
    sfobj <- data_sf()
    cutoff_fixed <- input$cutoff
    nbin_range <- 5:30
    
    out <- data.frame(nbins = integer(), width = numeric(),
                      SSErr = numeric(), min_np = integer())
    
    for (n in nbin_range) {
      w <- cutoff_fixed / n
      v <- tryCatch(
        variogram(z ~ 1, sfobj, cutoff = cutoff_fixed, width = w),
        error = function(e) NULL
      )
      if (is.null(v) || nrow(v) < 3) next
      
      min_np <- min(v$np)
      psill_guess <- max(v$gamma, na.rm = TRUE)
      range_guess  <- cutoff_fixed / 2
      nugget_guess <- min(v$gamma, na.rm = TRUE)
      
      fit <- tryCatch(
        fit.variogram(v, vgm(psill_guess, "Sph", range_guess, nugget_guess)),
        error = function(e) NULL
      )
      sserr <- if (!is.null(fit)) attr(fit, "SSErr") else NA
      
      out <- rbind(out, data.frame(nbins = n, width = w,
                                   SSErr = sserr, min_np = min_np))
    }
    out
  })
  
  best_row <- reactive({
    res <- optimal_results()
    req(nrow(res) > 0)
    valid <- res[!is.na(res$SSErr) & res$min_np >= 30, ]
    if (nrow(valid) == 0) valid <- res[!is.na(res$SSErr), ]
    validate(need(nrow(valid) > 0, "No valid model fit found across bin counts tried."))
    valid[which.min(valid$SSErr), ]
  })
  
  output$optimalPlot <- renderPlot({
    res <- optimal_results()
    req(nrow(res) > 0)
    best <- best_row()
    ggplot(res, aes(x = nbins, y = SSErr)) +
      geom_line(color = "gray30") +
      geom_point(color = "gray30") +
      geom_vline(xintercept = best$nbins, color = "red", linetype = "dashed") +
      annotate("point", x = best$nbins, y = best$SSErr, color = "red", size = 3) +
      labs(title = "Model-fit error vs. number of bins",
           subtitle = "Lower SSErr = better fit of a theoretical model to the binned variogram",
           x = "Number of bins", y = "Weighted SSE (model fit)") +
      theme_minimal(base_size = 13)
  })
  
  output$optimalText <- renderText({
    best <- best_row()
    paste0(
      "Suggested optimal binning: ", best$nbins, " bins  \u2192  width \u2248 ",
      round(best$width, 1), "  (cutoff fixed at ", round(input$cutoff, 1),
      "). This minimizes the weighted SSE between the fitted spherical ",
      "model and the binned variogram, while keeping at least ~30 pairs ",
      "in the smallest bin."
    )
  })
  
  output$optimalTable <- renderTable({
    res <- optimal_results()
    req(nrow(res) > 0)
    res %>%
      arrange(SSErr) %>%
      head(5) %>%
      mutate(width = round(width, 1), SSErr = signif(SSErr, 4)) %>%
      rename(`# bins` = nbins, `width` = width,
             `SSErr` = SSErr, `min pairs/bin` = min_np)
  })
  
  observeEvent(input$applyOptimal, {
    best <- best_row()
    updateNumericInput(session, "width", value = round(best$width, 2))
  })
  
  # ---- Tab 3: Step-by-step binning visualization --------------------------
  step_variogram <- reactive({
    req(input$cutoffStep, input$widthStep)
    variogram(z ~ 1, data_sf(), cutoff = input$cutoffStep, width = input$widthStep)
  })
  
  output$stepPlot <- renderPlot({
    req(input$cutoffStep, input$widthStep)
    vc <- cloud_data()
    v  <- step_variogram()
    breaks <- seq(0, input$cutoffStep, by = input$widthStep)
    
    p <- ggplot() +
      geom_point(data = vc, aes(x = dist, y = gamma),
                 color = "blue", alpha = 0.25, size = 1) +
      xlim(0, input$cutoffStep * 1.05) +
      labs(x = "Distance, h", y = expression(gamma(h))) +
      theme_minimal(base_size = 14)
    
    if (input$step >= 2) {
      p <- p + geom_vline(xintercept = breaks, linetype = "dashed",
                          color = "gray40", linewidth = 0.4)
    }
    if (input$step >= 3) {
      p <- p +
        geom_point(data = v, aes(x = dist, y = gamma),
                   shape = 1, color = "red", size = 3.5, stroke = 1.2) +
        geom_text(data = v, aes(x = dist, y = gamma, label = np),
                  hjust = -0.35, vjust = -0.35, size = 3.3, color = "black")
    }
    if (input$step >= 4) {
      p <- p + geom_line(data = v, aes(x = dist, y = gamma),
                         color = "red", linewidth = 1)
    }
    
    title <- switch(as.character(input$step),
                    "1" = "Step 1 — Variogram cloud (every point-pair's distance & semi-variance)",
                    "2" = "Step 2 — Points grouped into distance bins (dashed lines = bin edges)",
                    "3" = "Step 3 — Average gamma computed within each bin (red points)",
                    "4" = "Step 4 — Bin averages connected \u2192 the Binned Empirical Variogram"
    )
    p + ggtitle(title)
  })
}

shinyApp(ui, server)
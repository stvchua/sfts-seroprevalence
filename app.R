# app.R
library(shiny)
library(readxl)
library(dplyr)
library(ggplot2)
library(janitor)
library(DT)
library(sf)
library(leaflet)
library(plotly)
library(jpinfect)

# ---- Data loading & cleaning ----

db_path <- "data/sfts_seroprevalence_japan.xlsx"

sero_raw <- readxl::read_excel(db_path, sheet = "sero_estimates") |>
  janitor::clean_names()

sero <- sero_raw |>
  mutate(
    n_tested   = as.numeric(n_tested),
    n_positive = as.numeric(n_positive),
    prev_prop_calc = if_else(
      !is.na(n_tested) & !is.na(n_positive) & n_tested > 0,
      n_positive / n_tested,
      NA_real_
    ),
    prev_percent = 100 * prev_prop_calc
  ) |>
  filter(country == "Japan")

available_hosts <- sort(unique(sero$host_type))
year_min <- min(sero$sampling_year_start, na.rm = TRUE)
year_max <- max(sero$sampling_year_end,   na.rm = TRUE)

# ---- Prefecture shapes ----
pref_sf <- sf::st_read("data/japan_prefectures.geojson", quiet = TRUE) |>
  janitor::clean_names() |>
  rename(admin1 = name) |>
  mutate(
    admin1 = if_else(admin1 == "Naoasaki", "Nagasaki", admin1)
  )

# ---- Weekly SFTS data ----

weekly_path_national <- "data/sfts_weekly_national.csv"
weekly_path_pref     <- "data/sfts_weekly_all_prefectures.csv"

sfts_weekly_national <- if (file.exists(weekly_path_national)) {
  read.csv(weekly_path_national) |>
    mutate(date = as.Date(date))
} else {
  data.frame(date = as.Date(character()), cases = integer(), source = character())
}

sfts_weekly_pref <- if (file.exists(weekly_path_pref)) {
  read.csv(weekly_path_pref) |>
    mutate(date = as.Date(date))
} else {
  data.frame(date = as.Date(character()), prefecture = character(),
             cases = integer(), source = character())
}

# FIX 1: use as.integer() so these are numeric from the start, not character
weekly_year_min <- if (nrow(sfts_weekly_national) > 0) {
  as.integer(min(format(sfts_weekly_national$date, "%Y")))
} else { 2013L }

weekly_year_max <- if (nrow(sfts_weekly_national) > 0) {
  as.integer(max(format(sfts_weekly_national$date, "%Y")))
} else { as.integer(format(Sys.Date(), "%Y")) }

# ---- UI ----

ui <- fluidPage(
  titlePanel("SFTSV seroprevalence by prefecture (Japan)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      checkboxGroupInput(
        "host_type", "Host type",
        choices  = available_hosts,
        selected = available_hosts
      ),
      uiOutput("species_ui"),
      sliderInput(
        "year_range", "Sampling years (start-end):",
        min   = year_min,
        max   = year_max,
        value = c(year_min, year_max),
        step  = 1, sep = ""
      ),
      checkboxInput(
        "only_with_counts",
        "Include only rows with valid counts (n_tested > 0)",
        value = TRUE
      ),
      helpText("Pooled prevalence = total positives / total tested per prefecture.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Prefecture plot",  br(), plotlyOutput("prefecture_plot", height = "500px")),
        tabPanel("Prefecture table", br(), DTOutput("prefecture_table")),
        tabPanel("Map",              br(), leafletOutput("prefecture_map", height = "600px")),
        tabPanel("Weekly cases", br(),
                 fluidRow(
                   column(4,
                          # FIX 2: guard against min == max (only one year of data) so
                          # sliderInput doesn't break; clamp max to at least min + 1
                          sliderInput("weekly_years", "Year range",
                                      min   = weekly_year_min,
                                      max   = max(weekly_year_max, weekly_year_min + 1L),
                                      value = c(weekly_year_min, weekly_year_max),
                                      step  = 1, sep = ""
                          )
                   ),
                   column(4,
                          selectInput("weekly_source", "Data source",
                                      choices  = c("All", "confirmed", "provisional"),
                                      selected = "All"
                          )
                   ),
                   column(4,
                          checkboxInput("weekly_smooth", "Show 4-week rolling average", value = FALSE)
                   )
                 ),
                 plotlyOutput("weekly_plot", height = "420px"),
                 br(),
                 DTOutput("weekly_table")
        )                          # closes tabPanel("Weekly cases")
      )                            # closes tabsetPanel
    )                              # closes mainPanel
  )                                # closes sidebarLayout
)                                  # closes fluidPage

# ---- Server ----

server <- function(input, output, session) {
  
  # Species selector
  output$species_ui <- renderUI({
    dat <- sero
    if (!is.null(input$host_type) && length(input$host_type) > 0) {
      dat <- dat |> filter(host_type %in% input$host_type)
    }
    species_choices <- sort(unique(dat$species_common))
    selectInput("species", "Species",
                choices  = c("All", species_choices),
                selected = "All"
    )
  })
  
  # Base filtered data
  filtered_data <- reactive({
    dat <- sero
    
    if (!is.null(input$host_type) && length(input$host_type) > 0) {
      dat <- dat |> filter(host_type %in% input$host_type)
    }
    if (!is.null(input$species) && input$species != "All") {
      dat <- dat |> filter(species_common == input$species)
    }
    
    yr_min <- input$year_range[1]
    yr_max <- input$year_range[2]
    dat <- dat |>
      filter(sampling_year_start <= yr_max, sampling_year_end >= yr_min)
    
    if (isTRUE(input$only_with_counts)) {
      dat <- dat |> filter(!is.na(n_tested), !is.na(n_positive), n_tested > 0)
    }
    
    dat
  })
  
  # ---- Prefecture summaries ----
  
  prefecture_summary_overall <- reactive({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    dat |>
      group_by(admin1) |>
      summarise(
        n_estimates         = n_distinct(estimate_id),
        n_refs              = n_distinct(ref_id),
        total_tested        = sum(n_tested, na.rm = TRUE),
        total_positive      = sum(n_positive, na.rm = TRUE),
        pooled_prev_percent = if_else(
          total_tested > 0, 100 * total_positive / total_tested, NA_real_
        ),
        first_year = min(sampling_year_start, na.rm = TRUE),
        last_year  = max(sampling_year_end,   na.rm = TRUE),
        .groups = "drop"
      ) |>
      arrange(desc(pooled_prev_percent))
  })
  
  prefecture_summary_by_host <- reactive({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    dat |>
      group_by(admin1, host_type) |>
      summarise(
        n_estimates         = n(),
        total_tested        = sum(n_tested, na.rm = TRUE),
        total_positive      = sum(n_positive, na.rm = TRUE),
        pooled_prev_percent = if_else(
          total_tested > 0, 100 * total_positive / total_tested, NA_real_
        ),
        .groups = "drop"
      ) |>
      arrange(admin1, host_type)
  })
  
  # ---- Prefecture plot ----
  
  output$prefecture_plot <- renderPlotly({
    dat <- prefecture_summary_by_host()
    req(nrow(dat) > 0)
    
    p <- ggplot(dat, aes(
      x    = pooled_prev_percent,
      y    = reorder(admin1, pooled_prev_percent, FUN = max),
      fill = host_type,
      text = paste0(
        "Prefecture: ", admin1, "<br>",
        "Host type: ", host_type, "<br>",
        "Pooled prev: ", sprintf("%.1f%%", pooled_prev_percent), "<br>",
        "Total tested: ", total_tested, "<br>",
        "Total positive: ", total_positive
      )
    )) +
      geom_col(width = 0.7, position = "stack") +
      scale_x_continuous(limits = c(0, 100)) +
      labs(x = "Pooled seroprevalence (%)", y = "Prefecture",
           fill = "Host type",
           title = "SFTSV pooled seroprevalence by prefecture and host type") +
      theme_minimal(base_size = 13)
    
    ggplotly(p, tooltip = "text")
  })
  
  # ---- Prefecture table ----
  
  output$prefecture_table <- renderDT({
    dat <- prefecture_summary_overall() |>
      mutate(pooled_prev_percent = round(pooled_prev_percent, 1))
    
    datatable(dat, options = list(pageLength = 20), rownames = FALSE)
  })
  
  # ---- Map data ----
  
  prefecture_map_data <- reactive({
    summary <- prefecture_summary_overall()
    
    pref_sf |>
      left_join(summary, by = "admin1") |>
      mutate(
        prev_cat = case_when(
          is.na(pooled_prev_percent)                             ~ "No data",
          pooled_prev_percent == 0                               ~ "0%",
          pooled_prev_percent > 0  & pooled_prev_percent <= 5   ~ "1-5%",
          pooled_prev_percent > 5  & pooled_prev_percent <= 10  ~ "5-10%",
          pooled_prev_percent > 10 & pooled_prev_percent <= 25  ~ "10-25%",
          pooled_prev_percent > 25 & pooled_prev_percent <= 50  ~ "25-50%",
          pooled_prev_percent > 50                              ~ ">50%"
        ),
        prev_cat = factor(prev_cat,
                          levels = c("0%","1-5%","5-10%","10-25%","25-50%",">50%","No data"))
      )
  })
  
  # ---- Map ----
  
  output$prefecture_map <- renderLeaflet({
    shp <- prefecture_map_data()
    req(nrow(shp) > 0)
    
    cat_cols <- c(
      "No data" = "#FFFFFF", "0%"   = "#B0B0B0", "1-5%"  = "#FFF7BC",
      "5-10%"   = "#FEE391", "10-25%" = "#FEC44F",
      "25-50%"  = "#FD8D3C", ">50%" = "#E31A1C"
    )
    
    shp$fill_col <- cat_cols[as.character(shp$prev_cat)]
    
    leaflet(shp) |>
      addProviderTiles("CartoDB.Positron") |>
      addPolygons(
        fillColor = ~fill_col, color = "#444444",
        weight = 1, opacity = 1, fillOpacity = 0.7,
        highlightOptions = highlightOptions(weight = 2, color = "#000000", bringToFront = TRUE),
        label = ~paste0(admin1, ": ",
                        ifelse(is.na(pooled_prev_percent), "no data",
                               sprintf("%.1f%%", pooled_prev_percent))),
        popup = ~paste0(
          "<strong>", admin1, "</strong><br/>",
          "Pooled seroprevalence: ",
          ifelse(is.na(pooled_prev_percent), "No data",
                 sprintf("%.1f%%", pooled_prev_percent)), "<br/>",
          "Total tested: ",   ifelse(is.na(total_tested),   "NA", total_tested),   "<br/>",
          "Total positive: ", ifelse(is.na(total_positive), "NA", total_positive), "<br/>",
          "Studies (papers): ",          ifelse(is.na(n_refs),      "NA", n_refs),      "<br/>",
          "Estimates (species/hosts): ", ifelse(is.na(n_estimates), "NA", n_estimates), "<br/>",
          "Data years: ",
          ifelse(is.na(first_year) | is.na(last_year), "NA",
                 paste0(first_year, "-", last_year))
        )
      ) |>
      addLegend(
        position = "bottomright", title = "Pooled seroprevalence",
        colors  = cat_cols[c("No data","0%","1-5%","5-10%","10-25%","25-50%",">50%")],
        labels  = c("No data","0%","1-5%","5-10%","10-25%","25-50%",">50%"),
        opacity = 1
      )
  })
  
  # ---- Weekly cases ----
  
  weekly_filtered <- reactive({
    dat <- sfts_weekly_national
    req(nrow(dat) > 0)
    
    # FIX 3: cast both sides to integer so year comparison is numeric not
    # character — avoids lexicographic ordering bugs ("9" > "10" etc.)
    dat <- dat |>
      filter(
        as.integer(format(date, "%Y")) >= input$weekly_years[1],
        as.integer(format(date, "%Y")) <= input$weekly_years[2]
      )
    
    if (!is.null(input$weekly_source) && input$weekly_source != "All") {
      dat <- dat |> filter(source == input$weekly_source)
    }
    
    dat
  })
  
  output$weekly_plot <- renderPlotly({
    dat <- weekly_filtered()
    req(nrow(dat) > 0)
    
    if (isTRUE(input$weekly_smooth)) {
      dat <- dat |>
        arrange(date) |>
        mutate(cases_smooth = stats::filter(cases, rep(1/4, 4), sides = 1))
    }
    
    p <- ggplot(dat, aes(x = date)) +
      geom_col(aes(y = cases, fill = source), width = 5, alpha = 0.7) +
      scale_fill_manual(
        values   = c(confirmed = "#2166AC", provisional = "#F4A582"),
        na.value = "#888888"
      ) +
      labs(x = "Week", y = "Cases",
           title = "SFTS weekly reported cases — Japan (national total)",
           fill  = "Source") +
      theme_minimal(base_size = 13)
    
    if (isTRUE(input$weekly_smooth)) {
      p <- p + geom_line(aes(y = cases_smooth), colour = "#B2182B",
                         linewidth = 0.8, na.rm = TRUE)
    }
    
    ggplotly(p, tooltip = c("x", "y", "fill"))
  })
  
  output$weekly_table <- renderDT({
    weekly_filtered() |>
      arrange(desc(date)) |>
      mutate(date = format(date, "%Y-%m-%d")) |>
      datatable(options = list(pageLength = 15), rownames = FALSE)
  })
  
}                                  # closes server

shinyApp(ui, server)

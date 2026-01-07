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

available_hosts   <- sort(unique(sero$host_type))
year_min <- min(sero$sampling_year_start, na.rm = TRUE)
year_max <- max(sero$sampling_year_end,   na.rm = TRUE)

# ---- Prefecture shapes ----
# Assumes geojson has an "admin1" column with prefecture names in English.
pref_sf <- sf::st_read("data/japan_prefectures.geojson", quiet = TRUE) |>
  janitor::clean_names() |>
  rename(admin1 = name) |>
  mutate(
    # Fix typo: "Naoasaki" -> "Nagasaki"
    admin1 = if_else(admin1 == "Naoasaki", "Nagasaki", admin1)
  )

# ---- UI ----

ui <- fluidPage(
  titlePanel("SFTSV seroprevalence by prefecture (Japan)"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      checkboxGroupInput(
        "host_type",
        "Host type",
        choices  = available_hosts,
        selected = available_hosts
      ),
      uiOutput("species_ui"),
      sliderInput(
        "year_range",
        "Sampling years (start–end):",
        min   = year_min,
        max   = year_max,
        value = c(year_min, year_max),
        step  = 1,
        sep   = ""
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
        tabPanel("Raw estimates",     br(), DTOutput("raw_table")),
        tabPanel("Map",              br(), leafletOutput("prefecture_map", height = "600px"))
      )
    )
  )
)

# ---- Server ----

server <- function(input, output, session) {
  
  # Species selector depends on current host-type filter
  output$species_ui <- renderUI({
    dat <- sero
    if (!is.null(input$host_type) && length(input$host_type) > 0) {
      dat <- dat |> filter(host_type %in% input$host_type)
    }
    species_choices <- sort(unique(dat$species_common))
    selectInput(
      "species",
      "Species (optional)",
      choices  = c("All", species_choices),
      selected = "All"
    )
  })
  
  # Base filtered data
  filtered_data <- reactive({
    dat <- sero
    
    # Host type filter
    if (!is.null(input$host_type) && length(input$host_type) > 0) {
      dat <- dat |> filter(host_type %in% input$host_type)
    }
    
    # Species filter
    if (!is.null(input$species) && input$species != "All") {
      dat <- dat |> filter(species_common == input$species)
    }
    
    # Year window: any overlap between study and slider
    yr_min <- input$year_range[1]
    yr_max <- input$year_range[2]
    dat <- dat |>
      filter(
        sampling_year_start <= yr_max,
        sampling_year_end   >= yr_min
      )
    
    # Valid counts only?
    if (isTRUE(input$only_with_counts)) {
      dat <- dat |>
        filter(!is.na(n_tested), !is.na(n_positive), n_tested > 0)
    }
    
    dat
  })
  
  # ---- Prefecture summaries ----
  
  # Overall pooled per prefecture (for table + map)
  prefecture_summary_overall <- reactive({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    dat |>
      group_by(admin1) |>
      summarise(
        n_estimates     = n_distinct(estimate_id),  # number of unique estimates
        n_refs          = n_distinct(ref_id),       # number of unique papers (studies)
        total_tested    = sum(n_tested, na.rm = TRUE),
        total_positive  = sum(n_positive, na.rm = TRUE),
        pooled_prev_percent = if_else(
          total_tested > 0,
          100 * total_positive / total_tested,
          NA_real_
        ),
        first_year = min(sampling_year_start, na.rm = TRUE),
        last_year  = max(sampling_year_end,   na.rm = TRUE),
        .groups = "drop"
      ) |>
      arrange(desc(pooled_prev_percent))
  })
  
  # Pooled per prefecture *by host type* (for coloured stacked bar plot)
  prefecture_summary_by_host <- reactive({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    dat |>
      group_by(admin1, host_type) |>
      summarise(
        n_estimates     = n(),
        total_tested    = sum(n_tested, na.rm = TRUE),
        total_positive  = sum(n_positive, na.rm = TRUE),
        pooled_prev_percent = if_else(
          total_tested > 0,
          100 * total_positive / total_tested,
          NA_real_
        ),
        .groups = "drop"
      ) |>
      arrange(admin1, host_type)
  })
  
  # ---- Prefecture plot: stacked by host_type, 0–100% axis ----
  
  output$prefecture_plot <- renderPlotly({
    dat <- prefecture_summary_by_host()
    req(nrow(dat) > 0)
    
    # ggplot stacked bar (one bar per prefecture, coloured by host type)
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
      labs(
        x = "Pooled seroprevalence (%)",
        y = "Prefecture",
        fill = "Host type",
        title = "SFTSV pooled seroprevalence by prefecture and host type (stacked)"
      ) +
      theme_minimal(base_size = 13)
    
    # Convert to interactive plotly object; tooltip uses the `text` aesthetic
    ggplotly(p, tooltip = "text")
  })
  
  # ---- Prefecture table: pooled_prev_percent rounded to 1 decimal ----
  
  output$prefecture_table <- renderDT({
    dat <- prefecture_summary_overall() |>
      mutate(
        pooled_prev_percent = round(pooled_prev_percent, 1)
      )
    
    datatable(
      dat,
      options = list(pageLength = 20),
      rownames = FALSE
    )
  })
  
  # ---- Raw estimates table ----
  
  output$raw_table <- renderDT({
    dat <- filtered_data()
    
    datatable(
      dat |>
        select(
          estimate_id, ref_id, country, admin1, admin2,
          host_type, species_common, population_group,
          n_tested, n_positive,
          prev_prop,         # original from Excel (percent)
          prev_percent,      # recomputed from counts
          sampling_year_start, sampling_year_end,
          assay_type, outcome_detail
        ),
      options = list(pageLength = 20),
      rownames = FALSE
    )
  })
  
  # ---- Map data: join summary to shapes, create prevalence bins ----
  
  prefecture_map_data <- reactive({
    summary <- prefecture_summary_overall()
    
    shp <- pref_sf |>
      left_join(summary, by = "admin1")
    
    shp |>
      mutate(
        prev_cat = case_when(
          is.na(pooled_prev_percent)               ~ "No data",
          pooled_prev_percent == 0                 ~ "0%",
          pooled_prev_percent > 0 & pooled_prev_percent <= 5   ~ "1–5%",
          pooled_prev_percent > 5 & pooled_prev_percent <= 10  ~ "5–10%",
          pooled_prev_percent > 10 & pooled_prev_percent <= 25 ~ "10–25%",
          pooled_prev_percent > 25 & pooled_prev_percent <= 50 ~ "25–50%",
          pooled_prev_percent > 50                 ~ ">50%"
        ),
        prev_cat = factor(
          prev_cat,
          levels = c("0%", "1–5%", "5–10%", "10–25%", "25–50%", ">50%", "No data")
        )
      )
  })
  
  # ---- Map: white for NA, grey for 0%, beige→yellow→orange→red for >0 ----
  
  output$prefecture_map <- renderLeaflet({
    shp <- prefecture_map_data()
    req(nrow(shp) > 0)
    
    # Explicit colours
    cat_cols <- c(
      "No data" = "#FFFFFF",  # white
      "0%"      = "#B0B0B0",  # grey
      "1–5%"    = "#FFF7BC",  # very light beige/yellow
      "5–10%"   = "#FEE391",  # pale yellow
      "10–25%"  = "#FEC44F",  # yellow-orange
      "25–50%"  = "#FD8D3C",  # orange
      ">50%"    = "#E31A1C"   # red
    )
    
    shp$fill_col <- cat_cols[as.character(shp$prev_cat)]
    
    leaflet(shp) |>
      addProviderTiles("CartoDB.Positron") |>
      addPolygons(
        fillColor   = ~fill_col,
        color       = "#444444",
        weight      = 1,
        opacity     = 1,
        fillOpacity = 0.7,
        highlightOptions = highlightOptions(
          weight = 2,
          color = "#000000",
          bringToFront = TRUE
        ),
        label = ~paste0(
          admin1, ": ",
          ifelse(is.na(pooled_prev_percent),
                 "no data",
                 sprintf("%.1f%%", pooled_prev_percent))
        ),
        popup = ~paste0(
          "<strong>", admin1, "</strong><br/>",
          "Pooled seroprevalence: ",
          ifelse(is.na(pooled_prev_percent),
                 "No data",
                 sprintf("%.1f%%", pooled_prev_percent)), "<br/>",
          "Total tested: ", ifelse(is.na(total_tested), "NA", total_tested), "<br/>",
          "Total positive: ", ifelse(is.na(total_positive), "NA", total_positive), "<br/>",
          "Studies (papers): ", ifelse(is.na(n_refs), "NA", n_refs), "<br/>",
          "Estimates (species/hosts): ", ifelse(is.na(n_estimates), "NA", n_estimates), "<br/>",
          "Data years: ",
          ifelse(is.na(first_year) | is.na(last_year),
                 "NA",
                 paste0(first_year, "–", last_year))
        )
      ) |>
      addLegend(
        position = "bottomright",
        title    = "Pooled seroprevalence",
        colors   = cat_cols[c("No data", "0%", "1–5%", "5–10%", "10–25%", "25–50%", ">50%")],
        labels   = c("No data", "0%", "1–5%", "5–10%", "10–25%", "25–50%", ">50%"),
        opacity  = 1
      )
  })
}

shinyApp(ui, server)
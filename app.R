# app.R
library(shiny)
library(readxl)
library(dplyr)
library(ggplot2)
library(janitor)
library(DT)

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
        tabPanel("Prefecture plot",  br(), plotOutput("prefecture_plot", height = "500px")),
        tabPanel("Prefecture table", br(), DTOutput("prefecture_table")),
        tabPanel("Raw estimates",     br(), DTOutput("raw_table"))
      )
    )
  )
)

# ---- Server ----

server <- function(input, output, session) {
  
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
      filter(
        sampling_year_start <= yr_max,
        sampling_year_end   >= yr_min
      )
    
    if (isTRUE(input$only_with_counts)) {
      dat <- dat |>
        filter(!is.na(n_tested), !is.na(n_positive), n_tested > 0)
    }
    
    dat
  })
  
  prefecture_summary <- reactive({
    dat <- filtered_data()
    req(nrow(dat) > 0)
    
    dat |>
      group_by(admin1) |>
      summarise(
        n_estimates   = n(),
        total_tested  = sum(n_tested, na.rm = TRUE),
        total_positive = sum(n_positive, na.rm = TRUE),
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
  
  output$prefecture_plot <- renderPlot({
    dat <- prefecture_summary()
    req(nrow(dat) > 0)
    
    ggplot(dat, aes(x = reorder(admin1, pooled_prev_percent),
                    y = pooled_prev_percent)) +
      geom_col() +
      coord_flip() +
      labs(
        x = "Prefecture",
        y = "Pooled seroprevalence (%)",
        title = "SFTSV pooled seroprevalence by prefecture"
      ) +
      theme_minimal(base_size = 13)
  })
  
  output$prefecture_table <- renderDT({
    dat <- prefecture_summary()
    datatable(dat, options = list(pageLength = 20), rownames = FALSE)
  })
  
  output$raw_table <- renderDT({
    dat <- filtered_data()
    datatable(
      dat |>
        select(
          estimate_id, ref_id, country, admin1, admin2,
          host_type, species_common, population_group,
          n_tested, n_positive,
          prev_prop,         # original column from Excel
          prev_percent,      # recomputed from counts
          sampling_year_start, sampling_year_end,
          assay_type, outcome_detail
        ),
      options = list(pageLength = 20),
      rownames = FALSE
    )
  })
}

shinyApp(ui, server)

server <- function(input, output, session) {

  ############################################################
  #Homepage
output$welcome_text <- renderUI({

  fluidPage(

    # HERO

div(

  class = "landing-hero",

  div(

    class = "landing-brand",

    tags$img(
      src = "logo.png",
      style = "
        height:90px;
        width:auto;
        flex-shrink:0;
      "
    ),

    div(

      style="
        display:flex;
        flex-direction:column;
        justify-content:center;
        align-items:flex-start;
        line-height:1;
      ",

      div(

        HTML(
          "
          <span style='
            color:#0f172a;
            font-weight:700;
            font-size:52px;
          '>TPP</span>

          <span style='
            color:#0EA5A5;
            font-weight:300;
            font-size:42px;
          '>Analyst</span>
          "
        )

      ),

      div(
        class = "landing-divider"
      ),

      div(
        "THERMAL PROTEOMICS PLATFORM",
        class = "landing-subtitle"
      )

    )

  ),

  p(
    "Comprehensive analysis of Thermal Proteome Profiling datasets from quality control to biological interpretation.",

    style="
      font-size:13px;
      color:#64748b;
      max-width:700px;
      margin:10px auto 0 auto;
      line-height:1.4;
    "
  )

),

div(

  class = "manual-card",

  style = "
    text-align:center;
    margin:15px auto 25px auto;
    max-width:500px;
  ",

  icon(
    "book-open",
    style="
      font-size:22px;
      color:#0EA5A5;
      margin-bottom:8px;
    "
  ),

  br(),

  tags$strong(
    style="
      font-size:18px;
      display:block;
      margin-bottom:6px;
    ",
    "Documentation & User Guide"
  ),

  tags$a(
    href = "https://github.com/knth16/TPP_Analyst",
    target = "_blank",

    "Open GitHub documentation",

    style="
      color:#0EA5A5;
      font-weight:500;
      text-decoration:none;
    "
  )

),

    # WORKFLOW

fluidRow(

  column(
    12,

    div(

      class = "workflow-row",


      div(
        class="workflow-node",
        icon("database"),
        span("Data Input")
      ),

      div(class="workflow-arrow", icon("angle-right")),

      div(
        class="workflow-node",
        icon("search"),
        span("Quality Control")
      ),

      div(class="workflow-arrow", icon("angle-right")),

      div(
        class="workflow-node",
        icon("sliders"),
        span("Normalization")
      ),

      div(class="workflow-arrow", icon("angle-right")),

      div(
        class="workflow-node",
        icon("diagram-project"),
        span("Clustering")
      ),

      div(class="workflow-arrow", icon("angle-right")),

      div(
        class="workflow-node",
        icon("chart-simple"),
        span("Model Fitting")
      ),

      div(class="workflow-arrow", icon("angle-right")),

      div(
        class="workflow-node",
        icon("chart-column"),
        span("Protein Plots")
      ),

      div(class="workflow-arrow", icon("angle-right")),

      div(
        class="workflow-node",
        icon("dna"),
        span("GO Terms")
      ),

      div(class="workflow-arrow", icon("angle-right")),

      div(
        class="workflow-node",
        icon("scale-balanced"),
        span("Comparison Models")
      )

    )
  )
),

    br(),

    fluidRow(

      column(
        4,

        div(
          class = "landing-card",

          img(
            src = "heatmap.png",
            width = "100%"
          ),

          h4("Quality control"),

          p(
            "Identify missing data, evaluate replicate consistency,
and assess dataset quality before downstream analysis."
          )
        )
      ),

      column(
        4,

        div(
          class = "landing-card",

          img(
            src = "clustering.png",
            width = "100%"
          ),

          h4("Clustering analysis"),

          p(
            "Group proteins with similar thermal stability profiles
to identify coordinated responses to treatment."
          )
        )
      ),

      column(
        4,

        div(
          class = "landing-card",

          img(
            src = "TPP.png",
            width = "100%"
          ),

          h4("Comparison models"),

          p(
            "Detect thermal stability shifts and compare complete
melting curves between experimental conditions."
          )
        )
      )

    )
  )
})

  ############################################################
  #File input and data display 
  # Reactive data loading
  # Load datasets via modules
#example dataset reactive values
example_data <- reactiveValues(
  proteins = NULL,
  design = NULL,
  loaded = FALSE
)
using_example_data <- reactive({

  example_data$loaded &&
    is.null(input$file2)

})
design_source <- reactive({

  if (using_example_data()) {

    return(example_data$design)

  }

  NULL

})




protein_source <- reactive({

  if (using_example_data()){

    return(example_data$proteins)

  }

  NULL

})
#example data buttom
observeEvent(input$load_example, {

  example_data$design <- read.csv(
    "Example_Data/experimental_design.csv",
    check.names = FALSE
  )

  example_data$proteins <- read.csv(
    "Example_Data/protein_matrix.csv",
    check.names = FALSE
  )

  example_data$loaded <- TRUE

  # Automatically set preprocessing flags for example dataset
updateCheckboxInput(
  session,
  "already_log2",
  value = TRUE
)

updateCheckboxInput(
  session,
  "Already_lowestnorm",
  value = TRUE
)

})

observeEvent(using_example_data(), {

  if (!using_example_data()) {

    updateCheckboxInput(session, "already_log2", value = FALSE)

    updateCheckboxInput(session, "Already_lowestnorm", value = FALSE)

  }

})

#mapping for example dataset
example_mapping <- reactive({

  if (!example_data$loaded) {
    return(NULL)
  }

  list(
    "Vehicle"   = "H.sapiens",
    "Treatment" = "H.sapiens"
  )

})


dataset1 <- proteindataModuleServer(
  "dataset1",
  fileInput = reactive(input$file),
  requireOrganismColumn = require_protein_organism_column,
  exampleData = protein_source
)

selected_organisms <- reactive({

  mapping <- condition_organism_mapping()

  mapping <- mapping[
    sapply(
      mapping,
      function(x)
        !is.null(x) &&
        length(x) > 0
    )
  ]

  mapping

})
dataset2 <- designdataModuleServer(
  "dataset2",
  fileInput = reactive(input$file2),
  proteinInput = dataset1$data,
  customOrgData = custom_orgdata_files,
  selected_organisms = selected_organisms,
  exampleData = design_source
)
#debugs
observe({

  if (!is.null(dataset1$data())) {

    cat("\nEXAMPLE PROTEINS AVAILABLE\n")

    print(dim(dataset1$data()))

  }

})

observe({

  if (!is.null(dataset2$data())) {

    cat("\nEXAMPLE DESIGN AVAILABLE\n")

    print(dim(dataset2$data()))

    cat("\nDESIGN COLUMN NAMES\n")

    print(
      names(
        dataset2$data()
      )
    )

    cat("\nFIRST 10 DESIGN ROWS\n")

    print(
      head(
        dataset2$data(),
        10
      )
    )

    write.csv(
      dataset2$data(),
      "design_table.csv",
      row.names = FALSE
    )

  }

})


#join datasets


processed <- processedDataModuleServer(
  "processed",
  proteinInput = dataset1$data,
  designInput  = dataset2$data,
  splitMode    = reactive(input$split_mode),
  organismMapping = selected_organisms,
  customProteinMapping = custom_protein_mapping
)

available_organisms <- c(
  "Sulfolobus acidocaldarius",
  "Parageobacillus thermoglucosidasius",
  "Caldimonas thermodepolymerans",
  "Haloferax volcanii",
  "Thermosynechococcus elongatus",
  "Thermus thermophilus",
  "E.coli",
  "H.sapiens",
  "Other"
)

output$condition_organism_mapping_ui <- renderUI({

  req(dataset2$data())

  design <- dataset2$data()

  req(design)

  conditions <- unique(
    as.character(design$Condition)
  )

  tagList(


  tags$label(

    class = "control-label",

    "Condition to organism mapping ",

    tags$span(
      `data-toggle` = "tooltip",
      `data-placement` = "right",
      `data-html` = "true",

      title = HTML(
        "<b>Condition to organism mapping</b><br><br>

        Select the organism(s) used in each condition.<br><br>

        • Multiple organisms per condition are allowed.<br>

        • If the organism is not listed, select <b>Other</b>.<br>

        • Additional organism databases can be uploaded afterwards."
      ),

      tags$i(
        class = "fa fa-question-circle",
        style = "cursor:pointer;"
      )
    )
  )
,


    lapply(conditions, function(cond){

selectizeInput(
  inputId = paste0(
    "organism_map_",
    make.names(cond)
  ),
  

 label = tags$span(
  style = "
    font-weight: normal;
    color: #333;
    font-size: 14px;
  ",
  cond
),
        choices = available_organisms,
        multiple = TRUE,
        selected = if (
  using_example_data()
) {
  example_mapping()[[cond]]
} else {
  NULL
}
      )

    })

  )

})
output$custom_org_loader_ui <- renderUI({


  req(
  dataset2$data(),
  input$file,
  dataset1$data()
)

  mapping <- condition_organism_mapping()

  conditions_with_other <- names(mapping)[
    sapply(
      mapping,
      function(x) "Other" %in% x
    )
  ]

  if (length(conditions_with_other) == 0) {
    return(NULL)
  }

  tagList(

tags$div(

  class = "form-group",

  tags$label(

    class = "control-label",

    "Custom organism databases ",

    tags$span(
      `data-toggle` = "tooltip",
      `data-placement` = "right",
      `data-html` = "true",

      title = HTML(
        "<b>Required columns</b><br><br>

        • Entry<br>
        • Reviewed<br>
        • Protein names<br>
        • Gene Names<br>
        • Organism<br>
        • Length<br>
        • Gene Ontology (GO)<br>
        • Gene Ontology (biological process)<br>
        • Gene Ontology (cellular component)<br>
        • Gene Ontology (molecular function)<br>
        • Gene Ontology IDs<br><br>

        These databases can be downloaded directly from <b>UniProt</b> for the organism of interest."
      ),

      tags$i(
        class = "fa fa-question-circle",
        style = "cursor:pointer;"
      )
    )
  )

),

    lapply(
      conditions_with_other,
      function(cond) {

tagList(

  fileInput(
    inputId = paste0(
      "custom_org_",
      make.names(cond)
    ),
label = tags$span(
  style = "
    font-weight: normal;
    color: #333;
    font-size: 14px;
  ",
  paste(
    cond,
    "- Upload Org_data"
  )
),
    accept = c(
      ".csv",
      ".tsv",
      ".txt",
      ".xlsx",
      ".xls"
    )
  )

)

      }
    )
  )

})

output$custom_org_mapping_ui <- renderUI({

  req(
    dataset2$data(),
    input$file,
    dataset1$data()
  )

  mapping <- condition_organism_mapping()

  conditions_with_other <- names(mapping)[
    sapply(
      mapping,
      function(x)
        !is.null(x) &&
        "Other" %in% x
    )
  ]

  if (length(conditions_with_other) == 0) {
    return(NULL)
  }

  tagList(
tags$label(
  class = "control-label",

  "Condition to organism mapping ",

  tags$span(
    `data-toggle` = "tooltip",
    `data-placement` = "right",
    `data-html` = "true",

    title = HTML(
      "<b>Condition to organism mapping</b><br><br>

      Select the organism(s) used in each condition.<br><br>

      • Multiple organisms per condition are allowed.<br>

      • If the organism is not listed, select <b>Other</b>.<br>

      • Additional organism databases can be uploaded afterwards."
    ),

    tags$i(
      class = "fa fa-question-circle",
      style = "color:#337ab7; cursor:pointer;"
    )
  )
),

    lapply(
      conditions_with_other,
      function(cond) {

        selectInput(
          inputId = paste0(
            "custom_protein_org_",
            make.names(cond)
          ),

label = tags$span(
  style = "
    font-weight: normal;
    color: #333;
    font-size: 14px;
  ",
  paste(
    cond,
    "- Organism in protein dataset"
  )
),

          choices = sort(
            unique(
              as.character(
                dataset1$data()$Organism
              )
            )
          )
        )

      }
    )

  )

})

condition_organism_mapping <- reactive({
if (using_example_data()) {
  return(example_mapping())
}
  req(dataset2$data())

  design <- dataset2$data()

  conditions <- unique(
    as.character(design$Condition)
  )

  out <- lapply(
    conditions,
    function(cond) {

      input[[paste0(
        "organism_map_",
        make.names(cond)
      )]]

    }
  )

  names(out) <- conditions

  out

})
custom_protein_mapping <- reactive({

  req(dataset2$data())

  mapping <- condition_organism_mapping()

  conditions_with_other <- names(mapping)[
    sapply(
      mapping,
      function(x)
        !is.null(x) &&
        "Other" %in% x
    )
  ]

  out <- list()

  for (cond in conditions_with_other) {

    out[[cond]] <- input[[paste0(
      "custom_protein_org_",
      make.names(cond)
    )]]

  }

  out

})

custom_orgdata_files <- reactive({

  req(dataset2$data())

  mapping <- condition_organism_mapping()


conditions_with_other <- names(mapping)[
  sapply(
    mapping,
    function(x)
      !is.null(x) &&
      "Other" %in% x
  )
]


  out <- list()

  for (cond in conditions_with_other) {

    file_obj <- input[[paste0(
      "custom_org_",
      make.names(cond)
    )]]


id <- paste0(
  "custom_org_",
  make.names(cond)
)


    if (is.null(file_obj))
      next

    ext <- tolower(
  tools::file_ext(file_obj$name)
)

df <- switch(

  ext,

  csv = read.csv(
    file_obj$datapath,
    check.names = FALSE
  ),

  tsv = read.delim(
    file_obj$datapath,
    check.names = FALSE
  ),

  txt = read.delim(
    file_obj$datapath,
    check.names = FALSE
  ),

  xlsx = as.data.frame(
    readxl::read_excel(
      file_obj$datapath
    )
  ),

  xls = as.data.frame(
    readxl::read_excel(
      file_obj$datapath
    )
  ),

  NULL

)


out[[cond]] <- df



  }

  out

})

custom_orgdata_for_conditions <- reactive({

  x <- custom_orgdata_files()

  if (length(x) == 0) {
    return(NULL)
  }

  x

})





require_protein_organism_column <- reactive({

  mapping <- selected_organisms()

  if (length(mapping) == 0) {
    return(FALSE)
  }

  all_orgs <- unique(
    unlist(mapping)
  )

  length(all_orgs) > 1

})
mapping_complete <- reactive({
  if (using_example_data()) {
  return(TRUE)
}

  mapping <- condition_organism_mapping()

  all(
    sapply(
      mapping,
      function(x)
        !is.null(x) && length(x) > 0
    )
  )

})
output$protein_upload_message <- renderUI({

  req(dataset2$data())

  if (!mapping_complete()) {

    return(

      tags$div(

        style = "padding:10px; color:#666;",

        strong(
          "Please assign organisms to all conditions before uploading proteins."
        )

      )

    )

  }

if (require_protein_organism_column()) {

  return(

    tagList(

      tags$div(

        style = "
          margin-bottom:10px;
          padding:10px;
          background-color:#f8f9fa;
          border-left:4px solid #337ab7;
        ",

        strong("Information: "),
        "Multiple organisms detected. Protein dataset must contain an Organism column."

      ),

      tags$div(

        style = "
          margin-bottom:10px;
          padding:10px;
          background-color:#f8f9fa;
          border-left:4px solid #5cb85c;
        ",

        strong("GO analysis recommendation: "),
        "UniProt Accession IDs are recommended in the Accession column if you plan to use GO ORA or GO GSEA. Other identifier formats may not support GO annotation matching."

      )

    )

  )

}


tagList(

  tags$div(

    style = "
      margin-bottom:10px;
      padding:10px;
      background-color:#f8f9fa;
      border-left:4px solid #337ab7;
    ",

    strong("Information: "),
    "Only one organism detected. Organism column in the protein dataset is optional."

  ),

  tags$div(

    style = "
      margin-bottom:10px;
      padding:10px;
      background-color:#f8f9fa;
      border-left:4px solid #5cb85c;
    ",

    strong("GO analysis recommendation: "),
    "UniProt Accession IDs are recommended in the Accession column if you plan to use GO ORA or GO GSEA. Other identifier formats may not support GO annotation matching."

  )

)

})
output$protein_upload_ui <- renderUI({
if (using_example_data()) {

  return(

    tags$div(

      style = "
        margin-bottom:10px;
        padding:10px;
        background-color:#dff0d8;
        border-left:4px solid #5cb85c;
      ",

      strong("Example dataset loaded. "),
      "The protein matrix has been loaded automatically."

    )

  )

}
  if (!mapping_complete()) {

    return(NULL)

  }

fileInput(
  "file",

  tagList(
    "Proteins ",

    tags$span(
      `data-toggle` = "tooltip",
      `data-placement` = "right",
      `data-html` = "true",

        title = HTML(
          "<b>Required columns</b><br><br>

          • <b>Accession</b> (always required).<br>

          • <b>Organism</b> (required when more than one organism is present in the experiment).<br><br>

          All remaining columns should contain quantitative protein abundance measurements."
        ),

      tags$i(
        class = "fa fa-question-circle",
        style = "color:#337ab7; cursor:pointer;"
      )
    )
  ),
    accept = c(
      ".csv",
      ".tsv",
      ".txt",
      ".xlsx",
      ".xls"
    )
  )

})

output$split_mode_ui <- renderUI({

  req(condition_organism_mapping())

  all_orgs <- unique(
    unlist(
      condition_organism_mapping()
    )
  )

  if (length(all_orgs) <= 1) {

    return(

      selectInput(
        "split_mode",
        "Split dataset by:",
        choices = c(
          "Condition only" = "condition"
        ),
        selected = "condition"
      )

    )

  }

  selectInput(
    "split_mode",
    "Split dataset by:",
    choices = c(
      "Condition only" = "condition",
      "Condition + Organism" = "condition_organism"
    ),
    selected = "condition"
  )

})
#select how to split the datasets
observe({

  req(processed$data())

  updateSelectInput(
    session,
    "selected_split",
    choices = names(processed$data()),
    selected = names(processed$data())[1]
  )

})

  # Dataset selectors


output$dataset_selector2 <- renderUI({

  choices <- c(
    "Proteins" = "data1"
  )

  if (
    !is.null(dataset2$data())
  ) {

    choices <- c(
      choices,
      "Experimental Design" = "data2"
    )

  }

  radioButtons(
    inputId = "dataset2",
    label = "Select dataset to preview:",
    choices = choices,
    selected = "data1",
    inline = TRUE
  )

})


  # File info
  output$value <- renderPrint({
    if (input$dataset == "data1") {
      str(input$file)
    } else {
      str(input$file2)
    }
  })

# Dataset preview
output$head <- DT::renderDT({
  req(input$dataset2)

  # determine df
df <- if (input$dataset2 == "data1") {

  dataset1$data()

} else {

  req(!is.null(dataset2$data()))

  dataset2$data()

}

  req(!is.null(df))
  validate(need(NCOL(df) > 0, "No columns in data."))

numeric_cols <- which(
  vapply(
    df,
    is.numeric,
    logical(1)
  )
)

if (length(numeric_cols) > 0) {

  tbl <- DT::datatable(
    df,
    options = list(
      pageLength = 5,
      scrollX = TRUE,
      autoWidth = FALSE
    ),
    width = "100%"
  )

  tbl <- DT::formatRound(
    tbl,
    columns = numeric_cols,
    digits = 2
  )

  return(tbl)

}

DT::datatable(
  df,
  options = list(
    pageLength = 5,
    scrollX = TRUE,
    autoWidth = FALSE
  ),
  width = "100%"
)
})
outputOptions(output, "head", suspendWhenHidden = FALSE) #this line fixes the glitch of tables suspending when collapsed and repoen
#might use this line later in other tables as well.


output$str <- renderPrint({
  req(dataset1$data())
  df <- if (is.null(input$dataset3) || input$dataset3 == "data1") {
    dataset1$data()
  } else {
    req(input$file2, dataset2$data())
    dataset2$data()
  }
  str(df)
})

#split dataset table


output$split_table <- renderDT({

  req(
    processed$data(),
    input$selected_split
  )

  df <- processed$data()[[input$selected_split]]

  numeric_cols <- which(
    vapply(
      df,
      is.numeric,
      logical(1)
    )
  )

  DT::datatable(
    df,
    options = list(
      pageLength = 5,
      scrollX = TRUE
    )
  ) %>%
    DT::formatRound(
      columns = numeric_cols,
      digits = 1
    )

})





############################################################
# === Optional: Show dataset tabs toggle ===
output$show_dataset_tabs <- reactive({
  !is.null(input$file) && !is.null(input$file2)
})
outputOptions(output, "show_dataset_tabs", suspendWhenHidden = FALSE)
already_log2 <- reactive({
  isTRUE(input$already_log2)
})

already_normalized <- reactive({
isTRUE(input$already_normalized)
})
#cache
analysis_cache <- reactiveValues()

#QC tab - missing values and zero values
output$qc_tabs_ui <- renderUI({
  req(processed$data())

  datasets <- processed$data()

  tabs <- lapply(names(datasets), function(name) {
    tabPanel(
      title = name,
      value = name,
      qcPanelModuleUI(name, name)
    )
  })

  do.call(tabsetPanel, c(
    list(id = "qc_dataset_tabs"),
    tabs
  ))
})


#looped server
protplot_modules <- reactiveValues()
cluster_version <- reactiveVal(0)
local({

  datasets <- reactive({
    req(processed$data())
    processed$data()
  })

  # ✅ version only for clustering
  


normalization_results <- reactiveValues()
clustering_results <- reactiveValues()
go_ora_inputs <- reactiveValues()
go_ora_results <- reactiveValues()
model_results <- reactiveValues()

tpp_results_rv <- reactiveVal(NULL)
tpp_summary_rv <- reactiveVal(NULL)
tpp_fitdata_rv <- reactiveVal(NULL)





observeEvent(input$run_tpp_import, {
vehicle_name <- input$tpp_vehicle
treatment_name <- input$tpp_treatment

qc_vehicle_done <- tryCatch(
  normalization_results[[vehicle_name]]$qc_applied(),
  error = function(e) FALSE
)

qc_treatment_done <- tryCatch(
  normalization_results[[treatment_name]]$qc_applied(),
  error = function(e) FALSE
)

if (!isTRUE(qc_vehicle_done) || !isTRUE(qc_treatment_done)) {

  showNotification(
    "Quality Control must be applied to both selected conditions before running TPP analysis.",
    type = "error",
    duration = 8
  )

  return()

}
  withProgress(

    message = "Running TPP analysis",

detail = "Preparing TPP analysis...",

    value = 0,

    {


req(
  normalization_results[[vehicle_name]],
  normalization_results[[treatment_name]]
)


vehicle_df <-
  normalization_results[[vehicle_name]]$plot_data()

treatment_df <-
  normalization_results[[treatment_name]]$plot_data()

vehicle_lowest_temp <- tryCatch(
  normalization_results[[vehicle_name]]$lowest_temp_normalized(),
  error = function(e) FALSE
)

treatment_lowest_temp <- tryCatch(
  normalization_results[[treatment_name]]$lowest_temp_normalized(),
  error = function(e) FALSE
)

if(
  !isTRUE(input$Already_lowestnorm) &&
  (
    !isTRUE(vehicle_lowest_temp) ||
    !isTRUE(treatment_lowest_temp)
  )
){

  showNotification(
    paste(
      "TPP analysis requires either:",
      "\n• 'Normalize to lowest temperature' to be applied in both datasets",
      "\nor",
      "\n• 'Dataset already normalized to lowest temperature' to be checked."
    ),
    type = "error",
    duration = 8
  )

  return()

}

n_proteins <- length(
  unique(
    vehicle_df[[1]]
  )
)
incProgress(
  0.02,
  detail = paste(
    "Detected",
    format(n_proteins, big.mark = ","),
    "proteins."
  )
)
convert_to_tpp_long <- function(df){

df %>%
tidyr::pivot_longer(
cols = -1,
names_to = "Replicate",
values_to = "RelAbundance"
) %>%
tidyr::separate(
Replicate,
into = c(
"Temperature",
"Replicate"
),
sep = "_"
) %>%
dplyr::mutate(
Temperature = as.numeric(
Temperature
)
)
}
incProgress(
  0.10,
  detail = "Normalizing vehicle condition"
)

vehicle_rel <- convert_to_tpp_long(vehicle_df)

incProgress(
  0.10,
  detail = "Normalizing treatment condition"
)

treatment_rel <- convert_to_tpp_long(treatment_df)


build_tpp_matrix <- function(
    x,
    rep_name,
    all_temps
){

  out <- x %>%

    dplyr::filter(
      Replicate == rep_name
    ) %>%

    dplyr::select(
      Accession,
      Temperature,
      RelAbundance
    ) %>%

    tidyr::pivot_wider(
      names_from = Temperature,
      values_from = RelAbundance
    )

  # -----------------------------------
  # FORCE ALL TEMPERATURES TO EXIST
  # -----------------------------------

  missing_cols <- setdiff(
    as.character(all_temps),
    names(out)
  )

  if(length(missing_cols) > 0){

    for(col in missing_cols){

      out[[col]] <- NA_real_

    }

  }

  temp_cols <- sort(
    all_temps,
    decreasing = TRUE
  )

  out <- out %>%
    dplyr::select(
      Accession,
      all_of(as.character(temp_cols))
    )

  names(out)[-1] <- paste0(
    "rel_fc_L",
    seq_along(temp_cols)
  )

  as.data.frame(out)

}

#setting the replicates 
replicates <- sort(
  unique(
    vehicle_rel$Replicate
  )
)

all_temps <- sort(
  unique(
    c(
      vehicle_rel$Temperature,
      treatment_rel$Temperature
    )
  )
)

vehicle_tpp_list <- lapply(
  replicates,
  function(rep){

build_tpp_matrix(
  vehicle_rel,
  rep,
  all_temps
)

  }
)

names(vehicle_tpp_list) <- paste0(
  "Vehicle_",
  replicates
)

treatment_tpp_list <- lapply(
  replicates,
  function(rep){

build_tpp_matrix(
  treatment_rel,
  rep,
  all_temps
)

  }
)

names(treatment_tpp_list) <- paste0(
  "Treatment_",
  replicates
)

#set the comparisons
all_experiments <- c(

  names(vehicle_tpp_list),

  names(treatment_tpp_list)

)

temps <- sort(
  unique(vehicle_rel$Temperature),
  decreasing = TRUE
)

config_test <- data.frame(

  Experiment = all_experiments,

  Condition = c(

    rep(
      "Vehicle",
      length(vehicle_tpp_list)
    ),

    rep(
      "Treatment",
      length(treatment_tpp_list)
    )

  ),

  stringsAsFactors = FALSE

)

# comparison columns
for(i in seq_along(replicates)){

  col_name <- paste0(
    "ComparisonVT",
    i
  )

  config_test[[col_name]] <- NA

  config_test[
    i,
    col_name
  ] <- "x"

  config_test[
    i + length(replicates),
    col_name
  ] <- "x"

}
#temperature columns
for(i in seq_along(temps)){

  config_test[
    paste0("L", i)
  ] <- temps[i]

}

tpp_data_list <- c(
  vehicle_tpp_list,
  treatment_tpp_list
)

incProgress(
  0.15,
  detail = "Preparing TPP input tables"
)
trData <- tpptrImport(
  configTable = config_test,
  data = tpp_data_list,
  idVar = "Accession",
  fcStr = "rel_fc_",
  qualColName = NA
)
incProgress(
  0.25,
  detail = paste(
    "Fitting melting curves for",
    format(n_proteins, big.mark = ","),
    "proteins. This step may take several minutes."
  )
)
fitData <- tpptrCurveFit(
  data = trData,
  nCores = 1
)

cat("\nFITDATA STRUCTURE\n")

str(
  fitData,
  max.level = 2
)
tpp_fitdata_rv(
fitData
)
pValFilters <- list(

  minR2 = input$tpp_min_r2,

  maxPlateau = input$tpp_max_plateau

)
incProgress(
  0.20,
  detail = "Calculating ΔTm statistics"
)

TRresults <- tpptrAnalyzeMeltingCurves(
  data = fitData,
  pValFilter = pValFilters
)
incProgress(
  0.03,
  detail = "Generating summary table"
)
summary_tbl <- as.data.frame(TRresults)
pass_cols <- grep(
  "^passed_filter_",
  names(summary_tbl),
  value = TRUE
)

proteins_passing_all <- 0

if(length(pass_cols) > 0){

  proteins_passing_all <- sum(
    apply(
      summary_tbl[, pass_cols, drop = FALSE],
      1,
      function(x)
        all(x, na.rm = TRUE)
    ),
    na.rm = TRUE
  )

}
diff_cols <- grep(
  "^diff_meltP_",
  names(summary_tbl),
  value = TRUE
)

pval_cols <- grep(
  "^pVal_adj_",
  names(summary_tbl),
  value = TRUE
)
proteins_passing_all_and_tm <- 0

if(
  length(pass_cols) > 0 &&
  length(pval_cols) > 0
){

  proteins_passing_all_and_tm <- sum(

    apply(
      summary_tbl[, pass_cols, drop = FALSE],
      1,
      function(x)
        all(x, na.rm = TRUE)
    )

    &

    apply(
      summary_tbl[, pval_cols, drop = FALSE],
      1,
      function(x)
        any(is.finite(x))
    ),

    na.rm = TRUE
  )

}



proteins_with_tm <- sum(
  apply(
    summary_tbl[, pval_cols, drop = FALSE],
    1,
    function(x)
      any(is.finite(x))
  ),
  na.rm = TRUE
)

summary_stats <- list(

  proteins_tested =
    nrow(summary_tbl),

proteins_passing_all =
  proteins_passing_all,

proteins_with_tm =
  proteins_with_tm,

proteins_passing_all_and_tm =
  proteins_passing_all_and_tm,

  significant_005 =
    sum(
      apply(
        summary_tbl[, pval_cols, drop = FALSE],
        1,
        min,
        na.rm = TRUE
      ) < 0.05,
      na.rm = TRUE
    ),

  significant_001 =
    sum(
      apply(
        summary_tbl[, pval_cols, drop = FALSE],
        1,
        min,
        na.rm = TRUE
      ) < 0.01,
      na.rm = TRUE
    ),

  median_dTm =
    median(
      apply(
        summary_tbl[, diff_cols, drop = FALSE],
        1,
        mean,
        na.rm = TRUE
      ),
      na.rm = TRUE
    ),

  percentile95_dTm =
    quantile(
      abs(
        apply(
          summary_tbl[, diff_cols, drop = FALSE],
          1,
          mean,
          na.rm = TRUE
        )
      ),
      0.95,
      na.rm = TRUE
    )

)
incProgress(
  0.02,
  detail = "Finalizing results"
)
tpp_summary_rv(summary_stats)
tpp_results_rv(TRresults)

  }

)
})

tpp_display_table <- reactive({

  req(
    tpp_results_rv()
  )
  

  tbl <- as.data.frame(
    tpp_results_rv()
  )
  pval_cols <- grep(
  "^pVal_adj_",
  names(tbl),
  value = TRUE
)

if(length(pval_cols) > 0){

  tbl$tm_calculated <-
    apply(
      tbl[, pval_cols, drop = FALSE],
      1,
      function(x)
        any(is.finite(x))
    )

}

  pass_cols <- grep(
    "^passed_filter_",
    names(tbl),
    value = TRUE
  )

  if(length(pass_cols) > 0){

    tbl$passes_all_replicates <-
      apply(
        tbl[, pass_cols, drop = FALSE],
        1,
        function(x)
          all(
            x,
            na.rm = TRUE
          )
      )

  }

if(
  isTRUE(
    input$tpp_show_tm_only
  )
){

  tbl <- tbl %>%
    dplyr::filter(
      tm_calculated
    )

}

if(
  isTRUE(
    input$tpp_show_hits_only
  )
){

  tbl <- tbl %>%
    dplyr::filter(
      passes_all_replicates
    )

}

  tbl

})
#result table
output$tpp_results_table <- DT::renderDataTable({

  req(
    tpp_results_rv()
  )

tbl <- tpp_display_table()

diff_cols <- grep(
  "^diff_meltP_",
  names(tbl),
  value = TRUE
)

pval_cols <- grep(
  "^pVal_adj_",
  names(tbl),
  value = TRUE
)

pass_cols <- grep(
  "^passed_filter_",
  names(tbl),
  value = TRUE
)
tbl$Mean_DeltaTm <- apply(
  tbl[, diff_cols, drop = FALSE],
  1,
  mean,
  na.rm = TRUE
)

tbl$SD_DeltaTm <- apply(
  tbl[, diff_cols, drop = FALSE],
  1,
  sd,
  na.rm = TRUE
)

tbl$Min_Adjusted_P <- apply(
  tbl[, pval_cols, drop = FALSE],
  1,
  min,
  na.rm = TRUE
)

tbl$Replicates_Passing <- rowSums(
  tbl[, pass_cols, drop = FALSE],
  na.rm = TRUE
)

pass_cols <- grep(
  "^passed_filter_",
  names(tbl),
  value = TRUE
)
if(length(pass_cols) > 0){

  tbl$passes_all_replicates <-
    apply(
      tbl[, pass_cols, drop = FALSE],
      1,
      function(x)
        all(
          x,
          na.rm = TRUE
        )
    )

}
if(
  isTRUE(
    input$tpp_show_hits_only
  )
){

  tbl <- tbl %>%
    dplyr::filter(
      passes_all_replicates
    )

}


  keep_cols <- c(

    "Protein_ID",

    grep(
      "^diff_meltP_",
      names(tbl),
      value = TRUE
    ),

    grep(
      "^pVal_adj_",
      names(tbl),
      value = TRUE
    ),

    grep(
      "^passed_filter_",
      names(tbl),
      value = TRUE
    )

  )

  keep_cols <- keep_cols[
    keep_cols %in% names(tbl)
  ]

summary_cols <- c(

  "Protein_ID",

  "Mean_DeltaTm",
  "SD_DeltaTm",

  "Min_Adjusted_P",

  "Replicates_Passing",

  diff_cols,
  pval_cols,
  pass_cols

)

summary_cols <- summary_cols[
  summary_cols %in% names(tbl)
]

DT::datatable(
  tbl[, summary_cols, drop = FALSE],
    options = list(
      pageLength = 5,
      scrollX = TRUE
    ),

    rownames = FALSE,
    selection = "single"

  )

})




selected_tpp_protein <- reactive({

  req(
    input$tpp_results_table_rows_selected
  )

  tbl <- tpp_display_table()

  idx <- input$tpp_results_table_rows_selected

  tbl$Protein_ID[idx]

})

make_tpp_curve_plot <- function() {

  req(
    selected_tpp_protein()
  )

  req(
    tpp_fitdata_rv()
  )

  protein_id <- selected_tpp_protein()

fd <- Biobase::pData(
  Biobase::featureData(
    tpp_fitdata_rv()[["Vehicle_REP1"]]
  )
)

cat(
  "\n======================\n"
)

cat(
  "SELECTED PROTEIN:",
  protein_id,
  "\n"
)

print(
  fd[protein_id, ]
)

plot.new()

title(
  main = protein_id
)
fit_list <- tpp_fitdata_rv()

cat("\n====================\n")
cat("FIT LIST CLASS\n")
cat("====================\n")

print(class(fit_list))
print(names(fit_list))
cat("\n====================\n")
cat("FIT LIST CLASS\n")
cat("====================\n")

print(class(fit_list))

print(names(fit_list))

cat("\nFIRST ELEMENT STRUCTURE\n")

str(
  fit_list[[1]],
  max.level = 1
)


obs_data <- list()
curve_data <- list()

for(exp_name in names(fit_list)) {

  eset <- fit_list[[exp_name]]


cat("\nPDATA COLNAMES:\n")
print(colnames(Biobase::pData(eset)))

cat("\nPDATA HEAD:\n")
print(head(Biobase::pData(eset)))

cat("\n=====================\n")
cat("EXPERIMENT:", exp_name, "\n")
cat("=====================\n")

print(
  colnames(
    Biobase::exprs(eset)
  )
)

cat("\nASSAYDATA ELEMENTS:\n")
print(Biobase::assayDataElementNames(eset))

cat("\nEXPRS DIM:\n")
print(dim(Biobase::exprs(eset)))

cat("\nALL SLOT NAMES:\n")
print(slotNames(eset))

  if(!(protein_id %in%
       Biobase::featureNames(eset))) {
    next
  }

  expr_vals <- as.numeric(
    Biobase::exprs(eset)[protein_id, ]
  )
condition <- sub(
  "_REP\\d+$",
  "",
  exp_name
)

replicate <- sub(
  "^.*_(REP\\d+)$",
  "\\1",
  exp_name
)
 temps <- Biobase::pData(eset)$temperature

obs_data[[exp_name]] <- data.frame(
  Temperature = temps,
  Abundance = expr_vals,
  Condition = condition,
  Replicate = replicate
)

  fit_params <- Biobase::pData(
    Biobase::featureData(eset)
  )

  fit_params <- fit_params[
    protein_id,
    ,
    drop = FALSE
  ]
cat("\n=====================\n")
cat("EXPERIMENT:", exp_name, "\n")
cat("=====================\n")

print(names(Biobase::fData(eset)))

print(colnames(Biobase::fData(eset)))

print(
  Biobase::fData(eset)[protein_id, , drop = FALSE]
)
if(
  is.na(fit_params$a) ||
  is.na(fit_params$b)
) {
  next
}

mTmp <- TPP:::fitSigmoidTR(
  xVec = temps,
  yVec = expr_vals,
  startPars = c(
    Pl = 0,
    a = 550,
    b = 10
  ),
  maxAttempts = 500,
  fixT0 = TRUE
)

if(inherits(mTmp, "try-error")) {
  next
}

x_grid <- seq(
  min(temps, na.rm = TRUE),
  max(temps, na.rm = TRUE),
  length.out = 100
)

y_pred <- TPP:::robustNlsPredict(
  model = mTmp,
  newdata = list(
    x = x_grid
  )
)

curve_data[[exp_name]] <- data.frame(
  Temperature = x_grid,
  Abundance = y_pred,
  Condition = condition,
  Replicate = replicate
)

}



obs_data <- dplyr::bind_rows(
  obs_data
)
curve_data <- dplyr::bind_rows(
  curve_data
)
p<- ggplot() +

geom_point(
  data = obs_data,
  aes(
    Temperature,
    Abundance,
    color = Condition,
    shape = Replicate
  ),
  size = 3
) +

geom_line(
  data = curve_data,
  aes(
    Temperature,
    Abundance,
    color = Condition,
    linetype = Replicate,
    group = interaction(
      Condition,
      Replicate
    )
  ),
  linewidth = 1
) +

  theme_classic() +

  labs(
    title = protein_id,
    x = "Temperature (°C)",
    y = "Relative abundance"
  )
return(p)
}
output$tpp_curve_plot <- renderPlot({

  make_tpp_curve_plot()

})
#tpp padjust plot

output$tpp_padj_hist <- renderPlot({

  req(tpp_display_table())

  tbl <- tpp_display_table()

  pval_cols <- grep(
    "^pVal_adj_",
    names(tbl),
    value = TRUE
  )

  vals <- unlist(
    tbl[, pval_cols, drop = FALSE]
  )

  vals <- vals[
    is.finite(vals)
  ]

  ggplot(
    data.frame(
      pAdj = vals
    ),
    aes(x = pAdj)
  ) +

    geom_histogram(
      bins = 40,
      fill = "steelblue",
      color = "black"
    ) +

    theme_minimal() +

    labs(
      title = "Adjusted P-value Distribution",
      x = "Adjusted P-value",
      y = "Count"
    )

})

#TPP delta tm
output$tpp_dtm_hist <- renderPlot({

  req(tpp_display_table())

  tbl <- tpp_display_table()

  diff_cols <- grep(
    "^diff_meltP_",
    names(tbl),
    value = TRUE
  )

  vals <- unlist(
    tbl[, diff_cols, drop = FALSE]
  )

  vals <- vals[
    is.finite(vals)
  ]

  ggplot(
    data.frame(
      DeltaTm = vals
    ),
    aes(x = DeltaTm)
  ) +

    geom_histogram(
      bins = 40,
      fill = "tomato",
      color = "black"
    ) +

    theme_minimal() +

    labs(
      title = expression(Delta * "Tm Distribution"),
      x = expression(Delta * "Tm"),
      y = "Count"
    )

})

#tpp summary text
output$tpp_summary <- renderText({

  req(tpp_summary_rv())

  s <- tpp_summary_rv()

paste0(

  "Proteins tested: ",
  s$proteins_tested,

  "\nProteins passing all replicate filters: ",
  s$proteins_passing_all,

"\nProteins with ΔTm calculated: ",
s$proteins_with_tm,

"\nProteins passing all filters AND with ΔTm calculated: ",
s$proteins_passing_all_and_tm,

    "\nSignificant (adj p < 0.05): ",
    s$significant_005,

    "\nSignificant (adj p < 0.01): ",
    s$significant_001,

    "\nMedian ΔTm: ",
    round(s$median_dTm, 2),

    "\n95th percentile |ΔTm|: ",
    round(s$percentile95_dTm, 2)

  )

})


observeEvent(datasets(), {

  withProgress(

    message = "Initializing analysis modules",
    value = 0,

  {
  t0 <- Sys.time()

  cat(
    "\n=============================\n",
    "DATASET INITIALIZATION START\n",
    "=============================\n"
  )

    #  bump version ONLY for clustering modules
    cluster_version(cluster_version() + 1)

    # reset normalization results
# clear all module stores

      for (n in names(reactiveValuesToList(normalization_results))) {
        normalization_results[[n]] <- NULL
      }

      for (n in names(reactiveValuesToList(clustering_results))) {
        clustering_results[[n]] <- NULL
      }


      for (n in names(reactiveValuesToList(model_results))) {
        model_results[[n]] <- NULL
      }

dataset_names <- names(datasets())

total_datasets <- length(dataset_names)
    for (name in names(datasets())) {

      local({
        dataset_name <- name
incProgress(
  1 / total_datasets,
  detail = paste(
    "Loading",
    dataset_name
  )
)

    safe_id <- gsub(
  "[^[:alnum:]_.+-]",
  "_",
  dataset_name
)
dataset_t0 <- Sys.time()


dataset_orgs <- reactive({

  message("ORGS_A before dataset_organisms()")

  tmp <- processed$dataset_organisms()

message("ORGS_B after dataset_organisms()")

str(tmp)

message("TMP NAMES")

print(names(tmp))

  message("DATASET NAME:")
  print(dataset_name)

  out <- tmp[[dataset_name]]

  message("ORGS_VALUE")
  str(out)

  out
})






        # ---- missing data ----
        filtered_data <- missingDataModuleServer(
          id = dataset_name,
          dataReactive = reactive({
            datasets()[[dataset_name]]
          }),
          choicesInput = reactive(input[[paste0(dataset_name, "-choices")]]),
          methodInput  = reactive(input[[paste0(dataset_name, "-impute_method")]])
        )

        # ---- stats ----
        missing_stats <- missingStatsModule(
          reactive({
            datasets()[[dataset_name]]
          })
        )

zero_stats <- zeroStatsModule(
  id = dataset_name,
  filtered_data = filtered_data,
  already_log2 = already_log2
)

        # ---- QC ----
        qcPanelModuleServer(
          id = dataset_name,
          filtered_data = filtered_data,
          missing_stats = missing_stats,
          zero_stats    = zero_stats,
          already_log2 = already_log2
        )

        # ---- NORMALIZATION ----
        
        norm_res <- normalizationPanelModuleServer(
          id = dataset_name,
          zero_stats_data = zero_stats,
          dataset_org_data = reactive({
            dataset2$Org_data()
          }),
          all_datasets = datasets
        )

        normalization_results[[dataset_name]] <- norm_res



        # ------------------------
        #  ONLY CLUSTERING IS VERSIONED
        # ------------------------

        cluster_id <- paste0(dataset_name, "_", cluster_version())

        cluster_res <- clusterPanelModuleServer(
          id = cluster_id,

          plot_data = reactive({
            req(normalization_results[[dataset_name]])
            normalization_results[[dataset_name]]$plot_data()
          }),

          dataset_names = reactive({
            setdiff(names(datasets()), dataset_name)
          }),

          normalization_results = normalization_results
        )

        # STORE RESULTS USING dataset_name (NOT cluster_id)
        clustering_results[[dataset_name]] <- cluster_res
        # ------------------------
        # MODEL FITTING
        # ------------------------

        model_results[[dataset_name]] <-
          modelPanelModuleServer(
            id = safe_id,

            normalized_data = reactive({
              req(normalization_results[[dataset_name]])
              normalization_results[[dataset_name]]$plot_data()
            }),

            dataclustered = reactive({
              req(clustering_results[[dataset_name]])
              clustering_results[[dataset_name]]$clustering$data()
            }),

            reclustered_data = reactive({
              req(clustering_results[[dataset_name]])
              clustering_results[[dataset_name]]$reclustering$data()
            })
          )

# ------------------------
# MELTING CURVE
# ------------------------

meltingCurveModuleServer(
  id = safe_id,

  normalized_data = reactive({
    req(normalization_results[[dataset_name]])
    normalization_results[[dataset_name]]$plot_data()
  }),

all_normalized_data = reactive({

  out <- list()

  for (nm in names(reactiveValuesToList(normalization_results))) {

    if (is.null(normalization_results[[nm]]))
      next

qc_done <- tryCatch(
  normalization_results[[nm]]$qc_applied(),
  error = function(e) FALSE
)

if (!isTRUE(qc_done))
  next


    out[[nm]] <- tryCatch(
      normalization_results[[nm]]$plot_data(),
      error = function(e) NULL
    )

  }

  out
}),
normalization_results = normalization_results,
already_lowestnorm = reactive(input$Already_lowestnorm)
)
# ------------------------
# ✅ GO ORA INPUTS
# ------------------------

go_ora_inputs[[dataset_name]] <-
  goORAClusterPanelServer(
    id = safe_id,
Org_data = reactive({

 

  orgs <- dataset_orgs()




  print(class(orgs))

      known_orgs <- c(
        "Sulfolobus acidocaldarius",
        "Parageobacillus thermoglucosidasius",
        "Caldimonas thermodepolymerans",
        "Haloferax volcanii",
        "Thermosynechococcus elongatus",
        "Thermus thermophilus",
        "E.coli",
        "Escherichia coli",
        "E coli",
        "H.sapiens",
        "Homo sapiens",
        "human"
      )

      if (!all(orgs %in% known_orgs)) {



        return(NULL)
      }

cat("\nORGS REQUESTED:\n")
print(orgs)

cat("\nAVAILABLE ORGANISMS:\n")
print(unique(dataset2$Org_data()$Organism))

      x <- dataset2$Org_data() %>%
  dplyr::filter(
    Organism %in% orgs
  )

cat(
  "\nFILTERED ORGANISMS:\n"
)

print(unique(x$Organism))


      x

    })
  )





# ------------------------
# GO ORA RESULTS
# ------------------------



go_ora_results[[dataset_name]] <-
  goORAPlotsPanelServer(

    id = safe_id,

    dataclustered = reactive({
      req(clustering_results[[dataset_name]])
      clustering_results[[dataset_name]]$clustering$data()
    }),

 Org_data = reactive({



  orgs <- dataset_orgs()



  known_orgs <- c(
    "Sulfolobus acidocaldarius",
    "Parageobacillus thermoglucosidasius",
    "Caldimonas thermodepolymerans",
    "Haloferax volcanii",
    "Thermosynechococcus elongatus",
    "Thermus thermophilus",
    "E.coli",
    "Escherichia coli",
    "E coli",
    "H.sapiens",
    "Homo sapiens",
    "human"
  )

  if (!all(orgs %in% known_orgs)) {



    return(NULL)
  }



  x <- dataset2$Org_data() %>%
  dplyr::filter(
    Organism %in% orgs
  )

cat(
  "\nFILTERED ORGANISMS:\n"
)

print(unique(x$Organism))


  x

}),

    num_clusters = reactive({
      req(clustering_results[[dataset_name]])
      clustering_results[[dataset_name]]$num_clusters()
    }),

    cluster_inputs = go_ora_inputs[[dataset_name]]

  )

# ------------------------
# ✅ GO GSEA
# ------------------------

goGSEAPanelModuleServer(

  id = safe_id,


  hypothesis_result = reactive({

    req(model_results[[dataset_name]])

    model_results[[dataset_name]]$hypothesis_result()

  }),

  term2gene_data = reactive({

    req(go_ora_results[[dataset_name]])

    go_ora_results[[dataset_name]]$term2gene_data()

  }),

  orgdata = reactive({

    req(go_ora_results[[dataset_name]])

    go_ora_results[[dataset_name]]$orgdata()

  })

)

#protein plots module

# ------------------------
#  PROTEIN PLOTS MODULE (NOW IN CORRECT SCOPE)
# ------------------------

if (is.null(protplot_modules[[dataset_name]])) {

  protplot_modules[[dataset_name]] <- protplotsPanelModuleServer(
    id = dataset_name,

    normalized_data = reactive({
      req(normalization_results[[dataset_name]])
      normalization_results[[dataset_name]]$plot_data()
    }),

    dataclustered = reactive({
      req(clustering_results[[dataset_name]])
      clustering_results[[dataset_name]]$clustering$data()
    })
  )

}

      })
    }
})
  }, ignoreInit = FALSE)
#protplot

observeEvent(
  {
    lapply(names(protplot_modules), function(n) {

      mod <- protplot_modules[[n]]

      if (
        !is.null(mod) &&
        is.list(mod) &&
        "plot_trigger" %in% names(mod) &&
        is.function(mod$plot_trigger)
      ) {
        mod$plot_trigger()
      }

    })
  },
 {

    withProgress(

      message = "Building protein plots",
      detail = "Preparing datasets...",
      value = 0,

      {

        cat("== Plot button clicked (multi-dataset) ==\n")

        combined_list <- list()



    # ✅ LOOP THROUGH DATASETS
total_datasets <- length(names(protplot_modules))
dataset_counter <- 0
    for (name in names(protplot_modules)) {
dataset_counter <- dataset_counter + 1

incProgress(
  1 / total_datasets,
  detail = paste(
    "Processing",
    name,
    "(",
    dataset_counter,
    "of",
    total_datasets,
    ")"
  )
)
      module <- protplot_modules[[name]]

      # ✅ selection
      selected <- isolate(tryCatch(module$selected_data(), error = function(e) NULL))

      cat("Selected for", name, ": ",
          if (is.null(selected)) "NULL" else paste("rows =", nrow(selected)), "\n")

      if (is.null(selected) || (is.data.frame(selected) && nrow(selected) == 0)) {
        next
      }

      # ✅ dataset (clustered -> normalized fallback)
df <- {

  cat(">> entering df construction for", name, "\n")

  if (is.null(normalization_results[[name]])) {
    cat("!! normalization missing for", name, "\n")
    NULL
  } else {

    d_norm <- tryCatch(
      normalization_results[[name]]$plot_data(),
      error = function(e) {
        cat("!! normalization failed for", name, "\n")
        NULL
      }
    )

    d_clust <- tryCatch(
      clustering_results[[name]]$clustering$data(),
      error = function(e) {
        cat("!! clustering failed for", name, "\n")
        NULL
      }
    )

    if (!is.null(d_clust) && is.data.frame(d_clust) && nrow(d_clust) > 0) {
      cat(".. using clustered data for", name, "\n")
      d_clust
    } else if (!is.null(d_norm) && is.data.frame(d_norm) && nrow(d_norm) > 0) {
      cat(".. using normalized data for", name, "\n")
      d_norm
    } else {
      cat("!! both datasets NULL for", name, "\n")
      NULL
    }
  }
}

cat("DF class for", name, ":", class(df), "\n")

if (!is.data.frame(df)) {
  cat("!! df is NOT a data.frame for", name, "\n")
  print(str(df))
}
      # ✅ Org data
      



org <- tryCatch({
  dataset2$Org_data()
}, error = function(e) {
  cat("!! Org_data failed for", name, "\n")
  NULL
})





      # PROCESS
df_long <- process_dataset(
  df,
  selected,
  org,
  dataset_label = name,
  line_type = name,
  debug_tag = name
)



      if (!is.null(df_long)) {
        combined_list[[name]] <- df_long
      }
    }

    # ✅ combine all datasets
    combined <- dplyr::bind_rows(combined_list)

    if (length(combined_list) == 0 ||
        is.null(combined) ||
        !is.data.frame(combined) ||
        nrow(combined) == 0) {


      plot_data("no_selection")

    } else {


      plot_data(combined)
    }

      }

    )

  },
  ignoreInit = TRUE
)


##############################################################################

# Model Comparison logic
fstats_empirical_rv   <- reactiveVal(NULL)

fstats_theoretical_rv <- reactiveVal(NULL)

metrics_rv <- reactiveVal(NULL)


#data for plotting
combined_data_rv <- reactiveVal(NULL)

id_col_rv <- reactiveVal(NULL)

#predictions
predictions_rv <- reactiveVal(NULL)




observeEvent(input$run_model_comparison, {
  shinyjs::html(
  id = "nparc_console",
  html = ""
)
cond1 <- input$comparison_condition_1
cond2 <- input$comparison_condition_2

qc1_done <- tryCatch(
  normalization_results[[cond1]]$qc_applied(),
  error = function(e) FALSE
)

qc2_done <- tryCatch(
  normalization_results[[cond2]]$qc_applied(),
  error = function(e) FALSE
)

if (!isTRUE(qc1_done) || !isTRUE(qc2_done)) {

  showNotification(
    "Quality Control must be applied to both selected conditions before running NPARC.",
    type = "error",
    duration = 8
  )

  return()

}

req(
  cond1,
  cond2,
  normalization_results[[cond1]],
  normalization_results[[cond2]]
)

lowest1 <- tryCatch(
  normalization_results[[cond1]]$lowest_temp_normalized(),
  error = function(e) FALSE
)

lowest2 <- tryCatch(
  normalization_results[[cond2]]$lowest_temp_normalized(),
  error = function(e) FALSE
)

if(
  !isTRUE(input$Already_lowestnorm) &&
  (
    !isTRUE(lowest1) ||
    !isTRUE(lowest2)
  )
){
  showNotification(
    paste(
      "NPARC requires either:",
      "\n• 'Normalize to lowest temperature' to be applied in both datasets",
      "\nor",
      "\n• 'Dataset already normalized to lowest temperature' to be checked."
    ),
    type = "error",
    duration = 8
  )

  return()
}


convert_to_nparc_long <- function(df){

df %>%
tidyr::pivot_longer(
cols = -1,
names_to = "Replicate",
values_to = "RelAbundance"
) %>%
tidyr::separate(
Replicate,
into = c(
"Temperature",
"Replicate"
),
sep = "_"
) %>%
dplyr::mutate(
Temperature = as.numeric(
Temperature
)
)

}

  # Apply normalization with error handling
compar1_rel <- tryCatch(
  convert_to_nparc_long(
    normalization_results[[cond1]]$plot_data()
  ),

    error = function(e) {
      showNotification(e$message, type = "error")
      return(NULL)
    } 
  )
tpp_test <- compar1_rel %>%
  tidyr::unite(
    Label,
    Temperature,
    Replicate,
    sep = "_"
  ) %>%
  tidyr::pivot_wider(
    names_from = Label,
    values_from = RelAbundance
  )

write.csv(
  tpp_test,
  "TPP_TEST_INPUT.csv",
  row.names = FALSE
)



compar2_rel <- tryCatch(
  convert_to_nparc_long(
    normalization_results[[cond2]]$plot_data()
  ),
  

    error = function(e) {
      showNotification(e$message, type = "error")
      return(NULL)
    }
  )


# --------------------------------------------------
# TEMPORARY EXPORT FOR TPP TESTING
# --------------------------------------------------


  # If either dataset failed, stop
  if (is.null(compar1_rel) || is.null(compar2_rel)) return()

  # Combine datasets
combined_data <- dplyr::bind_rows(
  compar1_rel %>% dplyr::mutate(Dataset = cond1),
  compar2_rel %>% dplyr::mutate(Dataset = cond2)
)
#filter out protein with NA in relabundance
combined_data <- combined_data %>%
  dplyr::filter(!is.na(RelAbundance))
#select only proteins with data in every replicate

id_col <- names(combined_data)[1]
combined_data <- combined_data %>%
  dplyr::group_by(Dataset, !!rlang::sym(id_col)) %>%  # Dynamically pick the first column
  dplyr::mutate(n = dplyr::n()) %>%
  dplyr::group_by(Dataset) %>%
  dplyr::mutate(max_n = max(n)) %>%
  dplyr::ungroup()
# Filter for full curves per protein:
combined_data <- combined_data %>% 
  filter(n == max_n) %>%
  dplyr::select(-n, -max_n)

combined_data <- combined_data %>%
  dplyr::mutate(Dataset = factor(Dataset))

id_col_rv(id_col)
combined_data_rv(combined_data)

#Fit NPARC models
set.seed(1234)
withProgress(message = "Fitting models", value = 0, {

  incProgress(0.2, detail = "Preparing data...")

  fit_results <- withCallingHandlers(

    {

      NPARCfit(
        x = combined_data$Temperature,
        y = combined_data$RelAbundance,
        id = combined_data[[1]],
        groupsNull = NULL,
        groupsAlt = data.frame(
          Dataset = combined_data$Dataset
        ),
        BPPARAM = BiocParallel::SerialParam(
          progressbar = FALSE
        ),
        returnModels = FALSE
      )

    },

    message = function(m) {

      shinyjs::html(
        id = "nparc_console",
        html = paste0(
          m$message,
          "<br>"
        ),
        add = TRUE
      )
      shinyjs::runjs("
  var el = document.getElementById('nparc_console');
  if(el){
    el.scrollTop = el.scrollHeight;
  }
")


    }

  )

  incProgress(
    0.8,
    detail = "Finalizing..."
  )

})

#perform ftest based on selected method
  safe_NPARCtest <- function(metrics, dfType = "empirical") {
    tryCatch(
      NPARCtest(metrics, dfType = dfType),
      error = function(e) {
        message("Empirical df failed (", e$message, "). Falling back to theoretical df.")
        NPARCtest(metrics, dfType = "theoretical")
      }
    )
  }

fStats_empirical <- NPARCtest(
  fit_results$metrics,
  dfType = "empirical"
)

fStats_theoretical <- NPARCtest(
  fit_results$metrics,
  dfType = "theoretical"
)





str(fit_results$metrics)
summary(fit_results$metrics)


  # Store results in reactive values


fstats_empirical_rv(
  fStats_empirical
)

fstats_theoretical_rv(
  fStats_theoretical
)

metrics_rv(
  fit_results$metrics
)

predictions_rv(
  fit_results$predictions
)

})

#selector for fstats based on user input
selected_fstats <- reactive({

  req(
    fstats_empirical_rv(),
    fstats_theoretical_rv()
  )

  if (input$nparc_df_type == "theoretical") {

    return(
      fstats_theoretical_rv()
    )

  }

  fstats_empirical_rv()

})


# Build the filtered table used for display and selection
model_table <- reactive({
validate(
  need(
    !is.null(selected_fstats()),
    "Click 'Run model comparison' to generate results."
  ))

 
  fStats_to_use <- selected_fstats()
  req(fStats_to_use)


  # Ensure numeric cutoff
  padj_cut <- suppressWarnings(as.numeric(input$fstatpadjust))

  rss_cut <- suppressWarnings(
  as.numeric(
    input$rssdiff_cutoff
  )
)

if (is.na(rss_cut))
  rss_cut <- 0


  if (is.na(padj_cut)) padj_cut <- 1

metrics <- metrics_rv()

rss_table <- metrics %>%
  dplyr::filter(
    modelType == "null"
  ) %>%
  dplyr::select(
    id,
    rss_null = rss
  ) %>%
  dplyr::left_join(

    metrics %>%
      dplyr::filter(
        modelType == "alternative"
      ) %>%
      dplyr::group_by(id) %>%
      dplyr::summarise(
        rss_alt = sum(rss, na.rm = TRUE),
        .groups = "drop"
      ),

    by = "id"
  ) %>%
dplyr::mutate(
  rss_diff = rss_null - rss_alt
)

id_name <- names(fStats_to_use)[1]

tbl <- fStats_to_use %>%

  dplyr::rename(
    id = !!rlang::sym(id_name)
  ) %>%

  dplyr::filter(
    !is.na(pAdj)
  ) %>%

  dplyr::left_join(
    rss_table,
    by = "id"
  ) %>%
dplyr::filter(
  rss_diff > 0
)%>%
dplyr::filter(
  pAdj <= padj_cut
) %>%

dplyr::filter(
  rss_diff >= rss_cut
) %>%

dplyr::select(
  1,
  rss_null,
  rss_alt,
  rss_diff,
  fStat,
  pVal,
  pAdj
) %>%

  dplyr::arrange(
    dplyr::desc(fStat)
  )

  validate(need(nrow(tbl) > 0,
    "No proteins pass the current threshold. 
     Try increasing the pAdj cutoff."
  ))

  tbl
})


output$model_comparison_table <- DT::renderDataTable({
  tbl <- model_table()
  DT::datatable(
    tbl,
    options = list(pageLength = 10, scrollX = TRUE),
    rownames = FALSE,
    selection = "single"
  ) %>%
    DT::formatRound(
  columns = intersect(
    c(
      "rss_null",
      "rss_alt",
      "rss_diff",
      "fStat",
      "pVal",
      "pAdj"
    ),
    names(tbl)
  ),
  digits = 4
)
})




selected_id <- reactive({
  tbl <- model_table()
  s <- input$model_comparison_table_rows_selected
  validate(need(length(s) == 1, "Select a protein from the table to plot."))
  # First column is the ID
  id_value <- tbl[s, 1][[1]]
  id_value
})

# Prepare prediction curves for plotting

pred_curves <- reactive({
  req(predictions_rv(), selected_id())

  predictions_rv() %>%
    dplyr::filter(id == selected_id()) %>%
    dplyr::mutate(
      Temperature = x  # <-- FIX: rename the groups column
    ) %>%
    dplyr::group_by(modelType, Dataset, Temperature) %>%
    dplyr::summarise(fitted = mean(.fitted, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(modelType, Dataset, Temperature)
})




# =====================================================================
# Reusable function to build the model comparison plot
# =====================================================================
make_model_comparison_plot <- function() {

  req(combined_data_rv(), id_col_rv(), selected_id(), pred_curves())

  combined_data <- combined_data_rv()
  id_col        <- id_col_rv()
  selected_protein <- selected_id()
  pc <- pred_curves()

  plot_data <- combined_data %>%
    dplyr::filter(.data[[id_col]] == selected_protein)

  validate(need(nrow(plot_data) > 0, "No data available for the selected protein."))

  max_y <- max(plot_data$RelAbundance, na.rm = TRUE) * 1.1

  # Base plot
  p <- ggplot(plot_data, aes(
    x = Temperature,
    y = RelAbundance
  )) +
    geom_point(aes(shape = factor(Replicate), color = Dataset), size = 2) +
    labs(
      title = paste(
      "Model Comparison for Protein:",
      selected_protein,
      "|",
      input$nparc_df_type,
      "F-test"
    ),
      x = "Temperature",
      y = "Relative Abundance",
      color = "Dataset"
    ) +
    theme_light(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
    ) +
    ylim(0, max_y)

  # --------------------------------------------------------------------
  # Alternative model curves: one per dataset
  # --------------------------------------------------------------------
  p <- p + geom_line(
    data = pc %>% dplyr::filter(modelType == "alternative"),
    aes(
      x = Temperature,
      y = fitted,
      color = Dataset,
      group = Dataset
    ),
    linewidth = 1.1
  )

  # --------------------------------------------------------------------
  # Null model curve: single dashed black line
  # --------------------------------------------------------------------
  p <- p + geom_line(
    data = pc %>% dplyr::filter(modelType == "null"),
    aes(
      x = Temperature,
      y = fitted,
      linetype = "Null model"
    ),
    color = "black",
    linewidth = 1
  )

  # --------------------------------------------------------------------
  # Customize legend for null model
  # --------------------------------------------------------------------
  p <- p +
    scale_linetype_manual(
      name = NULL,
      values = c("Null model" = "dashed")
    ) +
    guides(
      color    = guide_legend(order = 1),
      linetype = guide_legend(
        order = 2,
        override.aes = list(color = "black")
      )
    )

  return(p)
}


# =====================================================================
# Render the plot in the UI
# =====================================================================
output$model_comparison_plot <- renderPlot({
  make_model_comparison_plot()
})

#Histogram pvalue
output$nparc_pvalue_hist <- renderPlot({

req(
  selected_fstats()
)

ggplot(
  selected_fstats(),
    aes(x = pVal)
  ) +

    geom_histogram(
      bins = 40,
      fill = "steelblue",
      color = "black"
    ) +

    theme_minimal() +

    labs(
      title = "NPARC p-value distribution",
      x = "p-value",
      y = "Count"
    )

})

output$nparc_fstat_hist <- renderPlot({

req(selected_fstats())

ggplot(
  selected_fstats(),
    aes(x = fStat)
  ) +

    geom_histogram(
      bins = 50,
      fill = "tomato",
      color = "black"
    ) +

    coord_cartesian(
      xlim = c(
        0,
        quantile(
          selected_fstats()$fStat,
          0.99,
          na.rm = TRUE
        )
      )
    ) +

    theme_minimal() +

    labs(
      title = "NPARC F statistic distribution",
      x = "F statistic",
      y = "Count"
    )

})
output$nparc_fdist_fit <- renderPlot({

  req(selected_fstats())

  fstats <- selected_fstats()
  fstats <- fstats %>%
  dplyr::filter(
    is.finite(fStat)
  )


  df1_use <- median(
    fstats$df1,
    na.rm = TRUE
  )

  df2_use <- median(
    fstats$df2,
    na.rm = TRUE
  )

xmax <- quantile(
  fstats$fStat[
    is.finite(fstats$fStat)
  ],
  0.99,
  na.rm = TRUE
)
  x_grid <- seq(
    0,
    xmax,
    length.out = 1000
  )

  theoretical_df <- data.frame(
    x = x_grid,
    density = df(
      x_grid,
      df1 = df1_use,
      df2 = df2_use
    )
  )

  ggplot(
    fstats,
    aes(
      x = fStat
    )
  ) +

    geom_density(
      aes(y = after_stat(density)),
      colour = "black",
      fill = "steelblue",
      alpha = 0.4,
      linewidth = 0.8
    ) +

    geom_line(
      data = theoretical_df,
      aes(
        x = x,
        y = density
      ),
      inherit.aes = FALSE,
      colour = "red",
      linewidth = 1.2
    ) +

    coord_cartesian(
      xlim = c(
        0,
        xmax
      )
    ) +

    theme_minimal() +

    labs(
      title = paste(
        "Observed vs Expected F Distribution",
        "(",
        input$nparc_df_type,
        ")"
      ),
      x = "F statistic",
      y = "Density"
    )

})
#RSS histogram
output$nparc_rssdiff_hist <- renderPlot({

  req(metrics_rv())

  rss_tbl <- metrics_rv() %>%

    dplyr::filter(
      modelType == "null"
    ) %>%

    dplyr::select(
      id,
      rss_null = rss
    ) %>%

    dplyr::left_join(

      metrics_rv() %>%

        dplyr::filter(
          modelType == "alternative"
        ) %>%

        dplyr::group_by(id) %>%

        dplyr::summarise(
          rss_alt = sum(rss),
          .groups = "drop"
        ),

      by = "id"
    ) %>%

    dplyr::mutate(
      rss_diff = rss_null - rss_alt
    )

  ggplot(
    rss_tbl,
    aes(x = rss_diff)
  ) +

    geom_histogram(
      bins = 60,
      fill = "steelblue",
      color = "black"
    ) +

    theme_minimal() +

    labs(
      title = "RSS difference distribution",
      x = "RSSnull - RSSalt",
      y = "Count"
    )

})

#nparc summary text
output$nparc_summary <- renderText({

  req(
    selected_fstats(),
    metrics_rv()
  )

  fstats <- selected_fstats()

  metrics <- metrics_rv()

  rss_tbl <- metrics %>%
    dplyr::filter(
      modelType == "null"
    ) %>%
    dplyr::select(
      id,
      rss_null = rss
    ) %>%
    dplyr::left_join(

      metrics %>%
        dplyr::filter(
          modelType == "alternative"
        ) %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(
          rss_alt = sum(rss),
          .groups = "drop"
        ),

      by = "id"
    ) %>%
    dplyr::mutate(
      rss_diff = rss_null - rss_alt
    )

  paste0(

    "Method: ",
    input$nparc_df_type,

    "\n\nProteins tested: ",
    nrow(fstats),
    "\nMedian df1: ",
round(
  median(
    fstats$df1,
    na.rm = TRUE
  ),
  2
),

"\nMedian df2: ",
round(
  median(
    fstats$df2,
    na.rm = TRUE
  ),
  2
),
"\nMaximum finite F statistic: ",
round(
  max(
    fstats$fStat[
      is.finite(fstats$fStat)
    ],
    na.rm = TRUE
  ),
  2
),

"\n99th percentile F statistic: ",
round(
quantile(
  fstats$fStat[
    is.finite(fstats$fStat)
  ],
  0.99,
  na.rm = TRUE
),
  2
),

    "\nSignificant (pAdj < 0.05): ",
    sum(
      fstats$pAdj < 0.05,
      na.rm = TRUE
    ),
    

    "\nSignificant (pAdj < 0.01): ",
    sum(
      fstats$pAdj < 0.01,
      na.rm = TRUE
    ),
    "\nMedian RSS difference: ",
      round(
        median(
          rss_tbl$rss_diff,
          na.rm = TRUE
        ),
        4
      ),

      "\n95th percentile RSS difference: ",
      round(
        quantile(
          rss_tbl$rss_diff,
          0.95,
          na.rm = TRUE
        ),
        4
      ),

    "\nNegative RSS differences: ",
    sum(
      rss_tbl$rss_diff < 0,
      na.rm = TRUE
    )

  )

})
observeEvent(input$nparc_help, {

  showModal(

    modalDialog(

      title = "NPARC Documentation",

HTML(
  '
  <b>NPARC diagnostic plots and summary</b><br><br>

  NPARC (Non-Parametric Analysis of Response Curves) compares complete protein melting curves between experimental conditions rather than relying solely on melting temperature (Tm) shifts.<br><br>

  The method evaluates whether fitting separate curves for each condition significantly improves the fit compared with a shared model, allowing the detection of treatment effects that may not produce large Tm changes.<br><br>

  <b>Plots included in this section</b><br><br>

  • <b>P-value distribution</b>: overall significance landscape across all tested proteins.<br>

  • <b>F-statistic distribution</b>: strength of improvement when fitting condition-specific models.<br>

  • <b>Observed vs expected F-distribution</b>: comparison between the observed test statistics and the theoretical expectation.<br>

  • <b>RSS difference distribution</b>: magnitude of fit improvement between the null and alternative models.<br>

  • <b>NPARC summary</b>: overall statistics and significant hit counts.<br><br>

  These diagnostics help assess model quality, significance thresholds, effect size distributions and the overall reliability of the analysis.<br><br>

  <b>Reference</b><br><br>

  Childs D, Bach K, Franken H, Anders S, Kurzawa N, Bantscheff M, Savitski MM (2019).<br>

  <i>Nonparametric analysis of thermal proteome profiles reveals novel drug-binding proteins.</i><br>

  Molecular & Cellular Proteomics 18(12):2506-2515.<br><br>

  <b>Documentation</b><br><br>

  Copy and paste this URL into your browser:<br><br>

  bioconductor.posit.co/packages/devel/bioc/vignettes/NPARC/inst/doc/NPARC.html
  '
),

      easyClose = TRUE,
      size = "l"

    )

  )

})
observeEvent(input$tpp_help, {

  showModal(

    modalDialog(

      title = "TPP Documentation",

      HTML(
        '
        <b>Thermal Proteome Profiling (TPP)</b><br><br>

        TPP identifies proteins whose thermal stability changes between vehicle and treatment conditions, indicating potential protein-ligand interactions or downstream biological responses.<br><br>

        This implementation uses the Bioconductor <b>TPP</b> package and follows the TPP-TR workflow.<br><br>

        <b>Plots included in this section</b><br><br>

        • Adjusted p-value distribution<br>

        • Delta Tm distribution<br>

        • Protein-specific melting curves<br>

        • Summary statistics<br><br>

        <b>Reference</b><br><br>

        Franken H et al. (2015).<br>

        Thermal Proteome Profiling for unbiased identification of direct and indirect drug targets using multiplexed quantitative mass spectrometry.<br>

        Nature Protocols 10:1567–1593.<br><br>

        <b>Documentation</b><br><br>

        Copy and paste this URL into your browser:<br><br>

        bioconductor.org/packages/TPP
        '
      ),

      easyClose = TRUE,
      size = "l"

    )

  )

})


# =====================================================================
# Download handler
# =====================================================================
output$download_model_comparison_plot <- downloadHandler(
  filename = function() {
    paste0("model_comparison_", selected_id(), "_", Sys.Date(), ".png")
  },
  content = function(file) {
    p <- make_model_comparison_plot()
    ggsave(
      filename = file,
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )
  }
)

output$download_nparc_table <- downloadHandler(

  filename = function() {

    paste0(
      "NPARC_results_",
      Sys.Date(),
      ".csv"
    )

  },

  content = function(file) {

    write.csv(
      model_table(),
      file,
      row.names = FALSE
    )

  }

)

output$download_tpp_table <- downloadHandler(

  filename = function() {

    paste0(
      "TPP_results_",
      Sys.Date(),
      ".csv"
    )

  },

  content = function(file) {

    write.csv(
      tpp_display_table(),
      file,
      row.names = FALSE
    )

  }

)
output$download_tpp_curve_plot <- downloadHandler(

  filename = function() {

    paste0(
      "TPP_curve_",
      selected_tpp_protein(),
      "_",
      Sys.Date(),
      ".png"
    )

  },

  content = function(file) {

    p <- make_tpp_curve_plot()

    ggsave(
      filename = file,
      plot = p,
      width = 10,
      height = 6,
      dpi = 300
    )

  }

)

})




############################################################

# #normalization
output$norm_tabs_ui <- renderUI({
  req(processed$data())

  datasets <- processed$data()

  tabs <- lapply(names(datasets), function(name) {
    tabPanel(
      title = name,
      value = name,
      normalizationPanelModuleUI(name)
    )
  })

  do.call(tabsetPanel, c(
    list(id = "norm_dataset_tabs"),
    tabs
  ))
})
# === Clustering tabs ===
output$clustering_tabs_ui <- renderUI({
  req(processed$data())

  datasets <- processed$data()

  tabs <- lapply(names(datasets), function(name) {

    cluster_id <- paste0(name, "_", isolate(cluster_version()))

    tabPanel(
      title = name,
      value = cluster_id,
      clusterPanelModuleUI(cluster_id, name)
    )
  })

  do.call(tabsetPanel, c(
    list(id = "clust_dataset_tabs"),
    tabs
  ))
})


#modelfitting UI
output$model_fitting_tabs_ui <- renderUI({
  req(processed$data())

  datasets <- processed$data()

tabs <- lapply(names(datasets), function(name) {

  safe_id <- gsub(
    "[^[:alnum:]_.+-]",
    "_",
    name
  )

  tabPanel(
    title = name,
    value = name,
    modelPanelModuleUI(
      safe_id,
      name
    )
  )
})

  do.call(tabsetPanel, c(
    list(id = "model_dataset_tabs"),
    tabs
  ))
})

#protein plots UI
output$protein_plots_tabs_ui <- renderUI({
  req(processed$data())

  datasets <- processed$data()

  tabs <- lapply(names(datasets), function(name) {
    tabPanel(
      title = name,
      value = name,
      protplotsPanelModuleUI(name)
    )
  })

  do.call(tabsetPanel, c(
    list(id = "protein_plots_tabs"),
    tabs
  ))
})


 #store the clicked data
 plot_data <- reactiveVal(NULL)
 plot_trigger <- reactiveVal(FALSE)


#protein plots module

output$go_terms_tabs_ui <- renderUI({

  req(processed$data())

  datasets <- processed$data()

tabs <- lapply(names(datasets), function(name) {

  safe_id <- gsub(
    "[^[:alnum:]_.+-]",
    "_",
    name
  )

tabPanel(
  title = name,
  value = name,

  tabsetPanel(

tabPanel(

  title = tagList(

    "ORA",

    tags$span(
      `data-toggle` = "tooltip",
      `data-placement` = "right",
      `data-html` = "true",

      title = HTML(
        "<b>GO Over-Representation Analysis (ORA)</b><br><br>

        GO ORA identifies Gene Ontology (GO) terms that occur more frequently than expected within a selected protein cluster.<br><br>

        In this application, ORA is performed using proteins assigned to each cluster during the clustering step.<br><br>

        Results can help identify biological processes, molecular functions, and cellular components that are characteristic of each cluster."
      ),

      tags$i(
        class = "fa fa-question-circle",
        style = "margin-left:8px; color:#337ab7; cursor:pointer;"
      )
    )
  ),

  goORAClusterPanelUI(safe_id),

  br(),

  goORAPlotsPanelUI(safe_id)

),

tabPanel(

  title = tagList(

    "GSEA",

    tags$span(
      `data-toggle` = "tooltip",
      `data-placement` = "right",
      `data-html` = "true",

      title = HTML(
        "<b>Gene Set Enrichment Analysis (GSEA)</b><br><br>

        GSEA identifies biological processes that are enriched across the full ranked protein list rather than only among significantly changing proteins.<br><br>

        Protein ranking is obtained directly from the differential analysis results generated in the Limma section of the application.<br><br>

        Depending on the selected ranking metric, proteins are ordered using either:<br>
        • Limma moderated t-statistics (recommended)<br>
        • Limma log fold-change values"
      ),

      tags$i(
        class = "fa fa-question-circle",
        style = "margin-left:8px; color:#337ab7; cursor:pointer;"
      )
    )
  ),

  goGSEAPanelModuleUI(safe_id)

)

  )
)

  })

  do.call(tabsetPanel, c(
    list(id = "go_terms_dataset_tabs"),
    tabs
  ))

})








process_dataset <- function(df, selected_proteins, Org_data,
                            dataset_label, line_type, debug_tag = "") {

  cat(">> process_dataset START for:", dataset_label, "\n")

  # --------------------------
  # Basic checks
  # --------------------------
  if (is.null(selected_proteins) || !is.data.frame(selected_proteins) || nrow(selected_proteins) == 0) {
    cat("!! selected_proteins empty for", dataset_label, "\n")
    return(NULL)
  }

  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
    cat("!! df NULL/empty for", dataset_label, "\n")
    return(NULL)
  }

  cat("Columns in df:", paste(names(df), collapse = ", "), "\n")

  # --------------------------
  # Ensure ProteinID exists
  # --------------------------
  if (!"ProteinID" %in% names(df)) {
    names(df)[1] <- "ProteinID"
  }

  # --------------------------
  # Ensure Cluster column
  # --------------------------
  if (!"Cluster" %in% names(df)) {
    df$Cluster <- NA_character_
  }

  # --------------------------
  # Extract selected IDs
  # --------------------------
  ids <- NULL
  if ("ProteinID" %in% names(selected_proteins)) {
    ids <- selected_proteins$ProteinID
  } else {
    ids <- selected_proteins[[1]]
  }

  ids <- unique(as.character(ids))
  cat("Selected IDs:", paste(head(ids, 5), collapse = ", "), "\n")

  # --------------------------
  # Org data
  # --------------------------
  if (!is.null(Org_data) && is.data.frame(Org_data) && nrow(Org_data) > 0) {
    organism_df <- Org_data %>%
      dplyr::rename(ProteinID = 1) %>%   #  explicit
      dplyr::mutate(
        ProteinID = sub("^.*\\|(.*)\\|.*$", "\\1", ProteinID)  #  strip UniProt formats
      ) %>%
      dplyr::select(ProteinID, Protein.names, Gene.Names)

  } else {
    organism_df <- tibble::tibble(
      ProteinID     = df$ProteinID,
      Protein.names = NA_character_,
      Gene.Names    = NA_character_
    )
  }

  df <- df %>% dplyr::left_join(organism_df, by = "ProteinID")
cat("Joined organism_df rows:", nrow(organism_df), "\n")
cat("Matched proteins:", sum(!is.na(df$Protein.names)), "\n")
  # --------------------------
  #  Detect value columns
  # --------------------------
  fixed_cols <- intersect(c("ProteinID", "Protein.names", "Gene.Names", "Cluster"), names(df))
  value_cols <- setdiff(names(df), fixed_cols)

  # Fallback if empty (IMPORTANT FIX)
  if (length(value_cols) == 0) {
    cat("!! No value_cols found → fallback numeric detection\n")

    value_cols <- names(df)[sapply(df, is.numeric)]

    if (length(value_cols) == 0) {
      cat("!! Still no usable columns → abort\n")
      return(NULL)
    } else {
      cat("Using numeric columns:", paste(value_cols, collapse = ", "), "\n")
    }
  }

  # --------------------------
  #  Build long data
  # --------------------------
  df_long_list <- lapply(ids, function(pid) {

    dfi <- df %>% dplyr::filter(.data$ProteinID == pid)

prot_name <- dfi$Protein.names[1]
gene_name <- dfi$Gene.Names[1]

    if (nrow(dfi) == 0) {
      cat("!! ID not found:", pid, "\n")
      return(NULL)
    }

    # ✅ pivot safely
    out <- tryCatch({
      dfi %>%
        tidyr::pivot_longer(
          cols = dplyr::all_of(value_cols),
          names_to = "Replicate",
          values_to = "Intensity"
        )
    }, error = function(e) {
      cat("!! pivot failed for", pid, "\n")
      return(NULL)
    })

    if (is.null(out)) return(NULL)

    out %>%
      tidyr::separate(Replicate, into = c("temperature", "Replicate"),
                      sep = "_", remove = FALSE, fill = "right") %>%
      dplyr::group_by(temperature) %>%
      dplyr::mutate(
        Mean     = mean(Intensity, na.rm = TRUE),
        SD       = sd(Intensity, na.rm = TRUE),
        Protein  = pid,
        Protein.names = prot_name,
        Gene.Names    = gene_name,
        Dataset  = dataset_label,
        LineType = line_type
      ) %>%
      dplyr::ungroup()
  })

  df_long <- dplyr::bind_rows(df_long_list)

  # --------------------------
  #  Final check
  # --------------------------
  if (is.null(df_long) || !is.data.frame(df_long) || nrow(df_long) == 0) {
    cat("!! df_long empty for", dataset_label, "\n")
    return(NULL)
  }

  # --------------------------
  # Tooltip
  # --------------------------
  df_long <- df_long %>%
    dplyr::mutate(
      tooltip = paste0(
        "Protein: ", Protein,
        "<br>Gene: ", ifelse(is.na(Gene.Names), "NA", Gene.Names),
        "<br>Name: ", ifelse(is.na(Protein.names), "NA", Protein.names),
        "<br>Dataset: ", Dataset,
        "<br>Temperature: ", temperature,
        "<br>Mean: ", round(Mean, 2),
        "<br>SD: ", round(SD, 2)
      )
    )

  cat(" process_dataset SUCCESS → rows:", nrow(df_long), "\n")

  return(df_long)
}





#plot
output$protein_plot <- plotly::renderPlotly({

  df_long <- plot_data()



  if (identical(df_long, "no_selection") ||
      is.null(df_long) ||
      !is.data.frame(df_long) ||
      nrow(df_long) == 0) {
    return(NULL)
  }

  req(df_long)
  if (nrow(df_long) == 0) return(NULL)

  max_y <- max(df_long$Mean + df_long$SD, na.rm = TRUE) * 1.1

  # detect number of datasets
  n_datasets <- length(unique(df_long$Dataset))

  # --------------------------
  # CONDITIONAL MAPPING
  # --------------------------
  if (n_datasets == 1) {

    #  ONE DATASET → color by protein, same line style
    p <- ggplot(df_long, aes(
      x = as.numeric(temperature),
      y = Mean,
      color = Protein,
      linetype = Dataset,
      group = Protein,
      text = tooltip
    ))

  } else {

    # MULTIPLE DATASETS → color + linetype by dataset
    p <- ggplot(df_long, aes(
      x = as.numeric(temperature),
      y = Mean,
      color = Dataset,
      linetype = Dataset,
      group = interaction(Protein, Dataset),
      text = tooltip
    ))

  }

  # --------------------------
  # BUILD PLOT
  # --------------------------
  p <- p +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), width = 1) +
    geom_line(size = 0.8) +
    scale_linetype_manual(
      values = setNames(
        c("solid", "dashed", "dotted", "dotdash", "longdash")[seq_along(unique(df_long$Dataset))],
        unique(df_long$Dataset)
      )
    ) +
    labs(
      x = "Temperature",
      y = "Mean Abundance ± SD"
    ) +
    theme_light(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(color = "grey80", fill = NA, linewidth = 0.5)
    ) +
    ylim(0, max_y)

  # --------------------------
  # LEGEND CONTROL
  # --------------------------
  if (n_datasets == 1) {

    p <- p +
      guides(
        color = guide_legend(title = "Protein"),
        linetype = "none"
      )

  } else {

    p <- p +
      guides(
        color = guide_legend(title = "Dataset"),
        linetype = guide_legend(title = "Dataset")
      )

  }

  # --------------------------
  # ✅ RETURN PLOTLY
  # --------------------------
  plotly::ggplotly(p, tooltip = "text")

})

observeEvent(plot_trigger(), {
  # Reset the trigger after a short delay
  later::later(function() {
    plot_trigger(FALSE)
  }, delay = 0.5)  # Adjust delay if needed
})




###########################
# Model Comparison tab
output$model_comparison_ui <- renderUI({

  req(processed$data())

  available_conditions <- names(processed$data())

  if (length(available_conditions) < 2) {

    return(
    tags$strong("At least two conditions are required.")
    )

  }

  tagList(
      fluidRow(

            box(
              title = tagList(

                "Model control",

                tags$span(
                  `data-toggle` = "tooltip",
                  `data-placement` = "right",
                  `data-html` = "true",

                  title = HTML(
                    "<b>NPARC: Non-Parametric Analysis of Response Curves</b><br><br>

                    NPARC compares complete melting curves between two experimental conditions instead of relying solely on melting temperature (Tm) shifts.<br><br>

                    The method evaluates whether fitting separate models for each condition significantly improves the fit compared with a shared model.<br><br>"
                  ),

                  tags$i(
                    class = "fa fa-question-circle",
                    style = "margin-left:8px; color:#337ab7; cursor:pointer;"
                  )
                )

              ),
            width = 4,
            status = "primary",
            solidHeader = TRUE,

            selectInput(
              "comparison_condition_1",
              "Condition 1",
              choices = available_conditions,
              selected = available_conditions[1]
            ),

            selectInput(
              "comparison_condition_2",
              "Condition 2",
              choices = available_conditions,
              selected = available_conditions[min(2, length(available_conditions))]
            ),

            br(),

            actionButton(
              "run_model_comparison",
              "Run model comparison"
            ),
            br(),
            br(),

            tags$pre(
              id = "nparc_console",
              style = "
                background-color:#111;
                color:#00ff00;
                height:150px;
                overflow-y:auto;
                padding:10px;
                font-size:12px;
              "
            ),

            tags$hr(),

            selectInput(
              inputId = "nparc_df_type",
              label = "F statistic method",
              choices = c(
                "empirical",
                "theoretical"
              ),
              selected = "empirical"
            ),

            numericInput(
              "fstatpadjust",
              "Filter table pAdj values <= to:",
              value = 0.01,
              min = 0.0001,
              step = 0.1
            ),

            numericInput(
              "rssdiff_cutoff",
              "Minimum RSS difference:",
              value = 0,
              min = 0,
              step = 0.05
            )
          ),

            box(
              title = "Model comparison results",
              width = 8,
              status = "primary",
              solidHeader = TRUE,

              DT::dataTableOutput(
                "model_comparison_table"
              ),
              br(),

              downloadButton(
                "download_nparc_table",
                "Download Results"
              )
            )
      ),
fluidRow(

  box(
    title = "Model comparison plot",
    width = 12,
    status = "info",
    solidHeader = TRUE,

    plotOutput(
      "model_comparison_plot",
      height = 650
    ) %>% withSpinner(
color = "#0EA5A5"
),

    downloadButton(
      "download_model_comparison_plot",
      "Download Plot"
    )
  )

),
fluidRow(

    box(
      title = tagList(

        "NPARC summary",

        actionLink(
          "nparc_help",

          label = NULL,

    icon = icon("book-open"),

    style = "
      margin-left:8px;
      color:var(--text-secondary);
      display:inline;
    "
        )

      ),
    width = 6,
    status = "success",
    solidHeader = TRUE,

    verbatimTextOutput(
      "nparc_summary"
    )
  ),

  column(

    width = 6,

    fluidRow(

      box(
        title = "P-value distribution",
        width = 6,
        status = "warning",
        solidHeader = TRUE,

      plotOutput(
        "nparc_pvalue_hist",
        height = 220
      ) %>% withSpinner(
color = "#0EA5A5"
)
      ),

      box(
        title = "F statistic distribution",
        width = 6,
        status = "warning",
        solidHeader = TRUE,

        plotOutput(
          "nparc_fstat_hist",
          height = 220
        ) %>% withSpinner(
color = "#0EA5A5"
)

      )

    ),

    fluidRow(

      box(
        title = "F distribution fit",
        width = 6,
        status = "warning",
        solidHeader = TRUE,

        plotOutput(
          "nparc_fdist_fit",
          height = 220
        )%>% withSpinner(
color = "#0EA5A5"
)
      ),

      box(
        title = "RSS difference distribution",
        width = 6,
        status = "warning",
        solidHeader = TRUE,

        plotOutput(
          "nparc_rssdiff_hist",
          height = 220
        )%>% withSpinner(
color = "#0EA5A5"
)
      )

    )

  )

)
      )


})

output$tpp_analysis_ui <- renderUI({

  req(processed$data())

  available_conditions <- names(processed$data())

  fluidPage(

    fluidRow(

      box(
        title = tagList(

          "TPP analysis",

          tags$span(
            `data-toggle` = "tooltip",
            `data-placement` = "right",
            `data-html` = "true",

            title = HTML(
              "<b>Thermal Proteome Profiling (TPP)</b><br><br>

              TPP compares protein thermal stability between vehicle and treatment conditions.

              Compares the melting temperature point (Tm) of individual replicates.

              Statistical analysis is performed using the Bioconductor TPP package."
            ),

      tags$i(
        class = "fa fa-question-circle",
        style = "margin-left:8px; color:#337ab7; cursor:pointer;"
      )
          )

        ),
        width = 4,
        status = "primary",
        solidHeader = TRUE,

        selectInput(
          "tpp_vehicle",
          "Vehicle condition",
          choices = available_conditions,
          selected = available_conditions[1]
        ),

        selectInput(
          "tpp_treatment",
          "Treatment condition",
          choices = available_conditions,
          selected = available_conditions[min(2,length(available_conditions))]
        ),
        numericInput(
        "tpp_min_r2",
        "Minimum R²",
        value = 0.8,
        min = 0,
        max = 1,
        step = 0.05
      ),

      numericInput(
        "tpp_max_plateau",
        "Maximum Plateau",
        value = 0.3,
        min = 0,
        max = 2,
        step = 0.05
      ),

checkboxInput(
  "tpp_show_hits_only",
  "Show only proteins passing all replicate comparisons",
  FALSE
),
checkboxInput(
  "tpp_show_tm_only",
  "Show only proteins with ΔTm calculated in at least one replicate",
  FALSE
),
        actionButton(
          "run_tpp_import",
          "Run TPP analysis"
        )

      ),

      box(
        title = "TPP output",
        width = 8,
        status = "info",
        solidHeader = TRUE,

        DT::dataTableOutput(
          "tpp_results_table"
        ),
        br(),

        downloadButton(
          "download_tpp_table",
          "Download Results"
        )

      )

    ),
fluidRow(

  box(
    title = "TPP melting curve",
    width = 12,
    status = "primary",
    solidHeader = TRUE,

    plotOutput(
      "tpp_curve_plot",
      height = 650
    ) %>% withSpinner(
      color = "#0EA5A5"
    ),

    br(),

    downloadButton(
      "download_tpp_curve_plot",
      "Download Plot"
    )
  )

),
fluidRow(

  box(
    title = tagList(

      "TPP summary",

      actionLink(
        "tpp_help",
        label = NULL,
    icon = icon("book-open"),

    style = "
      margin-left:8px;
      color:var(--text-secondary);
      display:inline;
    "
      )

    ),

    width = 6,
    status = "success",
    solidHeader = TRUE,

    verbatimTextOutput(
      "tpp_summary"
    )
  ),

  column(

    width = 6,

    fluidRow(

      box(
        title = "Adjusted p-value distribution",
        width = 12,
        status = "warning",
        solidHeader = TRUE,

plotOutput(
  "tpp_padj_hist",
  height = 220
) %>% withSpinner(
color = "#0EA5A5"
)
      )

    ),

    fluidRow(

      box(
        title = "Delta Tm distribution",
        width = 12,
        status = "warning",
        solidHeader = TRUE,

      plotOutput(
        "tpp_dtm_hist",
        height = 220
      ) %>% withSpinner(
color = "#0EA5A5"
)
      )

    )

  )

)
  )
})



}



# ###########################################################################






# }


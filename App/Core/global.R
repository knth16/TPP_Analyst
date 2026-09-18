
#packages
library(shiny)
library(shinydashboard)
library(readr)
library(bslib)
library(bsicons)
library(tidyverse)
library(dplyr)
library(DT)
library(ComplexHeatmap)
library(grid)
library(pheatmap)
library(DEP)
library(shinyjs)
library(limma)
library(tibble)
library(minpack.lm)
library(pracma)
library(UniprotR)
library(purrr)
library(clusterProfiler)
library(enrichplot)
library(GO.db)
library(ggtext)
library(plotly)
library(Amelia)
library(ggnewscale)
library(VIM)
library(shinycssloaders)
library(shinyjs)
library(later)
library(shinybusy)
library(NPARC)
library(broom)
library(knitr)
library(callr)
library(processx)
library(TPP)

#####################################################################################
#functions


clean_names_safe <- function(df) {

  nm <- names(df)

  # replace NA or "" with safe names
  bad <- is.na(nm) | nm == ""
  nm[bad] <- paste0("V", seq_len(sum(bad)))

  # trim spaces
  nm <- trimws(nm)

  names(df) <- nm

  df
}

#ora launch

run_ORA <- function(
  dataclustered,
  term2g,
  pval,
  qval,
  padj,
  minGS,
  maxGS
) {
  library(clusterProfiler)
  library(GO.db)

  Test_data <- dataclustered[, c("Accession", "Cluster")]
  background_genes <- unique(Test_data$Accession)

  cluster_list <- split(Test_data$Accession, Test_data$Cluster)

  GOMF <- as.list(GOMFANCESTOR)
  GOCC <- as.list(GOCCANCESTOR)
  GOBP <- as.list(GOBPANCESTOR)

  
get_GO_term_name <- function(go_id) {

  obj <- GOTERM[[go_id]]

  if (is.null(obj)) {
    return(NA_character_)
  }

  tryCatch(
    Term(obj),
    error = function(e) {
      NA_character_
    }
  )
}


  get_ontology <- function(go_id) {
    if (!is.null(GOMF[[go_id]])) "MF"
    else if (!is.null(GOCC[[go_id]])) "CC"
    else if (!is.null(GOBP[[go_id]])) "BP"
    else NA
  }

  lapply(names(cluster_list), function(cluster_name) {
    result <- enricher(
      gene          = cluster_list[[cluster_name]],
      universe      = background_genes,
      TERM2GENE     = term2g,
      pvalueCutoff  = pval,
      qvalueCutoff  = qval,
      pAdjustMethod = padj,
      minGSSize     = minGS,
      maxGSSize     = maxGS
    )

    if (!is.null(result) && nrow(result@result) > 0) {
      result@result$Cluster     <- cluster_name
      result@result$Description <- sapply(result@result$ID, get_GO_term_name)
      result@result$ONTOLOGY    <- sapply(result@result$ID, get_ontology)
      result
    } else {
      NULL
    }
  })
}





ensure_protein_df <- function(df) {
  if (is.null(df) || !is.data.frame(df) || nrow(df) == 0)
    return(NULL)

  df <- as.data.frame(df)

  protein_id_from <- NULL

  # CASE 1 — normalized dataset: use the Accession column
  if ("Accession" %in% names(df)) {
    df$ProteinID <- as.character(df$Accession)
    protein_id_from <- "Accession"
  }
  # CASE 2 — if no Accession, use rownames
  else if (!is.null(rownames(df))) {
    df$ProteinID <- as.character(rownames(df))
    protein_id_from <- "rownames"
  }
  # CASE 3 — fallback: first column
  else {
    df$ProteinID <- as.character(df[[1]])
    protein_id_from <- names(df)[1]
  }

  # Put ProteinID FIRST
  other_cols <- setdiff(names(df), "ProteinID")
  df <- df[, c("ProteinID", other_cols), drop = FALSE]

  # Remove Accession if it was used as ID
  if (protein_id_from == "Accession") {
    df$Accession <- NULL
  }

  # Ensure Cluster exists
  if (!"Cluster" %in% names(df)) {
    df$Cluster <- NA_character_
  }

  df
}

clean_protein_ids <- function(ids) {

  ids <- as.character(ids)

  # Remove legacy IPI suffixes
  ids <- sub("_IPI.*", "", ids)

  ids

}
#QC functions
summarize_zeros <- function(df) {
  df %>%
    summarise(across(everything(), ~ sum(. == 0, na.rm = TRUE))) %>%
    pivot_longer(cols = everything(), names_to = "Column", values_to = "Zero_Count") %>%
    filter(row_number() != 1) %>%
    separate(Column, into = c("temperature", "Replicate"), sep = "_") %>%
    mutate(ID = str_c(temperature, Replicate, sep = "_"))
}

plot_detect_custom_df <- function(df, colors = c("#1f77b4", "#ff7f0e"), font_family = "Arial") {
  stopifnot("Accession" %in% colnames(df))
  numeric_data <- df[, sapply(df, is.numeric), drop = FALSE]
  if (!any(numeric_data == 0, na.rm = TRUE)) stop("No zero values in the input dataframe")


  df_long <- numeric_data %>%
    mutate(rowname = row_number()) %>%
    tidyr::pivot_longer(-rowname, names_to = "ID", values_to = "val")


  stat <- df_long %>%
    group_by(rowname) %>%
    summarise(mean = mean(val[val != 0], na.rm = TRUE), zeroval = any(val == 0, na.rm = TRUE))
cat("\n====================\n")
cat("DENSITY DEBUG\n")
cat("====================\n")

print(table(stat$zeroval))
  ggplot(stat, aes(mean, color = zeroval)) +
    geom_density(na.rm = TRUE) +
    scale_color_manual(values = colors) +
    labs(x = expression(log[2] ~ "Intensity"), y = "Density") +
    guides(color = guide_legend(title = "Zero values")) +
    theme_DEP1() +
    theme(text = element_text(family = font_family))
}
#GO
# Function to fetch the data
GetProteinGOInfoM <- function(ProteinAccList, directorypath = NULL, progress_callback = NULL) {
  if (!curl::has_internet()) {
    message("Please connect to the internet.")
    return()
  }

  options(timeout = 10000)
  ProteinInfoParsed_total <- data.frame()
  baseUrl <- "https://rest.uniprot.org/uniprotkb/search?query=accession:"
  columns <- "go_id,go,go_p,go_f,go_c"

  for (i in seq_along(ProteinAccList)) {
    ProteinAcc <- ProteinAccList[i]
    ProteinName_url <- paste0(ProteinAcc, "&format=tsv&fields=", columns)
    RequestUrl <- paste0(baseUrl, ProteinName_url)
    RequestUrl <- URLencode(RequestUrl)

    Request <- tryCatch({
      httr::GET(RequestUrl, httr::timeout(7))
    }, error = function(cond) {
      message("Internet connection problem.")
      message(cond)
      return(NULL)
    })

    if (is.null(Request) || Request$status_code != 200) {
      next
    }

    ProteinDataTable <- tryCatch({
      read.csv(RequestUrl, header = TRUE, sep = "\t")
    }, error = function(e) NULL)

    if (!is.null(ProteinDataTable)) {
      ProteinDataTable <- ProteinDataTable[1, ]
      ProteinInfoParsed <- as.data.frame(ProteinDataTable, row.names = ProteinAcc)
      ProteinInfoParsed_total <- rbind(ProteinInfoParsed_total, ProteinInfoParsed)
    }

    # Update Shiny progress bar
    if (!is.null(progress_callback)) {
      progress_callback(i / length(ProteinAccList))
    }
  }

  if (!is.null(directorypath)) {
    write.csv(ProteinInfoParsed_total, file.path(directorypath, "Protein GO Info.csv"))
  }

  return(ProteinInfoParsed_total)
}


#########################################################################################################
#modules
#data loading
#server
proteindataModuleServer <- function(
  id,
  fileInput,
  requireOrganismColumn,
  exampleData = NULL
) {
  moduleServer(id, function(input, output, session) {


data <- reactive({

  # ---------------------------------
  # Example dataset
  # ---------------------------------

  if (
    !is.null(exampleData) &&
    !is.null(exampleData())
  ) {

    return(
      exampleData()
    )

  }

  # ---------------------------------
  # Uploaded dataset
  # ---------------------------------

  req(fileInput())

  path <- fileInput()$datapath
  ext  <- tolower(tools::file_ext(fileInput()$name))
  supported_ext <- c("csv", "tsv", "txt", "xlsx", "xls")

  # Guard: supported extension
  if (!ext %in% supported_ext) {
    showModal(
      modalDialog(
        title = "❌ Unsupported file",
        paste("Unsupported file format:", ext,
              "\nPlease upload CSV, TSV, TXT, or Excel files."),
        easyClose = TRUE
      )
    )
    req(NULL)   # <- abort silently; no crash
  }

  # --- read file with error guard
  df <- tryCatch(
    switch(
      ext,
      csv  = read.csv(path, header = TRUE, check.names = FALSE),
      tsv  = read.delim(path, header = TRUE, check.names = FALSE),
      txt  = read.delim(path, header = TRUE, check.names = FALSE),
      xlsx = as.data.frame(readxl::read_excel(path)),
      xls  = as.data.frame(readxl::read_excel(path))
    ),
    error = function(e) {
      showModal(
        modalDialog(
          title = "❌ Could not read file",
          paste("Error while reading the file:", conditionMessage(e)),
          easyClose = TRUE
        )
      )
      return(NULL)
    }
  )

  # If read failed, abort reactive
  if (is.null(df)) req(FALSE)
cat("\n====================\n")
cat("RAW IMPORT CHECK\n")
cat("====================\n")

cat("NA count:\n")
print(sum(is.na(df)))

cat("NaN count:\n")
print(sum(is.nan(as.matrix(df)), na.rm = TRUE))

cat("Infinite count:\n")
print(sum(is.infinite(as.matrix(df)), na.rm = TRUE))
  # Drop accidental leading empty column (your testdata.csv had this)
  if (identical(names(df)[1], "")) {
    df <- df[, -1, drop = FALSE]
  }

# ✅ Normalize names
names(df) <- trimws(names(df))

# ✅ REMOVE EMPTY / NA COLUMN NAMES (CRITICAL FIX)
valid_cols <- !(is.na(names(df)) | names(df) == "")
df <- df[, valid_cols, drop = FALSE]

# ✅ REMOVE COMPLETELY EMPTY COLUMNS
df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]


  # ---- STRUCTURE CHECKS ----
  required_cols <- c("Accession")

if (requireOrganismColumn()) {

  required_cols <- c(
    required_cols,
    "Organism"
  )

}
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    showModal(
      modalDialog(
        title = "❌ Incorrect data structure",
if (requireOrganismColumn()) {

  msg <- paste0(
    "Multiple organisms were selected.\n\n",
    "Required columns:\n",
    "• Accession\n",
    "• Organism\n\n",
    "Missing columns:\n• ",
    paste(missing_cols, collapse = "\n• ")
  )

} else {

  msg <- paste0(
    "The protein file must contain:\n",
    "• Accession\n\n",
    "Missing columns:\n• ",
    paste(missing_cols, collapse = "\n• ")
  )

},
        easyClose = TRUE
      )
    )
    
    req(NULL)
  }

# ----------------------------
# Validate unique Accession IDs
# ----------------------------

dup_ids <- unique(
  df$Accession[
    duplicated(df$Accession)
  ]
)

if (length(dup_ids) > 0) {

  showModal(
    modalDialog(
      title = "❌ Duplicate protein IDs detected",

      paste0(
        "The Accession column must contain unique identifiers.\n\n",
        "Found ",
        length(dup_ids),
        " duplicated accession(s).\n\n",
        "Examples:\n",
        paste(
          head(dup_ids, 10),
          collapse = "\n"
        )
      ),

      easyClose = TRUE
    )
  )

  req(FALSE)

}

  # ---- reorder columns ----
  base_cols <- intersect(
  c("Accession", "Organism"),
  names(df)
)

df <- df[
  ,
  c(
    base_cols,
    setdiff(names(df), base_cols)
  ),
  drop = FALSE
]
  df
})


    return(list(data = data))
  })
}



designdataModuleServer <- function(
  id,
  fileInput,
  proteinInput,
  customOrgData,
  selected_organisms,
  exampleData = NULL
) {
  moduleServer(id, function(input, output, session) {
org_file_cache <- reactiveValues()
data <- reactive({

  # ---------------------------------
  # Example dataset
  # ---------------------------------

  if (
    !is.null(exampleData) &&
    !is.null(exampleData())
  ) {

    df <- exampleData()

  } else {

    req(fileInput())

    path <- fileInput()$datapath
    ext  <- tolower(tools::file_ext(fileInput()$name))

    supported_ext <- c(
      "csv",
      "tsv",
      "txt",
      "xlsx",
      "xls"
    )

    if (!ext %in% supported_ext) {

      showModal(
        modalDialog(
          title = "❌ Unsupported file",
          paste(
            "Unsupported file format:",
            ext
          ),
          easyClose = TRUE
        )
      )

      req(FALSE)

    }

    df <- tryCatch(
      switch(
        ext,
        csv  = read.csv(path, header = TRUE, check.names = FALSE),
        tsv  = read.delim(path, header = TRUE, check.names = FALSE),
        txt  = read.delim(path, header = TRUE, check.names = FALSE),
        xlsx = as.data.frame(readxl::read_excel(path)),
        xls  = as.data.frame(readxl::read_excel(path))
      ),
      error = function(e) {

        showModal(
          modalDialog(
            title = "❌ Could not read file",
            paste("Error:", conditionMessage(e)),
            easyClose = TRUE
          )
        )

        return(NULL)

      }
    )

    if (is.null(df))
      req(FALSE)

  }
     
      if (identical(names(df)[1], "")) {
        df <- df[, -1, drop = FALSE]
      }

      names(df) <- trimws(names(df))

      # ---- STRUCTURE CHECKS ----
      required_cols <- c("File Name", "Condition", "Temperature", "Replicate")
      missing_cols <- setdiff(required_cols, names(df))

      if (length(missing_cols) > 0) {
        showModal(
          modalDialog(
            title = "❌ Incorrect data structure",
            paste("Missing columns:", paste(missing_cols, collapse = ", ")),
            easyClose = TRUE
          )
        )
        req(FALSE)
      }

      # ---- type coercion ----
      df$Temperature <- suppressWarnings(as.numeric(df$Temperature))
      df$Replicate <- as.character(
        df$Replicate
      )

      if (any(trimws(df$Replicate) == "")) {

        showModal(
          modalDialog(
            "Replicate contains empty values",
            easyClose = TRUE
          )
        )

        return(NULL)

      }

      if (any(is.na(df$Replicate))) {
        showModal(modalDialog("Replicate must be integer", easyClose = TRUE))
        req(NULL)
      }

      if (any(is.na(df[, required_cols]))) {
        showModal(modalDialog("Missing values in required columns", easyClose = TRUE))
        req(FALSE)
      }

      if (any(duplicated(df$`File Name`))) {
        showModal(modalDialog("Duplicate File Names detected", easyClose = TRUE))
        req(NULL)
      }

      # ---- organism parsing ----


      df
    })

    # -----------------------------
    # Org_data reactive (clean)
    # -----------------------------
Org_data <- reactive({

  cat("\n====================\n")
  cat("ENTERING ORG_DATA()\n")
  cat("====================\n")

  message("[ORG_DATA] entered")

  req(data())

  message("[ORG_DATA] entered")

  req(data())

  message("[ORG_DATA] after data()")

 orgs <- unique(
  unlist(
    selected_organisms()
  )
)

orgs <- setdiff(
  orgs,
  "Other"
)

  message("[ORG_DATA] organisms found:")
  print(orgs)

  if (length(orgs) == 0) {
    return(NULL)
  }


    org_files <- list(
      list(
        names = c("Sulfolobus acidocaldarius"),
        file  = "Data/SACIuniprotkb_taxonomy_id_2285_2025_09_30.tsv"
      ),
      list(
        names = c("Parageobacillus thermoglucosidasius"),
        file  = "Data/uniprotkb_taxonomy_id_1426_2026_01_20.tsv"
      ),
      list(
        names = c("Caldimonas thermodepolymerans"),
        file  = "Data/uniprotkb_proteome_UP000239406_2026_01_20 (1).tsv"
      ),
      list(
        names = c("Haloferax volcanii"),
        file  = "Data/uniprotkb_taxonomy_id_2246_2026_01_20.tsv"
      ),
      list(
        names = c("Thermosynechococcus elongatus"),
        file  = "Data/uniprotkb_taxonomy_id_197221_2026_01_21.tsv"
      ),
      list(
        names = c("Thermus thermophilus"),
        file  = "Data/uniprotkb_proteome_UP000000532_2026_01_20.tsv"
      ),
      list(
        names = c("E.coli", "Escherichia coli", "E coli"),
        file  = "Data/uniprotkb_E_coli_AND_model_organism_833_2026_01_20.tsv"
      ),
      list(
        names = c("H.sapiens", "Homo sapiens", "human"),
        file  = "Data/uniprotkb_taxonomy_id_9606_2026_01_20.tsv"
      )
    )
    data_list <- lapply(orgs, function(org) {
message("[ORG_DATA] processing: ", org)
      # normalize input
      org_clean <- tolower(trimws(org))

      # find matching entry
      match_entry <- NULL

      for (entry in org_files) {
        aliases <- tolower(entry$names)

        if (any(org_clean == aliases)) {
          match_entry <- entry
          break
        }
      }

      if (!is.null(match_entry)) {

cache_key <- match_entry$file

if (!is.null(org_file_cache[[cache_key]])) {

  message("[ORG_DATA] CACHE HIT: ", cache_key)

  df_org <- org_file_cache[[cache_key]]

} else {

  message("[ORG_DATA] reading file: ", cache_key)

  df_org <- read.delim(
    cache_key,
    header = TRUE,
    stringsAsFactors = FALSE
  )

  org_file_cache[[cache_key]] <- df_org

  message("[ORG_DATA] cached file: ", cache_key)
}

        df_org$Organism <- org

      } else {
        df_org <- NULL
      }

      df_org
    })

  data_list <- Filter(
  Negate(is.null),
  data_list
)

# --------------------------------------------------
# Add custom uploaded organism databases
# --------------------------------------------------

custom_list <- customOrgData()

if (!is.null(custom_list)) {

  custom_dfs <- unname(custom_list)

  data_list <- c(
    data_list,
    custom_dfs
  )

}

if (length(data_list) == 0) {
  return(NULL)
}

final_df <- do.call(
  rbind,
  data_list
)
cat("\nORG_DATA DIM:\n")
print(dim(final_df))
message(
  "[ORG_DATA] final rows: ",
  nrow(final_df)
)

final_df
    })

    return(list(data = data, Org_data = Org_data))
  })
}


#module to bind datasets and split into conditions
processedDataModuleServer <- function(
  id,
  proteinInput,
  designInput,
  splitMode,
  organismMapping,
  customProteinMapping
) {
  moduleServer(id, function(input, output, session) {

    # ✅ GLOBAL SAFE NAME CLEANER
    clean_names_safe <- function(df) {
      nm <- names(df)
      bad <- is.na(nm) | nm == ""
      nm[bad] <- paste0("V", seq_len(sum(bad)))
      nm <- trimws(nm)
      names(df) <- nm
      df
    }
organism_aliases <- list(

  "H.sapiens" = c(
    "H.sapiens",
    "Homo sapiens",
    "human"
  ),

  "E.coli" = c(
    "E.coli",
    "Escherichia coli",
    "E coli"
  ),

  "Sulfolobus acidocaldarius" = c(
    "Sulfolobus acidocaldarius"
  ),

  "Parageobacillus thermoglucosidasius" = c(
    "Parageobacillus thermoglucosidasius"
  ),

  "Caldimonas thermodepolymerans" = c(
    "Caldimonas thermodepolymerans"
  ),

  "Haloferax volcanii" = c(
    "Haloferax volcanii"
  ),

  "Thermosynechococcus elongatus" = c(
    "Thermosynechococcus elongatus"
  ),

  "Thermus thermophilus" = c(
    "Thermus thermophilus"
  )

)
    data <- reactive({



      req(proteinInput(), designInput(), splitMode())
progress <- shiny::Progress$new()

on.exit(
  progress$close()
)

progress$set(
  value = 0,
  message = "Preparing datasets",
  detail = "Loading inputs"
)

      proteins <- proteinInput()
      design   <- designInput()

  progress$set(
  value = 0.10,
  detail = "Loading protein and design tables"
)
  cat("\n====================\n")
  cat("CUSTOM PROTEIN MAPPING RECEIVED\n")
  cat("====================\n")

  print(
    customProteinMapping()
  )
      # ✅ ALWAYS CLEAN INPUTS
      proteins <- clean_names_safe(proteins)
      design   <- clean_names_safe(design)
progress$set(
  value = 0.25,
  detail = "Validating dataset structure"
)
      # ----------------------------------------------------
# Ensure Organism column exists for single-organism datasets
# ----------------------------------------------------

if (!"Organism" %in% names(proteins)) {

  mapping <- organismMapping()

  all_orgs <- unique(unlist(mapping))

  if (length(all_orgs) == 1) {

    proteins$Organism <- all_orgs[1]

  } else {

    validate(
      need(
        FALSE,
        paste(
          "Multiple organisms detected.",
          "Protein dataset must contain an Organism column."
        )
      )
    )

  }

}

# ----------------------------------------------------
# Validate selected organisms are present
# in protein dataset
# ----------------------------------------------------

if ("Organism" %in% names(proteins)) {

  mapping <- organismMapping()

  selected_orgs <- unique(
    unlist(mapping)
  )

  selected_orgs <- setdiff(
    selected_orgs,
    "Other"
  )

  protein_orgs <- unique(
    as.character(
      proteins$Organism
    )
  )

  missing_orgs <- c()

for(org in selected_orgs) {

  aliases <- organism_aliases[[org]]

  if(is.null(aliases)) {
    aliases <- org
  }

  found <- FALSE

  for(alias in aliases) {

    matches <- grepl(
      pattern = tolower(alias),
      x = tolower(protein_orgs),
      fixed = TRUE
    )

    if(any(matches)) {

      found <- TRUE
      break

    }

  }

  if(!found) {

    missing_orgs <- c(
      missing_orgs,
      org
    )

  }

}


  if(length(missing_orgs) > 0) {

    showModal(

      modalDialog(

        title = "❌ Organism validation failed",

        paste0(
          "The following organism(s) were selected in the Condition to Organism Mapping section but were not found in the Organism column of the protein dataset:\n\n",
          paste(
            missing_orgs,
            collapse = "\n"
          )
        ),

        easyClose = TRUE

      )

    )

    req(FALSE)

  }

}

progress$set(
  value = 0.40,
  detail = "Validating organism mappings"
)


      # ---- 1. Build lookup ----
lookup <- design %>%
  mutate(
    new_name = paste0(
      Condition,
      "_",
      Temperature,
      "_REP",
      Replicate
    )
  ) %>%
  select(
    `File Name`,
    new_name,
    Condition
  )


      # ---- 2. Rename columns in protein table ----
      rename_columns_by_lookup <- function(prot, lookup) {
        cn <- colnames(prot)

        for (i in seq_len(nrow(lookup))) {
          old <- lookup$`File Name`[i]
          new <- lookup$new_name[i]

          match_idx <- grep(old, cn, fixed = TRUE)

          if (length(match_idx) > 0) {
            cn[match_idx] <- new
          }
        }

        colnames(prot) <- cn
        prot
      }

      proteins <- rename_columns_by_lookup(proteins, lookup)

progress$set(
  value = 0.55,
  detail = "Matching protein data to experimental design"
)
      # ✅ CLEAN AGAIN AFTER RENAMING (VERY IMPORTANT)
      proteins <- clean_names_safe(proteins)

      # ---- 3. Filter multiple accessions ----
      proteins <- proteins %>%
        filter(!stringr::str_detect(Accession, ";"))

      # ---- 4. Keep only relevant columns ----
      cols_to_keep <- base::intersect(
        c("Accession", "Organism", lookup$new_name),
        colnames(proteins)
      )
      proteins <- proteins[, cols_to_keep, drop = FALSE]

      # ✅ FINAL CLEAN BEFORE PIVOT (CRITICAL)
      proteins <- clean_names_safe(proteins)
progress$set(
  value = 0.70,
  detail = "Generating analysis datasets"
)
      # ============================================================
      #  MODE 1: CONDITION ONLY
      # ============================================================
if (splitMode() == "condition") {

  split_list <- list()

  conditions <- unique(lookup$Condition)

  for (cond in conditions) {

    cond_cols <- lookup %>%
      dplyr::filter(Condition == cond) %>%
      dplyr::pull(new_name)

    keep_cols <- intersect(
      c("Accession", "Organism", cond_cols),
      names(proteins)
    )

split_list[[cond]] <- proteins[
  ,
  keep_cols,
  drop = FALSE
] %>%
  dplyr::select(
    -Organism
  )

  }

}


      # ============================================================
      # MODE 2: CONDITION + ORGANISM
      # ============================================================
else if (splitMode() == "condition_organism") {

  t_pivot <- Sys.time()

  split_list <- list()

  mapping <- organismMapping()

  for (cond in names(mapping)) {

    orgs <- mapping[[cond]]

    cond_cols <- lookup %>%
      dplyr::filter(
        Condition == cond
      ) %>%
      dplyr::pull(new_name)

for (org in orgs) {

  # ----------------------------------
  # Standard organism
  # ----------------------------------

  if (org != "Other") {

    aliases <- organism_aliases[[org]]

    if (is.null(aliases)) {
      aliases <- org
    }

  } else {

    # ----------------------------------
    # Custom organism
    # ----------------------------------

    aliases <- customProteinMapping()[[cond]]

  }

  cat("\n---------------------\n")
  cat("CONDITION:\n")
  print(cond)

  cat("\nORG:\n")
  print(org)

  cat("\nALIASES USED:\n")
  print(aliases)

  keep_cols <- intersect(
    c(
      "Accession",
      "Organism",
      cond_cols
    ),
    names(proteins)
  )


dataset_name <- if (

  org == "Other"

) {

  customProteinMapping()[[cond]]

} else {

  org

}

split_list[[paste0(
  cond,
  "_",
  dataset_name
)]] <-
  proteins %>%
  dplyr::filter(
    tolower(trimws(Organism)) %in%
      tolower(trimws(aliases))
  ) %>%
  dplyr::select(
    all_of(keep_cols)
  ) %>%
  dplyr::select(
    -Organism
  )


    }

  }



}



if (
  splitMode() %in%
  c(
    "condition",
    "condition_organism"
  )
) {

cat("\n=====================\n")
cat("PROCESSED OUTPUT\n")
cat("=====================\n")

for (nm in names(split_list)) {

  cat("\nDATASET:", nm, "\n")

  print(names(split_list[[nm]]))

  cat("\nHEAD:\n")
  print(head(split_list[[nm]][, 1:min(6, ncol(split_list[[nm]]))]))

}

split_list <- lapply(
  split_list,
  function(df) {

    nm <- names(df)

    nm[-1] <- sub(
      "^[^_]+_",
      "",
      nm[-1]
    )

    names(df) <- nm

    # -----------------------
    # Convert string NaN -> NA
    # -----------------------

    value_cols <- names(df)[-1]

    df[value_cols] <- lapply(
      df[value_cols],
      function(x) {

        x <- as.character(x)

        x[trimws(x) %in% c(
          "NaN",
          "nan",
          "NA",
          ""
        )] <- NA

        suppressWarnings(
          as.numeric(x)
        )

      }
    )

    df
  }
)
progress$set(
  value = 1,
  detail = "Finished"
)
  return(split_list)

}


t_wide <- Sys.time()

      result <- lapply(split_list, function(df) {

        df <- clean_names_safe(df)  # ✅ IMPORTANT BEFORE PIVOT


        out <- df %>%
          select(-Condition) %>%
          pivot_wider(
            names_from = rest,
            values_from = abundance
          )

        colnames(out)[1] <- "Accession"

        # FIX numeric conversion safely
        value_cols <- setdiff(colnames(out), c("Accession", "Organism"))

        out[, value_cols] <- lapply(out[, value_cols], function(x) {
          x <- trimws(x)
          x[x == "NaN"] <- NA
          x[x == ""] <- NA
          as.numeric(x)
        })


        out
      })


removeNotification(progress_id)
return(result)
    })

dataset_organisms <- reactive({

  req(organismMapping())

  mapping <- organismMapping()

  # -------------------------
  # CONDITION ONLY
  # -------------------------
if (splitMode() == "condition") {

  out <- list()

  conditions <- names(mapping)

  all_orgs <- unique(
    setdiff(
      unlist(mapping),
      "Other"
    )
  )

  if (length(all_orgs) == 1) {

    for (cond in conditions) {

      out[[cond]] <- all_orgs[1]

    }

    cat("\nDATASET_ORGANISMS OUTPUT:\n")
    print(out)

    return(out)

  }

  return(NULL)

}

  # -------------------------
  # CONDITION + ORGANISM
  # -------------------------
  if (splitMode() == "condition_organism") {

    out <- list()

    for (cond in names(mapping)) {

      orgs <- mapping[[cond]]

      for (org in orgs) {

        out[[paste0(cond, "_", org)]] <- org

      }

    }

    return(out)

  }

  NULL

})
return(
  list(
    data = data,
    dataset_organisms = dataset_organisms
  )
)
  })
}




###############################################################################################
#QC
qcPanelModuleUI <- function(id, label) {

  ns <- NS(id)

  tagList(

    # ======================
    # Missing values section
    # ======================

    fluidRow(

      box(
        title = "Handle missing values",
        width = 12,
        status = "primary",
        solidHeader = TRUE,

fluidRow(

  uiOutput(
    ns("missing_info_card")
  )

),
fluidRow(

  column(

    width = 6,

    radioButtons(
      inputId = ns("choices"),
      label = "Would you like to:",
      choices = c(
        "Remove all rows with any missing values" =
          "Remove all",
        "Impute all missing values" =
          "Impute",
        "Remove all rows containing only missing values and Impute the rest" =
          "Remove and Impute"
      ),
      selected = "Impute",
      inline = FALSE
    )

  ),

  column(

    width = 6,

    conditionalPanel(
      condition = sprintf(
        "input['%s'] == 'Impute' || input['%s'] == 'Remove and Impute'",
        ns("choices"),
        ns("choices")
      ),

      radioButtons(
        inputId = ns("impute_method"),
        label = "Choose imputation method:",
        choices = c(
          "knn",
          "MLE",
          "min",
          "All to 0"
        ),
        selected = "All to 0"
      )

    )

  )

),

        actionButton(
          inputId = ns("apply"),
          label = "Apply transformation"
        )

      )

    ),

    # ======================
    # 0-value section
    # ======================

    fluidRow(

      box(
        title = "Handle zero values",
        width = 12,
        status = "primary",
        solidHeader = TRUE,

fluidRow(

  uiOutput(
    ns("zero_info_card")
  )

),


        uiOutput(
          ns("filter_0_ui")
        )

      )

    ),

    # ======================
    # Heatmap
    # ======================

    fluidRow(

      box(
        title = "Distribution of zero values",
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        radioButtons(
          inputId = ns("pattern_mode"),
          label = "Display mode:",
          choices = c(
            "Clustered" = "normal",
            "Per replicate" = "clustered"
          ),
          selected = "normal",
          inline = TRUE
        ),

        radioButtons(
          ns("row_sort_mode"),
          "Order proteins by:",
          choices = c(
            "Clustering" = "cluster",
            "Average intensity" = "intensity"
          ),
          selected = "intensity",
          inline = TRUE
        ),

        plotOutput(
          ns("missing_values_pattern"),
          height = 700
        ) %>% withSpinner(
          color = "#0EA5A5"
        ),

        tags$div(
          style = "height:8px;"
        ),

        downloadButton(
          ns("download_missing_values_pattern"),
          "Download heatmap"
        )

      )

    ),

    # ======================
    # Density plot
    # ======================

    fluidRow(

      box(
        title = "Density plots of average log2 intensity",
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        plotOutput(
          ns("missing_values_density"),
          height = 350
        ) %>% withSpinner(
color = "#0EA5A5"
),

        downloadButton(
          ns("download_missing_values_density"),
          "Download density plot"
        )

      )

    ),

    # ======================
    # Bottom plots
    # ======================

    fluidRow(

      box(
        title = "Average number of zero values per temperature",
        width = 6,
        status = "primary",
        solidHeader = TRUE,

        plotOutput(
          ns("zero_values_per_temp"),
          height = 300
        ) %>% withSpinner(
color = "#0EA5A5"
),

        downloadButton(
          ns("download_zero_values_per_temp"),
          "Download plot"
        )

      ),

      box(
        title = "Total abundance plot",
        width = 6,
        status = "primary",
        solidHeader = TRUE,

        plotOutput(
          ns("total_abundance_plot"),
          height = 300
        ) %>% withSpinner(
color = "#0EA5A5"
),

        downloadButton(
          ns("download_total_abundance_plot"),
          "Download plot"
        )

      )

    )

  )

}



#server


# ---- Missing data handling (apply-on-click) ----

missingDataModuleServer <- function(id, dataReactive, choicesInput, methodInput) {
  moduleServer(id, function(input, output, session) {

  observe({
    req(dataReactive())

    df <- dataReactive()

  })

    filteredData <- eventReactive(input$apply, {
      req(dataReactive())
      df <- dataReactive()


      choice <- choicesInput()
      method <- methodInput()
withProgress(
  message = "Applying transformation",
  value = 0,
  {
      # -------------------------
      # REMOVE ALL ROWS WITH ANY NA (consider columns 2:n
      # -------------------------
incProgress(
  0.10,
  detail = "Checking selected transformation"
)
      if ("Remove all" %in% choice) {

        df <- df[
          complete.cases(df[, 2:ncol(df), drop = FALSE]),
          ,
          drop = FALSE
        ]
incProgress(
  1,
  detail = "Finished"
)
      # -------------------------
      # IMPUTE ONLY (NO ROW REMOVAL)
      # -------------------------
      } else if ("Impute" %in% choice) {
        req(method)
incProgress(
  0.20,
  detail = paste("Preparing", method, "imputation")
)

        if (method == "All to 0") {
          df[, 2:ncol(df)] <- lapply(df[, 2:ncol(df), drop = FALSE], function(col) {
            col[is.na(col)] <- 0
            col
          })
incProgress(
  1,
  detail = "Finished"
)
        } else if (method == "knn") {
          # ---- Fast & robust KNN using Bioconductor impute::impute.knn ----
          shiny::withProgress(message = "KNN (Bioconductor) — impute.knn", value = 0, {
            cat("\n[impute.knn] start: nrow=", nrow(df), " ncol=", ncol(df), "\n", sep = "")
            incProgress(0.1)

            data_cols <- which(sapply(df, is.numeric))
            mat <- as.matrix(df[, data_cols, drop = FALSE])
            if (!is.numeric(mat)) mode(mat) <- "numeric"

            # Remove rows that are only NA across data columns
            only_na <- rowSums(!is.na(mat)) == 0
            if (any(only_na)) {
              cat("[impute.knn] removing rows with only-NA across data cols: ", sum(only_na), "\n", sep = "")
              df  <- df[!only_na, , drop = FALSE]
              mat <- as.matrix(df[, data_cols, drop = FALSE])
              if (!is.numeric(mat)) mode(mat) <- "numeric"
            }
            incProgress(0.25)

            donors_n <- sum(rowSums(!is.na(mat)) > 0)
            k_eff    <- max(1, min(10, floor(donors_n / 3)))
            cat("[impute.knn] donors=", donors_n, " k=", k_eff, " maxp=3000\n", sep = "")
            incProgress(0.5)

            imp <- impute::impute.knn(
              data   = mat,
              k      = k_eff,
              rowmax = 0.95,
              colmax = 0.95,
              maxp   = 3000
            )
            incProgress(0.9)

            df[, data_cols] <- imp$data
            cat("[impute.knn] done\n")
          })
incProgress(
  1,
  detail = "Finished"
)
        } else if (method == "MLE") {

          numeric_cols <- names(df)[sapply(df, is.numeric)]
          df_numeric <- df[, numeric_cols, drop = FALSE]

          # If no NA, return as-is
          if (!any(is.na(df_numeric))) {
            return(df)
          }

          # ---- SIZE GUARD: disable Amelia for large data ----
          row_threshold <- 3000L  # adjust if needed
          if (nrow(df_numeric) > row_threshold) {
            # NOTE: validate/need doesn't show here; use a notification + stop
            shiny::showNotification(
              "MLE (Amelia) is disabled for large data. Please select other method",
              type = "warning", duration = 6
            )
            req(FALSE)  # halt this eventReactive evaluation
          }

          rownames(df_numeric) <- NULL
          imputed_try <- tryCatch(
            Amelia::amelia(df_numeric, m = 1, quiet = TRUE),
            error = function(e) e
          )

          if (inherits(imputed_try, "error")) {
            stop(imputed_try$message)
          }

          imputed <- imputed_try$imputations[[1]]
          df[, numeric_cols] <- imputed[, numeric_cols]
incProgress(
  0.90,
  detail = "Applying imputed values"
)
          df[, numeric_cols] <- lapply(df[, numeric_cols], function(col) {
            col[col < 0] <- 0
            col
          })
incProgress(
  1,
  detail = "Finished"
)
        } else if (method == "min") {
         
          data_cols <- which(sapply(df, is.numeric))

          df[, data_cols] <- lapply(df[, data_cols, drop = FALSE], function(col) {
            col[is.na(col)] <- min(col, na.rm = TRUE)
            col
          })
        }

      # -------------------------
      # REMOVE ONLY-NA ROWS & IMPUTE REST
      # -------------------------
      } else if ("Remove and Impute" %in% choice) {
        req(method)

       
        data_cols <- which(sapply(df, is.numeric))

        # Step 1: remove rows that contain ONLY NA in columns 2:n
        only_na <- rowSums(!is.na(df[, data_cols, drop = FALSE])) == 0
        df <- df[!only_na, , drop = FALSE]

        # Step 2: impute remaining NAs
        if (method == "All to 0") {
          df[, data_cols] <- lapply(df[, data_cols, drop = FALSE], function(col) {
            col[is.na(col)] <- 0
            col
          })

        } else if (method == "knn") {
          shiny::withProgress(message = "KNN (Bioconductor) — impute.knn", value = 0, {
            cat("\n[impute.knn] start: nrow=", nrow(df), " ncol=", ncol(df), "\n", sep = "")
            incProgress(0.1)

            mat <- as.matrix(df[, data_cols, drop = FALSE])
            if (!is.numeric(mat)) mode(mat) <- "numeric"

            # (Rows with only NA already removed above)
            donors_n <- sum(rowSums(!is.na(mat)) > 0)
            k_eff    <- max(1, min(10, floor(donors_n / 3)))
            cat("[impute.knn] donors=", donors_n, " k=", k_eff, " maxp=3000\n", sep = "")
            incProgress(0.5)

            imp <- impute::impute.knn(
              data   = mat,
              k      = k_eff,
              rowmax = 0.95,
              colmax = 0.95,
              maxp   = 3000
            )
            incProgress(0.9)

            df[, data_cols] <- imp$data
            cat("[impute.knn] done\n")
          })

        } else if (method == "MLE") {

          numeric_cols <- names(df)[sapply(df, is.numeric)]
          df_numeric <- df[, numeric_cols, drop = FALSE]

          # If no NA in numeric part, nothing to do
          if (!any(is.na(df_numeric))) {
            return(df)
          }

          # ---- SIZE GUARD: disable Amelia for large data ----
          row_threshold <- 3000L  # adjust if needed
          if (nrow(df_numeric) > row_threshold) {
            shiny::showNotification(
              "MLE (Amelia) is disabled for large data. Please select other method",
              type = "warning", duration = 6
            )
            req(FALSE)  # halt this eventReactive evaluation
          }

          rownames(df_numeric) <- NULL
          imputed_try <- tryCatch(
            Amelia::amelia(df_numeric, m = 1, quiet = TRUE),
            error = function(e) e
          )

          if (!inherits(imputed_try, "error")) {
            imputed <- imputed_try$imputations[[1]]
            df[, numeric_cols] <- imputed[, numeric_cols]

            df[, numeric_cols] <- lapply(df[, numeric_cols], function(col) {
              col[col < 0] <- 0
              col
            })
          }

        } else if (method == "min") {
          df[-1] <- lapply(df[-1], function(col) {
            col[is.na(col)] <- min(col, na.rm = TRUE)
            col
          })
        }
      }
}
)
      df
    }, ignoreInit = TRUE)

    return(filteredData)
  })
}


#missing percentage
missingStatsModule <- function(dataReactive) {
  missing_percentage <- reactive({
    df <- dataReactive()[, -1, drop = FALSE]
    total_values <- prod(dim(df))
    missing_values <- sum(is.na(df) | is.infinite(as.matrix(df)) | is.nan(as.matrix(df)))
    (missing_values / total_values) * 100
  })



  return(list(
    percentage = missing_percentage
  ))
}
zeroStatsModule <- function(
  id,
  filtered_data,
  already_log2
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive: Filtered data based on 0-value threshold
    filtered_0_data <- reactive({
      
      df <- filtered_data()


numeric_df <- df[, sapply(df, is.numeric), drop = FALSE]

all_zero <- apply(
  numeric_df,
  1,
  function(x) all(x == 0, na.rm = TRUE)
)

print(sum(all_zero))

cat("Rows with at least one positive value:\n")
print(sum(rowSums(numeric_df > 0, na.rm = TRUE) > 0))
req(
  input$mean_intensity_filter
)


threshold <- input$filter_0

if (
  !is.null(threshold) &&
  threshold > 0
) {
        df <- df %>%
          filter(rowSums(across(-1, ~ . == 0)) <= threshold)
      }
            numeric_df <- df[, sapply(df, is.numeric), drop = FALSE]

      row_means <- apply(
        numeric_df,
        1,
        function(x) {

          x <- x[x > 0]

          if (length(x) == 0) {
            return(0)
          }

          if (already_log2()) {

            mean(x)

          } else {

            mean(log2(x + 1))

          }


        }
      )


        range_selected <- input$mean_intensity_filter

        if (
          !is.null(range_selected)
        ) {

          df <- df[
            row_means >= range_selected[1] &
              row_means <= range_selected[2],
            ,
            drop = FALSE
          ]

        }
        cat("\nRows AFTER filtering:\n")
print(nrow(df))

numeric_df2 <- df[, sapply(df, is.numeric), drop = FALSE]

all_zero2 <- apply(
  numeric_df2,
  1,
  function(x) all(x == 0, na.rm = TRUE)
)

cat("All-zero rows after filtering:\n")
print(sum(all_zero2))
      df
    })

    # Reactive: Percentage of 0 values
zero_percentage <- reactive({

  df <- filtered_0_data()

  if (nrow(df) == 0) {
    return(NA_real_)
  }

  df <- df[, -1, drop = FALSE]

  total_values <- prod(dim(df))

  if (total_values == 0) {
    return(NA_real_)
  }

  zero_values <- sum(df == 0, na.rm = TRUE)

  (zero_values / total_values) * 100

})

    # Reactive, icon based on percentage
    zero_icon <- reactive({
      pct <- zero_percentage()
      if (pct < 5) {
        icon("check-circle")   # positive icon
      } else {
        icon("exclamation-triangle")  # warning icon
      }
    })

    return(list(
      filtered = filtered_0_data,
      percentage = zero_percentage,
      icon = zero_icon
    ))
  })
}


qcPanelModuleServer <- function(
  id,
  filtered_data,
  missing_stats,
  zero_stats,
  already_log2
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
last_heatmap <- reactiveVal(NULL)
last_tick_labels <- reactiveVal(NULL)
#progress bar
withProgress(
  message = paste("Initializing QC module:", id),
  value = 0,
  {

    incProgress(
      0.2,
      detail = "Preparing filtered dataset"
    )

    # Helper: safely obtain filtered data (NULL if not available yet)
    safe_df <- reactive({
      tryCatch({
        df <- filtered_data()
        if (is.null(df) || !is.data.frame(df) || ncol(df) < 2) return(NULL)
        df
      }, error = function(e) NULL)
    })
    incProgress(
      1,
      detail = "QC module ready"
    )

  })
    # Helper: message string
    msg <- "Please choose how to Handle the missing data"

    #filter summary output
    output$filter_summary <- renderText({

      df_original <- safe_df()

      req(df_original)

      n_original <- nrow(df_original)

      n_remaining <- nrow(
        zero_stats$filtered()
      )

      threshold_0 <- input$filter_0

      intensity_range <- input$mean_intensity_filter

      paste0(

        "Proteins retained: ",
        format(n_remaining, big.mark = ","),
        " / ",
        format(n_original, big.mark = ","),

        "\n",

        "Proteins removed: ",
        format(
          n_original - n_remaining,
          big.mark = ","
        )

      )

    })
    # -------------------------------
    # Info box for missing values
    # -------------------------------
output$missing_info_card <- renderUI({

  div(
    class = "qc-info-card",

    div(
      class = "qc-info-title",
      "Initial Percentage of Missing Values"
    ),

    div(
      class = "qc-info-value",
      paste0(
        round(
          missing_stats$percentage(),
          2
        ),
        "%"
      )
    ),

    div(
      class = "qc-info-subtitle",
      "NA, NaN and Inf values"
    )
  )

})

    # -------------------------------
    # UI for filtering 0 values
    # -------------------------------

    output$filter_0_ui <- renderUI({
      df <- safe_df()
      if (is.null(df)) {
        return(helpText("Apply a transformation to enable zero-value filtering."))
      }



      max_val <- ncol(df) - 1

numeric_df <- df[, sapply(df, is.numeric), drop = FALSE]

row_means <- apply(
  numeric_df,
  1,
  function(x) {

    x <- x[x > 0]

    if (length(x) == 0) {
      return(0)
    }

    if (already_log2()) {

      mean(x)

    } else {

      mean(log2(x + 1))

    }

  }
)

min_intensity <- floor(
  min(row_means, na.rm = TRUE)
)

max_intensity <- ceiling(
  max(row_means, na.rm = TRUE)
)



has_zeros <- any(df == 0, na.rm = TRUE)

controls <- list()

# ---------------------------------
# Zero-value filter only if needed
# ---------------------------------

if (has_zeros) {

  controls <- c(
    controls,

    list(

      sliderInput(
        inputId = ns("filter_0"),
        label = "Filter out proteins with number of zero values >= to:",
        min = 0,
        max = max(0, max_val),
        value = max(0, max_val),
        step = 1
      )

    )

  )

}

# ---------------------------------
# Intensity filter ALWAYS available
# ---------------------------------

controls <- c(
  controls,

  list(

    sliderInput(
      inputId = ns("mean_intensity_filter"),
      label = "Average log2 intensity range:",
      min = min_intensity,
      max = max_intensity,
      value = c(
        min_intensity,
        max_intensity
      ),
      step = 0.1
    ),

    verbatimTextOutput(
      ns("filter_summary")
    )

  )

)

do.call(tagList, controls)

    })


    # -------------------------------
    # Info box for 0 values
    # -------------------------------
output$zero_info_card <- renderUI({

  pct <- zero_stats$percentage()

  value_text <- if (
    is.na(pct)
  ) {
    "No proteins remain after filtering"
  } else {
    paste0(round(pct, 2), "%")
  }

  div(
    class = "qc-info-card",

    div(
      class = "qc-info-title",
      "Percentage of Zero Values"
    ),

    div(
      class = "qc-info-value",
      value_text
    ),

    div(
      class = "qc-info-subtitle",
      "Current filtered dataset"
    )
  )

})

    # -------------------------------
    # Plot: Average number of 0 values per temperature
    # -------------------------------

    output$zero_values_per_temp <- renderPlot({
      dfm <- safe_df()
      validate(
        need(!is.null(dfm), msg),
        need(any(dfm == 0, na.rm = TRUE), "No zero values in the data to load this plot.")
      )

      df0 <- zero_stats$filtered()

      zero_df0 <- summarize_zeros(df0)
      zero_dfm <- summarize_zeros(dfm)

      max_mean <- zero_dfm %>%
        dplyr::group_by(temperature) %>%
        dplyr::summarise(Mean = mean(Zero_Count), SD = sd(Zero_Count), .groups = "drop") %>%
        dplyr::summarise(Max_Mean = max(Mean + SD), .groups = "drop") %>%
        dplyr::pull(Max_Mean)

      ylim_limit <- if (is.finite(max_mean)) max_mean + 5 else NA

      zero_df0 %>%
        dplyr::group_by(temperature) %>%
        dplyr::summarise(Mean = mean(Zero_Count), SD = sd(Zero_Count), .groups = "drop") %>%
        ggplot2::ggplot(ggplot2::aes(temperature, Mean)) +
        ggplot2::geom_col() +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = Mean - SD, ymax = Mean + SD), width = .2
        ) +
        ggplot2::labs(y = "Average number of zero values", x = "Temperature") +
        ggplot2::theme_minimal() +
        ggplot2::coord_cartesian(ylim = c(0, ylim_limit))
    })
output$download_zero_values_per_temp <- downloadHandler(

  filename = function() {
    "zero_values_per_temperature.png"
  },

  content = function(file) {

    png(
      file,
      width = 1800,
      height = 1200,
      res = 200
    )

    dfm <- safe_df()
    df0 <- zero_stats$filtered()

    zero_df0 <- summarize_zeros(df0)
    zero_dfm <- summarize_zeros(dfm)

    max_mean <- zero_dfm %>%
      dplyr::group_by(temperature) %>%
      dplyr::summarise(
        Mean = mean(Zero_Count),
        SD = sd(Zero_Count),
        .groups = "drop"
      ) %>%
      dplyr::summarise(
        Max_Mean = max(Mean + SD),
        .groups = "drop"
      ) %>%
      dplyr::pull(Max_Mean)

    ylim_limit <- if (is.finite(max_mean)) max_mean + 5 else NA

    p <- zero_df0 %>%
      dplyr::group_by(temperature) %>%
      dplyr::summarise(
        Mean = mean(Zero_Count),
        SD = sd(Zero_Count),
        .groups = "drop"
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(temperature, Mean)
      ) +
      ggplot2::geom_col() +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = Mean - SD,
          ymax = Mean + SD
        ),
        width = .2
      ) +
      ggplot2::coord_cartesian(
        ylim = c(0, ylim_limit)
      ) +
      ggplot2::theme_minimal()

    print(p)

    dev.off()

  }

)

    # -------------------------------
    # Plot: Density of average log2 intensity
    # -------------------------------

  output$missing_values_density <- renderPlot({

  withProgress(
    message = "Generating density plot",
    value = 0,
    {
      df <- safe_df()
      validate(
        need(!is.null(df), msg),
        need(any(df == 0, na.rm = TRUE), "No zero values in the data to load this plot.")
      )

      df0 <- zero_stats$filtered()
      validate(
  need(
    any(df0 == 0, na.rm = TRUE),
    "No proteins with zero values remain after filtering."
  )
)
if (already_log2()) {

  dflog2 <- df0

} else {

  dflog2 <- df0 %>%
    dplyr::mutate(
      dplyr::across(
        where(is.numeric),
        log2
      )
    )

  dflog2[dflog2 == -Inf] <- 0

}
      incProgress(
        0.8,
        detail = "Calculating density"
      )

p <- plot_detect_custom_df(dflog2)

incProgress(
  1,
  detail = "Done"
)

p

    }
  )

})
output$download_missing_values_density <- downloadHandler(

  filename = function() {
    "density_plot.png"
  },

  content = function(file) {

    png(
      file,
      width = 1800,
      height = 1200,
      res = 200
    )

    df0 <- zero_stats$filtered()

    dflog2 <- df0 %>%
      dplyr::mutate(
        dplyr::across(where(is.numeric), log2)
      )

    dflog2[dflog2 == -Inf] <- 0

    print(
      plot_detect_custom_df(dflog2)
    )

    dev.off()

  }

)
    # -------------------------------
    # Plot: Missing values pattern heatmap
    # -------------------------------
    heatmap_data <- reactive({

  req(
    input$pattern_mode,
    input$row_sort_mode
  )

  df <- safe_df()
  validate(
    need(!is.null(df), msg)
  )

  df0 <- zero_stats$filtered()

  list(
    df = df,
    df0 = df0
  )

})
output$missing_values_pattern <- renderPlot({

  withProgress(
    message = "Building heatmap",
    value = 0,
    {

req(
  input$pattern_mode,
  input$row_sort_mode
)
      df <- safe_df()

      validate(
        need(!is.null(df), msg)
      )

      df0 <- zero_stats$filtered()

      incProgress(
        0.25,
        detail = "Preparing matrix"
      )

      # keep numeric columns only
      df_numeric <- df0[, sapply(df0, is.numeric), drop = FALSE]

      validate(
        need(ncol(df_numeric) > 0, "No numeric columns to display."),
        need(any(df_numeric == 0, na.rm = TRUE), "No zero values in the data to display.")
      )

      # keep only rows with zeros
      filtered_mat <- as.matrix(df_numeric)
      filtered_mat <- filtered_mat[
        apply(filtered_mat == 0, 1, any),
        ,
        drop = FALSE
      ]

      validate(
        need(nrow(filtered_mat) > 0, "No zero values to display.")
      )

      # keep original intensities
      if (already_log2()) {

  missval <- filtered_mat

} else {

  missval <- log2(filtered_mat + 1)

}
cat("\n====================\n")
cat("HEATMAP DEBUG\n")
cat("====================\n")

cat("Rows:\n")
print(nrow(missval))

cat("Columns:\n")
print(ncol(missval))

cat("Zeros:\n")
print(sum(missval == 0, na.rm = TRUE))

cat("Positive values:\n")
print(sum(missval > 0, na.rm = TRUE))

cat("Range:\n")
print(range(missval, na.rm = TRUE))

      incProgress(
        0.50,
        detail = "Calculating ordering"
      )
      # ---------------------------------
      # Mean intensity using ONLY non-zero values
      # ---------------------------------

      row_means <- apply(
        missval,
        1,
        function(x) {

          x <- x[x > 0]

          if (length(x) == 0) {
            return(0)
          }

          mean(x)

        }
      )
      #rw order
      row_order_intensity <- order(
        row_means,
        decreasing = TRUE
      )
      row_means_sorted <- row_means[
      row_order_intensity
      ]
validate(
  need(
    length(row_means_sorted) > 0,
    "No proteins remain after filtering."
  )
)

      # intensity color scale
non_zero_vals <- missval[
  missval > 0
]
validate(
  need(
    length(non_zero_vals) > 0,
    "No proteins with non-zero intensities remain after filtering."
  )
)
vals <- c(
  0,
  median(non_zero_vals, na.rm = TRUE),
  max(non_zero_vals, na.rm = TRUE)
)

      if (length(unique(vals)) < 3) {

        vals <- seq(
          vals[1],
          vals[1] + 0.01,
          length.out = 3
        )

      }

      col_fun <- circlize::colorRamp2(
        vals,
        c(
          "white",
          "skyblue",
          "darkblue"
        )
      )
      # ---------------------------------
      # Reorder columns:
      # Temp1_REP1 Temp1_REP2 Temp1_REP3 ...
      # ---------------------------------

      col_info <- data.frame(
        colname = colnames(missval),
        stringsAsFactors = FALSE
      )

      col_info$temp <- as.numeric(
        sub(
          "_REP.*$",
          "",
          col_info$colname
        )
      )

      col_info$rep <- as.numeric(
        sub(
          ".*_REP",
          "",
          col_info$colname
        )
      )

        new_order <- order(
          col_info$temp,
          col_info$rep
        )

      missval_ordered <- missval[
        ,
        new_order,
        drop = FALSE
      ]

      missval_ordered_intensity <- missval_ordered[
        row_order_intensity,
        ,
        drop = FALSE
      ]
        missval_clustered <- missval[
          ,
          order(
            col_info$rep,
            col_info$temp
          ),
          drop = FALSE
        ]

      missval_clustered_intensity <- missval_clustered[
        row_order_intensity,
        ,
        drop = FALSE
      ]




annot_vals <- c(
  min(row_means_sorted, na.rm = TRUE),
  median(row_means_sorted, na.rm = TRUE),
  max(row_means_sorted, na.rm = TRUE)
)

if (length(unique(annot_vals)) < 3) {

  annot_vals <- seq(
    annot_vals[1],
    annot_vals[1] + 0.01,
    length.out = 3
  )

}

annot_col_fun <- circlize::colorRamp2(
  annot_vals,
  c(
    "white",
    "skyblue",
    "darkblue"
  )
)
intensity_annotation <- ComplexHeatmap::rowAnnotation(

MeanIntensity =
  ComplexHeatmap::anno_empty(
    width = grid::unit(18, "mm"),
    border = FALSE
  ),

show_annotation_name = FALSE
)



#intensity annotation only if intensity sorting is selected
left_annot <- NULL

if (
  input$row_sort_mode == "intensity"
) {

  left_annot <- intensity_annotation

}
      # switch depending on UI
      if (input$pattern_mode == "clustered") {

          plot_mat <- if (
            input$row_sort_mode == "cluster"
          ) {
            missval_clustered
          } else {
            missval_clustered_intensity
          }

          ht <- ComplexHeatmap::Heatmap(
            plot_mat,
            col = col_fun,
            left_annotation = left_annot,
            row_names_max_width = grid::unit(20, "mm"),
          column_names_side = "top",
          column_names_gp = grid::gpar(fontsize = 8),
          show_row_names = FALSE,
          show_column_names = TRUE,
          name = "Missing values pattern",
          heatmap_legend_param = list(
            title = "log2(Intensity)",
            direction = "horizontal",
            title_position = "leftcenter"
          ),
          cluster_rows =
  input$row_sort_mode == "cluster", # key difference
          cluster_columns = FALSE,   # keep temperature order
          use_raster = FALSE
        )
last_heatmap(ht)
      } else {

        plot_mat <- if (
          input$row_sort_mode == "cluster"
        ) {
          missval_ordered
        } else {
          missval_ordered_intensity
        }

        ht <- ComplexHeatmap::Heatmap(
          plot_mat,
          col = col_fun,
          left_annotation = left_annot,
          row_names_max_width = grid::unit(20, "mm"),
          heatmap_legend_param = list(
            title = "log2(Intensity)",
            direction = "horizontal",
            title_position = "leftcenter"
          ),

          column_names_side = "top",
          show_row_names = FALSE,
          show_column_names = TRUE,
          name = "Missing values pattern",
          column_names_gp = grid::gpar(fontsize = 8),
          cluster_columns = FALSE, 
          use_raster = FALSE,
          cluster_rows =
  input$row_sort_mode == "cluster",
        )
        last_heatmap(ht)
      }

#tick labels for intensity annotation
tick_positions <- seq(
  0,
  1,
  length.out = 11
)

tick_labels <- round(

  seq(
    min(row_means_sorted),
    max(row_means_sorted),
    length.out = 11
  ),

  1

)
last_tick_labels(tick_labels)
ComplexHeatmap::draw(
  ht,
  heatmap_legend_side = "top"
)
incProgress(
  0.90,
  detail = "Rendering heatmap"
)
if (input$row_sort_mode == "intensity") {

ComplexHeatmap::decorate_annotation(
  "MeanIntensity",
  {
grid::grid.text(

  label = "Average protein intensity",

  x = grid::unit(
    0,
    "mm"
  ),

  y = grid::unit(
    0.5,
    "npc"
  ),

  rot = 90,

  just = "centre",

  gp = grid::gpar(
    fontsize = 12,
    fontface = "bold"
  )

)
for (i in seq_along(tick_positions)) {

  y <- tick_positions[i]

grid::grid.lines(

  x = grid::unit(
    c(14, 17),
    "mm"
  ),

  y = grid::unit(
    c(y, y),
    "npc"
  ),

  gp = grid::gpar(
    col = "grey50",
    lwd = 0.3
  )

)

  grid::grid.text(

    label = sprintf("%.1f", tick_labels[i]),

    x = grid::unit(
      12,
      "mm"
    ),

    y = grid::unit(
      y,
      "npc"
    ),

    just = "right",

    gp = grid::gpar(
      fontsize = 10,
      fontface = "bold",
      col = "black"
    )

  )

}

  }
)
}
incProgress(
  1,
  detail = "Done"
)

    }
  )

})
output$download_missing_values_pattern <- downloadHandler(

  filename = function() {
    "missing_values_heatmap.png"
  },

  content = function(file) {

    req(last_heatmap())

    png(
      file,
      width = 2200,
      height = 1400,
      res = 250
    )
tick_positions <- seq(
  0,
  1,
  length.out = 11
)

tick_labels <- last_tick_labels()
    ComplexHeatmap::draw(
      last_heatmap(),
      heatmap_legend_side = "top"
    )
if (input$row_sort_mode == "intensity") {
  ComplexHeatmap::decorate_annotation(
  "MeanIntensity",
  {
    grid::grid.text(

      label = "Average protein intensity",

      x = grid::unit(
        0,
        "mm"
      ),

      y = grid::unit(
        0.5,
        "npc"
      ),

      rot = 90,

      just = "centre",

      gp = grid::gpar(
        fontsize = 12,
        fontface = "bold"
      )

    )

    for (i in seq_along(tick_positions)) {

      y <- tick_positions[i]

      grid::grid.lines(

        x = grid::unit(
          c(14, 17),
          "mm"
        ),

        y = grid::unit(
          c(y, y),
          "npc"
        ),

        gp = grid::gpar(
          col = "grey50",
          lwd = 0.3
        )

      )

      grid::grid.text(

        label = sprintf("%.1f", tick_labels[i]),

        x = grid::unit(
          12,
          "mm"
        ),

        y = grid::unit(
          y,
          "npc"
        ),

        just = "right",

        gp = grid::gpar(
          fontsize = 10,
          fontface = "bold",
          col = "black"
        )

      )

    }

  }
)}
    dev.off()

  }

)

    # -------------------------------
    # Plot: Total abundance per temperature
    # -------------------------------
output$total_abundance_plot <- renderPlot({

  df <- zero_stats$filtered()

  validate(
    need(!is.null(df), msg),
    need(nrow(df) > 0, "No proteins remain after filtering.")
  )


withProgress(
  message = "Calculating total abundance",
  value = 0,
  {

    incProgress(
      0.2,
      detail = "Reshaping dataset"
    )

    tmp <- df %>%
      tidyr::pivot_longer(
        cols = where(is.numeric),
        names_to = "Replicate",
        values_to = "Intensity"
      )

    incProgress(
      0.5,
      detail = "Extracting temperatures"
    )

    tmp <- tmp %>%
      tidyr::separate(
        Replicate,
        into = c(
          "temperature",
          "Replicate"
        ),
        sep = "_"
      )

    incProgress(
      0.8,
      detail = "Summarising abundance"
    )

    total_abundance_df <- tmp %>%
      dplyr::group_by(
        temperature
      ) %>%
      dplyr::summarise(
        Total_Abundance =
          sum(
            Intensity,
            na.rm = TRUE
          ),
        .groups = "drop"
      )

    incProgress(
      1,
      detail = "Done"
    )

  }
)


  # ---- Plot ----
  ggplot2::ggplot(total_abundance_df, ggplot2::aes(x = temperature, y = Total_Abundance)) +
    ggplot2::geom_line(group = 1, color = "steelblue4", linewidth = 1) +
    ggplot2::geom_point(color = "steelblue4") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Temperature",
      y = "Total abundance"
    )

    
})
output$download_total_abundance_plot <- downloadHandler(

  filename = function() {
    "total_abundance_plot.png"
  },

  content = function(file) {

    png(
      file,
      width = 1800,
      height = 1200,
      res = 200
    )

    df <- zero_stats$filtered()

    tmp <- df %>%
      tidyr::pivot_longer(
        cols = where(is.numeric),
        names_to = "Replicate",
        values_to = "Intensity"
      ) %>%
      tidyr::separate(
        Replicate,
        into = c("temperature", "Replicate"),
        sep = "_"
      )

    total_abundance_df <- tmp %>%
      dplyr::group_by(temperature) %>%
      dplyr::summarise(
        Total_Abundance = sum(Intensity, na.rm = TRUE),
        .groups = "drop"
      )

    p <- ggplot(
      total_abundance_df,
      aes(
        x = temperature,
        y = Total_Abundance
      )
    ) +
      geom_line(
        group = 1,
        color = "steelblue4"
      ) +
      geom_point(
        color = "steelblue4"
      ) +
      theme_minimal()

    print(p)

    dev.off()

  }

)
  })
}


####################################################################################
#normalization moduless
normalizationPanelModuleUI <- function(id, label) {

  ns <- NS(id)

  tagList(

    # ==============================
    # NORMALIZATION SETTINGS
    # ==============================

    fluidRow(
      box(
        title = tagList(
          "Normalization settings",

          tags$span(
            `data-toggle` = "tooltip",
            `data-placement` = "right",
            `data-html` = "true",

            title = HTML(
              "<b>Normalization methods</b><br><br>

              Select one or more normalization methods to apply before downstream analysis.<br><br>

              • <b>log2 transformation</b><br>
              Converts intensities to log2 scale and stabilizes variance.<br><br>

              • <b>Within replicate median normalization</b><br>
              Adjusts each replicate relative to the median signal at the same temperature.<br><br>

              • <b>Reference protein normalization</b><br>
              Normalizes all proteins relative to a reference protein.<br><br>

              • <b>Normalize to lowest temperature</b><br>
              For each replicate, all temperature points are divided by the abundance measured at the lowest temperature.<br><br>


              If no normalization method is selected, the raw filtered dataset is used unchanged."
            ),

            tags$i(
              class = "fa fa-question-circle",
              style = "margin-left:8px; color:#337ab7; cursor:pointer;"
            )
          )
        ),
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        uiOutput(ns("normalization_warning")),

        checkboxGroupInput(
          inputId = ns("normalization_methods"),
          label = "Select normalization methods to apply:",
          choices = c(
          "log2 transformation" = "log2",
          "Within replicate median normalization" = "rep",
          "Reference protein normalization" = "ref",
          "Normalize to lowest temperature" = "lowest_temp"
          ),
          selected = NULL,
          inline = TRUE
        ),

        conditionalPanel(
          condition = sprintf(
            "input['%s'].includes('ref')",
            ns("normalization_methods")
          ),

          tagList(

            fluidRow(

              column(
                width = 6,

                selectInput(
                  inputId = ns("ref_option"),
                  label = "Reference protein normalization method:",
                  choices = c(
                    "Manually select" = "manual",
                    "Program estimation" = "auto"
                  ),
                  selected = "auto"
                )

              ),

                column(
                  width = 6,

                  div(
                    class = "ref-protein-card",

                    tags$div(
                      class = "ref-protein-title",

                      icon("circle-info"),
                      " Reference protein"
                    ),

                    uiOutput(ns("ref_protein_used"))

                  )

                )

            ),

            conditionalPanel(
              condition = sprintf(
                "input['%s'] == 'manual'",
                ns("ref_option")
              ),

              fluidRow(
                column(
                  width = 6,

                  textInput(
                    inputId = ns("manual_ref_protein"),
                    label = "Enter reference protein name:",
                    placeholder = "e.g., P12345"
                  )

                )
              )

            )

          )

        ),

        actionButton(
          inputId = ns("apply_normalization"),
          label = "Apply selected normalization"
        )

      )
    ),

    # ==============================
    # BOXPLOT
    # ==============================

    fluidRow(

      box(
        title = "Intensity distribution across samples",
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        plotOutput(
          ns("normalized_data_preview"),
          height = 500
        ) %>% withSpinner(
          color = "#0EA5A5"
        ),

        tags$div(
          style = "height:8px;"
        ),

        downloadButton(
          ns("download_normalized_data_preview"),
          "Download plot"
        ),

        br(),
        br(),

        checkboxInput(
          inputId = ns("include_0_values"),
          label = "Include zero values in the plot",
          value = FALSE
        )

      )

    ),

    # ==============================
    # DENSITY PLOT
    # ==============================

    fluidRow(

      box(
        title = tagList(
          "Global intensity distribution",

          tags$span(
            `data-toggle` = "tooltip",
            `data-placement` = "right",
            `data-html` = "true",
            title = HTML(
              "<b>Global intensity distribution</b><br><br>

              This plot displays compares each protein thermal profile to the global average profile.<br><br>

              • Each colored line represents a single protein.<br>

              • The thick black line represents the average abundance trend across all proteins.<br><br>


              • Line color represents the distance between the protein profile and the global average profile.<br>

              • Distance is calculated using the squared Euclidean distance across all temperatures.<br><br>

              • Darker lines are closer to the global average profile.<br>"
            ),

            tags$i(
              class = "fa fa-question-circle",
              style = "margin-left:8px; color:#337ab7; cursor:pointer;"
            )
          )
        ),
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        plotOutput(
          ns("Total_density_plot"),
          height = 450
        ) %>% withSpinner(
          color = "#0EA5A5"
        ),

        tags$div(
          style = "height:8px;"
        ),

        downloadButton(
          ns("download_total_density_plot"),
          "Download plot"
        ),

        br(),
        br(),

        checkboxInput(
          ns("include_zeros"),
          "Include zero values in mean calculation",
          value = FALSE
        )

      )

    ),

    # ==============================
    # PCA
    # ==============================

    fluidRow(

      box(
      title = tagList(
        "Sample clustering and outlier detection",

        tags$span(
          `data-toggle` = "tooltip",
          `data-placement` = "right",
          `data-html` = "true",

          title = HTML(
            "<b>PCA and outlier detection</b><br><br>

            PCA is used to visualize sample similarity and identify potential outliers.<br><br>

            • <b>By sample</b><br>
            Displays individual samples and replicates.<br><br>

            • <b>By condition</b><br>
            Displays relationships between experimental conditions.<br><br>

            Samples can be selected directly from the PCA plot and moved to the removal list.<br>

            Removed samples can later be restored if needed.<br><br>

            <b>Important:</b><br>
            Whenever samples are removed or restored, normalization should be applied again so all downstream analyses use the updated dataset."
          ),

          tags$i(
            class = "fa fa-question-circle",
            style = "margin-left:8px; color:#337ab7; cursor:pointer;"
          )
        )
      ),
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        fluidRow(

          column(
            width = 6,

            selectInput(
              ns("pca_mode"),
              "PCA mode",
              choices = c(
                "By sample" = "sample",
                "By condition" = "condition"
              )
            )

          ),

          column(
            width = 6,

            br(),

            checkboxInput(
              ns("show_variance_plot"),
              "Show variance decomposition",
              value = FALSE
            )

          )

        ),

        plotlyOutput(
          ns("normalization_pca_plot"),
          height = "550px"
        ) %>% withSpinner(
color = "#0EA5A5"
),

        br(),

fluidRow(

  box(
    title = "Samples to remove",
    width = 6,
    status = "warning",
    solidHeader = TRUE,

    DT::dataTableOutput(
      ns("removed_samples_table_selected")
    ),

    br(),

    actionButton(
      ns("remove_selected"),
      "Remove selected",
      width = "100%"
    )

  ),

  box(
    title = "Removed samples",
    width = 6,
    status = "primary",
    solidHeader = TRUE,

    DT::dataTableOutput(
      ns("removed_samples_table_removed")
    ),

    br(),

    actionButton(
      ns("restore_selected"),
      "Restore selected",
      width = "100%"
    )

  )

)

      )

    ),

    # ==============================
    # CV PLOT
    # ==============================

    fluidRow(

      box(
      title = tagList(
        "Coefficient of variation plot",

        tags$span(
          `data-toggle` = "tooltip",
          `data-placement` = "right",
          `data-html` = "true",

          title = HTML(
            "<b>Coefficient of Variation (CV)</b><br><br>

            This plot summarizes the variability between replicates at each temperature.<br><br>

            The coefficient of variation is calculated as the standard deviation divided by the mean abundance for each protein.<br><br>

            Lower CV values indicate better agreement between replicates and generally reflect higher data quality.<br><br>

            Comparing CV distributions across temperatures can help identify experimental conditions associated with increased variability."
          ),

          tags$i(
            class = "fa fa-question-circle",
            style = "margin-left:8px; color:#337ab7; cursor:pointer;"
          )
        )
      ),
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        plotOutput(
          ns("CV_plot"),
          height = 450
        ) %>% withSpinner(
          color = "#0EA5A5"
        ),

        tags$div(
          style = "height:8px;"
        ),

        downloadButton(
          ns("download_cv_plot"),
          "Download plot"
        )

      )

    )

  )

}

# Helper: null coalescing
`%||%` <- function(a, b) if (!is.null(a)) a else b

normalizationPanelModuleServer <- function(id, zero_stats_data, dataset_org_data, all_datasets) {

  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    ref_protein_name <- reactiveVal(NULL)

    #for removing samples using PCA
selected_methods <- reactiveVal(character(0))
selected_samples <- reactiveValues(
  data = data.frame(
    Sample = character(),
    stringsAsFactors = FALSE
  )
)



removed_samples <- reactiveValues(
  data = data.frame(
    Sample = character(),
    stringsAsFactors = FALSE
  )
)
#removed sample message for renormalize
normalization_invalidated <- reactiveVal(FALSE)
# Tracks whether lowest temperature normalization has actually been applied
lowest_temp_normalized <- reactiveVal(FALSE)

observe({
  selected_samples$data
  removed_samples$data
})

    # 1) RAW DATA (canonical types & original order preserved)
    raw_data <- reactive({
      df <- req(zero_stats_data$filtered())


        # Remove selected samples
        
removed <- tryCatch(
  removed_samples$data,
  error = function(e) {
    cat("\n[DEBUG] removed_samples failed, resetting\n")
    data.frame(Sample = character())
  }
)


        if (nrow(removed) > 0) {

          sample_cols <- names(df)

          for (s in removed$Sample) {
            sample_cols <- sample_cols[!grepl(s, sample_cols)]
          }

          df <- df[, sample_cols, drop = FALSE]
        }
      df <- tibble::as_tibble(df)

      # enforce canonical types: first col as character, rest numeric
      first_col <- names(df)[1]

      # Remove Organism if present
      if ("Organism" %in% names(df)) {
        df <- df[, names(df) != "Organism", drop = FALSE]
      }

      #  Recompute columns AFTER removing Organism
      other_cols <- names(df)[-1]

      df[[first_col]] <- as.character(df[[first_col]])

      df <- dplyr::mutate(
        df,
        dplyr::across(
          dplyr::all_of(other_cols),
          ~ suppressWarnings(as.numeric(.))
        )
      )

      # keep original order explicitly
      attr(df, "orig_order") <- names(df)
      df
    })

    # 2) NORMALIZED DATA (only when button is clicked)
normalized_data <- eventReactive(input$apply_normalization, {

  df <- raw_data()
  methods <- isolate(input$normalization_methods) %||% character(0)

  withProgress(
    message = "Applying normalization",
    value = 0,
    {

      incProgress(
        0.05,
        detail = "Preparing dataset"
      )

      # If NO method selected, return EXACTLY the raw snapshot (types + order preserved)

        if (length(methods) == 0) {
          ref_protein_name(NULL)
          return(raw_data())
        }

      orig_order <- attr(df, "orig_order") %||% names(df)
      first_col  <- orig_order[1]
      temp_cols <- setdiff(orig_order[-1], "Organism")
normalize_to_lowest_temperature <- function(df, first_col) {

temp_cols <- names(df)[-1]

temps_str <- stringr::str_extract(
temp_cols,
"^\\d+\\.?\\d*"
)

temps_num <- as.numeric(temps_str)

replicates <- unique(
stringr::str_extract(
temp_cols,
"REP\\d+"
)
)

lowest_temp_str <- temps_str[
which.min(temps_num)
]

out <- df

for(rep in replicates){

rep_cols <- grep(
rep,
names(out),
value = TRUE
)

ref_col <- grep(
paste0(
"^",
lowest_temp_str,
"_",
rep,
"$"
),
names(out),
value = TRUE
)

if(length(ref_col) == 0)
next

ref_vals <- out[[ref_col]]

out[rep_cols] <- sweep(
out[rep_cols],
1,
ref_vals,
"/"
)

bad_rows <- is.na(ref_vals) | ref_vals == 0

if(any(bad_rows)){
out[rep_cols][bad_rows, ] <- NA
}
}
out
}

      # ---- LOG2 normalization ----
if ("log2" %in% methods) {

  incProgress(
    0.25,
    detail = "Applying log2 transformation"
  )
        df <- df %>%
          dplyr::mutate(
            dplyr::across(
              dplyr::all_of(temp_cols),
              ~ {
                v <- suppressWarnings(log2(.))
                v[!is.finite(v)] <- 0  # replace -Inf/Inf/NaN with 0
                v
              }
            )
          )
      }

      # ---- Replicate-based normalization ----
if ("rep" %in% methods) {

  incProgress(
    0.40,
    detail = "Applying replicate normalization"
  )
        # Long format -> compute normalization -> wide -> restore order
        df_long <- df %>%
          tidyr::pivot_longer(
            cols = dplyr::all_of(temp_cols),
            names_to = "Replicate",
            values_to = "Intensity"
          ) %>%
          tidyr::separate(Replicate, into = c("Temperature", "Replicate"), sep = "_", remove = FALSE)

        # Temperature kept as character to avoid accidental factor ordering issues
        df_long <- df_long %>%
          dplyr::mutate(Temperature = as.character(Temperature))

        # Group medians (skip zeros and non-finite)
        df_long <- df_long %>%
          dplyr::group_by(Temperature) %>%
          dplyr::mutate(
            Group_average = dplyr::coalesce(stats::median(Intensity[is.finite(Intensity) & Intensity != 0], na.rm = TRUE), 0)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::group_by(Temperature, Replicate) %>%
          dplyr::mutate(
            Rep_average = dplyr::coalesce(stats::median(Intensity[is.finite(Intensity) & Intensity != 0], na.rm = TRUE), 0)
          ) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(
            M_difference = Rep_average - Group_average,
            Normalized   = dplyr::case_when(
              Intensity == 0 ~ 0,
              is.finite(Intensity) ~ Intensity - M_difference,
              TRUE ~ NA_real_
            )
          )

        df <- df_long %>%
          dplyr::mutate(Replicate_full = paste0(Temperature, "_", Replicate)) %>%  # original name already contains temp + rep
          dplyr::select(dplyr::all_of(first_col), Replicate_full, Normalized) %>%
          tidyr::pivot_wider(names_from = Replicate_full, values_from = Normalized)

        # Restore original order of columns exactly
        df <- df[, intersect(orig_order, names(df)), drop = FALSE]
        # Ensure numerics
        df <- dplyr::mutate(df, dplyr::across(dplyr::all_of(temp_cols), ~ suppressWarnings(as.numeric(.))))
      }
# ---- Lowest temperature normalization ----

if ("lowest_temp" %in% methods) {

incProgress(
0.65,
detail = "Normalizing to lowest temperature"
)

df <- normalize_to_lowest_temperature(
df,
first_col
)

}
      # ---- Reference-based normalization ----
if ("ref" %in% methods) {

  incProgress(
    0.85,
    detail = "Applying reference protein normalization"
  )



        if (input$ref_option == "manual") {
          ref_protein <- input$manual_ref_protein
          if (!is.null(ref_protein) && ref_protein %in% df[[first_col]]) {
            reference_row <- df %>% dplyr::filter(.data[[first_col]] == ref_protein)
            reference_values <- reference_row %>% dplyr::select(dplyr::all_of(temp_cols))
            ref_vec <- as.numeric(reference_values[1, ])
            names(ref_vec) <- temp_cols

            # Divide each temp column by its reference value (safe divide)
            df <- df %>%
              dplyr::mutate(
                dplyr::across(
                  dplyr::all_of(temp_cols),
                  ~ {
                    denom <- unname(ref_vec[[cur_column()]])
                    if (isTRUE(denom == 0) || is.na(denom)) return(.x)
                    .x / denom
                  }
                )
              )

            ref_protein_name(ref_protein)
            # Remove the reference protein from data
            df <- df %>% dplyr::filter(.data[[first_col]] != ref_protein)
          } else {
            showNotification("Reference protein not found in data.", type = "error")
            ref_protein_name(NULL)
          }
        } else if (input$ref_option == "auto") {
          # Find most stable (lowest SD across temp cols)
          df_sd <- df %>%
            dplyr::rowwise() %>%
            dplyr::mutate(std_dev = stats::sd(c(dplyr::c_across(dplyr::all_of(temp_cols))), na.rm = TRUE)) %>%
            dplyr::ungroup() %>%
            dplyr::arrange(std_dev)

          most_stable <- df_sd %>%
            dplyr::slice(1)

          reference_values <- most_stable %>%
            dplyr::select(dplyr::all_of(temp_cols))

          ref_vec <- as.numeric(reference_values[1, ])
          names(ref_vec) <- temp_cols

          df <- df %>%
            dplyr::mutate(
              dplyr::across(
                dplyr::all_of(temp_cols),
                ~ {
                  denom <- unname(ref_vec[[cur_column()]])
                  if (isTRUE(denom == 0) || is.na(denom)) return(.x)
                  .x / denom
                }
              )
            )

          ref_protein_name(most_stable[[first_col]])
          df <- df %>% dplyr::filter(.data[[first_col]] != most_stable[[first_col]])
        }

        # After ref normalization, preserve types & order
        df <- df[, intersect(orig_order, names(df)), drop = FALSE]
        valid_cols <- intersect(setdiff(orig_order, first_col), names(df))

          df <- dplyr::mutate(
            df,
            dplyr::across(
              dplyr::all_of(valid_cols),
              ~ suppressWarnings(as.numeric(.))
            )
          )
      }

      # Final: restore original order (id + temp columns)
      df <- df[, intersect(orig_order, names(df)), drop = FALSE]

      # Keep canonical type for id column
      df[[first_col]] <- as.character(df[[first_col]])

      incProgress(
        1,
        detail = "Finished"
      )

      df

    }
  )

})

    # 3) CURRENT DATA used by downstream modules
   




# --- CURRENT DATA STATE ---
current_data <- reactiveVal(NULL)

observe({
  current_data(raw_data())
})

observeEvent(input$apply_normalization, {

  methods <- input$normalization_methods %||% character(0)

  if (length(methods) == 0) {
    current_data(raw_data())
  } else {
    current_data(normalized_data())
  }

  lowest_temp_normalized(
    "lowest_temp" %in% methods
  )

  normalization_invalidated(FALSE)

})

plot_data <- reactive({
  df <- current_data()
  validate(
    need(!is.null(df), "Data not available yet."),
    need(nrow(df) > 0, "No rows in data."),
    need(ncol(df) > 1, "No measurement columns in data.")
  )
  df
})

# --- PCA click handling ---
pca_click <- reactive({
  plotly::event_data(
    "plotly_click",
    source = ns("pca")
  )
})

clicked_sample <- reactive({
  d <- pca_click()
  req(d)
  req(d$key)
  as.character(d$key)
})

observeEvent(pca_click(), {

  d <- pca_click()
  if (is.null(d) || is.null(d$key)) return()

  sample_id <- as.character(d$key)

  current <- selected_samples$data

  if (!(sample_id %in% current$Sample)) {

    selected_samples$data <- rbind(
      current,
      data.frame(Sample = sample_id, stringsAsFactors = FALSE)
    )

    cat("\n[DEBUG CLICK ADDED]:", sample_id, "\n")
  }

}, ignoreNULL = TRUE)


# --- RESTORE ---
observeEvent(input$restore_selected, {

  isolate({

    idx <- input$removed_samples_table_removed_rows_selected
    if (is.null(idx)) return(NULL)

    df <- removed_samples$data

    if (length(idx) > 0) {
      df <- df[-idx, , drop = FALSE]
    }

    removed_samples$data <- df
    normalization_invalidated(TRUE)
    lowest_temp_normalized(FALSE)
  })

})

# --- REMOVE ---
observeEvent(input$remove_selected, {

 

  isolate({

    

    current <- selected_samples$data

   

    if (is.null(current) || nrow(current) == 0) {
      cat("No samples selected → exit\n")
      return(NULL)
    }

    removed <- removed_samples$data

    removed_samples$data <- unique(rbind(removed, current))

    normalization_invalidated(TRUE)
    lowest_temp_normalized(FALSE)
    selected_samples$data <- data.frame(
      Sample = character(),
      stringsAsFactors = FALSE
    )

   

  })

}, priority = 1)

# --- Tables ---

output$removed_samples_table_selected <- DT::renderDataTable({

  # ✅ force reactive dependency
  nrow(selected_samples$data)

  selected_samples$data
})

output$removed_samples_table_removed <- DT::renderDataTable({

  # ✅ same trick here
  nrow(removed_samples$data)

  removed_samples$data
})



output$normalization_warning <- renderUI({

  if (!normalization_invalidated()) {
    return(NULL)
  }

  div(
    style = paste(
      "padding:10px;",
      "margin-bottom:15px;",
      "background-color:#fff3cd;",
      "border:1px solid #ffeeba;",
      "border-radius:4px;",
      "color:#856404;",
      "font-weight:bold;"
    ),

    "Samples were removed or restored in the PCA section. ",
    "Please click 'Apply selected normalization' again to recompute all normalization steps."
  )

})

# 5) Reference protein info (only meaningful if ref normalization ran)
output$ref_protein_used <- renderUI({
  req(ref_protein_name())  # Only require the ref protein name

  # Try to get organism data; it may be NULL
  org <- dataset_org_data()

  # Base line: always show the reference protein
  base_line <- paste0("<strong>Reference protein:</strong> ",
                      htmltools::htmlEscape(ref_protein_name()), "<br>")

  # If organism data is NULL, stop here and just show the base line
  if (is.null(org)) {
    return(htmltools::HTML(base_line))
  }

  # Defensive: ensure the expected columns exist before filtering/selecting
  needed_cols <- c("Entry", "Protein.names", "Gene.Names")
  if (!all(needed_cols %in% colnames(org))) {
    # Columns missing => show only base line
    return(htmltools::HTML(base_line))
  }

  # Filter for the reference protein
  organism_df <- org %>%
    dplyr::filter(.data$Entry == ref_protein_name()) %>%
    dplyr::select(.data$Protein.names, .data$Gene.Names)

  # If no match found, return only the base line
  if (nrow(organism_df) == 0) {
    return(htmltools::HTML(base_line))
  }

  # If multiple matches, use the first (or customize as needed)
  protein_name <- organism_df$Protein.names[[1]]
  gene_name    <- organism_df$Gene.Names[[1]]

  # Build the full HTML (escape user-facing strings)
  htmltools::HTML(paste0(
    base_line,
    "<strong>Protein Name:</strong> ", htmltools::htmlEscape(protein_name), "<br>",
    "<strong>Gene Name:</strong> ",    htmltools::htmlEscape(gene_name)
  ))
})




# Boxplot
output$normalized_data_preview <- renderPlot({

  withProgress(
    message = "Generating boxplot",
    value = 0,
    {

      incProgress(
        0.20,
        detail = "Loading normalized data"
      )

      df <- plot_data()

      boxinput <- input$include_0_values
      methods <- selected_methods()

incProgress(
  0.50,
  detail = "Reshaping dataset"
)

      df_long <- df %>%
        pivot_longer(cols = 2:last_col(), names_to = "Replicate", values_to = "Intensity") %>%
        tidyr::extract(Replicate, into = c("temperature", "Replicate"), regex = "([^_]+)_([^_]+)") %>%
        mutate(
          temperature_num = as.numeric(str_extract(temperature, "\\d+")),
          ID = paste(temperature, Replicate, sep = "_")
        )


      if (!boxinput) {
        df_long <- df_long %>% filter(Intensity != 0)
      }

      y_label <- if ("log2" %in% methods) "Log2(Intensity)" else "Intensity"
incProgress(
  0.90,
  detail = "Rendering plot"
)
      p <- ggplot(df_long, aes(x = ID, y = Intensity, fill = temperature_num)) +
        geom_boxplot() +
        scale_fill_gradient(low = "blue", high = "red", name = "Temperature") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
        ylab(y_label) +
        xlab("Samples")


      incProgress(
        1,
        detail = "Done"
      )

      p

    }
  )
     
    })


output$download_normalized_data_preview <- downloadHandler(

  filename = function() {
    paste0(
      "Intensity_Distribution_",
      Sys.Date(),
      ".png"
    )
  },

  content = function(file) {

    df <- plot_data()

    boxinput <- input$include_0_values
    methods <- selected_methods()

    df_long <- df %>%
      pivot_longer(
        cols = 2:last_col(),
        names_to = "Replicate",
        values_to = "Intensity"
      ) %>%
      tidyr::extract(
        Replicate,
        into = c("temperature", "Replicate"),
        regex = "([^_]+)_([^_]+)"
      ) %>%
      mutate(
        temperature_num = as.numeric(
          stringr::str_extract(
            temperature,
            "\\d+"
          )
        ),
        ID = paste(
          temperature,
          Replicate,
          sep = "_"
        )
      )

    if (!boxinput) {

      df_long <- df_long %>%
        filter(Intensity != 0)

    }

    y_label <- if (
      "log2" %in% methods
    ) {
      "Log2(Intensity)"
    } else {
      "Intensity"
    }

    p <- ggplot(
      df_long,
      aes(
        x = ID,
        y = Intensity,
        fill = temperature_num
      )
    ) +
      geom_boxplot() +
      scale_fill_gradient(
        low = "blue",
        high = "red",
        name = "Temperature"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )
      ) +
      ylab(y_label) +
      xlab("Samples")

    ggsave(
      file,
      p,
      width = 12,
      height = 8,
      dpi = 300
    )

  }

)

    # PCA
    run_tpp_pca_multi <- function(dataset_list) {

  # dataset_list: named list of dataframes

  # ---- add condition ----
    S_all <- dplyr::bind_rows(
      lapply(names(dataset_list), function(cond) {

        df <- dataset_list[[cond]]

        # ✅ Remove Organism BEFORE anything else
        if ("Organism" %in% colnames(df)) {
          df <- df[, colnames(df) != "Organism", drop = FALSE]
        }

        first_col <- names(df)[1]

        df %>%
          dplyr::rename(Accession = !!first_col) %>%
          dplyr::mutate(Condition = cond)
      })
    )

  # ---- replace NA ----
  S_all[is.na(S_all)] <- 0



  # ---- reshape ----
  S_sample <- S_all %>%
    tidyr::pivot_longer(
      -c(Accession, Condition),
      names_to = c("Temperature", "Replicate"),
      names_sep = "_",
      values_to = "Intensity"
    ) %>%
    tidyr::unite(Sample, Condition, Temperature, Replicate, sep = "__") %>%
    tidyr::pivot_wider(
      names_from = Accession,
      values_from = Intensity,
      values_fill = 0
    )

  # ---- matrix ----
mat <- S_sample %>%
  dplyr::select(-Sample) %>%
  as.matrix()

rownames(mat) <- S_sample$Sample

mode(mat) <- "numeric"

# ✅ REMOVE CONSTANT COLUMNS (CRITICAL)
keep_cols <- apply(mat, 2, function(x) {
  sd_val <- sd(x, na.rm = TRUE)
  !is.na(sd_val) && sd_val > 0
})

keep_cols <- keep_cols[seq_len(ncol(mat))]

mat <- mat[, keep_cols, drop = FALSE]

# ✅ SAFETY CHECK
if (ncol(mat) < 2) {
  stop("Not enough variable proteins for PCA (after filtering)")
}

# ---- PCA ----
pca <- stats::prcomp(mat, center = TRUE, scale. = TRUE)

  # ---- metadata ----
  meta <- S_sample %>%
    tidyr::separate(Sample, into = c("Condition","Temperature","Replicate"), sep = "__") %>%
    dplyr::mutate(Temperature = as.numeric(Temperature))

  pca_var <- summary(pca)$importance[2, ]

  pc1_var <- round(pca_var[1] * 100, 1)
  pc2_var <- round(pca_var[2] * 100, 1)


  df_plot <- data.frame(
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    Condition = meta$Condition,
    Replicate = meta$Replicate,
    Temperature = meta$Temperature
  )


p <- ggplot(
  df_plot,
  aes(
    PC1,
    PC2,
    color = Condition,
    shape = Replicate,
    text = paste0(
      "Condition: ", Condition,
      "<br>Replicate: ", Replicate,
      "<br>Temperature: ", Temperature, " °C"
    ),
    key = paste(Temperature, Replicate, sep = "_")
  )
) +

  geom_point(size = 3) +

scale_color_discrete(
  name = "Condition"
) +

guides(
  color = guide_legend(order = 1),
  shape = "none"
) +

  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    title = "PCA by condition"
  ) +

  theme_classic()



#  variance explained
pca_var <- summary(pca)$importance[2, ]

#  batch effect
batch_var <- sapply(seq_along(pca_var), function(i) {
  fit <- lm(pca$x[, i] ~ meta$Replicate)
  summary(fit)$r.squared
})

# condition effect
condition_var <- sapply(seq_along(pca_var), function(i) {
  fit <- lm(pca$x[, i] ~ meta$Condition)
  summary(fit)$r.squared
})


n_show <- min(10, length(pca_var))

var_df <- data.frame(
  PC = factor(paste0("PC", 1:n_show), levels = paste0("PC", 1:n_show)),
  TotalVariance = pca_var[1:n_show],
  BatchVariance = (batch_var * pca_var)[1:n_show],
  ConditionVariance = (condition_var * pca_var)[1:n_show]
)


var_df$ResidualVariance <- var_df$TotalVariance -
  var_df$BatchVariance -
  var_df$ConditionVariance

var_df_long <- var_df %>%
  tidyr::pivot_longer(
    cols = c(BatchVariance, ConditionVariance, ResidualVariance),
    names_to = "Component",
    values_to = "Variance"
  )

p_var <- ggplot(var_df_long, aes(x = PC, y = Variance, fill = Component)) +
  geom_col() +
  theme_classic() +
  scale_fill_manual(
    values = c(
      BatchVariance = "steelblue4",
      ConditionVariance = "firebrick4",
      ResidualVariance = "grey70"
    ),
    name = "Variance source"
  ) +
  labs(title = "Variance decomposition per PC")
#return both plots
  return(list(
    plot = p,
    variance_plot = p_var  # placeholder for now
  ))

}

    run_tpp_pca <- function(df, scale_data = TRUE, impute_zero = TRUE) {

      # ---- Prepare data ----
      S <- df
      first_col <- names(S)[1]

      if (impute_zero) {
        S[is.na(S)] <- 0
      }

      # ---- Reshape to sample-level ----
      S_sample <- S %>%
        tidyr::pivot_longer(
          -all_of(first_col),
          names_to = c("Temperature", "Replicate"),
          names_sep = "_",
          values_to = "Intensity"
        ) %>%
        tidyr::unite(Sample, Temperature, Replicate, sep = "_") %>%
        tidyr::pivot_wider(
          names_from = all_of(first_col),
          values_from = Intensity
        )

      # ---- Matrix for PCA ----

        mat <- S_sample %>%
          dplyr::select(-Sample) %>%
          as.matrix()

        rownames(mat) <- S_sample$Sample

        # Step 1: ensure matrix is numeric (important)
        mode(mat) <- "numeric"

        # Step 2: ensure deterministic column order
        cn <- colnames(mat)
        mat <- mat[, order(cn), drop = FALSE]

        # Step 3: recompute safe keep_cols AFTER ordering
        keep_cols <- apply(mat, 2, function(x) {
          sd_val <- sd(x, na.rm = TRUE)
          !is.na(sd_val) && sd_val > 0
        })

        # Step 4: guard against mismatch
        keep_cols <- keep_cols[seq_len(ncol(mat))]

        # Step 5: subset safely
        mat <- mat[, keep_cols, drop = FALSE]

        # safety check
        validate(
          need(ncol(mat) > 2, "Not enough variable proteins for PCA")
        )


      # ---- PCA ----
      pca <- stats::prcomp(mat, center = TRUE, scale. = scale_data)

      # ---- Metadata ----
      meta <- S_sample %>%
        tidyr::separate(Sample, into = c("Temperature", "Replicate"),
         sep = "_", extra = "merge", fill = "right") %>%
        dplyr::mutate(
          Temperature = as.numeric(Temperature)
        )

      # ---- Variance explained ----
      pca_var <- summary(pca)$importance[2, ]

      pc1_var <- round(pca_var[1] * 100, 1)
      pc2_var <- round(pca_var[2] * 100, 1)

      # ---- PCA plot ----
      df_plot <- data.frame(
        PC1 = pca$x[, 1],
        PC2 = pca$x[, 2],
        Temperature = meta$Temperature,
        Replicate = meta$Replicate
      )






      # Compute batch effect 
      batch_var <- sapply(seq_along(pca_var), function(i) {
        fit <- lm(pca$x[, i] ~ meta$Replicate)
        summary(fit)$r.squared
      })
      total_batch_effect <- sum(batch_var * pca_var)
      
        p <- ggplot(df_plot, aes(
          PC1, PC2,
          color = Temperature,
          shape = Replicate,
          text = paste0(
            "Replicate: ", Replicate,
            "<br>Temperature: ", Temperature
          ),
          key = paste(Temperature, Replicate, sep = "_")
        ))+
        geom_point(size = 3) +
        scale_color_viridis_c(name = "Temperature") +
        labs(
          x = paste0("PC1 (", pc1_var, "%)"),
          y = paste0("PC2 (", pc2_var, "%)"),
          title = "PCA by sample"
        ) +
        theme_classic()

          #  Variance decomposition

        n_show <- min(10, length(pca_var))

        var_df <- data.frame(
          PC = factor(paste0("PC", 1:n_show), levels = paste0("PC", 1:n_show)),
          TotalVariance = pca_var[1:n_show],
          BatchVariance = (batch_var * pca_var)[1:n_show]
        )


          var_df$ResidualVariance <- var_df$TotalVariance - var_df$BatchVariance

          var_df_long <- var_df %>%
            tidyr::pivot_longer(
              cols = c(BatchVariance, ResidualVariance),
              names_to = "Component",
              values_to = "Variance"
            )

          p_var <- ggplot(var_df_long, aes(x = PC, y = Variance, fill = Component)) +
            geom_col() +
            theme_classic() +
            scale_fill_manual(
              values = c(BatchVariance = "steelblue4", ResidualVariance = "grey70"),
              name = "Variance source",
              labels = c(BatchVariance = "Batch", ResidualVariance = "Residual")
            ) +
            labs(
              y = "Proportion of total variance",
              title = "Variance decomposition per PC"
            ) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))


          # ✅ RETURN BOTH
          return(list(
            plot = p,
            variance_plot = p_var,
            total_batch_effect = total_batch_effect
          ))

    }

    #plot PCA
    output$normalization_pca_plot <- renderPlotly({
withProgress(
  message = "Generating PCA",
  value = 0,
  {

      mode <- input$pca_mode
      show_var <- input$show_variance_plot
incProgress(
  0.10,
  detail = "Loading normalized data"
)
      if (mode == "sample") {

        df <- plot_data()



        if ("Organism" %in% colnames(df)) {
          df <- df[, colnames(df) != "Organism", drop = FALSE]
        }

        validate(
          need(ncol(df) > 2, "Not enough data for PCA"),
          need(nrow(df) > 2, "Too few proteins for PCA")
        )

    incProgress(
  0.40,
  detail = "Computing PCA by sample"
)    
res <- tryCatch({

  run_tpp_pca(df)

}, error = function(e) {



  return(NULL)
})

validate(
  need(!is.null(res), "PCA failed — check debug log.")
)
incProgress(
  0.75,
  detail = "Preparing visualization"
)

        if (show_var) {

          validate(
            need(!is.null(res$variance_plot), "Variance plot not available")
          )

          # ✅ return ggplot directly
          return(res$variance_plot)

        } else {

          
incProgress(
  0.95,
  detail = "Rendering PCA"
)
p <- ggplotly(
  res$plot,
  tooltip = "text",
  source = ns("pca")
)


p <- plotly::event_register(p, "plotly_click")
incProgress(
  1,
  detail = "Done"
)
return(p)

        }

      } else {

datasets <- all_datasets()

dataset_list <- lapply(names(datasets), function(name) {

  df <- as.data.frame(datasets[[name]])

  # ✅ apply same preprocessing as raw_data
  if ("Organism" %in% colnames(df)) {
    df <- df[, colnames(df) != "Organism", drop = FALSE]
  }

  # ✅ apply normalization if active
  methods <- input$normalization_methods %||% character(0)

  if (length(methods) > 0) {

    df_tmp <- df
    first_col <- names(df_tmp)[1]
    temp_cols <- names(df_tmp)[-1]

    # ---- LOG2 ----
    if ("log2" %in% methods) {
      df_tmp[temp_cols] <- lapply(df_tmp[temp_cols], function(x) {
        v <- suppressWarnings(log2(as.numeric(x)))
        v[!is.finite(v)] <- 0
        v
      })
    }

    # 👉 (optional: extend later for rep/ref)

    return(df_tmp)
  }

  df
})


        names(dataset_list) <- names(datasets)

        validate(
          need(length(dataset_list) > 1, "Need multiple datasets for condition PCA")
        )
incProgress(
  0.40,
  detail = "Computing PCA across conditions"
)
        res <- run_tpp_pca_multi(dataset_list)
incProgress(
  0.75,
  detail = "Preparing visualization"
)
        if (show_var) {

          validate(
            need(!is.null(res$variance_plot), "Variance plot not available")
          )

          return(res$variance_plot)

        } else {

       incProgress(
  0.95,
  detail = "Rendering PCA"
)   
p <- plotly::ggplotly(
  res$plot,
  tooltip = "text",
  source = ns("pca")
)

p <- plotly::event_register(p, "plotly_click")
incProgress(
  1,
  detail = "Done"
)
return(p)

        }

      }
    }
  )

})




# Selected (pending removal)
output$removed_samples_table_selected <- DT::renderDataTable({
  selected_samples$data
})

# Removed (applied)
output$removed_samples_table_removed <- DT::renderDataTable({
  removed_samples$data
})




    # Total Density Plot
output$Total_density_plot <- renderPlot({

  withProgress(
    message = "Generating global intensity distribution",
    value = 0,
    {

      incProgress(
        0.20,
        detail = "Loading normalized data"
      )

      df <- plot_data()
incProgress(
  0.40,
  detail = "Reshaping dataset"
)
      df_long <- df %>%
        pivot_longer(cols = -1, names_to = "Replicate", values_to = "Intensity") %>%
        separate(Replicate, into = c("temperature", "Replicate"), sep = "_")

      df_long$temperature <- as.numeric(df_long$temperature)

      if (!input$include_zeros) {
        df_long <- df_long %>%
          group_by(across(1), temperature) %>%
          filter(!(all(Intensity == 0))) %>%
          filter(Intensity != 0) %>%
          ungroup()
      }
incProgress(
  0.70,
  detail = "Calculating protein trends"
)

      df2 <- df_long %>%
        group_by(across(1), temperature) %>%
        summarise(Mean_Intensity = mean(Intensity, na.rm = TRUE), .groups = "drop")

      average_trend <- df2 %>%
        group_by(temperature) %>%
        summarise(Avg_Intensity = mean(Mean_Intensity, na.rm = TRUE), .groups = "drop")

      df_wide <- df2 %>%
        pivot_wider(names_from = temperature, values_from = Mean_Intensity)

      avg_vector <- colMeans(df_wide[,-1], na.rm = TRUE)

      df_wide$distance <- apply(df_wide[,-1], 1, function(x) sum((x - avg_vector)^2, na.rm = TRUE))
      
      df_wide$distance <- log(pmax(df_wide$distance, 1e-10))
      
      df_wide$distance <- df_wide$distance - min(df_wide$distance, na.rm = TRUE)
      df_wide$distance <- df_wide$distance / max(df_wide$distance, na.rm = TRUE) * 0.1
      df_wide$inverse_distance <- 0.1 - df_wide$distance

      df2 <- df2 %>%
        left_join(df_wide %>% select(1, distance, inverse_distance), by = names(df2)[1])
incProgress(
  0.95,
  detail = "Rendering plot"
)
     p <- ggplot(df2, aes(x = temperature, y = Mean_Intensity, group = !!sym(names(df2)[1]))) +
        geom_line(aes(color = distance, size = inverse_distance), alpha = 0.3) +
        geom_line(data = average_trend, aes(x = temperature, y = Avg_Intensity, group = 1),
                  color = "black", size = 2, linetype = "solid", alpha = 0.5) +
        scale_color_viridis_c(option = "D") +
        scale_size_continuous(range = c(0.2, 1.5)) +
        labs(
          title = "",
          x = "Temperature (°C)",
          y = "Relative Abundance",
          color = "Distance",
          size = "Inverse Distance"
        ) +
        theme_minimal(base_size = 14) +
        theme(legend.position = "none")


incProgress(
  1,
  detail = "Done"
)

p
   
    }
  )

})
output$download_total_density_plot <- downloadHandler(

  filename = function() {
    paste0("Global_Intensity_Distribution_", Sys.Date(), ".png")
  },

  content = function(file) {

    df <- plot_data()

    df_long <- df %>%
      pivot_longer(
        cols = -1,
        names_to = "Replicate",
        values_to = "Intensity"
      ) %>%
      separate(
        Replicate,
        into = c("temperature", "Replicate"),
        sep = "_"
      )

    df_long$temperature <- as.numeric(df_long$temperature)

    if (!input$include_zeros) {

      df_long <- df_long %>%
        group_by(across(1), temperature) %>%
        filter(!(all(Intensity == 0))) %>%
        filter(Intensity != 0) %>%
        ungroup()

    }

    df2 <- df_long %>%
      group_by(across(1), temperature) %>%
      summarise(
        Mean_Intensity = mean(Intensity, na.rm = TRUE),
        .groups = "drop"
      )

    average_trend <- df2 %>%
      group_by(temperature) %>%
      summarise(
        Avg_Intensity = mean(Mean_Intensity, na.rm = TRUE),
        .groups = "drop"
      )

    df_wide <- df2 %>%
      pivot_wider(
        names_from = temperature,
        values_from = Mean_Intensity
      )

    avg_vector <- colMeans(
      df_wide[, -1],
      na.rm = TRUE
    )

    df_wide$distance <- apply(
      df_wide[, -1],
      1,
      function(x)
        sum((x - avg_vector)^2, na.rm = TRUE)
    )

    df_wide$distance <- log(pmax(df_wide$distance, 1e-10))

    df_wide$distance <- df_wide$distance -
      min(df_wide$distance, na.rm = TRUE)

    df_wide$distance <- df_wide$distance /
      max(df_wide$distance, na.rm = TRUE) *
      0.1

    df_wide$inverse_distance <- 0.1 -
      df_wide$distance

    df2 <- df2 %>%
      left_join(
        df_wide %>%
          select(1, distance, inverse_distance),
        by = names(df2)[1]
      )

    p <- ggplot(
      df2,
      aes(
        x = temperature,
        y = Mean_Intensity,
        group = !!sym(names(df2)[1])
      )
    ) +
      geom_line(
        aes(
          color = distance,
          size = inverse_distance
        ),
        alpha = 0.3
      ) +
      geom_line(
        data = average_trend,
        aes(
          x = temperature,
          y = Avg_Intensity,
          group = 1
        ),
        color = "black",
        size = 2,
        alpha = 0.5
      ) +
      scale_color_viridis_c(option = "D") +
      scale_size_continuous(
        range = c(0.2, 1.5)
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")

    ggsave(
      file,
      p,
      width = 12,
      height = 7,
      dpi = 300
    )

  }

)

    # CV Plot
output$CV_plot <- renderPlot({

  withProgress(
    message = "Generating CV plot",
    value = 0,
    {

      incProgress(
        0.20,
        detail = "Loading normalized data"
      )

      df <- plot_data()
      method <- selected_methods()
incProgress(
  0.40,
  detail = "Preparing replicate matrix"
)
      df_long <- df %>%
        pivot_longer(cols = -1, names_to = "Replicate", values_to = "Intensity") %>%
        separate(Replicate, into = c("temperature", "Replicate"), sep = "_")

      df_split <- df_long %>%
        group_by(temperature) %>%
        group_split()

      calc_cv_by_temp <- function(df, method = "raw") {
        wide <- df %>%
          select(1, Replicate, Intensity) %>%
          pivot_wider(names_from = Replicate, values_from = Intensity)

        id_col <- names(wide)[1]
        mat <- wide %>% select(-all_of(id_col))
        mat[mat == 0] <- NA
        mat <- as.matrix(mat)

        cv <- apply(mat, 1, function(x) {
          x_clean <- x[!is.na(x)]
          if (length(x_clean) >= 2) {
            if (method == "log2") {
              x_log <- log2(x_clean)
              sqrt(exp(sd(x_log)^2) - 1) * 100
            } else {
              (sd(x_clean) / mean(x_clean)) * 100
            }
          } else {
            NA
          }
        })

        tibble(
          !!id_col := wide[[id_col]],
          temperature = unique(df$temperature),
          CV = cv
        )
      }
incProgress(
  0.75,
  detail = "Calculating coefficient of variation"
)
      cv_all <- map_dfr(df_split, calc_cv_by_temp)
incProgress(
  0.95,
  detail = "Rendering plot"
)
      p <- ggplot(cv_all, aes(x = temperature, y = CV)) +
        geom_boxplot() +
        labs(
          title = "",
          x = "Temperature (°C)",
          y = ifelse(method == "log2", "Log2 CV", "CV")
        ) +
        theme_minimal(base_size = 14)


incProgress(
  1,
  detail = "Done"
)

p

    }
  )

})

output$download_cv_plot <- downloadHandler(

  filename = function() {
    paste0("CV_Plot_", Sys.Date(), ".png")
  },

  content = function(file) {

    p <- isolate({

      df <- plot_data()
      method <- selected_methods()

      df_long <- df %>%
        pivot_longer(
          cols = -1,
          names_to = "Replicate",
          values_to = "Intensity"
        ) %>%
        separate(
          Replicate,
          into = c("temperature", "Replicate"),
          sep = "_"
        )

      df_split <- df_long %>%
        group_by(temperature) %>%
        group_split()

      calc_cv_by_temp <- function(df) {

        wide <- df %>%
          select(1, Replicate, Intensity) %>%
          pivot_wider(
            names_from = Replicate,
            values_from = Intensity
          )

        id_col <- names(wide)[1]

        mat <- wide %>%
          select(-all_of(id_col))

        mat[mat == 0] <- NA

        mat <- as.matrix(mat)

        cv <- apply(mat, 1, function(x) {

          x_clean <- x[!is.na(x)]

          if (length(x_clean) >= 2) {

            (sd(x_clean) / mean(x_clean)) * 100

          } else {

            NA

          }

        })

        tibble(
          !!id_col := wide[[id_col]],
          temperature = unique(df$temperature),
          CV = cv
        )

      }

      cv_all <- purrr::map_dfr(
        df_split,
        calc_cv_by_temp
      )

      ggplot(
        cv_all,
        aes(
          x = temperature,
          y = CV
        )
      ) +
        geom_boxplot() +
        theme_minimal(base_size = 14)

    })

    ggsave(
      file,
      p,
      width = 10,
      height = 6,
      dpi = 300
    )

  }

)


    # Return reactives if needed


qc_applied <- reactive({

  !is.null(
    tryCatch(
      zero_stats_data$filtered(),
      error = function(e) NULL
    )
  )

})

return(
list(
plot_data = plot_data,
current_data = current_data,
raw_data = raw_data,
normalized_data = normalized_data,
ref_protein_name = ref_protein_name,
qc_applied = qc_applied,
selected_methods = reactive(
input$normalization_methods %||%
character(0)
),
lowest_temp_normalized = lowest_temp_normalized
)
)
  })
}
####################################################################################
#clustering modules
clusterPanelModuleUI <- function(id, label) {
  ns <- NS(id)

  tagList(

    # CONTROL PANEL (FULL WIDTH)
    fluidRow(
      box(
        title = tagList(
        "Clustering controls",
        tags$span(
          `data-toggle` = "tooltip",
          `data-placement` = "right",
          `data-html` = "true",
          title = HTML(
            "<b>Clustering analysis</b><br><br>

            Proteins are grouped according to the similarity of their thermal profiles across the temperature gradient.<br><br>

            Proteins displaying similar melting behaviour are assigned to the same cluster.<br><br>

            Clustering can be performed on a single condition or across multiple conditions using combined protein profiles."
          ),
          tags$i(
            class = "fa fa-question-circle",
            style = "margin-left:8px; color:#337ab7; cursor:pointer;"
          )
        )
      ),

        width = 12,
        status = "primary",
        solidHeader = TRUE,

        selectInput(
          inputId = ns("correlation_type"),
          label = tagList(
          "Select type of correlation:",
          tags$span(
            `data-toggle` = "tooltip",
            `data-placement` = "right",
            `data-html` = "true",
            title = HTML(
              "<b>Correlation modes</b><br><br>

              <b>Individual condition</b><br>

              Correlations are calculated using protein thermal profiles from a single condition only.<br><br>

              <b>Condition comparison</b><br>

              At least two conditions must be selected.<br><br>

              The thermal profiles of the same protein across all selected conditions are combined into a single profile and correlated together.<br><br>

              This allows clustering based on similarities and differences observed between conditions."
            ),
            tags$i(
              class = "fa fa-question-circle",
              style = "margin-left:8px; color:#337ab7; cursor:pointer;"
            )
          )
        ),
          choices = c(
            "Individual condition" = "individual",
            "Condition comparison" = "comparison"
          ),
          selected = "individual"
        ),

        conditionalPanel(
          condition = sprintf("input['%s'] == 'comparison'", ns("correlation_type")),

        uiOutput(ns("comparison_dataset_ui"))
        ),

        numericInput(
          inputId = ns("num_clusters"),
          label = "Select number of clusters (k):",
          value = 4,
          min = 2,
          max = 10,
          step = 1
        ),

      selectInput(
        inputId = ns("Correlation_method"),
        label = "Select correlation method:",
        choices = c(
          "Pearson" = "pearson",
          "Spearman" = "spearman"
        ),
        selected = "spearman"
      ),

        actionButton(ns("run_clustering"), "Run clustering"),

        br(),
        br(),
        tags$b(
          "Cluster recommendation ",
          tags$span(
            `data-toggle` = "tooltip",
            `data-placement` = "right",
            `data-html` = "true",
            title = HTML(
              "<b>Silhouette analysis</b><br><br>

              The recommended number of clusters is determined using silhouette analysis.<br><br>

              The silhouette score evaluates how similar each protein is to its own cluster compared to neighbouring clusters.<br><br>

              Higher average silhouette scores indicate stronger separation between clusters and generally represent a better clustering solution."
            ),
            tags$i(
              class = "fa fa-question-circle",
              style = "margin-left:5px; color:#337ab7; cursor:pointer;"
            )
          )
        ),
        br(),
        uiOutput(ns("cluster_recommendation")),
        plotOutput(
          ns("cluster_recommendation_plot"),
          height = "250px"
        ) %>% withSpinner(
color = "#0EA5A5"
)
      )
    ),

    # CLUSTERED DATA PLOT
    fluidRow(
      box(
        title = tagList(
        "Clustered data",
        tags$span(
          `data-toggle` = "tooltip",
          `data-placement` = "right",
          `data-html` = "true",
          title = HTML(
            "<b>Cluster colours</b><br><br>

            Cluster colours are shared between the cluster trend plot and the cluster heatmap.<br><br>"
          ),
          tags$i(
            class = "fa fa-question-circle",
            style = "margin-left:8px; color:#337ab7; cursor:pointer;"
          )
        )
      ),
        width = 12,
        status = "primary",
        solidHeader = TRUE,
        plotOutput(ns("clustering_plot"), height = "500px") %>% withSpinner(
color = "#0EA5A5"
),

        checkboxInput(
          inputId = ns("show_missing_lines"),
          label = "Show missing values as zero",
          value = FALSE
        ),


        downloadButton(ns("download_clustering_plot"), "Download cluster trend plot")
      )
    ),

    # HEATMAP (BELOW)
    fluidRow(
      box(
        title = "Heatmap with clusters",
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        div(
          style = "overflow-y: auto; max-height: 700px;",
          plotOutput(ns("clustering_heatmap"), height = "600px") %>% withSpinner(
color = "#0EA5A5"
)
        ),

        downloadButton(ns("download_clustering_heatmap"), "Download heatmap")
      )
    ),

    # RECLUSTER SECTION (UNCHANGED STRUCTURE)
fluidRow(
  box(
    title = "Sub-cluster data",
    width = 12,
    status = "primary",
    solidHeader = TRUE,
    collapsible = TRUE,
    collapsed = TRUE,

    # OPTIONS ROW (unchanged)
    fluidRow(
      box(
        title = tagList(
        "Sub-clustering options",
        tags$span(
          `data-toggle` = "tooltip",
          `data-placement` = "right",
          `data-html` = "true",
          title = HTML(
            "<b>Sub-clustering</b><br><br>

            Sub-clustering performs a second clustering analysis within an existing cluster.<br><br>

            This can reveal smaller groups of proteins with more specific melting behaviours that may not be visible in the initial clustering result."
          ),
          tags$i(
            class = "fa fa-question-circle",
            style = "margin-left:8px; color:#337ab7; cursor:pointer;"
          )
        )
      ),
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        numericInput(
          inputId = ns("recluster_num"),
          label = "Select cluster to sub-cluster (1 to k):",
          value = 1,
          min = 1,
          max = 10,
          step = 1
        ),

        numericInput(
          inputId = ns("recluster_k"),
          label = "Select number of clusters for sub-clustering:",
          value = 4,
          min = 2,
          max = 10,
          step = 1
        ),

        actionButton(
        ns("run_reclustering"),
        "Run Sub-lustering"
      ),

      br(),
      br(),
      tags$b(
        "Sub-cluster recommendation ",
        tags$span(
          `data-toggle` = "tooltip",
          `data-placement` = "right",
          `data-html` = "true",
          title = HTML(
            "<b>Silhouette analysis</b><br><br>

            The recommended number of sub-clusters is determined using silhouette analysis.<br><br>

            Higher silhouette values indicate that proteins are well matched to their assigned sub-cluster and poorly matched to neighbouring sub-clusters."
          ),
          tags$i(
            class = "fa fa-question-circle",
            style = "margin-left:5px; color:#337ab7; cursor:pointer;"
          )
        )
      ),
      br(),
      uiOutput(ns("recluster_recommendation")),

      plotOutput(
        ns("recluster_recommendation_plot"),
        height = "250px"
      ) %>% withSpinner(
color = "#0EA5A5"
)
      )
    ),

    #  RECLUSTERED PLOT (FULL WIDTH)
    fluidRow(
      box(
        title = tagList(
        "Sub-clustered plot",
        tags$span(
          `data-toggle` = "tooltip",
          `data-placement` = "right",
          `data-html` = "true",
          title = HTML(
            "<b>Sub-cluster colours</b><br><br>

            Sub-cluster colours are shared between the Sub-cluster Trend Plot and the Sub-cluster Heatmap.<br><br>

            A given colour always represents the same sub-cluster across both visualizations."
          ),
          tags$i(
            class = "fa fa-question-circle",
            style = "margin-left:8px; color:#337ab7; cursor:pointer;"
          )
        )
      ),
        width = 12,
        status = "success",
        solidHeader = TRUE,

        plotOutput(ns("reclustered_plot"), height = "500px") %>% withSpinner(
color = "#0EA5A5"
),
        checkboxInput(
          inputId = ns("show_missing_recluster"),
          label = "Show missing values as zero (sub-clustering plot)",
          value = FALSE
        ),
        downloadButton(
          ns("download_reclustered_plot"),
          "Download sub-ustered trend plot"
        )
      )
    ),

    #  RECLUSTERED HEATMAP (FULL WIDTH BELOW)
    fluidRow(
      box(
        title = "Sub-clustered heatmap",
        width = 12,
        status = "info",
        solidHeader = TRUE,

        div(
          style = "overflow-y: auto; max-height: 700px;",
          plotOutput(ns("reclustered_heatmap"), height = "600px") %>% withSpinner(
color = "#0EA5A5"
)
        ),

        downloadButton(
          ns("download_reclustered_heatmap"),
          "Download sub-clustered heatmap"
        )
      )
    )
  )
)
  )
}

#server
clusterPanelModuleServer <- function(id, plot_data, dataset_names, normalization_results) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    current_dataset_name <- sub(
  "_[0-9]+$",
  "",
  id
)
    print(current_dataset_name)
      #read dataset names and update the selectInput choices

output$comparison_dataset_ui <- renderUI({

  req(input$correlation_type == "comparison")

  choices <- dataset_names()

  req(length(choices) > 0)

  selectizeInput(
    inputId = ns("comparison_dataset"),
    label = "Select datasets for comparison:",
    choices = choices,
    selected = choices,
    multiple = TRUE,
    options = list(placeholder = "Select one or more datasets")
  )
})

#numb of clusters recommendation

cluster_recommendation <- reactive({

  withProgress(
    message = "Assessing cluster quality",
    value = 0,
    {

      incProgress(
        0.10,
        detail = "Preparing clustered matrix"
      )

      req(clustering_result())

      method <- input$Correlation_method

      mat <- clustering_result()$mat
incProgress(
  0.30,
  detail = "Calculating correlation matrix"
)
  cor_mat <- cor(
    t(mat),
    method = method
  )

  cor_mat[is.na(cor_mat)] <- 0

  diss_mat <- as.dist(
    (1 - cor_mat) / 2
  )
incProgress(
  0.50,
  detail = "Generating clustering tree"
)
  hc <- hclust(
    diss_mat,
    method = "complete"
  )

  max_k <- min(
    10,
    nrow(mat) - 1
  )

  req(max_k >= 2)

  k_values <- 2:max_k
incProgress(
  0.70,
  detail = "Calculating silhouette scores"
)
  scores <- sapply(
    k_values,
    function(k) {

      cl <- cutree(
        hc,
        k = k
      )

      sil <- cluster::silhouette(
        cl,
        diss_mat
      )

      mean(
        sil[, "sil_width"],
        na.rm = TRUE
      )

    }
  )
incProgress(
  1,
  detail = "Done"
)
  list(
    best_k = k_values[which.max(scores)],
    scores = scores,
    k_values = k_values
  )

    }
  )

})

output$cluster_recommendation <- renderUI({

  if (is.null(tryCatch(clustering_result(), error = function(e) NULL))) {

    return(
      div(
        style = "
          padding:15px;
          background-color:#f8f9fa;
          border:1px solid #dee2e6;
          border-radius:4px;
          text-align:center;
          color:#6c757d;
        ",

        tags$h4(
          icon("project-diagram"),
          " Clustering recommendation"
        ),

        tags$p(
          "Select a clustering mode and click ",
          tags$b("Run Clustering"),
          " to generate the optimal number of clusters recommendation."
        )
      )
    )
  }

  rec <- cluster_recommendation()

  div(

    style = "
      padding:10px;
      background-color:#f8f9fa;
      border:1px solid #dee2e6;
      border-radius:4px;
    ",

    tags$b(
      paste(
        "Recommended number of clusters:",
        rec$best_k
      )
    ),

    tags$br(),

    paste(
      "Average silhouette score:",
      round(
        max(rec$scores),
        3
      )
    )

  )

})
output$cluster_recommendation_plot <- renderPlot({

  withProgress(
    message = "Evaluating cluster quality",
    value = 0,
    {


  if (is.null(tryCatch(clustering_result(), error = function(e) NULL))) {

    plot.new()

    text(
      0.5,
      0.5,
      "Run clustering to generate\ncluster quality assessment",
      cex = 1.4,
      col = "grey50"
    )

    return()
  }
incProgress(
  0.50,
  detail = "Calculating silhouette scores"
)
  rec <- cluster_recommendation()
incProgress(
  0.80,
  detail = "Building recommendation plot"
)
  df_plot <- data.frame(
    k = rec$k_values,
    silhouette = rec$scores
  )

  p<- ggplot(
    df_plot,
    aes(
      x = k,
      y = silhouette
    )
  ) +

    geom_line(
      color = "steelblue4",
      linewidth = 1
    ) +

    geom_point(
      color = "steelblue4",
      size = 3
    ) +

    geom_vline(
      xintercept = rec$best_k,
      linetype = "dashed",
      color = "red"
    ) +

    annotate(
      "text",
      x = rec$best_k,
      y = max(rec$scores) * 1.05,
      label = paste(
        "Recommended k =",
        rec$best_k
      ),
      color = "red",
      hjust = 0,
      size = 4
    )+

scale_x_continuous(
  breaks = rec$k_values
) +

expand_limits(
  y = max(rec$scores) * 1.15
) +

    labs(
      x = "Number of clusters (k)",
      y = "Average silhouette score",
      title = "Cluster quality assessment"
    ) +

    theme_minimal() +
theme(
  plot.margin = margin(
    t = 10,
    r = 25,
    b = 10,
    l = 25
  )
)
incProgress(
  1,
  detail = "Done"
)

p

    }
  )

})


# Clustering logic
clustering_result <- eventReactive(input$run_clustering, {

  req(input$num_clusters, input$Correlation_method, input$correlation_type)

  k <- input$num_clusters
  method <- input$Correlation_method

  withProgress(
    message = "Running clustering",
    value = 0,
    {
      incProgress(
  0.05,
  detail = "Preparing dataset"
)


  # =====================================================
  # INDIVIDUAL MODE (UNCHANGED)
  # =====================================================
  if (input$correlation_type == "individual") {
incProgress(
  0.20,
  detail = "Loading normalized data"
)
    req(plot_data())

    
df <- plot_data()



df_clean <- df %>%
  mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))



mat <- as.matrix(df_clean[, -1])
rownames(mat) <- df_clean[[1]]


incProgress(
  0.50,
  detail = "Calculating correlation matrix"
)
cor_mat <- cor(t(mat), method = method)

    cor_mat[is.na(cor_mat)] <- 0

    diss_mat <- as.dist((1 - cor_mat) / 2)
    incProgress(
  0.75,
  detail = "Generating hierarchical clustering"
)
    hc <- hclust(diss_mat, method = "complete")
    clusters <- cutree(hc, k = k)
incProgress(
  1,
  detail = "Done"
)
    return(list(
      mode = "individual",
      mat = mat,
      hc = hc,
      clusters = clusters,
      df = df
    ))
  }

  # =====================================================
  # COMPARISON MODE (FIXED PIPELINE)
  # =====================================================
  req(input$comparison_dataset)
  req(length(input$comparison_dataset) >= 1)
incProgress(
  0.20,
  detail = "Loading comparison datasets"
)
  df_main <- plot_data()
  selected_names <- input$comparison_dataset
selected_names <- setdiff(
  selected_names,
  current_dataset_name
)
  req(all(selected_names %in% names(normalization_results)))

  # -----------------------------
  # 1. COLLECT DATASETS
  # -----------------------------
  datasets_list <- c(
    setNames(list(df_main), current_dataset_name),
    setNames(
      lapply(selected_names, function(name) {
        normalization_results[[name]]$plot_data()
      }),
      selected_names
    )
  )

  # -----------------------------
  # 2. RENAME COLUMNS (add dataset suffix)
  # -----------------------------
  datasets_list <- lapply(names(datasets_list), function(nm) {

    df <- datasets_list[[nm]]

    meas_cols <- names(df)[-1]

    names(df) <- c(
      names(df)[1],
      paste0(meas_cols, ".", nm)
    )

    df
  })

  names(datasets_list) <- c(current_dataset_name, selected_names)

  # -----------------------------
  # 3. NORMALIZATION (PER DATASET)
  # -----------------------------


datasets_norm <- datasets_list
incProgress(
  0.50,
  detail = "Clustering each dataset"
)
  # -----------------------------
  # 4. CLUSTER EACH DATASET INDEPENDENTLY
  # -----------------------------
  cluster_results <- lapply(datasets_norm, function(df) {

    df_clean <- df %>%
      mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))

    mat <- as.matrix(df_clean[, -1])
    rownames(mat) <- df_clean[[1]]



    cor_matrix <- cor(t(mat), method = method)
    cor_matrix[is.na(cor_matrix)] <- 0

    diss <- as.dist((1 - cor_matrix) / 2)
    hc <- hclust(diss, method = "complete")
    clusters <- cutree(hc, k = k)

    list(
      df = df,
      clusters = clusters,
      hc = hc
    )
  })

  # -----------------------------
  # 5. USE CURRENT DATASET AS REFERENCE
  # -----------------------------
  current_res <- cluster_results[[current_dataset_name]]
  clusters <- current_res$clusters

  # -----------------------------
  # 6. APPLY SAME CLUSTERS TO ALL DATASETS
  # -----------------------------
  clustered_data <- lapply(datasets_norm, function(df) {

    annotation <- data.frame(
      Accession = names(clusters),
      Cluster   = factor(clusters)
    )

    df %>% left_join(annotation, by = "Accession")
  })

  # -----------------------------
  # 7. JOIN FOR PLOTTING ONLY
  # -----------------------------
# keep cluster from reference dataset ONLY
cluster_ref <- clustered_data[[current_dataset_name]] %>%
  dplyr::select(Accession, Cluster)

# join all datasets WITHOUT cluster
data_only <- lapply(clustered_data, function(df) {
  dplyr::select(df, -Cluster)
})

df_joined <- Reduce(function(x, y) {
  dplyr::left_join(x, y, by = "Accession")
}, data_only)

# reattach cluster
df_joined <- dplyr::left_join(cluster_ref, df_joined, by = "Accession")

# ✅ REBUILD MATRIX FOR HEATMAP

incProgress(
  0.85,
  detail = "Building combined clustering structure"
)
df_clean <- df_joined %>%
  dplyr::mutate(across(where(is.numeric), ~ tidyr::replace_na(., 0)))

mat <- df_clean %>%
  dplyr::select(-Accession, -Cluster) %>%
  as.matrix()

rownames(mat) <- df_clean$Accession

# ✅ IMPORTANT: DO NOT TOUCH CLUSTERS LENGTH
# just reuse original clustering order
cor_matrix <- cor(t(mat), method = method)
cor_matrix[is.na(cor_matrix)] <- 0

diss <- as.dist((1 - cor_matrix) / 2)
hc <- hclust(diss, method = "complete")





cor_matrix <- cor(t(mat), method = method)
cor_matrix[is.na(cor_matrix)] <- 0

diss <- as.dist((1 - cor_matrix) / 2)
hc <- hclust(diss, method = "complete")
incProgress(
  1,
  detail = "Done"
)
return(list(
  mode = "comparison",
  mat = mat,
  hc = hc,
  clusters = clusters,
  df = df_joined
))

    }
  )

})



    # Heatmap

# =====================================================================
# Helper function that builds the heatmap (for both renderPlot & download)
# =====================================================================
make_clustering_heatmap <- function() {

  res_obj <- clustering_result()
  incProgress(
  0.10,
  detail = "Loading clustered matrix"
)

  req(res_obj)

  mat <- res_obj$mat
  hc  <- res_obj$hc
  clusters <- res_obj$clusters

  # APPLY ROW ORDER FROM HCLUST (CRITICAL)
  mat <- mat[hc$order, , drop = FALSE]

  # ORDER COLUMNS EXACTLY AS BEFORE


col_info <- tibble::tibble(
  name = colnames(mat)
)


  col_info$temperature <- as.numeric(stringr::str_extract(col_info$name, "^[0-9.]+"))
  col_info$replicate   <- as.integer(stringr::str_extract(col_info$name, "(?<=REP)\\d+"))

  #  extract dataset name
col_info$dataset <- stringr::str_remove(
  col_info$name,
  "^\\d+(?:\\.\\d+)?_REP\\d+\\."
)

col_info$dataset[col_info$dataset == col_info$name] <-
  current_dataset_name


  col_info$dataset <- factor(
    col_info$dataset,
    levels = c(
      current_dataset_name,
      setdiff(unique(col_info$dataset), current_dataset_name)
    )
  )

  incProgress(
  0.30,
  detail = "Ordering datasets and samples"
)
  ordered_cols <- col_info %>%
    dplyr::arrange(dataset, temperature, replicate) %>%
    dplyr::pull(name)
dataset_positions <- col_info %>%
  dplyr::arrange(dataset, temperature, replicate) %>%
  dplyr::mutate(pos = seq_len(n())) %>%
  dplyr::group_by(dataset) %>%
  dplyr::summarise(
    midpoint = mean(pos),
    .groups = "drop"
  )
# get cluster per protein
cluster_vec <- clusters[rownames(mat)]

# order: cluster FIRST, then dendrogram order
order_idx <- order(cluster_vec, hc$order)

mat <- mat[order_idx, , drop = FALSE]


  # COLOR PALETTE (SAME AS PHEATMAP DEFAULT)
  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(7, "RdYlBu")))(100)

  # GENERATE BREAKS (LIKE PHEATMAP)
  breaks <- seq(min(mat, na.rm = TRUE),
                max(mat, na.rm = TRUE),
                length.out = length(colors) + 1)
  midpoints <- (breaks[-1] + breaks[-length(breaks)]) / 2



  # LONG FORMAT
incProgress(
  0.50,
  detail = "Preparing heatmap matrix"
)
df_long <- as.data.frame(mat)
df_long$Protein <- rownames(df_long)


df_long <- tidyr::pivot_longer(
  df_long,
  cols = -Protein,
  names_to = "Sample",
  values_to = "Intensity"
)

df_long$Intensity <- df_long$Intensity




  # FACTORS FOR ORDER
  
df_long$Protein <- factor(df_long$Protein, levels = rownames(mat))
df_long$y_index <- as.numeric(df_long$Protein)

  df_long$Sample  <- factor(df_long$Sample, levels = ordered_cols)
sample_labels <- stringr::str_extract(
  ordered_cols,
  "^\\d+(?:\\.\\d+)?_REP\\d+"
)

  #  Build cluster annotation bar (left side)

cluster_levels <- sort(unique(clusters))

incProgress(
  0.70,
  detail = "Building cluster annotations"
)
cluster_df <- data.frame(
  Protein = factor(rownames(mat), levels = rownames(mat)),
  Cluster = factor(clusters[rownames(mat)], levels = sort(unique(clusters))),
  x = "ClusterBar"
)



cluster_df$Protein <- factor(cluster_df$Protein, levels = rownames(mat))
cluster_df$y_index <- as.numeric(cluster_df$Protein)
dataset_label_df <- data.frame(
  x = dataset_positions$midpoint,
  y = max(cluster_df$y_index) + 90,
  label = as.character(dataset_positions$dataset)
)

incProgress(
  0.90,
  detail = "Rendering heatmap"
)
p <- ggplot() +
geom_text(
  data = dataset_label_df,
  aes(
    x = x,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  fontface = "bold",
  size = 5
) +
  geom_tile(
    data = df_long,
    aes(x = Sample, y = y_index, fill = Intensity)
  ) +


scale_fill_gradientn(
  colours = colors,
  breaks = pretty(range(df_long$Intensity, na.rm = TRUE), n = 5),
  name = "Intensity",
  limits = range(df_long$Intensity, na.rm = TRUE)
) +

  ggnewscale::new_scale_fill() +

  geom_tile(
    data = cluster_df,
    aes(x = " ", y = y_index, fill = Cluster)
  ) +

  scale_fill_brewer(palette = "Set2", name = "Cluster") +

  theme_minimal() +
scale_x_discrete(
  labels = function(x) {
    ifelse(
      x == " ",
      "",
      stringr::str_extract(x, "^\\d+(?:\\.\\d+)?_REP\\d+")
    )
  }
)+
scale_y_continuous(
  expand = expansion(mult = c(0, 0))
)+
coord_cartesian(
  clip = "off"
) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0
    ),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    axis.title.x = element_blank(),
    plot.margin = margin(
  t = 20,
  r = 5,
  b = 5,
  l = 5
)
  )


incProgress(
  1,
  detail = "Done"
)

  return(p)
}


# =====================================================================
# Original plot output
# =====================================================================

output$clustering_heatmap <- renderPlot({

  withProgress(
    message = "Generating clustering heatmap",
    value = 0,
    {

      print(
        make_clustering_heatmap()
      )

    }
  )

})




# =====================================================================
# Download handler for the heatmap
# =====================================================================
output$download_clustering_heatmap <- downloadHandler(
  filename = function() {
    paste0("clustering_heatmap_", Sys.Date(), ".png")
  },
  content = function(file) {
    # Save heatmap directly to file using pheatmap's filename argument
    make_clustering_heatmap(filename = file)
  }
)


    # Clustered data
dataclustered <- reactive({
  res <- clustering_result()
  req(res)

  df <- res$df
  clusters <- res$clusters

  # ✅ SAFE ALIGNMENT (NO FILTERING REQUIRED)
  df$Cluster <- factor(clusters[df$Accession])

  df
})

    # Cluster counts
    output$cluster_counts <- renderTable({
      req(dataclustered())
      df <- dataclustered()
      cluster_counts <- as.data.frame(table(df$Cluster))
      colnames(cluster_counts) <- c("Cluster", "Number of Proteins")
      cluster_counts
    })

    # Clustered plot


# =====================================================================
# Helper function to build the clustering trend plot
# =====================================================================
make_clustering_plot <- function() {

  res <- clustering_result()
  req(res)

  df <- dataclustered()
incProgress(
  0.10,
  detail = "Loading clustered proteins"
)
  # =====================================================
  # CASE 1 — INDIVIDUAL CONDITION 
  # =====================================================
  if (res$mode == "individual") {
incProgress(
  0.25,
  detail = "Preparing cluster annotations"
)
    df$Cluster <- as.factor(df$Cluster)

    counts <- df %>%
      count(Cluster, name = "n")
cluster_levels <- sort(unique(df$Cluster))

cluster_colors <- setNames(
  RColorBrewer::brewer.pal(
    max(3, length(cluster_levels)),
    "Set2"
  )[seq_along(cluster_levels)],
  as.character(cluster_levels)
)
labeller_fn <- labeller(
  Cluster = function(x) {
    paste0(
      x,
      " (n=",
      counts$n[match(x, counts$Cluster)],
      ")"
    )
  }
)
incProgress(
  0.45,
  detail = "Reshaping protein profiles"
)
    # ---- reshape ----
    df_long <- df %>%
      pivot_longer(
        cols = -c(1, "Cluster"),
        names_to   = c("temperature", "Replicate", "Dataset"),
        names_pattern = "^(\\d+(?:\\.\\d+)?)_(REP\\d+)(?:\\.(.+))?$",
        values_to  = "Intensity",
        values_drop_na = FALSE
      ) %>%
      mutate(
        temperature = as.numeric(temperature),

        # ✅ MISSING VALUE CONTROL (KEY FIX)
        Intensity = if (isTRUE(input$show_missing_lines)) {
          Intensity
        } else {
          ifelse(Intensity == 0, NA, Intensity)
        }
      )
incProgress(
  0.65,
  detail = "Calculating protein trends"
)
    # ---- protein-level mean ----
    df2 <- df_long %>%
      group_by(across(1), Cluster, temperature) %>%
      summarise(
        Mean_Intensity = mean(Intensity, na.rm = TRUE),
        .groups = "drop"
      )

    # ---- cluster mean ----
    avg <- df_long %>%
      group_by(Cluster, temperature) %>%
      summarise(
        Mean_Intensity = mean(Intensity, na.rm = TRUE),
        .groups = "drop"
      )


facet_color_df <- data.frame(
  Cluster = factor(cluster_levels, levels = cluster_levels),
  x = min(df2$temperature, na.rm = TRUE),
  y = max(df2$Mean_Intensity, na.rm = TRUE),
  label = paste0(
    "Cluster ",
    cluster_levels,
    "\n(n=",
    counts$n[match(cluster_levels, counts$Cluster)],
    ")"
  )
)
cluster_levels <- sort(unique(df$Cluster))

cluster_colors <- setNames(
  RColorBrewer::brewer.pal(
    max(3, length(cluster_levels)),
    "Set2"
  )[seq_along(cluster_levels)],
  as.character(cluster_levels)
)

incProgress(
  0.90,
  detail = "Rendering cluster plot"
)
p <- ggplot(df2, aes(

  x = temperature,
  y = Mean_Intensity,
  group = !!sym(names(df2)[1])
)) +

geom_line(color = "steelblue2", alpha = 0.1) +

geom_line(
  data = avg,
  aes(group = 1),
  color = "steelblue4",
  linewidth = 1.5
) +



geom_point(
  data = facet_color_df,
  aes(
    x = x,
    y = y,
    colour = Cluster
  ),
  inherit.aes = FALSE,
  shape = 15,
  size = 5
) +
geom_text(
  data = facet_color_df,
  aes(
    x = x,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  hjust = -0.2,
  vjust = 1,
  size = 4,
  fontface = "bold",
  colour = "black"
) +

scale_colour_manual(
  values = cluster_colors,
  guide = "none"
) +

facet_wrap(~ Cluster, nrow = 1) +
      theme_light(base_size = 14) +

      theme(
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),

        panel.border = element_rect(
          color = "grey80",
          fill = NA
        ),

        strip.background = element_rect(
          fill = "white",
          color = NA
        ),

strip.text = element_text(
  face = "bold",
  size = 12
)
      ) +
      labs(
        title = "",
        x = "Temperature (°C)",
        y = "Relative abundance"
      )
incProgress(
  1,
  detail = "Done"
)
return(p)
  }
incProgress(
  0.35,
  detail = "Preparing comparison datasets"
)
  # =====================================================
  # CASE 2 — CONDITION COMPARISON 
  # =====================================================

  df_long <- df %>%
    pivot_longer(
      cols = -c(1, "Cluster"),
      names_to   = c("temperature", "Replicate", "Dataset"),
      names_pattern = "^(\\d+(?:\\.\\d+)?)_(REP\\d+)(?:\\.(.+))?$",
      values_to  = "Intensity",
      values_drop_na = FALSE
    ) %>%
    mutate(
      temperature = as.numeric(temperature),

      Dataset = ifelse(is.na(Dataset), current_dataset_name, Dataset),

      Dataset = factor(Dataset, levels = unique(Dataset)),

      # ✅ SAME FIX HERE (important!)
      Intensity = if (isTRUE(input$show_missing_lines)) {
        Intensity
      } else {
        ifelse(Intensity == 0, NA, Intensity)
      }
    )
incProgress(
  0.60,
  detail = "Calculating protein trends"
)
  # ---- protein curves ----
  df_protein <- df_long %>%
    group_by(across(1), Cluster, Dataset, temperature) %>%
    summarise(
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      .groups = "drop"
    )

  # ---- cluster trends ----
  avg <- df_long %>%
    group_by(Cluster, Dataset, temperature) %>%
    summarise(
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      .groups = "drop"
    )

  # ---- cluster sizes ----
  counts <- df %>%
    count(Cluster, name = "n")
cluster_levels <- sort(unique(df$Cluster))

cluster_colors <- setNames(
  RColorBrewer::brewer.pal(
    max(3, length(cluster_levels)),
    "Set2"
  )[seq_along(cluster_levels)],
  as.character(cluster_levels)
)
facet_color_df <- data.frame(
  Cluster = factor(cluster_levels, levels = cluster_levels),
  x = min(df_protein$temperature, na.rm = TRUE),
  y = max(df_protein$Mean_Intensity, na.rm = TRUE)
)
  labeller_fn <- labeller(
    Cluster = function(x) {
      paste0(
        x, " (n=",
        counts$n[match(x, counts$Cluster)],
        ")"
      )
    }
  )

  datasets_unique <- unique(df_long$Dataset)

  curve_colors <- setNames(
    scales::hue_pal()(length(datasets_unique)),
    datasets_unique
  )

  trend_colors <- setNames(
    scales::hue_pal(l = 40)(length(datasets_unique)),
    datasets_unique
  )
incProgress(
  0.90,
  detail = "Rendering comparison plot"
)
  p <- ggplot(
    df_protein,
    aes(
      x = temperature,
      y = Mean_Intensity,
      group = interaction(Dataset, !!sym(names(df)[1])),
      color = Dataset
    )
  ) +
    geom_line(alpha = 0.08, linewidth = 0.4) +

    geom_line(
      data = avg,
      aes(
        x = temperature,
        y = Mean_Intensity,
        group = Dataset,
        color = Dataset
      ),
      linewidth = 2.2
    ) +

    scale_color_manual(values = curve_colors, guide = "none") +
    ggnewscale::new_scale_color() +

    geom_line(
      data = avg,
      aes(
        x = temperature,
        y = Mean_Intensity,
        group = Dataset,
        color = Dataset
      ),
      linewidth = 1
    ) +

    scale_color_manual(
      values = trend_colors,
      name = "Dataset"
    ) +
    ggnewscale::new_scale_colour() +

geom_point(
  data = facet_color_df,
  aes(
    x = x,
    y = y,
    colour = Cluster
  ),
  inherit.aes = FALSE,
  shape = 15,
  size = 5
) +

scale_colour_manual(
  values = cluster_colors,
  guide = "none"
) +
    facet_wrap(~ Cluster, nrow = 1, labeller = labeller_fn) +

    labs(
      title = "Cluster comparison across datasets",
      x = "Temperature (°C)",
      y = "Relative abundance"
    ) +

    theme_light(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(
        color = "grey80",
        fill = NA
      )
    )
incProgress(
  1,
  detail = "Done"
)

return(p)
}



# =====================================================================
# Render plot in UI
# =====================================================================
output$clustering_plot <- renderPlot({

  withProgress(
    message = "Generating cluster trend plot",
    value = 0,
    {

      p <- make_clustering_plot()

      p

    }
  )

})


# =====================================================================
# Download handler
# =====================================================================
output$download_clustering_plot <- downloadHandler(
  filename = function() {
    paste0("clustering_plot_", Sys.Date(), ".png")
  },
  content = function(file) {
    p <- make_clustering_plot()
    ggsave(
      filename = file,
      plot = p,
      width = 10,
      height = 5,
      dpi = 300
    )
  }
)
observe({

  cat(
    "\nRECLUSTER INPUT CHANGED:",
    input$recluster_num,
    "\n"
  )

})

recluster_recommendation <- eventReactive(

  list(
    input$recluster_num,
    input$Correlation_method,
    clustering_result()
  ),

  {

withProgress(
  message = "Assessing recluster quality",
  value = 0,
  {

    incProgress(
      0.10,
      detail = "Preparing cluster subset"
    )

    res <- clustering_result()

    df <- res$df

    clusters <- res$clusters

    df$Cluster <- factor(
      clusters[df$Accession]
    )

    target_cluster <- input$recluster_num

    df_subset <- df %>%
      dplyr::filter(
        Cluster == target_cluster
      ) %>%
      dplyr::select(-Cluster)



    req(nrow(df_subset) >= 5)

    mat <- as.matrix(
      df_subset[, -1]
    )

    rownames(mat) <- df_subset[[1]]
incProgress(
  0.30,
  detail = "Calculating correlation matrix"
)
    cor_mat <- cor(
      t(mat),
      method = input$Correlation_method
    )

    cor_mat[is.na(cor_mat)] <- 0

    diss_mat <- as.dist(
      (1 - cor_mat) / 2
    )
incProgress(
  0.50,
  detail = "Generating clustering tree"
)
    hc <- hclust(
      diss_mat,
      method = "complete"
    )

    max_k <- min(
      10,
      nrow(df_subset) - 1
    )

    k_values <- 2:max_k
incProgress(
  0.70,
  detail = "Calculating silhouette scores"
)
    scores <- sapply(
      k_values,
      function(k) {

        cl <- cutree(
          hc,
          k = k
        )

        sil <- cluster::silhouette(
          cl,
          diss_mat
        )

        mean(
          sil[, "sil_width"],
          na.rm = TRUE
        )

      }
    )
incProgress(
  1,
  detail = "Done"
)
    list(
      best_k = k_values[which.max(scores)],
      scores = scores,
      k_values = k_values
    )

    }

  )

  }

)

output$recluster_recommendation <- renderUI({

  if (is.null(tryCatch(clustering_result(), error = function(e) NULL))) {

    return(
      div(
        style = "
          padding:10px;
          background-color:#f8f9fa;
          border:1px solid #dee2e6;
          border-radius:4px;
        ",
        "Run clustering first."
      )
    )

  }

  rec <- recluster_recommendation()

  div(

    style = "
      padding:10px;
      background-color:#f8f9fa;
      border:1px solid #dee2e6;
      border-radius:4px;
    ",

    tags$b(
      paste(
        "Recommended number of clusters:",
        rec$best_k
      )
    ),

    tags$br(),

    paste(
      "Average silhouette score:",
      round(
        max(rec$scores),
        3
      )
    )

  )

})

output$recluster_recommendation_plot <- renderPlot({

  if (is.null(tryCatch(clustering_result(), error = function(e) NULL))) {
    return(NULL)
  }

  rec <- recluster_recommendation()

  df_plot <- data.frame(
    k = rec$k_values,
    silhouette = rec$scores
  )

  ggplot(
    df_plot,
    aes(
      x = k,
      y = silhouette
    )
  ) +

    geom_line(
      color = "steelblue4",
      linewidth = 1
    ) +

    geom_point(
      color = "steelblue4",
      size = 3
    ) +

    geom_vline(
      xintercept = rec$best_k,
      linetype = "dashed",
      color = "red"
    ) +

    annotate(
      "text",
      x = rec$best_k - 0.2,
      y = max(rec$scores) * 1.05,

      label = paste(
        "Recommended k =",
        rec$best_k
      ),

      color = "red",
      hjust = 1,
      size = 4
    ) +

    scale_x_continuous(
      breaks = rec$k_values
    ) +

    expand_limits(
      y = max(rec$scores) * 1.15
    ) +

    labs(
      x = "Number of clusters (k)",
      y = "Average silhouette score",
      title = "Recluster quality assessment"
    ) +

    theme_minimal() +

    theme(
      plot.margin = margin(
        t = 10,
        r = 25,
        b = 10,
        l = 25
      )
    )

})

    # Reclustering logic
reclustering_result <- eventReactive(input$run_reclustering, {

  req(dataclustered(), input$recluster_num, input$recluster_k, input$Correlation_method)

  res_main <- clustering_result()     
  df <- dataclustered()

  recluster_num <- input$recluster_num
  recluster_k <- input$recluster_k
  method <- input$Correlation_method

  if (is.null(recluster_num) || recluster_num < 1 || recluster_num > input$num_clusters) return(NULL)
  if (is.null(recluster_k) || recluster_k < 2) return(NULL)

  # subset cluster
  df_subset <- df %>%
    dplyr::filter(Cluster == recluster_num) %>%
    dplyr::select(-Cluster)

  if (nrow(df_subset) == 0 || nrow(df_subset) < recluster_k) return(NULL)

  # matrix creation (works for both modes)
  mat <- as.matrix(df_subset[, -1])
  rownames(mat) <- df_subset[[1]]

  # correlation (same logic for both modes)
  cor_mat <- cor(t(mat), method = method)
  cor_mat[is.na(cor_mat)] <- 0

  diss_mat <- as.dist((1 - cor_mat) / 2)
  hc <- hclust(diss_mat, method = "complete")
  clusters <- cutree(hc, k = recluster_k)

  # IMPORTANT: keep "mode"
  list(
    mode = res_main$mode,   
    mat = mat,
    hc = hc,
    clusters = clusters,
    df_subset = df_subset
  )
})

    # Heatmap

# =====================================================================
# Helper function for reclustered heatmap
# =====================================================================
make_reclustered_heatmap <- function(res) {


  validate(
    need(!is.null(res), "No data to recluster. Try a different cluster or reduce recluster_k.")
  )

  mat <- res$mat
  hc <- res$hc
  clusters <- res$clusters
incProgress(
  0.10,
  detail = "Loading reclustered matrix"
)
  # ✅ ORDER ROWS: cluster first, then dendrogram
  cluster_vec <- clusters[rownames(mat)]
  order_idx <- order(cluster_vec, hc$order)
  mat <- mat[order_idx, , drop = FALSE]

  # ============================
  # COLUMN ORDERING
  # ============================

  col_info <- tibble::tibble(
    name = colnames(mat)
  )

  col_info$temperature <- as.numeric(stringr::str_extract(col_info$name, "^[0-9.]+"))
  col_info$replicate   <- as.integer(stringr::str_extract(col_info$name, "(?<=REP)\\d+"))

col_info$dataset <- stringr::str_remove(
  col_info$name,
  "^\\d+(?:\\.\\d+)?_REP\\d+\\."
)

# individual mode detection
col_info$dataset[
  col_info$dataset == col_info$name
] <- current_dataset_name



  col_info$dataset <- factor(
    col_info$dataset,
    levels = c(
      current_dataset_name,
      setdiff(unique(col_info$dataset), current_dataset_name)
    )
  )
incProgress(
  0.30,
  detail = "Ordering datasets and samples"
)
  ordered_cols <- col_info %>%
    dplyr::arrange(dataset, temperature, replicate) %>%
    dplyr::pull(name)

  mat <- mat[, ordered_cols, drop = FALSE]

  # ✅ FULL X-AXIS (Cluster + Samples)
  all_cols <- c("Cluster", ordered_cols)

  # ============================
  # LONG FORMAT
  # ============================
incProgress(
  0.50,
  detail = "Preparing heatmap matrix"
)

  df_long <- as.data.frame(mat)
  df_long$Protein <- rownames(mat)

  df_long <- tidyr::pivot_longer(
    df_long,
    cols = -Protein,
    names_to = "Sample",
    values_to = "Intensity"
  )
df_long$Intensity[is.na(df_long$Intensity)] <- 0
  df_long$Protein <- factor(df_long$Protein, levels = rownames(mat))
  df_long$y_index <- as.numeric(df_long$Protein)

  # ✅ FIXED factor levels
  df_long$Sample <- factor(df_long$Sample, levels = all_cols)

  # ============================
  # CLUSTER SIDEBAR
  # ============================

  cluster_levels <- sort(unique(clusters))
incProgress(
  0.70,
  detail = "Building cluster annotations"
)
  cluster_df <- data.frame(
    Protein = factor(rownames(mat), levels = rownames(mat)),
    Cluster = factor(clusters[rownames(mat)], levels = cluster_levels)
  )

  cluster_df$y_index <- as.numeric(cluster_df$Protein)

  # ✅ FIX: proper x column
  cluster_df$x <- factor("Cluster", levels = all_cols)


  dataset_positions <- col_info %>%
  dplyr::arrange(dataset, temperature, replicate) %>%
  dplyr::mutate(pos = seq_len(n())) %>%
  dplyr::group_by(dataset) %>%
  dplyr::summarise(
    midpoint = mean(pos) + 1, # +1 because Cluster bar occupies first column
    .groups = "drop"
  )

dataset_label_df <- data.frame(
  x = dataset_positions$midpoint,
  y = max(cluster_df$y_index) + 90,
  label = as.character(dataset_positions$dataset)
)

  # ============================
  # COLORS
  # ============================

  colors <- colorRampPalette(
    rev(RColorBrewer::brewer.pal(7, "RdYlBu"))
  )(100)

  # ============================
  # PLOT
  # ============================
incProgress(
  0.90,
  detail = "Rendering heatmap"
)
p <- ggplot() +

    geom_text(
      data = dataset_label_df,
      aes(
        x = x,
        y = y,
        label = label
      ),
      inherit.aes = FALSE,
      fontface = "bold",
      size = 5
    ) +

    geom_tile(
      data = df_long,
      aes(x = Sample, y = y_index, fill = Intensity)
    ) +

    scale_fill_gradientn(
      colours = colors,
      breaks = pretty(range(df_long$Intensity, na.rm = TRUE), n = 5),
      name = "Intensity",
      limits = range(df_long$Intensity, na.rm = TRUE)
    ) +

    ggnewscale::new_scale_fill() +

    # cluster bar
    geom_tile(
      data = cluster_df,
      aes(x = x, y = y_index, fill = Cluster)
    ) +

    scale_fill_brewer(palette = "Set2", name = "Cluster") +
    scale_x_discrete(
  labels = function(x) {
    ifelse(
      x == "Cluster",
      "",
      stringr::str_extract(
        x,
        "^\\d+(?:\\.\\d+)?_REP\\d+"
      )
    )
  }
) +

theme_minimal() +

scale_y_continuous(
  expand = expansion(mult = c(0, 0))
) +

coord_cartesian(
  clip = "off"
) +

theme(
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),

  axis.text.x = element_text(
    angle = 90,
    hjust = 1,
    vjust = 0
  ),

  axis.title.y = element_blank(),
  axis.title.x = element_blank(),

  panel.grid = element_blank(),

  plot.margin = margin(
    t = 20,
    r = 5,
    b = 5,
    l = 5
  )
)
incProgress(
  1,
  detail = "Done"
)
  return(p)
}


# =====================================================================
# Render the plot in the UI
# =====================================================================

output$reclustered_heatmap <- renderPlot({

  withProgress(
    message = "Generating reclustered heatmap",
    value = 0,
    {

      res <- reclustering_result()

      if (is.null(res)) return(NULL)

      print(
        make_reclustered_heatmap(res)
      )

    }
  )

})


# =====================================================================
# Download handler
# =====================================================================
output$download_reclustered_heatmap <- downloadHandler(
  filename = function() {
    paste0("reclustered_heatmap_", Sys.Date(), ".png")
  },
  content = function(file) {

    res <- isolate(reclustering_result())

    if (is.null(res)) {
      stop("Run reclustering first")
    }

    p <- isolate(make_reclustered_heatmap(res))

    ggplot2::ggsave(
      filename = file,
      plot = p,
      device = "png",
      width = 10,
      height = 8,
      dpi = 300
    )
  }
)



    # Reclustered data
    reclustered_data <- reactive({
      res <- reclustering_result()
      req(res)
      df <- res$df_subset
      clusters <- res$clusters
      df$Recluster <- factor(clusters)
      df
    })


    # Update recluster_num max when num_clusters changes
    observeEvent(input$num_clusters, {
      updateNumericInput(session, "recluster_num", max = input$num_clusters)
    })

observeEvent(input$run_reclustering, {

  # Trigger resize after plot is rendered
  session$sendCustomMessage("triggerResize", list())

})



    # Reclustered plot

# =====================================================================
# Helper function to build the reclustered trend plot
# =====================================================================
make_reclustered_plot <- function(res) {



  df <- reclustered_data()
incProgress(
  0.10,
  detail = "Loading reclustered proteins"
)
  # =====================================================
  # CASE 1 — INDIVIDUAL MODE
  # =====================================================
  if (res$mode == "individual") {
incProgress(
  0.25,
  detail = "Preparing cluster annotations"
)
    df$Recluster <- as.factor(df$Recluster)

    # cluster sizes
    counts <- df %>%
      dplyr::count(Recluster, name = "n")

    labeller_fn <- labeller(
      Recluster = function(x) {
        paste0(
          x, " (n=",
          counts$n[match(x, counts$Recluster)],
          ")"
        )
      }
    )

    # reshape
incProgress(
  0.45,
  detail = "Reshaping protein profiles"
)
    df_long <- df %>%
      pivot_longer(
        cols = -c(1, "Recluster"),
        names_to = "Replicate",
        values_to = "Intensity"
      ) %>%
      separate(Replicate, into = c("temperature", "Replicate"), sep = "_") %>%
      mutate(
        temperature = as.numeric(temperature),

        # NEW SEPARATE CONTROL
        Intensity = if (isTRUE(input$show_missing_recluster)) {
          Intensity
        } else {
          ifelse(Intensity == 0, NA, Intensity)
        }
      )

incProgress(
  0.65,
  detail = "Calculating protein trends"
)
    # protein-level curves
    df2 <- df_long %>%
      group_by(across(1), Recluster, temperature) %>%
      summarise(
        Mean_Intensity = mean(Intensity, na.rm = TRUE),
        .groups = "drop"
      )

    # average trends
    avg <- df_long %>%
      group_by(Recluster, temperature) %>%
      summarise(
        Mean_Intensity = mean(Intensity, na.rm = TRUE),
        .groups = "drop"
      )

      cluster_levels <- sort(unique(df$Recluster))

cluster_colors <- setNames(
  RColorBrewer::brewer.pal(
    max(3, length(cluster_levels)),
    "Set2"
  )[seq_along(cluster_levels)],
  as.character(cluster_levels)
)

facet_color_df <- data.frame(
  Recluster = factor(cluster_levels, levels = cluster_levels),
  x = min(df2$temperature, na.rm = TRUE),
  y = max(df2$Mean_Intensity, na.rm = TRUE),
  label = paste0(
    "Cluster ",
    cluster_levels,
    "\n(n=",
    counts$n[match(cluster_levels, counts$Recluster)],
    ")"
  )
)
incProgress(
  0.90,
  detail = "Rendering recluster plot"
)
p <- ggplot(df2, aes(
  x = temperature,
  y = Mean_Intensity,
  group = !!sym(names(df2)[1])
)) +

geom_line(
  color = "steelblue2",
  alpha = 0.1
) +

geom_line(
  data = avg,
  aes(group = 1),
  color = "steelblue4",
  linewidth = 1.5
) +

geom_point(
  data = facet_color_df,
  aes(
    x = x,
    y = y,
    colour = Recluster
  ),
  inherit.aes = FALSE,
  shape = 15,
  size = 5
) +

geom_text(
  data = facet_color_df,
  aes(
    x = x,
    y = y,
    label = label
  ),
  inherit.aes = FALSE,
  hjust = -0.2,
  vjust = 1,
  size = 4,
  fontface = "bold",
  colour = "black"
) +

scale_colour_manual(
  values = cluster_colors,
  guide = "none"
) +

facet_wrap(~ Recluster, nrow = 1) +

theme_minimal(base_size = 14) +

labs(
  title = "Reclustered trends (individual condition)",
  x = "Temperature (°C)",
  y = "Relative abundance"
)

incProgress(
  1,
  detail = "Done"
)
    return(p)
  }

  # =====================================================
  # CASE 2 — COMPARISON MODE
  # =====================================================
incProgress(
  0.35,
  detail = "Preparing comparison datasets"
)
df_long <- df %>%
  pivot_longer(
    cols = -c(1, "Recluster"),
    names_to   = c("temperature", "Replicate", "Dataset"),
    names_pattern = "^(\\d+(?:\\.\\d+)?)_(REP\\d+)(?:\\.(.+))?$",
    values_to  = "Intensity",
    values_drop_na = FALSE
  ) %>%
  mutate(
    temperature = as.numeric(temperature),
    Dataset = ifelse(is.na(Dataset), current_dataset_name, Dataset),
    Dataset = factor(Dataset, levels = unique(Dataset)),

    # NEW SEPARATE CONTROL (important)
    Intensity = if (isTRUE(input$show_missing_recluster)) {
      Intensity
    } else {
      ifelse(Intensity == 0, NA, Intensity)
    }
  )
  #  cluster sizes
  counts <- df %>%
    dplyr::count(Recluster, name = "n")

  labeller_fn <- labeller(
    Recluster = function(x) {
      paste0(
        x, " (n=",
        counts$n[match(x, counts$Recluster)],
        ")"
      )
    }
  )
incProgress(
  0.60,
  detail = "Calculating protein trends"
)

  # protein-level curves
  df_protein <- df_long %>%
    group_by(across(1), Recluster, Dataset, temperature) %>%
    summarise(
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      .groups = "drop"
    )

  # average trends
  avg <- df_long %>%
    group_by(Recluster, Dataset, temperature) %>%
    summarise(
      Mean_Intensity = mean(Intensity, na.rm = TRUE),
      .groups = "drop"
    )

    cluster_levels <- sort(unique(df$Recluster))

cluster_colors <- setNames(
  RColorBrewer::brewer.pal(
    max(3, length(cluster_levels)),
    "Set2"
  )[seq_along(cluster_levels)],
  as.character(cluster_levels)
)

facet_color_df <- data.frame(
  Recluster = factor(cluster_levels, levels = cluster_levels),
  x = min(df_protein$temperature, na.rm = TRUE),
  y = max(df_protein$Mean_Intensity, na.rm = TRUE)
)

  # colors
  datasets_unique <- unique(df_long$Dataset)

  curve_colors <- setNames(
    scales::hue_pal()(length(datasets_unique)),
    datasets_unique
  )

  trend_colors <- setNames(
    scales::hue_pal(l = 40)(length(datasets_unique)),
    datasets_unique
  )
incProgress(
  0.90,
  detail = "Rendering comparison plot"
)
  p <- ggplot(
    df_protein,
    aes(
      x = temperature,
      y = Mean_Intensity,
      group = interaction(Dataset, !!sym(names(df)[1])),
      color = Dataset
    )
  ) +

    geom_line(alpha = 0.08, linewidth = 0.4) +

    geom_line(
      data = avg,
      aes(
        x = temperature,
        y = Mean_Intensity,
        group = Dataset,
        color = Dataset
      ),
      linewidth = 2.2
    ) +

    scale_color_manual(values = curve_colors, guide = "none") +
    ggnewscale::new_scale_color() +

    geom_line(
      data = avg,
      aes(
        x = temperature,
        y = Mean_Intensity,
        group = Dataset,
        color = Dataset
      ),
      linewidth = 1
    ) +

scale_color_manual(
  values = trend_colors,
  name = "Dataset"
) +

ggnewscale::new_scale_colour() +

geom_point(
  data = facet_color_df,
  aes(
    x = x,
    y = y,
    colour = Recluster
  ),
  inherit.aes = FALSE,
  shape = 15,
  size = 5
) +


scale_colour_manual(
  values = cluster_colors,
  guide = "none"
) +

facet_wrap(
  ~ Recluster,
  nrow = 1,
  labeller = labeller_fn
) +

    labs(
      title = "Reclustered comparison across datasets",
      x = "Temperature (°C)",
      y = "Relative abundance"
    ) +

    theme_light(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.border = element_rect(
        color = "grey80",
        fill = NA,
        linewidth = 0.5
      )
    )
incProgress(
  1,
  detail = "Done"
)

  return(p)
}


# =====================================================================
# Render plot in UI
# =====================================================================


output$reclustered_plot <- renderPlot({

  withProgress(
    message = "Generating reclustered trend plot",
    value = 0,
    {

      res <- reclustering_result()

      if (is.null(res)) return(NULL)

      p <- make_reclustered_plot(res)

      p

    }
  )

})



# this is so collapsed things still work when the tab is not visible
outputOptions(output, "reclustered_plot", suspendWhenHidden = FALSE)
outputOptions(output, "reclustered_heatmap", suspendWhenHidden = FALSE)

outputOptions(output, "recluster_recommendation", suspendWhenHidden = FALSE)
outputOptions(output, "recluster_recommendation_plot", suspendWhenHidden = FALSE)

# =====================================================================
# Download handler (PNG)
# =====================================================================
output$download_reclustered_plot <- downloadHandler(
  filename = function() {
    paste0("reclustered_plot_", Sys.Date(), ".png")
  },
  content = function(file) {

    res <- isolate(reclustering_result())

    if (is.null(res)) {
      stop("Run reclustering first")
    }

    p <- isolate(make_reclustered_plot(res))

    ggplot2::ggsave(
      filename = file,
      plot = p,
      device = "png",
      width = 10,
      height = 5,
      dpi = 300
    )
  }
)


outputOptions(output, "download_reclustered_plot", suspendWhenHidden = FALSE)
outputOptions(output, "download_reclustered_heatmap", suspendWhenHidden = FALSE)

    # Return clustered data if needed

return(list(

  # =========================
  # CLUSTERING
  # =========================
  clustering = list(
    data = dataclustered,
    result = clustering_result,
    mode = reactive({
      req(clustering_result())
      clustering_result()$mode
    })
  ),

  # =========================
  # RECLUSTERING
  # =========================
  reclustering = list(
    data = reclustered_data,
    result = reclustering_result,
    mode = reactive({
      req(reclustering_result())
      reclustering_result()$mode
    })
  ),

  # =========================
  # PARAMETERS
  # =========================
  num_clusters = reactive(input$num_clusters)

))


  })
}


modelPanelModuleUI <- function(id, label) {
  ns <- NS(id)

tagList(

  tabsetPanel(

        id = ns("model_fitting"),


        # =========================
        # LIMMA TAB
        # =========================
        tabPanel(
          "Limma",
          fluidRow(
              box(
              title = tagList(
              "Test parameters",
              tags$span(
                `data-toggle` = "tooltip",
                `data-placement` = "right",
                `data-html` = "true",
                title = HTML(
                  "<b>Limma analysis of TPP data</b><br><br>

                  Limma is applied to protein thermal profiles to identify proteins with unusual melting behaviour.<br><br>

                  <b>Slope different from 0</b><br>
                  Identifies proteins whose abundance changes significantly across the temperature gradient.<br><br>

                  <b>Slope different from average trend</b><br>
                  Identifies proteins whose melting behaviour deviates significantly from the overall proteome trend.<br><br>

                  Proteins with large positive or negative logFC values may represent proteins with altered thermal stability."
                ),
                tags$i(
                  class = "fa fa-question-circle",
                  style = "margin-left:8px; color:#337ab7; cursor:pointer;"
                )
              )
            ),
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              radioButtons(
                inputId = ns("Null_hypothesis_options"),
                label = "Select the Null hypothesis to test:",
                choices = c(
                  "Is the slope of the proteins significantly different from 0?" = 1,
                  "Is the slope of the proteins significantly different from the average trend?" = 2
                ),
                selected = 2,
                inline = FALSE
              )
            ),

          box(
            title = "Limma results",
            width = 6,
            status = "primary",
            solidHeader = TRUE,

            withSpinner(
              DTOutput(ns("limma_results_table"))
            )
          ),



            box(
              title = "Plot options",
              width = 6,
              status = "primary",
              solidHeader = TRUE,
              uiOutput(ns("cluster_options_ui")),
              numericInput(
                inputId = ns("quantile_cutoff"),
                label = "Enter quantile cutoff percentage (e.g., 5 for 5%)",
                value = 5,
                min = 0,
                max = 50,
                step = 0.1
              ),
              numericInput(
                inputId = ns("label_top_n"),
                label = "Number of top proteins to label (by logFC)",
                value = 10,
                min = 1,
                step = 1
              )
            ),

            box(
              title = "Volcano plot",
              width = 12,
              status = "primary",
              solidHeader = TRUE,
              plotOutput(ns("limma_volcano_plot"), height = 450) %>% withSpinner(
color = "#0EA5A5"
),
              downloadButton(ns("download_limma_volcano_plot"), "Download volcano plot")
            )
          )
        ),

        # =========================
        # MELTING CURVE TAB
        # =========================
        tabPanel(
          "Melting curve",

          fluidRow(
            box(
              title = tagList(
                "Melting curve fitting results",
                tags$span(
                  `data-toggle` = "tooltip",
                  `data-placement` = "right",
                  `data-html` = "true",
                  title = HTML(
                    "<b>Melting curve fitting</b><br><br>

                    A sigmoidal melting model is fitted independently to every protein thermal profile.<br><br>

                    The fitted parameters are used to estimate protein melting temperatures (Tm), curve quality metrics (R²), slope and plateau values.<br><br>

                    Only proteins passing the selected quality filters are displayed."
                  ),
                  tags$i(
                    class = "fa fa-question-circle",
                    style = "margin-left:8px; color:#337ab7; cursor:pointer;"
                  )
                )
              ),
              width = 8,
              status = "primary",
              solidHeader = TRUE,
              div(
                style = "overflow-x: auto;",
                DTOutput(ns("melting_curve_results_table")) %>% withSpinner(
color = "#0EA5A5"
)
              )
            ),

            box(
              title = tagList(
                "Control parameters",

                tags$span(
                  `data-toggle` = "tooltip",
                  `data-placement` = "right",
                  `data-html` = "true",
                  title = HTML(
                    "<b>Curve fitting parameters</b><br><br>

                    Melting curves are fitted using a sigmoidal model.<br><br>

                    <b>Program estimation</b><br>
                    Starting values are estimated automatically:

                    <br>• Tm = temperature closest to the median abundance

                    <br>• k = 1 + SD(abundance)

                    <br>• Plateau = minimum abundance value

                    <br><br>

                    <b>Manual values</b><br>

                    User-defined values are used only as starting guesses for the nonlinear fit.<br><br>

                    R², slope and plateau filters determine which fitted proteins are retained for downstream analyses."
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
                inputId = ns("starting_values"),
                label = "Select how to determine initial values for nonlinear fitting:",
                choices = c("Program estimation", "Manually specified"),
                selected = "Program estimation"
              ),

              conditionalPanel(
                condition = sprintf(
                  "input['%s'] == 'Manually specified'",
                  ns("starting_values")
                ),
                box(
                  title = "Manual parameters",
                  width = 12,
                  status = "warning",
                  solidHeader = TRUE,
                  numericInput(ns("param_tm"), "Tm:", value = 85),
                  numericInput(ns("param_k"), "k:", value = 2),
                  numericInput(ns("param_p"), "p (plateau):", value = 0.05)
                )
              ),




              numericInput(
                inputId = ns("R2_threshold"),
                label = "Filter: R² ≥",
                value = 0.7,
                min = 0,
                max = 1,
                step = 0.05
              ),

              numericInput(
                inputId = ns("slope_threshold"),
                label = "Filter: slope ≤",
                value = -0.06,
                step = 0.005
              ),

              numericInput(
                ns("plateau_threshold"),
                "Filter p ≤",
                value = 0.2,
                min = 0,
                max = 1
              ),
              numericInput(
              ns("tm_min"),
              "Filter: Tm ≥",
              value = 30,
              step = 1
            ),

            numericInput(
              ns("tm_max"),
              "Filter: Tm ≤",
              value = 95,
              step = 1
            ),
            actionButton(
              ns("run_curve_fitting"),
              "Run melting curve fitting",
              icon = icon("play"),
              class = "btn-success"
            ),


            )
          ),

          fluidRow(
            box(
              title = "Melting curve plot",
              width = 12,
              status = "info",
              solidHeader = TRUE,
              plotOutput(ns("melting_curve_example_plot"), height = 400) %>% withSpinner(
color = "#0EA5A5"
),
              downloadButton(ns("download_melting_curve_example_plot"), "Download melting curve plot")
            )
          ),
          fluidRow(
              box(
                title = tagList(
                "Melting temperature distribution",
                tags$span(
                  `data-toggle` = "tooltip",
                  `data-placement` = "right",
                  `data-html` = "true",
                  title = HTML(
                    "<b>Melting temperature (Tm) distribution</b><br><br>

                    This plot summarizes the melting temperature estimates obtained from all proteins whose curve fits pass the currently selected quality filters.<br><br>

                    Each condition is displayed using their corresponding set of accepted protein fits.<br><br>

                    <b>Show only shared proteins</b><br>
                    Restricts the analysis to proteins with valid curve fits shared between conditions.<br><br>"
                  ),
                  tags$i(
                    class = "fa fa-question-circle",
                    style = "margin-left:8px; color:#337ab7; cursor:pointer;"
                  )
                )
              ),
                width = 12,
                status = "info",
                solidHeader = TRUE,

                checkboxInput(
                  ns("shared_proteins_only"),
                  "Show only shared proteins",
                  value = FALSE
                ),

                br(),

                plotOutput(
                  ns("tm_distribution_plot"),
                  height = 300
                ) %>% withSpinner(
                color = "#0EA5A5"
                ),
                br(),
                downloadButton(
                  ns("download_tm_distribution"),
                  "Download Tm distribution plot"
                ),
              )
          )
        )
      )
    )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

modelPanelModuleServer <- function(id, normalized_data, dataclustered = reactive(NULL), reclustered_data = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # UI for cluster highlighting options


output$cluster_options_ui <- renderUI({

  # Require normalized data only
  req(normalized_data())
  # ---- BUILD CHOICES ----
  choices <- c("None")

  # dataclustered
  dc <- tryCatch(dataclustered(), error = function(e) NULL)
  if (!is.null(dc) && (is.data.frame(dc) || is.matrix(dc)) && NROW(dc) > 0) {
    choices <- c(choices, "Clusters")
  }

  # reclustered_data
  rdc <- tryCatch(reclustered_data(), error = function(e) NULL)
  if (!is.null(rdc) && (is.data.frame(rdc) || is.matrix(rdc)) && NROW(rdc) > 0) {
    choices <- c(choices, "Sub-clusters")
  }

  # ---- UI ----
  radioButtons(
    inputId = ns("highlight_clusters"),
    label   = "Highlight clusters in volcano plot:",
    choices = choices,
    selected = "None",
    inline   = TRUE
  )
})



    # Limma results as a reactive for reuse
limma_results <- reactive({

  withProgress(
    message = "Running Limma analysis",
    value = 0,
    {

      req(normalized_data())

      incProgress(
        0.1,
        detail = "Loading normalized data"
      )

      df <- normalized_data()
      mat <- as.matrix(df[,-1])
      incProgress(
        0.25,
        detail = "Preparing expression matrix"
      )
      rownames(mat) <- df[[1]]

      temperature <- as.numeric(sub("_.*", "", colnames(mat)))
            incProgress(
        0.45,
        detail = "Building design matrix"
      )
      if (any(is.na(temperature))) stop("Temperature extraction failed. Check column name format.")

      design <- model.matrix(~ temperature)
            incProgress(
        0.60,
        detail = "Fitting linear models"
      )
      fit <- lmFit(mat, design)
            incProgress(
        0.80,
        detail = "Applying empirical Bayes moderation"
      )
      fit <- eBayes(fit, robust = TRUE)

      results <- topTable(fit, coef = "temperature", number = Inf, adjust.method = "BH", sort.by = "P")
            incProgress(
        1,
        detail = "Limma analysis complete"
      )

results$Accession <- rownames(results)

results

    }
  )

})

    # Limma results for centered data
    centered_limma_results <- reactive({
      req(normalized_data())
      df <- normalized_data()
      mat <- as.matrix(df[,-1])
      rownames(mat) <- df[[1]]

      temperature <- as.numeric(sub("_.*", "", colnames(mat)))
      if (any(is.na(temperature))) stop("Temperature extraction failed. Check column name format.")

      average_profile <- colMeans(mat, na.rm = TRUE)
      centered_matrix <- sweep(mat, 2, average_profile)

      design <- model.matrix(~ temperature)
      fit <- lmFit(centered_matrix, design)
      fit <- eBayes(fit)

      results <- topTable(fit, coef = "temperature", number = Inf, adjust.method = "BH", sort.by = "P")
      results$Accession <- rownames(results)
      results
    })

    # Reactive expression to compute results based on selected null hypothesis
    hypothesis_result <- reactive({
      req(input$Null_hypothesis_options)
      if (input$Null_hypothesis_options == 1) {
        limma_results()
      } else if (input$Null_hypothesis_options == 2) {
        centered_limma_results()
      } else {
        NULL
      }
    })


# =========================
# Helpers / reactives
# =========================

# Base results as a reactive (keeps row order)
base_results <- reactive({
  req(hypothesis_result())
  hypothesis_result()
})

# Row selection from the DT table -> protein IDs (rownames)
selected_proteins <- reactive({
  idx <- input$limma_results_table_rows_selected
  if (is.null(idx) || length(idx) == 0) return(character(0))

  # The table is built from hypothesis_result() with the same row order,
  # so we can map the selected row indices to rownames (protein IDs)
  rownames(base_results())[idx]
})

# =========================
# Data table (unchanged behavior)
# =========================

output$limma_results_table <- renderDT({
  df <- base_results() %>% dplyr::select(-Accession)

  numeric_cols <- names(df)[sapply(df, is.numeric)]
  cols_to_exclude <- names(df)[c(4, 5)]                 # keep your original rule
  cols_to_round <- setdiff(numeric_cols, cols_to_exclude)

  DT::datatable(
    df,
    options = list(pageLength = 5, scrollX = TRUE),
    selection = "multiple"
  ) %>%
    DT::formatRound(columns = cols_to_round, digits = 2) %>%
    DT::formatSignif(columns = cols_to_exclude, digits = 3)
})


    # Reactive quantile cutoff
    cutoff_percent <- reactive({
      req(input$quantile_cutoff)
      input$quantile_cutoff / 100
    })

    # Limma volcano plot
    
# ===== Function that builds the volcano plot (reused for both display & download) =====

make_volcano_plot <- function() {
withProgress(
  message = "Building volcano plot",
  value = 0,
  {
  # Required inputs
  req(base_results(), cutoff_percent(), input$label_top_n)

  # Optional sources
  clustered   <- tryCatch(dataclustered(),     error = function(e) NULL)
  reclustered <- tryCatch(reclustered_data(),  error = function(e) NULL)

  # Main results with protein IDs
  results <- base_results()
  results <- tibble::rownames_to_column(results, var = "proteins")
  results$negLog10AdjP <- -log10(results$adj.P.Val)
incProgress(
  0.15,
  detail = "Loading Limma results"
)

  top_n <- max(1, round(input$label_top_n))

  # Cutoffs
  lower_cutoff <- quantile(results$logFC, cutoff_percent(), na.rm = TRUE)
  upper_cutoff <- quantile(results$logFC, 1 - cutoff_percent(), na.rm = TRUE)
incProgress(
  0.30,
  detail = "Calculating significance cutoffs"
)
  # Expression classification
  results <- results %>%
    dplyr::mutate(
      Expression = dplyr::case_when(
        logFC < lower_cutoff & adj.P.Val < 0.05 ~ "Downregulated",
        logFC > upper_cutoff & adj.P.Val < 0.05 ~ "Upregulated",
        TRUE ~ "No significant"
      )
    )
incProgress(
  0.45,
  detail = "Classifying proteins"
)
  # Default color scale
  default_colors <- scale_color_manual(values = c(
    "Downregulated" = "#00AFBB",
    "No significant" = "grey",
    "Upregulated"   = "#bb0c00"
  ))

  boxinput <- input$highlight_clusters
if (is.null(boxinput) ||
    length(boxinput) == 0 ||
    is.na(boxinput)) {

  boxinput <- "None"
}

  color_scale <- default_colors  # fallback
incProgress(
  0.60,
  detail = "Applying cluster annotations"
)
  # ---- CLUSTER MODE ----
  if (!is.null(boxinput) &&
    !is.na(boxinput) &&
    boxinput == "Clusters" &&
      !is.null(clustered) &&
      (is.data.frame(clustered) || is.matrix(clustered)) &&
      NROW(clustered) > 0) {

    colnames(clustered)[1] <- "proteins"
    clusters <- clustered %>% dplyr::select(proteins, Cluster)
    results <- dplyr::left_join(results, clusters, by = "proteins")

    color_scale <- NULL  # let ggplot auto-assign colors for clusters

  # ---- SUBCLUSTER MODE ----
  } else if (!is.null(boxinput) &&
           !is.na(boxinput) &&
           boxinput == "Sub-clusters" &&
             !is.null(reclustered) &&
             (is.data.frame(reclustered) || is.matrix(reclustered)) &&
             NROW(reclustered) > 0) {

    colnames(reclustered)[1] <- "proteins"
    subcl <- reclustered %>% dplyr::select(proteins, Recluster)
    results <- dplyr::left_join(results, subcl, by = "proteins")

    results$Cluster <- ifelse(is.na(results$Recluster), "all", results$Recluster)

    # Subcluster colors
    cluster_levels <- unique(results$Cluster)
    subcluster_colors <- setNames(
      scales::hue_pal()(length(cluster_levels)),
      cluster_levels
    )
    subcluster_colors["all"] <- "lightgray"

    color_scale <- scale_color_manual(values = subcluster_colors)

  # ---- NONE MODE (default) ----
  } else {
    results$Cluster <- results$Expression
    color_scale <- default_colors
  }

  # Labels (top up/down)
  label_data <- results %>%
    dplyr::filter(!is.na(logFC), adj.P.Val < 0.05) %>%
    dplyr::arrange(dplyr::desc(logFC)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::bind_rows(
      results %>%
        dplyr::filter(!is.na(logFC), adj.P.Val < 0.05) %>%
        dplyr::arrange(logFC) %>%
        dplyr::slice_head(n = top_n)
    )
incProgress(
  0.75,
  detail = "Preparing labels"
)
  # Base plot
  p <- ggplot2::ggplot(results, ggplot2::aes(x = logFC, y = negLog10AdjP, color = Cluster)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    ggplot2::geom_vline(
  xintercept = c(lower_cutoff, upper_cutoff),
  color = "black",
  linetype = "dashed"
) +
    ggplot2::geom_text(data = label_data, ggplot2::aes(label = proteins),
                       vjust = -1, size = 3, show.legend = FALSE) +
    ggplot2::labs(x = "logFC", y = "-log10(adj.P.Val)") +
    ggplot2::theme_minimal()
incProgress(
  0.90,
  detail = "Generating volcano plot"
)
  if (!is.null(color_scale)) {
    p <- p + color_scale
  }

  # ---- Highlight selection from table ----

# ---- Highlight selection from table ----
sel <- selected_proteins()

if (length(sel) > 0) {

  highlight_df <- results %>% dplyr::filter(proteins %in% sel)

  # Orange points
  p <- p +
    ggplot2::geom_point(
      data = highlight_df,
      color = "#FF8C00", fill = "#FF8C00",
      size = 4, alpha = 1, shape = 21, stroke = 0
    ) +
    # Labels using FIRST COLUMN ("proteins")
    ggplot2::geom_text(
      data = highlight_df,
      ggplot2::aes(label = proteins),
      vjust = -1,
      color = "#FF8C00",
      fontface = "bold",
      size = 4,
      show.legend = FALSE
    )
}

incProgress(
  1,
  detail = "Done"
)
  return(p)

  }
)

}



# ===== Render the plot in the UI =====
output$limma_volcano_plot <- renderPlot({
  make_volcano_plot()
})


# ===== Download handler =====
output$download_limma_volcano_plot <- downloadHandler(
  filename = function() {
    paste0("volcano_plot_", Sys.Date(), ".png")
  },
  content = function(file) {
    ggsave(
      filename = file,
      plot = make_volcano_plot(),
      width = 8,
      height = 6,
      dpi = 300
    )
  }
)

    return(list(
        hypothesis_result = hypothesis_result,
        limma_results = limma_results,
        centered_limma_results = centered_limma_results
    ))
  })
}



meltingCurveModuleServer <- function(
  id,
  normalized_data,
  all_normalized_data = NULL,
  normalization_results,
  already_lowestnorm
){
  moduleServer(id, function(input, output, session) {
    ns <- session$ns


#check all conditions
all_datasets_lowest_temp <- reactive({

  if (isTRUE(already_lowestnorm())) {
    return(TRUE)
  }

  status <- sapply(
    names(reactiveValuesToList(normalization_results)),
    function(nm){

      tryCatch(
        isTRUE(
          normalization_results[[nm]]$lowest_temp_normalized()
        ),
        error = function(e) FALSE
      )

    }
  )

  print(status)

  all(status)

})
    # ------------------------------------------------------------
    # Models
    # ------------------------------------------------------------

    # Standard TPP logistic
    logistic_model <- function(T, Tm, k, p) {
      (1 - p) / (1 + exp((T - Tm) / k)) + p
    }


    compute_rmse <- function(obs, pred) sqrt(mean((obs - pred)^2))

    # ------------------------------------------------------------
    # Extract temperatures
    # ------------------------------------------------------------
    temperature_vector <- reactive({
      col_names <- colnames(normalized_data())[-1]
      temperature_values <- as.numeric(sub("^([^_]+)_.*$", "\\1", col_names))
      if (any(is.na(temperature_values))) stop("Temperature extraction failed")
      temperature_values
    })

    # ------------------------------------------------------------
    # protein_df
    # ------------------------------------------------------------
protein_df <- reactive({

df <- normalized_data()

df <- df %>%
dplyr::filter(complete.cases(df[, -1]))

      df <- as.data.frame(df)
      rownames(df) <- df[, 1]
      df <- df %>% dplyr::select(-1)
      df[is.nan(as.matrix(df))] <- 0
      df
    })

    # ------------------------------------------------------------
    # STANDARD LOGISTIC FIT
    # ------------------------------------------------------------
    run_standard_logistic_fit <- function() {

      df <- protein_df()

      result <- data.frame(
        Protein = rownames(df),
        Tm = NA, k = NA, p = NA,
        R2 = NA, AUC = NA,
        slope = NA
      )

      numeric_cols <- sapply(df, is.numeric)
      col_names <- names(df)[numeric_cols]
      temp_raw_all <- as.numeric(sub("_.*", "", col_names))

      withProgress(

  message = "Fitting melting curves",
  detail = "Starting...",
  value = 0,

  {

    n_proteins <- nrow(df)

    for (i in seq_len(n_proteins)) {

      incProgress(
        1 / n_proteins,
        detail = paste(
          "Protein",
          i,
          "of",
          n_proteins
        )
      )


        raw_y <- as.numeric(df[i, numeric_cols])
        df_temp <- data.frame(temp = temp_raw_all, y = raw_y)
        df_temp <- df_temp[order(df_temp$temp), ]

        temp <- df_temp$temp
        y <- df_temp$y

        if (input$starting_values == "Manually specified") {
          start_Tm <- input$param_tm
          start_k  <- input$param_k
          start_p  <- input$param_p
        } else {
          start_Tm <- temp[which.min(abs(y - median(y)))]
          start_k  <- 1 + sd(y)
          start_p  <- min(y)
        }

        fit <- tryCatch({
          nlsLM(
            y ~ logistic_model(temp, Tm, k, p),
            start = list(Tm = start_Tm, k = start_k, p = start_p),
            control = nls.lm.control(maxiter = 500)
          )
        }, error = function(e) NULL)

        if (is.null(fit)) next

        coef_vals <- coef(fit)
        y_pred <- predict(fit)

        R2 <- 1 - sum((y - y_pred)^2) / sum((y - mean(y))^2)
        slope <- -(1 - coef_vals["p"]) / (4 * coef_vals["k"])
        auc_val <- trapz(temp, y)

      result[i, c("Tm", "k", "p", "R2", "AUC", "slope")] <-
        c(
          coef_vals["Tm"],
          coef_vals["k"],
          coef_vals["p"],
          R2,
          auc_val,
          slope
        )
    }

    result

  }
)
    }

    # ------------------------------------------------------------
    # Main reactive: run selected model
    # ------------------------------------------------------------
# ------------------------------------------------------------
# Main reactive: run selected model
# ------------------------------------------------------------

raw_fit_results <- eventReactive(
  input$run_curve_fitting,
  {

    if (!all_datasets_lowest_temp()) {

      showNotification(
        paste(
          "Melting Curve analysis requires",
          "'Normalize to lowest temperature'",
          "to be selected in all datasets."
        ),
        type = "error",
        duration = 8
      )

      return(NULL)
    }

    run_standard_logistic_fit()

  },
  ignoreInit = TRUE
)


    # ------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------
delayed_filters <- reactive({

  list(
    R2 = input$R2_threshold,
    slope = input$slope_threshold,
    p = input$plateau_threshold,
    tm_min = input$tm_min,
    tm_max = input$tm_max
  )

}) %>%
  debounce(2000)

fit_results <- reactive({
req(raw_fit_results())
  filt <- delayed_filters()

  raw_fit_results() %>%
    filter(
      !is.na(R2),
      R2 >= filt$R2,
      slope <= filt$slope,
      p <= filt$p,
      Tm >= filt$tm_min,
      Tm <= filt$tm_max
    )
})



tm_boxplot_data <- eventReactive(
  input$run_curve_fitting,
  {

  req(all_normalized_data)
if (!all_datasets_lowest_temp()) {

showNotification(

paste(

"Tm distribution analysis requires",
"'Normalize to lowest temperature'",
"to be selected in all datasets."

),

type = "error",

duration = 8

)
return(NULL)
}
  datasets <- all_normalized_data()

  withProgress(

    message = "Calculating Tm distributions",
    detail = "Preparing datasets...",
    value = 0,

    {

      out_list <- list()
total_proteins <- sum(
  sapply(
    datasets,
    function(x) if (is.null(x)) 0 else nrow(x)
  )
)

processed_proteins <- 0

  for (nm in names(datasets)) {

    df <- datasets[[nm]]
    
    if (is.null(df))
    next
    
    df <- df %>%
    dplyr::filter(
    complete.cases(df[, -1])
    )

    df <- as.data.frame(df)

    rownames(df) <- df[,1]

    df <- df %>%
      dplyr::select(-1)

    df[is.nan(as.matrix(df))] <- 0

    numeric_cols <- sapply(df, is.numeric)

    result <- data.frame(
      Protein = rownames(df),
      Tm = NA,
      k = NA,
      p = NA,
      R2 = NA,
      AUC = NA,
      slope = NA
    )

    col_names <- names(df)[numeric_cols]
    temp_raw_all <- as.numeric(sub("_.*", "", col_names))

    for (i in seq_len(nrow(df))) {
processed_proteins <- processed_proteins + 1

incProgress(
  1 / total_proteins,
  detail = paste(
    nm,
    "- protein",
    processed_proteins,
    "of",
    total_proteins
  )
)
      raw_y <- as.numeric(df[i, numeric_cols])

      df_temp <- data.frame(
        temp = temp_raw_all,
        y = raw_y
      )

      df_temp <- df_temp[order(df_temp$temp), ]

      temp <- df_temp$temp
      y <- df_temp$y

      start_Tm <- temp[which.min(abs(y - median(y)))]
      start_k  <- 1 + sd(y)
      start_p  <- min(y)

      fit <- tryCatch({

        nlsLM(
          y ~ logistic_model(temp, Tm, k, p),
          start = list(
            Tm = start_Tm,
            k  = start_k,
            p  = start_p
          ),
          control = nls.lm.control(maxiter = 500)
        )

      }, error = function(e) NULL)

      if (is.null(fit))
        next

      coef_vals <- coef(fit)

      y_pred <- predict(fit)

      R2 <- 1 -
        sum((y - y_pred)^2) /
        sum((y - mean(y))^2)

      slope <- -(1 - coef_vals["p"]) /
        (4 * coef_vals["k"])

      auc_val <- trapz(temp, y)

      result[i, c(
        "Tm",
        "k",
        "p",
        "R2",
        "AUC",
        "slope"
      )] <- c(
        coef_vals["Tm"],
        coef_vals["k"],
        coef_vals["p"],
        R2,
        auc_val,
        slope
      )
    }

filt <- delayed_filters()

result <- result %>%
  dplyr::filter(
    !is.na(R2),
    R2 >= filt$R2,
    slope <= filt$slope,
    p <= filt$p,
    Tm >= filt$tm_min,
    Tm <= filt$tm_max
  )


    result$Dataset <- nm

    out_list[[nm]] <- result
  }
combined <- dplyr::bind_rows(out_list)

if (isTRUE(input$shared_proteins_only)) {

  # ------------------------------------------------
  # Infer organism from dataset name
  # ------------------------------------------------

  combined$OrganismGroup <- dplyr::case_when(

    grepl("H.sapiens|Homo sapiens|human",
          combined$Dataset,
          ignore.case = TRUE) ~ "Human",

    TRUE ~ "Sulfolobus"
  )

  # ------------------------------------------------
  # Shared proteins INSIDE EACH ORGANISM
  # ------------------------------------------------

  shared_by_organism <- combined %>%

    dplyr::group_by(OrganismGroup) %>%

    dplyr::group_modify(~{

      required_sets <- dplyr::n_distinct(.x$Dataset)

      shared_ids <- .x %>%

        dplyr::group_by(Protein) %>%

        dplyr::summarise(
          n_sets = dplyr::n_distinct(Dataset),
          .groups = "drop"
        ) %>%

        dplyr::filter(
          n_sets == required_sets
        ) %>%

        dplyr::pull(Protein)

      .x %>%
        dplyr::filter(
          Protein %in% shared_ids
        )
    }) %>%

    dplyr::ungroup()

  combined <- shared_by_organism
}

combined

    }
  )
},
ignoreInit = TRUE
)

    # ------------------------------------------------------------
    # Table
    # ------------------------------------------------------------
    output$melting_curve_results_table <- DT::renderDataTable({
      df <- fit_results()
      DT::datatable(df,
        options = list(pageLength = 5, scrollX = TRUE),
        selection = "single"
      ) %>%
        DT::formatRound(columns = names(df)[sapply(df, is.numeric)], digits = 2)
    })

    selected_protein <- reactive({
      req(input$melting_curve_results_table_rows_selected)
      fit_results()[input$melting_curve_results_table_rows_selected, "Protein"]
    })

    # ------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------

make_melting_curve_example_plot <- function() {

  req(selected_protein(), protein_df())

  df <- protein_df()
  protein_name <- selected_protein()
  protein_index <- which(rownames(df) == protein_name)

  numeric_cols <- sapply(df, is.numeric)
  raw_y <- as.numeric(df[protein_index, numeric_cols])
  col_names <- names(df)[numeric_cols]

  temp_raw <- as.numeric(sub("_.*$", "", col_names))
  rep_raw  <- sub("^.*_(.*)$", "\\1", col_names)

  df_temp <- data.frame(
    temp = temp_raw,
    y = raw_y,
    rep = factor(rep_raw)
  )

  df_temp <- df_temp[order(df_temp$temp), ]

  temp <- df_temp$temp
  y <- df_temp$y

  # --------------------------------------------------------
  # Fit curve based on selected model
  # --------------------------------------------------------

    # Starting values
    if (input$starting_values == "Manually specified") {
      start_Tm <- input$param_tm
      start_k  <- input$param_k
      start_p  <- input$param_p
    } else {
      start_Tm <- temp[which.min(abs(y - median(y)))]
      start_k  <- 1 + sd(y)
      start_p  <- min(y)
    }

    fit <- tryCatch({
      nlsLM(
        y ~ logistic_model(temp, Tm, k, p),
        start = list(Tm = start_Tm, k = start_k, p = start_p),
        control = nls.lm.control(maxiter = 500)
      )
    }, error = function(e) NULL)

    if (is.null(fit)) 
      return(ggplot() + ggtitle("Fit failed"))

    y_pred <- predict(fit)


  fit_df <- data.frame(temp = temp, y_pred = y_pred)

  # --------------------------------------------------------
  # Final Plot
  # --------------------------------------------------------
  p <- ggplot(df_temp, aes(x = temp, y = y)) +
    geom_point(aes(shape = rep, color = rep), size = 3, alpha = 0.8) +
    geom_line(data = fit_df, aes(x = temp, y = y_pred), size = 1.3) +
    labs(
      title = protein_name,
      x = "Temperature (°C)",
      y = "Intensity",
      shape = "Replicate"
    ) +
    theme_minimal(base_size = 14)

  return(p)
}


# =====================================================================
# Render the plot
# =====================================================================
output$melting_curve_example_plot <- renderPlot({
  make_melting_curve_example_plot()
})
# =====================================================================
# Tm distribution plot
# =====================================================================

output$tm_distribution_plot <- renderPlot({

  df <- tm_boxplot_data()

  validate(
    need(
      nrow(df) > 0,
      "No proteins pass current filters."
    )
  )

  ggplot(
    df,
    aes(
      x = Dataset,
      y = Tm,
      fill = Dataset
    )
  ) +

    geom_jitter(
      width = 0.15,
      alpha = 0.15,
      color = "grey40"
    ) +

    geom_boxplot(
      width = 0.5,
      alpha = 0.8,
      outlier.shape = NA
    ) +

    theme_minimal(base_size = 14) +

    labs(
      title = "Tm Distribution Across Datasets",
      x = NULL,
      y = "Melting Temperature (Tm)"
    )
})

# =====================================================================
# Download handler
# =====================================================================
output$download_melting_curve_example_plot <- downloadHandler(
  filename = function() {
    paste0("melting_curve_", selected_protein(), "_", Sys.Date(), ".png")
  },
  content = function(file) {
    p <- make_melting_curve_example_plot()
    ggsave(
      filename = file,
      plot = p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
)
output$download_tm_distribution <- downloadHandler(

  filename = function() {
    paste0(
      "Tm_distribution_",
      Sys.Date(),
      ".png"
    )
  },

  content = function(file) {

    df <- tm_boxplot_data()

    p <- ggplot(
      df,
      aes(
        x = Dataset,
        y = Tm,
        fill = Dataset
      )
    ) +

      geom_jitter(
        width = 0.15,
        alpha = 0.15,
        color = "grey40"
      ) +

      geom_boxplot(
        width = 0.5,
        alpha = 0.8,
        outlier.shape = NA
      ) +

      theme_minimal(base_size = 14) +

      labs(
        title = "Tm Distribution Across Datasets",
        x = NULL,
        y = "Melting Temperature (Tm)"
      )

ggsave(
  filename = file,
  plot = p,
  width = 16,
  height = 8,
  dpi = 300
)
  }
)

  })
}


#####################################################################################
#protein plots module
protplotsPanelModuleUI <- function(id) {

  ns <- NS(id)

  tagList(

    div(

      style = "padding-top: 0px; margin-top: 0px;",

      tags$div(

        style = "margin-bottom: 5px;",

        tags$span(
          "Protein Selection",
          tags$span(
            `data-toggle` = "tooltip",
            `data-placement` = "right",
            `data-html` = "true",
            title = HTML(
              "<b>Select proteins for visualization</b><br><br>

              Select one or more proteins from the current condition table and click
              <b>Plot selected proteins</b> to visualize their melting profiles.<br><br>

              • Proteins can be selected independently within each condition.<br>

              • Switching to another condition will show that condition's protein table.<br>

              • <b>Clear selection</b> only removes the currently selected proteins from the active condition table and does not affect selections made in other conditions."
            ),
            tags$i(
              class = "fa fa-question-circle",
              style = "margin-left:8px; color:#337ab7; cursor:pointer;"
            )
          )
        )

      ),

      DTOutput(
        ns("protein_list_table"),
        height = "auto"
      )
    ),

    actionButton(
      ns("clear_selection"),
      "Clear selection"
    ),

    actionButton(
      ns("plot_button"),
      "Plot selected proteins"
    ),

    downloadButton(
      ns("download_table"),
      "Download table"
    ),

    tags$span(
      `data-toggle` = "tooltip",
      `data-placement` = "right",
      `data-html` = "true",
      title = HTML(
        "<b>Download table</b><br><br>

        Export the protein table for the currently active condition.<br><br>

        The export only contains data from the currently displayed condition."
      ),
      tags$i(
        class = "fa fa-question-circle",
        style = "margin-left:8px; color:#337ab7; cursor:pointer;"
      )
    )

  )
}

protplotsPanelModuleServer <- function(id, normalized_data, dataclustered = reactive(NULL)) {
  moduleServer(id, function(input, output, session) {

    ns <- session$ns

#helper function to download filtered data
apply_dt_filters <- function(df, input, table_id) {

  state <- input[[paste0(table_id, "_state")]]

  # no filters applied
  if (is.null(state) || is.null(state$columns)) {
    return(df)
  }

  for (i in seq_along(state$columns)) {

    col_name <- names(df)[i]
    filter_val <- state$columns[[i]]$search$value

    if (!is.null(filter_val) && filter_val != "") {

      if (is.numeric(df[[col_name]])) {

        val <- filter_val

        df <- tryCatch({
          if (grepl("^>=", val)) {
            df[df[[col_name]] >= as.numeric(sub(">=", "", val)), ]
          } else if (grepl("^<=", val)) {
            df[df[[col_name]] <= as.numeric(sub("<=", "", val)), ]
          } else if (grepl("^>", val)) {
            df[df[[col_name]] > as.numeric(sub(">", "", val)), ]
          } else if (grepl("^<", val)) {
            df[df[[col_name]] < as.numeric(sub("<", "", val)), ]
          } else if (grepl("^!=", val)) {
            df[df[[col_name]] != as.numeric(sub("!=", "", val)), ]
          } else {
            df[df[[col_name]] == as.numeric(val), ]
          }
        }, error = function(e) df)

      } else {

        df <- df[grepl(filter_val, df[[col_name]], ignore.case = TRUE), ]

      }
    }

  }

  return(df)
}


    # Choose clustered if available, else normalized


active_df <- reactive({
  dc <- tryCatch(dataclustered(), error = function(e) NULL)
  if (!is.null(dc) && is.data.frame(dc) && nrow(dc) > 0) {

    df <- ensure_protein_df(dc)

    #  detect if we are in comparison mode (presence of suffix)
    if (any(grepl("\\.", names(df)))) {

      #  get suffix from first data column (SAFE and deterministic)
      # col1 = Accession, col2 = Cluster, col3 = first measurement



      first_data_col <- names(df)[3]
      current_suffix <- sub(".*\\.", "", first_data_col)


matched_cols <- names(df)[
  endsWith(names(df), paste0(".", current_suffix))
]

df <- df %>%
  dplyr::select(
    names(df)[1],
    Cluster,
    all_of(matched_cols)
  )

    } else {

      # individual mode (no suffix)
      df <- df %>%
        dplyr::select(
          names(df)[1],
          Cluster,
          dplyr::everything()
        )
    }

    # ensure Cluster is correctly placed
    if ("Cluster" %in% names(df)) {
      df <- df %>% dplyr::relocate(Cluster, .after = 1)
    }

    # remove dataset suffix for display
# remove dataset suffix for display
if (!any(grepl("\\.", names(df)))) {

  df <- df %>%
    dplyr::rename_with(
      ~ sub("\\..*$", "", .x),
      -c(names(df)[1], "Cluster")
    )
}

# ==================================================
# Reorder temperature columns
# ==================================================

id_col <- names(df)[1]

meta_cols <- id_col

if ("Cluster" %in% names(df)) {
  meta_cols <- c(meta_cols, "Cluster")
}

temp_cols <- setdiff(names(df), meta_cols)

temp_info <- data.frame(
  col = temp_cols,
  temp = as.numeric(sub("_.*", "", temp_cols)),
  rep = as.numeric(sub(".*REP", "", temp_cols))
)

temp_info <- temp_info %>%
  dplyr::arrange(temp, rep)

ordered_cols <- c(
  meta_cols,
  temp_info$col
)

df <- df[, ordered_cols, drop = FALSE]

return(df)
  }

  # ✅ fallback to normalized data (individual mode, no clustering yet)
  nd <- tryCatch(normalized_data(), error = function(e) NULL)
  if (!is.null(nd) && is.data.frame(nd) && nrow(nd) > 0) {
    return(ensure_protein_df(nd))
  }

  NULL
})




    # Identify which dataset is active (key used to reset selections when switching)
    dataset_key <- reactive({
      dc <- tryCatch(dataclustered(), error = function(e) NULL)
      if (!is.null(dc) && is.data.frame(dc) && nrow(dc) > 0) "clustered" else "normalized"
    })

    # Store selected rows
    selected_rows <- reactiveVal(NULL)

    # Render table
output$protein_list_table <- DT::renderDT({

  df <- active_df()
  req(df)

  df <- df %>%
    dplyr::mutate(across(where(is.numeric), round, 3))

DT::datatable(
  df,
  filter = "top",
  selection = "multiple",
  options = list(
    pageLength = 5,
    scrollX = TRUE,
    autoWidth = TRUE,
    deferRender = TRUE,
    columnDefs = list(),


    # ✅ THIS IS THE CRITICAL FIX
    stateSave = FALSE,
    dom = "tip"
  ),
  rownames = FALSE
) %>%
    DT::formatRound(columns = names(df)[sapply(df, is.numeric)], digits = 2)

}, server = TRUE)   

    proxy <- DT::dataTableProxy("protein_list_table")

    # When data source switches, clear selection to avoid mismatched row indices
    observeEvent(dataset_key(), {
      selected_rows(NULL)
      selectRows(proxy, NULL)
    }, ignoreInit = TRUE)

    # Restore selection on re-render (if any)
    observeEvent(active_df(), {
      if (!is.null(selected_rows())) selectRows(proxy, selected_rows())
    })


    # Capture new selections
    observeEvent(input$protein_list_table_rows_selected, {
      selected_rows(input$protein_list_table_rows_selected)
      selectRows(proxy, selected_rows()) # avoid blinking
    })

    # Clear button
    observeEvent(input$clear_selection, {
      selected_rows(NULL)
      selectRows(proxy, NULL)
    })

    # Selected data (rows)
    selected_data <- reactive({
      df  <- active_df()
      sel <- selected_rows()
      if (is.null(df) || is.null(sel) || length(sel) == 0) return(NULL)
      df[sel, , drop = FALSE]
    })

    #donwload option
output$download_table <- downloadHandler(

  filename = function() {
    paste0("protein_table_", Sys.Date(), ".csv")
  },

  content = function(file) {

    df <- active_df()

    # ✅ get properly namespaced table id
    table_id <- ns("protein_list_table")

    # ✅ apply correct filtering
    df_filtered <- apply_dt_filters(df, input, table_id)

    write.csv(df_filtered, file, row.names = FALSE)
  }
)
    return(list(
      selected_data = selected_data,
      plot_trigger  = reactive(input$plot_button),
      # (Optional) expose which source is active if you want to show a banner
      dataset_key   = dataset_key
    ))
  })
}



#################################################################################
#Enrichment analysis module
goORAClusterPanelUI <- function(id) {

  ns <- NS(id)

  div(
    style = "padding-top:0px; margin-top:0px;",
    uiOutput(ns("Cluster_OR"))
  )

}
goORAPlotsPanelUI <- function(id) {
  ns <- NS(id)
  fluidRow(
    box(
      title = "Dotplot",
      width = 6,
      status = "info",
      solidHeader = TRUE,

      plotOutput(
        ns("GO_dotplot"),
        height = "600px"
      ) %>% withSpinner(
        color = "#0EA5A5"
      ),

      downloadButton(
        ns("download_go_dotplot"),
        "Download plot"
      )
    ),
box(
  title = "Enrichment plot",
  width = 6,
  status = "success",
  solidHeader = TRUE,

  plotOutput(
    ns("GO_enrichment_plot"),
    height = "600px"
  ) %>% withSpinner(
    color = "#0EA5A5"
  ),

  downloadButton(
    ns("download_go_enrichment_plot"),
    "Download plot"
  )
),
box(
  title = "Upset plot",
  width = 12,
  status = "primary",
  solidHeader = TRUE,

  plotOutput(
    ns("GO_upset_plot"),
    height = "600px"
  ) %>% withSpinner(
    color = "#0EA5A5"
  ),

  downloadButton(
    ns("download_go_upset_plot"),
    "Download plot"
  )
)

  )
}
 

#server
goORAClusterPanelServer <- function(id, Org_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
observeEvent(input$ora_advanced_info, {

cat("\nORA ADVANCED INFO CLICKED\n")
  showModal(

    modalDialog(

      title = "ORA Advanced Options",

      HTML(
        "
        <b>ORA settings</b><br><br>

        These parameters control the Gene Ontology Over-Representation Analysis (ORA)
        performed using the <b>clusterProfiler</b> R package.<br><br>

        <b>Available settings</b><br><br>

        • Multiple-testing correction method.<br>

        • Minimum GO term size.<br>

        • Maximum GO term size.<br>

        • P-value significance threshold.<br>

        • Q-value significance threshold.<br>

        • Optional GO hierarchy expansion using buildGOmap().<br><br>

        GO hierarchy expansion propagates gene annotations to parent GO terms,
        increasing coverage while potentially producing broader enrichment results.<br><br>

        <b>Reference</b><br><br>

        Yu G, Wang LG, Han Y, He QY (2012).<br>

        <i>clusterProfiler: an R package for comparing biological themes among gene clusters.</i><br>

        OMICS 16(5):284-287.<br><br>

        <b>Documentation</b><br><br>

        Copy this URL into your browser:<br><br>

        bioconductor.org/packages/clusterProfiler/
        "
      ),

      easyClose = TRUE,
      size = "l"

    )

  )

})

   output$Cluster_OR <- renderUI({

  tmp_org <- Org_data()



  if (!is.null(tmp_org)) {


  }
advancedOptionsBox <- box(
title = tagList(

  "Advanced options",

  actionLink(
    ns("ora_advanced_info"),

    label = NULL,

    icon = icon("book-open"),

    style = "
      margin-left:8px;
      color:var(--text-secondary);
      display:inline;
    "
  )

),
  width = 12,
  status = "warning",
  solidHeader = TRUE,
  collapsible = TRUE,
  collapsed = TRUE,

  selectInput(
    ns("PAdjust_method"),
    "p-value adjustment method:",
    choices = c(
      "BH",
      "holm",
      "hochberg",
      "hommel",
      "bonferroni",
      "BY",
      "fdr",
      "none"
    )
  ),

  numericInput(
    ns("minGSSize"),
    "Minimum gene set size:",
    value = 10
  ),

  numericInput(
    ns("maxGSSize"),
    "Maximum gene set size:",
    value = 500
  ),

  numericInput(
    ns("P_value_cutoff"),
    "p-value cutoff",
    value = 0.05,
    min = 0,
    max = 1,
    step = 0.01
  ),

  numericInput(
    ns("Qvalue_cutoff"),
    "q-value cutoff:",
    value = 0.2,
    min = 0,
    max = 1,
    step = 0.01
  ),

  checkboxInput(
    ns("use_buildGOmap"),
    "Expand GO hierarchy (buildGOmap)",
    value = TRUE
  )
)

        leftBox <- if (is.null(Org_data())) {

          box(
            title = "Fetch Uniprot data",
            width = 4,
            status = "info",
            solidHeader = TRUE,

            actionButton(
              ns("fetch_data"),
              label = "Fetch data"
            ),

            uiOutput(ns("dynamic_cluster_selector")),

            advancedOptionsBox
          )

        } else {

            leftBox <- box(
              title = "Select cluster for ORA",
              width = 4,
              status = "info",
              solidHeader = TRUE,

              uiOutput(ns("dynamic_cluster_selector")),

              br(),

              div(
                style = "
                  margin-bottom:10px;
                  padding:10px;
                  background-color:#f8f9fa;
                  border-left:4px solid #337ab7;
                ",

                strong("Required: "),
                "Run the Clustering section first to generate protein clusters for enrichment analysis."
              ),
              div(
            style = "
              margin-bottom:10px;
              padding:10px;
              background-color:#f8f9fa;
              border-left:4px solid #5cb85c;
            ",

            strong("Recommended: "),
            "The protein Accession column should contain UniProt accession IDs. ORA relies on UniProt GO annotations and may return incomplete results for other identifier formats."
          ),


              actionButton(
                ns("run_ora"),
                "Run ORA"
              ),

              tags$div(
                style = "height:40px;"
              ),

              tags$hr(
              style = "margin-top:15px; margin-bottom:15px;"
            ),

              advancedOptionsBox
            )
        }
        rightBox <- box(
          title = "ORA result table",
          width = 8,
          status = "info",
          solidHeader = TRUE,

          DT::dataTableOutput(ns("GO_results_table")),

          br(),

          downloadButton(
            ns("download_go_results_table"),
            "Download table"
          ),

          br(),
          br(),

          textOutput(ns("GO_results_note")),


tags$div(
  style = "height:100%;"
)

        )

fluidRow(
  leftBox,
  rightBox
)

    })

    # Return reactive values for use in other modules
ora_settings <- eventReactive(input$run_ora, {

  list(
    P_value_cutoff = input$P_value_cutoff,
    PAdjust_method = input$PAdjust_method,
    minGSSize = input$minGSSize,
    maxGSSize = input$maxGSSize,
    Qvalue_cutoff = input$Qvalue_cutoff,
    use_buildGOmap = input$use_buildGOmap
  )

})
current_go_settings <- reactive({

  list(
    use_buildGOmap = isTRUE(input$use_buildGOmap)
  )

})


return(list(

  settings = ora_settings,

  current_go_settings = current_go_settings,

  run_ora = reactive(input$run_ora)

))
  })
}
.term2gene_cache <- new.env(parent = emptyenv())
goORAPlotsPanelServer <- function(
  id,
  dataclustered,
  Org_data,
  num_clusters,
  cluster_inputs
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

# Cache object
go_cache <- reactiveValues()

# Fetch the data
observeEvent(input$fetch_data, {
  if (is.null(go_cache[["all_clusters"]])) {
    withProgress(message = "Fetching GO data for all clusters", value = 0, {
      data <- dataclustered()
      accessions <- dplyr::pull(data, 1)

      get_go_info_in_batches <- function(accessions, batch_size = 50) {
        all_results <- list()
        total_batches <- ceiling(length(accessions) / batch_size)

        for (i in seq(1, length(accessions), by = batch_size)) {
          batch <- accessions[i:min(i + batch_size - 1, length(accessions))]
          batch_number <- ceiling(i / batch_size)

          withProgress(message = paste("Processing batch", batch_number), value = 0, {
            tryCatch({
              result <- GetProteinGOInfoM(
                batch,
                progress_callback = function(p) incProgress(p)
              )
              all_results <- c(all_results, list(result))
            }, error = function(e) {
              cat("Error in batch", batch_number, ":", conditionMessage(e), "\n")
            })
          })

          incProgress(1 / total_batches)
          Sys.sleep(2)
        }

        do.call(rbind, all_results)
      }

      go_result <- get_go_info_in_batches(accessions)
      go_cache[["all_clusters"]] <- go_result
    })
  }
})



# --- Prepare orgdata ---
orgdata <- reactive({



tmp_org <- Org_data()



  if (!is.null(tmp_org)) {



    d <- tmp_org
cat("\nORGDATA INPUT DIM:\n")
print(dim(d))

cat("\nORGDATA INPUT COLS:\n")
print(names(d))


  } else {



    req(go_cache[["all_clusters"]])

    d <- go_cache[["all_clusters"]]
    d$Entry <- rownames(d)
  }





  x <- d %>%
    dplyr::select(
      Entry,
      Gene.Ontology..biological.process.,
      Gene.Ontology..cellular.component.,
      Gene.Ontology..molecular.function.
    )

cat("\n====================\n")
cat("RAW UNIPROT GO DATA\n")
cat("====================\n")

print(head(x, 5))

  x <- x %>%
    setNames(c("gene ID", "BP", "CC", "MF"))


  x <- x %>%
    pivot_longer(
      cols = c("BP", "CC", "MF"),
      names_to = "Group",
      values_to = "GO"
    )



  x <- x %>%
    filter(!is.na(GO), GO != "")



  x <- x %>%
    separate_rows(GO, sep = ";")



  x <- x %>%
    mutate(
      `GO ID` =
        str_extract(
          GO,
          "\\[GO:\\d+\\]"
        ) %>%
        str_remove_all("\\[|\\]"),

      GO_Description =
        str_trim(
          str_remove(
            GO,
            "\\[GO:\\d+\\]"
          )
        )
    )



  x <- x %>%
    filter(!is.na(`GO ID`))

cat("\n====================\n")
cat("ORGDATA DEBUG\n")
cat("====================\n")

cat("Rows after GO extraction:\n")
print(nrow(x))

cat("\nFirst rows:\n")
print(head(x))

cat("\nFirst GO column values:\n")
print(head(x$GO))

cat("\nFirst extracted GO IDs:\n")
print(head(x$`GO ID`))

x <- x %>%
  dplyr::select(
    `gene ID`,
    `GO ID`,
    GO_Description
  )


  return(x)

})

    # --- Prepare term2gene data ---

term2gene_data <- reactive({





tmp <- orgdata()

org_key <- paste(
  sort(unique(
    if ("Organism" %in% names(Org_data()))
      Org_data()$Organism
    else
      "UNKNOWN"
  )),
  collapse = "_"
)

print(cluster_inputs$current_go_settings())

use_buildGOmap <- isTRUE(
  cluster_inputs$current_go_settings()$use_buildGOmap
)

cache_key <- paste0(
  org_key,
  "_",
  use_buildGOmap
)



if (exists(cache_key, envir = .term2gene_cache)) {



  return(
get(
  cache_key,
  envir = .term2gene_cache
)
  )
}

df <- tmp %>%
  rename(
    GeneID = `gene ID`,
    GO_ID = `GO ID`
  )



df <- df %>%
  mutate(
    GO_ID = str_trim(as.character(GO_ID))
  )



df <- df %>%
  dplyr::select(GO_ID, GeneID)



df <- df %>%
  distinct()



df_base <- as.data.frame(df)


use_buildGOmap <- isTRUE(
  cluster_inputs$current_go_settings()$use_buildGOmap
)

if (use_buildGOmap) {



 print("ABOUT TO CALL buildGOmap")

print(dim(df_base))

print(head(df_base))

res <- buildGOmap(df_base)

} else {


  res <- df_base

}



assign(
  cache_key,
  res,
  envir = .term2gene_cache
)

res

})


    # --- Enrichment results ---
enrichment_results <- reactive({

  req(cluster_inputs$run_ora() > 0)

  withProgress(

    message = "Running GO enrichment analysis",
    detail = "Preparing clusters...",
    value = 0,

    {

  req(
        term2gene_data(),
        dataclustered(),
        cluster_inputs$settings()$P_value_cutoff,
        cluster_inputs$settings()$Qvalue_cutoff,
        cluster_inputs$settings()$PAdjust_method,
        cluster_inputs$settings()$minGSSize,
        cluster_inputs$settings()$maxGSSize
      )
      message("DATASET:")
print(head(dataclustered()))

message("COLNAMES:")
print(colnames(dataclustered()))

      Test_data <- dataclustered() %>% dplyr::select(Accession = 1, Cluster)
      Test_data$Accession <- clean_protein_ids(
  Test_data$Accession
)

cat(
  "\nBEFORE CLEANING EXAMPLE:\n"
)
print(
  head(dataclustered()[[1]], 10)
)

cat(
  "\nAFTER CLEANING EXAMPLE:\n"
)
print(
  head(Test_data$Accession, 10)
)

print(names(Test_data))

str(Test_data)
      message("ENRICHMENT A")
      message("ENRICHMENT B")

cat(
  "\nBackground genes:\n",
  length(unique(Test_data$Accession)),
  "\n"
)
 term2g <- term2gene_data()
cat(
  "Overlap with TERM2GENE:\n",
  sum(
    unique(Test_data$Accession) %in%
      unique(term2g$GeneID)
  ),
  "\n"
)

     
      message("ENRICHMENT C")
      background_genes <- unique(Test_data$Accession)
message("ENRICHMENT D")

cluster_list <- Test_data %>%
  group_by(Cluster) %>%
  summarise(
    Genes = list(Accession),
    .groups = "drop"
  ) %>%
  deframe()
total_clusters <- length(cluster_list)
message("ENRICHMENT E")

cat(
  "Clusters:",
  length(cluster_list),
  "\n"
)

      GOMF <- as.list(GOMFANCESTOR)
      GOCC <- as.list(GOCCANCESTOR)
      GOBP <- as.list(GOBPANCESTOR)

      get_GO_term_name <- function(go_id) {

        obj <- GOTERM[[go_id]]

        if (is.null(obj)) {
          return(NA_character_)
        }

        Term(obj)

      }

get_ontology <- function(go_id) {

  if (!is.null(GOMF[[go_id]]))
    return("MF")

  if (!is.null(GOCC[[go_id]]))
    return("CC")

  if (!is.null(GOBP[[go_id]]))
    return("BP")

  NA_character_

}

message("ENRICHMENT F")

results <- lapply(
  seq_along(cluster_list),
  function(idx) {

    cluster_name <- names(cluster_list)[idx]

    incProgress(
      1 / total_clusters,
      detail = paste(
        "Cluster",
        idx,
        "of",
        total_clusters
      )
    )

  message(
    paste(
      "Running cluster",
      cluster_name
    )
  )
        result <- tryCatch({
          enricher(
            gene = cluster_list[[cluster_name]],
            universe = background_genes,
            TERM2GENE = term2g,
            pvalueCutoff = cluster_inputs$settings()$P_value_cutoff,
            qvalueCutoff = cluster_inputs$settings()$Qvalue_cutoff,
            pAdjustMethod = cluster_inputs$settings()$PAdjust_method,
            minGSSize = cluster_inputs$settings()$minGSSize,
            maxGSSize = cluster_inputs$settings()$maxGSSize
          )
        }, error = function(e) {
          message("Error in enrichment for cluster ", cluster_name, ": ", e$message)
          return(NULL)
        })

        if (!is.null(result) && nrow(result@result) > 0) {
          result@result$Cluster <- cluster_name
          cat(
  "\nUnique GO terms:",
  length(unique(result@result$ID)),
  "\n"
)

cat(
  "GO terms missing from GOTERM:",
  sum(
    sapply(
      result@result$ID,
      function(x) is.null(GOTERM[[x]])
    )
  ),
  "\n"
)

          result@result$Description <- sapply(result@result$ID, get_GO_term_name)
          result@result$ONTOLOGY <- sapply(result@result$ID, get_ontology)
          return(result)
        } else {
          message("No enrichment for cluster ", cluster_name)
          return(NULL)
        }
}
)
results
}
)
})






    # --- Outputs ---
output$GO_results_table <- DT::renderDataTable({

  req(
    enrichment_results(),
    input$Select_cluster_ORA
  )

  Cluster_selected <- as.numeric(
    gsub(
      "Cluster ",
      "",
      input$Select_cluster_ORA
    )
  )

  result <- enrichment_results()[[Cluster_selected]]

  req(
    !is.null(result),
    !is.null(result@result)
  )

  tbl <- result@result %>%
    arrange(p.adjust) %>%
    dplyr::select(
      -c(ID, Cluster)
    )

  numeric_cols <- names(tbl)[
    sapply(tbl, is.numeric)
  ]

  DT::datatable(
    tbl,

    filter = "top",

    options = list(
      pageLength = 3,
      scrollX = TRUE
    )
  ) %>%
    DT::formatRound(
      columns = numeric_cols,
      digits = 1
    )

})

dotplot_reactive <- reactive({

  req(
    enrichment_results(),
    input$Select_cluster_ORA
  )

  Cluster_selected <- as.numeric(
    gsub(
      "Cluster ",
      "",
      input$Select_cluster_ORA
    )
  )

  result <- enrichment_results()[[Cluster_selected]]

  req(
    !is.null(result),
    !is.null(result@result)
  )

  if (
    !any(result@result$p.adjust < 0.05)
  ) {
    return(NULL)
  }

  dotplot(result)

})
    output$GO_dotplot <- renderPlot({
      req(enrichment_results(), input$Select_cluster_ORA)
      Cluster_selected <- as.numeric(gsub("Cluster ", "", input$Select_cluster_ORA))
      result <- enrichment_results()[[Cluster_selected]]
      if (is.null(result) || is.null(result@result) || !any(result@result$p.adjust < 0.05)) {
        plot.new(); text(0.5, 0.5, "No enriched terms found", cex = 1.5)
      } else {
        print(dotplot_reactive())
      }
    })
emapplot_reactive <- reactive({

  req(
    enrichment_results(),
    input$Select_cluster_ORA
  )

  Cluster_selected <- as.numeric(
    gsub(
      "Cluster ",
      "",
      input$Select_cluster_ORA
    )
  )

  result <- enrichment_results()[[Cluster_selected]]

  req(
    !is.null(result),
    !is.null(result@result)
  )

  sig_terms <- result@result %>%
    dplyr::filter(p.adjust < 0.05)

  if (nrow(sig_terms) < 3) {

    return(NULL)

  }

  result <- pairwise_termsim(result)

  tryCatch({

    emapplot(result)

  }, error = function(e) {

    NULL

  })

})
output$GO_enrichment_plot <- renderPlot({

  req(
    enrichment_results(),
    input$Select_cluster_ORA
  )

  Cluster_selected <- as.numeric(
    gsub(
      "Cluster ",
      "",
      input$Select_cluster_ORA
    )
  )

  result <- enrichment_results()[[Cluster_selected]]

  if (
    is.null(result) ||
    is.null(result@result) ||
    !any(result@result$p.adjust < 0.05)
  ) {

    plot.new()

    text(
      0.5,
      0.5,
      "No enriched terms found",
      cex = 1.5
    )

  } else {

    p <- emapplot_reactive()

    if (is.null(p)) {

      plot.new()

      text(
        0.5,
        0.5,
        "Enrichment map requires at least 3 significant GO terms.",
        cex = 1.2
      )

    } else {

      print(p)

    }

  }

})

upsetplot_reactive <- reactive({

  req(
    enrichment_results(),
    input$Select_cluster_ORA
  )

  Cluster_selected <- as.numeric(
    gsub(
      "Cluster ",
      "",
      input$Select_cluster_ORA
    )
  )

  result <- enrichment_results()[[Cluster_selected]]

  req(
    !is.null(result),
    !is.null(result@result)
  )

  if (
    !any(result@result$p.adjust < 0.05)
  ) {
    return(NULL)
  }

  upsetplot(result)


})
    output$GO_upset_plot <- renderPlot({
      req(enrichment_results(), input$Select_cluster_ORA)
      Cluster_selected <- as.numeric(gsub("Cluster ", "", input$Select_cluster_ORA))
      result <- enrichment_results()[[Cluster_selected]]
      if (is.null(result) || is.null(result@result) || !any(result@result$p.adjust < 0.05)) {
        plot.new(); text(0.5, 0.5, "No enriched terms found", cex = 1.5)
      } else {
        upsetplot(result)
      }
    })

    output$GO_results_note <- renderText({
      req(enrichment_results(), input$Select_cluster_ORA)
      Cluster_selected <- as.numeric(gsub("Cluster ", "", input$Select_cluster_ORA))
      result <- enrichment_results()[[Cluster_selected]]
      summary_line <- capture.output(show(result))
      clean_line <- gsub("^#\\.\\.\\.", "", summary_line[grep("enriched terms found", summary_line)])
      return(clean_line)
    })

    output$dynamic_cluster_selector <- renderUI({
      cluster_choices <- if (!is.null(num_clusters())) {
        paste("Cluster", seq_len(num_clusters()))
      } else {
        "Cluster 1"
      }

      selectInput(ns("Select_cluster_ORA"),
                  "Select cluster to analyze",
                  choices = cluster_choices,
                  selected = cluster_choices[1])
    })
output$download_go_results_table <- downloadHandler(

  filename = function() {
    paste0(
      "ORA_Table_",
      Sys.Date(),
      ".csv"
    )
  },

  content = function(file) {

    Cluster_selected <- as.numeric(
      gsub(
        "Cluster ",
        "",
        input$Select_cluster_ORA
      )
    )

    result <- enrichment_results()[[Cluster_selected]]

    write.csv(
      result@result,
      file,
      row.names = FALSE
    )

  }

)
output$download_go_dotplot <- downloadHandler(

  filename = function() {
    paste0(
      "ORA_Dotplot_",
      Sys.Date(),
      ".png"
    )
  },

  content = function(file) {

ggsave(
  file,
  plot = dotplot_reactive(),
  width = 12,
  height = 8,
  dpi = 300
)


  }

)
output$download_go_enrichment_plot <- downloadHandler(

  filename = function() {
    paste0(
      "ORA_Enrichment_Map_",
      Sys.Date(),
      ".png"
    )
  },

  content = function(file) {

ggsave(
  file,
  plot = emapplot_reactive() +
    theme(
      panel.background = element_rect(
        fill = "white",
        colour = "white"
      ),
      plot.background = element_rect(
        fill = "white",
        colour = "white"
      )
    ),
  width = 12,
  height = 8,
  dpi = 300,
  bg = "white"
)

  }

)
output$download_go_upset_plot <- downloadHandler(

  filename = function() {
    paste0(
      "ORA_Upset_Plot_",
      Sys.Date(),
      ".png"
    )
  },

  content = function(file) {

png(
  file,
  width = 4000,
  height = 3000,
  res = 300
)

print(
  upsetplot_reactive()
)

dev.off()

    dev.off()

  }

)
    return(list(
  term2gene_data = term2gene_data,
  orgdata = orgdata
))
  })
}





#GSEA module
goGSEAPanelModuleUI <- function(id) {
  ns <- NS(id)
tagList(
  div(
    style = "padding-top:0px; margin-top:0px;",
    fluidRow(
      box(title = "Control parameters", width = 6, status = "primary", solidHeader = TRUE,
          selectInput(
          ns("rank_by"),
          "Rank genes by:",
          choices = c(
          "t-statistic (recommended)" = "t",
          "logFC (requires normalized/log-transformed data)" = "logFC"
              )
            ),
          div(
            style = "
              margin-bottom:10px;
              padding:10px;
              background-color:#f8f9fa;
              border-left:4px solid #337ab7;
            ",

            strong("Required: "),
            "Run the Limma section first. GSEA uses the differential expression statistics generated by Limma to rank proteins."
          ),
          div(
          style = "
            margin-bottom:10px;
            padding:10px;
            background-color:#f8f9fa;
            border-left:4px solid #5cb85c;
          ",

          strong("Recommended: "),
          "The protein Accession column should contain UniProt accession IDs. GSEA relies on UniProt GO annotations and may return incomplete results when alternative protein identifiers are used."
        ),


          actionButton(
            ns("run_gsea"),
            "Run GSEA"
          ),

            tags$hr(
              style = "margin-top:15px; margin-bottom:15px;"
            ),

box(

  title = tagList(

    "Advanced options",

actionLink(
  ns("gsea_advanced_info"),

  label = NULL,

    icon = icon("book-open"),

    style = "
      margin-left:8px;
      color:var(--text-secondary);
      display:inline;
    "
)
  ),

  width = 12,
  status = "warning",
  solidHeader = TRUE,
  collapsible = TRUE,
  collapsed = TRUE,
            selectInput(
              inputId = ns("gsea_pAdjust_method"),
              label = "p-value adjustment method:",
              choices = c("BH", "holm", "hochberg", "hommel", "bonferroni", "BY", "fdr", "none")
            ),
            numericInput(ns("gsea_minGSSize"), "Minimum gene set size:", value = 5),
            numericInput(ns("gsea_maxGSSize"), "Maximum gene set size:", value = 500),
            numericInput(ns("gsea_pvalue_cutoff"), "p-value cutoff", value = 0.05, min = 0, max = 1, step = 0.01),
            numericInput(ns("gsea_exponent"), "Exponent:", value = 1, min = 0, max = 5, step = 0.1),
            numericInput(ns("permutations"), "Number of permutations:", value = 10000, min = 10000, max = 500000, step = 10000),
            checkboxInput(
            ns("gsea_buildGOmap"),
            "Expand GO hierarchy (buildGOmap)",
            value = TRUE
            ) 
          )
      )
    ),
    fluidRow(
    box(title = "Upset plot", width = 12, status = "info", solidHeader = TRUE,

        plotOutput(
          ns("upsetplot2"),
          height = "600px"
        ) %>% withSpinner(
          color = "#0EA5A5"
        ),

        downloadButton(
          ns("download_upsetplot2"),
          "Download plot"
        )

    )
,
box(title = "Ridge line plot", width = 12, status = "success", solidHeader = TRUE,

    plotOutput(
      ns("ridgeplot"),
      height = "800px"
    ) %>% withSpinner(
      color = "#0EA5A5"
    ),

    downloadButton(
      ns("download_ridgeplot"),
      "Download plot"
    )

),
      box(
        title = "Running score and preranked list plot",
        width = 12,
        status = "primary",
        solidHeader = TRUE,

        uiOutput(ns("gsea_term_ui")),

        br(),

plotOutput(
  ns("gseaplot2"),
  height = "600px"
) %>% withSpinner(
  color = "#0EA5A5"
),

downloadButton(
  ns("download_gseaplot2"),
  "Download plot"
)
      )
    )
  )
)
}


goGSEAPanelModuleServer <- function(
  id,
  hypothesis_result,
  term2gene_data,
  orgdata
) {
  moduleServer(id, function(input, output, session) {
        message(sprintf("[DEBUG][%s] goGSEAPanelModuleServer initialized; ns prefix = %s",
                    id, session$ns("")))
                        observe({
      # show whether key inputs exist in this module instance
      message(sprintf("[DEBUG][%s] Inputs available: run_gsea present? %s, rank_by present? %s",
                      id,
                      !is.null(input$run_gsea),
                      !is.null(input$rank_by)))
      # run once
      isolate(NULL)
      invalidateLater(Inf) # prevents repeated firing
    })
    ns <- session$ns

    observeEvent(input$gsea_advanced_info, {

  showModal(

    modalDialog(

      title = "GSEA Advanced Options",

      HTML(
        "
        <b>Gene Set Enrichment Analysis (GSEA)</b><br><br>

        GSEA evaluates whether predefined Gene Ontology gene sets are enriched toward the top or bottom of a ranked protein list rather than relying on significance thresholds.<br><br>

        In this application, proteins are ranked using the differential analysis results generated in the Limma section.<br><br>

        <b>Advanced parameters</b><br><br>

        • P-value adjustment method.<br>

        • Minimum gene set size.<br>

        • Maximum gene set size.<br>

        • P-value cutoff.<br>

        • Enrichment exponent (weighting parameter).<br>

        • Number of permutations.<br>

        • Optional GO hierarchy expansion using buildGOmap().<br><br>

        Increasing the number of permutations generally improves statistical stability but increases computation time.<br><br>

        <b>Reference</b><br><br>

        Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu X, Liu S, Bo X, Yu G (2021).<br><br>

        <i>clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.</i><br>

        The Innovation 2(3):100141.<br><br>

        "
      ),

      easyClose = TRUE,
      size = "l"

    )

  )

})

    #####term2gene
gsea_term2gene <- reactive({

  withProgress(

    message = "Preparing GO term mappings",
    value = 0,

    {

      df <- orgdata()

  df <- df %>%
    rename(
      GeneID = `gene ID`,
      GO_ID = `GO ID`
    ) %>%
    distinct()


gsea_term2name <- df %>%
  dplyr::select(
    GO_ID,
    GO_Description
  ) %>%
  distinct() %>%
  dplyr::filter(
    !is.na(GO_ID),
    !is.na(GO_Description),
    nzchar(GO_ID),
    nzchar(GO_Description)
  )

df_term2gene <- df %>%
  dplyr::select(
    GO_ID,
    GeneID
  ) %>%
  distinct()

attr(df, "TERM2NAME") <- gsea_term2name
if (isTRUE(input$gsea_buildGOmap)) {

  incProgress(
    0.5,
    detail = "Expanding GO hierarchy"
  )

res <- buildGOmap(
  as.data.frame(df_term2gene)
)
attr(res, "TERM2NAME") <- gsea_term2name
incProgress(
  1,
  detail = "GO hierarchy ready"
)

res

  } else {

res <- as.data.frame(df_term2gene)

incProgress(
  1,
  detail = "GO mappings ready"
)

res


  }

    }
  )

})
    # --- Ranked list for GSEA ---
    ranked_list <- reactive({
      req(hypothesis_result(), input$rank_by)
      message("[DEBUG] Building ranked list...")

      results <- hypothesis_result()
      cat(
  "\n=== VOLCANO RESULTS DEBUG ===\n",
  "Rows:", nrow(results),
  "\nColumns:\n"
)



cat("\nlogFC summary:\n")
print(summary(results$logFC))

cat("\nt summary:\n")
print(summary(results$t))
      rank_col <- input$rank_by
if (
  rank_col == "logFC" &&
  max(abs(results$logFC), na.rm = TRUE) > 100
) {

  showNotification(
    "logFC values appear to be on the raw intensity scale. Use t-statistic ranking or normalize/log-transform the data first.",
    type = "warning",
    duration = 10
  )

}
      if (!rank_col %in% colnames(results)) {
        message("[ERROR] rank_by column not found: ", rank_col)
        return(NULL)
      }

      results <- tibble::rownames_to_column(results, var = "proteins")
      ranked_vector <- results[[rank_col]]
cat(
  "\nGSEA DEBUG:",
  id,
  "\nrank_col =", rank_col,
  "\nrange =",
  paste(range(ranked_vector, na.rm = TRUE), collapse = " -> "),
  "\nNAs =",
  sum(is.na(ranked_vector)),
  "\n"
)
      # Validate
      if (!is.numeric(ranked_vector)) {
        message("[ERROR] rank column must be numeric.")
        return(NULL)
      }

      names(ranked_vector) <- results$proteins
      ranked_vector <- ranked_vector[!is.na(ranked_vector)]

      message("[DEBUG] Ranked list ready: ", length(ranked_vector), " entries.")
      sort(ranked_vector, decreasing = TRUE)
    })

    # --- Run GSEA ---
gsea_result <- eventReactive(input$run_gsea, {

  withProgress(

    message = "Running Gene Set Enrichment Analysis",
    value = 0,

    {
  message("[DEBUG][", id, "] GSEA event triggered")

  # Safely isolate reactivity
  ranked_vals <- isolate(ranked_list())
  term2gene_vals <- isolate(gsea_term2gene())

incProgress(
  0.1,
  detail = "Preparing ranked gene list"
)

  # --- Print diagnostics safely ---
  message("[DEBUG][", id, "] Checking inputs before GSEA...")
  if (is.null(ranked_vals)) {
    message("[ERROR][", id, "] ranked_list() returned NULL.")
  } else {

    print(utils::head(ranked_vals))
  }

  if (is.null(term2gene_vals)) {
    message("[ERROR][", id, "] term2gene_data() returned NULL.")
  } else {
    message("[DEBUG][", id, "] term2gene_data() rows: ", nrow(term2gene_vals))
    message("[DEBUG][", id, "] term2gene_data() columns: ", paste(names(term2gene_vals), collapse = ", "))
  }

  # --- Run GSEA only if both are valid ---
req(ranked_vals, term2gene_vals)

incProgress(
  0.3,
  detail = "Preparing GO term mappings"
)

message("[DEBUG][", id, "] Starting GSEA...")

  gse <- tryCatch({
    message("[DEBUG][", id, "] Entered GSEA call")
    incProgress(
  0.6,
  detail = "Running enrichment permutations"
)
    sum(is.na(term2gene_vals$GO_ID))

sum(term2gene_vals$GO_ID == "NA", na.rm = TRUE)

head(term2gene_vals)
tail(term2gene_vals)
term2name_vals <- attr(
  term2gene_vals,
  "TERM2NAME"
)
res <- GSEA(
  geneList      = ranked_vals,
  TERM2GENE     = term2gene_vals,
  TERM2NAME     = term2name_vals,
  pvalueCutoff  = as.numeric(input$gsea_pvalue_cutoff),
  pAdjustMethod = input$gsea_pAdjust_method,
  minGSSize     = as.numeric(input$gsea_minGSSize),
  maxGSSize     = as.numeric(input$gsea_maxGSSize),
  verbose       = TRUE,
  exponent      = as.numeric(input$gsea_exponent),
  nPermSimple   = as.numeric(input$permutations)
)

incProgress(
  0.95,
  detail = "Finalizing results"
)
message("[DEBUG][", id, "] Finished GSEA call")
if (!is.null(res@result)) {

  cat("\n====================\n")
  cat("GSEA DESCRIPTION CHECK\n")
  cat("====================\n")

  print(
    head(
      res@result[, c("ID", "Description")]
    )
  )

}
message(
  "[DEBUG][", id, "] Class = ",
  paste(class(res), collapse = ", ")
)

message(
  "[DEBUG][", id, "] Is NULL = ",
  is.null(res)
)

if (!is.null(res@result)) {

  message(
    "[DEBUG][", id, "] Result rows = ",
    nrow(res@result)
  )

}

res


  }, error = function(e) {
    message("[ERROR][", id, "] GSEA failed: ", e$message)
    return(NULL)
  })

  if (is.null(gse)) {
    message("[DEBUG][", id, "] GSEA returned NULL.")
    return(NULL)
  }

  message("[DEBUG][", id, "] GSEA output class: ", class(gse))
  if (!is.null(gse@result))
if (!is.null(gse@result)) {

  go_lookup <- orgdata() %>%
    dplyr::select(
      `GO ID`,
      GO_Description
    ) %>%
    distinct()

  idx <- match(
    gse@result$ID,
    go_lookup$`GO ID`
  )

  desc <- go_lookup$GO_Description[idx]

  gse@result$Description <- ifelse(
    is.na(desc) | desc == "",
    gse@result$ID,
    desc
  )

}
    message("[DEBUG][", id, "] Rows in gse@result: ", nrow(gse@result))

  gse

    }
  )
})

    # --- Helper: safe plotting wrapper ---
    safePlot <- function(plot_expr, error_label) {
      tryCatch({
        plot_expr()
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5,
             paste(error_label, e$message),
             cex = 1.2)
      })
    }

    # --- upsetplot ---
    output$upsetplot2 <- renderPlot({
      result <- gsea_result()

      safePlot(function() {
        if (is.null(result) || is.null(result@result) ||
            !any(result@result$p.adjust < 0.05)) {
          plot.new(); text(0.5, 0.5, "No enriched terms found", cex = 1.5)
        } else {
          message("[DEBUG] Rendering upsetplot...")
          print(upsetplot(result))
        }
      }, "Error in upsetplot:")
    })

    # --- ridgeplot ---
    output$ridgeplot <- renderPlot({
      result <- gsea_result()
      

      safePlot(function() {
        if (is.null(result) || is.null(result@result) ||
            !any(result@result$p.adjust < 0.05)) {
          plot.new(); text(0.5, 0.5, "No enriched terms found", cex = 1.5)
        } else {
          message("[DEBUG] Rendering ridgeplot...")
          print(
          ridgeplot(result) +
            theme(
              axis.text.y = element_text(size = 8)
            )
        )
        }
      }, "Error in ridgeplot:")
    })

    # --- gseaplot2 ---

output$gsea_term_ui <- renderUI({

  result <- gsea_result()

  req(
    !is.null(result),
    !is.null(result@result),
    nrow(result@result) > 0
  )

sig_terms <- result@result[
  result@result$p.adjust < 0.05,
]

sig_terms <- sig_terms |>
  dplyr::filter(
    !is.na(ID),
    nzchar(ID)
  )

term_choices <- setNames(
  sig_terms$ID,
  sig_terms$Description
)

term_choices <- term_choices[
  !is.na(names(term_choices))
]

term_choices <- term_choices[
  nzchar(names(term_choices))
]



  term_choices <- term_choices[
    !is.na(term_choices) &
    nzchar(term_choices)
  ]

selectizeInput(
  inputId = ns("gsea_term"),
  label = "Select term(s) to display:",
  choices = term_choices,
  selected = sig_terms$ID[1],
  multiple = TRUE,
  options = list(
    placeholder = "Search GO terms...",
    maxOptions = 1000
  )
)

})
output$gseaplot2 <- renderPlot({

  req(input$gsea_term)

  result <- gsea_result()

  req(
    !is.null(result),
    !is.null(result@result)
  )

  safePlot(function() {

    valid_terms <- input$gsea_term[
      input$gsea_term %in% result@result$ID
    ]

    req(length(valid_terms) > 0)

    gseaplot2(
      result,
      geneSetID = valid_terms
    )

  }, "Error in gseaplot2:")

})
output$download_upsetplot2 <- downloadHandler(

  filename = function() {
    paste0(
      "GSEA_UpsetPlot_",
      Sys.Date(),
      ".png"
    )
  },

  content = function(file) {

    png(
      file,
      width = 2000,
      height = 1600,
      res = 300
    )

    result <- gsea_result()

    print(
      upsetplot(result)
    )

    dev.off()

  }

)
output$download_ridgeplot <- downloadHandler(

  filename = function() {
    paste0(
      "GSEA_RidgePlot_",
      Sys.Date(),
      ".png"
    )
  },

  content = function(file) {

    png(
      file,
      width = 2000,
      height = 1600,
      res = 300
    )

    result <- gsea_result()

    print(
      ridgeplot(result)
    )

    dev.off()

  }

)
output$download_gseaplot2 <- downloadHandler(

  filename = function() {
    paste0(
      "GSEA_RunningScore_",
      Sys.Date(),
      ".png"
    )
  },

  content = function(file) {

    png(
      file,
      width = 2200,
      height = 1800,
      res = 300
    )

    result <- gsea_result()

    valid_terms <- input$gsea_term[
      input$gsea_term %in%
        result@result$ID
    ]

    print(
      gseaplot2(
        result,
        geneSetID = valid_terms
      )
    )

    dev.off()

  }

)


  }) 
}



cat("========================================\n")
cat("TPP Analyst Dependency Installer\n")
cat("========================================\n\n")
#Update the R version requirement here if needed
required_r <- package_version("4.6.1")

if (getRversion() < required_r) {
  stop(
    sprintf(
    "TPP Analyst requires R >= %s. Current version: %s",
    as.character(required_r),
    as.character(getRversion())
    )
  )
}


cat("[OK] R version:", as.character(getRversion()), "\n")

cran_repo <- "https://cloud.r-project.org"

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = cran_repo)

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = cran_repo)

cran_packages <- c(
  "Amelia",
  "broom",
  "bsicons",
  "bslib",
  "callr",
  "curl",
  "dplyr",
  "DT",
  "ggnewscale",
  "ggplot2",
  "ggtext",
  "ggupset",
  "httr",
  "knitr",
  "later",
  "minpack.lm",
  "pheatmap",
  "plotly",
  "pracma",
  "processx",
  "purrr",
  "RColorBrewer",
  "readr",
  "readxl",
  "shiny",
  "shinybusy",
  "shinycssloaders",
  "shinydashboard",
  "shinyjs",
  "stringr",
  "tibble",
  "tidyr",
  "tidyverse",
  "VIM",
  "ggridges",
  "UniprotR"
)

missing_cran <- cran_packages[
  !sapply(cran_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_cran) > 0) {

  cat("\nInstalling CRAN packages...\n")

  install.packages(
    missing_cran,
    repos = cran_repo,
    dependencies = TRUE
  )

} else {

  cat("\n[OK] All CRAN packages already installed.\n")

}

cat("\nConfiguring Bioconductor...\n")

BiocManager::install(
  version = "3.23",
  ask = FALSE
)

bioc_packages <- c(
  "Biobase",
  "BiocParallel",
  "Biostrings",
  "IRanges",
  "GenomicAlignments",
  "circlize",
  "enrichplot",
  "clusterProfiler",
  "ComplexHeatmap",
  "GO.db",
  "impute",
  "limma",
  "TPP",
  "NPARC"
)

missing_bioc <- bioc_packages[
  !sapply(bioc_packages, requireNamespace, quietly = TRUE)
]

if (length(missing_bioc) > 0) {

  cat("\nInstalling Bioconductor packages...\n")

  BiocManager::install(
    missing_bioc,
    ask = FALSE,
    update = FALSE
  )

} else {

  cat("\n[OK] All Bioconductor packages already installed.\n")

}
#DEP
if (!requireNamespace("DEP", quietly = TRUE)) {

  cat("\nInstalling DEP...\n")

  remotes::install_github(
    "arnesmits/DEP",
    ref = "b425d8d0db67b15df4b8bcf87729ef0bf5800256",
    dependencies = TRUE
  )

} else {

  cat("\n[OK] DEP already installed.\n")

}

test_packages <- c(
  "DEP",
  "TPP",
  "NPARC",
  "ComplexHeatmap",
  "clusterProfiler",
  "shiny"
)

cat("\nValidating installation...\n")

for(pkg in test_packages){

  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )

  cat("[OK]", pkg, "\n")
}


# ------------------------------------------------------------------
# Environment validation
# ------------------------------------------------------------------

cat("\n========================================\n")
cat("Validating software versions\n")
cat("========================================\n")

tested_r <- package_version("4.6.1")
current_r <- getRversion()

if (current_r != tested_r) {

  message(
    sprintf(
      paste(
        "TPP Analyst was developed and tested using R %s.",
        "You are running R %s.",
        "The application may still work correctly,",
        "but some features could behave differently."
      ),
      tested_r,
      current_r
    )
  )

} else {

  cat("[OK] R version:", as.character(current_r), "\n")

}

cat(
  "[OK] Bioconductor version:",
  as.character(BiocManager::version()),
  "\n"
)

required_versions <- c(
  "DEP"               = "1.7.1",
  "TPP"               = "3.40.0",
  "NPARC"             = "1.24.0",
  "ComplexHeatmap"    = "2.28.0",
  "clusterProfiler"   = "4.20.0",
  "shiny"             = "1.14.0",
  "shinydashboard"    = "0.7.3",
  "shinyjs"           = "2.1.1",
  "shinybusy"         = "0.3.3",
  "shinycssloaders"   = "1.1.0",
  "DT"                = "0.34.0",
  "plotly"            = "4.12.1",
  "ggplot2"           = "4.0.3",
  "dplyr"             = "1.2.1",
  "tidyr"             = "1.3.2",
  "tibble"            = "3.3.1",
  "purrr"             = "1.2.2",
  "stringr"           = "1.6.0",
  "readr"             = "2.2.0",
  "readxl"            = "1.4.5",
  "bslib"             = "0.12.0",
  "bsicons"           = "0.1.2",
  "broom"             = "1.0.13",
  "callr"             = "3.8.0",
  "processx"          = "3.9.0",
  "curl"              = "8.0.0",
  "later"             = "1.4.8",
  "pheatmap"          = "1.0.13",
  "pracma"            = "2.4.6",
  "minpack.lm"        = "1.2-4",
  "Amelia"            = "1.8.3",
  "ggnewscale"        = "0.5.2",
  "ggtext"            = "0.1.2",
  "UniprotR"          = "2.5.1",
  "VIM"               = "7.0.0",
  "Biobase"           = "2.72.0",
  "BiocParallel"      = "1.46.0",
  "Biostrings"        = "2.80.1",
  "IRanges"           = "2.46.0",
  "GenomicAlignments" = "1.48.0",
  "circlize"          = "0.4.18",
  "enrichplot"        = "1.32.0",
  "GO.db"             = "3.23.1",
  "impute"            = "1.86.0",
  "limma"             = "3.68.5"
)

cat("\nChecking package versions...\n")

for (pkg in names(required_versions)) {

  if (!requireNamespace(pkg, quietly = TRUE)) {

    warning(
      sprintf(
        "%s is not installed.",
        pkg
      )
    )

    next

  }

  installed_version <- as.character(packageVersion(pkg))
  expected_version  <- required_versions[[pkg]]

  if (installed_version == expected_version) {

    cat(
      sprintf(
        "[OK] %-20s %s\n",
        pkg,
        installed_version
      )
    )

  } else {

    warning(
      sprintf(
        "%s version mismatch. Installed: %s | Tested: %s",
        pkg,
        installed_version,
        expected_version
      )
    )

  }

}

cat("\n")
cat("========================================\n")
cat("TPP Analyst dependencies installed.\n")
cat("You can now launch the app using:\n\n")
cat("shiny::runApp('App/Core')\n")
cat("========================================\n")

# Installation

TPP Analyst can be run in three different ways depending on your needs and experience level.

## Option 1: Run the Pre-Built Docker Image (Recommended)

This is the easiest and most reliable option. The Docker image contains all required R packages, Bioconductor packages, and system dependencies.

```bash
docker run -p 3838:3838 <your-docker-image>
```

Then open:

```text
http://localhost:3838
```

This option guarantees that TPP Analyst runs in the exact environment in which it was developed and tested.

---

## Option 2: Build the Docker Image Locally

If you would like to inspect, customize, or modify the Docker environment, you can build the image directly from the source code.

Clone the repository:

```bash
git clone https://github.com/<your-username>/TPP_Analyst.git
cd TPP_Analyst
```

Build the image:

```bash
docker build -t tppanalyst .
```

Run the container:

```bash
docker run -p 3838:3838 tppanalyst
```

The application will then be available at:

```text
http://localhost:3838
```

This option is recommended for developers and advanced users who wish to customize the application or Docker environment.

---

## Option 3: Run Locally in R

TPP Analyst can also be executed directly from an R session without Docker.

### Requirements

- R ≥ 4.6.1
- Bioconductor 3.22
- Internet connection during the initial package installation

### Installation

Clone the repository:

```bash
git clone https://github.com/<your-username>/TPP_Analyst.git
```

Open R and navigate to the project directory:

```r
setwd("TPP_Analyst")
```

Install all required dependencies:

```r
source("dependencies.R")
```

Launch the application:

```r
shiny::runApp("App/Core")
```

### What the Installation Script Does

The installation script automatically:

- Checks that the installed R version meets the minimum requirement.
- Installs missing CRAN packages.
- Installs required Bioconductor packages.
- Installs the validated DEP version used during development.
- Verifies that critical packages load successfully.
- Reports package version differences relative to the tested environment.

### Notes

Some operating systems may require additional system libraries to compile certain packages. These dependencies are installed automatically when using Docker, but may need to be installed manually for local installations.

Examples include:

- OpenSSL
- libcurl
- libxml2
- HDF5
- NetCDF
- GLPK
- ImageMagick

If you encounter package installation issues, we strongly recommend using one of the Docker-based installation methods above.

---

# Reproducibility

TPP Analyst was developed and validated using the following software environment:

```text
R 4.6.1
Bioconductor 3.22

DEP               1.7.1
TPP               3.40.0
NPARC             1.24.0
ComplexHeatmap    2.28.0
clusterProfiler   4.20.0

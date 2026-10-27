# Installation

TPP Analyst can be run in three different ways depending on your needs and experience level.

## Option 1: Run the Pre-Built Docker Image (Recommended)

TPP Analyst is available as a pre-built Docker image containing all required R packages, Bioconductor packages, and system dependencies.

This is the simplest and most reliable installation method and is recommended for most users.

### Step 1: Install Docker Desktop

Download and install Docker Desktop for your operating system:

https://www.docker.com/products/docker-desktop/

After installation, launch Docker Desktop and ensure it is running before continuing.

### Step 2: Run TPP Analyst

Open a terminal (Command Prompt, PowerShell, Terminal, or Bash) and run:

```bash
docker run -p 3838:3838 knth16/tppanalyst:latest
```

The first time you run this command, Docker will automatically download the image from Docker Hub. The image is approximately 5-6 GB and may take several minutes to download depending on your internet connection.

Docker automatically selects the correct image for your system, including:

- Windows (Intel/AMD)
- Linux (Intel/AMD)
- Intel Macs
- Apple Silicon Macs (M1/M2/M3/M4)

### Step 3: Open the Application

Open your web browser and navigate to:

```text
http://localhost:3838
```

TPP Analyst should now be available.

### Alternative: Run Through Docker Desktop

If you prefer using Docker Desktop instead of the command line:

1. Open Docker Desktop.
2. Navigate to **Images**.
3. Search for:

```text
knth16/tppanalyst:latest
```

4. Pull the image.
5. Click **Run**.
6. Configure the port mapping:

```text
Host Port: 3838
Container Port: 3838
```

7. Start the container.
8. Open:

```text
http://localhost:3838
```

in your browser.

### Notes

- No additional configuration is required.
- The Docker image contains all required software dependencies.
- This option guarantees that TPP Analyst runs in the same environment in which it was developed and tested.
- Docker automatically downloads the appropriate image for Intel/AMD and Apple Silicon processors.

---

## Option 2: Build the Docker Image Locally

If you would like to inspect, customize, or modify the Docker environment, you can build the image directly from the source code.

Clone the repository:

```bash
git clone https://github.com/knth16/TPP_Analyst.git
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
- Internet connection during the initial package installation

### Installation

Clone the repository:

```bash
git clone https://github.com/knth16/TPP_Analyst.git
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
```
# User Manual
## 📖 Documentation

A complete step-by-step user manual for TPP Analyst is available online:

[**Open the TPP Analyst User Manual**](https://knth16.github.io/TPP_Analyst/)



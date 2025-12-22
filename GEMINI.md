# GEMINI.md - sdm-runner-conteiner

## Project Overview

This project provides a containerized R environment for performing Species Distribution Modeling (SDM). It is designed to train SDM models, generate predictions, and evaluate model performance in a reproducible and isolated Docker container.

The core logic is implemented in the `utils.r` file, which contains a suite of functions for the entire SDM workflow, including:

*   **Data Acquisition:** Fetching species occurrence data from GBIF.
*   **Model Training:** Training SDM models using various algorithms (e.g., Random Forest).
*   **Prediction:** Generating prediction rasters from trained models and environmental data.
*   **Evaluation:** Calculating model evaluation statistics and generating reports.
*   **Utilities:** Helper functions for tasks like VIF calculation and data resampling.

The project uses the `rocker/geospatial` Docker image, which comes with many common geospatial libraries pre-installed. The `Dockerfile` adds further R packages required for the SDM analysis, such as `sdm`, `raster`, and `rgbif`.

This project appears to be a template or library, as the main execution scripts (`run.sh`, `r/train_model.r`, and `r/predict_sdm.r`) are empty. The user is expected to write the necessary R code to call the functions in `utils.r` to perform a specific analysis.

## Building and Running

### Building the Docker Image

To build the Docker image, run the following command in the project's root directory:

```bash
docker build -t sdm-runner .
```

### Running the Container

Since `run.sh` (the container's entrypoint) is empty, you need to provide a command to run when starting the container. Here are two common ways to do this:

**1. Run an Interactive R Session**

This is useful for developing and testing your analysis interactively. This command starts the container, mounts the current directory into `/app` in the container, and opens an R session.

```bash
docker run -it --rm -v "$(pwd):/app" sdm-runner R
```

Inside the R session, you can then source the utility functions and run them:

```R
source("utils.r")
# Now you can call functions like generate_model(), predict_sdm(), etc.
```

**2. Run a Custom R Script**

For a non-interactive workflow, you can write an R script (e.g., `my_analysis.r`) and run it with Docker.

First, create your script, for example `my_analysis.r`:

```R
# my_analysis.r
source("utils.r")

# Example: Train a model
generate_model(
  data_path   = "/app/data",
  output_path = "/app/output",
  model_name  = "my_species_rf",
  csv_filename = "predictors.csv"
)

# Example: Make a prediction
predict_sdm(
  model_path      = "/app/output",
  predictors_path = "/app/rasters",
  output_file     = "/app/output/my_species_prediction.tif",
  csv_path        = "/app/data/predictors.csv",
  model_name      = "my_species_rf"
)
```

Then, run the script with Docker. This command mounts the necessary data, output, and raster directories.

```bash
docker run --rm \
  -v "$(pwd)/data:/app/data" \
  -v "$(pwd)/output:/app/output" \
  -v "$(pwd)/rasters:/app/rasters" \
  sdm-runner Rscript my_analysis.r
```

**TODO:** The `run.sh` script should be implemented to define the main execution workflow of the project. This script should call the necessary R scripts to perform the analysis.

## Development Conventions

*   **Language:** The project is written in R.
*   **Code Style:** The code in `utils.r` follows a functional programming paradigm, with clear separation of concerns.
*   **Documentation:** Function comments and documentation are written in Portuguese.
*   **File Structure:**
    *   `utils.r`: Contains reusable utility functions.
    *   `r/`: Intended to hold the main R scripts for training and prediction (currently empty).
    *   `run.sh`: The main entrypoint for the Docker container (currently empty).
*   **Dependencies:** R packages are managed in the `Dockerfile`.

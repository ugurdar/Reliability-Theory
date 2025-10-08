# Script to initialize renv for dependency management in the Reliability-Theory project # nolint

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv") # Install renv if not already installed # nolint

if (!file.exists("renv.lock")) {
  renv::init()
} else {
  message("renv is already initialized for this project.")
}

# Initialize project-local R environment with renv
# renv::init()  # This line is commented out since it's handled conditionally above # nolint
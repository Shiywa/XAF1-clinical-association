# Shared path helpers for the XAF1 figure-code repository.
# Run scripts from the repository root or set XAF1_PROJECT_ROOT explicitly.

xaf1_project_root <- Sys.getenv("XAF1_PROJECT_ROOT", unset = getwd())
xaf1_project_root <- normalizePath(xaf1_project_root, mustWork = FALSE)

xaf1_path <- function(...) {
    file.path(xaf1_project_root, ...)
}

xaf1_output_dir <- function(...) {
    path <- xaf1_path("results", ...)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    path
}

xaf1_env_path <- function(name, default = NULL, must_exist = FALSE) {
    value <- Sys.getenv(name, unset = "")
    if (!nzchar(value)) {
        value <- default
    }
    if (is.null(value) || !nzchar(value)) {
        stop("Environment variable ", name, " is required.")
    }
    value <- normalizePath(value, mustWork = FALSE)
    if (must_exist && !file.exists(value)) {
        stop(name, " does not exist: ", value)
    }
    value
}

xaf1_threads <- function(default = 4L) {
    value <- suppressWarnings(as.integer(Sys.getenv("XAF1_THREADS", unset = default)))
    if (is.na(value) || value < 1L) default else value
}

xaf1_require_file <- function(path, label = "input") {
    if (!file.exists(path)) {
        stop("Missing ", label, ": ", path)
    }
    path
}

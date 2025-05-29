# Helper to load or install a CRAN package
load_or_install_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

# Helper to load or install a GitHub package
load_or_install_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github(repo, quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
}
#' Load or install a CRAN package
#'
#' @param pkg Character. Package name.
#' @return Loads the package. Installs it if not already installed.
#' @export
load_or_install_cran <- function(pkg) { ... }

#' Load or install a GitHub package
#'
#' @param pkg Character. Package name.
#' @param repo Character. GitHub repo string in the form "user/repo".
#' @return Loads the package. Installs it from GitHub if needed.
#' @export
load_or_install_github <- function(pkg, repo) { ... }


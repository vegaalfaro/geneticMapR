#' Load or Install a CRAN Package
#'
#' This function checks if a specified package is installed. If the package is
#' not installed, it installs it from CRAN, then loads it into the R session.
#'
#' @param pkg A character string specifying the name of the package to be loaded or installed.
#'
#' @details
#' The function first checks if the package is available in the current session using `requireNamespace`.
#' If the package is not installed, it installs it using `install.packages` and then loads it using `library`.
#'
#' @return This function does not return a value. It is called for its side effect of loading or installing the package.
#'
#' @examples
#' \dontrun{
#' load_or_install_cran("ggplot2")
#' load_or_install_cran("dplyr")
#' }
#'
#' @export
load_or_install_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}

load_or_install_cran <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
  library(pkg, character.only = TRUE)
}




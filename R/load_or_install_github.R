#' Load or Install a GitHub Package
#'
#' This function checks if a specified package is installed. If the package is not installed,
#' it installs the package from GitHub using the `remotes` package, and then loads it into the R session.
#'
#' @param pkg A character string specifying the name of the package to be loaded or installed.
#' @param repo A character string specifying the GitHub repository in the form "user/repo" from which to install the package.
#'
#' @details
#' The function first checks if the package is available in the current session using `requireNamespace`.
#' If the package is not installed, it installs it from the specified GitHub repository using `remotes::install_github`,
#' and then loads the package using `library`.
#'
#' @return This function does not return a value. It is called for its side effect of loading or installing the package.
#'
#' @examples
#' \dontrun{
#' load_or_install_github("geneticMapR", "vegaalfaro/geneticMapR")
#' load_or_install_github("dplyr", "tidyverse/dplyr")
#' }
#'
#' @export
load_or_install_github <- function(pkg, repo) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      install.packages("devtools")
    }
    devtools::install_github(repo, quiet = TRUE)
  }
  library(pkg, character.only = TRUE)
}

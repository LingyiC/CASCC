#' Check package dependencies
#'
#' @examples
#' checkDependencies()
#' @export





checkDependencies <- function(){
  # if(!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
  # if(!requireNamespace("remotes", quietly = TRUE)){install.packages("remotes")}

  checkVersion <- function(x) {
    tryCatch({
      packageVersion(x)
    }, error = function(e) {
      NA
    })}

  packages <- c("cafr", "Seurat", "SeuratObject")
  versions <- c(checkVersion("cafr"), checkVersion("Seurat"), checkVersion("SeuratObject"))
  # repos <- c(NA, "https://satijalab.r-universe.dev", "https://satijalab.r-universe.dev")

  for (i in seq_along(packages)) {
    package <- packages[i]
    version <- versions[i]
    # repo <- repos[i]

      if (!is.na(version)) {
        cat(paste0(package, " is already installed."), "\n")
      } else {
         cat(paste0(package, " is not installed."), "\n")
      }
  }

  if((packageVersion("Seurat") != "4.3.0")){
  cat(paste0("You are using Seurat(", packageVersion("Seurat"), "). A more stable version of Seurat (4.3.0) is highly recommended to find seeds."), "\n")
  # cat('yes - remotes::install_version("Seurat", "4.3.0")', "\n") Do you want to downgrade? (yes/no)
  remotes::install_version("Seurat", "4.3.0")
  }
}



# checkDependencies <- function(){
#   if(!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
#   if(!requireNamespace("remotes", quietly = TRUE)){install.packages("remotes")}

#   packages <- c("cafr", "Seurat", "SeuratObject", "NbClust", "ggplot2", "stats", "Matrix", "dplyr", "RColorBrewer", "gridExtra", "tidyverse", "parallel")
#   versions <- c(NA, "4.3", "4.1.3", NA, NA, NA, NA, NA, NA, NA, NA, NA)
#   repos <- c(NA, "https://satijalab.r-universe.dev", "https://satijalab.r-universe.dev", NA, NA, NA, NA, NA, NA, NA, NA, NA)

#   for (i in seq_along(packages)) {
#     package <- packages[i]
#     version <- versions[i]
#     repo <- repos[i]

#     if (!requireNamespace(package, quietly = TRUE)) {
#       if (!is.na(version)) {
#         if(package == "Seurat" || package == "SeuratObject"){
#           remotes::install_version(package, version, repos = c(repo, getOption("repos")))
#         }
#       } else {
#         install.packages(package)
#       }
#       cat(paste0(package, " has been successfully installed just now."), "\n")
#     } else {
#       cat(paste0(package, " is already installed."), "\n")
#     }
#   }

#   # Note: The 'cafr' package installation requires a special command.
#   if (!requireNamespace("cafr", quietly = TRUE)) {
#     devtools::install_github("weiyi-bitw/cafr")
#     cat("cafr was installed\n")
#   } else {
#     cat("cafr is already installed\n")
#   }
# }

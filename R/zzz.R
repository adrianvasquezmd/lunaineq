#' @keywords internal
.onLoad <- function(libname, pkgname) {
  conflicted::conflict_prefer("filter",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("select",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("mutate",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("rename",    "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("summarise", "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("arrange",   "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("group_by",  "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("ungroup",   "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("slice",     "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("pull",      "dplyr", quiet = TRUE)
  conflicted::conflict_prefer("lag",       "dplyr", quiet = TRUE)
}

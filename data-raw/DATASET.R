## code to prepare `grembi_parameters` dataset goes here

grembi_parameters <- read.csv("data-raw/grembi_parameters.csv")
usethis::use_data(grembi_parameters, overwrite = TRUE)

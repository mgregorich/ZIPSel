rm(list=ls())

library(quarto)
plan(multisession, workers=5)
vec_propnz <- c(0.9, 0.7, 0.5, 0.3, 0.1)
res <- future_sapply(vec_propnz, function(x) quarto_render(here::here("example", "example_methods.qmd"), 
                                                           output_file = paste0("example_methods_pnz", x, ".html"),
                                                    execute_params = list(propnz=x)), future.seed=TRUE)
plan(sequential)

library(tidymodels)

rand_forest() |> 
  set_engine("ranger") |>
  set_mode("regression") |>
  translate()
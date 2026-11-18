library(m6APrediction)
ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction",
                             mustWork = TRUE))
input <- read.csv(system.file("extdata", "m6A_input_example.csv",
                              package = "m6APrediction", mustWork = TRUE))
out <- prediction_multiple(ml_fit, input, positive_threshold = 0.6)
stopifnot(nrow(out) == 100L, ncol(out) == 9L,
          identical(round(out$predicted_m6A_prob[1:3], 3), c(0.876, 0, 0.498)),
          identical(out$predicted_m6A_status[1:3], c("Positive", "Negative", "Negative")),
          all(out$predicted_m6A_prob >= 0 & out$predicted_m6A_prob <= 1),
          setequal(getNamespaceExports("m6APrediction"),
                   c("prediction_multiple", "prediction_single")))
single <- prediction_single(ml_fit, 0.6, "mRNA", "CDS", 12, 5, 0.8, "GGACA")
stopifnot(is.character(single), length(single) == 2L,
          identical(names(single), c("predicted_m6A_prob", "predicted_m6A_status")),
          identical(unname(single), c("0.51", "Positive")))
# Equality at the cutoff is Negative, with unchanged scores.
equal <- prediction_multiple(ml_fit, input[1, ], out$predicted_m6A_prob[1])
stopifnot(equal$predicted_m6A_status == "Negative")
# Single and batch interfaces must agree for the same independent inputs.
for (i in 1:5) {
  row <- as.list(input[i, ])
  one <- do.call(prediction_single, c(list(ml_fit = ml_fit), row,
                                     list(positive_threshold = 0.6)))
  stopifnot(as.numeric(one[1]) == out$predicted_m6A_prob[i],
            one[2] == out$predicted_m6A_status[i])
}
with_id <- input[1:3, ]; with_id$site_id <- c("a", "b", "c")
stopifnot(identical(prediction_multiple(ml_fit, with_id)$site_id, with_id$site_id))
fails <- function(expr) inherits(try(force(expr), silent = TRUE), "try-error")
stopifnot(fails(prediction_multiple(ml_fit, input[0, ])),
          fails(prediction_multiple(ml_fit, input[, -1])),
          fails(prediction_multiple(ml_fit, input, positive_threshold = NA_real_)))
bad <- input[1, ]; bad$DNA_5mer <- "GGACN"
stopifnot(fails(prediction_multiple(ml_fit, bad)))
bad <- input[1, ]; bad$RNA_type <- "unknown"
stopifnot(fails(prediction_multiple(ml_fit, bad)))
bad <- input[1, ]; bad$gc_content <- Inf
stopifnot(fails(prediction_multiple(ml_fit, bad)))
bad <- input[1, ]; bad$evolutionary_conservation <- 1.1
stopifnot(fails(prediction_multiple(ml_fit, bad)))

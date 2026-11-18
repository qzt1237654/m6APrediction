#' Encode five nucleotide positions as categorical predictors
#'
#' Internal helper used by the two exported prediction functions. Each position
#' is a factor, not a one-hot indicator. The supplied model uses positions
#' 1, 2, 4 and 5; its training contexts have A at position 3.
#'
#' @param dna_strings A nonempty character vector of five-base uppercase DNA
#'   contexts containing only A, T, C and G. Missing values are not allowed.
#' @return A data frame with five factor columns named `nt_pos1` to `nt_pos5`,
#'   each with levels `c("A", "T", "C", "G")`, and one row per context.
#' @keywords internal
dna_encoding <- function(dna_strings) {
  if (!is.character(dna_strings) || !length(dna_strings) ||
      anyNA(dna_strings) || any(!grepl("^[ATCG]{5}$", dna_strings))) {
    stop("DNA_5mer must contain five uppercase A/T/C/G bases per row.",
         call. = FALSE)
  }
  seq_m <- matrix(unlist(strsplit(dna_strings, "", fixed = TRUE)),
                  ncol = 5L, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:5)
  seq_df <- as.data.frame(seq_m, stringsAsFactors = FALSE)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  seq_df
}

#' Predict m6A probability and status for multiple sites
#'
#' Applies a fitted classification random forest to a table of site features.
#' The supplied model and example table are available in `extdata`; use
#' `system.file()` to find them after installation.
#'
#' @param ml_fit A fitted `randomForest` classification model with a
#'   `Positive` response class, such as the bundled `rf_fit.rds`.
#' @param feature_df A nonempty data frame containing these seven columns:
#'   `gc_content` (numeric, 0 to 1), `RNA_type` (mRNA, lincRNA, lncRNA or
#'   pseudogene), `RNA_region` (CDS, intron, 3'UTR or 5'UTR), `exon_length`
#'   (nonnegative numeric on the course log2 scale), `distance_to_junction`
#'   (nonnegative numeric on the course log2 scale),
#'   `evolutionary_conservation` (numeric, 0 to 1), and `DNA_5mer` (five
#'   uppercase A/T/C/G bases). Numeric inputs must be finite, with no missing
#'   values. Extra columns are retained. Intended contexts have central A.
#' @param positive_threshold One finite number between 0 and 1. A score strictly
#'   greater than this cutoff is Positive; a score equal to it is Negative.
#'
#' @return A data frame with the original rows and input columns, followed by
#'   numeric `predicted_m6A_prob` and character `predicted_m6A_status` columns.
#'   `RNA_type` and `RNA_region` are converted to factors with the training
#'   levels. Existing prediction columns are replaced. Probabilities are the
#'   fraction of trees voting Positive, not measured methylation levels.
#' @details No model is trained by this function. Length and distance inputs
#'   must already use the transformations in the teaching data; raw base-pair
#'   counts should not be substituted. The model was developed for educational
#'   exploration of the supplied dataset and is not a clinical diagnostic tool.
#' @import randomForest
#' @importFrom stats predict
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv",
#'                                    package = "m6APrediction", mustWork = TRUE))
#' prediction_multiple(ml_fit, feature_df[1:3, ], positive_threshold = 0.6)
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5) {
  required <- c("gc_content", "RNA_type", "RNA_region", "exon_length",
                "distance_to_junction", "evolutionary_conservation", "DNA_5mer")
  if (!inherits(ml_fit, "randomForest") || ml_fit$type != "classification") {
    stop("ml_fit must be a classification randomForest model.", call. = FALSE)
  }
  if (!is.data.frame(feature_df) || !nrow(feature_df)) {
    stop("feature_df must be a nonempty data frame.", call. = FALSE)
  }
  if (anyDuplicated(names(feature_df))) {
    stop("feature_df must have unique column names.", call. = FALSE)
  }
  missing <- setdiff(required, names(feature_df))
  if (length(missing)) {
    stop(paste("Missing required columns:", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  if (!is.numeric(positive_threshold) || length(positive_threshold) != 1L ||
      !is.finite(positive_threshold) || positive_threshold < 0 ||
      positive_threshold > 1) {
    stop("positive_threshold must be one finite number between 0 and 1.",
         call. = FALSE)
  }
  numeric_fields <- c("gc_content", "exon_length", "distance_to_junction",
                      "evolutionary_conservation")
  for (nm in numeric_fields) {
    x <- feature_df[[nm]]
    if (!is.numeric(x) || any(!is.finite(x)) || any(x < 0)) {
      stop(paste(nm, "must contain finite nonnegative numbers."), call. = FALSE)
    }
  }
  for (nm in c("gc_content", "evolutionary_conservation")) {
    if (any(feature_df[[nm]] > 1)) {
      stop(paste(nm, "must be between 0 and 1."), call. = FALSE)
    }
  }
  allowed <- list(RNA_type = c("mRNA", "lincRNA", "lncRNA", "pseudogene"),
                  RNA_region = c("CDS", "intron", "3'UTR", "5'UTR"))
  for (nm in names(allowed)) {
    x <- as.character(feature_df[[nm]])
    if (anyNA(x) || any(!x %in% allowed[[nm]])) {
      stop(paste(nm, "has missing or unsupported categories."), call. = FALSE)
    }
    feature_df[[nm]] <- factor(x, levels = allowed[[nm]])
  }
  encoded <- dna_encoding(as.character(feature_df$DNA_5mer))
  # Restrict the prediction table to the intended features. Preserve all input
  # columns in the returned table, including identifiers supplied by the user.
  newdata <- cbind(feature_df[setdiff(required, "DNA_5mer")], encoded)
  inferred_prob <- predict(ml_fit, newdata = newdata, type = "prob")
  if (!"Positive" %in% colnames(inferred_prob)) {
    stop("ml_fit must predict a class named Positive.", call. = FALSE)
  }
  feature_df$predicted_m6A_prob <- as.numeric(inferred_prob[, "Positive"])
  feature_df$predicted_m6A_status <- ifelse(
    feature_df$predicted_m6A_prob > positive_threshold, "Positive", "Negative")
  feature_df
}

#' Predict m6A probability and status for one site
#'
#' Constructs one feature row and calls [prediction_multiple()]. The model,
#' feature scales and strict threshold rule are the same in both functions.
#'
#' @param ml_fit A fitted classification `randomForest` with a Positive class.
#' @param gc_content One finite numeric GC fraction between 0 and 1.
#' @param RNA_type One of mRNA, lincRNA, lncRNA or pseudogene.
#' @param RNA_region One of CDS, intron, 3'UTR or 5'UTR.
#' @param exon_length One finite nonnegative number on the course log2 scale.
#' @param distance_to_junction One finite nonnegative number on the course
#'   log2 scale, not a raw base-pair distance.
#' @param evolutionary_conservation One finite numeric value between 0 and 1.
#' @param DNA_5mer One five-base uppercase A/T/C/G context. The supplied model
#'   was trained on central-A contexts, such as GGACA.
#' @param positive_threshold One finite number between 0 and 1; default 0.5.
#'   Scores strictly greater than this value are Positive; equality is Negative.
#'
#' @return A named character vector of length two, with names
#'   `predicted_m6A_prob` and `predicted_m6A_status`. Combining a number and a
#'   class label coerces the score to character. Use `as.numeric()` to recover
#'   the numeric score when needed.
#' @export
#' @examples
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds",
#'                              package = "m6APrediction", mustWork = TRUE))
#' prediction_single(ml_fit, gc_content = 0.6, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 5,
#'                   evolutionary_conservation = 0.8, DNA_5mer = "GGACA",
#'                   positive_threshold = 0.5)
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region,
                              exon_length, distance_to_junction,
                              evolutionary_conservation, DNA_5mer,
                              positive_threshold = 0.5) {
  features <- list(gc_content = gc_content, RNA_type = RNA_type,
                   RNA_region = RNA_region, exon_length = exon_length,
                   distance_to_junction = distance_to_junction,
                   evolutionary_conservation = evolutionary_conservation,
                   DNA_5mer = DNA_5mer)
  if (any(lengths(features) != 1L)) {
    stop("Each feature must have length one for single-site prediction.",
         call. = FALSE)
  }
  pred_outcome <- prediction_multiple(
    ml_fit, as.data.frame(features, stringsAsFactors = FALSE),
    positive_threshold = positive_threshold)
  c(predicted_m6A_prob = as.character(pred_outcome$predicted_m6A_prob),
    predicted_m6A_status = pred_outcome$predicted_m6A_status)
}

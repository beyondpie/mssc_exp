## Ref: https://en.wikipedia.org/wiki/Receiver_operating_characteristic

tpr <- function(tp, fp, tn, fn) {
  ## * true positve rate / recall / sensitivity
  return(tp / (tp + fn))
}

fnr <- function(tp, fp, tn, fn) {
  ## * false negative (type II error) rate / miss rate / 1- tpr
  return(fn / (tp + fn))
}

fpr <- function(tp, fp, tn, fn) {
  ## * false positive (type I error) rate, fall-out
  return(fp / (tn + fp))
}

tnr <- function(tp, fp, tn, fn) {
  ## * true negative rate / specificity / selectivity / 1 - fpr
  return(tn / (tn / (tn + fp)))
}

acc <- function(tp, fp, tn, fn) {
  ## * accuracy
  return((tp + tn) / (tp + fp + tn + fn))
}

ppv <- function(tp, fp, tn, fn) {
  ## * positive predictive value / precision
  return(tp / (tp + fp))
}

fdr <- function(tp, fp, tn, fn) {
  ## * false discovery rate  /  1 - ppv
  return(fp / (tp + fp))
}

f1 <- function(tp, fp, tn, fn) {
  ## * harmonic mean of precision and sensitivity
  return(2 * tp / (2 * tp + fp + fn))
}

mcc <- function(tp, fp, tn, fn) {
  ## * Matthew correlation coefficient
  return((tp * tn - fp * fn) /
         sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)))
}

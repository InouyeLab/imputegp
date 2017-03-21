# Filter out samples with abnormal measurement values
#
# Look up a table to check whether input measurements / imputed glycoprotein
# measurements are within their range of values present in DILGOM, the training
# data. Abnormal concentrations are set to NA, since the glycoprotein
# concentrations could not/cannot be reliably imputed.
#
# @param measurement vector of measurements
# @param name of the measurement
#
# @return the measurement vector where entries outside of the predefined range
# are set to 'NA'.
check_range <- function(measurement, name) {
  before_n <- sum(!is.na(measurement))
  if (name == "Sex") {
    measurement[!(measurement %in% c(1L, 2L))] <- NA
  } else {
    # measurement_ranges is an internal data structure
    measurement[measurement < measurement_ranges[name, "min_val"]] <- NA
    measurement[measurement > measurement_ranges[name, "max_val"]] <- NA
  }
  after_n <- sum(!is.na(measurement))
  filtered_n <- before_n - after_n
  if (filtered_n > 0) {
    if (name == "Sex") {
      warning(filtered_n, " samples with unrecognisable sex coding (Male == 1, Female == 2)",
              " removed")
    } else if (name %in% c("AGP", "A1AT", "HP", "TF")) {
      message(filtered_n, " imputed \"", name, "\" measurements outside of acceptable range (",
              measurement_ranges[name, "min_val"], "--", measurement_ranges[name, "max_val"],
              " ", measurement_ranges[name, "units"], ") set to 'NA'")
    } else {
      message(filtered_n, " \"", name, "\" measurements outside of acceptable range (",
              measurement_ranges[name, "min_val"], "--", measurement_ranges[name, "max_val"],
              " ", measurement_ranges[name, "units"], ") set to 'NA'")
    }
  }
  return(measurement)
}





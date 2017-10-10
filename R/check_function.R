# Filter out samples with abnormal measurement values
#
# Look up a table to check whether input measurements / imputed glycoprotein
# measurements are within their range of values present in DILGOM, the training
# data. Abnormal concentrations are set to NA, since the glycoprotein
# concentrations could not/cannot be reliably imputed.
#
# Will also set 0 NMR values to the minim value in DILGOM
#
# @param measurement vector of measurements
# @param name of the measurement
# @param range_check should range checking be performed?
# @param na.omit should samples with missing values be excluded or should
#   missing values be imputed as the model training dataset median?
#
# @return the measurement vector where entries outside of the predefined range
# are set to 'NA'.
check_range <- function(measurement, name, range_check, na.omit) {
  # Set 0 values to DILGOM lower detection limits
  if (!(name %in% c("Sex", "Age", "BMI"))) {
    measurement[measurement == 0] <- measurement_ranges[name, "min_val"]
  }

  if (!na.omit) {
    measurement[is.na(measurement)] <- measurement_ranges[name, "median_val"]
  }

  before_n <- sum(!is.na(measurement))
  if (name == "Sex") {
    measurement[!(measurement %in% c(1L, 2L))] <- NA
  } else {
    if (range_check) {
      # measurement_ranges is an internal data structure
      measurement[measurement < measurement_ranges[name, "min_val"]] <- NA
      measurement[measurement > measurement_ranges[name, "max_val"]] <- NA
    }
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

# Handle missing measurements
#
# @param values vector of measurement values
# @param name name of the measurement
# @param standardised logical; has the measurement been standardised?
# @param na.omit logical; should NAs be kept as NA or imputed to the
#   model training dataset median?
#
# @return a vector of values
handle_nas <- function(values, name, standardised, na.omit) {
  if (na.omit) {
    return(values)
  } else {
    if (standardised) {
      values[is.na(values)] <- 0
    } else {
      values[is.na(values)] <- measurement_ranges[name, "median_val"]
    }
  }
}






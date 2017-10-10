#' Impute Alpha-1 antitrypsin (AAT)
#'
#' Imputes the concentrations of AAT from serum NMR measurements in healthy
#' population-based samples.
#'
#' @details
#'  The imputation models will only return a concentration where no input
#'  measurements were missing, and where all input measurements were within
#'  their acceptable predefined range of values. These correspond to their
#'  range of values in the training dataset. Similarly, the imputed measurement
#'  will be set to missing if the imputation model returns a concentration
#'  outside the range of values that were found in the training dataset.
#'
#'  Standardised measurements may also be used by setting \code{standardised = 'TRUE'},
#'  which is useful for cases where measurements must be adjusted for
#'  technical effects. In this case, all NMR measurements and BMI should
#'  be log transformed before scaling. Age and sex must also be standardised.
#'
#' @param GlycA NMR measurement for glycoprotein acetyls. Measurements should
#'  be between 0.869--2.24 mmol/L.
#' @param FAw3 NMR measurement for Omega-3 fatty acids. Measurements should be
#'  between 0.196--1.18 mmol/L.
#' @param VLDL.D NMR measurement for the mean diameter of VLDL particles.
#'  Measurements should be between 33--41 nm.
#' @param HDL3.C NMR measurement for cholesterol in the HDL-3 fraction.
#'  Measurements should be between 0.421--0.647 mmol/L.
#' @param LDL.D NMR measurement for the mean diameter of LDL particles.
#'  Measurements should be between 23.2--24.1 nm.
#' @param Phe NMR measurement for phenylalanine. Measurements should be between
#'  0.0555--0.118 mmol/L.
#' @param Leu NMR measurement for leucine. Measurements should be between
#'  0.0295--0.138 mmol/L.
#' @param ApoB NMR measurement for Apolipoprotein B. Measurements should be
#'  between 0.365--1.48 g/L.
#' @param Alb NMR measurement for albumin. Measurements should have a signal
#'  area ranging between 0.0733--0.101.
#' @param Tyr NMR measurement for tyrosine. Measurements should be between
#'  0.0279--0.124 mmol/L.
#' @param bOHBut NMR measurement for 3-hydroxybutyrate. Measurements should be
#'  between 0.0331--0.872 mmol/L.
#' @param BMI Body Mass Index. Measurements should be between 16.17--47.44
#'  kg/m^2.
#' @param Ala NMR measurement for alanine. Measurements should be between
#'  0.274--0.557 mmol/L.
#' @param L.HDL.TG NMR measurement for triglycerides within large HDL
#'  particles. Measurements should be between 0.000462--0.0933 mmol/L.
#' @param Ile NMR measurement for isoleucine. Measurements should be between
#'  0.024--0.103 mmol/L.
#' @param Ace NMR measurement for acetate. Measurements should be between
#'  0.029--1.32 mmol/L.
#' @param His NMR measurement for histidine. Measurements should be between
#'  0.0446--0.0913 mmol/L.
#' @param HDL.TG NMR measurement for triglycerides within all HDL particles.
#'  Measurements should be between 0.0698--0.315 mmol/L.
#' @param range_check logical; if \code{TRUE} discard measurements that
#'  are not in the accepted range of values (see Details). If \code{FALSE},
#'  no checking of input measurements or predicted concentrations will
#'  be performed.
#' @param standardised logical; have measurements been standardised
#'  (\emph{i.e.} using the \code{scale} function.)
#' @param na.omit logical; should samples with missing values be omited?
#'  If \code{FALSE} missing values are set to the measurement's median
#'  in the model training dataset. Alternatively consider imputing
#'  missing values using \code{\link[impute]{impute.knn}}.
#'
#' @return A vector of AAT measurements ranging between 0.64--2.58 mg/L or
#'  measurements standardised to the population if \code{standardised = 'TRUE'}.
#'
#' @export
impute_AAT <- function(
  GlycA, FAw3, VLDL.D, HDL3.C, LDL.D, Phe, Leu, ApoB, Alb, Tyr, bOHBut, BMI,
  Ala, L.HDL.TG, Ile, Ace, His, HDL.TG, range_check=TRUE, standardised=FALSE,
  na.omit=TRUE
) {
  if (standardised) {
    AAT <- AAT_coef["intercept", "standardised"] +
      AAT_coef["GlycA", "standardised"] * handle_nas(GlycA, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["FAw3", "standardised"] * handle_nas(FAw3, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["VLDL.D", "standardised"] * handle_nas(VLDL.D, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["HDL3.C", "standardised"] * handle_nas(HDL3.C, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["LDL.D", "standardised"] * handle_nas(LDL.D, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["Phe", "standardised"] * handle_nas(Phe, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["Leu", "standardised"] * handle_nas(Leu, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["ApoB", "standardised"] * handle_nas(ApoB, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["Alb", "standardised"] * handle_nas(Alb, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["Tyr", "standardised"] * handle_nas(Tyr, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["bOHBut", "standardised"] * handle_nas(bOHBut, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["BMI", "standardised"] * handle_nas(BMI, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["Ala", "standardised"] * handle_nas(Ala, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["L.HDL.TG", "standardised"] * handle_nas(L.HDL.TG, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["Ile", "standardised"] * handle_nas(Ile, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["Ace", "standardised"] * handle_nas(Ace, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["His", "standardised"] * handle_nas(His, standardised=TRUE, na.omit=na.omit) +
      AAT_coef["HDL.TG", "standardised"] * handle_nas(HDL.TG, standardised=TRUE, na.omit=na.omit)
    AAT <- scale(AAT)
  } else {
    log_AAT <- AAT_coef["intercept", "raw"] +
      AAT_coef["GlycA", "raw"] * log(check_range(GlycA, "GlycA", range_check, na.omit)) +
      AAT_coef["FAw3", "raw"] * log(check_range(FAw3, "FAw3", range_check, na.omit)) +
      AAT_coef["VLDL.D", "raw"] * log(check_range(VLDL.D, "VLDL.D", range_check, na.omit)) +
      AAT_coef["HDL3.C", "raw"] * log(check_range(HDL3.C, "HDL3.C", range_check, na.omit)) +
      AAT_coef["LDL.D", "raw"] * log(check_range(LDL.D, "LDL.D", range_check, na.omit)) +
      AAT_coef["Phe", "raw"] * log(check_range(Phe, "Phe", range_check, na.omit)) +
      AAT_coef["Leu", "raw"] * log(check_range(Leu, "Leu", range_check, na.omit)) +
      AAT_coef["ApoB", "raw"] * log(check_range(ApoB, "ApoB", range_check, na.omit)) +
      AAT_coef["Alb", "raw"] * log(check_range(Alb, "Alb", range_check, na.omit)) +
      AAT_coef["Tyr", "raw"] * log(check_range(Tyr, "Tyr", range_check, na.omit)) +
      AAT_coef["bOHBut", "raw"] * log(check_range(bOHBut, "bOHBut", range_check, na.omit)) +
      AAT_coef["BMI", "raw"] * log(check_range(BMI, "BMI", range_check, na.omit)) +
      AAT_coef["Ala", "raw"] * log(check_range(Ala, "Ala", range_check, na.omit)) +
      AAT_coef["L.HDL.TG", "raw"] * log(check_range(L.HDL.TG, "L.HDL.TG", range_check, na.omit)) +
      AAT_coef["Ile", "raw"] * log(check_range(Ile, "Ile", range_check, na.omit)) +
      AAT_coef["Ace", "raw"] * log(check_range(Ace, "Ace", range_check, na.omit)) +
      AAT_coef["His", "raw"] * log(check_range(His, "His", range_check, na.omit)) +
      AAT_coef["HDL.TG", "raw"] * log(check_range(HDL.TG, "HDL.TG", range_check, na.omit))

    # Transform back to raw concentration units
    AAT <- exp(log_AAT)

    # Remove imputed concentrations that are outside the range of concentration
    # values observed in the model training data.
    AAT <- check_range(AAT, "AAT", range_check, na.omit)
  }

  message("Successfully imputed AAT for ", sum(!is.na(AAT)), " samples")

  return(AAT)
}

#' Impute Alpha-1-acid glcyoprotein (AGP)
#'
#' Imputes the concentrations of AGP from serum NMR measurements in healthy
#' population-based samples.
#'
#' @details
#'  The imputation models will only return a concentration where no input
#'  measurements were missing, and where all input measurements were within
#'  their acceptable predefined range of values. These correspond to their
#'  range of values in the training dataset. Similarly, the imputed measurement
#'  will be set to missing if the imputation model returns a concentration
#'  outside the range of values that were found in the training dataset.
#'
#'  Standardised measurements may also be used by setting \code{standardised = 'TRUE'},
#'  which is useful for cases where measurements must be adjusted for
#'  technical effects. In this case, all NMR measurements and BMI should
#'  be log transformed before scaling. Age and sex must also be standardised.
#'
#' @param GlycA NMR measurement for glycoprotein acetyls. Measurements should
#'  be between 0.869--2.24 mmol/L.
#' @param TotFA NMR measurement for total fatty acids. Measurements should be
#'  between 5.58--17.1 mmol/L
#' @param IDL.FC NMR measurement for free cholesterol within IDL particles.
#'  Measurements should be between 0.0515--0.368 mmol/L.
#' @param L.HDL.FC NMR measurement for free cholesterol within large HDL
#'  particles. Measurements should be between 0.00134--0.208 mmol/L.
#' @param His NMR measurement for histidine. Measurements should be between
#'  0.0446--0.0913 mmol/L.
#' @param HDL.TG NMR measurement for triglycerides within all HDL particles.
#'  Measurements should be between 0.0698--0.315 mmol/L.
#' @param BMI Body Mass Index. Measurements should be between 16.17--47.44
#'  kg/m^2.
#' @param S.HDL.FC NMR measurement for free cholesterol within small HDL
#'  particles. Measurements should be between 0.0674--0.168 mmol/L.
#' @param S.LDL.TG NMR measurement for triglycerides within small LDL
#'  particles. Measurements should be between 0.00478--0.076 mmol/L.
#' @param bOHBut NMR measurement for 3-hydroxybutyrate. Measurements should be
#'  between 0.0331--0.872 mmol/L.
#' @param LA NMR measurement for 18:2 linoleic acid. Measurements should be
#'  between 0.449--5.15 mmol/L.
#' @param S.HDL.CE NMR measurement for cholesterol esters in small HDL
#'  particles. Measurements should be between 0.061--0.466 mmol/L.
#' @param Lac NMR measurement for lactate. Measurements should be between
#'  0.69--4.93 mmol/L.
#' @param S.VLDL.TG NMR measurement for triglycerides within small VLDL
#'  particles. Measurements should be between 0.0135--0.76 mmol/L.
#' @param Ace NMR measurement for acetate. Measurements should be between
#'  0.029--1.32 mmol/L.
#' @param Cit NMR measurement for citrate. Measurements should be between
#'  0.0739--0.186 mmol/L.
#' @param SFA NMR measurement for saturated fatty acids. Measurements should
#'  be between 2.07--7.02 mmol/L.
#' @param Ala NMR measurement for alanine. Measurements should be between
#'  0.274--0.557 mmol/L.
#' @param XXL.VLDL.CE NMR measurement for cholesterol esters within extremely
#'  large VLDL particles. Measurements should be between 0.000000133--0.0169
#'  mmol/L.
#' @param Glol NMR measurement for glycerol. Measurements should be between
#'  0.0415--0.312 mmol/L.
#' @param Age The study participant's age in years. Measurements should be
#'  between 25--74 years old.
#' @param Crea NMR measurement for creatinine. Measurements should be between
#'  0.0182--0.347 mmol/L.
#' @param Gly NMR measurement for glycine. Measurements should be between
#'  0.182--0.632 mmol/L.
#' @param range_check logical; if \code{TRUE} discard measurements that
#'  are not in the accepted range of values (see Details). If \code{FALSE},
#'  no checking of input measurements or predicted concentrations will
#'  be performed.
#' @param standardised logical; have measurements been standardised
#'  (\emph{i.e.} using the \code{scale} function.)
#' @param na.omit logical; should samples with missing values be omited?
#'  If \code{FALSE} missing values are set to the measurement's median
#'  in the model training dataset. Alternatively consider imputing
#'  missing values using \code{\link[impute]{impute.knn}}.
#'
#' @return A vector of AGP measurements ranging between 362--1,880 mg/L or
#'  measurements standardised to the population if \code{standardised = 'TRUE'}.
#'
#' @export
impute_AGP <- function(
  GlycA, TotFA, IDL.FC, L.HDL.FC, His, HDL.TG, BMI, S.HDL.FC, S.LDL.TG, bOHBut,
  LA, S.HDL.CE, Lac, S.VLDL.TG, Ace, Cit, SFA, Ala, XXL.VLDL.CE, Glol, Age,
  Crea, Gly, range_check=TRUE, standardised=FALSE, na.omit=TRUE
) {
  if (standardised) {
    AGP <- AGP_coef["intercept", "standardised"] +
      AGP_coef["GlycA", "standardised"] * handle_nas(GlycA, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["TotFA", "standardised"] * handle_nas(TotFA, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["IDL.FC", "standardised"] * handle_nas(IDL.FC, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["L.HDL.FC", "standardised"] * handle_nas(L.HDL.FC, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["His", "standardised"] * handle_nas(His, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["HDL.TG", "standardised"] * handle_nas(HDL.TG, standardised=TRUE, na.omit=na.omit)+
      AGP_coef["BMI", "standardised"] * handle_nas(BMI, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["S.HDL.FC", "standardised"] * handle_nas(S.HDL.FC, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["S.LDL.TG", "standardised"] * handle_nas(S.LDL.TG, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["bOHBut", "standardised"] * handle_nas(bOHBut, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["LA", "standardised"] * handle_nas(LA, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["S.HDL.CE", "standardised"] * handle_nas(S.HDL.CE, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["Lac", "standardised"] * handle_nas(Lac, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["S.VLDL.TG", "standardised"] * handle_nas(S.VLDL.TG, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["Ace", "standardised"] * handle_nas(Ace, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["Cit", "standardised"] * handle_nas(Cit, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["SFA", "standardised"] * handle_nas(SFA, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["Ala", "standardised"] * handle_nas(Ala, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["XXL.VLDL.CE", "standardised"] * handle_nas(XXL.VLDL.CE, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["Glol", "standardised"] * handle_nas(Glol, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["Age", "standardised"] * handle_nas(Age, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["Crea", "standardised"] * handle_nas(Crea, standardised=TRUE, na.omit=na.omit) +
      AGP_coef["Gly", "standardised"] * handle_nas(Gly, standardised=TRUE, na.omit=na.omit)
    AGP <- scale(AGP)
  } else {
    log_AGP <- AGP_coef["intercept", "raw"] +
      AGP_coef["GlycA", "raw"] * log(check_range(GlycA, "GlycA", range_check, na.omit)) +
      AGP_coef["TotFA", "raw"] * log(check_range(TotFA, "TotFA", range_check, na.omit)) +
      AGP_coef["IDL.FC", "raw"] * log(check_range(IDL.FC, "IDL.FC", range_check, na.omit)) +
      AGP_coef["L.HDL.FC", "raw"] * log(check_range(L.HDL.FC, "L.HDL.FC", range_check, na.omit)) +
      AGP_coef["His", "raw"] * log(check_range(His, "His", range_check, na.omit)) +
      AGP_coef["HDL.TG", "raw"] * log(check_range(HDL.TG, "HDL.TG", range_check, na.omit)) +
      AGP_coef["BMI", "raw"] * log(check_range(BMI, "BMI", range_check, na.omit)) +
      AGP_coef["S.HDL.FC", "raw"] * log(check_range(S.HDL.FC, "S.HDL.FC", range_check, na.omit)) +
      AGP_coef["S.LDL.TG", "raw"] * log(check_range(S.LDL.TG, "S.LDL.TG", range_check, na.omit)) +
      AGP_coef["bOHBut", "raw"] * log(check_range(bOHBut, "bOHBut", range_check, na.omit)) +
      AGP_coef["LA", "raw"] * log(check_range(LA, "LA", range_check, na.omit)) +
      AGP_coef["S.HDL.CE", "raw"] * log(check_range(S.HDL.CE, "S.HDL.CE", range_check, na.omit)) +
      AGP_coef["Lac", "raw"] * log(check_range(Lac, "Lac", range_check, na.omit)) +
      AGP_coef["S.VLDL.TG", "raw"] * log(check_range(S.VLDL.TG, "S.VLDL.TG", range_check, na.omit)) +
      AGP_coef["Ace", "raw"] * log(check_range(Ace, "Ace", range_check, na.omit)) +
      AGP_coef["Cit", "raw"] * log(check_range(Cit, "Cit", range_check, na.omit)) +
      AGP_coef["SFA", "raw"] * log(check_range(SFA, "SFA", range_check, na.omit)) +
      AGP_coef["Ala", "raw"] * log(check_range(Ala, "Ala", range_check, na.omit)) +
      AGP_coef["XXL.VLDL.CE", "raw"] * log(check_range(XXL.VLDL.CE, "XXL.VLDL.CE", range_check, na.omit)) +
      AGP_coef["Glol", "raw"] * log(check_range(Glol, "Glol", range_check, na.omit)) +
      AGP_coef["Age", "raw"] * check_range(Age, "Age", range_check, na.omit) +
      AGP_coef["Crea", "raw"] * log(check_range(Crea, "Crea", range_check, na.omit)) +
      AGP_coef["Gly", "raw"] * log(check_range(Gly, "Gly", range_check, na.omit))

    # Transform back to raw concentration units
    AGP <- exp(log_AGP)

    # Remove imputed concentrations that are outside the range of concentration
    # values observed in the model training data.
    AGP <- check_range(AGP, "AGP", range_check, na.omit)
  }

  message("Successfully imputed AGP for ", sum(!is.na(AGP)), " samples")

  return(AGP)
}

#' Impute Haptoglobin (HP)
#'
#' Imputes the concentrations ofHP from serum NMR measurements in healthy
#' population-based samples.
#'
#' @details
#'  The imputation models will only return a concentration where no input
#'  measurements were missing, and where all input measurements were within
#'  their acceptable predefined range of values. These correspond to their
#'  range of values in the training dataset. Similarly, the imputed measurement
#'  will be set to missing if the imputation model returns a concentration
#'  outside the range of values that were found in the training dataset.
#'
#'  Standardised measurements may also be used by setting \code{standardised = 'TRUE'},
#'  which is useful for cases where measurements must be adjusted for
#'  technical effects. In this case, all NMR measurements and BMI should
#'  be log transformed before scaling. Age and sex must also be standardised.
#'
#' @param GlycA NMR measurement for glycoprotein acetyls. Measurements should
#'  be between 0.869--2.24 mmol/L.
#' @param LA NMR measurement for 18:2 linoleic acid. Measurements should be
#'  between 0.449--5.15 mmol/L.
#' @param IDL.FC NMR measurement for free cholesterol within IDL particles.
#'  Measurements should be between 0.0515--0.368 mmol/L.
#' @param SM NMR measurement for sphingomyelins. Measurements should be between
#'  0.103--0.733 mmol/L.
#' @param FAw3 NMR measurement for Omega-3 fatty acids. Measurements should be
#'  between 0.196--1.18 mmol/L.
#' @param HDL.TG NMR measurement for triglycerides within all HDL particles.
#'  Measurements should be between 0.0698--0.315 mmol/L.
#' @param S.VLDL.CE NMR measurement for cholesterol esters within small VLDL
#'  particles. Measurements should be between 0.000221--0.288 mmol/L.
#' @param Age The study participant's age in years. Measurements should be
#'  between 25--74 years old.
#' @param Alb NMR measurement for albumin. Measurements should have a signal
#'  area ranging between 0.0733--0.101.
#' @param Ile NMR measurement for isoleucine. Measurements should be between
#'  0.024--0.103 mmol/L.
#' @param Cit NMR measurement for citrate. Measurements should be between
#'  0.0739--0.186 mmol/L.
#' @param VLDL.D NMR measurement for the mean diameter of VLDL particles.
#'  Measurements should be between 33--41 nm.
#' @param Leu NMR measurement for leucine. Measurements should be between
#'  0.0295--0.138 mmol/L.
#' @param Val NMR measurement for valine. Measurements should be between
#'  0.0728--0.297 mmol/L.
#' @param L.VLDL.CE NMR measurement for cholesterol esters within large VLDL
#'  particles. Measurements should be between 0.00000152--0.161 mmol/L.
#' @param Pyr NMR measurement for pyruvate. Measurements should be between
#'  0.0512--0.179 mmol/L.
#' @param Lac NMR measurement for lactate. Measurements should be between
#'  0.69--4.93 mmol/L.
#' @param Gln NMR measurement for glutamine. Measurements should be between
#'  0.331--1.36 mmol/L.
#' @param M.HDL.FC NMR measurement for free cholesterol within medium HDL
#'  particles. Measurements should be between 0.0208--0.145 mmol/L.
#' @param XL.HDL.TG NMR measurement for triglycerides within very large HDL
#'  particles. Measurements should be between 0.000287--0.055 mmol/L.
#' @param XL.HDL.PL NMR measurement for phospholipids within very large HDL
#'  particles. Measurements should be between 0.00176--0.639 mmol/L.
#' @param His NMR measurement for histidine. Measurements should be between
#'  0.0446--0.0913 mmol/L.
#' @param Tyr NMR measurement for tyrosine. Measurements should be between
#'  0.0279--0.124 mmol/L.
#' @param BMI Body Mass Index. Measurements should be between 16.17--47.44
#'  kg/m^2.
#' @param L.HDL.TG NMR measurement for triglycerides within large HDL
#'  particles. Measurements should be between 0.000462--0.0993 mmol/L.
#' @param PUFA NMR measurement for polyunsaturated fatty acids. Measurements
#'  should be between 2.29--6.62 mmol/L.
#' @param S.LDL.FC NMR measurement for free cholesterol within small LDL
#'  particles. Measurements should be between 0.0215--0.133 mmol/L.
#' @param range_check logical; if \code{TRUE} discard measurements that
#'  are not in the accepted range of values (see Details). If \code{FALSE},
#'  no checking of input measurements or predicted concentrations will
#'  be performed.
#' @param standardised logical; have measurements been standardised
#'  (\emph{i.e.} using the \code{scale} function.)
#' @param na.omit logical; should samples with missing values be omited?
#'  If \code{FALSE} missing values are set to the measurement's median
#'  in the model training dataset. Alternatively consider imputing
#'  missing values using \code{\link[impute]{impute.knn}}.
#'
#' @return A vector of HP measurements ranging between 0.14--3.95 mg/L or
#'  measurements standardised to the population if \code{standardised = 'TRUE'}.
#'
#' @export
impute_HP <- function(
  GlycA, LA, IDL.FC, SM, FAw3, HDL.TG, S.VLDL.CE, Age, Alb, Ile, Cit, VLDL.D,
  Leu, Val, L.VLDL.CE, Pyr, Lac, Gln, M.HDL.FC, XL.HDL.TG, XL.HDL.PL, His, Tyr,
  BMI, L.HDL.TG, PUFA, S.LDL.FC, range_check=TRUE, standardised=FALSE, na.omit=TRUE
) {
  if (standardised) {
    HP <- HP_coef["intercept", "standardised"] +
      HP_coef["GlycA", "standardised"] * handle_nas(GlycA, standardised=TRUE, na.omit=na.omit) +
      HP_coef["LA", "standardised"] * handle_nas(LA, standardised=TRUE, na.omit=na.omit) +
      HP_coef["IDL.FC", "standardised"] * handle_nas(IDL.FC, standardised=TRUE, na.omit=na.omit) +
      HP_coef["SM", "standardised"] * handle_nas(SM, standardised=TRUE, na.omit=na.omit) +
      HP_coef["FAw3", "standardised"] * handle_nas(FAw3, standardised=TRUE, na.omit=na.omit) +
      HP_coef["HDL.TG", "standardised"] * handle_nas(HDL.TG, standardised=TRUE, na.omit=na.omit) +
      HP_coef["S.VLDL.CE", "standardised"] * handle_nas(S.VLDL.CE, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Age", "standardised"] * handle_nas(Age, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Alb", "standardised"] * handle_nas(Alb, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Ile", "standardised"] * handle_nas(Ile, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Cit", "standardised"] * handle_nas(Cit, standardised=TRUE, na.omit=na.omit) +
      HP_coef["VLDL.D", "standardised"] * handle_nas(VLDL.D, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Leu", "standardised"] * handle_nas(Leu, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Val", "standardised"] * handle_nas(Val, standardised=TRUE, na.omit=na.omit) +
      HP_coef["L.VLDL.CE", "standardised"] * handle_nas(L.VLDL.CE, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Pyr", "standardised"] * handle_nas(Pyr, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Lac", "standardised"] * handle_nas(Lac, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Gln", "standardised"] * handle_nas(Gln, standardised=TRUE, na.omit=na.omit) +
      HP_coef["M.HDL.FC", "standardised"] * handle_nas(M.HDL.FC, standardised=TRUE, na.omit=na.omit) +
      HP_coef["XL.HDL.TG", "standardised"] * handle_nas(XL.HDL.TG, standardised=TRUE, na.omit=na.omit) +
      HP_coef["XL.HDL.PL", "standardised"] * handle_nas(XL.HDL.PL, standardised=TRUE, na.omit=na.omit) +
      HP_coef["His", "standardised"] * handle_nas(His, standardised=TRUE, na.omit=na.omit) +
      HP_coef["Tyr", "standardised"] * handle_nas(Tyr, standardised=TRUE, na.omit=na.omit) +
      HP_coef["BMI", "standardised"] * handle_nas(BMI, standardised=TRUE, na.omit=na.omit) +
      HP_coef["L.HDL.TG", "standardised"] * handle_nas(L.HDL.TG, standardised=TRUE, na.omit=na.omit) +
      HP_coef["PUFA", "standardised"] * handle_nas(PUFA, standardised=TRUE, na.omit=na.omit) +
      HP_coef["S.LDL.FC", "standardised"] * handle_nas(S.LDL.FC, standardised=TRUE, na.omit=na.omit)
    HP <- scale(HP)
  } else {
    log_HP <- HP_coef["intercept", "raw"] +
      HP_coef["GlycA", "raw"] * log(check_range(GlycA, "GlycA", range_check, na.omit)) +
      HP_coef["LA", "raw"] * log(check_range(LA, "LA", range_check, na.omit)) +
      HP_coef["IDL.FC", "raw"] * log(check_range(IDL.FC, "IDL.FC", range_check, na.omit)) +
      HP_coef["SM", "raw"] * log(check_range(SM, "SM", range_check, na.omit)) +
      HP_coef["FAw3", "raw"] * log(check_range(FAw3, "FAw3", range_check, na.omit)) +
      HP_coef["HDL.TG", "raw"] * log(check_range(HDL.TG, "HDL.TG", range_check, na.omit)) +
      HP_coef["S.VLDL.CE", "raw"] * log(check_range(S.VLDL.CE, "S.VLDL.CE", range_check, na.omit)) +
      HP_coef["Age", "raw"] * check_range(Age, "Age", range_check, na.omit) +
      HP_coef["Alb", "raw"] * log(check_range(Alb, "Alb", range_check, na.omit)) +
      HP_coef["Ile", "raw"] * log(check_range(Ile, "Ile", range_check, na.omit)) +
      HP_coef["Cit", "raw"] * log(check_range(Cit, "Cit", range_check, na.omit)) +
      HP_coef["VLDL.D", "raw"] * log(check_range(VLDL.D, "VLDL.D", range_check, na.omit)) +
      HP_coef["Leu", "raw"] * log(check_range(Leu, "Leu", range_check, na.omit)) +
      HP_coef["Val", "raw"] * log(check_range(Val, "Val", range_check, na.omit)) +
      HP_coef["L.VLDL.CE", "raw"] * log(check_range(L.VLDL.CE, "L.VLDL.CE", range_check, na.omit)) +
      HP_coef["Pyr", "raw"] * log(check_range(Pyr, "Pyr", range_check, na.omit)) +
      HP_coef["Lac", "raw"] * log(check_range(Lac, "Lac", range_check, na.omit)) +
      HP_coef["Gln", "raw"] * log(check_range(Gln, "Gln", range_check, na.omit)) +
      HP_coef["M.HDL.FC", "raw"] * log(check_range(M.HDL.FC, "M.HDL.FC", range_check, na.omit)) +
      HP_coef["XL.HDL.TG", "raw"] * log(check_range(XL.HDL.TG, "XL.HDL.TG", range_check, na.omit)) +
      HP_coef["XL.HDL.PL", "raw"] * log(check_range(XL.HDL.PL, "XL.HDL.PL", range_check, na.omit)) +
      HP_coef["His", "raw"] * log(check_range(His, "His", range_check, na.omit)) +
      HP_coef["Tyr", "raw"] * log(check_range(Tyr, "Tyr", range_check, na.omit)) +
      HP_coef["BMI", "raw"] * log(check_range(BMI, "BMI", range_check, na.omit)) +
      HP_coef["L.HDL.TG", "raw"] * log(check_range(L.HDL.TG, "L.HDL.TG", range_check, na.omit)) +
      HP_coef["PUFA", "raw"] * log(check_range(PUFA, "PUFA", range_check, na.omit)) +
      HP_coef["S.LDL.FC", "raw"] * log(check_range(S.LDL.FC, "S.LDL.FC", range_check, na.omit))

    # Transform back to raw concentration units
    HP <- exp(log_HP)

    # Remove imputed concentrations that are outside the range of concentration
    # values observed in the model training data.
    HP <- check_range(HP, "HP", range_check, na.omit)
  }

  message("Successfully imputed HP for ", sum(!is.na(HP)), " samples")

  return(HP)
}

#' Impute Transferrin (TF)
#'
#' Imputes the concentrations of AGP from serum NMR measurements in healthy
#' population-based samples.
#'
#' @details
#'  The imputation models will only return a concentration where no input
#'  measurements were missing, and where all input measurements were within
#'  their acceptable predefined range of values. These correspond to their
#'  range of values in the training dataset. Similarly, the imputed measurement
#'  will be set to missing if the imputation model returns a concentration
#'  outside the range of values that were found in the training dataset.
#'
#'  Standardised measurements may also be used by setting \code{standardised = 'TRUE'},
#'  which is useful for cases where measurements must be adjusted for
#'  technical effects. In this case, all NMR measurements and BMI should
#'  be log transformed before scaling. Age and sex must also be standardised.
#'
#' @param GlycA NMR measurement for glycoprotein acetyls. Measurements should
#'  be between 0.869--2.24 mmol/L.
#' @param Sex The study participant's sex. The coding should be 1 for males, 2 for
#'  females.
#' @param Age The study participant's age in years. Measurements should be
#'  between 25--74 years old.
#' @param S.HDL.FC NMR measurement for free cholesterol within small HDL
#'  particles. Measurements should be between 0.0674--0.168 mmol/L.
#' @param Ace NMR measurement for acetate. Measurements should be between
#'  0.029--1.32 mmol/L.
#' @param Ala NMR measurement for alanine. Measurements should be between
#'  0.274--0.557 mmol/L.
#' @param SFA NMR measurement for saturated fatty acids. Measurements should
#'  be between 2.07--7.02 mmol/L.
#' @param His NMR measurement for histidine. Measurements should be between
#'  0.0446--0.0913 mmol/L.
#' @param Gln NMR measurement for glutamine. Measurements should be between
#'  0.331--1.36 mmol/L.
#' @param range_check logical; if \code{TRUE} discard measurements that
#'  are not in the accepted range of values (see Details). If \code{FALSE},
#'  no checking of input measurements or predicted concentrations will
#'  be performed.
#' @param standardised logical; have measurements been standardised
#'  (\emph{i.e.} using the \code{scale} function.)
#' @param na.omit logical; should samples with missing values be omited?
#'  If \code{FALSE} missing values are set to the measurement's median
#'  in the model training dataset. Alternatively consider imputing
#'  missing values using \code{\link[impute]{impute.knn}}.
#'
#' @return A vector of TF measurements ranging between 1.39--4.38 mg/L or
#'  measurements standardised to the population if \code{standardised = 'TRUE'}.
#'
impute_TF <- function(GlycA, Sex, Age, S.HDL.FC, Ace, Ala, SFA, His, Gln,
                      range_check=TRUE, standardised=FALSE, na.omit=TRUE) {

  if (standardised) {
    TF <- TF_coef["intercept", "standardised"] +
      TF_coef["GlycA", "standardised"] * handle_nas(GlycA, standardised=TRUE, na.omit=na.omit) +
      TF_coef["Sex", "standardised"] * handle_nas(Sex, standardised=TRUE, na.omit=na.omit) +
      TF_coef["Age", "standardised"] * handle_nas(Age, standardised=TRUE, na.omit=na.omit) +
      TF_coef["S.HDL.FC", "standardised"] * handle_nas(S.HDL.FC, standardised=TRUE, na.omit=na.omit) +
      TF_coef["Ace", "standardised"] * handle_nas(Ace, standardised=TRUE, na.omit=na.omit) +
      TF_coef["Ala", "standardised"] * handle_nas(Ala, standardised=TRUE, na.omit=na.omit) +
      TF_coef["SFA", "standardised"] * handle_nas(SFA, standardised=TRUE, na.omit=na.omit) +
      TF_coef["His", "standardised"] * handle_nas(His, standardised=TRUE, na.omit=na.omit) +
      TF_coef["Gln", "standardised"] * handle_nas(Gln, standardised=TRUE, na.omit=na.omit)
    TF <- scale(TF)
  } else {
    log_TF <- TF_coef["intercept", "raw"] +
      TF_coef["GlycA", "raw"] * log(check_range(GlycA, "GlycA", range_check, na.omit)) +
      TF_coef["Sex", "raw"] * check_range(Sex, "Sex", range_check, na.omit) +
      TF_coef["Age", "raw"] * check_range(Age, "Age", range_check, na.omit) +
      TF_coef["S.HDL.FC", "raw"] * log(check_range(S.HDL.FC, "S.HDL.FC", range_check, na.omit)) +
      TF_coef["Ace", "raw"] * log(check_range(Ace, "Ace", range_check, na.omit)) +
      TF_coef["Ala", "raw"] * log(check_range(Ala, "Ala", range_check, na.omit)) +
      TF_coef["SFA", "raw"] * log(check_range(SFA, "SFA", range_check, na.omit)) +
      TF_coef["His", "raw"] * log(check_range(His, "His", range_check, na.omit)) +
      TF_coef["Gln", "raw"] * log(check_range(Gln, "Gln", range_check, na.omit))

    # Transform back to raw concentration units
    TF <- exp(log_TF)

    # Remove imputed concentrations that are outside the range of concentration
    # values observed in the model training data.
    TF <- check_range(TF, "TF", range_check, na.omit)
  }

  message("Successfully imputed TF for ", sum(!is.na(TF)), " samples")
  warning("The concentrations predicted by this function are unlikely to reflect true TF levels", immediate.=TRUE)
  return(TF)
}

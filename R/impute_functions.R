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
#'
#' @return A vector of AAT measurements ranging between 0.64--2.58 mg/L or
#'  measurements standardised to the population if \code{standardised = 'TRUE'}.
#'
#' @export
impute_AAT <- function(
  GlycA, FAw3, VLDL.D, HDL3.C, LDL.D, Phe, Leu, ApoB, Alb, Tyr, bOHBut, BMI,
  Ala, L.HDL.TG, Ile, Ace, His, HDL.TG, range_check=TRUE, standardised=FALSE
) {
  if (standardised) {
    AAT <- AAT_coef["intercept", "standardised"] +
      AAT_coef["GlycA", "standardised"] * GlycA +
      AAT_coef["FAw3", "standardised"] * FAw3 +
      AAT_coef["VLDL.D", "standardised"] * VLDL.D +
      AAT_coef["HDL3.C", "standardised"] * HDL3.C +
      AAT_coef["LDL.D", "standardised"] * LDL.D +
      AAT_coef["Phe", "standardised"] * Phe +
      AAT_coef["Leu", "standardised"] * Leu +
      AAT_coef["ApoB", "standardised"] * ApoB +
      AAT_coef["Alb", "standardised"] * Alb +
      AAT_coef["Tyr", "standardised"] * Tyr +
      AAT_coef["bOHBut", "standardised"] * bOHBut +
      AAT_coef["BMI", "standardised"] * BMI +
      AAT_coef["Ala", "standardised"] * Ala +
      AAT_coef["L.HDL.TG", "standardised"] * L.HDL.TG +
      AAT_coef["Ile", "standardised"] * Ile +
      AAT_coef["Ace", "standardised"] * Ace +
      AAT_coef["His", "standardised"] * His +
      AAT_coef["HDL.TG", "standardised"] * HDL.TG
    AAT <- scale(AAT)
  } else {
    log_AAT <- AAT_coef["intercept", "raw"] +
      AAT_coef["GlycA", "raw"] * log(check_range(GlycA, "GlycA", range_check)) +
      AAT_coef["FAw3", "raw"] * log(check_range(FAw3, "FAw3", range_check)) +
      AAT_coef["VLDL.D", "raw"] * log(check_range(VLDL.D, "VLDL.D", range_check)) +
      AAT_coef["HDL3.C", "raw"] * log(check_range(HDL3.C, "HDL3.C", range_check)) +
      AAT_coef["LDL.D", "raw"] * log(check_range(LDL.D, "LDL.D", range_check)) +
      AAT_coef["Phe", "raw"] * log(check_range(Phe, "Phe", range_check)) +
      AAT_coef["Leu", "raw"] * log(check_range(Leu, "Leu", range_check)) +
      AAT_coef["ApoB", "raw"] * log(check_range(ApoB, "ApoB", range_check)) +
      AAT_coef["Alb", "raw"] * log(check_range(Alb, "Alb", range_check)) +
      AAT_coef["Tyr", "raw"] * log(check_range(Tyr, "Tyr", range_check)) +
      AAT_coef["bOHBut", "raw"] * log(check_range(bOHBut, "bOHBut", range_check)) +
      AAT_coef["BMI", "raw"] * log(check_range(BMI, "BMI", range_check)) +
      AAT_coef["Ala", "raw"] * log(check_range(Ala, "Ala", range_check)) +
      AAT_coef["L.HDL.TG", "raw"] * log(check_range(L.HDL.TG, "L.HDL.TG", range_check)) +
      AAT_coef["Ile", "raw"] * log(check_range(Ile, "Ile", range_check)) +
      AAT_coef["Ace", "raw"] * log(check_range(Ace, "Ace", range_check)) +
      AAT_coef["His", "raw"] * log(check_range(His, "His", range_check)) +
      AAT_coef["HDL.TG", "raw"] * log(check_range(HDL.TG, "HDL.TG", range_check))

    # Transform back to raw concentration units
    AAT <- exp(log_AAT)

    # Remove imputed concentrations that are outside the range of concentration
    # values observed in the model training data.
    AAT <- check_range(AAT, "AAT", range_check)
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
#'
#' @return A vector of AGP measurements ranging between 362--1,880 mg/L or
#'  measurements standardised to the population if \code{standardised = 'TRUE'}.
#'
#' @export
impute_AGP <- function(
  GlycA, TotFA, IDL.FC, L.HDL.FC, His, HDL.TG, BMI, S.HDL.FC, S.LDL.TG, bOHBut,
  LA, S.HDL.CE, Lac, S.VLDL.TG, Ace, Cit, SFA, Ala, XXL.VLDL.CE, Glol, Age,
  Crea, Gly, range_check=TRUE, standardised=FALSE
) {
  if (standardised) {
    AGP <- AGP_coef["intercept", "standardised"] +
      AGP_coef["GlycA", "standardised"] * GlycA +
      AGP_coef["TotFA", "standardised"] * TotFA +
      AGP_coef["IDL.FC", "standardised"] * IDL.FC +
      AGP_coef["L.HDL.FC", "standardised"] * L.HDL.FC +
      AGP_coef["His", "standardised"] * His +
      AGP_coef["HDL.TG", "standardised"] * HDL.TG +
      AGP_coef["BMI", "standardised"] * BMI +
      AGP_coef["S.HDL.FC", "standardised"] * S.HDL.FC +
      AGP_coef["S.LDL.TG", "standardised"] * S.LDL.TG +
      AGP_coef["bOHBut", "standardised"] * bOHBut +
      AGP_coef["LA", "standardised"] * LA +
      AGP_coef["S.HDL.CE", "standardised"] * S.HDL.CE +
      AGP_coef["Lac", "standardised"] * Lac +
      AGP_coef["S.VLDL.TG", "standardised"] * S.VLDL.TG +
      AGP_coef["Ace", "standardised"] * Ace +
      AGP_coef["Cit", "standardised"] * Cit +
      AGP_coef["SFA", "standardised"] * SFA +
      AGP_coef["Ala", "standardised"] * Ala +
      AGP_coef["XXL.VLDL.CE", "standardised"] * XXL.VLDL.CE +
      AGP_coef["Glol", "standardised"] * Glol +
      AGP_coef["Age", "standardised"] * Age +
      AGP_coef["Crea", "standardised"] * Crea +
      AGP_coef["Gly", "standardised"] * Gly
    AGP <- scale(AGP)
  } else {
    log_AGP <- AGP_coef["intercept", "raw"] +
      AGP_coef["GlycA", "raw"] * log(check_range(GlycA, "GlycA", range_check)) +
      AGP_coef["TotFA", "raw"] * log(check_range(TotFA, "TotFA", range_check)) +
      AGP_coef["IDL.FC", "raw"] * log(check_range(IDL.FC, "IDL.FC", range_check)) +
      AGP_coef["L.HDL.FC", "raw"] * log(check_range(L.HDL.FC, "L.HDL.FC", range_check)) +
      AGP_coef["His", "raw"] * log(check_range(His, "His", range_check)) +
      AGP_coef["HDL.TG", "raw"] * log(check_range(HDL.TG, "HDL.TG", range_check)) +
      AGP_coef["BMI", "raw"] * log(check_range(BMI, "BMI", range_check)) +
      AGP_coef["S.HDL.FC", "raw"] * log(check_range(S.HDL.FC, "S.HDL.FC", range_check)) +
      AGP_coef["S.LDL.TG", "raw"] * log(check_range(S.LDL.TG, "S.LDL.TG", range_check)) +
      AGP_coef["bOHBut", "raw"] * log(check_range(bOHBut, "bOHBut", range_check)) +
      AGP_coef["LA", "raw"] * log(check_range(LA, "LA", range_check)) +
      AGP_coef["S.HDL.CE", "raw"] * log(check_range(S.HDL.CE, "S.HDL.CE", range_check)) +
      AGP_coef["Lac", "raw"] * log(check_range(Lac, "Lac", range_check)) +
      AGP_coef["S.VLDL.TG", "raw"] * log(check_range(S.VLDL.TG, "S.VLDL.TG", range_check)) +
      AGP_coef["Ace", "raw"] * log(check_range(Ace, "Ace", range_check)) +
      AGP_coef["Cit", "raw"] * log(check_range(Cit, "Cit", range_check)) +
      AGP_coef["SFA", "raw"] * log(check_range(SFA, "SFA", range_check)) +
      AGP_coef["Ala", "raw"] * log(check_range(Ala, "Ala", range_check)) +
      AGP_coef["XXL.VLDL.CE", "raw"] * log(check_range(XXL.VLDL.CE, "XXL.VLDL.CE", range_check)) +
      AGP_coef["Glol", "raw"] * log(check_range(Glol, "Glol", range_check)) +
      AGP_coef["Age", "raw"] * check_range(Age, "Age", range_check) +
      AGP_coef["Crea", "raw"] * log(check_range(Crea, "Crea", range_check)) +
      AGP_coef["Gly", "raw"] * log(check_range(Gly, "Gly", range_check))

    # Transform back to raw concentration units
    AGP <- exp(log_AGP)

    # Remove imputed concentrations that are outside the range of concentration
    # values observed in the model training data.
    AGP <- check_range(AGP, "AGP", range_check)
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
#'
#' @return A vector of HP measurements ranging between 0.14--3.95 mg/L or
#'  measurements standardised to the population if \code{standardised = 'TRUE'}.
#'
#' @export
impute_HP <- function(
  GlycA, LA, IDL.FC, SM, FAw3, HDL.TG, S.VLDL.CE, Age, Alb, Ile, Cit, VLDL.D,
  Leu, Val, L.VLDL.CE, Pyr, Lac, Gln, M.HDL.FC, XL.HDL.TG, XL.HDL.PL, His, Tyr,
  BMI, L.HDL.TG, PUFA, S.LDL.FC, range_check=TRUE, standardised=FALSE
) {
  if (standardised) {
    HP <- HP_coef["intercept", "standardised"] +
      HP_coef["GlycA", "standardised"] * GlycA +
      HP_coef["LA", "standardised"] * LA +
      HP_coef["IDL.FC", "standardised"] * IDL.FC +
      HP_coef["SM", "standardised"] * SM +
      HP_coef["FAw3", "standardised"] * FAw3 +
      HP_coef["HDL.TG", "standardised"] * HDL.TG +
      HP_coef["S.VLDL.CE", "standardised"] * S.VLDL.CE +
      HP_coef["Age", "standardised"] * Age +
      HP_coef["Alb", "standardised"] * Alb +
      HP_coef["Ile", "standardised"] * Ile +
      HP_coef["Cit", "standardised"] * Cit +
      HP_coef["VLDL.D", "standardised"] * VLDL.D +
      HP_coef["Leu", "standardised"] * Leu +
      HP_coef["Val", "standardised"] * Val +
      HP_coef["L.VLDL.CE", "standardised"] * L.VLDL.CE +
      HP_coef["Pyr", "standardised"] * Pyr +
      HP_coef["Lac", "standardised"] * Lac +
      HP_coef["Gln", "standardised"] * Gln +
      HP_coef["M.HDL.FC", "standardised"] * M.HDL.FC +
      HP_coef["XL.HDL.TG", "standardised"] * XL.HDL.TG +
      HP_coef["XL.HDL.PL", "standardised"] * XL.HDL.PL +
      HP_coef["His", "standardised"] * His +
      HP_coef["Tyr", "standardised"] * Tyr +
      HP_coef["BMI", "standardised"] * BMI +
      HP_coef["L.HDL.TG", "standardised"] * L.HDL.TG +
      HP_coef["PUFA", "standardised"] * PUFA +
      HP_coef["S.LDL.FC", "standardised"] * S.LDL.FC
    HP <- scale(HP)
  } else {
    log_HP <- HP_coef["intercept", "raw"] +
      HP_coef["GlycA", "raw"] * log(check_range(GlycA, "GlycA", range_check)) +
      HP_coef["LA", "raw"] * log(check_range(LA, "LA", range_check)) +
      HP_coef["IDL.FC", "raw"] * log(check_range(IDL.FC, "IDL.FC", range_check)) +
      HP_coef["SM", "raw"] * log(check_range(SM, "SM", range_check)) +
      HP_coef["FAw3", "raw"] * log(check_range(FAw3, "FAw3", range_check)) +
      HP_coef["HDL.TG", "raw"] * log(check_range(HDL.TG, "HDL.TG", range_check)) +
      HP_coef["S.VLDL.CE", "raw"] * log(check_range(S.VLDL.CE, "S.VLDL.CE", range_check)) +
      HP_coef["Age", "raw"] * check_range(Age, "Age", range_check) +
      HP_coef["Alb", "raw"] * log(check_range(Alb, "Alb", range_check)) +
      HP_coef["Ile", "raw"] * log(check_range(Ile, "Ile", range_check)) +
      HP_coef["Cit", "raw"] * log(check_range(Cit, "Cit", range_check)) +
      HP_coef["VLDL.D", "raw"] * log(check_range(VLDL.D, "VLDL.D", range_check)) +
      HP_coef["Leu", "raw"] * log(check_range(Leu, "Leu", range_check)) +
      HP_coef["Val", "raw"] * log(check_range(Val, "Val", range_check)) +
      HP_coef["L.VLDL.CE", "raw"] * log(check_range(L.VLDL.CE, "L.VLDL.CE", range_check)) +
      HP_coef["Pyr", "raw"] * log(check_range(Pyr, "Pyr", range_check)) +
      HP_coef["Lac", "raw"] * log(check_range(Lac, "Lac", range_check)) +
      HP_coef["Gln", "raw"] * log(check_range(Gln, "Gln", range_check)) +
      HP_coef["M.HDL.FC", "raw"] * log(check_range(M.HDL.FC, "M.HDL.FC", range_check)) +
      HP_coef["XL.HDL.TG", "raw"] * log(check_range(XL.HDL.TG, "XL.HDL.TG", range_check)) +
      HP_coef["XL.HDL.PL", "raw"] * log(check_range(XL.HDL.PL, "XL.HDL.PL", range_check)) +
      HP_coef["His", "raw"] * log(check_range(His, "His", range_check)) +
      HP_coef["Tyr", "raw"] * log(check_range(Tyr, "Tyr", range_check)) +
      HP_coef["BMI", "raw"] * log(check_range(BMI, "BMI", range_check)) +
      HP_coef["L.HDL.TG", "raw"] * log(check_range(L.HDL.TG, "L.HDL.TG", range_check)) +
      HP_coef["PUFA", "raw"] * log(check_range(PUFA, "PUFA", range_check)) +
      HP_coef["S.LDL.FC", "raw"] * log(check_range(S.LDL.FC, "S.LDL.FC", range_check))

    # Transform back to raw concentration units
    HP <- exp(log_HP)

    # Remove imputed concentrations that are outside the range of concentration
    # values observed in the model training data.
    HP <- check_range(HP, "HP", range_check)
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
#'
#' @return A vector of TF measurements ranging between 1.39--4.38 mg/L or
#'  measurements standardised to the population if \code{standardised = 'TRUE'}.
#'
impute_TF <- function(GlycA, Sex, Age, S.HDL.FC, Ace, Ala, SFA, His, Gln,
                      range_check=TRUE, standardised=FALSE) {

  if (standardised) {
    TF <- TF_coef["intercept", "standardised"] +
      TF_coef["GlycA", "standardised"] * GlycA +
      TF_coef["Sex", "standardised"] * Sex +
      TF_coef["Age", "standardised"] * Age +
      TF_coef["S.HDL.FC", "standardised"] * S.HDL.FC +
      TF_coef["Ace", "standardised"] * Ace +
      TF_coef["Ala", "standardised"] * Ala +
      TF_coef["SFA", "standardised"] * SFA +
      TF_coef["His", "standardised"] * His +
      TF_coef["Gln", "standardised"] * Gln
    TF <- scale(TF)
  } else {
    log_TF <- TF_coef["intercept", "raw"] +
      TF_coef["GlycA", "raw"] * log(check_range(GlycA, "GlycA", range_check)) +
      TF_coef["Sex", "raw"] * check_range(Sex, "Sex", range_check) +
      TF_coef["Age", "raw"] * check_range(Age, "Age", range_check) +
      TF_coef["S.HDL.FC", "raw"] * log(check_range(S.HDL.FC, "S.HDL.FC", range_check)) +
      TF_coef["Ace", "raw"] * log(check_range(Ace, "Ace", range_check)) +
      TF_coef["Ala", "raw"] * log(check_range(Ala, "Ala", range_check)) +
      TF_coef["SFA", "raw"] * log(check_range(SFA, "SFA", range_check)) +
      TF_coef["His", "raw"] * log(check_range(His, "His", range_check)) +
      TF_coef["Gln", "raw"] * log(check_range(Gln, "Gln", range_check))

    # Transform back to raw concentration units
    TF <- exp(log_TF)

    # Remove imputed concentrations that are outside the range of concentration
    # values observed in the model training data.
    TF <- check_range(TF, "TF", range_check)
  }

  message("Successfully imputed TF for ", sum(!is.na(TF)), " samples")
  warning("The concentrations predicted by this function are unlikely to reflect true TF levels", immediate.=TRUE)
  return(TF)
}

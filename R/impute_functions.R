#' Impute Alpha-1 antitrypsin (A1AT)
#'
#' Imputes the concentrations of A1AT from serum NMR measurements in healthy
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
#'
#' @return A vector of A1AT measurements ranging between 0.64--2.58 mg/L.
#'
#' @export
impute_A1AT <- function(
  GlycA, FAw3, VLDL.D, HDL3.C, LDL.D, Phe, Leu, ApoB, Alb, Tyr, bOHBut, BMI,
  Ala, L.HDL.TG, Ile, Ace, His, HDL.TG
) {
  log_A1AT <- A1AT_coef["intercept"] +
    A1AT_coef["GlycA"] * log(check_range(GlycA, "GlycA")) +
    A1AT_coef["FAw3"] * log(check_range(FAw3, "FAw3")) +
    A1AT_coef["VLDL.D"] * log(check_range(VLDL.D, "VLDL.D")) +
    A1AT_coef["HDL3.C"] * log(check_range(HDL3.C, "HDL3.C")) +
    A1AT_coef["LDL.D"] * log(check_range(LDL.D, "LDL.D")) +
    A1AT_coef["Phe"] * log(check_range(Phe, "Phe")) +
    A1AT_coef["Leu"] * log(check_range(Leu, "Leu")) +
    A1AT_coef["ApoB"] * log(check_range(ApoB, "ApoB")) +
    A1AT_coef["Alb"] * log(check_range(Alb, "Alb")) +
    A1AT_coef["Tyr"] * log(check_range(Tyr, "Tyr")) +
    A1AT_coef["bOHBut"] * log(check_range(bOHBut, "bOHBut")) +
    A1AT_coef["BMI"] * log(check_range(BMI, "BMI")) +
    A1AT_coef["Ala"] * log(check_range(Ala, "Ala")) +
    A1AT_coef["L.HDL.TG"] * log(check_range(L.HDL.TG, "L.HDL.TG")) +
    A1AT_coef["Ile"] * log(check_range(Ile, "Ile")) +
    A1AT_coef["Ace"] * log(check_range(Ace, "Ace")) +
    A1AT_coef["His"] * log(check_range(His, "His")) +
    A1AT_coef["HDL.TG"] * log(check_range(HDL.TG, "HDL.TG"))

  # Transform back to raw concentration units
  A1AT <- exp(log_A1AT)

  # Remove imputed concentrations that are outside the range of concentration
  # values observed in the model training data.
  A1AT <- check_range(A1AT, "A1AT")

  message("Successfully imputed A1AT for ", sum(!is.na(A1AT)), " samples")

  return(A1AT)
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
#'
#' @return A vector of AGP measurements ranging between 362--1,880 mg/L.
#'
#' @export
impute_AGP <- function(
  GlycA, TotFA, IDL.FC, L.HDL.FC, His, HDL.TG, BMI, S.HDL.FC, S.LDL.TG, bOHBut,
  LA, S.HDL.CE, Lac, S.VLDL.TG, Ace, Cit, SFA, Ala, XXL.VLDL.CE, Glol, Age,
  Crea, Gly
) {
  log_AGP <- AGP_coef["intercept"] +
    AGP_coef["GlycA"] * log(check_range(GlycA, "GlycA")) +
    AGP_coef["TotFA"] * log(check_range(TotFA, "TotFA")) +
    AGP_coef["IDL.FC"] * log(check_range(IDL.FC, "IDL.FC")) +
    AGP_coef["L.HDL.FC"] * log(check_range(L.HDL.FC, "L.HDL.FC")) +
    AGP_coef["His"] * log(check_range(His, "His")) +
    AGP_coef["HDL.TG"] * log(check_range(HDL.TG, "HDL.TG")) +
    AGP_coef["BMI"] * log(check_range(BMI, "BMI")) +
    AGP_coef["S.HDL.FC"] * log(check_range(S.HDL.FC, "S.HDL.FC")) +
    AGP_coef["S.LDL.TG"] * log(check_range(S.LDL.TG, "S.LDL.TG")) +
    AGP_coef["bOHBut"] * log(check_range(bOHBut, "bOHBut")) +
    AGP_coef["LA"] * log(check_range(LA, "LA")) +
    AGP_coef["S.HDL.CE"] * log(check_range(S.HDL.CE, "S.HDL.CE")) +
    AGP_coef["Lac"] * log(check_range(Lac, "Lac")) +
    AGP_coef["S.VLDL.TG"] * log(check_range(S.VLDL.TG, "S.VLDL.TG")) +
    AGP_coef["Ace"] * log(check_range(Ace, "Ace")) +
    AGP_coef["Cit"] * log(check_range(Cit, "Cit")) +
    AGP_coef["SFA"] * log(check_range(SFA, "SFA")) +
    AGP_coef["Ala"] * log(check_range(Ala, "Ala")) +
    AGP_coef["XXL.VLDL.CE"] * log(check_range(XXL.VLDL.CE, "XXL.VLDL.CE")) +
    AGP_coef["Glol"] * log(check_range(Glol, "Glol")) +
    AGP_coef["Age"] * check_range(Age, "Age") +
    AGP_coef["Crea"] * log(check_range(Crea, "Crea")) +
    AGP_coef["Gly"] * log(check_range(Gly, "Gly"))

  # Transform back to raw concentration units
  AGP <- exp(log_AGP)

  # Remove imputed concentrations that are outside the range of concentration
  # values observed in the model training data.
  AGP <- check_range(AGP, "AGP")

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
#'
#' @return A vector of HP measurements ranging between 0.14--3.95 mg/L.
#'
#' @export
impute_HP <- function(
  GlycA, LA, IDL.FC, SM, FAw3, HDL.TG, S.VLDL.CE, Age, Alb, Ile, Cit, VLDL.D,
  Leu, Val, L.VLDL.CE, Pyr, Lac, Gln, M.HDL.FC, XL.HDL.TG, XL.HDL.PL, His, Tyr,
  BMI, L.HDL.TG, PUFA, S.LDL.FC
) {
  log_HP <- HP_coef["intercept"] +
    HP_coef["GlycA"] * log(check_range(GlycA, "GlycA")) +
    HP_coef["LA"] * log(check_range(LA, "LA")) +
    HP_coef["IDL.FC"] * log(check_range(IDL.FC, "IDL.FC")) +
    HP_coef["SM"] * log(check_range(SM, "SM")) +
    HP_coef["FAw3"] * log(check_range(FAw3, "FAw3")) +
    HP_coef["HDL.TG"] * log(check_range(HDL.TG, "HDL.TG")) +
    HP_coef["S.VLDL.CE"] * log(check_range(S.VLDL.CE, "S.VLDL.CE")) +
    HP_coef["Age"] * check_range(Age, "Age") +
    HP_coef["Alb"] * log(check_range(Alb, "Alb")) +
    HP_coef["Ile"] * log(check_range(Ile, "Ile")) +
    HP_coef["Cit"] * log(check_range(Cit, "Cit")) +
    HP_coef["VLDL.D"] * log(check_range(VLDL.D, "VLDL.D")) +
    HP_coef["Leu"] * log(check_range(Leu, "Leu")) +
    HP_coef["Val"] * log(check_range(Val, "Val")) +
    HP_coef["L.VLDL.CE"] * log(check_range(L.VLDL.CE, "L.VLDL.CE")) +
    HP_coef["Pyr"] * log(check_range(Pyr, "Pyr")) +
    HP_coef["Lac"] * log(check_range(Lac, "Lac")) +
    HP_coef["Gln"] * log(check_range(Gln, "Gln")) +
    HP_coef["M.HDL.FC"] * log(check_range(M.HDL.FC, "M.HDL.FC")) +
    HP_coef["XL.HDL.TG"] * log(check_range(XL.HDL.TG, "XL.HDL.TG")) +
    HP_coef["XL.HDL.PL"] * log(check_range(XL.HDL.PL, "XL.HDL.PL")) +
    HP_coef["His"] * log(check_range(His, "His")) +
    HP_coef["Tyr"] * log(check_range(Tyr, "Tyr")) +
    HP_coef["BMI"] * log(check_range(BMI, "BMI")) +
    HP_coef["L.HDL.TG"] * log(check_range(L.HDL.TG, "L.HDL.TG")) +
    HP_coef["PUFA"] * log(check_range(PUFA, "PUFA")) +
    HP_coef["S.LDL.FC"] * log(check_range(S.LDL.FC, "S.LDL.FC"))

  # Transform back to raw concentration units
  HP <- exp(log_HP)

  # Remove imputed concentrations that are outside the range of concentration
  # values observed in the model training data.
  HP <- check_range(HP, "HP")

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
#'
#' @return A vector of TF measurements ranging between 1.39--4.38 mg/L.
#'
#' @export
impute_TF <- function(GlycA, Sex, Age, S.HDL.FC, Ace, Ala, SFA, His, Gln) {
  log_TF <- TF_coef["intercept"] +
    TF_coef["GlycA"] * log(check_range(GlycA, "GlycA")) +
    TF_coef["Sex"] * check_range(Sex, "Sex") +
    TF_coef["Age"] * check_range(Age, "Age") +
    TF_coef["S.HDL.FC"] * log(check_range(S.HDL.FC, "S.HDL.FC")) +
    TF_coef["Ace"] * log(check_range(Ace, "Ace")) +
    TF_coef["Ala"] * log(check_range(Ala, "Ala")) +
    TF_coef["SFA"] * log(check_range(SFA, "SFA")) +
    TF_coef["His"] * log(check_range(His, "His")) +
    TF_coef["Gln"] * log(check_range(Gln, "Gln"))

  # Transform back to raw concentration units
  TF <- exp(log_TF)

  # Remove imputed concentrations that are outside the range of concentration
  # values observed in the model training data.
  TF <- check_range(TF, "TF")

  message("Successfully imputed TF for ", sum(!is.na(TF)), " samples")

  return(TF)
}

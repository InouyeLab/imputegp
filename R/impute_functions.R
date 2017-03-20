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
  # Model in the paper (coefficients rounded):
  # A1AT = −11 + 0.68 GlycA − 0.094 FAw3 − 0.8 VLDL-D + 0.37 HDL3-C + 4.2 LDL-D
  #            + 0.14 Phe − 0.08 Leu − 0.064 ApoB − 0.19 Alb − 0.045 Tyr
  #            + 0.019 bOHBut − 0.045 BMI - 0.043 Ala + 0.011 L-HDL-TG
  #            − 0.022 Ile − 0.025 Ace − 0.038 His + 0.0024 HDL-TG
  log_A1AT <- -10.76239 +
    0.676036663 * log(check_range(GlycA, "GlycA")) -
    0.093856499 * log(check_range(FAw3, "FAw3")) -
    0.801330272 * log(check_range(VLDL.D, "VLDL.D")) +
    0.369888369 * log(check_range(HDL3.C, "HDL3.C")) +
    4.209011261 * log(check_range(LDL.D, "LDL.D")) +
    0.142793364 * log(check_range(Phe, "Phe")) -
    0.080290581 * log(check_range(Leu, "Leu")) -
    0.063905074 * log(check_range(ApoB, "ApoB")) -
    0.193478006 * log(check_range(Alb, "Alb")) -
    0.044871373 * log(check_range(Tyr, "Tyr")) +
    0.018540450 * log(check_range(bOHBut, "bOHBut")) -
    0.044926761 * log(check_range(BMI, "BMI")) -
    0.042894831 * log(check_range(Ala, "Ala")) +
    0.011385043 * log(check_range(L.HDL.TG, "L.HDL.TG")) -
    0.021748967 * log(check_range(Ile, "Ile")) -
    0.025021140 * log(check_range(Ace, "Ace")) -
    0.038414793 * log(check_range(His, "His")) +
    0.002385282 * log(check_range(HDL.TG, "HDL.TG"))

  # Transform back to raw concentration units
  A1AT <- exp(log_A1AT)

  # Remove imputed concentrations that are outside the range of concentration
  # values observed in the model training data.
  A1AT <- check_range(A1AT, "A1AT")

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
  # Model in the paper (coefficients rounded):
  # AGP = 5.7 + 1.7 GlycA − 0.34 TotFA + 0.22 IDL-FC + 0.043 L-HDL-FC − 0.22 His
  #           − 0.092 HDL-TG + 0.11 BMI − 0.16 S-HDL-FC − 0.057 S-LDL-TG
  #           − 0.023 bOHBut − 0.062 LA + 0.066 S-HDL-CE − 0.055 Lac
  #           − 0.021 S-VLDL-TG − 0.042 Ace − 0.052 Cit − 0.038 SFA
  #           − 0.049 Ala + 0.0013 XXL-VLDL-CE + 0.014 Glol + 0.00017 Age
  #           + 0.0084 Crea + 0.0061 Gly
  log_AGP <- 5.686557 +
    1.6835689376 * log(check_range(GlycA, "GlycA")) -
    0.3423080143 * log(check_range(TotFA, "TotFA")) +
    0.2171134779 * log(check_range(IDL.FC, "IDL.FC")) +
    0.0425565289 * log(check_range(L.HDL.FC, "L.HDL.FC")) -
    0.2194217249 * log(check_range(His, "His")) -
    0.0915754435 * log(check_range(HDL.TG, "HDL.TG")) +
    0.1141298188 * log(check_range(BMI, "BMI")) -
    0.1559315573 * log(check_range(S.HDL.FC, "S.HDL.FC")) -
    0.0568963001 * log(check_range(S.LDL.TG, "S.LDL.TG")) -
    0.0232790321 * log(check_range(bOHBut, "bOHBut")) -
    0.0623795538 * log(check_range(LA, "LA")) +
    0.0661691313 * log(check_range(S.HDL.CE, "S.HDL.CE")) -
    0.0548688667 * log(check_range(Lac, "Lac")) -
    0.0208627226 * log(check_range(S.LDL.TG, "S.LDL.TG")) -
    0.0416815172 * log(check_range(Ace, "Ace")) -
    0.0520885532 * log(check_range(Cit, "Cit")) -
    0.0378144874 * log(check_range(SFA, "SFA")) -
    0.0487374784 * log(check_range(Ala, "Ala")) +
    0.0013157870 * log(check_range(XXL.VLDL.CE, "XXL.VLDL.CE")) +
    0.0141442945 * log(check_range(Glol, "Glol")) +
    0.0001718169 * check_range(Age, "Age") +
    0.0084351085 * log(check_range(Crea, "Crea")) +
    0.0060648743 * log(check_range(Gly, "Gly"))

  # Transform back to raw concentration units
  AGP <- exp(log_AGP)

  # Remove imputed concentrations that are outside the range of concentration
  # values observed in the model training data.
  AGP <- check_range(AGP, "AGP")

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
  # Model in the paper (coefficients rounded):
  # HP = 0.53 + 3.0 GlycA − 0.82 LA + 0.36 IDL-FC + 0.48 SM − 0.26 FAw3
  #           − 0.29 HDL-TG − 0.068 S-VLDL-CE + 0.003 Age − 0.75 Alb − 0.13 Ile
  #           − 0.24 Cit − 0.89 VLDL-D − 0.13 Leu + 0.13 Val + 0.0066 L-VLDL-CE
  #           + 0.084 Pyr − 0.084 Lac − 0.094 Gln + 0.042 M-HDL-FC
  #           − 0.011 XL-HDL-TG + 0.0089 XL-HDL-PL − 0.048 His − 0.026 Tyr
  #           − 0.024 BMI − 0.0037 L-HDL-TG − 0.004 PUFA + 0.00056 S-LDL-FC
  log_HP <- 0.5276596 +
    3.039616949 * log(check_range(GlycA, "GlycA")) -
    0.824141412 * log(check_range(LA, "LA")) +
    0.364675335 * log(check_range(IDL.FC, "IDL.FC")) +
    0.479815330 * log(check_range(SM, "SM")) -
    0.256070355 * log(check_range(FAw3, "FAw3")) -
    0.290408926 * log(check_range(HDL.TG, "HDL.TG")) -
    0.068466341 * log(check_range(S.VLDL.CE, "S.VLDL.CE")) +
    0.003019050 * check_range(Age, "Age") -
    0.749877375 * log(check_range(Alb, "Alb")) -
    0.134308479 * log(check_range(Ile, "Ile")) -
    0.235930832 * log(check_range(Cit, "Cit")) -
    0.887777970 * log(check_range(VLDL.D, "VLDL.D")) -
    0.132830119 * log(check_range(Leu, "Leu")) +
    0.132092539 * log(check_range(Val, "Val")) +
    0.006579148 * log(check_range(L.VLDL.CE, "L.VLDL.CE")) +
    0.083861053 * log(check_range(Pyr, "Pyr")) -
    0.084077503 * log(check_range(Lac, "Lac")) -
    0.093876184 * log(check_range(Gln, "Gln")) +
    0.041966889 * log(check_range(M.HDL.FC, "M.HDL.FC")) -
    0.011330501 * log(check_range(XL.HDL.TG, "XL.HDL.TG")) +
    0.008918302 * log(check_range(XL.HDL.PL, "XL.HDL.PL")) -
    0.048117394 * log(check_range(His, "His")) -
    0.025676721 * log(check_range(Tyr, "Tyr")) -
    0.024479223 * log(check_range(BMI, "BMI")) -
    0.003680147 * log(check_range(L.HDL.TG, "L.HDL.TG")) -
    0.003989682 * log(check_range(PUFA, "PUFA")) +
    0.000558750 * log(check_range(S.LDL.FC, "S.LDL.FC"))

  # Transform back to raw concentration units
  HP <- exp(log_HP)

  # Remove imputed concentrations that are outside the range of concentration
  # values observed in the model training data.
  HP <- check_range(HP, "HP")

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
  # Model in the paper (coefficients rounded):
  # TF = 1.1 + 0.14 GlycA + 0.032 Sex − 0.0010 Age + 0.090 S-HDL-FC − 0.037 Ace
  #          + 0.039 Ala + 0.024 SFA + 0.013 His − 0.0097 Gln
  log_TF <- 1.058351 +
    0.139028225 * log(check_range(GlycA, "GlycA")) +
    0.032289249 * check_range(Sex, "Sex") -
    0.001008589 * check_range(Age, "Age") +
    0.089914967 * log(check_range(S.HDL.FC, "S.HDL.FC")) -
    0.037428622 * log(check_range(Ace, "Ace")) +
    0.038566477 * log(check_range(Ala, "Ala")) +
    0.023512582 * log(check_range(SFA, "SFA")) +
    0.012685158 * log(check_range(His, "His")) -
    0.009667614 * log(check_range(Gln, "Gln"))

  # Transform back to raw concentration units
  TF <- exp(log_TF)

  # Remove imputed concentrations that are outside the range of concentration
  # values observed in the model training data.
  TF <- check_range(TF, "TF")

  return(TF)
}

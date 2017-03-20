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

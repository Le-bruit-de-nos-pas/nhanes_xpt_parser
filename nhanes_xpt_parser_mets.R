library(tidyverse)
library(data.table)
library(survey)
options(scipen = 999)
options(survey.lonely.psu='adjust')

cat("R package versions:\n")

for (p in c("base", "survey","dplyr", "tidyverse", "data.table")) { 
  cat(p, ": ", as.character(packageVersion(p)), "\n")
}

# DEMO ***************************************
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEMO_L.XPT", tf <- tempfile(), mode="wb")
Demographcics_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_20211_2022 %>% select(SEQN, RIDAGEYR) %>% filter(RIDAGEYR>=18)


 # How to Import all data -------------------

# EXAM ***************************************
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BMX_L.XPT", tf <- tempfile(), mode="wb")
Body_Measures_20211_2022 <- foreign::read.xport(tf)[,]

Body_Measures_20211_2022 %>% select(BMXBMI) %>% summarise(mean=mean(BMXBMI,na.rm=T))

Body_Measures_20211_2022 %>% select(SEQN, BMXBMI) %>%
  inner_join(Adults) %>%
  drop_na() %>% mutate(BMXBMI=ifelse(BMXBMI>=30,1,0)) %>% group_by(BMXBMI) %>% count()


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BPXO_L.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/LUX_L.XPT", tf <- tempfile(), mode="wb")
Elastography_20211_2022 <- foreign::read.xport(tf)[,]


# LABS ***************************************
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/AGP_L.XPT", tf <- tempfile(), mode="wb")
alpha_Acid_Glycoprotein_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HDL_L.XPT", tf <- tempfile(), mode="wb")
Cholesterol_High_Density_Lipoprotein_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/TCHOL_L.XPT", tf <- tempfile(), mode="wb")
Cholesterol_Total_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/CBC_L.XPT", tf <- tempfile(), mode="wb")
Complete_Blood_Count_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/FASTQX_L.XPT", tf <- tempfile(), mode="wb")
Fasting_Questionnaire_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/FERTIN_L.XPT", tf <- tempfile(), mode="wb")
Ferritin_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/FOLATE_L.XPT", tf <- tempfile(), mode="wb")
Folate_RBC_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/GHB_L.XPT", tf <- tempfile(), mode="wb")
Glycohemoglobin_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HEPA_L.XPT", tf <- tempfile(), mode="wb")
Hepatitis_A_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HEPB_S_L.XPT", tf <- tempfile(), mode="wb")
Hepatitis_B_Surface_Antibody_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HSCRP_L.XPT", tf <- tempfile(), mode="wb")
High_Sensitivity_C_Reactive_Protein_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/INS_L.XPT", tf <- tempfile(), mode="wb")
Insulin_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/PBCD_L.XPT", tf <- tempfile(), mode="wb")
Lead_Cadmium_Total_Mercury_Selenium_Manganese_Blood_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/IHGEM_L.XPT", tf <- tempfile(), mode="wb")
Mercury_Inorganic_Ethyl_and_Methyl_Blood_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/GLU_L.XPT", tf <- tempfile(), mode="wb")
Plasma_Fasting_Glucose_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/FOLFMS_L.XPT", tf <- tempfile(), mode="wb")
Serum_Folate_Forms_Total_Individual_Serum_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/TFR_L.XPT", tf <- tempfile(), mode="wb")
Transferrin_Receptor_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/UCPREG_L.XPT", tf <- tempfile(), mode="wb")
Urine_Pregnancy_Test_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/VID_L.XPT", tf <- tempfile(), mode="wb")
Vitamin_D_20211_2022 <- foreign::read.xport(tf)[,]


# Questionnaire  ***************************************
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/AUQ_L.XPT", tf <- tempfile(), mode="wb")
Audiometry_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/ACQ_L.XPT", tf <- tempfile(), mode="wb")
Acculturation_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/ALQ_L.XPT", tf <- tempfile(), mode="wb")
Alcohol_Use_Total_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BPQ_L.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_Cholesterol_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HSQ_L.XPT", tf <- tempfile(), mode="wb")
Current_Health_Status_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEQ_L.XPT", tf <- tempfile(), mode="wb")
Dermatology_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DIQ_L.XPT", tf <- tempfile(), mode="wb")
Diabetes_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DBQ_L.XPT", tf <- tempfile(), mode="wb")
Diet_Behavior_Nutrition_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/ECQ_L.XPT", tf <- tempfile(), mode="wb")
Early_Childhood_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/FNQ_L.XPT", tf <- tempfile(), mode="wb")
Functioning_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HIQ_L.XPT", tf <- tempfile(), mode="wb")
Health_Insurance_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HEQ_L.XPT", tf <- tempfile(), mode="wb")
Hepatitis_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HUQ_L.XPT", tf <- tempfile(), mode="wb")
Hospital_Utilization_Access_to_Care_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HOQ_L.XPT", tf <- tempfile(), mode="wb")
Housing_Characteristics_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/IMQ_L.XPT", tf <- tempfile(), mode="wb")
Immunization_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/INQ_L.XPT", tf <- tempfile(), mode="wb")
Income_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/KIQ_U_L.XPT", tf <- tempfile(), mode="wb")
Kidney_Conditions_Urology_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/MCQ_L.XPT", tf <- tempfile(), mode="wb")
Medical_Conditions_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DPQ_L.XPT", tf <- tempfile(), mode="wb")
Mental_Health_Depression_Screener_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/OCQ_L.XPT", tf <- tempfile(), mode="wb")
Occupation_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/OHQ_L.XPT", tf <- tempfile(), mode="wb")
Oral_Health_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/PUQMEC_L.XPT", tf <- tempfile(), mode="wb")
Pesticide_Use_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/PAQ_L.XPT", tf <- tempfile(), mode="wb")
Physical_Activity_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/PAQY_L.XPT", tf <- tempfile(), mode="wb")
Physical_Activity_Youth_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/RXQ_RX_L.XPT", tf <- tempfile(), mode="wb")
Prescription_Medications_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/RXQASA_L.XPT", tf <- tempfile(), mode="wb")
Preventive_Aspirin_Use_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/RHQ_L.XPT", tf <- tempfile(), mode="wb")
Reproductive_Health_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/SLQ_L.XPT", tf <- tempfile(), mode="wb")
Sleep_Disorders_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/SMQ_L.XPT", tf <- tempfile(), mode="wb")
Smoking_Cigarette_Use_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/SMQFAM_L.XPT", tf <- tempfile(), mode="wb")
Smoking_Household_Smokers_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/SMQRTU_L.XPT", tf <- tempfile(), mode="wb")
Smoking_Recent_Tobacco_Use_20211_2022 <- foreign::read.xport(tf)[,]

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/WHQ_L.XPT", tf <- tempfile(), mode="wb")
Weight_History_20211_2022 <- foreign::read.xport(tf)[,]
# -----------------

# Comorbidity per BMI group -----------------------

# Patientid
# Age and gender
# BMI (<25 ,25-27, 27-30, 30-35, 35-40, >40)
# HbA1c


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEMO_L.XPT", tf <- tempfile(), mode="wb")
Demographcics_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_20211_2022 %>% select(SEQN, RIDAGEYR, RIAGENDR) %>% filter(RIDAGEYR>=18) %>%
  mutate(RIAGENDR=ifelse(RIAGENDR==1,"M","F"))

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BMX_L.XPT", tf <- tempfile(), mode="wb")
Body_Measures_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Adults %>% left_join(Body_Measures_20211_2022 %>% select(SEQN, BMXBMI))

Adults <- Adults %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25",
                                           ifelse(BMXBMI<27,"25-27",
                                 ifelse(BMXBMI<30,"27-30",
                                        ifelse(BMXBMI<35,"30-35",
                                               ifelse(BMXBMI<40,"35-40",">40"))))))




download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/GHB_L.XPT", tf <- tempfile(), mode="wb")
Glycohemoglobin_20211_2022 <- foreign::read.xport(tf)[,]
Glycohemoglobin_20211_2022 <- Glycohemoglobin_20211_2022 %>% select(SEQN, LBXGH )
Adults <- Adults %>% left_join(Glycohemoglobin_20211_2022)



# Comorbs
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BPQ_L.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_Cholesterol_20211_2022 <- foreign::read.xport(tf)[,]
Blood_Pressure_Cholesterol_20211_2022 <- Blood_Pressure_Cholesterol_20211_2022 %>% select(SEQN,BPQ020,BPQ080) %>%
  rename("HTN"="BPQ020", "CHOL"="BPQ080")
Adults <- Adults %>% left_join(Blood_Pressure_Cholesterol_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/KIQ_U_L.XPT", tf <- tempfile(), mode="wb")
Kidney_Conditions_Urology_20211_2022 <- foreign::read.xport(tf)[,]
Kidney_Conditions_Urology_20211_2022 <- Kidney_Conditions_Urology_20211_2022 %>% select(SEQN,KIQ022) %>%
  rename("CKD"="KIQ022")
Adults <- Adults %>% left_join(Kidney_Conditions_Urology_20211_2022)

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DIQ_L.XPT", tf <- tempfile(), mode="wb")
Diabetes_20211_2022 <- foreign::read.xport(tf)[,]
Diabetes_20211_2022 <- Diabetes_20211_2022 %>% select(SEQN, DIQ010,DIQ160 ) %>%
     rename("T2D"="DIQ010", "PRET2D"="DIQ160")
Adults <- Adults %>% left_join(Diabetes_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/MCQ_L.XPT", tf <- tempfile(), mode="wb")
Medical_Conditions_20211_2022 <- foreign::read.xport(tf)[,]
Medical_Conditions_20211_2022 <- Medical_Conditions_20211_2022 %>% 
  select(SEQN, MCQ160A, MCQ160B,MCQ160C ,MCQ160E,MCQ160F ,MCQ160L ,MCQ510D ,MCQ510E ,MCQ510F ,MCQ220    ) %>%
     rename("ARTH"="MCQ160A", "CHF"="MCQ160B", "CAD"="MCQ160C", "HEARTATTACK"="MCQ160E", 
            "STROKE"="MCQ160F", "Liver"="MCQ160L", "ViralHepat"="MCQ510D", "AutoHepat"="MCQ510E",
            "OtherHepat"="MCQ510F", "Cancer"="MCQ220" )

Adults <- Adults %>% left_join(Medical_Conditions_20211_2022)


# BMI buckets
Adults %>% group_by(BMI_GRP) %>% count()
# HTN
Adults %>% group_by(BMI_GRP,HTN ) %>% count() %>% filter(HTN==1|HTN==2) %>% spread(key=HTN, value=n)
# CHOL
Adults %>% group_by(BMI_GRP,CHOL ) %>% count() %>% filter(CHOL==1|CHOL==2) %>% spread(key=CHOL, value=n)
# CKD
Adults %>% group_by(BMI_GRP,CKD ) %>% count() %>% filter(CKD==1|CKD==2) %>% spread(key=CKD, value=n)
# T2D
Adults %>% group_by(BMI_GRP,T2D ) %>% mutate(T2D=ifelse(T2D==3,1,T2D)) %>% count() %>% filter(T2D==1|T2D==2) %>% spread(key=T2D, value=n)
# Pre
Adults %>% group_by(BMI_GRP,PRET2D ) %>% count() %>% filter(PRET2D==1|PRET2D==2) %>% spread(key=PRET2D, value=n)
# Arth
Adults %>% group_by(BMI_GRP,ARTH ) %>% count() %>% filter(ARTH==1|ARTH==2) %>% spread(key=ARTH, value=n)
# CHF
Adults %>% group_by(BMI_GRP,CHF ) %>% count() %>% filter(CHF==1|CHF==2) %>% spread(key=CHF, value=n)
# CAD
Adults %>% group_by(BMI_GRP,CAD ) %>% count() %>% filter(CAD==1|CAD==2) %>% spread(key=CAD, value=n)
# HEARTATTACK
Adults %>% group_by(BMI_GRP,HEARTATTACK ) %>% count() %>% filter(HEARTATTACK==1|HEARTATTACK==2) %>% spread(key=HEARTATTACK, value=n)
# STROKE
Adults %>% group_by(BMI_GRP,STROKE ) %>% count() %>% filter(STROKE==1|STROKE==2) %>% spread(key=STROKE, value=n)
# Cancer
Adults %>% group_by(BMI_GRP,Cancer ) %>% count() %>% filter(Cancer==1|Cancer==2) %>% spread(key=Cancer, value=n)
# Liver
Adults %>% group_by(BMI_GRP,Liver ) %>% count() %>% filter(Liver==1|Liver==2) %>% spread(key=Liver, value=n)


Adults %>% filter(HEARTATTACK%in%c(1,2) & CAD%in%c(1,2) & STROKE%in%c(1,2) & CHF%in%c(1,2)) %>% # 7730
  filter(HEARTATTACK%in%c(1) | CAD%in%c(1) | STROKE%in%c(1) | CHF%in%c(1))

# Heart Failure, CAD, Stroke, PAD, CKD
# HYpertention, Cholestereol, Diabetes /Pre
# GERD, Sleep Apnea, NASH/NAFLD, RA, Osteoarthritis

# Cancer

CVD <- Adults %>% filter(HEARTATTACK%in%c(1) | CAD%in%c(1) | STROKE%in%c(1) | CHF%in%c(1) | CKD%in%c(1)) %>% select(SEQN)
RISK <- Adults %>% filter(HTN%in%c(1) | CHOL%in%c(1) | T2D%in%c(1) | PRET2D%in%c(1,3)) %>% select(SEQN)
OTHER <- Adults %>% filter(Liver%in%c(1) | ARTH%in%c(1) ) %>% select(SEQN)
CANCER <- Adults %>% filter(Cancer %in%c(1)) %>% select(SEQN)

length(unique(CVD$SEQN)) /  length(unique(Adults$SEQN))
length(unique(RISK$SEQN)) /  length(unique(Adults$SEQN))
length(unique(CANCER$SEQN)) /  length(unique(Adults$SEQN))
length(unique(OTHER$SEQN)) /  length(unique(Adults$SEQN))

RISK %>% inner_join(Adults) %>% group_by(BMI_GRP) %>% count()

# ---------------------
# Number of Comorbidities per BMI group -----------------------

# Patientid
# Age and gender
# BMI (<25 ,25-27, 27-30, 30-35, 35-40, >40)
# HbA1c


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEMO_L.XPT", tf <- tempfile(), mode="wb")
Demographcics_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_20211_2022 %>% select(SEQN, RIDAGEYR, RIAGENDR) %>% filter(RIDAGEYR>=18) %>%
  mutate(RIAGENDR=ifelse(RIAGENDR==1,"M","F"))

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BMX_L.XPT", tf <- tempfile(), mode="wb")
Body_Measures_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Adults %>% left_join(Body_Measures_20211_2022 %>% select(SEQN, BMXBMI))

Adults <- Adults %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25",
                                           ifelse(BMXBMI<27,"25-27",
                                 ifelse(BMXBMI<30,"27-30",
                                        ifelse(BMXBMI<35,"30-35",
                                               ifelse(BMXBMI<40,"35-40",">40"))))))




download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/GHB_L.XPT", tf <- tempfile(), mode="wb")
Glycohemoglobin_20211_2022 <- foreign::read.xport(tf)[,]
Glycohemoglobin_20211_2022 <- Glycohemoglobin_20211_2022 %>% select(SEQN, LBXGH )
Adults <- Adults %>% left_join(Glycohemoglobin_20211_2022)



# Comorbs
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BPQ_L.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_Cholesterol_20211_2022 <- foreign::read.xport(tf)[,]
Blood_Pressure_Cholesterol_20211_2022 <- Blood_Pressure_Cholesterol_20211_2022 %>% select(SEQN,BPQ020,BPQ080) %>%
  rename("HTN"="BPQ020", "CHOL"="BPQ080")
Adults <- Adults %>% left_join(Blood_Pressure_Cholesterol_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/KIQ_U_L.XPT", tf <- tempfile(), mode="wb")
Kidney_Conditions_Urology_20211_2022 <- foreign::read.xport(tf)[,]
Kidney_Conditions_Urology_20211_2022 <- Kidney_Conditions_Urology_20211_2022 %>% select(SEQN,KIQ022) %>%
  rename("CKD"="KIQ022")
Adults <- Adults %>% left_join(Kidney_Conditions_Urology_20211_2022)

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DIQ_L.XPT", tf <- tempfile(), mode="wb")
Diabetes_20211_2022 <- foreign::read.xport(tf)[,]
Diabetes_20211_2022 <- Diabetes_20211_2022 %>% select(SEQN, DIQ010,DIQ160 ) %>%
     rename("T2D"="DIQ010", "PRET2D"="DIQ160")
Adults <- Adults %>% left_join(Diabetes_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/MCQ_L.XPT", tf <- tempfile(), mode="wb")
Medical_Conditions_20211_2022 <- foreign::read.xport(tf)[,]
Medical_Conditions_20211_2022 <- Medical_Conditions_20211_2022 %>% 
  select(SEQN, MCQ160A, MCQ160B,MCQ160C ,MCQ160E,MCQ160F ,MCQ160L ,MCQ510D ,MCQ510E ,MCQ510F ,MCQ220    ) %>%
     rename("ARTH"="MCQ160A", "CHF"="MCQ160B", "CAD"="MCQ160C", "HEARTATTACK"="MCQ160E", 
            "STROKE"="MCQ160F", "Liver"="MCQ160L", "ViralHepat"="MCQ510D", "AutoHepat"="MCQ510E",
            "OtherHepat"="MCQ510F", "Cancer"="MCQ220" )

Adults <- Adults %>% left_join(Medical_Conditions_20211_2022)


# BMI buckets
Adults %>% group_by(BMI_GRP) %>% count()
# HTN
Adults %>% group_by(BMI_GRP,HTN ) %>% count() %>% filter(HTN==1|HTN==2) %>% spread(key=HTN, value=n)
# CHOL
Adults %>% group_by(BMI_GRP,CHOL ) %>% count() %>% filter(CHOL==1|CHOL==2) %>% spread(key=CHOL, value=n)
# CKD
Adults %>% group_by(BMI_GRP,CKD ) %>% count() %>% filter(CKD==1|CKD==2) %>% spread(key=CKD, value=n)
# T2D
Adults %>% group_by(BMI_GRP,T2D ) %>% mutate(T2D=ifelse(T2D==3,1,T2D)) %>% count() %>% filter(T2D==1|T2D==2) %>% spread(key=T2D, value=n)
# Pre
Adults %>% group_by(BMI_GRP,PRET2D ) %>% count() %>% filter(PRET2D==1|PRET2D==2) %>% spread(key=PRET2D, value=n)
# Arth
Adults %>% group_by(BMI_GRP,ARTH ) %>% count() %>% filter(ARTH==1|ARTH==2) %>% spread(key=ARTH, value=n)
# CHF
Adults %>% group_by(BMI_GRP,CHF ) %>% count() %>% filter(CHF==1|CHF==2) %>% spread(key=CHF, value=n)
# CAD
Adults %>% group_by(BMI_GRP,CAD ) %>% count() %>% filter(CAD==1|CAD==2) %>% spread(key=CAD, value=n)
# HEARTATTACK
Adults %>% group_by(BMI_GRP,HEARTATTACK ) %>% count() %>% filter(HEARTATTACK==1|HEARTATTACK==2) %>% spread(key=HEARTATTACK, value=n)
# STROKE
Adults %>% group_by(BMI_GRP,STROKE ) %>% count() %>% filter(STROKE==1|STROKE==2) %>% spread(key=STROKE, value=n)
# Cancer
Adults %>% group_by(BMI_GRP,Cancer ) %>% count() %>% filter(Cancer==1|Cancer==2) %>% spread(key=Cancer, value=n)
# Liver
Adults %>% group_by(BMI_GRP,Liver ) %>% count() %>% filter(Liver==1|Liver==2) %>% spread(key=Liver, value=n)


Adults %>% filter(HEARTATTACK%in%c(1,2) & CAD%in%c(1,2) & STROKE%in%c(1,2) & CHF%in%c(1,2)) %>% # 7730
  filter(HEARTATTACK%in%c(1) | CAD%in%c(1) | STROKE%in%c(1) | CHF%in%c(1))

# Heart Failure, CAD, Stroke, PAD, CKD
# HYpertention, Cholestereol, Diabetes /Pre
# GERD, Sleep Apnea, NASH/NAFLD, RA, Osteoarthritis

# Cancer

CVD <- Adults %>% filter(HEARTATTACK%in%c(1) | CAD%in%c(1) | STROKE%in%c(1) | CHF%in%c(1) | CKD%in%c(1)) %>% select(SEQN)
RISK <- Adults %>% filter(HTN%in%c(1) | CHOL%in%c(1) | T2D%in%c(1) | PRET2D%in%c(1,3)) %>% select(SEQN)
OTHER <- Adults %>% filter(Liver%in%c(1) | ARTH%in%c(1) ) %>% select(SEQN)
CANCER <- Adults %>% filter(Cancer %in%c(1)) %>% select(SEQN)

length(unique(CVD$SEQN)) /  length(unique(Adults$SEQN))
length(unique(RISK$SEQN)) /  length(unique(Adults$SEQN))
length(unique(CANCER$SEQN)) /  length(unique(Adults$SEQN))
length(unique(OTHER$SEQN)) /  length(unique(Adults$SEQN))

RISK %>% inner_join(Adults) %>% group_by(BMI_GRP) %>% count()


Adults <- Adults %>% select(-c(RIDAGEYR, RIAGENDR, LBXGH, ViralHepat, AutoHepat, OtherHepat, Cancer))


Adults <- Adults %>% 
  mutate(HTN=ifelse(HTN==1,1,0)) %>%
  mutate(CHOL=ifelse(CHOL==1,1,0)) %>%
  mutate(CKD =ifelse(CKD ==1,1,0)) %>%
  mutate(T2D =ifelse(T2D==1|T2D==3,1,0)) %>%
  mutate(PRET2D =ifelse(PRET2D ==1,1,0)) %>%
  mutate(ARTH =ifelse(ARTH ==1,1,0)) %>%
  mutate(CHF =ifelse(CHF ==1,1,0)) %>%
  mutate(CAD =ifelse(CAD ==1,1,0)) %>%
  mutate(HEARTATTACK  =ifelse(HEARTATTACK  ==1,1,0)) %>%
  mutate(STROKE  =ifelse(STROKE  ==1,1,0)) %>%
  mutate(Liver =ifelse(Liver ==1,1,0)) 

Adults <- Adults %>% filter(!is.na(BMI_GRP))


Adults[is.na(Adults)] <- 0


Adults <- Adults %>% mutate(Total=HTN+CHOL+CKD+T2D+PRET2D+ARTH+CHF+CAD+HEARTATTACK+STROKE+Liver)

Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  group_by(GROUP) %>% count()

Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  group_by(BMI_GRP, GROUP) %>% count() %>%
  spread(key=BMI_GRP, value=n)



# ---------------------
# BMI breakdown by AGe and Race -----------



download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEMO_L.XPT", tf <- tempfile(), mode="wb")
Demographcics_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_20211_2022 %>% select(SEQN, RIDAGEYR, RIAGENDR,RIDRETH3 ,RIDRETH1) %>% filter(RIDAGEYR>=18) %>%
  mutate(RIAGENDR=ifelse(RIAGENDR==1,"M","F"))



download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BMX_L.XPT", tf <- tempfile(), mode="wb")
Body_Measures_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Adults %>% left_join(Body_Measures_20211_2022 %>% select(SEQN, BMXBMI))

Adults <- Adults %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25",
                                           ifelse(BMXBMI<27,"25-27",
                                 ifelse(BMXBMI<30,"27-30",
                                        ifelse(BMXBMI<35,"30-35",
                                               ifelse(BMXBMI<40,"35-40",">40"))))))


Adults <- Adults %>% filter(!is.na(BMI_GRP)) %>% filter(RIDAGEYR>=18)
Adults <- Adults %>% select(-BMXBMI)

unique(Adults$RIDRETH1)
unique(Adults$RIDRETH3)


range(Adults$RIDAGEYR)


Adults <- Adults %>% mutate(RIDAGEYR=ifelse(RIDAGEYR<=30, "<30",
                                  ifelse(RIDAGEYR<=40,"<40",
                                         ifelse(RIDAGEYR<=50,"<50",
                                                ifelse(RIDAGEYR<=60,"<60",
                                                       ifelse(RIDAGEYR<=70,"<70",
                                                              ifelse(RIDAGEYR<=80,"80", NA)))))))


Adults %>% group_by(RIDAGEYR, BMI_GRP) %>% count() %>%
  spread(key=BMI_GRP, value=n)


Adults %>% group_by(RIDRETH3, BMI_GRP) %>% count() %>%
  spread(key=BMI_GRP, value=n)

# -----------
# BMIs in 18-19 y/o ------


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEMO_L.XPT", tf <- tempfile(), mode="wb")
Demographcics_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_20211_2022 %>% select(SEQN, RIDAGEYR) %>% filter(RIDAGEYR>=18&RIDAGEYR<=19) 


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BMX_L.XPT", tf <- tempfile(), mode="wb")
Body_Measures_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Adults %>% left_join(Body_Measures_20211_2022 %>% select(SEQN, BMXBMI))

Adults <- Adults %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25", ifelse(BMXBMI<30,"25-30", ">30")))

Adults <- Adults %>% filter(!is.na(BMI_GRP))


Adults %>% group_by(BMI_GRP) %>% count() %>% mutate(n=n/265)


# ----

# 2017- 2018 data  Check Sleep Apnea-------------
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/SLQ_J.XPT", tf <- tempfile(), mode="wb")
Sleep_20217_2018 <- foreign::read.xport(tf)[,]

Sleep_20217_2018 <- Sleep_20217_2018 %>% select(SEQN, SLQ030, SLQ040)



download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BMX_J.XPT", tf <- tempfile(), mode="wb")
BMI_20217_2018 <- foreign::read.xport(tf)[,]
BMI_20217_2018 <- BMI_20217_2018 %>% select(SEQN, BMXBMI) 
BMI_20217_2018 <- BMI_20217_2018 %>% filter(!is.na(BMXBMI))


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT", tf <- tempfile(), mode="wb")
Demographcics_20217_2018 <- foreign::read.xport(tf)[,]
Demographcics_20217_2018 <- Demographcics_20217_2018 %>% select(SEQN, RIAGENDR , RIDAGEYR) %>% filter(RIDAGEYR>=18) 


Sleep_20217_2018 <- Demographcics_20217_2018 %>% inner_join(BMI_20217_2018) %>% left_join(Sleep_20217_2018)

Sleep_20217_2018 <- Sleep_20217_2018 %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25",
                                           ifelse(BMXBMI<27,"25-27",
                                 ifelse(BMXBMI<30,"27-30",
                                        ifelse(BMXBMI<35,"30-35",
                                               ifelse(BMXBMI<40,"35-40",">40"))))))


unique(Sleep_20217_2018$SLQ040)
Sleep_20217_2018 <- Sleep_20217_2018 %>% mutate(OSA=ifelse(SLQ040>=3,1,0))

Sleep_20217_2018 %>% group_by(BMI_GRP, OSA) %>% count() %>% spread(key=OSA, value=n) %>%
  mutate(perc=`1`/(`1`+`0`)) %>% arrange(-perc)


Sleep_20217_2018 %>% group_by(OSA) %>% count()





# Comorbs
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPQ_J.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_Cholesterol_2017_2018 <- foreign::read.xport(tf)[,]
Blood_Pressure_Cholesterol_2017_2018 <- Blood_Pressure_Cholesterol_2017_2018 %>% select(SEQN,BPQ020,BPQ080) %>%
  rename("HTN"="BPQ020", "CHOL"="BPQ080")


Sleep_20217_2018 <- Sleep_20217_2018 %>% left_join(Blood_Pressure_Cholesterol_2017_2018)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/KIQ_U_J.XPT", tf <- tempfile(), mode="wb")
Kidney_Conditions_Urology_2017_2018 <- foreign::read.xport(tf)[,]
Kidney_Conditions_Urology_2017_2018 <- Kidney_Conditions_Urology_2017_2018 %>% select(SEQN,KIQ022) %>%
  rename("CKD"="KIQ022")
Sleep_20217_2018 <- Sleep_20217_2018 %>% left_join(Kidney_Conditions_Urology_2017_2018)

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DIQ_J.XPT", tf <- tempfile(), mode="wb")
Diabetes_2017_2018 <- foreign::read.xport(tf)[,]
Diabetes_2017_2018 <- Diabetes_2017_2018 %>% select(SEQN, DIQ010,DIQ160 ) %>%
     rename("T2D"="DIQ010", "PRET2D"="DIQ160")
Sleep_20217_2018 <- Sleep_20217_2018 %>% left_join(Diabetes_2017_2018)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/MCQ_J.XPT", tf <- tempfile(), mode="wb")
Medical_Conditions_2017_2018 <- foreign::read.xport(tf)[,]
Medical_Conditions_2017_2018 <- Medical_Conditions_2017_2018 %>% 
  select(SEQN, MCQ160A, MCQ160B,MCQ160C ,MCQ160E,MCQ160F ,MCQ160L ,MCQ220    ) %>%
     rename("ARTH"="MCQ160A", "CHF"="MCQ160B", "CAD"="MCQ160C", "HEARTATTACK"="MCQ160E", 
            "STROKE"="MCQ160F", "Liver"="MCQ160L",  "Cancer"="MCQ220" )

Sleep_20217_2018 <- Sleep_20217_2018 %>% left_join(Medical_Conditions_2017_2018)



Sleep_20217_2018 <- Sleep_20217_2018 %>% 
  mutate(HTN=ifelse(HTN==1,1,0)) %>%
  mutate(CHOL=ifelse(CHOL==1,1,0)) %>%
  mutate(CKD =ifelse(CKD ==1,1,0)) %>%
  mutate(T2D =ifelse(T2D==1|T2D==3,1,0)) %>%
  mutate(PRET2D =ifelse(PRET2D ==1,1,0)) %>%
  mutate(ARTH =ifelse(ARTH ==1,1,0)) %>%
  mutate(CHF =ifelse(CHF ==1,1,0)) %>%
  mutate(CAD =ifelse(CAD ==1,1,0)) %>%
  mutate(HEARTATTACK  =ifelse(HEARTATTACK  ==1,1,0)) %>%
  mutate(STROKE  =ifelse(STROKE  ==1,1,0)) %>%
  mutate(Liver =ifelse(Liver ==1,1,0))  %>%
  mutate(Cancer =ifelse(Cancer ==1,1,0))  

Sleep_20217_2018[is.na(Sleep_20217_2018)] <- 0


temp <- data.frame(Sleep_20217_2018 %>% group_by(BMI_GRP, HTN, CHOL, CKD, T2D, PRET2D, ARTH, CHF, CAD, HEARTATTACK, STROKE, Liver, OSA) %>% 
  count() %>%
  spread(key=OSA, value=n) %>% mutate(OSA=`1`/(`1`+`0`)) %>%
  select(-c(`1`, `0`)))





Sleep_20217_2018 <- Sleep_20217_2018 %>% mutate(RIDAGEYR=ifelse(RIDAGEYR<=30, "<30",
                                  ifelse(RIDAGEYR<=40,"<40",
                                         ifelse(RIDAGEYR<=50,"<50",
                                                ifelse(RIDAGEYR<=60,"<60",
                                                       ifelse(RIDAGEYR<=70,"<70",
                                                              ifelse(RIDAGEYR<=80,"80", NA)))))))


Sleep_20217_2018 %>% group_by(RIDAGEYR, BMI_GRP, OSA) %>%  count()  %>%
  spread(key=OSA, value=n) %>% mutate(OSA=100*`1`/(`1`+`0`)) %>%
  select(-c(`1`, `0`)) %>%
  spread(key=BMI_GRP,value=OSA) %>% 
  select(RIDAGEYR  , `<25`, `25-27` ,`27-30` ,`30-35`, `35-40`, `>40`)


#   RIDAGEYR `<25` `25-27` `27-30` `30-35` `35-40` `>40`
# 1 <30       2.56    4.17    4.58    8.25    7.62  6.98
# 2 <40       3.03    9.09    2.08   11.1    11.8  16.5 
# 3 <50       4.14    9.57   10.7    16.1    17.2  21.6 
# 4 <60       8.37    6.36   13.1    12.3    14.3  23   
# 5 <70       9.31   10.3     9.84   14.0    18.0  26.7 
# 6 80       13.2     9.23   16.7    14.5    16.7  28.6 
# -----------------


# Number of Comorbidities per BMI group with OSA -----------------------


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEMO_L.XPT", tf <- tempfile(), mode="wb")
Demographcics_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_20211_2022 %>% select(SEQN, RIDAGEYR, RIAGENDR) %>% filter(RIDAGEYR>=18) %>%
  mutate(RIAGENDR=ifelse(RIAGENDR==1,"M","F"))

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BMX_L.XPT", tf <- tempfile(), mode="wb")
Body_Measures_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Adults %>% left_join(Body_Measures_20211_2022 %>% select(SEQN, BMXBMI))

Adults <- Adults %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25",
                                           ifelse(BMXBMI<27,"25-27",
                                 ifelse(BMXBMI<30,"27-30",
                                        ifelse(BMXBMI<35,"30-35",
                                               ifelse(BMXBMI<40,"35-40",">40"))))))



Adults <- Adults %>% mutate(RIDAGEYR=ifelse(RIDAGEYR<=30, "<30",
                                  ifelse(RIDAGEYR<=40,"<40",
                                         ifelse(RIDAGEYR<=50,"<50",
                                                ifelse(RIDAGEYR<=60,"<60",
                                                       ifelse(RIDAGEYR<=70,"<70",
                                                              ifelse(RIDAGEYR<=80,"80", NA)))))))



# Comorbs
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BPQ_L.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_Cholesterol_20211_2022 <- foreign::read.xport(tf)[,]
Blood_Pressure_Cholesterol_20211_2022 <- Blood_Pressure_Cholesterol_20211_2022 %>% select(SEQN,BPQ020,BPQ080) %>%
  rename("HTN"="BPQ020", "CHOL"="BPQ080")
Adults <- Adults %>% left_join(Blood_Pressure_Cholesterol_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/KIQ_U_L.XPT", tf <- tempfile(), mode="wb")
Kidney_Conditions_Urology_20211_2022 <- foreign::read.xport(tf)[,]
Kidney_Conditions_Urology_20211_2022 <- Kidney_Conditions_Urology_20211_2022 %>% select(SEQN,KIQ022) %>%
  rename("CKD"="KIQ022")
Adults <- Adults %>% left_join(Kidney_Conditions_Urology_20211_2022)

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DIQ_L.XPT", tf <- tempfile(), mode="wb")
Diabetes_20211_2022 <- foreign::read.xport(tf)[,]
Diabetes_20211_2022 <- Diabetes_20211_2022 %>% select(SEQN, DIQ010,DIQ160 ) %>%
     rename("T2D"="DIQ010", "PRET2D"="DIQ160")
Adults <- Adults %>% left_join(Diabetes_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/MCQ_L.XPT", tf <- tempfile(), mode="wb")
Medical_Conditions_20211_2022 <- foreign::read.xport(tf)[,]
Medical_Conditions_20211_2022 <- Medical_Conditions_20211_2022 %>% 
  select(SEQN, MCQ160A, MCQ160B,MCQ160C ,MCQ160E,MCQ160F ,MCQ160L ,MCQ510D ,MCQ510E ,MCQ510F ,MCQ220    ) %>%
     rename("ARTH"="MCQ160A", "CHF"="MCQ160B", "CAD"="MCQ160C", "HEARTATTACK"="MCQ160E", 
            "STROKE"="MCQ160F", "Liver"="MCQ160L", "ViralHepat"="MCQ510D", "AutoHepat"="MCQ510E",
            "OtherHepat"="MCQ510F", "Cancer"="MCQ220" )

Adults <- Adults %>% left_join(Medical_Conditions_20211_2022)



Adults <- Adults %>% 
  mutate(HTN=ifelse(HTN==1,1,0)) %>%
  mutate(CHOL=ifelse(CHOL==1,1,0)) %>%
  mutate(CKD =ifelse(CKD ==1,1,0)) %>%
  mutate(T2D =ifelse(T2D==1|T2D==3,1,0)) %>%
  mutate(PRET2D =ifelse(PRET2D ==1,1,0)) %>%
  mutate(ARTH =ifelse(ARTH ==1,1,0)) %>%
  mutate(CHF =ifelse(CHF ==1,1,0)) %>%
  mutate(CAD =ifelse(CAD ==1,1,0)) %>%
  mutate(HEARTATTACK  =ifelse(HEARTATTACK  ==1,1,0)) %>%
  mutate(STROKE  =ifelse(STROKE  ==1,1,0)) %>%
  mutate(Liver =ifelse(Liver ==1,1,0)) 

Adults <- Adults %>% filter(!is.na(BMI_GRP))


names(Adults)


Adults <- Adults %>% select(-c(ViralHepat, AutoHepat, OtherHepat, Cancer))

Adults[is.na(Adults)] <- 0

# 
#   RIDAGEYR `<25` `25-27` `27-30` `30-35` `35-40` `>40`
# 1 <30       2.56    4.17    4.58    8.25    7.62  6.98
# 2 <40       3.03    9.09    2.08   11.1    11.8  16.5
# 3 <50       4.14    9.57   10.7    16.1    17.2  21.6
# 4 <60       8.37    6.36   13.1    12.3    14.3  23
# 5 <70       9.31   10.3     9.84   14.0    18.0  26.7
# 6 80       13.2     9.23   16.7    14.5    16.7  28.6



unique(Adults$RIDAGEYR)
unique(Adults$BMI_GRP)


Adults <- Adults %>%
  mutate(
    percent_ones = case_when(
      RIDAGEYR == "<30" & BMI_GRP  == "<25" ~ 2.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "<25" ~ 3.0,  
      RIDAGEYR == "<50" & BMI_GRP  == "<25" ~ 4.1,  
      RIDAGEYR == "<60" & BMI_GRP  == "<25" ~ 8.4,  
      RIDAGEYR == "<70" & BMI_GRP  == "<25" ~ 9.3,  
      RIDAGEYR == "80" & BMI_GRP  == "<25" ~ 13.2,  
       RIDAGEYR == "<30" & BMI_GRP  == "25-27" ~ 4.2,  
      RIDAGEYR == "<40" & BMI_GRP  == "25-27" ~ 9.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "25-27" ~ 9.6,  
      RIDAGEYR == "<60" & BMI_GRP  == "25-27" ~ 6.4,  
      RIDAGEYR == "<70" & BMI_GRP  == "25-27" ~ 10.3,  
      RIDAGEYR == "80" & BMI_GRP  == "25-27" ~ 9.2,  
       RIDAGEYR == "<30" & BMI_GRP  == "27-30" ~ 4.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "27-30" ~ 2.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "27-30" ~ 10.7,  
      RIDAGEYR == "<60" & BMI_GRP  == "27-30" ~ 13.1,  
      RIDAGEYR == "<70" & BMI_GRP  == "27-30" ~ 9.8,  
      RIDAGEYR == "80" & BMI_GRP  == "27-30" ~ 16.7,  
       RIDAGEYR == "<30" & BMI_GRP  == "30-35" ~ 8.3,  
      RIDAGEYR == "<40" & BMI_GRP  == "30-35" ~ 11.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "30-35" ~ 16.1,  
      RIDAGEYR == "<60" & BMI_GRP  == "30-35" ~ 12.3,  
      RIDAGEYR == "<70" & BMI_GRP  == "30-35" ~ 14.0,  
      RIDAGEYR == "80" & BMI_GRP  == "30-35" ~ 14.5,  
      RIDAGEYR == "<30" & BMI_GRP  == "35-40" ~ 7.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "35-40" ~ 11.8,  
      RIDAGEYR == "<50" & BMI_GRP  == "35-40" ~ 17.2,  
      RIDAGEYR == "<60" & BMI_GRP  == "35-40" ~ 14.3,  
      RIDAGEYR == "<70" & BMI_GRP  == "35-40" ~ 18.0,  
      RIDAGEYR == "80" & BMI_GRP  == "35-40" ~ 16.7,  
        RIDAGEYR == "<30" & BMI_GRP  == ">40" ~ 7.0,  
      RIDAGEYR == "<40" & BMI_GRP  == ">40" ~ 16.5,  
      RIDAGEYR == "<50" & BMI_GRP  == ">40" ~ 21.6,  
      RIDAGEYR == "<60" & BMI_GRP  == ">40" ~ 23,  
      RIDAGEYR == "<70" & BMI_GRP  == ">40" ~ 26.7,  
      RIDAGEYR == "80" & BMI_GRP  == ">40" ~ 14
    ))


temp <- data.frame(Adults %>% arrange(BMI_GRP, RIDAGEYR, SEQN) %>%
  group_by(RIDAGEYR, BMI_GRP) %>% mutate(pop=n()) %>%
  select(SEQN, RIDAGEYR, BMI_GRP, percent_ones, pop) %>%
  mutate(lines= round(151*percent_ones/100) ) %>%
  mutate(line=row_number()))


temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>% group_by(RIDAGEYR, BMI_GRP, OSA) %>% count() %>%
   spread(key=OSA, value=n) %>% mutate(OSA=100*`1`/(`1`+`0`)) %>%
  select(-c(`1`, `0`)) %>%
  spread(key=BMI_GRP,value=OSA) %>% 
  select(RIDAGEYR  , `<25`, `25-27` ,`27-30` ,`30-35`, `35-40`, `>40`)



temp %>% mutate(OSA=ifelse(line<=lines,1,0))  %>% group_by(OSA) %>% count()


temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>%
  inner_join(Adults %>% select(SEQN, Liver)) %>%
  group_by(Liver, OSA) %>% count() 


Adults <- temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>% inner_join(Adults) %>%
  select(-c(pop, percent_ones, lines, line))

Adults <- Adults %>% mutate(Total=HTN+CHOL+CKD+T2D+PRET2D+ARTH+CHF+CAD+HEARTATTACK+STROKE+Liver+OSA)

Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  group_by(GROUP) %>% count()

Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  group_by(BMI_GRP, GROUP) %>% count() %>%
  spread(key=BMI_GRP, value=n)


unique(Adults$BMI_GRP)

names(Adults)

Adults %>% filter(BMI_GRP!="25-27"&BMI_GRP!="<25") %>%
  summarise(HTN=mean(OSA))


# ---------


# Income be Number of Comorbidities per BMI group 2017 2018 -----------------------


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DEMO_J.XPT", tf <- tempfile(), mode="wb")
Demographcics_2017_2018 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_2017_2018 %>% select(SEQN, RIDAGEYR, RIAGENDR) %>% filter(RIDAGEYR>=18) %>%
  mutate(RIAGENDR=ifelse(RIAGENDR==1,"M","F"))


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BMX_J.XPT", tf <- tempfile(), mode="wb")
BMI_20217_2018 <- foreign::read.xport(tf)[,]
BMI_20217_2018 <- BMI_20217_2018 %>% select(SEQN, BMXBMI) 
BMI_20217_2018 <- BMI_20217_2018 %>% filter(!is.na(BMXBMI))

Adults <- Adults %>% inner_join(BMI_20217_2018)

Adults <- Adults %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25",
                                           ifelse(BMXBMI<27,"25-27",
                                 ifelse(BMXBMI<30,"27-30",
                                        ifelse(BMXBMI<35,"30-35",
                                               ifelse(BMXBMI<40,"35-40",">40"))))))



# Comorbs
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/BPQ_J.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_Cholesterol_2017_2018 <- foreign::read.xport(tf)[,]
Blood_Pressure_Cholesterol_2017_2018 <- Blood_Pressure_Cholesterol_2017_2018 %>% select(SEQN,BPQ020,BPQ080) %>%
  rename("HTN"="BPQ020", "CHOL"="BPQ080")
Adults <- Adults %>% left_join(Blood_Pressure_Cholesterol_2017_2018)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/KIQ_U_J.XPT", tf <- tempfile(), mode="wb")
Kidney_Conditions_Urology_2017_2018 <- foreign::read.xport(tf)[,]
Kidney_Conditions_Urology_2017_2018 <- Kidney_Conditions_Urology_2017_2018 %>% select(SEQN,KIQ022) %>%
  rename("CKD"="KIQ022")
Adults <- Adults %>% left_join(Kidney_Conditions_Urology_2017_2018)

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/DIQ_J.XPT", tf <- tempfile(), mode="wb")
Diabetes_2017_2018 <- foreign::read.xport(tf)[,]
Diabetes_2017_2018 <- Diabetes_2017_2018 %>% select(SEQN, DIQ010,DIQ160 ) %>%
     rename("T2D"="DIQ010", "PRET2D"="DIQ160")
Adults <- Adults %>% left_join(Diabetes_2017_2018)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/MCQ_J.XPT", tf <- tempfile(), mode="wb")
Medical_Conditions_2017_2018 <- foreign::read.xport(tf)[,]
Medical_Conditions_2017_2018 <- Medical_Conditions_2017_2018 %>% 
  select(SEQN, MCQ160A, MCQ160B,MCQ160C ,MCQ160E,MCQ160F ,MCQ160L ,MCQ510D ,MCQ510E ,MCQ510F ,MCQ220    ) %>%
     rename("ARTH"="MCQ160A", "CHF"="MCQ160B", "CAD"="MCQ160C", "HEARTATTACK"="MCQ160E", 
            "STROKE"="MCQ160F", "Liver"="MCQ160L", "ViralHepat"="MCQ510D", "AutoHepat"="MCQ510E",
            "OtherHepat"="MCQ510F", "Cancer"="MCQ220" )

Adults <- Adults %>% left_join(Medical_Conditions_2017_2018)



Adults <- Adults %>% 
  mutate(HTN=ifelse(HTN==1,1,0)) %>%
  mutate(CHOL=ifelse(CHOL==1,1,0)) %>%
  mutate(CKD =ifelse(CKD ==1,1,0)) %>%
  mutate(T2D =ifelse(T2D==1|T2D==3,1,0)) %>%
  mutate(PRET2D =ifelse(PRET2D ==1,1,0)) %>%
  mutate(ARTH =ifelse(ARTH ==1,1,0)) %>%
  mutate(CHF =ifelse(CHF ==1,1,0)) %>%
  mutate(CAD =ifelse(CAD ==1,1,0)) %>%
  mutate(HEARTATTACK  =ifelse(HEARTATTACK  ==1,1,0)) %>%
  mutate(STROKE  =ifelse(STROKE  ==1,1,0)) %>%
  mutate(Liver =ifelse(Liver ==1,1,0)) 

Adults <- Adults %>% filter(!is.na(BMI_GRP))


names(Adults)

Adults <- Adults %>% select(-c(ViralHepat, AutoHepat, OtherHepat, Cancer))

Adults[is.na(Adults)] <- 0

Adults <- Adults %>% mutate(Total=HTN+CHOL+CKD+T2D+PRET2D+ARTH+CHF+CAD+HEARTATTACK+STROKE+Liver)

Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  group_by(GROUP) %>% count()

Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  group_by(BMI_GRP, GROUP) %>% count() %>%
  spread(key=BMI_GRP, value=n)



download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2017-2018/INQ_J.XPT", tf <- tempfile(), mode="wb")
Income_2017_2018 <- foreign::read.xport(tf)[,]

Income_2017_2018 <- Income_2017_2018 %>% select(SEQN, IND235)
unique(Income_2017_2018$IND235)

Income_2017_2018 <- Income_2017_2018 %>% filter(IND235!=77&IND235!=99)


data.frame(Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  select(SEQN, BMI_GRP, GROUP) %>% 
  group_by(SEQN) %>%
  mutate(group2=paste(BMI_GRP, GROUP)) %>%
  inner_join(Income_2017_2018) %>%
  group_by(group2, IND235) %>% count() %>%
  spread(key=IND235, value=n))



# ---------


# Insurances --------

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEMO_L.XPT", tf <- tempfile(), mode="wb")
Demographcics_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_20211_2022 %>% select(SEQN, RIDAGEYR, RIAGENDR) %>% filter(RIDAGEYR>=18) %>%
  mutate(RIAGENDR=ifelse(RIAGENDR==1,"M","F"))

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BMX_L.XPT", tf <- tempfile(), mode="wb")
Body_Measures_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Adults %>% left_join(Body_Measures_20211_2022 %>% select(SEQN, BMXBMI))

Adults <- Adults %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25",
                                           ifelse(BMXBMI<27,"25-27",
                                 ifelse(BMXBMI<30,"27-30",
                                        ifelse(BMXBMI<35,"30-35",
                                               ifelse(BMXBMI<40,"35-40",">40"))))))



Adults <- Adults %>% mutate(RIDAGEYR=ifelse(RIDAGEYR<=30, "<30",
                                  ifelse(RIDAGEYR<=40,"<40",
                                         ifelse(RIDAGEYR<=50,"<50",
                                                ifelse(RIDAGEYR<=60,"<60",
                                                       ifelse(RIDAGEYR<=70,"<70",
                                                              ifelse(RIDAGEYR<=80,"80", NA)))))))



# Comorbs
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BPQ_L.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_Cholesterol_20211_2022 <- foreign::read.xport(tf)[,]
Blood_Pressure_Cholesterol_20211_2022 <- Blood_Pressure_Cholesterol_20211_2022 %>% select(SEQN,BPQ020,BPQ080) %>%
  rename("HTN"="BPQ020", "CHOL"="BPQ080")
Adults <- Adults %>% left_join(Blood_Pressure_Cholesterol_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/KIQ_U_L.XPT", tf <- tempfile(), mode="wb")
Kidney_Conditions_Urology_20211_2022 <- foreign::read.xport(tf)[,]
Kidney_Conditions_Urology_20211_2022 <- Kidney_Conditions_Urology_20211_2022 %>% select(SEQN,KIQ022) %>%
  rename("CKD"="KIQ022")
Adults <- Adults %>% left_join(Kidney_Conditions_Urology_20211_2022)

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DIQ_L.XPT", tf <- tempfile(), mode="wb")
Diabetes_20211_2022 <- foreign::read.xport(tf)[,]
Diabetes_20211_2022 <- Diabetes_20211_2022 %>% select(SEQN, DIQ010,DIQ160 ) %>%
     rename("T2D"="DIQ010", "PRET2D"="DIQ160")
Adults <- Adults %>% left_join(Diabetes_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/MCQ_L.XPT", tf <- tempfile(), mode="wb")
Medical_Conditions_20211_2022 <- foreign::read.xport(tf)[,]
Medical_Conditions_20211_2022 <- Medical_Conditions_20211_2022 %>% 
  select(SEQN, MCQ160A, MCQ160B,MCQ160C ,MCQ160E,MCQ160F ,MCQ160L ,MCQ510D ,MCQ510E ,MCQ510F ,MCQ220    ) %>%
     rename("ARTH"="MCQ160A", "CHF"="MCQ160B", "CAD"="MCQ160C", "HEARTATTACK"="MCQ160E", 
            "STROKE"="MCQ160F", "Liver"="MCQ160L", "ViralHepat"="MCQ510D", "AutoHepat"="MCQ510E",
            "OtherHepat"="MCQ510F", "Cancer"="MCQ220" )

Adults <- Adults %>% left_join(Medical_Conditions_20211_2022)



Adults <- Adults %>% 
  mutate(HTN=ifelse(HTN==1,1,0)) %>%
  mutate(CHOL=ifelse(CHOL==1,1,0)) %>%
  mutate(CKD =ifelse(CKD ==1,1,0)) %>%
  mutate(T2D =ifelse(T2D==1|T2D==3,1,0)) %>%
  mutate(PRET2D =ifelse(PRET2D ==1,1,0)) %>%
  mutate(ARTH =ifelse(ARTH ==1,1,0)) %>%
  mutate(CHF =ifelse(CHF ==1,1,0)) %>%
  mutate(CAD =ifelse(CAD ==1,1,0)) %>%
  mutate(HEARTATTACK  =ifelse(HEARTATTACK  ==1,1,0)) %>%
  mutate(STROKE  =ifelse(STROKE  ==1,1,0)) %>%
  mutate(Liver =ifelse(Liver ==1,1,0)) 

Adults <- Adults %>% filter(!is.na(BMI_GRP))


names(Adults)


Adults <- Adults %>% select(-c(ViralHepat, AutoHepat, OtherHepat, Cancer))

Adults[is.na(Adults)] <- 0

# 
#   RIDAGEYR `<25` `25-27` `27-30` `30-35` `35-40` `>40`
# 1 <30       2.56    4.17    4.58    8.25    7.62  6.98
# 2 <40       3.03    9.09    2.08   11.1    11.8  16.5
# 3 <50       4.14    9.57   10.7    16.1    17.2  21.6
# 4 <60       8.37    6.36   13.1    12.3    14.3  23
# 5 <70       9.31   10.3     9.84   14.0    18.0  26.7
# 6 80       13.2     9.23   16.7    14.5    16.7  28.6



unique(Adults$RIDAGEYR)
unique(Adults$BMI_GRP)


Adults <- Adults %>%
  mutate(
    percent_ones = case_when(
      RIDAGEYR == "<30" & BMI_GRP  == "<25" ~ 2.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "<25" ~ 3.0,  
      RIDAGEYR == "<50" & BMI_GRP  == "<25" ~ 4.1,  
      RIDAGEYR == "<60" & BMI_GRP  == "<25" ~ 8.4,  
      RIDAGEYR == "<70" & BMI_GRP  == "<25" ~ 9.3,  
      RIDAGEYR == "80" & BMI_GRP  == "<25" ~ 13.2,  
       RIDAGEYR == "<30" & BMI_GRP  == "25-27" ~ 4.2,  
      RIDAGEYR == "<40" & BMI_GRP  == "25-27" ~ 9.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "25-27" ~ 9.6,  
      RIDAGEYR == "<60" & BMI_GRP  == "25-27" ~ 6.4,  
      RIDAGEYR == "<70" & BMI_GRP  == "25-27" ~ 10.3,  
      RIDAGEYR == "80" & BMI_GRP  == "25-27" ~ 9.2,  
       RIDAGEYR == "<30" & BMI_GRP  == "27-30" ~ 4.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "27-30" ~ 2.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "27-30" ~ 10.7,  
      RIDAGEYR == "<60" & BMI_GRP  == "27-30" ~ 13.1,  
      RIDAGEYR == "<70" & BMI_GRP  == "27-30" ~ 9.8,  
      RIDAGEYR == "80" & BMI_GRP  == "27-30" ~ 16.7,  
       RIDAGEYR == "<30" & BMI_GRP  == "30-35" ~ 8.3,  
      RIDAGEYR == "<40" & BMI_GRP  == "30-35" ~ 11.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "30-35" ~ 16.1,  
      RIDAGEYR == "<60" & BMI_GRP  == "30-35" ~ 12.3,  
      RIDAGEYR == "<70" & BMI_GRP  == "30-35" ~ 14.0,  
      RIDAGEYR == "80" & BMI_GRP  == "30-35" ~ 14.5,  
      RIDAGEYR == "<30" & BMI_GRP  == "35-40" ~ 7.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "35-40" ~ 11.8,  
      RIDAGEYR == "<50" & BMI_GRP  == "35-40" ~ 17.2,  
      RIDAGEYR == "<60" & BMI_GRP  == "35-40" ~ 14.3,  
      RIDAGEYR == "<70" & BMI_GRP  == "35-40" ~ 18.0,  
      RIDAGEYR == "80" & BMI_GRP  == "35-40" ~ 16.7,  
        RIDAGEYR == "<30" & BMI_GRP  == ">40" ~ 7.0,  
      RIDAGEYR == "<40" & BMI_GRP  == ">40" ~ 16.5,  
      RIDAGEYR == "<50" & BMI_GRP  == ">40" ~ 21.6,  
      RIDAGEYR == "<60" & BMI_GRP  == ">40" ~ 23,  
      RIDAGEYR == "<70" & BMI_GRP  == ">40" ~ 26.7,  
      RIDAGEYR == "80" & BMI_GRP  == ">40" ~ 14
    ))


temp <- data.frame(Adults %>% arrange(BMI_GRP, RIDAGEYR, SEQN) %>%
  group_by(RIDAGEYR, BMI_GRP) %>% mutate(pop=n()) %>%
  select(SEQN, RIDAGEYR, BMI_GRP, percent_ones, pop) %>%
  mutate(lines= round(151*percent_ones/100) ) %>%
  mutate(line=row_number()))


temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>% group_by(RIDAGEYR, BMI_GRP, OSA) %>% count() %>%
   spread(key=OSA, value=n) %>% mutate(OSA=100*`1`/(`1`+`0`)) %>%
  select(-c(`1`, `0`)) %>%
  spread(key=BMI_GRP,value=OSA) %>% 
  select(RIDAGEYR  , `<25`, `25-27` ,`27-30` ,`30-35`, `35-40`, `>40`)



temp %>% mutate(OSA=ifelse(line<=lines,1,0))  %>% group_by(OSA) %>% count()


temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>%
  inner_join(Adults %>% select(SEQN, Liver)) %>%
  group_by(Liver, OSA) %>% count() 


Adults <- temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>% inner_join(Adults) %>%
  select(-c(pop, percent_ones, lines, line))

Adults <- Adults %>% mutate(Total=HTN+CHOL+CKD+T2D+PRET2D+ARTH+CHF+CAD+HEARTATTACK+STROKE+Liver+OSA)

Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  group_by(GROUP) %>% count()

Adults %>%
    mutate(Total=ifelse(Total>=3,"3+", Total)) %>%
  mutate(GROUP=ifelse(Total=="0","0",
                               ifelse(T2D=="1","T2D", Total))) %>%
  group_by(BMI_GRP, GROUP) %>% count() %>%
  spread(key=BMI_GRP, value=n)


unique(Adults$BMI_GRP)

names(Adults)

Adults %>% filter(BMI_GRP!="25-27"&BMI_GRP!="<25") %>%
  summarise(HTN=mean(OSA))


Adults <- Adults %>% select(SEQN, BMI_GRP, Total, T2D)

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/HIQ_L.XPT", tf <- tempfile(), mode="wb")
Health_Insurance_20211_2022 <- foreign::read.xport(tf)[,]

Adults %>% inner_join(Health_Insurance_20211_2022)  %>%
  select(SEQN, BMI_GRP, Total, T2D, HIQ011, HIQ032A , HIQ032B, HIQ032C , HIQ032D, HIQ032F, HIQ032H,HIQ032I   ) %>%
  filter(BMI_GRP !="25-27" & BMI_GRP!="<25") %>%
  mutate(Total=ifelse(Total==0,0,
                      ifelse(Total<=2,1,3))) %>%
  mutate(Total=ifelse(T2D==1, "T2D", Total)) %>%
  mutate(BMI_GRP=paste(BMI_GRP, Total)) %>% select(-Total) %>%
  mutate(HIQ032A=ifelse(is.na(HIQ032A), 9999, HIQ032A)) %>%
  mutate(HIQ032B=ifelse(is.na(HIQ032B), 9999, HIQ032B)) %>%
  mutate(HIQ011 =ifelse(is.na(HIQ011 ), 9999, HIQ011 )) %>%
    mutate(HIQ032C =ifelse(is.na(HIQ032C ), 9999, HIQ032C )) %>%
  mutate(HIQ032D =ifelse(is.na(HIQ032D ), 9999, HIQ032D )) %>%
  mutate(HIQ032F =ifelse(is.na(HIQ032F ), 9999, HIQ032F )) %>%
  mutate(HIQ032H =ifelse(is.na(HIQ032H ), 9999, HIQ032H )) %>%
  mutate(HIQ032I =ifelse(is.na(HIQ032I ), 9999, HIQ032I )) %>%
  mutate(Insurance=ifelse(HIQ032A==1,"Private",
                          ifelse(HIQ032B==2,"Public", 
                                 ifelse(HIQ032C==3,"Public", 
                                        ifelse(HIQ032D==4,"Public", 
                                               ifelse(HIQ032F==6,"Public", 
                                                      ifelse(HIQ032H==8,"Public", 
                                                             ifelse(HIQ032I==9,"Public", "Other/None")))))))) %>%
  group_by(BMI_GRP, Insurance) %>% count() %>% spread(key=Insurance, value=n)


# ----------

# Number of Comorbidities for each comorbdity group -----------------------


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DEMO_L.XPT", tf <- tempfile(), mode="wb")
Demographcics_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Demographcics_20211_2022 %>% select(SEQN, RIDAGEYR, RIAGENDR) %>% filter(RIDAGEYR>=18) %>%
  mutate(RIAGENDR=ifelse(RIAGENDR==1,"M","F"))

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BMX_L.XPT", tf <- tempfile(), mode="wb")
Body_Measures_20211_2022 <- foreign::read.xport(tf)[,]
Adults <- Adults %>% left_join(Body_Measures_20211_2022 %>% select(SEQN, BMXBMI))

Adults <- Adults %>% mutate(BMI_GRP=ifelse(BMXBMI<25, "<25",
                                           ifelse(BMXBMI<27,"25-27",
                                 ifelse(BMXBMI<30,"27-30",
                                        ifelse(BMXBMI<35,"30-35",
                                               ifelse(BMXBMI<40,"35-40",">40"))))))



Adults <- Adults %>% mutate(RIDAGEYR=ifelse(RIDAGEYR<=30, "<30",
                                  ifelse(RIDAGEYR<=40,"<40",
                                         ifelse(RIDAGEYR<=50,"<50",
                                                ifelse(RIDAGEYR<=60,"<60",
                                                       ifelse(RIDAGEYR<=70,"<70",
                                                              ifelse(RIDAGEYR<=80,"80", NA)))))))



# Comorbs
download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/BPQ_L.XPT", tf <- tempfile(), mode="wb")
Blood_Pressure_Cholesterol_20211_2022 <- foreign::read.xport(tf)[,]
Blood_Pressure_Cholesterol_20211_2022 <- Blood_Pressure_Cholesterol_20211_2022 %>% select(SEQN,BPQ020,BPQ080) %>%
  rename("HTN"="BPQ020", "CHOL"="BPQ080")
Adults <- Adults %>% left_join(Blood_Pressure_Cholesterol_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/KIQ_U_L.XPT", tf <- tempfile(), mode="wb")
Kidney_Conditions_Urology_20211_2022 <- foreign::read.xport(tf)[,]
Kidney_Conditions_Urology_20211_2022 <- Kidney_Conditions_Urology_20211_2022 %>% select(SEQN,KIQ022) %>%
  rename("CKD"="KIQ022")
Adults <- Adults %>% left_join(Kidney_Conditions_Urology_20211_2022)

download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/DIQ_L.XPT", tf <- tempfile(), mode="wb")
Diabetes_20211_2022 <- foreign::read.xport(tf)[,]
Diabetes_20211_2022 <- Diabetes_20211_2022 %>% select(SEQN, DIQ010,DIQ160 ) %>%
     rename("T2D"="DIQ010", "PRET2D"="DIQ160")
Adults <- Adults %>% left_join(Diabetes_20211_2022)


download.file("https://wwwn.cdc.gov/Nchs/Nhanes/2021-2022/MCQ_L.XPT", tf <- tempfile(), mode="wb")
Medical_Conditions_20211_2022 <- foreign::read.xport(tf)[,]
Medical_Conditions_20211_2022 <- Medical_Conditions_20211_2022 %>% 
  select(SEQN, MCQ160A, MCQ160B,MCQ160C ,MCQ160E,MCQ160F ,MCQ160L ,MCQ510D ,MCQ510E ,MCQ510F ,MCQ220    ) %>%
     rename("ARTH"="MCQ160A", "CHF"="MCQ160B", "CAD"="MCQ160C", "HEARTATTACK"="MCQ160E", 
            "STROKE"="MCQ160F", "Liver"="MCQ160L", "ViralHepat"="MCQ510D", "AutoHepat"="MCQ510E",
            "OtherHepat"="MCQ510F", "Cancer"="MCQ220" )

Adults <- Adults %>% left_join(Medical_Conditions_20211_2022)



Adults <- Adults %>% 
  mutate(HTN=ifelse(HTN==1,1,0)) %>%
  mutate(CHOL=ifelse(CHOL==1,1,0)) %>%
  mutate(CKD =ifelse(CKD ==1,1,0)) %>%
  mutate(T2D =ifelse(T2D==1|T2D==3,1,0)) %>%
  mutate(PRET2D =ifelse(PRET2D ==1,1,0)) %>%
  mutate(ARTH =ifelse(ARTH ==1,1,0)) %>%
  mutate(CHF =ifelse(CHF ==1,1,0)) %>%
  mutate(CAD =ifelse(CAD ==1,1,0)) %>%
  mutate(HEARTATTACK  =ifelse(HEARTATTACK  ==1,1,0)) %>%
  mutate(STROKE  =ifelse(STROKE  ==1,1,0)) %>%
  mutate(Liver =ifelse(Liver ==1,1,0)) 

Adults <- Adults %>% filter(!is.na(BMI_GRP))


names(Adults)


Adults <- Adults %>% select(-c(ViralHepat, AutoHepat, OtherHepat, Cancer))

Adults[is.na(Adults)] <- 0

# 
#   RIDAGEYR `<25` `25-27` `27-30` `30-35` `35-40` `>40`
# 1 <30       2.56    4.17    4.58    8.25    7.62  6.98
# 2 <40       3.03    9.09    2.08   11.1    11.8  16.5
# 3 <50       4.14    9.57   10.7    16.1    17.2  21.6
# 4 <60       8.37    6.36   13.1    12.3    14.3  23
# 5 <70       9.31   10.3     9.84   14.0    18.0  26.7
# 6 80       13.2     9.23   16.7    14.5    16.7  28.6



unique(Adults$RIDAGEYR)
unique(Adults$BMI_GRP)


Adults <- Adults %>%
  mutate(
    percent_ones = case_when(
      RIDAGEYR == "<30" & BMI_GRP  == "<25" ~ 2.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "<25" ~ 3.0,  
      RIDAGEYR == "<50" & BMI_GRP  == "<25" ~ 4.1,  
      RIDAGEYR == "<60" & BMI_GRP  == "<25" ~ 8.4,  
      RIDAGEYR == "<70" & BMI_GRP  == "<25" ~ 9.3,  
      RIDAGEYR == "80" & BMI_GRP  == "<25" ~ 13.2,  
       RIDAGEYR == "<30" & BMI_GRP  == "25-27" ~ 4.2,  
      RIDAGEYR == "<40" & BMI_GRP  == "25-27" ~ 9.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "25-27" ~ 9.6,  
      RIDAGEYR == "<60" & BMI_GRP  == "25-27" ~ 6.4,  
      RIDAGEYR == "<70" & BMI_GRP  == "25-27" ~ 10.3,  
      RIDAGEYR == "80" & BMI_GRP  == "25-27" ~ 9.2,  
       RIDAGEYR == "<30" & BMI_GRP  == "27-30" ~ 4.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "27-30" ~ 2.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "27-30" ~ 10.7,  
      RIDAGEYR == "<60" & BMI_GRP  == "27-30" ~ 13.1,  
      RIDAGEYR == "<70" & BMI_GRP  == "27-30" ~ 9.8,  
      RIDAGEYR == "80" & BMI_GRP  == "27-30" ~ 16.7,  
       RIDAGEYR == "<30" & BMI_GRP  == "30-35" ~ 8.3,  
      RIDAGEYR == "<40" & BMI_GRP  == "30-35" ~ 11.1,  
      RIDAGEYR == "<50" & BMI_GRP  == "30-35" ~ 16.1,  
      RIDAGEYR == "<60" & BMI_GRP  == "30-35" ~ 12.3,  
      RIDAGEYR == "<70" & BMI_GRP  == "30-35" ~ 14.0,  
      RIDAGEYR == "80" & BMI_GRP  == "30-35" ~ 14.5,  
      RIDAGEYR == "<30" & BMI_GRP  == "35-40" ~ 7.6,  
      RIDAGEYR == "<40" & BMI_GRP  == "35-40" ~ 11.8,  
      RIDAGEYR == "<50" & BMI_GRP  == "35-40" ~ 17.2,  
      RIDAGEYR == "<60" & BMI_GRP  == "35-40" ~ 14.3,  
      RIDAGEYR == "<70" & BMI_GRP  == "35-40" ~ 18.0,  
      RIDAGEYR == "80" & BMI_GRP  == "35-40" ~ 16.7,  
        RIDAGEYR == "<30" & BMI_GRP  == ">40" ~ 7.0,  
      RIDAGEYR == "<40" & BMI_GRP  == ">40" ~ 16.5,  
      RIDAGEYR == "<50" & BMI_GRP  == ">40" ~ 21.6,  
      RIDAGEYR == "<60" & BMI_GRP  == ">40" ~ 23,  
      RIDAGEYR == "<70" & BMI_GRP  == ">40" ~ 26.7,  
      RIDAGEYR == "80" & BMI_GRP  == ">40" ~ 14
    ))


temp <- data.frame(Adults %>% arrange(BMI_GRP, RIDAGEYR, SEQN) %>%
  group_by(RIDAGEYR, BMI_GRP) %>% mutate(pop=n()) %>%
  select(SEQN, RIDAGEYR, BMI_GRP, percent_ones, pop) %>%
  mutate(lines= round(151*percent_ones/100) ) %>%
  mutate(line=row_number()))


temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>% group_by(RIDAGEYR, BMI_GRP, OSA) %>% count() %>%
   spread(key=OSA, value=n) %>% mutate(OSA=100*`1`/(`1`+`0`)) %>%
  select(-c(`1`, `0`)) %>%
  spread(key=BMI_GRP,value=OSA) %>% 
  select(RIDAGEYR  , `<25`, `25-27` ,`27-30` ,`30-35`, `35-40`, `>40`)



temp %>% mutate(OSA=ifelse(line<=lines,1,0))  %>% group_by(OSA) %>% count()


temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>%
  inner_join(Adults %>% select(SEQN, Liver)) %>%
  group_by(Liver, OSA) %>% count() 


Adults <- temp %>% mutate(OSA=ifelse(line<=lines,1,0)) %>% inner_join(Adults) %>%
  select(-c(pop, percent_ones, lines, line))

Adults <- Adults %>% mutate(Total=HTN+CHOL+CKD+T2D+PRET2D+ARTH+CHF+CAD+HEARTATTACK+STROKE+Liver+OSA)


Adults <- Adults %>% select(-c(BMI_GRP, RIDAGEYR, RIAGENDR, BMXBMI))

Adults <- Adults %>% select(SEQN, Total, OSA:Liver)



comorbidities <- c("OSA", "HTN", "CHOL", "CKD", "T2D", "PRET2D", "ARTH", "CHF", "CAD", "HEARTATTACK", "STROKE", "Liver")

# Initialize a result dataframe
results <- data.frame(Comorbidity = character(),
                      Avg_Total = numeric(),
                      stringsAsFactors = FALSE)

# Loop through each comorbidity and calculate the average Total for patients with that comorbidity
for (comorbidity in comorbidities) {
  # Filter patients with the comorbidity
  patients_with_comorbidity <- Adults[Adults[[comorbidity]] == 1, ]
  
  # Calculate the average Total
  avg_total <- mean(patients_with_comorbidity$Total, na.rm = TRUE)
  
  # Add the result to the dataframe
  results <- rbind(results, data.frame(Comorbidity = comorbidity, Avg_Total = avg_total))
}

# View results
print(results)

   Comorbidity Avg_Total
1          OSA  2.941456
2          HTN  3.245051
3         CHOL  3.085415
4          CKD  4.429224
5          T2D  3.667606
6       PRET2D  3.031073
7         ARTH  3.169500
8          CHF  5.240816
9          CAD  4.964968
10 HEARTATTACK  5.023715
11      STROKE  4.297398
12       Liver  3.846591


comorbidity_matrix <- matrix(0, 
                             nrow = length(comorbidities), 
                             ncol = length(comorbidities),
                             dimnames = list(comorbidities, comorbidities))

# Populate the matrix
for (row_comorbidity in comorbidities) {
  for (col_comorbidity in comorbidities) {
    # Count patients with both comorbidities
    count <- sum(Adults[[row_comorbidity]] == 1 & Adults[[col_comorbidity]] == 1, na.rm = TRUE)
    # Assign the count to the matrix
    comorbidity_matrix[row_comorbidity, col_comorbidity] <- count
  }
}

# View the resulting matrix
print(comorbidity_matrix)

# ---------


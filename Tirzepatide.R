
# Tirzepatide ---------------------------------------------
# ---------------------------------------------
# Packages ----------------------------------------
# install.packages(c("DBI", "dplyr","dbplyr","odbc"))

library(data.table)
library(DBI)
library(dplyr)
library(dbplyr)
library(odbc)
library(tidyverse)
# ---------------------------------------------------------
# Establish connection ------------------------------

myconn <- DBI::dbConnect(odbc::odbc(), 
                         "RWD PROD AMER on Snowflake", 
                         uid="********",
                         pwd='********')

# -------------------------------------

# Tirzepatide patients -------------------------



# CLAIMS_RX 
# CLAIMS_LAB_RESULTS
# CLAIMS_MEDICAL
# CLAIMS_INPATIENT_CONFINEMENT
# NLP_MEASUREMENTCLAIMS_CONT_ENROLL
# NLP_MEASUREMENT
# NLP_SDS
# NLP_DRUG_RATIONALE
# OBSERVATIONS
# PROCEDURES


start <- Sys.time()
tirzepatide_pats <- DBI::dbGetQuery(myconn,"SELECT DISTINCT PTID FROM RWD_PROD.MARKET_CLARITY.CLAIMS_RX 
                                            WHERE BRND_NM='ZEPBOUND' AND TRANSACTION_TYPE='PAID';") 
end <- Sys.time()
print(end-start)
names(tirzepatide_pats)
fwrite(tirzepatide_pats, "Source/Tirzepatide_patients.txt")

# -----------------------
# CLAIMS_RX Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

claims_rx <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
 
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.CLAIMS_RX WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  claims_rx <- rbind(claims_rx, data_df)

  if(i %% 1000 ==0){
    fwrite(claims_rx, paste0("Source/Temp/claims_rx_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(claims_rx)))
}

end <- Sys.time()

print(end-start)

length(unique(claims_rx$PTID))

fwrite(claims_rx, "Source/tirzepatide_claims_rx.txt")

# -----------------------
# CLAIMS_LAB_RESULTS Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

claims_lab_results <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.CLAIMS_LAB_RESULTS WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  claims_lab_results <- rbind(claims_lab_results, data_df)
  
  if(i %% 1000 ==0){
    fwrite(claims_lab_results, paste0("Source/Temp/claims_lab_results_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(claims_lab_results)))
}

end <- Sys.time()

print(end-start)

length(unique(claims_lab_results$PTID))

fwrite(claims_lab_results, "Source/claims_lab_results.txt")

# -----------------------
# CLAIMS_MEDICAL Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

claims_medical <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.CLAIMS_MEDICAL WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  claims_medical <- rbind(claims_medical, data_df)
  
  if(i %% 1000 ==0){
    fwrite(claims_medical, paste0("Source/Temp/claims_medical_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(claims_medical)))
}

end <- Sys.time()

print(end-start)

length(unique(claims_medical$PTID))

fwrite(claims_medical, "Source/claims_medical.txt")

# -----------------------
# CLAIMS_INPATIENT_CONFINEMENT Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

claims_inpatient_confinement <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.CLAIMS_INPATIENT_CONFINEMENT WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  claims_inpatient_confinement <- rbind(claims_inpatient_confinement, data_df)
  
  if(i %% 1000 ==0){
    fwrite(claims_inpatient_confinement, paste0("Source/Temp/claims_inpatient_confinement_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(claims_inpatient_confinement)))
}

end <- Sys.time()

print(end-start)

length(unique(claims_inpatient_confinement$PTID))

fwrite(claims_inpatient_confinement, "Source/claims_inpatient_confinement.txt")

# -----------------------
# CLAIMS_CONT_ENROLL Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

claims_cont_enroll <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.CLAIMS_CONT_ENROLL WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  claims_cont_enroll <- rbind(claims_cont_enroll, data_df)
  
  if(i %% 1000 ==0){
    fwrite(claims_cont_enroll, paste0("Source/Temp/claims_cont_enroll_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(claims_cont_enroll)))
}

end <- Sys.time()

print(end-start)

length(unique(claims_cont_enroll$PTID))

fwrite(claims_cont_enroll, "Source/claims_cont_enroll.txt")

# -----------------------
# NLP_MEASUREMENT Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

nlp_measurement <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.NLP_MEASUREMENT WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  nlp_measurement <- rbind(nlp_measurement, data_df)
  
  if(i %% 1000 ==0){
    fwrite(nlp_measurement, paste0("Source/Temp/nlp_measurement_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(nlp_measurement)))
}

end <- Sys.time()

print(end-start)

length(unique(nlp_measurement$PTID))

fwrite(nlp_measurement, "Source/nlp_measurement.txt")

# -----------------------

# NLP_SDS Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

nlp_sds <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.NLP_SDS WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  nlp_sds <- rbind(nlp_sds, data_df)
  
  if(i %% 1000 ==0){
    fwrite(nlp_sds, paste0("Source/Temp/nlp_sds_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(nlp_sds)))
}

end <- Sys.time()

print(end-start)

length(unique(nlp_sds$PTID))

fwrite(nlp_sds, "Source/nlp_sds.txt")

# -----------------------
# NLP_DRUG_RATIONALE Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

nlp_drug_rationale <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.NLP_DRUG_RATIONALE WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  nlp_drug_rationale <- rbind(nlp_drug_rationale, data_df)
  
  if(i %% 1000 ==0){
    fwrite(nlp_drug_rationale, paste0("Source/Temp/nlp_drug_rationale_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(nlp_drug_rationale)))
}

end <- Sys.time()

print(end-start)

length(unique(nlp_drug_rationale$PTID))

fwrite(nlp_drug_rationale, "Source/nlp_drug_rationale")

# -----------------------

# OBSERVATIONS Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

observations <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.OBSERVATIONS WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  observations <- rbind(observations, data_df)
  
  if(i %% 1000 ==0){
    fwrite(observations, paste0("Source/Temp/observations_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(observations)))
}

end <- Sys.time()

print(end-start)

length(unique(observations$PTID))

fwrite(observations, "Source/observations.txt")

# -----------------------



# PROCEDURES Tirzepatide patients -----------------------------

pagify <- function(data = NULL, by = 1000){
  pagemin <- seq(1, length(data), by = by)
  pagemax <- pagemin - 1 + by
  pagemax[length(pagemax)] <- length(data)
  pages <- list(min = pagemin, max = pagemax)
}


tirzepatide_pats <- fread("Source/Tirzepatide_patients.txt")


pages <- pagify(tirzepatide_pats$PTID, 14069)

procedures <- data.table()

start <- Sys.time()

for(i in length(pages$max)){
  
  pts <- paste0(tirzepatide_pats$PTID[pages$min[i]:pages$max[i]], collapse = "','")
  query <- paste0("SELECT * FROM RWD_PROD.MARKET_CLARITY.PROCEDURES WHERE PTID IN ('",pts,"');")
  data_df <- setDT(DBI::dbGetQuery(myconn, query))
  procedures <- rbind(procedures, data_df)
  
  if(i %% 1000 ==0){
    fwrite(procedures, paste0("Source/Temp/procedures_",i,".txt"))
  }
  
  print(paste0(i, " of ", length(pages$max), " rows: ", nrow(procedures)))
}

end <- Sys.time()

print(end-start)

length(unique(procedures$PTID))

fwrite(procedures, "Source/procedures.txt")

# -----------------------
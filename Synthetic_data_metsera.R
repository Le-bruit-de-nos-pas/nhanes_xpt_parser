DANU_Measures <- fread("DANU Measures 1.1/DANU Measures.txt",  integer64 = "character", stringsAsFactors = F)

unique(DANU_Measures$test)
DANU_Measures <- DANU_Measures %>% filter(test=="BMI") %>% select(patid, weight, value, claimed) %>% distinct()

DANU_Measures <- DANU_Measures %>% group_by(patid, weight) %>% summarise(value=mean(value)) %>% distinct() %>% ungroup()
DANU_Measures <- DANU_Measures %>% filter(value>15&value<75)

DANU_Demographics <- fread("DANU Demographics 1.1/DANU Demographics.txt")
DANU_Demographics <- DANU_Demographics %>% select(patid, age, gender) 

DANU_Measures <- DANU_Demographics %>% inner_join(DANU_Measures) 

mean(DANU_Measures$value) # 30.64611
median(DANU_Measures$value)  # 29.3368
sd(DANU_Measures$value)  # 6.676264

e1071::kurtosis(DANU_Measures$value)  # 2.744577
e1071::skewness(DANU_Measures$value)  # 1.251627

DANU_Measures %>% ggplot(aes(value)) + geom_density()
range(DANU_Measures$value) # 15.00931 74.99776

DANU_Measures %>% group_by(gender) %>% summarise(mean=mean(value))

# 1 F       30.8
# 2 M       30.4

DANU_Measures %>% group_by(gender) %>% summarise(sd=sd(value))

# 1 F       7.25
# 2 M       5.90

DANU_Measures %>% group_by(gender) %>% summarise(kurtosis=e1071::kurtosis(value))

# 1 F          2.07
# 2 M          3.71

DANU_Measures %>% group_by(gender) %>% summarise(kurtosis=e1071::skewness(value))
# 1 F          2.07
# 2 M          3.71

DANU_Measures %>% group_by(gender) %>% summarise(age=mean(age))

summary(DANU_Measures$age)

DANU_Measures %>%
  ggplot(aes(age, colour=gender, fill=gender)) +
  geom_density(alpha=0.5) + theme_minimal()


range(DANU_Measures$age)

DANU_Measures %>% mutate(age=ifelse(age<=25,"25",
                                    ifelse(age<=40,"40",
                                           ifelse(age<=65,"65","99")))) %>%
  group_by(age) %>% count() %>% ungroup() %>% mutate(tot=sum(n)) %>% mutate(perc=n/tot)






n_patients <- 165000  # Number of rows/patients
mean_bmi <- 29  # Mean on the original scale
sd_bmi <- 6     # Standard deviation on the original scale
# Convert to log-normal parameters
mu_log <- log(mean_bmi^2 / sqrt(sd_bmi^2 + mean_bmi^2))
sigma_log <- sqrt(log(1 + (sd_bmi^2 / mean_bmi^2)))
# Generate log-normal BMI values
set.seed(123)  # For reproducibility
bmi_values <- rlnorm(n_patients, meanlog = mu_log, sdlog = sigma_log)
# Summary statistics to check
summary(bmi_values)
BMI_males <- data.frame(BMI = bmi_values)
BMI_males$gender <- "M"

n_patients <- 185000  # Number of rows/patients
mean_bmi <- 30  # Mean on the original scale
sd_bmi <- 7     # Standard deviation on the original scale
# Convert to log-normal parameters
mu_log <- log(mean_bmi^2 / sqrt(sd_bmi^2 + mean_bmi^2))
sigma_log <- sqrt(log(1 + (sd_bmi^2 / mean_bmi^2)))
# Generate log-normal BMI values
set.seed(123)  # For reproducibility
bmi_values <- rlnorm(n_patients, meanlog = mu_log, sdlog = sigma_log)
# Summary statistics to check
summary(bmi_values)
BMI_females <- data.frame(BMI = bmi_values)
BMI_females$gender <- "F"



age_buckets <- c(18:25, 26:40, 41:65, 66:80)

age_probs <- c(rep(0.09 / length(18:25), length(18:25)),   # 9% for ages 18-25
               rep(0.15 / length(26:40), length(26:40)),   # 15% for ages 26-40
               rep(0.51 / length(41:65), length(41:65)),   # 51% for ages 41-65
               rep(0.25 / length(66:80), length(66:80)))   # 25% for ages 66-80

# Generate age values for males
age_males <- sample(age_buckets, size = nrow(BMI_males), replace = TRUE, prob = age_probs)

# Generate age values for females
age_females <- sample(age_buckets, size = nrow(BMI_females), replace = TRUE, prob = age_probs)

# Add age to the datasets
BMI_males$age <- age_males
BMI_females$age <- age_females




BMI_females %>%
  bind_rows(BMI_males) %>%
  ggplot(aes(BMI, colour=gender, fill=gender)) +
  geom_density(alpha=0.5) + theme_minimal()




BMI_females %>%
  bind_rows(BMI_males) %>%
  ggplot(aes(age, colour=gender, fill=gender)) +
  geom_density(alpha=0.5) + theme_minimal()

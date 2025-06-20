library(tidyverse)

select <- dplyr::select

################################################################################
# Load, preprocess data
################################################################################
setwd("supporting_code")

save_dir <- paste(getwd(), "data", sep="/")

demog <- read.csv("data/Participant/Primary/participant.sept23.d090723.csv")
asa24 <- read.csv("data/ASA24/Daily Total Nutrients/asa24_tn.d071919.csv")

################################################################################
# Exclusions
################################################################################
demog <-
  filter(demog, (cv1_bmi_cat != 4) & (cv2_bmi_cat != 4) & (cv3_bmi_cat != 4))

################################################################################
# How much missingness in there in BMI?
################################################################################
dat <- demog %>% select(iid, age=AGE, sex,
                        cv1_bmi, cv2_bmi, cv3_bmi)
dat$cv1_bmi[dat$cv1_bmi=="N"] <- NA
dat$cv2_bmi[dat$cv2_bmi=="N"] <- NA
dat$cv3_bmi[dat$cv3_bmi=="N"] <- NA

dat$cv1_bmi <- as.numeric(dat$cv1_bmi)
dat$cv2_bmi <- as.numeric(dat$cv2_bmi)
dat$cv3_bmi <- as.numeric(dat$cv3_bmi)

################################################################################
# Merge demographics with ASA24 after filtering only task 1 and 2
################################################################################
# asa24 <- filter(asa24, task <= 2)
merge_dat <- merge(demog, asa24, by="iid", suffixes=c(".demog", ".asa24"))

################################################################################
# flag nutrient fields (needs to be verified if any changes to the merged data
# frame)
################################################################################
nutrient_cols <- 130:194

################################################################################
# Any missingness in this data set? Note: this line of code works with the 
# assumption that the fields are numeric, which seems to be true.
################################################################################
sum(complete.cases(merge_dat[, nutrient_cols]))

################################################################################
# Include only the wanted nutrient columns
################################################################################
nutrients_inc <- c("tn_acar_asa24", "tn_alc_asa24", "tn_atoc_asa24", "tn_b12_add_asa24",
                   "tn_bcar_asa24", "tn_caff_asa24", "tn_calc_asa24", "tn_carb_asa24",
                   "tn_chole_asa24", "tn_choln_asa24", "tn_copp_asa24", "tn_cryp_asa24",
                   "tn_fdfe_asa24", "tn_fibe_asa24", "tn_iron_asa24", "tn_lyco_asa24",
                   "tn_lz_asa24", "tn_magn_asa24", "tn_mfat_asa24",
                   "tn_niac_asa24", "tn_pfat_asa24", "tn_phos_asa24", "tn_pota_asa24",
                   "tn_prot_asa24", "tn_sele_asa24", "tn_sfat_asa24", "tn_sodi_asa24",
                   "tn_theo_asa24", "tn_vara_asa24", "tn_vb12_asa24",
                   "tn_vb1_asa24", "tn_vb2_asa24", "tn_vb6_asa24", "tn_vc_asa24",
                   "tn_vitd_asa24", "tn_vite_add_asa24", "tn_vk_asa24", "tn_zinc_asa24")

################################################################################
# Create newer merged data with only included nutrients
################################################################################
primary_inc <- c("iid", "AGE", "sex", "cv1_bmi", "cv2_bmi", "cv3_bmi", "task")
inc_col <- c(primary_inc, nutrients_inc)
analysis_data <- merge_dat[, inc_col]

################################################################################
# Check classical error homoscedasticity assumption
################################################################################
log_analysis_data <- analysis_data
log_analysis_data[, nutrients_inc] <-
  log( analysis_data[, nutrients_inc] + 1 )

#' to look at raw data additive error assumption, use analysis_data instead of
#' log_analysis_data in the line below:
summ_dat <- log_analysis_data %>% group_by(iid) %>%
  summarise(across( .cols=all_of(nutrients_inc),
                    .fns=list(mean=mean, sd=sd),
                    .names="{.col}.{.fn}")
            )
summ_dat <- summ_dat[complete.cases(summ_dat), ]

summ_dat_grouped <- summ_dat %>% pivot_longer(
  cols=-iid,
  names_to=c("nutrient", "summary"),
  names_sep="\\.",
  values_to="level"
  ) %>%
  pivot_wider(names_from=summary, values_from=level) %>%
  arrange(iid, nutrient)

ggplot(
  summ_dat_grouped,
  aes(x=mean, y=sd)
  ) +
  geom_point() +
  facet_wrap(~nutrient, scales="free") +
  theme_minimal()

nutrients_exc_additive <- c("tn_acar_asa24", "tn_alc_asa24", "tn_b12_add_asa24",
                            "tn_caff_asa24", "tn_cryp_asa24", "tn_lyco_asa24",
                            "tn_sodi_asa24", "tn_theo_asa24", "tn_vitd_asa24",
                            "tn_vite_add_asa24")

nutrients_inc_final <- setdiff(nutrients_inc, nutrients_exc_additive)

inc_col <- c(primary_inc, nutrients_inc_final)
log_analysis_data <- log_analysis_data[, inc_col]

################################################################################
# Check classical error assumption on correlation with age, sex
################################################################################
summ_dat_grouped <- summ_dat_grouped %>%
  filter(nutrient %in% nutrients_inc_final)

merge_summ_grouped <- merge(summ_dat_grouped, dat, by="iid") %>%
  select(-c("mean", "cv1_bmi", "cv2_bmi", "cv3_bmi"))

ggplot(
  merge_summ_grouped,
  aes(x=age, y=sd)
) +
  geom_point() +
  facet_wrap(~nutrient, scales="free") +
  theme_minimal()

ggplot(
  merge_summ_grouped,
  aes(x=sex, y=sd)
) +
  geom_point() +
  facet_wrap(~nutrient, scales="free") +
  theme_minimal()

################################################################################
# log_analysis_data will be used for analysis. Continue parsing
################################################################################
log_analysis_data$cv1_bmi[log_analysis_data$cv1_bmi=="N"] <- NA
log_analysis_data$cv2_bmi[log_analysis_data$cv2_bmi=="N"] <- NA
log_analysis_data$cv3_bmi[log_analysis_data$cv3_bmi=="N"] <- NA

log_analysis_data$cv1_bmi <- as.numeric(log_analysis_data$cv1_bmi)
log_analysis_data$cv2_bmi <- as.numeric(log_analysis_data$cv2_bmi)
log_analysis_data$cv3_bmi <- as.numeric(log_analysis_data$cv3_bmi)

log_analysis_data <- log_analysis_data %>% 
  pivot_longer(
    cols = c("cv1_bmi", "cv2_bmi", "cv3_bmi"),
    names_to = "month",
    values_to = "bmi"
  )

log_analysis_data$month[log_analysis_data$month=="cv1_bmi"] <- 0
log_analysis_data$month[log_analysis_data$month=="cv2_bmi"] <- 5
log_analysis_data$month[log_analysis_data$month=="cv3_bmi"] <- 11
log_analysis_data$month <- as.numeric(log_analysis_data$month)
log_analysis_data <- log_analysis_data[complete.cases(log_analysis_data),]

log_analysis_data$sex <- log_analysis_data$sex - 1

log_analysis_data <- log_analysis_data %>%
  mutate(age=AGE) %>%
  select(-c(AGE))

write.csv(log_analysis_data,
          file=paste(save_dir, "log_analysis_data.csv", sep="/"),
          row.names=FALSE)

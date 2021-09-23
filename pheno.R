library('data.table')
library('dplyr')
source("funs.R")
set.seed(1234)

# load phenotypes
disease_id <- paste0("41270-0.", seq(0, 212))
disease_name <- paste0("diagnoses_icd10_41270.0.", seq(0, 212))

f <- "/tmp/tmp.FB5J2pN0Te/data.33352.csv"
pheno <- fread(f, select=c(
        "eid",
        "31-0.0",
        "21022-0.0",
        "21001-0.0",
        "30620-0.0",
        "30600-0.0",
        "30610-0.0",
        "30630-0.0",
        "30640-0.0",
        "30650-0.0",
        "30710-0.0",
        "30680-0.0",
        "30690-0.0",
        "30700-0.0",
        "30720-0.0",
        "30660-0.0",
        "30730-0.0",
        "30740-0.0",
        "30750-0.0",
        "30760-0.0",
        "30770-0.0",
        "30780-0.0",
        "30790-0.0",
        "30800-0.0",
        "30810-0.0",
        "30820-0.0",
        "30830-0.0",
        "30850-0.0",
        "30840-0.0",
        "30860-0.0",
        "30870-0.0",
        "30880-0.0",
        "30670-0.0",
        "30890-0.0",
        "30140-0.0",
        "30120-0.0",
        "3581-0.0",
        "2814-0.0",
        "20116-0.0",
        "22034-0.0",
        "1558-0.0",
        "100004-0.0",
        "100008-0.0",
        "6150-0.0",
        "6150-0.1",
        "6150-0.2",
        "6150-0.3",
        "129-0.0",
        "130-0.0",
        disease_id
    ),
    col.names=c(
        "eid", 
        "sex.31.0.0",
        "age_at_recruitment.21022.0.0",
        "body_mass_index.21001.0.0",
        "alanine_aminotransferase.30620.0.0",
        "albumin.30600.0.0",
        "alkaline_phosphatase.30610.0.0",
        "apolipoprotein_a.30630.0.0",
        "apolipoprotein_b.30640.0.0",
        "aspartate_aminotransferase.30650.0.0",
        "c_reactive_protein.30710.0.0",
        "calcium.30680.0.0",
        "cholesterol.30690.0.0",
        "creatinine.30700.0.0",
        "cystatin_c.30720.0.0",
        "direct_bilirubin.30660.0.0",
        "gamma_glutamyltransferase.30730.0.0",
        "glucose.30740.0.0",
        "glycated_haemoglobin.30750.0.0",
        "hdl_cholesterol.30760.0.0",
        "igf_1.30770.0.0",
        "ldl_direct.30780.0.0",
        "lipoprotein_a.30790.0.0",
        "oestradiol.30800.0.0",
        "phosphate.30810.0.0",
        "rheumatoid_factor.30820.0.0",
        "shbg.30830.0.0",
        "testosterone.30850.0.0",
        "total_bilirubin.30840.0.0",
        "total_protein.30860.0.0",
        "triglycerides.30870.0.0",
        "urate.30880.0.0",
        "urea.30670.0.0",
        "vitamin_d.30890.0.0",
        "neutrophill_count.30140.0.0",
        "lymphocyte_count.30120.0.0",
        "age_at_menopause.3581.0.0",
        "ever_used_hormone_replacement_therapy.2814.0.0",
        "smoking_status.20116.0.0",
        "summed_minutes_activity.22034.0.0",
        "alcohol_intake_frequency.1558.0.0",
        "estimated_fat_yesterday.100004.0.0",
        "estimated_total_sugars_yesterday.100008.0.0",
        "vascular_heart_problems_diagnosed_by_doctor.6150.0.0",
        "vascular_heart_problems_diagnosed_by_doctor.6150.0.1",
        "vascular_heart_problems_diagnosed_by_doctor.6150.0.2",
        "vascular_heart_problems_diagnosed_by_doctor.6150.0.3",
        "place_of_birth_in_UK_north_co_ordinate.129.0.0",
        "place_of_birth_in_UK_east_co_ordinate.130.0.0",
        disease_name
    )
)
unlink(f)

# process phenotypes
pheno$neutrophill_to_lymphocyte_count_ratio <- pheno$neutrophill_count.30140.0.0 / pheno$lymphocyte_count.30120.0.0
pheno <- pheno %>% mutate_at(c('age_at_menopause.3581.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('age_at_menopause.3581.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('ever_used_hormone_replacement_therapy.2814.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('ever_used_hormone_replacement_therapy.2814.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('smoking_status.20116.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('alcohol_intake_frequency.1558.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('vascular_heart_problems_diagnosed_by_doctor.6150.0.0'), na_if, -3)
pheno <- pheno %>% mutate_at(c('vascular_heart_problems_diagnosed_by_doctor.6150.0.1'), na_if, -3)
pheno <- pheno %>% mutate_at(c('vascular_heart_problems_diagnosed_by_doctor.6150.0.2'), na_if, -3)
pheno <- pheno %>% mutate_at(c('vascular_heart_problems_diagnosed_by_doctor.6150.0.3'), na_if, -3)
pheno <- pheno %>% mutate_at(c('place_of_birth_in_UK_north_co_ordinate.129.0.0'), na_if, -1)
pheno <- pheno %>% mutate_at(c('place_of_birth_in_UK_east_co_ordinate.130.0.0'), na_if, -1)

# 6150
pheno$heart_attack.6150 <- apply(pheno[,c('vascular_heart_problems_diagnosed_by_doctor.6150.0.0', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.1', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.2', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.3'),with=F], 1, function(x) {sum(x==1, na.rm=T)>0})
pheno$angina.6150 <- apply(pheno[,c('vascular_heart_problems_diagnosed_by_doctor.6150.0.0', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.1', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.2', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.3'),with=F], 1, function(x) {sum(x==2, na.rm=T)>0})
pheno$stroke.6150 <- apply(pheno[,c('vascular_heart_problems_diagnosed_by_doctor.6150.0.0', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.1', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.2', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.3'),with=F], 1, function(x) {sum(x==3, na.rm=T)>0})
pheno$hypertension.6150 <- apply(pheno[,c('vascular_heart_problems_diagnosed_by_doctor.6150.0.0', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.1', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.2', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.3'),with=F], 1, function(x) {sum(x==4, na.rm=T)>0})
pheno$none_of_the_above.6150 <- apply(pheno[,c('vascular_heart_problems_diagnosed_by_doctor.6150.0.0', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.1', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.2', 'vascular_heart_problems_diagnosed_by_doctor.6150.0.3'),with=F], 1, function(x) {sum(x==-7, na.rm=T)>0})
pheno[is.na(pheno$vascular_heart_problems_diagnosed_by_doctor.6150.0.0)]$heart_attack.6150 <- NA
pheno[is.na(pheno$vascular_heart_problems_diagnosed_by_doctor.6150.0.0)]$angina.6150 <- NA
pheno[is.na(pheno$vascular_heart_problems_diagnosed_by_doctor.6150.0.0)]$stroke.6150 <- NA
pheno[is.na(pheno$vascular_heart_problems_diagnosed_by_doctor.6150.0.0)]$hypertension.6150 <- NA
pheno[is.na(pheno$vascular_heart_problems_diagnosed_by_doctor.6150.0.0)]$none_of_the_above.6150 <- NA
pheno$vascular_problems.6150 <- !pheno$none_of_the_above.6150

# derive binary disease
pheno$liver_disease <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "K7"), na.rm=T)>0})
pheno$alcoholic_liver_disease <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "K70"), na.rm=T)>0})
pheno$fibrosis_liver_disease <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "K74"), na.rm=T)>0})
pheno$fatty_liver_disease <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "K760"), na.rm=T)>0})

pheno$CKD <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "N18"), na.rm=T)>0})
pheno$gout <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "M10"), na.rm=T)>0})
pheno$T2DM <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "E11"), na.rm=T)>0})
pheno$HF <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "I50"), na.rm=T)>0})
pheno$MI <- apply(pheno %>% select(all_of(disease_name)), 1, function(x) {sum(startsWith(x, "I21"), na.rm=T)>0})

# save data
save.image(file = "data/pheno.RData")
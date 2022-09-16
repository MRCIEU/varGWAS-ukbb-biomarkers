# recode filename
westerman <- westerman %>% 
    dplyr::mutate(outcome=dplyr::recode(
        f,
        "albumin.csv"="ALB",
        "ALP.csv"="ALP",
        "ALT.csv"="ALT",
        "apoa.csv"="ApoA",
        "apob.csv"="ApoB",
        "AST.csv"="AST",
        "Billirubin_direct.csv"="Direct BR",
        "Billirubin_total.csv"="Total BR",
        "Creatinine.csv"="Creatinine",
        "CRP.csv"="CRP",
        "Cystatin C.csv"="Cystatin C",
        "GGT.csv"="GGT",
        "Glucose.csv"="Glucose",
        "HbA1C.csv"="HbA1C",
        "HDL.csv"="HDL",
        "LDL.csv"="LDL",
        "LIPA.csv"="LipoA",
        "TC.csv"="TC",
        "TG.csv"="TG",
        "Urate.csv"="Urate",
        .default = NA_character_
    ))

# add key
westerman$key <- paste0(westerman$SNP, "_", westerman$outcome)
lyon$key <- paste0(lyon$rsid, "_", lyon$outcome)

# filter lyon
lyon <- lyon %>% dplyr::filter(outcome %in% westerman$outcome)

# count lyon in westerman
merge(lyon, westerman, "key")
library('data.table')
library('dplyr')
source("funs.R")
set.seed(1234)

# load phenotypes
f <- "/tmp/tmp.785ctU0Nrl/data.33352.csv"
pheno <- fread(f, select=c(
        "eid",
        "31-0.0",
        "21022-0.0",
        "21001-0.0"
    ),
    col.names=c(
        "eid", 
        "sex.31.0.0",
        "age_at_recruitment.21022.0.0",
        "body_mass_index.21001.0.0"
    )
)
unlink(f)

# save data
save.image(file = "data/pheno.RData")

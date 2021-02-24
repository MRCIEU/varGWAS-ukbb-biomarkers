# GWAS of trait variance

## Phenotypes

Copy pheno file to ramdisk

```sh
salloc --nodes=1 --cpus-per-task=21 --mem=80G --time=06:00:00 --partition=mrcieu
d=$(mktemp -d)
echo "copying pheno file to $d"
cp /mnt/storage/private/mrcieu/data/ukbiobank/phenotypic/applications/15825/2019-05-02/data/data.33352.csv "$d"/
```

Extract phenotypes

```sh
module load languages/r/3.6.0
Rscript pheno.R
```

Save variables to csv for vGWAS

```sh
module load languages/r/3.6.0
Rscript extract.R
```

## vGWAS

Run C++ implementation of B-P

```sh
~/projects/jlst_cpp/build/bin/jlst_cpp \
-v data/ukb_bmi.txt \
-s , \
-o data/ukb_bmi.vgwas.txt \
-b /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen \
-p body_mass_index.21001.0.0 \
-i appieu \
-t 1
```
-c sex.31.0.0,age_at_recruitment.21022.0.0,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
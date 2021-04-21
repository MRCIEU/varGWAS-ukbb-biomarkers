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
for trait in alanine_aminotransferase.30620 albumin.30600 alkaline_phosphatase.30610 apolipoprotein_a.30630 apolipoprotein_b.30640 aspartate_aminotransferase.30650 c-reactive_protein.30710 calcium.30680 cholesterol.30690 creatinine.30700 cystatin_c.30720 direct_bilirubin.30660 gamma_glutamyltransferase.30730 glucose.30740 glycated_haemoglobin.30750 hdl_cholesterol.30760 igf-1.30770 ldl_direct.30780 lipoprotein_a.30790 oestradiol.30800 phosphate.30810 rheumatoid_factor.30820 shbg.30830 testosterone.30850 total_bilirubin.30840 total_protein.30860 triglycerides.30870 urate.30880 urea.30670 vitamin_d.30890; do
    Rscript extract.R --trait "$trait"
done
```

## vGWAS

Run C++ implementation of B-P

```sh
# more efficient to run the larger chromosomes first
for chr in  $(seq -f "%02g" 1 22); do
    for trait in alanine_aminotransferase.30620 albumin.30600 alkaline_phosphatase.30610 apolipoprotein_a.30630 apolipoprotein_b.30640 aspartate_aminotransferase.30650 c-reactive_protein.30710 calcium.30680 cholesterol.30690 creatinine.30700 cystatin_c.30720 direct_bilirubin.30660 gamma_glutamyltransferase.30730 glucose.30740 glycated_haemoglobin.30750 hdl_cholesterol.30760 igf-1.30770 ldl_direct.30780 lipoprotein_a.30790 oestradiol.30800 phosphate.30810 rheumatoid_factor.30820 shbg.30830 testosterone.30850 total_bilirubin.30840 total_protein.30860 triglycerides.30870 urate.30880 urea.30670 vitamin_d.30890; do
        sbatch run_cpp.sh "$trait" "$chr"
    done
done
```

Run R implementation of B-P

```sh
# randomly select 10000 SNPs for analysis
module load apps/bgen/1.1.6
echo -e "chromosome\tposition\tfirst_allele\talternative_alleles" > snps.txt
bgenix \
-g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen \
-incl-range 22:0- \
-list | \
awk 'NR > 2 {print $3"\t"$4"\t"$6"\t"$7}' | \
shuf | \
head -n 10000 >> data/snps.txt

# run vGWAS on subset of SNPs
sbatch run_R.sh
```

## QC

```sh
Rscript qc.R
```

## Validation

Vaidate tophits using B-F model

```sh
for chr in  $(seq -f "%02g" 1 22); do
    for trait in albumin.30600 alkaline_phosphatase.30610 apolipoprotein_a.30630 apolipoprotein_b.30640 aspartate_aminotransferase.30650 c-reactive_protein.30710 cholesterol.30690 creatinine.30700 cystatin_c.30720 direct_bilirubin.30660; do
        sbatch run_validate.sh "$trait" "$chr"
    done
done
```

## MR

Find disease outcomes associated with vQTLs

```sh
sbatch run_mr.sh "glycated_haemoglobin.30750" "ukb-d-30750_irnt"
```
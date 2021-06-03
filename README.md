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
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    Rscript extract.R --trait "$trait"
done
```

## Histogram && QQ plot of trait

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    Rscript dist.R -t "$trait"
done
```

## vGWAS

```sh
# more efficient to run the larger chromosomes first
for chr in  $(seq -f "%02g" 1 22); do
    for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
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
sbatch runR.sh Rscript bp.R -p data/ukb_bmi.txt -t body_mass_index.21001 -s data/snps.txt -o data/ukb_bmi.vgwas.r_subsample.txt
```

## QC

Map data to LDSC format (mean effect only)

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh ldsc.R -t "$trait"
done
```

Estimate LD intercept for mean effect

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch ldsc.sh data/"$trait".ldsc
done
```

Plot QQ & Manhattan

```sh
# TODO add LDSC intercept to plots
sbatch runR.sh qc-fig.R
```

Failed QC due to T1E inflation of variance test:

- alanine_aminotransferase.30620.0.0
- alkaline_phosphatase.30610.0.0
- aspartate_aminotransferase.30650.0.0
- creatinine.30700.0.0
- cystatin_c.30720.0.0
- direct_bilirubin.30660.0.0
- glycated_haemoglobin.30750.0.0
- igf_1.30770.0.0
- oestradiol.30800.0.0
- phosphate.30810.0.0
- total_bilirubin.30840.0.0
- urea.30670.0.0

## Clump vQTL

Clump tophits & filter out vQTLs without mean effect

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh clump.R -t "$trait"
done
```

Not vQTLs detected:

- albumin.30600.0.0

## GxG interaction analysis

Test each vQTL for interaction with all other vQTLs to find GxG interaction effects

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxg.R \
    -p data/"$trait".txt \
    -t "$trait" \
    -s data/"$trait".clump.txt \
    -o data/"$trait".gxg.txt
done
```

Combine GxG analyses

```sh
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxg.txt
grep -v term data/*.0.0.gxg.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxg.txt:/\t/g' >> data/gxg.txt
```

## Colocalization with eQTLs and pQTLs

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh coloc.R -t "$trait"
done
```

Combine coloc analyses

```sh
echo -e "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tgene\tregion\ttrait" > data/coloc.txt
grep -hv nsnps data/*.0.0.coloc.txt >> data/coloc.txt
```

## MR of expression on biomarker concentration

Estimate casual effect of gene product on biomaker concentration & bidirectional effect
Estimate reverse causation effect using Stiger filtering

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh mr.R -t "$trait" -o $(echo "$trait" | cut -d. -f2 | sed 's/^/ukb-d-/g' | sed 's/$/_irnt/g')
done
```

Combine MR analyses

```sh
echo -e "id.exposure\tid.outcome\toutcome\texposure\tmethod\tnsnp\tb\tse\tpval\tsnp_r2.exposure\tsnp_r2.outcome\tcorrect_causal_direction\tsteiger_pval" > data/mr.txt
cat data/*.0.0.mr.txt | grep -v ^id >> data/mr.txt
```

Report only unidirectional effects of gene product on biomarker conc where coloc evidence is in cis-region of coding gene

```sh
Rscript mr_analysis.R
```

## Subset mvQTLs with drug target loci

```sh
Rscript chembl.R
```





## GxE interaction analysis

TODO
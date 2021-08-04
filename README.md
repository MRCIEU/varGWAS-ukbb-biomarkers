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
    sbatch runR.sh extract.R --trait "$trait"
    sbatch runR.sh extract.R --trait "$trait" --log
done
```

The following have low sample size ~ 50k
- oestradiol.30800.0.0
- rheumatoid_factor.30820

## Histogram && QQ plot of trait

```sh
sbatch runR.sh dist.R
```

## vGWAS

TODO - if re-running add 20 genetic PCs
TODO - ? add snp x covariate in model

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
echo -e "chromosome\tposition\tfirst_allele\talternative_alleles" > data/alkaline_phosphatase.30610.0.0.30k_snps.txt
bgenix \
-g /mnt/storage/private/mrcieu/data/ukbiobank/genetic/variants/arrays/imputed/released/2018-09-18/data/dosage_bgen/data.chr22.bgen \
-incl-range 22:0- \
-list | \
awk 'NR > 2 {print $3"\t"$4"\t"$6"\t"$7}' | \
shuf | \
head -n 30000 >> data/alkaline_phosphatase.30610.0.0.30k_snps.txt

# run vGWAS on subset of SNPs
sbatch runR.sh \
validate_app.R \
-t alkaline_phosphatase.30610.0.0 \
-m BP

sbatch runR.sh \
validate_app.R \
-t alkaline_phosphatase.30610.0.0 \
-m JLSSC

sbatch runR.sh \
validate_app.R \
-t alkaline_phosphatase.30610.0.0 \
-m MOM

sbatch runR.sh \
validate_app.R \
-t alkaline_phosphatase.30610.0.0 \
-m GxS
```

## vGWAS QC

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

Combine results

```sh
grep Intercept data/*ldsc.log | sed 's/.ldsc.log:/\t/g' | sed 's/data\///g' | sed 's/Intercept: //g' | sed 's/ (/\t/g' | tr -d ')' > data/ldsc.txt
```

Plot QQ & Manhattan

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh qc.R -t "$trait"
done
```

Combine plots

```sh
Rscript qc-fig.R
```

## Clump vQTLs

Clump tophits & test for mean effect using heteroscedaticity robust model

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh clump.R -t "$trait"
    sbatch runR.sh clump.R -t "$trait" -p 5e-5
done
```

## Annotate vQTLs with nearest gene

bedtools closest

```sh
module load apps/bedtools/2.3.0

for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    bedtools \
    closest \
    -g /mnt/storage/home/ml18692/db/reference_genomes/released/2019-08-30/data/2.8/b37/human_g1k_v37.fasta.fai \
    -a <(awk -F"," 'NR>1 {print $1"\t"$2-1"\t"$2}' "data/""$trait"".clump.txt" | sort -k1,1n -k2,2n) \
    -b data/Homo_sapiens.GRCh37.82.sorted.bed \
    > "data/""$trait"".nearest-gene.txt"
done
```

Combine nearest gene

```sh
awk '{print $1"\t"$3"\t"$8}' data/*nearest* | uniq > data/nearest.txt
```

## Colocalize vQTLs with eQTLs and pQTLs

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh run_coloc.R -t "$trait" -o "data/""$trait""_coloc.txt" -s "data/""$trait"".clump.txt"
done
```

Combine coloc analyses

```sh
echo -e "nsnps\tPP.H0.abf\tPP.H1.abf\tPP.H2.abf\tPP.H3.abf\tPP.H4.abf\tgene\tregion\ttrait" > data/coloc.txt
grep -hv nsnps data/*.0.0_coloc.txt | grep -v ^$ >> data/coloc.txt
```

## Produce vQTL table

Combine all vQTLs

```sh
echo -n "trait," > data/vqtls.txt
head -n1 data/body_mass_index.21001.0.0.clump.txt >> data/vqtls.txt
grep -v chr data/*clump* | sed 's/data\///g' | sed 's/.clump.txt:/,/g' >> data/vqtls.txt
```

## Filter colocalization results to encoding gene cis region

Select mvQTL that colocalize with gene / protein expression in blood that are in the cis region of gene coding sequence

```sh
Rscript filter_coloc_cis.R
```

| Biomarker                          | Target        |
|------------------------------------|---------------|
| c_reactive_protein.30710.0.0       | CRP           |
| apolipoprotein_a.30630.0.0         | DDX28         |
| triglycerides.30870.0.0            | DOCK7         |
| apolipoprotein_a.30630.0.0         | DPEP3         |
| apolipoprotein_a.30630.0.0         | DUS2          |
| alanine_aminotransferase.30620.0.0 | HSD17B13      |
| shbg.30830.0.0                     | MRPL45P2      |
| vitamin_d.30890.0.0                | NADSYN1       |
| apolipoprotein_a.30630.0.0         | NFATC3        |
| hdl_cholesterol.30760.0.0          | NUP160        |
| apolipoprotein_b.30640.0.0         | PSRC1         |
| cholesterol.30690.0.0              | PSRC1         |
| ldl_direct.30780.0.0               | PSRC1         |
| hdl_cholesterol.30760.0.0          | PTPRJ         |
| vitamin_d.30890.0.0                | RP11-660L16.2 |
| hdl_cholesterol.30760.0.0          | SLC12A3       |
| apolipoprotein_b.30640.0.0         | SLC12A3       |
| cholesterol.30690.0.0              | TMEM258       |
| shbg.30830.0.0                     | ZNF554        |

- No interaction effects reported for these genes except for HSD17B13

## Filter mvQTLs at drug target loci with MR evidence for uni-directional target effect on biomarker conc

Select mvQTLs at drug target loci and perform bidirectional MR & Stieger filtering to retain targets with evidence of a mean and variance effect on biomarker concentration

```sh
Rscript filter_mvqtl_chembl.R
```

| Biomarker                 | Target|
|---------------------------|-------|
| hdl_cholesterol.30760.0.0 | NR1H3 |
| triglycerides.30870.0.0   | KIF11 |
| vitamin_d.30890.0.0       | PDE3B |

- NR1H3 agonists increase LXR-alpha activity which reduces cholesterol. rs60515486:A increases expression and reduces mean and variance of cholesterol. Not other mvQTLs detected. No mvQTLs detected at NR1H2.
- Increased KIF11 is associated with lower non-HDL and increased HDL. Existing drugs inhibit KIF11. Not progressed further.
- Increased PDE3B is associated with reduced adiposity and increased BP and reduced Vit D. Lead SNP in high LD with CYP2R1 synonymous variant which is probably the casual gene. Not progressed furhter due to lack of data.

## All targets where coloc evidence exists & MR evidence for effect of gene / protein expression on biomarker concentration

Estimate casual effect of gene product on biomaker concentration & bidirectional effect
Estimate reverse causation effect using Stiger filtering

```sh
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh run_coloc_mr.R -t "$trait" -o $(echo "$trait" | cut -d. -f2 | sed 's/^/ukb-d-/g' | sed 's/$/_irnt/g')
done
```

Combine MR analyses

```sh
echo -e "id.exposure\tid.outcome\toutcome\texposure\tmethod\tnsnp\tb\tse\tpval\tsnp_r2.exposure\tsnp_r2.outcome\tcorrect_causal_direction\tsteiger_pval" > data/mr.txt
cat data/*.0.0.mr.txt | grep -v ^id >> data/mr.txt
```

Report only unidirectional effects of gene product on biomarker conc where coloc evidence is in cis-region of coding gene

```sh
Rscript filter_coloc_mr.R
```

Results

Previously reported targets
- APOA5
- APOE
- HSD17B13

## Test for vQTL effect on biomarker using pQTLs

Extract Tier1 pQTLs from Zheng et al

```cypher
MATCH (e:Exposure)<-[r:INST_EXP { tier: 'Tier1', trans_cis:'cis' }]-(i:Instruments) RETURN r.rs_ID as rsid, r.pval_exp as pval, r.se_exp as se, r.nea as nea, r.ea as ea, r.eaf_exp as eaf, r.beta_exp as beta, r.samplesize_exp as n, r.expID as gene, r.author_exp as study, r.units_exp as units
```

```sh
Rscript pqlt.R
```

- This analysis identified effects of ATP1B2 & PLG on variance of biomarker conc. but these genes are close to coding region of biomakers so do not represent a casual effect of the target but LD

## Select mvQTLs at LDL-c drug target loci

Select instruments for LDL-c at drug target loci and test for variance effect on LDL-c & compare with RCT evidence

```sh
Rscript lipids_drug_targets.R
```

## GxG interaction analysis

Test each vQTL for interaction with all other vQTLs to find GxG interaction effects

```sh
# main analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxg.R -t "$trait"
done

# robust SE analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxg_robust.R -t "$trait"
done

# WF analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxg_wf.R -t "$trait"
done

# disease analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxg_disease.R -t "$trait"
done
```

Combine GxG analyses

```sh
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxg.txt
grep -v term data/*.0.0.gxg.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxg.txt:/\t/g' >> data/gxg.txt
echo -e "trait\tbeta\tse\tpvalue\tsample_size\tterm" > data/gxg-wf.txt
grep -v term data/*.0.0.gxg-wf.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxg-wf.txt:/\t/g' >> data/gxg-wf.txt
```

## GxE interaction analysis

```sh
# main analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxe.R -t "$trait"
done

# WF analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxe_wf.R -t "$trait"
done
```

Combine GxE analyses

```sh
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxe.txt
grep -v term data/*.0.0.gxe.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxe.txt:/\t/g' >> data/gxe.txt
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxe-log.txt
grep -v term data/*.0.0.gxe-log.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxe-log.txt:/\t/g' >> data/gxe-log.txt
echo -e "trait\tbeta\tse\tpvalue\tsample_size\tterm" > data/gxe-wf.txt
grep -v term data/*.0.0.gxe-wf.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxe-wf.txt:/\t/g' >> data/gxe-wf.txt
```
# GWAS of trait variance

All outcomes are SD normalised + log after where relevant

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

The following have low sample size
- oestradiol.30800.0.0
- rheumatoid_factor.30820

## Phenotype distribution plots

```sh
sbatch runR.sh dist.R
```

## vGWAS

```sh
# NOTE more efficient to run the larger chromosomes first
for chr in  $(seq -f "%02g" 1 22); do
    for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
        sbatch run_vargwas.sh "$trait" "$chr"
    done
done
```

## QC

Q-Q Plot & Manhattan

```sh
for trait in alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh qc.R -t "$trait"
done
```

Combine plots

```sh
Rscript qc-fig.R
```

## Clump vQTLs

Clump tophits & flip effect allele

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
    -b data/Homo_sapiens.GRCh37.104.hgnc.bed \
    > "data/""$trait"".nearest-gene.txt"
done
```

Combine nearest gene

```sh
awk '{print $1"\t"$3"\t"$7}' data/*nearest* | uniq > data/nearest.txt
```

## Top hits table

Combine all vQTLs

```sh
echo -n "trait," > data/vqtls.txt
head -n1 data/body_mass_index.21001.0.0.clump.txt >> data/vqtls.txt
grep -v chr data/*.0.0.clump.txt | sed 's/data\///g' | sed 's/.clump.txt:/,/g' >> data/vqtls.txt
sbatch runR.sh vqtls.R
```

## GxG interaction analysis

Test each vQTL for interaction with all other vQTLs to find GxG interaction effects

```sh
# main analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxg.R -t "$trait"
done

Combine GxG analyses

```sh
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxg.txt
grep -v term data/*.0.0.gxg.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxg.txt:/\t/g' >> data/gxg.txt
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxg-log.txt
grep -v term data/*.0.0.gxg-log.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxg-log.txt:/\t/g' >> data/gxg-log.txt
```

# finemap adjusted analysis
sbatch runR.sh gxg_finemap.R -t alanine_aminotransferase.30620.0.0

# qual analysis
sbatch runR.sh gxg-qual.R -t alanine_aminotransferase.30620.0.0
```

Plots

```
Rscript gxg_plot.R
Rscript gxg_qual_plot.R
```

## GxE interaction analysis

```sh
# main analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxe.R -t "$trait"
done

# qualitative analysis
for trait in body_mass_index.21001.0.0 alanine_aminotransferase.30620.0.0 albumin.30600.0.0 alkaline_phosphatase.30610.0.0 apolipoprotein_a.30630.0.0 apolipoprotein_b.30640.0.0 aspartate_aminotransferase.30650.0.0 c_reactive_protein.30710.0.0 calcium.30680.0.0 cholesterol.30690.0.0 creatinine.30700.0.0 cystatin_c.30720.0.0 direct_bilirubin.30660.0.0 gamma_glutamyltransferase.30730.0.0 glucose.30740.0.0 glycated_haemoglobin.30750.0.0 hdl_cholesterol.30760.0.0 igf_1.30770.0.0 ldl_direct.30780.0.0 lipoprotein_a.30790.0.0 oestradiol.30800.0.0 phosphate.30810.0.0 rheumatoid_factor.30820.0.0 shbg.30830.0.0 testosterone.30850.0.0 total_bilirubin.30840.0.0 total_protein.30860.0.0 triglycerides.30870.0.0 urate.30880.0.0 urea.30670.0.0 vitamin_d.30890.0.0; do
    sbatch runR.sh gxe-qual.R -t "$trait"
done
```

Combine GxE analyses

```sh
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxe.txt
grep -v term data/*.0.0.gxe.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxe.txt:/\t/g' >> data/gxe.txt
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxe-main.txt
grep -v term data/*.0.0.gxe-main.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxe-main.txt:/\t/g' >> data/gxe-main.txt
echo -e "trait\tterm\testimate\tstd.error\tstatistic\tp.value" > data/gxe-log.txt
grep -v term data/*.0.0.gxe-log.txt | grep -v :$ | sed 's/data\///g' | sed 's/.gxe-log.txt:/\t/g' >> data/gxe-log.txt
head -n1 data/alanine_aminotransferase.30620.0.0.gxe-qual.txt > data/gxe-qual.txt
cat data/*gxe-qual.txt | grep -v ^term >> data/gxe-qual.txt
```

Plots

```
Rscript gxe_plot.R
Rscript gxe_qual_plot.R
```
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
for chr in  $(seq -f "%02g" 1 22); do
    sbatch run_cpp.sh "$chr"
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
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
# run only chr22
sbatch run_R.sh 22
```
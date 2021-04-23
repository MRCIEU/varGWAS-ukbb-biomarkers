# script to extract rare LoF variants in SEC31B from UK Biobank
set -euo pipefail
module load apps/plink-1.90

# extract exome calls SEC31B (+/- 100kb)
echo "##fileformat=VCFv4.2" > ~/ukbb.SEC31B.exome.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ~/ukbb.SEC31B.exome.vcf
awk '$4 >= 100486546 && $4 <= 100519964 {print "chr"$1"\t"$4"\t"$2"\t"$6"\t"$5"\t.\tPASS\t."}' \
/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_snp_chr10.bim \
>> ~/ukbb.SEC31B.exome.vcf

# annotate calls using VEP
/newhome/ml18692/apps/ensembl-vep/vep \
-i ~/ukbb.SEC31B.exome.vcf \
--cache \
--force_overwrite \
--sift p \
--polyphen p \
--tab \
--no_stats \
-o ~/ukbb.SEC31B.exome.vep.txt

# extract exome calls to text file
plink \
--bim /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_snp_chr10.bim \
--fam /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/raw/exome_download/data.exome.ukbapp.15825.fam \
--bed /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_cal_chr10_b0_v1.bed \
--extract <(grep -v '#' ~/ukbb.SEC31B.exome.vep.txt | cut -s -f1 | sort -u) \
--recode A \
--out ~/ukbb.SEC31B.exome.txt

# copy to BC4
scp ~/ukbb.SEC31B.exome.txt.raw bc4:/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data
scp ~/ukbb.SEC31B.exome.vep.txt bc4:/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data

# clean up
rm ~/ukbb.SEC31B.exome.txt.*
rm ~/ukbb.SEC31B.exome.vcf
rm ~/ukbb.SEC31B.exome.vep.txt
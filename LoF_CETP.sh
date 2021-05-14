# script to extract rare LoF variants in CETP from UK Biobank
set -euo pipefail
module load apps/plink-1.90

# extract exome calls CETP (+/- 100kb)
echo "##fileformat=VCFv4.2" > ~/ukbb.CETP.exome.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ~/ukbb.CETP.exome.vcf
awk '$4 >= 56961923 && $4 <= 56983845 {print "chr"$1"\t"$4"\t"$2"\t"$6"\t"$5"\t.\tPASS\t."}' \
/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_snp_chr16.bim \
>> ~/ukbb.CETP.exome.vcf

# annotate calls using VEP
/newhome/ml18692/apps/ensembl-vep/vep \
-i ~/ukbb.CETP.exome.vcf \
--cache \
--force_overwrite \
--sift p \
--polyphen p \
--tab \
--no_stats \
-o ~/ukbb.CETP.exome.vep.txt

# extract exome calls to text file
plink \
--bim /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_snp_chr16.bim \
--fam /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/raw/exome_download/data.exome.ukbapp.15825.fam \
--bed /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_cal_chr16_b0_v1.bed \
--extract <(grep -v '#' ~/ukbb.CETP.exome.vep.txt | cut -s -f1 | sort -u) \
--recode A \
--out ~/ukbb.CETP.exome.txt

# copy to BC4
scp ~/ukbb.CETP.exome.txt.raw bc4:/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data
scp ~/ukbb.CETP.exome.vep.txt bc4:/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data

# clean up
rm ~/ukbb.CETP.exome.txt.*
rm ~/ukbb.CETP.exome.vcf
rm ~/ukbb.CETP.exome.vep.txt
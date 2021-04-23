# script to extract rare LoF variants in HSD17B13 from UK Biobank
set -euo pipefail
module load apps/plink-1.90

# extract exome calls HSD17B13 (+/- 100kb)
echo "##fileformat=VCFv4.2" > ~/ukbb.HSD17B13.exome.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" >> ~/ukbb.HSD17B13.exome.vcf
awk '$4 >= 87303689 && $4 <= 87322986 {print "chr"$1"\t"$4"\t"$2"\t"$6"\t"$5"\t.\tPASS\t."}' \
/projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_snp_chr4.bim \
>> ~/ukbb.HSD17B13.exome.vcf

# annotate calls using VEP
/newhome/ml18692/apps/ensembl-vep/vep \
-i ~/ukbb.HSD17B13.exome.vcf \
--cache \
--force_overwrite \
--sift p \
--polyphen p \
--tab \
--no_stats \
-o ~/ukbb.HSD17B13.exome.vep.txt

# extract exome calls to text file
plink \
--bim /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_snp_chr4.bim \
--fam /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/raw/exome_download/data.exome.ukbapp.15825.fam \
--bed /projects/MRC-IEU/research/data/ukbiobank/genetic/variants/exome/dev/release_candidate/data/raw_downloaded/ukb_cal_chr4_b0_v1.bed \
--extract <(grep -v '#' ~/ukbb.HSD17B13.exome.vep.txt | cut -s -f1 | sort -u) \
--recode A \
--out ~/ukbb.HSD17B13.exome.txt

# copy to BC4
scp ~/ukbb.HSD17B13.exome.txt.raw bc4:/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data
scp ~/ukbb.HSD17B13.exome.vep.txt bc4:/mnt/storage/home/ml18692/projects/jlst-cpp-vgwas/data

# clean up
rm ~/ukbb.HSD17B13.exome.txt.*
rm ~/ukbb.HSD17B13.exome.vcf
rm ~/ukbb.HSD17B13.exome.vep.txt
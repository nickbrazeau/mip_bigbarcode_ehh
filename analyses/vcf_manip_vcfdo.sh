#!/usr/bin/env bash
PATH="/Users/NFB/Documents/GitHub/mip_bigbarcode_ehh/data" # user can define path 
FASTA="/Users/NFB/Documents/MountPoints/mountIDEEL/julianog/users/apm/p_reichenowi/ancestral.fa" # apm downloaded from p_reich MS and imputed alleles
gunzip -c $PATH/mipbi_bigbarcodepanel.vcf.gz | bgzip >  $PATH/mipbi_bigbarcodepanel.vcf.bgz #written out as gzip compress, need bgzip
gunzip -c $PATH/mipbi_drugrespanel.vcf.gz | bgzip > $PATH/mipbi_drugrespanel.vcf.bgz
bcftools index $PATH/mipbi_bigbarcodepanel.vcf.bgz
bcftools index $PATH/mipbi_drugrespanel.vcf.bgz
bcftools concat $PATH/mipbi_bigbarcodepanel.vcf $PATH/mipbi_drugrespanel.vcf  --allow-overlaps | \
bcftools sort | \
vcfdo polarize -f $FASTA | bgzip >  $PATH/polarized_mipbi_bigbarcode.vcf.gz

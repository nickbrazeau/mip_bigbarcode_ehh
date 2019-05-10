#!/usr/bin/env bash

path=~/Documents/GitHub/mip_bigbarcode_ehh/data/derived/vcfs/ # user can define path
FASTA=~/Documents/MountPoints/mountIDEEL/julianog/users/apm/p_reichenowi/ancestral.fa # apm downloaded from p_reich MS and imputed alleles

cd $path

gunzip -c mipbi_bigbarcodepanel.vcf.gz | bgzip >  mipbi_bigbarcodepanel.vcf.bgz #written out as gzip compress, need bgzip
gunzip -c mipbi_drugrespanel.vcf.gz | bgzip > mipbi_drugrespanel.vcf.bgz
bcftools index mipbi_bigbarcodepanel.vcf.bgz
bcftools index mipbi_drugrespanel.vcf.bgz
bcftools concat mipbi_bigbarcodepanel.vcf.bgz mipbi_drugrespanel.vcf.bgz  --allow-overlaps | \
bcftools sort | \
vcfdo polarize -f $FASTA | bgzip >  polarized_mipbi_drugres_bigbarcode_panels.vcf.bgz

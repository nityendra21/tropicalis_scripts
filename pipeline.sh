#!/usr/bin/env bash
set -e
source config.sh

mkdir -p ${OUTDIR}/{fastp,bam,gvcf,vcf,logs}

bash scripts/02_bwa_index.sh

while read SAMPLE; do
    bash scripts/01_fastp.sh ${SAMPLE}
    bash scripts/03_bwa_align.sh ${SAMPLE}
    bash scripts/04_markduplicates.sh ${SAMPLE}
    bash scripts/05_haplotypecaller.sh ${SAMPLE}
done < ${SAMPLES}

bash scripts/06_genotypegvcfs.sh
bash scripts/07_variant_filtering.sh

echo "Pipeline completed successfully."

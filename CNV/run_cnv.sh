#!/usr/bin/env bash
set -e
source config.sh

mkdir -p ${CNV_OUTDIR}/{depth,windows,normalized,calls}

while read SAMPLE; do
    bash depth_calculation.sh ${SAMPLE}
    bash window_average.sh ${SAMPLE}
    bash normalize.sh ${SAMPLE}
    bash call_cnv.sh ${SAMPLE}
done < ../${SAMPLES}

echo "CNV calling completed."

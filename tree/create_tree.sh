#!/usr/bin/env bash
set -e
mkdir -p ${PHYLO_OUTDIR}

PHYLIP=$1
BASENAME=$(basename ${PHYLIP} .phy)

${IQTREE} \
  -s ${PHYLIP} \
  -m "GTR+G" -mrate CAT \
  -bb 1000 \
  -nt 40 --mem 300G \
  -st DNA \
  -pre ${PHYLO_OUTDIR}/${BASENAME} \
  -safe

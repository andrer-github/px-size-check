#!/usr/bin/env bash
# Download the EM map (EMD-17992) and the E. coli DPS X-ray models used in Dickerson et al. 2024.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p data models

curl -L -o data/emd_17992.map.gz \
  https://ftp.ebi.ac.uk/pub/databases/emdb/structures/EMD-17992/map/emd_17992.map.gz
gunzip -f data/emd_17992.map.gz

for id in 1DPS 1JRE 1JTS 1L8H 1L8I 6GCM 7AQS 8OUC; do
  curl -L -o "models/${id}.cif" "https://files.rcsb.org/download/${id}-assembly1.cif"
done

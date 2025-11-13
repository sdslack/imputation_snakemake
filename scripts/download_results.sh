#!/bin/bash

# TO NOTE: Need to use different host paths to retreive job/file information
# then to make downloadable links.

# Parse args
while getopts "i:c:o:j:" opt; do
  case $opt in
    i) imp="$OPTARG" ;;
    c) code_dir="$OPTARG" ;;
    o) out_dir="$OPTARG" ;;
    j) job_id="$OPTARG" ;;
    \?) usage ;;
  esac
done

# Setup
cd "$code_dir"

# Based on server, set host name and get API key
host_i=""
host_d=""
key=""
if [ "$imp" = "topmed" ]; then
  host_i="https://imputation.biodatacatalyst.nhlbi.nih.gov/api/v2"
  host_d="https://imputation.biodatacatalyst.nhlbi.nih.gov"
  key=$(python -c "import config.key as k; print(k.TOPMED_API)")

else
  host_i="https://imputationserver.sph.umich.edu/api/v2"
  host_d="https://imputationserver.sph.umich.edu"
  key=$(python -c "import config.key as k; print(k.MICH_API)")
fi

# Get initial json with file information
curl -s -H "X-Auth-Token: $key" \
  "${host_i}/jobs/${job_id}" > ${out_dir}/job_metadata.json

# Use json to make downloadable links
  #TODO: this pattern works for TOPMed, need to confirm for Mich!
jq -r --arg host "$host_d" '.outputParams[].files[] | "\($host)/share/results/\(.hash)/\(.name)"' \
  ${out_dir}/job_metadata.json > ${out_dir}/job_download_links.txt

# Download with aria2 (recommended by server)
# Raise --min-split-size to avoid hitting imputation server download limit
aria2c \
  -j 4 -x 4 --split 4 --min-split-size=1G --continue=true \
  -i "${out_dir}/job_download_links.txt" -d "${out_dir}"

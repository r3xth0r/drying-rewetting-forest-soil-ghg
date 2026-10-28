#!/usr/bin/env bash

base_url="https://public.hub.geosphere.at/datahub/resources/winfore-v2-1d-1km/filelisting/SPEI30"
output_dir="dat/raw"

mkdir -p "$output_dir"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting download of SPEI data..."
echo "» See https://data.hub.geosphere.at/dataset/winfore-v2-1d-1km for details"

for year in {2016..2025}; do
    remote_file="WINFORE2_SPEI30_${year}.nc"

    local_file="${remote_file,,}"
    local_file="${local_file//_/-}"

    output_path="${output_dir}/${local_file}"

    if [[ -s "$output_path" ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Skipping existing file: $output_path"
        continue
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Downloading $local_file..."

    curl --fail --location --retry 3 \
         --output "$output_path" \
         "${base_url}/${remote_file}"

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Finished downloading $local_file."
done


echo "[$(date '+%Y-%m-%d %H:%M:%S')] SPEI data download completed \\o/"
echo "    » Downloaded files are located in: $output_dir"

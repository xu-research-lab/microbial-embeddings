#!/usr/bin/env bash
# Fetch the genome_mapping data files that exceed GitHub's 100 MB limit from the
# GitHub Release and place them in ./data. Safe to re-run: finished files are skipped.
set -euo pipefail

TAG="genome-mapping-data-v1"
BASE="https://github.com/xu-research-lab/microbial-embeddings/releases/download/${TAG}"

cd "$(dirname "$0")/data"

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }

# asset  sha256  final file in data/
while read -r asset sum target; do
  if [[ -f "$target" ]]; then echo "skip $target (exists)"; continue; fi
  if [[ -f "$asset" ]]; then
    echo "use local $asset"; mv "$asset" "$asset.part"
  else
    echo "download $asset"
    curl -fL --retry 3 --connect-timeout 30 -o "$asset.part" "$BASE/$asset" || {
      rm -f "$asset.part"
      echo "Download failed (no access to github.com?). Get $asset from" >&2
      echo "  https://github.com/xu-research-lab/microbial-embeddings/releases/tag/${TAG}" >&2
      echo "and place it in $(pwd), then re-run this script." >&2
      exit 1
    }
  fi
  [[ "$(sha256 "$asset.part")" == "$sum" ]] || { echo "checksum mismatch: $asset" >&2; rm -f "$asset.part"; exit 1; }
  if [[ "$target" == "$asset" ]]; then
    mv "$asset.part" "$target"
  else
    gzip -dc "$asset.part" > "$target.part"
    mv "$target.part" "$target"
    rm "$asset.part"
  fi
done <<'LIST'
bac120_metadata_r220.tsv.gz f133b46c54a7b7df3ee09e92d1930d84bb9fa987baa55c817dc509ccde49fd81 bac120_metadata_r220.tsv.gz
barrnap.fna.gz 10a1138e3d888b2e58ce05e7408f7be1ccf4a82ccb76add0d6bde85be870cb5e barrnap.fna
LIST
echo "done"

#!/usr/bin/env bash
# Fetch the function_phylogeny_hgt inputs that exceed GitHub's 100 MB limit from
# the GitHub Release and place them under ./data. Safe to re-run: finished files
# are skipped.
set -euo pipefail

TAG="function-phylogeny-hgt-data-v1"
BASE="https://github.com/xu-research-lab/microbial-embeddings/releases/download/${TAG}"

cd "$(dirname "$0")/data"

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }

# asset  sha256  final file under data/
while read -r asset sum target; do
  if [[ -f "$target" ]]; then echo "skip $target (exists)"; continue; fi
  mkdir -p "$(dirname "$target")"
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
  gzip -dc "$asset.part" > "$target.part"
  mv "$target.part" "$target"
  rm "$asset.part"
done <<'LIST'
table.co.gz a67a8e4b9ba81e6d69a72fbb9ef4b94572f8485e5711b79e4c7da8c984c43b75 cooccur_otuembedding/table.co
bac_EC_predicted.tsv.gz f9f8ea7bb1ab8186aa02ae30a868280b680687e5d68fd5bbb7a97a4db23e02b4 picrust/bac_EC_predicted.tsv
bac_KO_predicted.tsv.gz 65e2ec127b71bbf5a8bbf89c2301ae391e6a6221accd750843ac32724687e842 picrust/bac_KO_predicted.tsv
LIST
echo "done"

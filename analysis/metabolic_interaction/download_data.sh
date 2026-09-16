#!/usr/bin/env bash
# Fetch the metabolic_interaction inputs that exceed GitHub's file size limits
# from the GitHub Release and unpack them under ./data. Safe to re-run: finished
# files are skipped. An asset already placed in ./data by hand is used as is and
# kept; downloaded assets are removed once unpacked.
set -euo pipefail

TAG="metabolic-interaction-data-v1"
BASE="https://github.com/xu-research-lab/microbial-embeddings/releases/download/${TAG}"

cd "$(dirname "$0")/data"

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }

# fetch ASSET SHA256: make sure a verified ASSET is in data/; prints "downloaded" if we fetched it
fetch() {
  local asset=$1 sum=$2
  if [[ -f "$asset" ]]; then
    echo "use local $asset" >&2
    [[ "$(sha256 "$asset")" == "$sum" ]] || { echo "checksum mismatch: local $asset" >&2; exit 1; }
    return
  fi
  echo "download $asset" >&2
  curl -fL --retry 3 --connect-timeout 30 -o "$asset.part" "$BASE/$asset" || {
    rm -f "$asset.part"
    echo "Download failed (no access to github.com?). Get $asset from" >&2
    echo "  https://github.com/xu-research-lab/microbial-embeddings/releases/tag/${TAG}" >&2
    echo "and place it in $(pwd), then re-run this script." >&2
    exit 1
  }
  [[ "$(sha256 "$asset.part")" == "$sum" ]] || { echo "checksum mismatch: $asset" >&2; rm -f "$asset.part"; exit 1; }
  mv "$asset.part" "$asset"
  echo downloaded
}

# Gzipped tables -> data/<table>
while read -r asset sum target; do
  if [[ -f "$target" ]]; then echo "skip $target (exists)"; continue; fi
  got=$(fetch "$asset" "$sum")
  gzip -dc "$asset" > "$target.part"
  mv "$target.part" "$target"
  [[ "$got" == downloaded ]] && rm "$asset"
done <<'LIST'
mapping_bigg_gene_table.tsv.gz 252612c56727a8191eddeb2bf6ab94a6ec97b3af44e7529e21c39874b41ccbfa mapping_bigg_gene_table.tsv
mapping_scores_table.tsv.gz 8a9bc663e46d55e4e3839db80883ab3b1020b57d9c339555fdb9e716d2a532d0 mapping_scores_table.tsv
bigg_gene_predicted.tsv.gz 27e1ddb08e32786563b62bb58a70ea47e415ad931b16b720a9c43d5401453245 bigg_gene_predicted.tsv
scores_predicted.tsv.gz 75d58b924945961d9d42a701193b570257853147285a90968de0197ff2e77a31 scores_predicted.tsv
LIST

# Tarballs -> data/<dir>/ ; the archive either comes whole or as split parts (asset:sha256 ...)
while read -r dir sum parts; do
  if [[ -d "$dir" ]]; then echo "skip $dir/ (exists)"; continue; fi
  tgz="$dir.tar.gz"
  if [[ -f "$tgz" || -z "$parts" ]]; then
    got=$(fetch "$tgz" "$sum")
  else
    got=downloaded   # the joined archive is ours to remove
    names=() pdl=()
    for p in $parts; do
      pgot=$(fetch "${p%%:*}" "${p#*:}")
      names+=("${p%%:*}")
      [[ "$pgot" == downloaded ]] && pdl+=("${p%%:*}")
    done
    cat "${names[@]}" > "$tgz.part"
    [[ "$(sha256 "$tgz.part")" == "$sum" ]] || { echo "checksum mismatch: joined $tgz" >&2; rm -f "$tgz.part"; exit 1; }
    mv "$tgz.part" "$tgz"
    [[ ${#pdl[@]} -gt 0 ]] && rm -f "${pdl[@]}"
  fi
  echo "unpack $tgz" >&2
  rm -rf "$dir.unpack"; mkdir "$dir.unpack"
  tar -xzf "$tgz" -C "$dir.unpack"
  mv "$dir.unpack/$dir" "$dir"
  rmdir "$dir.unpack"
  [[ "$got" == downloaded ]] && rm "$tgz"
done <<'LIST'
blast_output_bigg 2089c574435c6028ae8d57db30ded5640be85fcd946137f2f0f2dbf18bf18b32
OTU_metabolic_model_M3 fe058a7b8935b21874854f4bdd54c03b2f40538443cc5c70cae87c5026a13178 OTU_metabolic_model_M3.tar.gz.part00:8fcc812f0d11467cb6645f1caab99120274748d6c02b46f81a30e07972c8a40f OTU_metabolic_model_M3.tar.gz.part01:b5fdb366efdc3def1f0d84b352459570a4f7fc14f0ce2cab9b1bfa27e65edf14
LIST
echo "done"

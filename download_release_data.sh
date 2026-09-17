#!/usr/bin/env bash
# Fetch the analysis inputs that exceed GitHub's file size limits from the
# GitHub Releases and put each one where its analysis expects it.
#
#   bash download_release_data.sh                         # all groups (~5.7 GB download)
#   bash download_release_data.sh metabolic_interaction   # only the listed groups
#
# Groups: genome_mapping, function_phylogeny_hgt, metabolic_interaction.
# Safe to re-run: targets that already exist are skipped. An asset placed by hand
# in its target's directory is used as is and kept; downloaded assets are
# removed once unpacked.
set -euo pipefail

REPO=https://github.com/xu-research-lab/microbial-embeddings
GROUPS_ALL=(genome_mapping function_phylogeny_hgt metabolic_interaction)

cd "$(dirname "$0")"
ROOT=$(pwd)

for g in "$@"; do
  [[ " ${GROUPS_ALL[*]} " == *" $g "* ]] || { echo "unknown group: $g (choose from: ${GROUPS_ALL[*]})" >&2; exit 2; }
done
SEL="$*"
selected() { [[ -z "$SEL" || " $SEL " == *" $1 "* ]]; }

sha256() { if command -v sha256sum >/dev/null; then sha256sum "$1"; else shasum -a 256 "$1"; fi | cut -d' ' -f1; }

# fetch TAG ASSET SHA256: make sure a verified ASSET is in the current directory;
# prints "downloaded" if we fetched it (so the caller may remove it afterwards)
fetch() {
  local tag=$1 asset=$2 sum=$3
  if [[ -f "$asset" ]]; then
    echo "use local $asset" >&2
    [[ "$(sha256 "$asset")" == "$sum" ]] || { echo "checksum mismatch: local $asset" >&2; exit 1; }
    return
  fi
  echo "download $asset ($tag)" >&2
  curl -#fL --retry 3 --connect-timeout 30 -o "$asset.part" "$REPO/releases/download/$tag/$asset" || {
    rm -f "$asset.part"
    echo "Download failed (no access to github.com?). Get $asset from" >&2
    echo "  $REPO/releases/tag/$tag" >&2
    echo "and place it in $(pwd), then re-run this script." >&2
    exit 1
  }
  [[ "$(sha256 "$asset.part")" == "$sum" ]] || { echo "checksum mismatch: $asset" >&2; rm -f "$asset.part"; exit 1; }
  mv "$asset.part" "$asset"
  echo downloaded
}

# group | release tag | asset | sha256 | target, relative to the repository root
#   target == asset name     -> kept as downloaded
#   *.gz                     -> gunzipped to target
#   *.tar.gz, target ends /  -> unpacked to that directory
#   optional 6th column      -> the asset ships as split parts (name:sha256,...)
while read -r group tag asset sum target parts; do
  selected "$group" || continue
  dir=$(dirname "$target"); name=$(basename "$target")
  mkdir -p "$ROOT/$dir"
  cd "$ROOT/$dir"
  if [[ -e "$name" ]]; then echo "skip $target (exists)"; continue; fi

  if [[ -f "$asset" || -z "${parts:-}" ]]; then
    got=$(fetch "$tag" "$asset" "$sum")
  else
    got=downloaded   # the joined archive is ours to remove
    names=() pdl=()
    for p in ${parts//,/ }; do
      pgot=$(fetch "$tag" "${p%%:*}" "${p#*:}")
      names+=("${p%%:*}")
      [[ "$pgot" == downloaded ]] && pdl+=("${p%%:*}")
    done
    cat "${names[@]}" > "$asset.part"
    [[ "$(sha256 "$asset.part")" == "$sum" ]] || { echo "checksum mismatch: joined $asset" >&2; rm -f "$asset.part"; exit 1; }
    mv "$asset.part" "$asset"
    [[ ${#pdl[@]} -gt 0 ]] && rm -f "${pdl[@]}"
  fi

  if [[ "$asset" == "$name" ]]; then
    echo "ready $target"; continue
  elif [[ "$target" == */ ]]; then
    echo "unpack $asset" >&2
    rm -rf "$name.unpack"; mkdir "$name.unpack"
    tar -xzf "$asset" -C "$name.unpack"
    mv "$name.unpack/$name" "$name"
    rmdir "$name.unpack"
  else
    gzip -dc "$asset" > "$name.part"
    mv "$name.part" "$name"
  fi
  [[ "$got" == downloaded ]] && rm "$asset"
  echo "ready $target"
done <<'LIST'
genome_mapping genome-mapping-data-v1 bac120_metadata_r220.tsv.gz f133b46c54a7b7df3ee09e92d1930d84bb9fa987baa55c817dc509ccde49fd81 analysis/resources/genome_mapping/data/bac120_metadata_r220.tsv.gz
genome_mapping genome-mapping-data-v1 barrnap.fna.gz 10a1138e3d888b2e58ce05e7408f7be1ccf4a82ccb76add0d6bde85be870cb5e analysis/resources/genome_mapping/data/barrnap.fna
function_phylogeny_hgt function-phylogeny-hgt-data-v1 table.co.gz a67a8e4b9ba81e6d69a72fbb9ef4b94572f8485e5711b79e4c7da8c984c43b75 analysis/function_phylogeny_hgt/data/cooccur_otuembedding/table.co
function_phylogeny_hgt function-phylogeny-hgt-data-v1 bac_EC_predicted.tsv.gz f9f8ea7bb1ab8186aa02ae30a868280b680687e5d68fd5bbb7a97a4db23e02b4 analysis/function_phylogeny_hgt/data/picrust/bac_EC_predicted.tsv
function_phylogeny_hgt function-phylogeny-hgt-data-v1 bac_KO_predicted.tsv.gz 65e2ec127b71bbf5a8bbf89c2301ae391e6a6221accd750843ac32724687e842 analysis/function_phylogeny_hgt/data/picrust/bac_KO_predicted.tsv
metabolic_interaction metabolic-interaction-data-v1 mapping_bigg_gene_table.tsv.gz 252612c56727a8191eddeb2bf6ab94a6ec97b3af44e7529e21c39874b41ccbfa analysis/metabolic_interaction/data/mapping_bigg_gene_table.tsv
metabolic_interaction metabolic-interaction-data-v1 mapping_scores_table.tsv.gz 8a9bc663e46d55e4e3839db80883ab3b1020b57d9c339555fdb9e716d2a532d0 analysis/metabolic_interaction/data/mapping_scores_table.tsv
metabolic_interaction metabolic-interaction-data-v1 bigg_gene_predicted.tsv.gz 27e1ddb08e32786563b62bb58a70ea47e415ad931b16b720a9c43d5401453245 analysis/metabolic_interaction/data/bigg_gene_predicted.tsv
metabolic_interaction metabolic-interaction-data-v1 scores_predicted.tsv.gz 75d58b924945961d9d42a701193b570257853147285a90968de0197ff2e77a31 analysis/metabolic_interaction/data/scores_predicted.tsv
metabolic_interaction metabolic-interaction-data-v1 blast_output_bigg.tar.gz 2089c574435c6028ae8d57db30ded5640be85fcd946137f2f0f2dbf18bf18b32 analysis/metabolic_interaction/data/blast_output_bigg/
metabolic_interaction metabolic-interaction-data-v1 OTU_metabolic_model_M3.tar.gz fe058a7b8935b21874854f4bdd54c03b2f40538443cc5c70cae87c5026a13178 analysis/metabolic_interaction/data/OTU_metabolic_model_M3/ OTU_metabolic_model_M3.tar.gz.part00:8fcc812f0d11467cb6645f1caab99120274748d6c02b46f81a30e07972c8a40f,OTU_metabolic_model_M3.tar.gz.part01:b5fdb366efdc3def1f0d84b352459570a4f7fc14f0ce2cab9b1bfa27e65edf14
LIST
echo "done"

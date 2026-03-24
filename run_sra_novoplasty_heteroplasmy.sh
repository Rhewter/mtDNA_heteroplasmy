#!/usr/bin/env bash
# run_sra_novoplasty_heteroplasmy.sh
# Pipeline:
# 1. SRA download (prefetch + fasterq-dump)
# 2. Trimmomatic QC slidingwindow 4:15
# 3. NOVOPlasty mtDNA assembly
# 4. heteroplasmy analysis using bwa/samtools/bcftools

set -euo pipefail
IFS=$'\n\t'

usage() {
  cat <<EOF
Usage: $0 -i sra_ids.txt [-o outdir] [-t threads] [--trim-params "SLIDINGWINDOW:4:15 MINLEN:36"] [--step all|download|trim|novo|hetero] [--novoplasty-config path]

Required:
  -i FILE             list of SRA run IDs (one per line, e.g., SRR31925970)

Options:
  -o DIR              output base directory (default: ./mtDNA_pipeline_out)
  -t N                threads (default: 8)
  -n PATH             path to Trimmomatic jar (default: trimmomatic in PATH)
  -a PATH             adapter fasta (default: adapters/TruSeq3-PE.fa or packaged)
  -p PATH             path to NOVOPlasty.pl script (default: NOVOPlasty.pl in PATH)
  -r PATH             mitochondrial reference fasta for mapping/heteroplasmy (default: auto from assembly)
  --trim-params STR   custom trimming params (default: SLIDINGWINDOW:4:15 MINLEN:36)
  --step NAME         run from one step to end: download, trim, novo, hetero, all (default: all)
  -h|--help           show this message
EOF
  exit 1
}

# command line parse
POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    -i) SRA_LIST="$2"; shift; shift ;; 
    -o) OUTDIR="$2"; shift; shift ;; 
    -t) THREADS="$2"; shift; shift ;; 
    -n) TRIMMOMATIC_JAR="$2"; shift; shift ;; 
    -a) ADAPTERS="$2"; shift; shift ;; 
    -p) NOVOPLASTY_SCRIPT="$2"; shift; shift ;; 
    -r) REF_FA="$2"; shift; shift ;; 
    --trim-params) TRIM_PARAMS="$2"; shift; shift ;; 
    --step) STEP="$2"; shift; shift ;; 
    -h|--help) usage ;; 
    *) POSITIONAL+=("$1"); shift ;; 
  esac
done
set -- "${POSITIONAL[@]}"

# defaults (set if not provided)
OUTDIR="${OUTDIR:-./mtDNA_pipeline_out}"
THREADS="${THREADS:-8}"
TRIM_PARAMS="${TRIM_PARAMS:-SLIDINGWINDOW:4:15 MINLEN:36}"
STEP="${STEP:-all}"
NOVOPLASTY_SCRIPT="${NOVOPLASTY_SCRIPT:-NOVOPlasty.pl}"
TRIMMOMATIC_JAR="${TRIMMOMATIC_JAR:-trimmomatic}"
ADAPTERS="${ADAPTERS:-TruSeq3-PE.fa}"
REF_FA="${REF_FA:-}"

# function to detect tools automatically
detect_tools() {
  # Detect Trimmomatic jar if not provided or not in PATH
  if [[ "$TRIMOMATIC_JAR" == "trimmomatic" ]] || ! command -v "$TRIMOMATIC_JAR" &>/dev/null; then
    # Try common locations for trimmomatic jar
    for path in /usr/share/java/trimmomatic*.jar /opt/trimmomatic/trimmomatic*.jar /usr/local/share/java/trimmomatic*.jar "$HOME"/trimmomatic/trimmomatic*.jar; do
      if [[ -f "$path" ]]; then
        TRIMOMATIC_JAR="$path"
        echo "Detected Trimmomatic jar at: $TRIMOMATIC_JAR"
        break
      fi
    done
  fi

  # Detect NOVOPlasty.pl if not provided or not in PATH
  if [[ "$NOVOPLASTY_SCRIPT" == "NOVOPlasty.pl" ]] || ! command -v "$NOVOPLASTY_SCRIPT" &>/dev/null; then
    # Try current directory and subdirectories
    for path in NOVOPlasty.pl ./*/NOVOPlasty.pl ./*/*/NOVOPlasty.pl; do
      if [[ -f "$path" ]]; then
        NOVOPLASTY_SCRIPT="$path"
        echo "Detected NOVOPlasty script at: $NOVOPLASTY_SCRIPT"
        break
      fi
    done
  fi

  # Detect adapters if not provided
  if [[ ! -f "$ADAPTERS" ]]; then
    # Try common locations
    for path in /usr/share/trimmomatic/adapters/TruSeq3-PE.fa /opt/trimmomatic/adapters/TruSeq3-PE.fa "$HOME"/adapters/TruSeq3-PE.fa; do
      if [[ -f "$path" ]]; then
        ADAPTERS="$path"
        echo "Detected adapters at: $ADAPTERS"
        break
      fi
    done
  fi
}

if [[ -z "${SRA_LIST:-}" ]]; then
  echo "ERROR: -i sra_ids.txt is required" >&2; usage
fi

# detect tools automatically
detect_tools

# check tools
for cmd in prefetch fasterq-dump bwa samtools bcftools; do
  if ! command -v "$cmd" &>/dev/null; then
    echo "ERROR: required command '$cmd' not found in PATH" >&2
    exit 1
  fi
done

if ! command -v "$TRIMOMATIC_JAR" &>/dev/null && [[ ! -f "$TRIMOMATIC_JAR" ]]; then
  echo "ERROR: Trimmomatic not found at '$TRIMOMATIC_JAR'" >&2; exit 1
fi

if ! command -v "$NOVOPLASTY_SCRIPT" &>/dev/null && [[ ! -f "$NOVOPLASTY_SCRIPT" ]]; then
  echo "ERROR: NOVOPlasty script not found at '$NOVOPLASTY_SCRIPT'" >&2; exit 1
fi

mkdir -p "$OUTDIR"
raw_dir="$OUTDIR/01_raw"
trim_dir="$OUTDIR/02_trimmed"
novo_dir="$OUTDIR/03_novoplasty"
hetero_dir="$OUTDIR/04_heteroplasmy"
mkdir -p "$raw_dir" "$trim_dir" "$novo_dir" "$hetero_dir"

# adapter file handling: prefer given path else installed local
if [[ ! -f "$ADAPTERS" ]]; then
  # default from Trimmomatic package path heuristics
  if [[ -f "/usr/share/trimmomatic/$ADAPTERS" ]]; then
    ADAPTERS="/usr/share/trimmomatic/$ADAPTERS"
  elif [[ -f "${HOME}/adapters/$ADAPTERS" ]]; then
    ADAPTERS="${HOME}/adapters/$ADAPTERS"
  else
    echo "ERROR: adapter file $ADAPTERS not found" >&2
    exit 1
  fi
fi

if [[ -n "$REF_FA" && ! -f "$REF_FA" ]]; then
  echo "ERROR: reference fasta '$REF_FA' not found" >&2
  exit 1
fi

# load ids
mapfile -t sra_ids < "$SRA_LIST"

run_download() {
  local sid="$1"
  local sid_dir="$raw_dir/$sid"
  mkdir -p "$sid_dir"
  local done_flag="$sid_dir/download.done"
  if [[ -f "$done_flag" ]]; then
    echo "[download] $sid already done, skipping"
    return
  fi

  echo "[download] $sid: prefetch"
  prefetch "$sid" --output-directory "$sid_dir"

  echo "[download] $sid: fasterq-dump"
  fasterq-dump "$sid" --split-files --outdir "$sid_dir" --threads "$THREADS"

  touch "$done_flag"
}

run_trim() {
  local sid="$1"
  local sid_dir="$raw_dir/$sid"
  local done_flag="$trim_dir/$sid.trim.done"
  if [[ -f "$done_flag" ]]; then
    echo "[trim] $sid already done, skipping"
    return
  fi

  local r1="$sid_dir/${sid}_1.fastq"
  local r2="$sid_dir/${sid}_2.fastq"
  if [[ ! -f "$r1" || ! -f "$r2" ]]; then
    echo "ERROR: read files for $sid not found in $sid_dir" >&2
    exit 1
  fi

  mkdir -p "$trim_dir/$sid"
  local outP="$trim_dir/$sid/${sid}_1P.fastq"
  local outU1="$trim_dir/$sid/${sid}_1U.fastq"
  local outP2="$trim_dir/$sid/${sid}_2P.fastq"
  local outU2="$trim_dir/$sid/${sid}_2U.fastq"

  echo "[trim] $sid: running Trimmomatic"
  java -jar "$TRIMOMATIC_JAR" PE -threads "$THREADS" -phred33 \
    "$r1" "$r2" \
    "$outP" "$outU1" "$outP2" "$outU2" \
    ILLUMINACLIP:"$ADAPTERS":2:30:10 $TRIM_PARAMS

  touch "$done_flag"
}

run_novoplasty() {
  local sid="$1"
  local done_flag="$novo_dir/$sid.novo.done"
  if [[ -f "$done_flag" ]]; then
    echo "[novo] $sid already done, skipping"
    return
  fi

  local merged_fa="$novo_dir/$sid/${sid}_novoplasty_config.txt"
  mkdir -p "$novo_dir/$sid"

  # locate trimmed pair files
  local p1="$trim_dir/$sid/${sid}_1P.fastq"
  local p2="$trim_dir/$sid/${sid}_2P.fastq"
  if [[ ! -f "$p1" || ! -f "$p2" ]]; then
    echo "ERROR: trimmed paired files for $sid not found" >&2; exit 1
  fi

  cat > "$merged_fa" <<EOF
Project name          = ${sid}_novoplasty
Insert size           = 350
Insert size aut        = yes
Read Length           = 150
Type                  = mito
Genome Range          = 14000-20000
K-mer                 = 39
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            =
Reference sequence    =
Variance detection    = no
Chloroplast sequence  =
Dataset 1             = $p1
Dataset 2             = $p2
Dataset 3             =
Dataset 4             =

Optional:  
Output path           = $novo_dir/$sid
EOF

  echo "[novo] $sid: running NOVOPlasty"
  perl "$NOVOPLASTY_SCRIPT" -c "$merged_fa" > "$novo_dir/$sid/${sid}_novoplasty.log" 2>&1

  touch "$done_flag"
}

run_heteroplasmy() {
  local sid="$1"
  local done_flag="$hetero_dir/$sid.hetero.done"
  if [[ -f "$done_flag" ]]; then
    echo "[hetero] $sid already done, skipping"
    return
  fi

  local assembly="$novo_dir/$sid/Final_assembly_${sid}_novoplasty.fasta"
  # fallback for common filename
  if [[ ! -f "$assembly" ]]; then
    assembly="$novo_dir/$sid/Final_assembly.fasta"
  fi
  if [[ -n "$REF_FA" ]]; then
    assembly="$REF_FA"
  fi

  if [[ ! -f "$assembly" ]]; then
    echo "ERROR: assembly reference for heteroplasmy not found ($assembly)" >&2
    exit 1
  fi

  mkdir -p "$hetero_dir/$sid"
  local p1="$trim_dir/$sid/${sid}_1P.fastq"
  local p2="$trim_dir/$sid/${sid}_2P.fastq"

  # index reference
  bwa index "$assembly"

  # mapping
  samtools faidx "$assembly"
  bwa mem -t "$THREADS" "$assembly" "$p1" "$p2" | samtools view -bS - | samtools sort -@ "$THREADS" -o "$hetero_dir/$sid/${sid}_mt.bam"
  samtools index "$hetero_dir/$sid/${sid}_mt.bam"

  # variant calling
  bcftools mpileup -f "$assembly" "$hetero_dir/$sid/${sid}_mt.bam" | bcftools call -mv -Oz -o "$hetero_dir/$sid/${sid}_mt.raw.vcf.gz"
  bcftools index "$hetero_dir/$sid/${sid}_mt.raw.vcf.gz"
  bcftools filter -i 'MAF>0.01 && DP>100' "$hetero_dir/$sid/${sid}_mt.raw.vcf.gz" -Oz -o "$hetero_dir/$sid/${sid}_mt.filtered.vcf.gz"
  bcftools index "$hetero_dir/$sid/${sid}_mt.filtered.vcf.gz"

  # heteroplasmy summarization
  bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%DP\t%AF\n' "$hetero_dir/$sid/${sid}_mt.filtered.vcf.gz" > "$hetero_dir/$sid/${sid}_heteroplasmy.tsv"

  touch "$done_flag"
}

# work by step
start_from=0
case "$STEP" in
  all) start_from=0 ;;
  download) start_from=0 ;;
  trim) start_from=1 ;;
  novo) start_from=2 ;;
  hetero) start_from=3 ;;
  *) echo "Unknown step $STEP" >&2; usage ;; 
esac

for sid in "${sra_ids[@]}"; do
  sid=$(echo "$sid" | tr -d '\r\n' | xargs)
  [[ -z "$sid" ]] && continue

  if (( start_from <= 0 )); then run_download "$sid"; fi
  if (( start_from <= 1 )); then run_trim "$sid"; fi
  if (( start_from <= 2 )); then run_novoplasty "$sid"; fi
  if (( start_from <= 3 )); then run_heteroplasmy "$sid"; fi

  # after first sample if starting later, reset baseline to full pipeline in case user wants full chain
  start_from=0
  echo "[sample] $sid pipeline completed"
done

echo "Pipeline completed all samples. Output in $OUTDIR"

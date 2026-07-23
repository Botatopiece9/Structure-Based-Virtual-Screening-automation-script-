#!/usr/bin/env bash
set -euo pipefail

# vs_screen.sh
#vs_screen.sh - General SBDD Virtual Screening Pipline

#Defaults: declaring the defaults of the program
RECEPTOR=""
LIGAND_DIR=""
CONFIG=""
CENTER_X=""
CENTER_Y=""
CENTER_Z=""
SIZE_X="20"
SIZE_Y="20"
SIZE_Z="20"
OUTPUT_DIR="results_$(date +%Y%m%d_%H%M%S)"
TOP_N=10
EXHAUSTIVENESS=8
NUM_MODES=9
ENERGY_RANGE=3

#Usage functions: the message that will show asking for the input and guide the user for correct arguments
usage() {
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo " -r, --receptor      Receptor .pdb file (required)"
    echo " -l, --library       Ligand library directory (required)"
    echo " -c, --config        Vina grid box config file (required)"
    echo " -o, --output        Output directory (default: results_<timestamp>)"
    echo " -n, --top           Top hits to report (default: 10)"
    echo " -e, --exhaustivness Vina exhaustiveness (default: 8)"
    echo " -h, --help          Show this help message"
    exit 1
}

#Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -r|--receptor)          RECEPTOR="$2";      shift 2 ;;
        -l|--library)           LIGAND_DIR="$2";    shift 2 ;;
        -c|--config)            CONFIG="$2";        shift 2 ;;
        --center-x)             CENTER_X="$2";      shift 2 ;;
        --center-y)             CENTER_Y="$2";      shift 2 ;;
        --center-z)             CENTER_Z="$2";      shift 2 ;;
        --size)                 SIZE_X="$2"; SIZE_Y="$2"; SIZE_Z="$2"; shift 2 ;;
        -o|--output)            OUTPUT_DIR="$2";     shift 2 ;;
        -n|--top)               TOP_N="$2";          shift 2 ;;
        -e|--exhaustiveness)    EXHAUSTIVENESS="$2"; shift 2 ;;
        -h|--help)              usage;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

#Arguments validation: validation for the required paramters needed to start the docking process
if [[ -z "$RECEPTOR" ]]; then
    echo "Error: receptor file not specified. Use -r"
    usage
fi

if [[ -z "$LIGAND_DIR" ]]; then
    echo "Error: Ligand library not specified. Use -l"
    usage
fi

if [[ -z "$CONFIG" && ( -z "$CENTER_X" || -z "$CENTER_Y" || -z "$CENTER_Z" ) ]]; then
    echo "Error: provide either -c config file or --center-x/y/z coordinates"
    usage
fi

#Files and directories validation
if [[ ! -f "$RECEPTOR" ]]; then
    echo "Error: receptor file not found: $RECEPTOR"
    exit 1
fi
if [[ ! -d "$LIGAND_DIR" ]]; then
    echo "Error: ligand library not found: $LIGAND_DIR"
    exit 1
fi
if [[ -n "$CONFIG" && ! -f "$CONFIG" ]]; then
    echo "Error: grid box config file not found: $CONFIG"
    exit 1
fi

#Receptor extension validation: to ensure the formats are correct for the docking
EXT="${RECEPTOR##*.}"
if [[ "$EXT" != "pdb" && "$EXT" != "pdbqt" ]]; then
    echo "Error: receptor must be .pdb or .pdbqt (got .$EXT)"
    exit 1
fi

#Creating output directory
mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

#Check required tools: to ensure the tools are available on the system before the start of the docking
check_dependency() {
    local tool=$1
    local install_hint=$2

    if ! command -v "$tool" &>/dev/null; then
        echo "Error: '$tool' not found."
        echo " Install hint: $install_hint"
        exit 1
    else 
        echo " [ok] $tool ($(command -v "$tool"))"
    fi
}

#MGLTools path hint
MGLTOOLS_PATH="/opt/mgltools/bin"
if [[ -d "$MGLTOOLS_PATH" ]]; then
    export PATH="$MGLTOOLS_PATH:$PATH"
fi

echo "Checking dependecies..."
check_dependency "vina"                   "https://vina.scripps.edu/downloads/"
check_dependency "python3"                "sudo apt install python3"
check_dependency "prepare_receptor4.py"   "Install MGLTools: https://ccsb.scripps.edu/mgltools/"
check_dependency "prepare_ligand4.py"     "Install MGLTools: https://ccsb.scripps.edu/mgltools/"
echo "All dependencies found!"
echo ""

#Directories: assigning the directories to which results will be stored and retrived from
PREP_DIR="${OUTPUT_DIR}/prepared"
DOCK_DIR="${OUTPUT_DIR}/docked"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "$PREP_DIR" "$DOCK_DIR" "$LOG_DIR"

#Prepare receptor
prepare_receptor() {
    local input=$1
    local ext="${input##*.}"
    local name
    name=$(basename "$input" ".$ext")
    local out="${PREP_DIR}/${name}.pdbqt"

    if [[ "$ext" == "pdbqt" ]]; then
        echo " Receptor already in .pdbqt format - copying..."
        cp "$input" "$out"
    else
        echo " Preparing receptor: $name"
        prepare_receptor4.py \
            -r "$input" \
            -o "$out" \
            -A hydrogens \
            -U nphs_lps_waters_deleteAltB \
            2>>"${LOG_DIR}/receptor_prep.log"
    fi
    if [[ ! -f "$out" ]]; then
        echo "Error: receptor preparation failed. Check ${LOG_DIR}/receptor_prep.log"
        exit 1
    fi

    echo " [ok] Receptor prepared -> $out"
    RECEPTOR_PDBQT="$out"
}

#Prepare ligand
prepare_ligand() {
    local ligand_dir=$1
    local prepared=0
    local failed=0

    echo " Preparing ligands from: $ligand_dir"

    for ligand in "$ligand_dir"/*.{sdf,mol2,pdbqt}; do
        [[ -e "$ligand" ]] || continue

        local ext="${ligand##*.}"
        local name
        name=$(basename "$ligand" ".$ext")
        local out="${PREP_DIR}/${name}.pdbqt"

        if [[ -f "$out" ]]; then
            echo " [SKIP] Already prepared: $name"
            (( prepared++ )) || true
            continue
        fi

        if [[ "$ext" == "pdbqt" ]]; then
            cp "$ligand" "$out"
        else
            prepare_ligand4.py \
                -l "$ligand" \
                -o "$out" \
                -A hydrogens \
                2>>"${LOG_DIR}/ligand_prep.log"
        fi

        if [[ -f "$out" ]]; then
            echo "[OK] Prepared: $name"
            (( prepared++ )) || true
        else
            echo " [WARN] Failed: $name - check ${LOG_DIR}/ligand_prep.log"
            (( failed++ )) || true
        fi
    done

    echo " Ligands prepared: $prepared | Failed: $failed"

    if [[ $prepared -eq 0 ]]; then
        echo " [ERROR] No ligands were successfully prepared. Exiting."
        exit 1
    fi        
            
}

#Config file handling
setup_config() {
    local final_config="${OUTPUT_DIR}/vina_config.txt"

    if [[ -n "$CONFIG" ]]; then
        if [[ ! -f "$CONFIG" ]]; then
        echo "[ERROR] Config file not found: $CONFIG"
        exit 1
        fi

        echo "[INFO] Using provided config: $CONFIG"
        cp "$CONFIG" "$final_config"

        if ! grep -q "exhaustiveness" "$final_config"; then
            echo "exhaustiveness = $EXHAUSTIVENESS" >> "$final_config"
        fi

    elif [[ -n "$CENTER_X" && -n "$CENTER_Y" && -n "$CENTER_Z" ]]; then
        echo "[INFO] Building config from provided coordinates..."
        cat > "$final_config" <<EOF
center_x = $CENTER_X
center_y = $CENTER_Y
center_z = $CENTER_Z
size_x   = $SIZE_X
size_y   = $SIZE_Y
size_z   = $SIZE_Z
exhaustiveness = $EXHAUSTIVENESS
num_modes      = $NUM_MODES
energy_range   = $ENERGY_RANGE
EOF

    else
        echo "[ERROR] No grid box defined."
        echo "        Provide a config file with -c, or coordinates with --center-x/y/z"
        exit 1
    fi

    echo "[INFO] Config ready: $final_config"
    VINA_CONFIG="$final_config"     
}

#Docking loop: the fuction that run the docking process and will alow the process to repeat for all the ligands in the library
run_docking() {
    local dock_dir="${OUTPUT_DIR}/docking_results"
    mkdir -p "$dock_dir"

    local total=0
    local success=0
    local failed=0

    echo ""
    echo "[INFO] Starting docking run..."
    echo "[INFO] Receptor:         $RECEPTOR_PDBQT"
    echo "[INFO] Config:           $VINA_CONFIG"
    echo "[INFO] Output dir:       $dock_dir"
    echo ""

   local receptor_name
   receptor_name=$(basename "$RECEPTOR_PDBQT")
    
    for ligand in "${PREP_DIR}"/*.pdbqt; do
        [[ -e "$ligand" ]] || continue

        local name
        name=$(basename "$ligand" .pdbqt)

	if [[ "$(basename "$ligand")" == "$receptor_name" ]]; then
        echo " [SKIP] Skipping receptor: $name"
        continue
        fi

        local out="${dock_dir}/${name}_docked.pdbqt"
        local log="${dock_dir}/${name}.log"

        (( total++ )) || true

        echo -n " Docking $name ..."

        if vina \
            --receptor "$RECEPTOR_PDBQT" \
            --config   "$VINA_CONFIG" \
            --ligand   "$ligand" \
            --out      "$out" \
            --log      "$log" \
            2>/dev/null; then

            echo "Done!"
            (( success++ )) || true
        else
            echo "FAILED"
            (( failed++ )) || true
        fi
    done

    echo ""
    echo "[INFO] Docking complete."
    echo "[INFO] Total: $total | Success: $success | Failed: $failed"
    echo ""

    if [[ $success -eq 0 ]]; then
        echo "[ERROR] All docking runs failed. Check your config and receptor."
        exit 1
    fi

    DOCK_DIR="$dock_dir"
}

#Parse docking results

parse_results() {
    local raw="${OUTPUT_DIR}/raw_scores.tsv"
    local parsed=0
    local skipped=0

    echo "[INFO] Parsing docking logs..."

    echo -e "Compound\tAffinity_kcal_mol\tRSMD_lb\tRSMD_ub" > "$raw"

    for log in "${DOCK_DIR}"/*.log; do
        [[ -e "$log" ]] || continue

        local name
        name=$(basename "$log" .log)

        local result 
        result=$(grep -A1 -- "-----+------------" "$log" \
                  | grep -v "^\-\-" \
                  | grep -v "mode" \
                  | awk 'NR==1 {print $3"\t"$5"\t"$7}')
        
        if [[ -z "$result" ]]; then
            echo " [WARN] Could not parse: $name"
            (( skipped++ )) || true
            continue
        fi

        echo -e "${name}\t${result}" >> "$raw"
        (( parsed++ )) || true
    done

    echo "[INFO] Parsed: $parsed | Skipped: $skipped"
    echo "[INFO] Raw scores written to: $raw"

    RAW_SCORES="$raw"
}

#Final report: function that will generate the final report with details such as best compoud, best score etc.
generate_report() {
    local sorted="${OUTPUT_DIR}/results_sorted.tsv"
    local report="${OUTPUT_DIR}/report.txt"
    local top_hits="${OUTPUT_DIR}/top_hits.tsv"

    echo "[INFO] Generating final report..."

    head -1 "$RAW_SCORES" > "$sorted"
    tail -n +2 "$RAW_SCORES" | sort -k2 -n >> "$sorted"

    head -1 "$sorted" > "$top_hits"
    tail -n +2 "$sorted" | head -n "$TOP_N" >> "$top_hits"

    local total
    total=$(tail -n +2 "$sorted" | wc -l)

    local best_compound best_score
    best_compound=$(tail -n +2 "$sorted" | awk 'NR==1 {print $1}')
    best_score=$(tail -n +2 "$sorted" | awk 'NR==1 {print $2}')

    local mean_score
    mean_score=$(tail -n +2 "$sorted" | awk '{sum+=$2} END {printf "%.2f", sum/NR}')

    cat > "$report" << EOF
═══════════════════════════════════════════════════
  Virtual Screening Report
  $(date "+%Y-%m-%d %H:%M:%S")
═══════════════════════════════════════════════════

  Receptor:        $RECEPTOR
  Ligand library:  $LIGAND_DIR
  Config:          $VINA_CONFIG
  Exhaustiveness:  $EXHAUSTIVENESS

───────────────────────────────────────────────────
  Summary
───────────────────────────────────────────────────

  Compounds screened:   $total
  Best compound:        $best_compound
  Best affinity:        $best_score kcal/mol
  Mean affinity:        $mean_score kcal/mol
  Top N reported:       $TOP_N

───────────────────────────────────────────────────
  Top $TOP_N Hits
───────────────────────────────────────────────────

$(column -t "$top_hits")

═══════════════════════════════════════════════════
EOF

    echo "[INFO] Report written to:       $report"
    echo "[INFO] Sorted results:          $sorted"
    echo "[INFO] Top hits TSV:            $top_hits"
    echo ""
    echo "════════════════════════════════════════"
    cat "$report"
    echo "════════════════════════════════════════"
}

#Execution: calling the functions
prepare_receptor "$RECEPTOR"
prepare_ligand "$LIGAND_DIR"
setup_config
run_docking
parse_results
generate_report

echo "[INFO] Pipeline finished."
echo "[INFO] All output saved in: $OUTPUT_DIR"


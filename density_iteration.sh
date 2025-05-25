#!/bin/bash

# Funkció a CPU magok számának lekérdezésére
get_cpu_cores() {
    local cores
    if command -v nproc >/dev/null 2>&1; then
        cores=$(nproc)
    elif command -v sysctl >/dev/null 2>&1 && sysctl -n hw.logicalcpu >/dev/null 2>&1; then
        cores=$(sysctl -n hw.logicalcpu)
    elif [ -f /proc/cpuinfo ]; then
        cores=$(grep -c ^processor /proc/cpuinfo)
    else
        echo "Figyelmeztetés: Nem sikerült automatikusan meghatározni a CPU magok számát. 1 mag lesz használva." >&2
        cores=1
    fi
    if ! [[ "$cores" =~ ^[0-9]+$ ]] || [ "$cores" -lt 1 ]; then
        echo "Figyelmeztetés: Érvénytelen magszám ($cores) detektálva. 1 mag lesz használva." >&2
        cores=1
    fi
    echo "$cores"
}

# Clean up previous results - ÓVATOSAN HASZNÁLD, HA MÁS FONTOS .txt IS VAN!
# Lehet, hogy specifikusabb törlést szeretnél, pl. csak a generáltakat.
# rm -f *.txt # Ez töröl minden .txt-t az aktuális könyvtárban
rm -f tails_ratio_Q*
rm -f extra_calc_Q*
rm -rf tails_ratios
mkdir -p tails_ratios

# --- Compilation ---
echo "Compiling code..."
if g++ -std=c++20 integrator.cpp -fopenmp -o iteration_solver \
    -lboost_math_c99 -lboost_math_c99f -lboost_math_c99l \
    -lboost_math_tr1 -lboost_math_tr1f -lboost_math_tr1l; then
    echo "Compilation successful."
else
    echo "ERROR: Compilation failed!"
    exit 1
fi
echo ""


# --- Configuration ---
SOLVER_PROGRAM="./iteration_solver"
RESULTS_DIR="simulation_results" # Directory to store detailed logs and CSVs

# Core Simulation Parameters (exported for the C++ program)
NUM_CORES=$(get_cpu_cores)
export OMP_NUM_THREADS="$NUM_CORES"

# --- Fixed Parameters ---
export CONST_SIGMA="0.05"
export CONST_MU=$(echo "scale=10; ${CONST_SIGMA} * ${CONST_SIGMA} / 2" | bc) # mu = sigma^2 / 2
export GAMMA_PLUS="0.003"  # gamma
export GAMMA_MINUS="0.003" # gamma
export P_PROBABILITY=$(echo "scale=10; 1/12" | bc) # p = 1/12
export NUM_ITERATIONS="100" # num_iterations
export N_POINTS_GRID="201"  # grid_points

# Grid Parameters (exported for the C++ program) - ezeket is rögzítheted, ha szükséges
export X_MIN="-0.02" # Maradhat ez, vagy állítsd be
export X_MAX="0.02" # Maradhat ez, vagy állítsd be

# --- Q Values for Iteration ---
Q_VALUES_TO_ITERATE=("0.0" "0.2" "0.4" "0.8")


# --- Script Start ---
echo "Starting simulation series with the following global settings:"
echo "CONST_SIGMA: ${CONST_SIGMA}, CONST_MU: ${CONST_MU}"
echo "GAMMA: ${GAMMA_PLUS} (GP & GM)"
echo "P_PROBABILITY: ${P_PROBABILITY}"
echo "NUM_ITERATIONS: ${NUM_ITERATIONS}"
echo "N_POINTS_GRID: ${N_POINTS_GRID}"
echo "X_MIN: ${X_MIN}, X_MAX: ${X_MAX}"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
echo "Q_FOR_H will iterate through: ${Q_VALUES_TO_ITERATE[*]}"
echo "---------------------------------------------------------------------"

mkdir -p "${RESULTS_DIR}"

if [ ! -f "${SOLVER_PROGRAM}" ]; then
    echo "Error: Solver program (${SOLVER_PROGRAM}) not found!"
    exit 1
fi
if [ ! -x "${SOLVER_PROGRAM}" ]; then
    echo "Error: Solver program (${SOLVER_PROGRAM}) not executable!"
    exit 1
fi

total_sims=${#Q_VALUES_TO_ITERATE[@]}
current_sim=0

# Iteration over Q_FOR_H values
for Q_H in "${Q_VALUES_TO_ITERATE[@]}"; do
    export Q_FOR_H="${Q_H}" # Export current Q_FOR_H

    current_sim=$((current_sim + 1))
    echo ""
    echo "--- Running Simulation ${current_sim}/${total_sims} ---"

    # Parameter string for filenames
    sigma_fn=$(echo "${CONST_SIGMA}" | tr '.' 'p' | sed 's/-/neg/g')
    mu_fn=$(echo "${CONST_MU}" | tr '.' 'p' | sed 's/-/neg/g')
    gamma_fn=$(echo "${GAMMA_PLUS}" | tr '.' 'p' | sed 's/-/neg/g') # GAMMA_PLUS-t használjuk, mert egyenlőek
    q_h_fn=$(echo "${Q_FOR_H}" | tr '.' 'p')
    p_prob_fn=$(echo "${P_PROBABILITY}" | tr '.' 'p' | sed 's/-/neg/g')
    xmin_fn=$(echo "${X_MIN}" | tr '.' 'p' | sed 's/-/neg/g')
    xmax_fn=$(echo "${X_MAX}" | tr '.' 'p' | sed 's/-/neg/g')
    npoints_fn="${N_POINTS_GRID}"
    niter_fn="${NUM_ITERATIONS}"

    # A fájlnév most főleg a Q_H-tól fog függeni, ha a többi rögzített
    # Ha a C++ program más paramétereket is használ a fájlnévben (pl. NPT, NIT a Q_FOR_H-n kívül),
    # akkor azokat is bele kell venni. Most feltételezem, hogy a param_suffix_log_csv elegendő.
    param_suffix_log_csv="S${sigma_fn}_M${mu_fn}_G${gamma_fn}_QH${q_h_fn}_P${p_prob_fn}_XMIN${xmin_fn}_XMAX${xmax_fn}_NPT${npoints_fn}_NIT${niter_fn}"

    echo "Parameters for this run:"
    echo "  CONST_SIGMA: ${CONST_SIGMA}, CONST_MU: ${CONST_MU}"
    echo "  GAMMA: ${GAMMA_PLUS} (GP: ${GAMMA_PLUS}, GM: ${GAMMA_MINUS})"
    echo "  Q_FOR_H: ${Q_FOR_H}, P_PROBABILITY: ${P_PROBABILITY}"
    echo "  X_MIN: ${X_MIN}, X_MAX: ${X_MAX}"
    echo "  N_POINTS_GRID: ${N_POINTS_GRID}, NUM_ITERATIONS: ${NUM_ITERATIONS}"
    echo "  OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
    echo "Log/CSV Suffix: ${param_suffix_log_csv}"

    # A C++ program ezt a nevet használja, töröljük az előző futásból, ha létezik
    rm -f f_functions_output.csv 

    log_filename="log_${param_suffix_log_csv}.txt"
    log_path="${RESULTS_DIR}/${log_filename}"
    rm -f "${log_path}" # Töröljük az előző futásból származó logot ehhez a paraméterkombinációhoz

    echo "Running solver... Saving log to: ${log_path}"
    if "${SOLVER_PROGRAM}" >"${log_path}" 2>&1; then
        echo "Solver finished successfully."
    else
        echo "ERROR: Solver program returned an error! Check log: ${log_path}"
        # Hiba esetén itt megállhat a szkript, ha szükséges:
        # exit 1
    fi

    if [ -f "f_functions_output.csv" ]; then
        output_csv_filename="f_out_${param_suffix_log_csv}.csv"
        destination_csv_path="${RESULTS_DIR}/${output_csv_filename}"
        mv "f_functions_output.csv" "${destination_csv_path}"
        echo "Detailed CSV result saved to: ${destination_csv_path}"
    else
        echo "WARNING: f_functions_output.csv was not created by the solver. Check log: ${log_path}"
    fi
done # End of Q_VALUES_TO_ITERATE loop

echo ""
echo "Moving any generated .txt files (tails_ratio and extra_calc) to tails_ratios/ directory..."
# Mivel a Q_FOR_H változik, a fájlnevekben a QH rész is változni fog.
# A 'tails_ratio_Q*' és 'extra_calc_Q*' minták általánosabbak, minden Q-val kezdődőt áthelyeznek.
find . -maxdepth 1 -type f -name 'tails_ratio_Q*.txt' -exec mv -t tails_ratios/ {} + 2>/dev/null || echo "No tails_ratio_Q* files found to move or error during move."
find . -maxdepth 1 -type f -name 'extra_calc_Q*.txt' -exec mv -t tails_ratios/ {} + 2>/dev/null || echo "No extra_calc_Q* files found to move or error during move."

echo "--- Simulation series finished ---"

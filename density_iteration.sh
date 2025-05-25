#!/bin/bash

# Funkció a CPU magok számának lekérdezésére
get_cpu_cores_physical() {
    local cores
    # macOS
    if command -v sysctl >/dev/null 2>&1 && sysctl -n hw.physicalcpu_max >/dev/null 2>&1; then # hw.physicalcpu_max a biztosabb újabb macOS-eken
        cores=$(sysctl -n hw.physicalcpu_max)
    elif command -v sysctl >/dev/null 2>&1 && sysctl -n hw.physicalcpu >/dev/null 2>&1; then # Régebbi macOS, vagy ha a _max nem elérhető
        cores=$(sysctl -n hw.physicalcpu)
    # Linux - lscpu a legmegbízhatóbb, ha elérhető
    elif command -v lscpu >/dev/null 2>&1; then
        # lscpu kimenetéből: Cores per socket * Sockets
        local cores_per_socket=$(lscpu | grep -E "^Core\(s\) per socket:" | awk '{print $4}')
        local sockets=$(lscpu | grep -E "^Socket\(s\):" | awk '{print $2}')
        if [[ -n "$cores_per_socket" && -n "$sockets" && "$cores_per_socket" =~ ^[0-9]+$ && "$sockets" =~ ^[0-9]+$ ]]; then
            cores=$((cores_per_socket * sockets))
        else # Fallback, ha a fenti nem sikerül, próbálkozás a CPU(s) sorból a Core(s) per socket-tel
              # Ez egy kicsit heurisztikusabb, egy socketes gépekre jó lehet
            cores=$(lscpu -p=CORE | grep -v '^#' | sort -u | wc -l)
            if ! [[ "$cores" =~ ^[0-9]+$ ]] || [ "$cores" -lt 1 ]; then
                 # Végső lscpu fallback, ha minden más kudarcot vall lscpu-val
                 # Ez a 'CPU(s)' sorból veszi a logikaiakat, és elosztja a 'Thread(s) per core'-ral
                 local logical_cpus=$(lscpu | grep -E "^CPU\(s\):" | awk '{print $2}')
                 local threads_per_core=$(lscpu | grep -E "^Thread\(s\) per core:" | awk '{print $4}')
                 if [[ "$logical_cpus" =~ ^[0-9]+$ && "$threads_per_core" =~ ^[0-9]+$ && "$threads_per_core" -gt 0 ]]; then
                     cores=$((logical_cpus / threads_per_core))
                 else
                     cores="" # Jelzi, hogy nem sikerült
                 fi
            fi
        fi
    # Linux - /sys/devices/system/cpu/cpu*/topology/core_id (másik megbízható módszer)
    elif [ -d /sys/devices/system/cpu ]; then
        # Számolja az egyedi core_id-kat a cpu topológiából
        cores=$(find /sys/devices/system/cpu/cpu[0-9]*/topology/core_id -print0 2>/dev/null | xargs -0 cat 2>/dev/null | sort -u | wc -l)
        if ! [[ "$cores" =~ ^[0-9]+$ ]] || [ "$cores" -lt 1 ]; then
            cores="" # Ha a find/xargs/cat sor hibát ad vagy nem számot, töröljük
        fi
    # Linux - /proc/cpuinfo (kevésbé megbízható a fizikai magok számának direkt lekérdezésére több socket esetén)
    elif [ -f /proc/cpuinfo ]; then
        # Ez a módszer feltételezi, hogy a "cpu cores" sor minden fizikai CPU leírásában egyszer szerepel
        # és az ott lévő érték az egy processzorban lévő magok száma. Több socket esetén össze kell adni.
        # Egyszerűsítésként az első "cpu cores" értéket vesszük, ami egyprocesszoros rendszereken jó.
        cores=$(grep -m1 "cpu cores" /proc/cpuinfo | awk '{print $4}')
        # Ha több processzor van (socket), ez bonyolultabb:
        # num_sockets=$(grep "physical id" /proc/cpuinfo | sort -u | wc -l)
        # if [ "$num_sockets" -gt 1 ]; then
        #    cores_per_socket_val=$(grep -m1 "cpu cores" /proc/cpuinfo | awk '{print $4}')
        #    cores=$(( num_sockets * cores_per_socket_val ))
        # fi
        if ! [[ "$cores" =~ ^[0-9]+$ ]] || [ "$cores" -lt 1 ]; then
            cores="" # Ha nem számot adott vissza
        fi
    fi

    # Végső ellenőrzés és alapértelmezett érték
    if ! [[ "$cores" =~ ^[0-9]+$ ]] || [ "$cores" -lt 1 ]; then
        echo "Figyelmeztetés: Nem sikerült meghatározni a FIZIKAI CPU magok számát. Fallback logikai magokra..." >&2
        # Itt visszatérhetünk az eredeti logikai magokat lekérdező függvényhez, vagy beállíthatunk egy fix értéket
        # Most az eredeti logikát használjuk fallbackként:
        if command -v nproc >/dev/null 2>&1; then
            cores=$(nproc)
        elif command -v sysctl >/dev/null 2>&1 && sysctl -n hw.logicalcpu >/dev/null 2>&1; then
            cores=$(sysctl -n hw.logicalcpu)
        elif [ -f /proc/cpuinfo ]; then
            cores=$(grep -c ^processor /proc/cpuinfo)
        else
            echo "Figyelmeztetés: A logikai magok száma sem határozható meg. 1 mag lesz használva." >&2
            cores=1
        fi
        
        # Újabb ellenőrzés a fallback után
        if ! [[ "$cores" =~ ^[0-9]+$ ]] || [ "$cores" -lt 1 ]; then
            echo "Figyelmeztetés: Érvénytelen magszám (fallback után is). 1 mag lesz használva." >&2
            cores=1
        fi
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
NUM_CORES=$(get_cpu_cores_physical)
export OMP_NUM_THREADS="$NUM_CORES"

# --- Fixed Parameters ---
export CONST_SIGMA="0.05"
export CONST_MU=$(echo "scale=10; ${CONST_SIGMA} * ${CONST_SIGMA} / 2" | bc) # mu = sigma^2 / 2
export GAMMA_PLUS="0.003"  # gamma
export GAMMA_MINUS="0.003" # gamma
export P_PROBABILITY=$(echo "scale=10; 1/12" | bc) # p = 1/12
export NUM_ITERATIONS="1000" # num_iterations
export N_POINTS_GRID="801"  # grid_points

# Grid Parameters (exported for the C++ program) - ezeket is rögzítheted, ha szükséges
export X_MIN="-0.01" # Maradhat ez, vagy állítsd be
export X_MAX="0.01" # Maradhat ez, vagy állítsd be

# --- Q Values for Iteration ---
Q_VALUES_TO_ITERATE=("0.0" "0.05" "0.1" "0.2")


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

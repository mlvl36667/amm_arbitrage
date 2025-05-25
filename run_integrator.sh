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


# Clean up previous results
rm -rf *.txt
rm -rf tails_ratio_Q*
rm -rf extra_calc_Q* # Also remove extra calculation files
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


# OMP_NUM_THREADS beállítása és exportálása
export OMP_NUM_THREADS="$NUM_CORES"
export CONST_SIGMA="0.05"
# MU is calculated based on CONST_SIGMA later

# Grid and Iteration Parameters (exported for the C++ program)
export X_MIN="-0.02" # Fixed X_MIN
export X_MAX="0.02" # Fixed X_MAX
export N_POINTS_GRID="801" # Fixed N_POINTS_GRID

# Fixed SIGMA value (used for MU calculation and parameter sweeps)
FIXED_SIGMA="${CONST_SIGMA}"

# --- Parameter Sweeps Configuration ---
# NUM_ITERATIONS values to iterate over
NUM_ITERATIONS_VALUES=("100" "500" "1000" "2000" "4000" "8000" "16000") # New iteration loop

# GAMMA values (in basis points, will be converted)
GAMMA_BP_VALUES=("1" "5" "10" "30" "100")

# Q_FOR_H value is fixed
Q_VALUES=("0.0") # Q is now fixed at 0.0

# P_PROBABILITY values (as expressions, calculated with bc)
P_EXPRESSIONS=(
    "1/(10*60)"       # approx. 0.001666
    "1/(2*60)"        # approx. 0.008333
    "1/12"            # approx. 0.083333
    "1/2"             # 0.5 (A te eredeti kommented "2" volt, de 1/2 = 0.5)
)


# --- Script Start ---
echo "Starting simulation series with the following global settings:"
echo "X_MIN: ${X_MIN}, X_MAX: ${X_MAX}, N_POINTS_GRID: ${N_POINTS_GRID}"
echo "OMP_NUM_THREADS: ${OMP_NUM_THREADS}, CONST_SIGMA: ${CONST_SIGMA}"
echo "Q_FOR_H will be fixed at: ${Q_VALUES[0]}"
echo "NUM_ITERATIONS will iterate through: ${NUM_ITERATIONS_VALUES[*]}"
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

# Calculate MU based on FIXED_SIGMA (which is CONST_SIGMA)
export CONST_MU=$(echo "scale=10; ${FIXED_SIGMA} * ${FIXED_SIGMA} / 2" | bc)

# Update total_sims calculation
total_sims=$(( ${#NUM_ITERATIONS_VALUES[@]} * ${#GAMMA_BP_VALUES[@]} * ${#Q_VALUES[@]} * ${#P_EXPRESSIONS[@]} ))
current_sim=0

for CURRENT_NUM_ITERATIONS in "${NUM_ITERATIONS_VALUES[@]}"; do
  export NUM_ITERATIONS="${CURRENT_NUM_ITERATIONS}" # Export current NUM_ITERATIONS

  for Q_H in "${Q_VALUES[@]}"; do # This loop will only run once for Q_H = "0.0"
    export Q_FOR_H="${Q_H}"
    for GAMMA_BP in "${GAMMA_BP_VALUES[@]}"; do
      CURRENT_GAMMA=$(echo "scale=10; ${GAMMA_BP} / 10000" | bc)
      export GAMMA_PLUS="${CURRENT_GAMMA}"
      export GAMMA_MINUS="${CURRENT_GAMMA}"
          for P_EXPR in "${P_EXPRESSIONS[@]}"; do
              CURRENT_P_PROB=$(echo "scale=10; ${P_EXPR}" | bc)
              export P_PROBABILITY="${CURRENT_P_PROB}"

              # A P_PROBABILITY ellenőrzése: 0 < P < 1 (a te eredeti kódod P < 1000-t ellenőrzött, de a valószínűségeknek 0 és 1 között kell lenniük)
              # Ha a P_PROBABILITY tényleg lehet 1-nél nagyobb a te modelledben, akkor módosítsd ezt az ellenőrzést.
              is_p_valid=$(echo "${CURRENT_P_PROB} > 0 && ${CURRENT_P_PROB} < 1" | bc -l)
              if [ "$is_p_valid" -ne 1 ]; then
                  echo ""
                  echo "--- Skipping simulation due to invalid P_PROBABILITY ---"
                  echo "Expression: ${P_EXPR} -> Calculated P: ${CURRENT_P_PROB}"
                  echo "P_PROBABILITY must be > 0 and < 1 (for standard probabilities)."
                  continue
              fi

              current_sim=$((current_sim + 1))
              echo ""
              echo "--- Running Simulation ${current_sim}/${total_sims} ---"

              sigma_fn=$(echo "${CONST_SIGMA}" | tr '.' 'p' | sed 's/-/neg/g')
              mu_fn=$(echo "${CONST_MU}" | tr '.' 'p' | sed 's/-/neg/g')
              gamma_fn=$(echo "${CURRENT_GAMMA}" | tr '.' 'p' | sed 's/-/neg/g')
              q_h_fn=$(echo "${Q_FOR_H}" | tr '.' 'p') # Nincs szükség sed-re, mert Q_FOR_H=0.0
              p_prob_fn=$(echo "${P_PROBABILITY}" | tr '.' 'p' | sed 's/-/neg/g')
              xmin_fn=$(echo "${X_MIN}" | tr '.' 'p' | sed 's/-/neg/g')
              xmax_fn=$(echo "${X_MAX}" | tr '.' 'p' | sed 's/-/neg/g')
              npoints_fn="${N_POINTS_GRID}"
              niter_fn="${NUM_ITERATIONS}" # Uses the currently exported NUM_ITERATIONS

              param_suffix_log_csv="S${sigma_fn}_M${mu_fn}_G${gamma_fn}_QH${q_h_fn}_P${p_prob_fn}_XMIN${xmin_fn}_XMAX${xmax_fn}_NPT${npoints_fn}_NIT${niter_fn}"

              echo "Parameters for this run:"
              echo "  CONST_SIGMA: ${CONST_SIGMA}, CONST_MU: ${CONST_MU}"
              echo "  GAMMA: ${CURRENT_GAMMA} (GP: ${GAMMA_PLUS}, GM: ${GAMMA_MINUS})"
              echo "  Q_FOR_H: ${Q_FOR_H}, P_PROBABILITY: ${P_PROBABILITY}"
              echo "  X_MIN: ${X_MIN}, X_MAX: ${X_MAX}"
              echo "  N_POINTS_GRID: ${N_POINTS_GRID}, NUM_ITERATIONS: ${NUM_ITERATIONS}"
              echo "  OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
              echo "Log/CSV Suffix: ${param_suffix_log_csv}"

              rm -f f_functions_output.csv # A C++ program ezt a nevet használja

              log_filename="log_${param_suffix_log_csv}.txt"
              log_path="${RESULTS_DIR}/${log_filename}"
              rm -f "${log_path}" # Töröljük az előző futásból származó logot

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
          done
      done
    done
done # End of NUM_ITERATIONS_VALUES loop

echo ""
echo "Moving all generated .txt files (tails_ratio and extra_calc) to tails_ratios/ directory..."
# Mivel a Q_FOR_H mindig 0.0, a fájlnevekben QH0p0 lesz
find . -maxdepth 1 -type f -name 'tails_ratio_Q0p0*.txt' -exec mv -t tails_ratios/ {} + 2>/dev/null || echo "No tails_ratio_Q0p0 files found to move or error during move."
find . -maxdepth 1 -type f -name 'extra_calc_Q0p0*.txt' -exec mv -t tails_ratios/ {} + 2>/dev/null || echo "No extra_calc_Q0p0 files found to move or error during move."

echo "--- Simulation series finished ---"

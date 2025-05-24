#!/bin/bash

# --- Parancssori argumentumok feldolgozása ---
USE_Q_FROM_PROCESSOR=true # Alapértelmezés: használjuk a Q-t a processor.py-ból
if [[ "$1" == "--no-q" ]]; then
    USE_Q_FROM_PROCESSOR=false
    echo "INFO: --no-q flag specified. Q parameter from processor.py will be ignored, and Q_FOR_H will be set to 0.0."
    shift # Eltávolítja a feldolgozott argumentumot
fi

# Clean up previous results
rm -f *.txt
rm -rf tails_ratio_Q*
rm -rf extra_calc_Q*
rm -rf tails_ratios
mkdir -p tails_ratios
rm -f estimated_sigma_mu_q.txt

# --- 0. Run processor.py to get SIGMA, MU, and Q_JUMP ---
echo "--- Running processor.py to estimate SIGMA, MU, and Q_JUMP ---"
if python3 processor.py; then
    echo "processor.py finished successfully."
else
    echo "ERROR: processor.py failed to run or encountered an error!"
    exit 1
fi

ESTIMATED_PARAMS_FILE="estimated_sigma_mu_q.txt"
if [ ! -f "${ESTIMATED_PARAMS_FILE}" ]; then
    echo "ERROR: Estimated parameters file '${ESTIMATED_PARAMS_FILE}' not found after running processor.py!"
    exit 1
fi

echo "--- Reading estimated parameters from ${ESTIMATED_PARAMS_FILE} ---"
source "${ESTIMATED_PARAMS_FILE}" # Beolvassa CONST_SIGMA, CONST_MU, Q_JUMP_PROBABILITY

# --- Paraméterek ellenőrzése és véglegesítése ---
# CONST_SIGMA (mindig a processor.py-ból jön)
if [ -z "${CONST_SIGMA}" ] || [ "${CONST_SIGMA}" == "N/A" ]; then
    echo "ERROR: CONST_SIGMA (from ${ESTIMATED_PARAMS_FILE}) not read or is N/A. Exiting."
    exit 1
fi

# CONST_MU (a processor.py-ból jövő érték, mielőtt esetleg felülírnánk)
PROCESSOR_MU="${CONST_MU}" # Mentsük el a processor.py MU-ját az üzenethez/ellenőrzéshez
if [ -z "${PROCESSOR_MU}" ] || [ "${PROCESSOR_MU}" == "N/A" ]; then
    echo "ERROR: CONST_MU (from ${ESTIMATED_PARAMS_FILE}) not read or is N/A. Exiting."
    exit 1
fi

# Q_JUMP_PROBABILITY kezelése
Q_FROM_PROCESSOR="${Q_JUMP_PROBABILITY}" # Mentsük el a processor.py Q-ját
EFFECTIVE_Q=""

if [ "${USE_Q_FROM_PROCESSOR}" = true ]; then
    if [ -z "${Q_FROM_PROCESSOR}" ] || [ "${Q_FROM_PROCESSOR}" == "N/A" ]; then
        echo "ERROR: Q_JUMP_PROBABILITY (from ${ESTIMATED_PARAMS_FILE}) not read or is N/A."
        echo "       Q usage from processor.py is ENABLED (run with --no-q to ignore processor's Q and use 0.0)."
        exit 1
    fi
    EFFECTIVE_Q="${Q_FROM_PROCESSOR}"
else
    # --no-q van megadva. Q_FROM_PROCESSOR lehet "N/A" vagy hiányozhat, nem probléma.
    EFFECTIVE_Q="0.0"
fi
export Q_FOR_H="${EFFECTIVE_Q}" # Exportáljuk a szimulációhoz használandó Q értéket

# CONST_MU felülírása a szkriptben definiált fix értékkel
# Az eredeti szkript ezt a felülírást tartalmazta.
# Ha a processor.py MU-ját kellene használni, ezt a következő két sort ki kell kommentezni/törölni.
FIXED_MU_VALUE="0.03712680"
export CONST_MU="${FIXED_MU_VALUE}" # Ez felülírja a source-olt CONST_MU-t

# A CONST_SIGMA már be van állítva a 'source' paranccsal és exportálva lesz, ha nem volt az.
# A biztonság kedvéért explicit export:
export CONST_SIGMA="${CONST_SIGMA}"


echo "--- Final parameters for simulation: ---"
echo "  CONST_SIGMA (from processor.py): ${CONST_SIGMA}"

if [ "${CONST_MU}" == "${PROCESSOR_MU}" ]; then
    # Ez az ág akkor fut, ha a FIXED_MU_VALUE megegyezik a PROCESSOR_MU-val, vagy ha a felülírás ki van kommentezve
    echo "  CONST_MU (from processor.py): ${CONST_MU}"
elif [ -n "${FIXED_MU_VALUE}" ] && [ "${CONST_MU}" == "${FIXED_MU_VALUE}" ]; then
    echo "  CONST_MU (OVERRIDDEN by script): ${CONST_MU} (original from processor.py: '${PROCESSOR_MU}')"
else
    # Ha a felülírás ki van kommentezve, és a PROCESSOR_MU nem üres
    echo "  CONST_MU (from processor.py): ${CONST_MU}"
fi


if [ "${USE_Q_FROM_PROCESSOR}" = true ]; then
    echo "  Q_FOR_H (using Q_JUMP_PROBABILITY from processor.py): ${Q_FOR_H}"
else
    echo "  Q_FOR_H (SET TO 0.0 due to --no-q): ${Q_FOR_H} (Q_JUMP_PROBABILITY from processor.py: '${Q_FROM_PROCESSOR:-'Not provided or N/A'}' was IGNORED)"
fi
echo "--------------------------------------------------------"


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

# --- Configuration for a single targeted run ---
SOLVER_PROGRAM="./iteration_solver"
RESULTS_DIR="simulation_results_targeted_run" # Külön könyvtár ennek a futásnak

# Core Simulation Parameters (CONST_SIGMA, CONST_MU, Q_FOR_H már exportálva vannak)
export OMP_NUM_THREADS=4 # Használj annyi szálat, amennyi hatékony

# Grid and Iteration Parameters
export X_MIN="-0.02"
export X_MAX="0.02"
export N_POINTS_GRID="200"
export NUM_ITERATIONS="5000"

# Fixed parameters for this specific run
TARGET_P_PROBABILITY_EXPR="1/12"
TARGET_GAMMA_VALUE="0.003" # 30 bp

# --- Script Start ---
echo "Starting targeted simulation run with the following parameters:"
echo "  (Exported) CONST_SIGMA: ${CONST_SIGMA}"
echo "  (Exported) CONST_MU: ${CONST_MU}"
echo "  (Exported) Q_FOR_H: ${Q_FOR_H}"
echo "  P_PROBABILITY (expression): ${TARGET_P_PROBABILITY_EXPR}"
echo "  GAMMA_PLUS / GAMMA_MINUS: ${TARGET_GAMMA_VALUE}"
echo "  X_MIN: ${X_MIN}, X_MAX: ${X_MAX}"
echo "  N_POINTS_GRID: ${N_POINTS_GRID}"
echo "  NUM_ITERATIONS: ${NUM_ITERATIONS}"
echo "  OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
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

# Beállítjuk a konkrét P és Gamma értékeket exportálásra
export P_PROBABILITY=$(echo "scale=10; ${TARGET_P_PROBABILITY_EXPR}" | bc)
export GAMMA_PLUS="${TARGET_GAMMA_VALUE}"
export GAMMA_MINUS="${TARGET_GAMMA_VALUE}"

# --- Futtatás ---
echo ""
echo "--- Running Single Targeted Simulation ---"

# Fájlnévbarát stringek készítése a paraméterekből
sigma_fn=$(echo "${CONST_SIGMA}" | tr '.' 'p' | sed 's/-/neg/g')
mu_fn=$(echo "${CONST_MU}" | tr '.' 'p' | sed 's/-/neg/g')
gamma_fn=$(echo "${GAMMA_PLUS}" | tr '.' 'p' | sed 's/-/neg/g')
q_h_fn=$(echo "${Q_FOR_H}" | tr '.' 'p' | sed 's/-/neg/g') # Ha Q_FOR_H=0.0, q_h_fn=0p0
p_prob_fn=$(echo "${P_PROBABILITY}" | tr '.' 'p' | sed 's/-/neg/g')
xmin_fn=$(echo "${X_MIN}" | tr '.' 'p' | sed 's/-/neg/g')
xmax_fn=$(echo "${X_MAX}" | tr '.' 'p' | sed 's/-/neg/g')
npoints_fn="${N_POINTS_GRID}"
niter_fn="${NUM_ITERATIONS}"

param_suffix_log_csv="TARGETED_S${sigma_fn}_M${mu_fn}_G${gamma_fn}_QH${q_h_fn}_P${p_prob_fn}_XMIN${xmin_fn}_XMAX${xmax_fn}_NPT${npoints_fn}_NIT${niter_fn}"

rm -f f_functions_output.csv

log_filename="log_${param_suffix_log_csv}.txt"
log_path="${RESULTS_DIR}/${log_filename}"
rm -f "${log_path}"

echo "Running solver... Saving log to: ${log_path}"
if "${SOLVER_PROGRAM}" >"${log_path}" 2>&1; then
    echo "Solver finished successfully."
else
    echo "ERROR: Solver program returned an error! Check log: ${log_path}"
    # exit 1 # Hiba esetén itt megállhat a szkript, ha szükséges
fi

if [ -f "f_functions_output.csv" ]; then
    output_csv_filename="f_out_${param_suffix_log_csv}.csv"
    destination_csv_path="${RESULTS_DIR}/${output_csv_filename}"
    mv "f_functions_output.csv" "${destination_csv_path}"
    echo "Detailed CSV result saved to: ${destination_csv_path}"
else
    echo "WARNING: f_functions_output.csv was not created by the solver. Check log: ${log_path}"
fi

echo ""
echo "Moving generated .txt files (tails_ratio and extra_calc) to ${RESULTS_DIR} and tails_ratios/ ..."
# Feltételezzük, hogy a C++ program a Q_FOR_H értékét (még ha 0 is) beépíti a fájlnévbe, pl. tails_ratio_Q0p0...txt
# Így a 'tails_ratio_Q*.txt' minta továbbra is működik.
find . -maxdepth 1 -type f -name 'tails_ratio_Q*.txt' -exec mv -t "${RESULTS_DIR}/" {} + 2>/dev/null
find . -maxdepth 1 -type f -name 'extra_calc_Q*.txt' -exec mv -t "${RESULTS_DIR}/" {} + 2>/dev/null

cp "${RESULTS_DIR}/tails_ratio_Q"*.txt tails_ratios/ 2>/dev/null || echo "INFO: No tails_ratio_Q*.txt files found in ${RESULTS_DIR} to copy to tails_ratios/, or error during copy."
cp "${RESULTS_DIR}/extra_calc_Q"*.txt tails_ratios/ 2>/dev/null || echo "INFO: No extra_calc_Q*.txt files found in ${RESULTS_DIR} to copy to tails_ratios/, or error during copy."

echo "--- Targeted simulation run finished ---"

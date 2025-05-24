#!/bin/bash

# Clean up previous results
rm -f *.txt
rm -rf tails_ratio_Q*
rm -rf extra_calc_Q*
rm -rf tails_ratios
mkdir -p tails_ratios
rm -f estimated_sigma_mu_q.txt # Új fájlnév a q-val

# --- 0. Run processor.py to get SIGMA, MU, and Q_JUMP ---
echo "--- Running processor.py to estimate SIGMA, MU, and Q_JUMP ---"
# Tegyük fel, hogy a processor.py az aktuális könyvtárban van
# és a szükséges bemeneti fájl (november_rates.csv) is ott van,
# vagy a processor.py tudja, hol keresse.
# Módosítsd a processor.py-t, hogy kiírja a q_jump_probability-t is!
if python3 processor.py; then # Feltételezzük, hogy processor.py kiírja a q-t is
    echo "processor.py finished successfully."
else
    echo "ERROR: processor.py failed to run or encountered an error!"
    exit 1
fi

ESTIMATED_PARAMS_FILE="estimated_sigma_mu_q.txt" # Új fájlnév
if [ ! -f "${ESTIMATED_PARAMS_FILE}" ]; then
    echo "ERROR: Estimated parameters file '${ESTIMATED_PARAMS_FILE}' not found!"
    exit 1
fi

echo "--- Reading estimated parameters from ${ESTIMATED_PARAMS_FILE} ---"
source "${ESTIMATED_PARAMS_FILE}" # Beolvassa CONST_SIGMA, CONST_MU, Q_JUMP_PROBABILITY

# Ellenőrzések (maradnak, de most a Q_JUMP_PROBABILITY-t is nézzük)
if [ -z "${CONST_SIGMA}" ] || [ "${CONST_SIGMA}" == "N/A" ]; then
    echo "ERROR: CONST_SIGMA not read or N/A. Exiting."
    exit 1
fi
if [ -z "${CONST_MU}" ] || [ "${CONST_MU}" == "N/A" ]; then
    echo "ERROR: CONST_MU not read or N/A. Exiting."
    exit 1
fi
if [ -z "${Q_JUMP_PROBABILITY}" ] || [ "${Q_JUMP_PROBABILITY}" == "N/A" ]; then
    echo "ERROR: Q_JUMP_PROBABILITY not read or N/A. Exiting."
    exit 1
fi

echo "Using parameters from processor.py:"
echo "  CONST_SIGMA=${CONST_SIGMA}"
echo "  CONST_MU=${CONST_MU}"
echo "  Q_FOR_H (from jump probability) = ${Q_JUMP_PROBABILITY}"
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

# Core Simulation Parameters (CONST_SIGMA és CONST_MU már exportálva van a `source` paranccsal)
export OMP_NUM_THREADS=4 # Használj annyi szálat, amennyi hatékony
export Q_FOR_H="${Q_JUMP_PROBABILITY}" # A processzorból jövő q érték

# Grid and Iteration Parameters
export X_MIN="-0.02"  # Ezt finomíthatod az adatok alapján
export X_MAX="0.02"   # Ezt finomíthatod az adatok alapján
export N_POINTS_GRID="200" # "Jó nagy" grid
export NUM_ITERATIONS="100" # Az 1000. függvényhez (f0...f999)

# Fixed parameters for this specific run
TARGET_P_PROBABILITY_EXPR="1/12"
TARGET_GAMMA_VALUE="0.003" # 30 bp

# --- Script Start ---
echo "Starting targeted simulation run with the following parameters:"
echo "  CONST_SIGMA: ${CONST_SIGMA}"
echo "  CONST_MU: ${CONST_MU}"
echo "  Q_FOR_H: ${Q_FOR_H}"
echo "  P_PROBABILITY (expression): ${TARGET_P_PROBABILITY_EXPR}"
echo "  GAMMA_PLUS / GAMMA_MINUS: ${TARGET_GAMMA_VALUE}"
echo "  X_MIN: ${X_MIN}, X_MAX: ${X_MAX}"
echo "  N_POINTS_GRID: ${N_POINTS_GRID}"
echo "  NUM_ITERATIONS: ${NUM_ITERATIONS}"
echo "  OMP_NUM_THREADS: ${OMP_NUM_THREADS}"
echo "---------------------------------------------------------------------"

mkdir -p "${RESULTS_DIR}" # Biztosítjuk, hogy a célkönyvtár létezzen

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
gamma_fn=$(echo "${GAMMA_PLUS}" | tr '.' 'p' | sed 's/-/neg/g') # Használhatjuk a GAMMA_PLUS-t
q_h_fn=$(echo "${Q_FOR_H}" | tr '.' 'p' | sed 's/-/neg/g') # Lehet negatív q? Valószínűleg nem.
p_prob_fn=$(echo "${P_PROBABILITY}" | tr '.' 'p' | sed 's/-/neg/g')
xmin_fn=$(echo "${X_MIN}" | tr '.' 'p' | sed 's/-/neg/g')
xmax_fn=$(echo "${X_MAX}" | tr '.' 'p' | sed 's/-/neg/g')
npoints_fn="${N_POINTS_GRID}"
niter_fn="${NUM_ITERATIONS}"

param_suffix_log_csv="TARGETED_S${sigma_fn}_M${mu_fn}_G${gamma_fn}_QH${q_h_fn}_P${p_prob_fn}_XMIN${xmin_fn}_XMAX${xmax_fn}_NPT${npoints_fn}_NIT${niter_fn}"

# Töröljük az előző futásból származó potenciális kimeneti fájlokat (ha ugyanaz a suffix)
rm -f f_functions_output.csv # A C++ ezt a nevet használja ideiglenesen
# A tails_ratio és extra_calc fájlok neveit a C++ generálja egyedien,
# így azokat nem kell itt törölni, ha a C++ kód a paramétereket beépíti a nevükbe.
# Ha a C++ fix neveket használna ezekre is, akkor azokat is törölni kellene itt.

log_filename="log_${param_suffix_log_csv}.txt"
log_path="${RESULTS_DIR}/${log_filename}"
rm -f "${log_path}" # Töröljük az előző logot, ha létezik

echo "Running solver... Saving log to: ${log_path}"
if "${SOLVER_PROGRAM}" >"${log_path}" 2>&1; then
    echo "Solver finished successfully."
else
    echo "ERROR: Solver program returned an error! Check log: ${log_path}"
    # Itt lehetne `exit 1` ha hiba esetén meg akarunk állni.
fi

# Az f_functions_output.csv átnevezése és mozgatása
if [ -f "f_functions_output.csv" ]; then
    output_csv_filename="f_out_${param_suffix_log_csv}.csv"
    destination_csv_path="${RESULTS_DIR}/${output_csv_filename}"
    mv "f_functions_output.csv" "${destination_csv_path}"
    echo "Detailed CSV result saved to: ${destination_csv_path}"
else
    echo "WARNING: f_functions_output.csv was not created by the solver. Check log: ${log_path}"
fi

# A C++ által generált tails_ratio és extra_calc fájlok mozgatása
# Ezeknek egyedi nevük kell legyen a C++ kódban a paraméterek alapján,
# hogy ne írják felül egymást, ha a szkriptet többször futtatod más célzott paraméterekkel.
echo ""
echo "Moving generated .txt files (tails_ratio and extra_calc) to ${RESULTS_DIR} and tails_ratios/ ..."
# Először a célkönyvtárba, hogy megmaradjanak a logok mellett
find . -maxdepth 1 -type f -name 'tails_ratio_Q*.txt' -exec mv -t "${RESULTS_DIR}/" {} + 2>/dev/null
find . -maxdepth 1 -type f -name 'extra_calc_Q*.txt' -exec mv -t "${RESULTS_DIR}/" {} + 2>/dev/null

# Majd másolatot (vagy az eredetit, ha csak ott kell) a tails_ratios-ba
# Ha a C++ a tails_ratios-ba ír, akkor ez a rész nem kell.
# Ha a C++ a gyökérbe ír (ahogy az előző példákban feltételeztük):
cp "${RESULTS_DIR}/tails_ratio_Q"*.txt tails_ratios/ 2>/dev/null || echo "No tails_ratio files to copy to tails_ratios/ or error during copy."
cp "${RESULTS_DIR}/extra_calc_Q"*.txt tails_ratios/ 2>/dev/null || echo "No extra_calc files to copy to tails_ratios/ or error during copy."


echo "--- Targeted simulation run finished ---"

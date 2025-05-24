#!/bin/bash

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

OMP_THREADS=12
NUM_ITER=1000 # As set in the C++ code, but this could also be an env var

# Fixed SIGMA value
FIXED_SIGMA="0.05"

# GAMMA values (in basis points, will be converted)
# 1 bp, 5 bp, 10 bp, 30 bp, 100 bp
GAMMA_BP_VALUES=("1" "5" "10" "30" "100")

# Q_FOR_H values (from your example, comma corrected to dot)
Q_VALUES=("0.0" "0.2" "0.4" "0.5" "0.6" "0.8")


# P_PROBABILITY values (as expressions, calculated with bc)
# Or perhaps P = 1/20 = 0.05.
# The expressions below are kept as provided, but ensure they result in valid probabilities (0 to 1).
# For `bc`, `scale` must be set for floating point division.
P_EXPRESSIONS=(
    "1/(10*60)"       # approx. 0.001666
    "1/(2*60)"        # approx. 0.008333
    "1/12"            # approx. 0.083333
    "1/2"             # 0.5
    "0.05"            # Explicitly P=0.05 for testing the python script's filter
    "1/20"            # This is also 0.05
    "1/100"           # P=0.01
)


# --- Script Start ---
echo "Starting simulation series..."
mkdir -p "${RESULTS_DIR}"

if [ ! -f "${SOLVER_PROGRAM}" ]; then
    echo "Error: Solver program (${SOLVER_PROGRAM}) not found!"
    exit 1
fi
if [ ! -x "${SOLVER_PROGRAM}" ]; then
    echo "Error: Solver program (${SOLVER_PROGRAM}) not executable!"
    exit 1
fi

# Calculate MU based on SIGMA (using bc for floating point calculation)
# scale=10 gives, for example, 10 decimal places of precision
CALCULATED_MU=$(echo "scale=10; ${FIXED_SIGMA} * ${FIXED_SIGMA} / 2" | bc)

total_sims=$(( ${#GAMMA_BP_VALUES[@]} * ${#Q_VALUES[@]} * ${#P_EXPRESSIONS[@]} ))
current_sim=0

for Q_H in "${Q_VALUES[@]}"; do
  for GAMMA_BP in "${GAMMA_BP_VALUES[@]}"; do
    # Calculate GAMMA value (bp / 10000)
    CURRENT_GAMMA=$(echo "scale=10; ${GAMMA_BP} / 10000" | bc)
        for P_EXPR in "${P_EXPRESSIONS[@]}"; do
            # Calculate P_PROBABILITY
            # Ensure `bc` is available and `scale` is set for division
            CURRENT_P_PROB=$(echo "scale=10; ${P_EXPR}" | bc)

            # Validate P_PROBABILITY (must be between 0 and 1, exclusive of 0 for 1/P)
            is_p_valid=$(echo "${CURRENT_P_PROB} > 0 && ${CURRENT_P_PROB} < 1" | bc -l)
            if [ "$is_p_valid" -ne 1 ]; then
                echo ""
                echo "--- Skipping simulation due to invalid P_PROBABILITY ---"
                echo "Expression: ${P_EXPR} -> Calculated P: ${CURRENT_P_PROB}"
                echo "P_PROBABILITY must be > 0 and < 1."
                continue
            fi


            current_sim=$((current_sim + 1))
            echo ""
            echo "--- Running Simulation ${current_sim}/${total_sims} ---"

            # Create filename-friendly strings from parameters
            sigma_fn=$(echo "${FIXED_SIGMA}" | tr '.' 'p' | sed 's/-/neg/g')
            mu_fn=$(echo "${CALCULATED_MU}" | tr '.' 'p' | sed 's/-/neg/g')
            gamma_fn=$(echo "${CURRENT_GAMMA}" | tr '.' 'p' | sed 's/-/neg/g')
            q_h_fn=$(echo "${Q_H}" | tr '.' 'p')
            p_prob_fn=$(echo "${CURRENT_P_PROB}" | tr '.' 'p' | sed 's/-/neg/g') # Replace . with p for filename

            # This param_suffix is for the log and detailed CSV, not for the `tails_ratio_...txt`
            param_suffix_log_csv="S${sigma_fn}_M${mu_fn}_G${gamma_fn}_QH${q_h_fn}_P${p_prob_fn}"
            echo "Parameters: SIGMA=${FIXED_SIGMA}, MU=${CALCULATED_MU}, GAMMA=${CURRENT_GAMMA}, Q_H=${Q_H}, P_PROB=${CURRENT_P_PROB}"
            echo "Log/CSV Suffix: ${param_suffix_log_csv}"


            export OMP_NUM_THREADS="${OMP_THREADS}"
            export CONST_SIGMA="${FIXED_SIGMA}"
            export CONST_MU="${CALCULATED_MU}"
            export Q_FOR_H="${Q_H}"
            export P_PROBABILITY="${CURRENT_P_PROB}"
            export GAMMA_PLUS="${CURRENT_GAMMA}"  # GAMMA_PLUS and GAMMA_MINUS are the same
            export GAMMA_MINUS="${CURRENT_GAMMA}"
            # NUM_ITERATIONS is set in C++ main, but can be exported if C++ reads it via get_double_from_env
            # export NUM_ITERATIONS="${NUM_ITER}"

            # The C++ program is expected to create 'f_functions_output.csv',
            # 'tails_ratio_Q...txt', and 'extra_calc_Q...txt' in the current directory.
            # Clean up specific output files from previous runs in the loop
            rm -f f_functions_output.csv
            rm -f tails_ratio_Q*.txt # remove any previous ratio files
            rm -f extra_calc_Q*.txt  # remove any previous extra calc files


            log_filename="log_${param_suffix_log_csv}.txt"
            log_path="${RESULTS_DIR}/${log_filename}"
            rm -f "${log_path}" # Clear previous log for this specific parameter set

            echo "Running solver... Saving log to: ${log_path}"
            if "${SOLVER_PROGRAM}" >"${log_path}" 2>&1; then
                echo "Solver finished successfully."
            else
                echo "ERROR: Solver program returned an error! Check log: ${log_path}"
                # Optionally continue to next simulation or exit
                # continue
            fi

            # Move the detailed CSV if it was created
            if [ -f "f_functions_output.csv" ]; then
                output_csv_filename="f_out_${param_suffix_log_csv}.csv"
                destination_csv_path="${RESULTS_DIR}/${output_csv_filename}"
                mv "f_functions_output.csv" "${destination_csv_path}"
                echo "Detailed CSV result saved to: ${destination_csv_path}"
            else
                echo "WARNING: f_functions_output.csv was not created by the solver. Check log: ${log_path}"
            fi

            # The C++ code directly creates files like "tails_ratio_Q<q>_P<p>_... .txt"
            # These will be moved at the end of all simulations.
            # We can check if they were created for this iteration for feedback
            found_ratio_file=$(ls tails_ratio_Q${q_h_fn}_P${p_prob_fn}_*.txt 2>/dev/null | wc -l)
            if [ "$found_ratio_file" -gt 0 ]; then
                echo "Tails ratio file(s) created for this iteration."
            else
                echo "WARNING: No tails_ratio_...txt file found for this iteration. Check log: ${log_path}"
            fi
             found_extra_calc_file=$(ls extra_calc_Q${q_h_fn}_P${p_prob_fn}_*.txt 2>/dev/null | wc -l)
            if [ "$found_extra_calc_file" -gt 0 ]; then
                echo "Extra calculation file(s) created for this iteration."
            else
                echo "WARNING: No extra_calc_...txt file found for this iteration. Check log: ${log_path}"
            fi

        done
    done
done

echo ""
echo "Moving all generated .txt files (tails_ratio and extra_calc) to tails_ratios/ directory..."
# Use find to avoid "argument list too long" if many files, and handle no-match gracefully
find . -maxdepth 1 -type f -name 'tails_ratio_Q*.txt' -exec mv {} tails_ratios/ \;
find . -maxdepth 1 -type f -name 'extra_calc_Q*.txt' -exec mv {} tails_ratios/ \;

echo "--- Simulation series finished ---"

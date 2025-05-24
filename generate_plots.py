import os
import re
import glob
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as mticker

# --- (Reference Table and Helper Functions remain the same) ---
# --- 1. Reference Table Data ---
REFERENCE_TABLE = {
    "10 min": {"1 bp": 96.7, "5 bp": 85.5, "10 bp": 74.7, "30 bp": 49.6, "100 bp": 22.8},
    "2 min":  {"1 bp": 92.9, "5 bp": 72.5, "10 bp": 56.9, "30 bp": 30.5, "100 bp": 11.6},
    "12 sec": {"1 bp": 80.7, "5 bp": 45.6, "10 bp": 29.5, "30 bp": 12.3, "100 bp": 4.0},
    "2 sec":  {"1 bp": 63.0, "5 bp": 25.4, "10 bp": 14.5, "30 bp": 5.4, "100 bp": 1.7}
}

# --- Helper function to map P_PROBABILITY to Δt string key ---
def p_to_delta_t_key(p_value_float):
    if p_value_float == 0: return None
    delta_t_seconds = 1.0 / p_value_float
    epsilon = 1e-3
    if math.isclose(delta_t_seconds, 10 * 60, rel_tol=epsilon): return "10 min"
    if math.isclose(delta_t_seconds, 2 * 60, rel_tol=epsilon): return "2 min"
    if math.isclose(delta_t_seconds, 12, rel_tol=epsilon): return "12 sec"
    if math.isclose(delta_t_seconds, 2, rel_tol=epsilon): return "2 sec"
    return None

# --- Helper function to map GP_float to gamma string key ---
def gp_to_gamma_key(gp_value_float):
    gamma_bp = gp_value_float * 10000
    epsilon = 1e-1
    if math.isclose(gamma_bp, 1, rel_tol=epsilon): return "1 bp"
    if math.isclose(gamma_bp, 5, rel_tol=epsilon): return "5 bp"
    if math.isclose(gamma_bp, 10, rel_tol=epsilon): return "10 bp"
    if math.isclose(gamma_bp, 30, rel_tol=epsilon): return "30 bp"
    if math.isclose(gamma_bp, 100, rel_tol=epsilon): return "100 bp"
    return None

# --- Function to parse filename ---
def parse_filename_and_get_ratio(filepath):
    filename = os.path.basename(filepath)
    match = re.match(
        r"tails_ratio_Q0p[0-9]+_"
        r"P(?P<p_str>[0-9p]+)_"
        r"GP(?P<gp_str>[0-9p]+)_"
        r"GM(?P<gm_str>[0-9p]+)_"
        r"CS(?P<cs_str>[0-9p]+)_"
        r"CM(?P<cm_str>[0-9p]+)_"
        r"NPT(?P<npt_str>[0-9]+)\.txt",
        filename
    )
    if not match: return None
    data = match.groupdict()
    try:
        p_float = float(data['p_str'].replace('p', '.'))
        gp_float = float(data['gp_str'].replace('p', '.'))
        npt_int = int(data['npt_str'])
        with open(filepath, 'r') as f:
            simulated_ratio_percent = float(f.read().strip())
        delta_t_key = p_to_delta_t_key(p_float)
        gamma_key = gp_to_gamma_key(gp_float)

        print(f"  Extracting from: {filename}")
        print(f"    P_float: {p_float:.6f} -> Δt key: {delta_t_key}")
        print(f"    GP_float: {gp_float:.6f} -> γ key: {gamma_key}")
        print(f"    NPT: {npt_int}")
        print(f"    Simulated Tails Ratio (%): {simulated_ratio_percent:.5f}")

        if delta_t_key and gamma_key:
            return {
                "filename": filename, "p_float": p_float, "gp_float": gp_float,
                "npt": npt_int, "simulated_ratio": simulated_ratio_percent,
                "delta_t_key": delta_t_key, "gamma_key": gamma_key,
            }
        else:
            print(f"    Could not map P/GP to reference table keys. Skipping.")
            return None
    except (ValueError, IOError) as e:
        print(f"  Error processing/reading {filename}: {e}. Skipping.")
        return None


# --- Main processing logic ---
def main():
    results_dir = "tails_ratios/"
    plot_output_dir = "analysis_plots"
    os.makedirs(plot_output_dir, exist_ok=True)
    consolidated_plot_filename = os.path.join(plot_output_dir, "all_series_absolute_errors_readable.pdf") # New filename

    all_results = []
    print(f"Searching for results in directory: {results_dir}\n" + "="*30)
    found_files_count = 0
    parsed_files_count = 0
    for filepath in glob.glob(os.path.join(results_dir, "tails_ratio_Q0p*.txt")):
        found_files_count += 1
        parsed_data = parse_filename_and_get_ratio(filepath)
        if parsed_data:
            all_results.append(parsed_data)
            parsed_files_count +=1
        print("-" * 20)

    print(f"\n{found_files_count} files found matching pattern.")
    print(f"{parsed_files_count} files successfully parsed and mapped for plotting.\n" + "="*30)

    if not all_results:
        print("No results to plot.")
        return

    grouped_results = {}
    for res in all_results:
        key = (res['delta_t_key'], res['gamma_key'])
        if key not in grouped_results: grouped_results[key] = []
        grouped_results[key].append(res)

    # --- 4. Plotting ---
    print("\nStarting plot generation...\n" + "="*30)

    # Increased font sizes for better readability
    label_fontsize = 14
    tick_fontsize = 12
    legend_fontsize = 10 # Increased legend font size
    marker_size = 6 # Slightly larger markers
    line_width = 1.8 # Slightly thicker lines

    delta_t_order = ["2 sec", "12 sec", "2 min", "10 min"]
    gamma_order = ["1 bp", "5 bp", "10 bp", "30 bp", "100 bp"]

    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']
    delta_t_colors = {key: colors[i % len(colors)] for i, key in enumerate(delta_t_order)}
    
    linestyles = ['-', '--', '-.', ':', (0, (3, 1, 1, 1))] # Simplified custom linestyle
    gamma_linestyles = {key: linestyles[i % len(linestyles)] for i, key in enumerate(gamma_order)}
    
    sorted_keys = sorted(
        grouped_results.keys(),
        key=lambda k: (
            delta_t_order.index(k[0]) if k[0] in delta_t_order else float('inf'),
            gamma_order.index(k[1]) if k[1] in gamma_order else float('inf')
        )
    )
    
    # Increased figure size for more space
    fig, ax = plt.subplots(figsize=(16, 10)) 
    all_npt_values = set()

    for key_tuple in sorted_keys:
        delta_t_key, gamma_key = key_tuple
        series_results = grouped_results[key_tuple]
        
        series_results.sort(key=lambda x: x['npt'])
        npt_values = [r['npt'] for r in series_results]
        simulated_ratios = [r['simulated_ratio'] for r in series_results]
        
        all_npt_values.update(npt_values)

        try:
            reference_value = REFERENCE_TABLE[delta_t_key][gamma_key]
            print(f"  Processing series for Δt={delta_t_key}, γ={gamma_key}: Ref={reference_value:.1f}%")
        except KeyError:
            print(f"  Warning: No ref value for Δt={delta_t_key}, γ={gamma_key}. Skipping.")
            continue
            
        absolute_errors = [abs(sr - reference_value) for sr in simulated_ratios]
        print(f"    NPT: {npt_values} -> AbsErrors (%): {[f'{ae:.2f}' for ae in absolute_errors]}")

        plot_label = f"Δt={delta_t_key}, γ={gamma_key}"
        color = delta_t_colors.get(delta_t_key, 'k') 
        linestyle = gamma_linestyles.get(gamma_key, '-')

        ax.plot(npt_values, absolute_errors, 
                color=color,
                linestyle=linestyle,
                marker='o', # Consider distinct markers if lines are too similar
                markersize=marker_size,
                linewidth=line_width,
                label=plot_label)
        
    ax.set_xlabel("NPT (Number of Grid Points)", fontsize=label_fontsize, labelpad=10) # Added labelpad
    ax.set_ylabel("Absolute Error from Reference (%)", fontsize=label_fontsize, labelpad=10) # Added labelpad
    ax.grid(True, linestyle=':', alpha=0.6) # Slightly more subtle grid
    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize, pad=7) # Added pad for tick labels

    sorted_unique_npt = sorted(list(all_npt_values))
    if sorted_unique_npt:
        ax.set_xticks(sorted_unique_npt)
    ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter()) # Keep standard number format
    
    # Adjust rotation if many NPT labels to prevent overlap
    if len(sorted_unique_npt) > 8 : # Adjusted heuristic
         plt.setp(ax.get_xticklabels(), rotation=30, ha="right", rotation_mode="anchor")

    # Log scale for x-axis if NPT range is very large
    if sorted_unique_npt and len(sorted_unique_npt) > 1 and (max(sorted_unique_npt) / min(sorted_unique_npt) > 10):
        ax.set_xscale('log')
        ax.get_xaxis().set_major_formatter(mticker.ScalarFormatter()) # Ensure formatter is good for log

    # --- Legend Adjustment ---
    # Place legend in the upper right corner INSIDE the plot
    # Adjust ncol for better fit (e.g., 2 or 3 columns)
    legend = ax.legend(loc='upper right', ncol=2, fontsize=legend_fontsize, 
                       fancybox=True, framealpha=0.9, title="Series (Δt, γ)")
    legend.get_title().set_fontsize(legend_fontsize + 1) # Make legend title slightly larger

    # Add padding around the plot to ensure labels/legend don't get cut off
    plt.subplots_adjust(left=0.1, right=0.95, top=0.95, bottom=0.1)


    try:
        plt.savefig(consolidated_plot_filename, format='pdf', bbox_inches='tight', dpi=300) # Added dpi
        print(f"\nConsolidated absolute error plot saved to: {consolidated_plot_filename}")
    except Exception as e:
        print(f"Error saving consolidated plot {consolidated_plot_filename}: {e}")
    
    plt.close(fig)

if __name__ == "__main__":
    main()

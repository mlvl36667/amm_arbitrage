import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import os
import re
import glob # <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< ADD THIS LINE

# --- Configuration ---
# Automatically find the newest CSV file in the specified directory
# or allow manual specification.

RESULTS_DIR = "simulation_results_targeted_run/"
PLOT_OUTPUT_DIR = "analysis_plots"
os.makedirs(PLOT_OUTPUT_DIR, exist_ok=True)

# Functions to plot (e.g., f0, f10, f50, f99).
# If you want all, it will be very cluttered.
# Use an empty list to try and plot all available f_k columns.
# Or specify a range or specific indices, e.g., [0, 9, 49, 99] for f0, f9, f49, f99
# For NUM_ITERATIONS=1000, you might plot f0, f100, f500, f999
FUNCTIONS_TO_PLOT_INDICES = [0, 19, 29, 49, 99] # Example: f0, f99, ... f999

# Plot styling
label_fontsize = 16
tick_fontsize = 14
legend_fontsize = 12
line_width = 2.0
figure_size = (15, 9)


def find_latest_csv_file(directory):
    """Finds the most recently modified CSV file in the given directory."""
    list_of_files = glob.glob(os.path.join(directory, '*.csv'))
    if not list_of_files:
        return None
    latest_file = max(list_of_files, key=os.path.getctime)
    return latest_file

def plot_functions_from_csv(csv_filepath, function_indices_to_plot):
    """
    Reads a CSV file and plots specified f_k functions against x_value.
    """
    if not os.path.exists(csv_filepath):
        print(f"Error: CSV file not found at {csv_filepath}")
        return

    print(f"--- Reading data from: {csv_filepath} ---")
    try:
        df = pd.read_csv(csv_filepath)
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return

    if 'x_value' not in df.columns:
        print("Error: 'x_value' column not found in the CSV.")
        return

    x_values = df['x_value']
    
    # Determine available f_k columns
    fk_columns = [col for col in df.columns if re.match(r"f\d+", col)]
    
    if not fk_columns:
        print("Error: No function columns (f0, f1, ...) found in the CSV.")
        return

    if not function_indices_to_plot: # If empty, attempt to plot all
        columns_to_plot = fk_columns
        print(f"Plotting all available {len(fk_columns)} functions. This might be cluttered.")
    else:
        columns_to_plot = []
        for idx in function_indices_to_plot:
            col_name = f"f{idx}"
            if col_name in df.columns:
                columns_to_plot.append(col_name)
            else:
                print(f"Warning: Column {col_name} not found in CSV. Skipping.")
    
    if not columns_to_plot:
        print("No valid function columns selected or available for plotting.")
        return

    print(f"Plotting functions: {', '.join(columns_to_plot)}")

    fig, ax = plt.subplots(figsize=figure_size)

    # Plotting each selected function
    for col_name in columns_to_plot:
        ax.plot(x_values, df[col_name], linewidth=line_width, label=col_name)

    ax.set_xlabel("x_value", fontsize=label_fontsize, labelpad=10)
    ax.set_ylabel("Function Value f_k(x)", fontsize=label_fontsize, labelpad=10)
    ax.grid(True, linestyle=':', alpha=0.7)
    ax.tick_params(axis='both', which='major', labelsize=tick_fontsize, pad=7)

    # Legend
    legend = ax.legend(loc='upper right', ncol=1, fontsize=legend_fontsize,
                       fancybox=True, framealpha=0.9)
    # If too many legend entries, consider loc='best' or placing outside
    if len(columns_to_plot) > 10 : # Heuristic for many legend entries
        legend = ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), ncol=1,
                           fontsize=legend_fontsize, fancybox=True, framealpha=0.9,
                           title="Functions")
        legend.get_title().set_fontsize(legend_fontsize)


    # X-axis ticks - make them readable, especially if x_values are dense
    # ax.xaxis.set_major_locator(mticker.MaxNLocator(nbins=10, prune='both')) # Example: at most 10 ticks
    # For y-axis, default locator is often fine, or use similar MaxNLocator
    # ax.yaxis.set_major_locator(mticker.MaxNLocator(nbins=8))


    # Try to make the plot look nice
    plt.subplots_adjust(left=0.1, right=0.85 if len(columns_to_plot) > 10 else 0.95,
                        top=0.95, bottom=0.1) # Adjust right margin for external legend


    # Sanitize filename for saving
    base_csv_filename = os.path.basename(csv_filepath)
    plot_filename_base = os.path.splitext(base_csv_filename)[0]
    pdf_filename = os.path.join(PLOT_OUTPUT_DIR, f"plot_{plot_filename_base}.pdf")

    try:
        plt.savefig(pdf_filename, format='pdf', bbox_inches='tight', dpi=300)
        print(f"--- Plot saved to: {pdf_filename} ---")
    except Exception as e:
        print(f"Error saving plot {pdf_filename}: {e}")
    
    plt.show() # Display the plot (optional, can be commented out if only saving)
    plt.close(fig)


if __name__ == "__main__":
    # Option 1: Manually specify the CSV file
    # csv_file_to_plot = os.path.join(RESULTS_DIR, "f_out_TARGETED_S0p02124475_M0p03712680_G0p003_QH0p02319134_Pp0833333333_XMINneg0p02_XMAX0p02_NPT200_NIT100.csv")
    
    # Option 2: Automatically find the latest CSV in RESULTS_DIR
    latest_csv = find_latest_csv_file(RESULTS_DIR)
    
    if latest_csv:
        plot_functions_from_csv(latest_csv, FUNCTIONS_TO_PLOT_INDICES)
    else:
        print(f"No CSV files found in {RESULTS_DIR}. Please specify a file manually or check the directory.")

    # Example of plotting a specific file if auto-find fails or you want to override:
    # specific_file = "path/to/your/file.csv"
    # if os.path.exists(specific_file):
    #     plot_functions_from_csv(specific_file, FUNCTIONS_TO_PLOT_INDICES)

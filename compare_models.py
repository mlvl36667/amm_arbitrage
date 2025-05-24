#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import seaborn as sns
from matplotlib.ticker import MaxNLocator, FuncFormatter

# --- Configuration ---
RESULTS_DIR = "simulation_results_targeted_run"

FILE_WITH_Q_NAME = "f_out_TARGETED_S0p02124475_M0p03712680_G0p003_QH0p02319134_Pp0833333333_XMINneg0p02_XMAX0p02_NPT200_NIT5000.csv"
FILE_WITHOUT_Q_NAME = "f_out_TARGETED_S0p02124475_M0p03712680_G0p003_QH0p0_Pp0833333333_XMINneg0p02_XMAX0p02_NPT200_NIT5000.csv"

X_COLUMN_NAME = 'x_value'

# Plot styling
AXIS_LABEL_FONTSIZE = 14
TICK_LABEL_FONTSIZE = 11
LEGEND_FONTSIZE = 12
LINE_WIDTH = 2.0
FIG_SIZE = (10, 6)

def extract_q_from_filename(filename):
    match = re.search(r"_QH([^_\/]+?)_", filename)
    if match:
        q_str_raw = match.group(1)
        q_str_clean = q_str_raw.replace('p', '.').replace('neg', '-')
        try:
            return float(q_str_clean)
        except ValueError:
            return q_str_raw
    return "Unknown Q"

def get_last_f_col_name(columns):
    f_cols = sorted([col for col in columns if re.match(r'^f\d+$', col)],
                    key=lambda x: int(re.search(r'\d+', x).group()))
    if f_cols:
        return f_cols[-1]
    return None

def load_data(filepath, filename_for_log):
    try:
        data = pd.read_csv(filepath, low_memory=False)
        actual_last_f_col_name = get_last_f_col_name(data.columns)
        
        if X_COLUMN_NAME in data.columns and actual_last_f_col_name:
            x_data = data[X_COLUMN_NAME]
            f_data = data[actual_last_f_col_name]
            return x_data, f_data, actual_last_f_col_name
        else:
            # ... (hibaüzenetek, mint korábban)
            raise pd.errors.ParserError("Required columns not found with headers, trying without header.")

    except (KeyError, ValueError, pd.errors.ParserError, FileNotFoundError) as e:
        try:
            data = pd.read_csv(filepath, header=None, low_memory=False)
            x_data = data.iloc[:, 0]
            f_data = data.iloc[:, -1]
            return x_data, f_data, "Last Column"
        except Exception as e_fallback:
            print(f"ERROR: Failed to load {filename_for_log} even with header=None. Error: {e_fallback}")
            return None, None, None

def plot_comparison():
    filepath_with_q = os.path.join(RESULTS_DIR, FILE_WITH_Q_NAME)
    filepath_without_q = os.path.join(RESULTS_DIR, FILE_WITHOUT_Q_NAME)

    if not os.path.isdir(RESULTS_DIR):
        print(f"ERROR: Results directory not found: {RESULTS_DIR}")
        return

    x_data_with_q, f_data_with_q, _ = load_data(filepath_with_q, FILE_WITH_Q_NAME)
    if x_data_with_q is None: return

    x_data_without_q, f_data_without_q, _ = load_data(filepath_without_q, FILE_WITHOUT_Q_NAME)
    if x_data_without_q is None: return

    label_with_q = 'With Q Jump'
    label_without_q = 'Without Q Jump'

    # === Stílus és rács kezelése ===
    chosen_style = "white" # Választható: "white", "whitegrid", "darkgrid", stb. vagy None az alapértelmezett Matplotlib stílushoz

    if chosen_style:
        try:
            sns.set_theme(style=chosen_style)
            # print(f"INFO: Seaborn style '{chosen_style}' applied.")
        except Exception as e_sns:
            # print(f"INFO: Could not apply Seaborn style '{chosen_style}' (Error: {e_sns}). Falling back.")
            try:
                plt.style.use(chosen_style if 'seaborn' not in chosen_style else 'ggplot') # Próbáljuk meg Matplotlib-ként, ha nem seaborn-specifikus
            except:
                # print("INFO: Fallback Matplotlib style also failed. Using library default.")
                pass # Matplotlib alapértelmezett stílusát használja
    # ==============================

    fig, ax = plt.subplots(figsize=FIG_SIZE)

    ax.plot(x_data_with_q, f_data_with_q, label=label_with_q, linewidth=LINE_WIDTH, alpha=0.8)
    ax.plot(x_data_without_q, f_data_without_q, label=label_without_q, linewidth=LINE_WIDTH, linestyle='--', alpha=0.8)

    ax.set_xlabel('Value of Logarithmic Misprice', fontsize=AXIS_LABEL_FONTSIZE)
    ax.set_ylabel('Value of Density', fontsize=AXIS_LABEL_FONTSIZE)

    ax.set_xlim(-0.015, 0.025)

    ax.yaxis.set_major_locator(MaxNLocator(nbins=6, prune='both'))
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5, prune='both'))

    def x_axis_formatter(x, pos):
        return f'{x:.3f}'
    ax.xaxis.set_major_formatter(FuncFormatter(x_axis_formatter))
    
    ax.tick_params(axis='both', which='major', labelsize=TICK_LABEL_FONTSIZE)
    
    legend = ax.legend(fontsize=LEGEND_FONTSIZE, frameon=True, loc='best')
    legend.get_frame().set_edgecolor('darkgray')
    legend.get_frame().set_alpha(0.9)

    # === Rács hozzáadása, ha a választott stílus nem tartalmazza (pl. "white") ===
    if chosen_style == "white": # Vagy más rács nélküli stílusok esetén is hozzáadhatod
        ax.grid(True, linestyle=':', alpha=0.5, color='lightgray')
    # Ha a chosen_style pl. "whitegrid", akkor a rácsot a stílus már tartalmazza.
    # ==========================================================================

    try:
        fig.tight_layout(pad=0.5)
    except Exception as e_layout:
        print(f"Warning: tight_layout() raised an error: {e_layout}. Plot may not be optimally arranged.")
    
    plot_filename_pdf = "final_f_comparison_q_vs_no_q.pdf"
    try:
        plt.savefig(plot_filename_pdf, format='pdf', dpi=300, bbox_inches='tight')
        print(f"Plot saved to {plot_filename_pdf}")
    except Exception as e:
        print(f"Error saving plot to PDF: {e}")
        
    plt.show()

if __name__ == "__main__":
    plot_comparison()

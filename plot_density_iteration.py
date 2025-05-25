import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.ticker as mticker # Opcionális a tengelyformázáshoz

# Stílusbeállítások a jobb olvashatóságért (opcionális, a dokumentum elejére)
plt.rcParams.update({
    'font.size': 12, # Alapértelmezett betűméret növelése
    'axes.labelsize': 14, # Tengelyfeliratok
    'axes.titlesize': 16, # Címek (subplot)
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 10,
    'figure.titlesize': 18 # Főcím (suptitle)
})


def parse_value_from_param_str(param_str_with_prefix, default_precision=6, sci_limit=1e-4):
    match_val = re.match(r"^[A-Za-z]+([0-9neg]+p?[0-9]*)$", param_str_with_prefix)
    if not match_val:
        try:
            float_val = float(param_str_with_prefix)
        except ValueError:
            return param_str_with_prefix
    else:
        value_part_str = match_val.group(1)
        is_negative = False
        if value_part_str.startswith("neg"):
            is_negative = True
            value_part_str = value_part_str[3:]
        numeric_str = value_part_str.replace('p', '.') if 'p' in value_part_str else value_part_str
        try:
            float_val = float(numeric_str)
            if is_negative: float_val *= -1
        except ValueError:
            return param_str_with_prefix

    if float_val.is_integer(): return str(int(float_val))
    if (abs(float_val) < sci_limit and float_val != 0) or abs(float_val) > 1/sci_limit:
        return f"{float_val:.2e}" # Kicsit több pontosság a tudományosnál
    formatted = f"{float_val:.{default_precision}f}".rstrip('0').rstrip('.')
    return formatted if formatted else "0"


def generate_fixed_param_title_and_q(q_str_for_qh):
    # NPT és NIT javítva a megadott fájlnevek alapján
    s_disp = parse_value_from_param_str("S0p05")
    m_disp = parse_value_from_param_str("Mp0012500000")
    g_disp = parse_value_from_param_str("G0p003")
    p_disp = parse_value_from_param_str("Pp0833333333")
    npt_disp = "801" # A megadott fájlnevekben NPT801 van
    nit_disp = "1000" # A megadott fájlnevekben NIT100 van
    
    common_title_parts = [
        f"σ={s_disp}", f"μ={m_disp}", f"γ={g_disp}", 
        f"P(jump)={p_disp}", f"NPT={npt_disp}", f"NIT={nit_disp}"
    ]
    # Két sorba tördeljük a jobb olvashatóságért, ha túl hosszú lenne
    if len(", ".join(common_title_parts)) > 60 : # Tetszőleges határ
        common_title = "Common Params: " + ", ".join(common_title_parts[:3]) + "\n" + ", ".join(common_title_parts[3:])
    else:
        common_title = "Common Params: " + ", ".join(common_title_parts)
        
    q_value_for_subplot = parse_value_from_param_str(f"QH{q_str_for_qh}")
    return common_title, q_value_for_subplot


def plot_specific_q_series_2x2(source_directory, q_str_parts, output_image_file="q_series_plot.pdf"):
    base_filename_prefix = "f_out_S0p05_Mp0012500000_G0p003_QH"
    # NPT és NIT javítva a fájlnévben
    base_filename_suffix = "_Pp0833333333_XMINneg0p02_XMAX0p02_NPT801_NIT1000.csv"

    if len(q_str_parts) != 4:
        print("Hiba: Pontosan 4 Q string részre van szükség a 2x2-es plot-hoz.")
        return

    # Nagyobb ábra, hogy a feliratok elférjenek
    fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharex=True, sharey=False) 
    axes = axes.flatten()

    common_suptitle, _ = generate_fixed_param_title_and_q(q_str_parts[0])
#    fig.suptitle(f"Value of Density Functions (f)\n{common_suptitle}", y=0.99) # y kicsit feljebb a suptitle-nek

    for i, q_str in enumerate(q_str_parts):
        ax = axes[i]
        csv_filename = f"{base_filename_prefix}{q_str}{base_filename_suffix}"
        file_path = os.path.join(source_directory, csv_filename)
        _, q_val_for_subplot_title = generate_fixed_param_title_and_q(q_str)

        print(f"Plotting: {csv_filename} (Q={q_val_for_subplot_title}) -> subplot {i+1}")

        if not os.path.exists(file_path):
            ax.text(0.5, 0.5, f"File not found:\n{csv_filename}", ha='center', va='center', 
                    transform=ax.transAxes, color='red', wrap=True, fontsize=10)
            print(f"  ERROR: File not found: '{file_path}'")
            continue
            
        try:
            df = pd.read_csv(file_path)
            if 'x_value' not in df.columns:
                ax.text(0.5, 0.5, "'x_value' column missing", ha='center', va='center', transform=ax.transAxes)
                continue

            x_values = df['x_value']
            f_columns_all = sorted([col for col in df.columns if re.match(r'^f\d+$', col)],
                                   key=lambda name: int(name[1:]))
            
            if not f_columns_all:
                ax.text(0.5, 0.5, "No 'f' columns found", ha='center', va='center', transform=ax.transAxes)
                continue

            f_columns_to_plot = [f_columns_all[j] for j in range(0, len(f_columns_all), 100)]
            if not f_columns_to_plot:
                ax.text(0.5, 0.5, "No functions to plot (after filtering)", ha='center', va='center', transform=ax.transAxes)
                continue

            # Színpaletta (opcionális, ha az alapértelmezett nem elég kontrasztos)
            # colors = plt.cm.get_cmap('tab10', len(f_columns_to_plot)) # Vagy 'viridis', 'plasma'

            for idx, f_col in enumerate(f_columns_to_plot):
                # ax.plot(x_values, df[f_col], label=f_col, linewidth=1.2, color=colors(idx))
                ax.plot(x_values, df[f_col], label=f_col, linewidth=1.2) # Vastagabb vonalak
            
            combined_f_values = df[f_columns_to_plot].copy()
            valid_x_indices = combined_f_values.notna().any(axis=1)
            if valid_x_indices.any():
                first_valid_idx = valid_x_indices.idxmax()
                last_valid_idx = valid_x_indices[::-1].idxmax()
                min_x_plot = df.loc[first_valid_idx, 'x_value']
                max_x_plot = df.loc[last_valid_idx, 'x_value']
                
                # Kis puffer hozzáadása, hogy a görbék ne érjenek a széléhez
                plot_range = max_x_plot - min_x_plot
                if plot_range > 0: # Csak ha van értelmes tartomány
                    ax.set_xlim([min_x_plot - 0.05 * plot_range, max_x_plot + 0.05 * plot_range])
                elif len(x_values) > 0 : # Ha csak egy pont van, vagy minden x ugyanaz
                    ax.set_xlim([x_values.min() - 0.1, x_values.max() + 0.1]) # Kis fix puffer
            
            ax.set_title(f"q = {q_val_for_subplot_title}") # A globális rcParams.update() már növelte a méretet
            # Az X tengely feliratát csak az alsó sor subplotjaira tesszük ki, ha sharex=True
            if i >= 2 : # Alsó sor (index 2 és 3)
                 ax.set_xlabel("Logarithmic Misprice")
            if i % 2 == 0:
                ax.set_ylabel("Value of Density")
            
            # Tengely tick formázás (opcionális, ha finomhangolni kell)
            # ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%.2f'))
            # ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%.1e'))
            
            ax.grid(True, linestyle=':', alpha=0.5) # Halványabb, pontozott rács
            
            if len(f_columns_to_plot) <= 12: # Kicsit több vonalat engedünk a legendában
                 ax.legend(loc='best', frameon=False, ncol=1 if len(f_columns_to_plot) <=6 else 2)

        except pd.errors.EmptyDataError:
            ax.text(0.5, 0.5, "File is empty", ha='center', va='center', transform=ax.transAxes)
        except Exception as e:
            ax.text(0.5, 0.5, f"Error:\n{e}", ha='center', va='center', transform=ax.transAxes, wrap=True, fontsize=9)
            print(f"  Error processing '{csv_filename}': {e}")

    plt.tight_layout(rect=[0, 0.02, 1, 0.95]) # Finomhangolás a suptitle és a tengelyfeliratok helyének
    plt.savefig(output_image_file, bbox_inches='tight', dpi=150) # dpi növelése a jobb felbontásért
    plt.close(fig)
    print(f"\nDone! 2x2 plot saved to '{output_image_file}'.")

if __name__ == "__main__":
    simulation_dir = "simulation_results/"
    q_identifier_strings = ["0p0", "0p2", "0p4", "0p8"] 
    output_file = "f_functions_Q_series_2x2_improved.pdf"
    plot_specific_q_series_2x2(simulation_dir, q_identifier_strings, output_file)

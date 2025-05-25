import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 16,
    'xtick.labelsize': 11,
    'ytick.labelsize': 11,
    'legend.fontsize': 9, # Kicsit kisebb a sűrűbb legendához
    'figure.titlesize': 18
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
            # Ha a numerikus rész önmagában nem szám (pl. üres string a prefix után),
            # akkor az eredeti teljes stringet adjuk vissza.
            return param_str_with_prefix
    
    # Float formázása, hibakezeléssel, ha float_val nem jött létre helyesen
    try:
        if float_val.is_integer(): return str(int(float_val))
        if (abs(float_val) < sci_limit and float_val != 0) or abs(float_val) > 1/sci_limit:
            return f"{float_val:.2e}"
        formatted = f"{float_val:.{default_precision}f}".rstrip('0').rstrip('.')
        return formatted if formatted else "0"
    except AttributeError: # Ha float_val nem rendelkezik is_integer stb. metódussal
        return param_str_with_prefix


def generate_fixed_param_title_and_q(q_str_for_qh):
    # Paraméterek a Bash log alapján
    s_disp = parse_value_from_param_str("S0p05")           # CONST_SIGMA: 0.05
    m_disp = parse_value_from_param_str("Mp0012500000")    # CONST_MU: .0012500000
    g_disp = parse_value_from_param_str("G0p003")          # GAMMA: 0.003
    p_disp = parse_value_from_param_str("Pp0833333333")    # P_PROBABILITY: .0833333333
    npt_disp = "801"                                      # N_POINTS_GRID: 801
    nit_disp = "1000"                                     # NUM_ITERATIONS: 1000
    # XMIN/XMAX a fájlnévből jön, de a suptitle-be is betehetjük, ha kell
    # xmin_disp = parse_value_from_param_str("XMINneg0p01")
    # xmax_disp = parse_value_from_param_str("XMAX0p01")

    common_title_parts = [
        f"σ={s_disp}", f"μ={m_disp}", f"γ={g_disp}", 
        f"P(jump)={p_disp}", f"NPT={npt_disp}", f"NIT={nit_disp}"
    ]
    
    title_line1 = ", ".join(common_title_parts[:3])
    title_line2 = ", ".join(common_title_parts[3:])
    common_title = f"{title_line1}\n{title_line2}"
        
    q_value_for_subplot = parse_value_from_param_str(f"QH{q_str_for_qh}")
    return common_title, q_value_for_subplot


def plot_specific_q_series_2x2(source_directory, q_str_parts, output_image_file="q_series_plot.pdf"):
    # Fájlnév részek a Bash log alapján
    base_filename_prefix = "f_out_S0p05_Mp0012500000_G0p003_QH"
    base_filename_suffix = "_Pp0833333333_XMINneg0p01_XMAX0p01_NPT801_NIT1000.csv"

    if len(q_str_parts) != 4:
        print("Hiba: Pontosan 4 Q string részre van szükség a 2x2-es plot-hoz.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharex=True, sharey=False) 
    axes = axes.flatten()

    common_suptitle, _ = generate_fixed_param_title_and_q(q_str_parts[0]) # Az első Q-val generáljuk a közös részt
#    fig.suptitle(f"Value of Density Functions (f) for Varying q\n{common_suptitle}", y=0.98)

    for i, q_str in enumerate(q_str_parts):
        ax = axes[i]
        csv_filename = f"{base_filename_prefix}{q_str}{base_filename_suffix}"
        file_path = os.path.join(source_directory, csv_filename)
        _, q_val_for_subplot_title = generate_fixed_param_title_and_q(q_str)

        print(f"Plotting: {csv_filename} (q={q_val_for_subplot_title}) -> subplot {i+1}")

        if not os.path.exists(file_path):
            ax.text(0.5, 0.5, f"File not found:\n{csv_filename}", ha='center', va='center', 
                    transform=ax.transAxes, color='red', wrap=True, fontsize=10)
            print(f"  ERROR: File not found: '{file_path}'")
            ax.set_title(f"q = {q_val_for_subplot_title}\n(FILE NOT FOUND)", fontsize=12, color='red')
            continue
            
        try:
            df = pd.read_csv(file_path)
            if 'x_value' not in df.columns:
                ax.text(0.5, 0.5, "'x_value' column missing", ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f"q = {q_val_for_subplot_title}\n(NO x_value)", fontsize=12)
                continue

            x_values = df['x_value']
            f_columns_all = sorted([col for col in df.columns if re.match(r'^f\d+$', col)],
                                   key=lambda name: int(name[1:]))
            
            if not f_columns_all:
                ax.text(0.5, 0.5, "No 'f' columns found", ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f"q = {q_val_for_subplot_title}\n(NO f columns)", fontsize=12)
                continue
            
            # NPT801 -> f0...f800. Ha minden 100-adikat, az 8-9 vonal.
            # Ha a fv-ek száma NPT, akkor f_columns_all hossza NPT lesz.
            # Ha f0..f99 van, akkor a range(0, len(f_columns_all), 10) 10 vonalat ad.
            # Mivel NPT=801, az f függvények f0-tól f800-ig mehetnek.
            # Ha 10-15 vonalnál nem akarunk többet:
            num_f_cols = len(f_columns_all)
            step = max(1, num_f_cols // 10) # Kb. 10 vonal legyen, de legalább 1-es lépésköz
            if num_f_cols > 150: # Ha nagyon sok fv van, ritkítsunk jobban
                step = max(1, num_f_cols // 8)


            f_columns_to_plot = [f_columns_all[j] for j in range(0, num_f_cols, step)]
            if not f_columns_to_plot: # Elvileg nem fordulhat elő, ha van f_columns_all
                ax.text(0.5, 0.5, "No functions to plot (after filtering)", ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f"q = {q_val_for_subplot_title}\n(FILTERING ERROR)", fontsize=12)
                continue

            for idx, f_col in enumerate(f_columns_to_plot):
                ax.plot(x_values, df[f_col], label=f_col, linewidth=1.2)
            
            combined_f_values = df[f_columns_to_plot].copy()
            valid_x_indices = combined_f_values.notna().any(axis=1)
            if valid_x_indices.any():
                first_valid_idx = valid_x_indices.idxmax()
                last_valid_idx = valid_x_indices[::-1].idxmax()
                min_x_plot = df.loc[first_valid_idx, 'x_value']
                max_x_plot = df.loc[last_valid_idx, 'x_value']
                plot_range = max_x_plot - min_x_plot
                if plot_range > 1e-9: # Csak ha van értelmes tartomány (nem egyetlen pont)
                    ax.set_xlim([min_x_plot - 0.05 * plot_range, max_x_plot + 0.05 * plot_range])
                elif len(x_values) > 0 : 
                    # Ha a tartomány nulla (pl. minden érték ugyanaz), vagy kevés pont van,
                    # adjunk egy kis fix margót az x_min és x_max köré, ha azok eltérnek.
                    # Vagy hagyjuk az alapértelmezett matplotlib határokat.
                    # Most egy kis fix puffert adunk, ha min_x_plot és max_x_plot ugyanaz
                    if abs(x_values.min() - x_values.max()) < 1e-9 and len(x_values)>1: #Minden x ugyanaz
                         ax.set_xlim([x_values.min() - 0.01, x_values.max() + 0.01]) # Fix puffer
                    # Ha a min_x_plot és max_x_plot ugyanaz volt, de az x_values változatosabb,
                    # akkor a fentebbi xlim beállítás már jó lehetett.
            
            ax.set_title(f"q = {q_val_for_subplot_title}")
            if i >= 2 : 
                 ax.set_xlabel("Logarithmic Misprice (x)")
            if i % 2 == 0:
                ax.set_ylabel("Value of Density (f)")
            
            ax.grid(True, linestyle=':', alpha=0.5)
            
            if len(f_columns_to_plot) <= 15: # Legenda, ha nem túl sok vonal van
                 ax.legend(fontsize='small', loc='best', frameon=False, ncol=1 if len(f_columns_to_plot) <=7 else 2)

        except pd.errors.EmptyDataError:
            ax.text(0.5, 0.5, "File is empty", ha='center', va='center', transform=ax.transAxes)
            ax.set_title(f"q = {q_val_for_subplot_title}\n(FILE EMPTY)", fontsize=12)
        except Exception as e:
            ax.text(0.5, 0.5, f"Error:\n{e}", ha='center', va='center', transform=ax.transAxes, wrap=True, fontsize=9)
            ax.set_title(f"q = {q_val_for_subplot_title}\n(PROCESSING ERROR)", fontsize=12)
            print(f"  Error processing '{csv_filename}': {e}")

    plt.tight_layout(rect=[0, 0.02, 1, 0.95])
    plt.savefig(output_image_file, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"\nDone! 2x2 plot saved to '{output_image_file}'.")

if __name__ == "__main__":
    simulation_dir = "simulation_results/"
    # A QH paraméter string részei a Bash log alapján
    q_identifier_strings = ["0p0", "0p05", "0p1", "0p2"] 
    
    output_file = "f_functions_Q_series_NIT1000_NPT801.pdf" # Fájlnév frissítve

    plot_specific_q_series_2x2(simulation_dir, q_identifier_strings, output_file)

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def parse_value_from_param_str(param_str_with_prefix, default_precision=6, sci_limit=1e-4):
    """
    Kinyeri a numerikus értéket a paraméter stringből (pl. 'S0p05' -> 0.05).
    A param_str_with_prefix formája pl. "S0p05", "QH0p2", "Pp0833333333".
    """
    # Elválasztja a betűs prefixet a numerikus/p-s résztől
    match_val = re.match(r"^[A-Za-z]+([0-9neg]+p?[0-9]*)$", param_str_with_prefix)
    if not match_val:
        # Ha nincs prefix, vagy nem illeszkedik, próbáljuk meg sima floatként
        try:
            float_val = float(param_str_with_prefix) # Pl. ha csak "0.05" jönne
        except ValueError:
            return param_str_with_prefix # Vissza az eredetivel
    else:
        value_part_str = match_val.group(1) # Pl. "0p05", "0p2", "0833333333"
        
        is_negative = False
        if value_part_str.startswith("neg"):
            is_negative = True
            value_part_str = value_part_str[3:]

        if 'p' in value_part_str:
            numeric_str = value_part_str.replace('p', '.')
        else:
            numeric_str = value_part_str
        
        try:
            float_val = float(numeric_str)
            if is_negative:
                float_val *= -1
        except ValueError:
            return param_str_with_prefix # Nem sikerült float-tá alakítani

    # Float formázása
    if float_val.is_integer():
        return str(int(float_val))
    if (abs(float_val) < sci_limit and float_val != 0) or abs(float_val) > 1/sci_limit:
        return f"{float_val:.1e}"
    else:
        formatted = f"{float_val:.{default_precision}f}".rstrip('0').rstrip('.')
        return formatted if formatted else "0"


def generate_fixed_param_title_and_q(filename_template, q_str_for_qh):
    """
    Generál egy címet a fix paraméterekből és visszaadja a Q értéket.
    A filename_template tartalmazza a fix részeket.
    A q_str_for_qh a QH részhez tartozó string (pl. "0p0").
    """
    # Itt feltételezzük, hogy a template-ben a QH helyén egy placeholder van,
    # vagy a paramétereket a template alapján "keménykódoljuk" a title-be.
    # A te esetedben a fix paraméterek:
    # sigma = 0.05, mu = 0.00125, gamma = 0.003, p = 1/12 (0.0833...)
    # NPT=201, NIT=100
    
    # Ezeket most manuálisan adjuk meg, mivel tudjuk, hogy fixek ehhez a 4 fájlhoz.
    s_disp = parse_value_from_param_str("S0p05")
    m_disp = parse_value_from_param_str("Mp0012500000") # Ez 0.00125 lesz
    g_disp = parse_value_from_param_str("G0p003")
    p_disp = parse_value_from_param_str("Pp0833333333") # Ez 0.083333 lesz
    npt_disp = "201" # parse_value_from_param_str("NPT201") de itt a prefix más
    nit_disp = "100" # parse_value_from_param_str("NIT100")
    
    common_title = f"σ={s_disp}, μ={m_disp}, γ={g_disp}, P(jump)={p_disp}\nNPT={npt_disp}, NIT={nit_disp}"
    
    # A Q értéket a q_str_for_qh-ból nyerjük ki
    q_value_for_subplot = parse_value_from_param_str(f"QH{q_str_for_qh}") # Hozzáadjuk a QH prefixet a parse-oláshoz
    
    return common_title, q_value_for_subplot


def plot_specific_q_series_2x2(source_directory, q_str_parts, output_image_file="q_series_plot.pdf"):
    """
    A megadott q_str_parts ('0p0', '0p2', '0p4', '0p8') alapján
    összeállítja a 4 CSV fájlnevet és egy 2x2-es subplotban ábrázolja őket.
    """
    
    # A fájlnév fix részei (a QH kivételével)
    # f_out_S0p05_Mp0012500000_G0p003_QH<Q_STR>_Pp0833333333_XMINneg0p02_XMAX0p02_NPT201_NIT100.csv
    base_filename_prefix = "f_out_S0p05_Mp0012500000_G0p003_QH"
    base_filename_suffix = "_Pp0833333333_XMINneg0p02_XMAX0p02_NPT201_NIT100.csv"

    if len(q_str_parts) != 4:
        print("Hiba: Pontosan 4 Q string részre van szükség a 2x2-es plot-hoz.")
        return

    fig, axes = plt.subplots(2, 2, figsize=(16, 11), sharex=False, sharey=False) # ShareY lehet False is
    axes = axes.flatten()

    # A közös címet az első fájlnév alapján (vagy manuálisan) generáljuk
    # A q_str_parts[0] az első Q érték stringje (pl. "0p0")
    common_suptitle, _ = generate_fixed_param_title_and_q(None, q_str_parts[0]) # A template most nem kell neki
    fig.suptitle(f"Szimulációs Eredmények (f függvények)\n{common_suptitle}", fontsize=14, y=0.98)


    for i, q_str in enumerate(q_str_parts):
        ax = axes[i]
        # Teljes fájlnév összeállítása
        csv_filename = f"{base_filename_prefix}{q_str}{base_filename_suffix}"
        file_path = os.path.join(source_directory, csv_filename)

        # Q érték kinyerése a subplot címéhez
        _, q_val_for_subplot_title = generate_fixed_param_title_and_q(None, q_str)


        print(f"Plot rajzolása: {csv_filename} (Q={q_val_for_subplot_title}) -> subplot {i+1}")

        if not os.path.exists(file_path):
            ax.text(0.5, 0.5, f"Fájl nem található:\n{csv_filename}", ha='center', va='center', transform=ax.transAxes, color='red', wrap=True)
            print(f"  HIBA: A '{file_path}' fájl nem található!")
            continue
            
        try:
            df = pd.read_csv(file_path)
            if 'x_value' not in df.columns:
                ax.text(0.5, 0.5, "'x_value' oszlop hiányzik", ha='center', va='center', transform=ax.transAxes)
                continue

            x_values = df['x_value']
            f_columns_all = sorted([col for col in df.columns if re.match(r'^f\d+$', col)],
                                   key=lambda name: int(name[1:]))
            
            if not f_columns_all:
                ax.text(0.5, 0.5, "Nincsenek 'f' oszlopok", ha='center', va='center', transform=ax.transAxes)
                continue

            f_columns_to_plot = [f_columns_all[j] for j in range(0, len(f_columns_all), 10)]
            if not f_columns_to_plot:
                ax.text(0.5, 0.5, "Nincs mit ábrázolni (szűrés után)", ha='center', va='center', transform=ax.transAxes)
                continue

            for f_col in f_columns_to_plot:
                ax.plot(x_values, df[f_col], label=f_col, linewidth=1.0)
            
            combined_f_values = df[f_columns_to_plot].copy()
            valid_x_indices = combined_f_values.notna().any(axis=1)
            if valid_x_indices.any():
                first_valid_idx = valid_x_indices.idxmax()
                last_valid_idx = valid_x_indices[::-1].idxmax()
                min_x_plot = df.loc[first_valid_idx, 'x_value']
                max_x_plot = df.loc[last_valid_idx, 'x_value']
                if min_x_plot < max_x_plot:
                    ax.set_xlim([min_x_plot, max_x_plot])
            
            ax.set_title(f"Q = {q_val_for_subplot_title}", fontsize=12)
            ax.set_xlabel("x_value", fontsize=10)
            if i % 2 == 0: # Csak a bal oldali subplotok kapnak y tengely címet
                ax.set_ylabel("f értékek", fontsize=10)
            ax.tick_params(axis='both', which='major', labelsize=9)
            ax.grid(True, linestyle='--', alpha=0.6)
            
            if len(f_columns_to_plot) <= 10:
                 ax.legend(fontsize='xx-small', loc='upper right') # Próbáljuk meg a jobb felső sarkot

        except pd.errors.EmptyDataError:
            ax.text(0.5, 0.5, "Fájl üres", ha='center', va='center', transform=ax.transAxes)
        except Exception as e:
            ax.text(0.5, 0.5, f"Hiba: {e}", ha='center', va='center', transform=ax.transAxes, wrap=True)
            print(f"  Hiba történt a '{csv_filename}' feldolgozása közben: {e}")

    plt.tight_layout(rect=[0, 0.03, 1, 0.94]) # rect y1 csökkentve a suptitle jobb elhelyezéséhez
    plt.savefig(output_image_file, bbox_inches='tight')
    plt.close(fig)
    print(f"\nKész! A 2x2 plot a '{output_image_file}' fájlba mentve.")

if __name__ == "__main__":
    simulation_dir = "simulation_results/"
    # A QH paraméter string részei, amelyek a fájlnevekben szerepelnek
    # Ezek PONTOSAN azok a stringek, amik a QH után és a következő '_' előtt vannak
    q_identifier_strings = ["0p0", "0p2", "0p4", "0p8"] 
    
    output_file = "f_functions_Q_series_2x2.pdf"

    plot_specific_q_series_2x2(simulation_dir, q_identifier_strings, output_file)

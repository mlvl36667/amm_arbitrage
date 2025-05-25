# --- START OF FILE processor.py --- (A függvénydefiníciók maradnak a régiek)

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime
import os # Szükséges a fájl létezésének ellenőrzéséhez
from scipy.stats import norm as gaussian_dist # Átnevezés, hogy ne ütközzön a norm modulnévvel

# --- Konfiguráció ---
INPUT_FILE_ORIGINAL = 'november_rates.csv' # Az eredeti, nagy fájl
PROCESSED_1SEC_FILE = 'november_rates_1sec.csv' # A másodperces felbontású fájl neve

TARGET_RESOLUTION_SECONDS = 1
PLOT_JUMP_HISTOGRAM = False # Ez vezérli, hogy megjelenjenek-e a hisztogramok

# Napi alapú számításokhoz
SECONDS_IN_DAY = 24 * 60 * 60

# --- A függvénydefiníciók (process_original_file_to_seconds_resolution, 
# load_or_process_data, calculate_sigma_jump_for_threshold, 
# estimate_jump_diffusion_params_daily_from_df, 
# find_optimal_jump_threshold_grid_search_from_df) itt helyezkednek el,
# ahogy az előző kódban megadtad. Ezek nem változnak.
# ... (feltételezzük, hogy a függvényeid itt vannak) ...
def process_original_file_to_seconds_resolution(original_filepath, output_1sec_filepath, target_seconds=1):
    """
    Beolvassa az eredeti CSV fájlt soronként, másodpercenként az utolsó árat gyűjti,
    és kiírja egy új CSV fájlba a másodperces adatokat.
    Visszatér egy Pandas DataFrame-mel, ami a másodperces adatokat tartalmazza.
    """
    print(f"--- Eredeti fájl feldolgozása ({original_filepath}) másodperces felbontásra ({target_seconds}s) ---")
    secondly_prices = {} 
    line_count = 0
    processed_lines = 0
    skipped_malformed_lines = 0
    last_printed_progress_second = None

    try:
        with open(original_filepath, 'r') as f:
            for line_idx, line in enumerate(f):
                line_count += 1
                # Egyszerűsített fejléc ellenőrzés: ha az első sor tartalmazza a "Date" szót
                if line_idx == 0 and "Date" in line and "Price" in line:
                    print("   - Fejléc sor átugorva.")
                    continue
                try:
                    timestamp_str, price_str = line.strip().split(',')
                    dt_obj = datetime.strptime(timestamp_str.split('.')[0], '%Y-%m-%d %H:%M:%S')
                    price = float(price_str)
                    current_second_dt = dt_obj.replace(microsecond=0) 
                    secondly_prices[current_second_dt] = price
                    processed_lines += 1

                    if last_printed_progress_second is None or \
                       (current_second_dt - last_printed_progress_second).total_seconds() >= 60:
                        print(f"   ...feldolgozva eddig: {line_count} sor, legutóbbi idő: {current_second_dt.strftime('%Y-%m-%d %H:%M:%S')}")
                        last_printed_progress_second = current_second_dt
                except ValueError:
                    skipped_malformed_lines +=1
                    continue
                except Exception: # Általánosabb hibakezelés itt
                    skipped_malformed_lines +=1
                    continue # Folytatjuk a következő sorral
        
        print(f"--- Eredeti fájl feldolgozása befejeződött ---")
        print(f"   Összesen olvasott sor: {line_count}")
        print(f"   Sikeresen feldolgozott (árképző) sor: {processed_lines}")
        print(f"   Átugrott hibás formátumú sor: {skipped_malformed_lines}")
        print(f"   Egyedi másodpercek száma (adatpontok): {len(secondly_prices)}")

        if not secondly_prices:
            print("   Hiba: Nem sikerült adatokat kinyerni az eredeti fájlból.")
            return None

        df_secondly = pd.DataFrame(list(secondly_prices.items()), columns=['timestamp', 'price'])
        df_secondly = df_secondly.sort_values('timestamp') # Indexelés előtt rendezés
        
        # Másodperces adatok kiírása CSV-be
        print(f"   Másodperces adatok kiírása a '{output_1sec_filepath}' fájlba...")
        df_secondly.to_csv(output_1sec_filepath, index=False, date_format='%Y-%m-%d %H:%M:%S')
        print(f"   Kiírás sikeres.")
        
        df_secondly = df_secondly.set_index('timestamp')
        print(f"   Másodperces DataFrame létrehozva és indexelve a memóriában.")
        return df_secondly

    except FileNotFoundError:
        print(f"Hiba: Az eredeti fájl ({original_filepath}) nem található.")
        return None
    except Exception as e:
        print(f"Kritikus hiba az eredeti fájl feldolgozása közben: {e}")
        import traceback
        traceback.print_exc()
        return None

def load_or_process_data(original_file, processed_file):
    if os.path.exists(processed_file):
        print(f"--- Meglévő másodperces fájl betöltése: {processed_file} ---")
        try:
            df_s = pd.read_csv(processed_file)
            df_s['timestamp'] = pd.to_datetime(df_s['timestamp'])
            df_s = df_s.set_index('timestamp')
            print(f"   '{processed_file}' sikeresen betöltve és indexelve ({len(df_s)} sor).")
            if df_s.empty:
                print("   Figyelmeztetés: A betöltött másodperces fájl üres. Újragenerálás megkísérlése...")
                try:
                    os.remove(processed_file)
                    print(f"   '{processed_file}' törölve az újrageneráláshoz.")
                except OSError as e_del:
                    print(f"   Hiba a hibás '{processed_file}' törlésekor: {e_del}")
                    return None
                return process_original_file_to_seconds_resolution(original_file, processed_file, TARGET_RESOLUTION_SECONDS)
            return df_s
        except Exception as e:
            print(f"Hiba a(z) '{processed_file}' betöltése közben: {e}. Megpróbáljuk újragenerálni.")
            return process_original_file_to_seconds_resolution(original_file, processed_file, TARGET_RESOLUTION_SECONDS)
    else:
        print(f"--- Másodperces fájl ({processed_file}) nem található. Feldolgozás az eredetiből ({original_file}). ---")
        return process_original_file_to_seconds_resolution(original_file, processed_file, TARGET_RESOLUTION_SECONDS)

def calculate_sigma_jump_for_threshold(log_returns_series, dt_days, sigma_diff_daily_val, threshold_stds_val):
    if sigma_diff_daily_val <= 0 or dt_days <= 0:
        return np.nan
    jump_abs_threshold = threshold_stds_val * sigma_diff_daily_val * np.sqrt(dt_days)
    jumps = log_returns_series[np.abs(log_returns_series) > jump_abs_threshold]
    if len(jumps) >= 2:
        return np.std(jumps)
    else:
        return np.nan

def estimate_jump_diffusion_params_daily_from_df(df_secondly, jump_threshold_stds=3.0):
    if df_secondly is None or len(df_secondly) < 2:
        return None, None
    log_returns = np.log(df_secondly['price'] / df_secondly['price'].shift(1)).dropna()
    if len(log_returns) < 2:
        return None, None
    dt_days_resample = TARGET_RESOLUTION_SECONDS / SECONDS_IN_DAY
    total_time_span_seconds = (log_returns.index[-1] - log_returns.index[0]).total_seconds()
    total_time_span_days = total_time_span_seconds / SECONDS_IN_DAY
    if total_time_span_days == 0 and len(log_returns) > 0:
        total_time_span_days = dt_days_resample * len(log_returns)
    jump_abs_threshold_value = np.nan
    if len(log_returns) >= 2 and total_time_span_days > 0:
        bpv_sum = (np.pi / 2) * np.sum(np.abs(log_returns.iloc[1:]) * np.abs(log_returns.shift(1).iloc[1:]))
        sigma_diff_daily_sq = bpv_sum / total_time_span_days
        sigma_diff_daily = np.sqrt(max(0, sigma_diff_daily_sq))
    elif len(log_returns) == 1 and dt_days_resample > 0:
        sigma_diff_daily = np.abs(log_returns.iloc[0]) / np.sqrt(dt_days_resample)
    else:
        sigma_diff_daily = 0 # Default to 0 if not enough data
    
    # Biztosítjuk, hogy sigma_diff_daily ne legyen NaN, ha a fenti feltételek nem teljesülnek
    if np.isnan(sigma_diff_daily): sigma_diff_daily = 0


    if sigma_diff_daily > 0 and dt_days_resample > 0:
        jump_abs_threshold_value = jump_threshold_stds * sigma_diff_daily * np.sqrt(dt_days_resample)
        jumps = log_returns[np.abs(log_returns) > jump_abs_threshold_value]
        non_jumps = log_returns[np.abs(log_returns) <= jump_abs_threshold_value]
    else:
        jumps = pd.Series(dtype=float) # Üres series
        non_jumps = log_returns.copy() # Minden hozam nem-ugrás, ha nincs diffúziós vola
        # jump_abs_threshold_value marad NaN, ami jelzi, hogy nem volt értelmes küszöb

    num_jumps = len(jumps)
    if num_jumps > 0 and total_time_span_days > 0:
        lambda_jump_daily = num_jumps / total_time_span_days
        mu_jump = np.mean(jumps)
        sigma_jump = np.std(jumps) if num_jumps >= 2 else np.nan
    else:
        lambda_jump_daily = 0.0
        mu_jump = np.nan
        sigma_jump = np.nan

    if len(non_jumps) > 0 and dt_days_resample > 0:
        mean_log_return_non_jump = np.mean(non_jumps)
        # sigma_diff_daily itt már definiált (akár 0 is lehet)
        mu_diff_daily = (mean_log_return_non_jump / dt_days_resample) + (sigma_diff_daily**2 / 2)
    elif len(log_returns) > 0 and dt_days_resample > 0: # Fallback, ha non_jumps üres, de log_returns nem
        mean_log_return_all = np.mean(log_returns)
        mu_diff_daily = (mean_log_return_all / dt_days_resample) + (sigma_diff_daily**2 / 2)
    else:
        mu_diff_daily = np.nan
        
    params_dict = {
        'mu_gbm_daily': mu_diff_daily, 'sigma_gbm_daily': sigma_diff_daily,
        'lambda_jump_daily': lambda_jump_daily, 'mu_jump': mu_jump, 'sigma_jump': sigma_jump,
        'num_total_log_returns': len(log_returns), 'num_jumps': num_jumps,
        'jump_threshold_value_abs': jump_abs_threshold_value,
        'jump_threshold_stds_used': jump_threshold_stds,
        'total_time_span_days': total_time_span_days, 'dt_days_resample': dt_days_resample,
    }
    return params_dict, jumps


if __name__ == '__main__':
    try:
        print(f"===== Kezdődik a paraméterbecslési folyamat =====")

        df_s = load_or_process_data(INPUT_FILE_ORIGINAL, PROCESSED_1SEC_FILE)

            # ... (hibaüzenet változatlan) ...
        if df_s is not None:
            # ... (adatmennyiség ellenőrzése változatlan) ...

            TARGET_THRESHOLD_FOR_PARAMS = 3.0
            output_params_filename = "estimated_sigma_mu_q.txt" # Új fájlnév

            print(f"\n--- Paraméterek becslése a {TARGET_THRESHOLD_FOR_PARAMS:.1f} küszöbértékkel a kiíráshoz ---")
            params, identified_jumps = estimate_jump_diffusion_params_daily_from_df(
                                                          df_s,
                                                          jump_threshold_stds=TARGET_THRESHOLD_FOR_PARAMS)

            if params:
                estimated_sigma_gbm_daily = params['sigma_gbm_daily']
                estimated_mu_gbm_daily = params['mu_gbm_daily']

                # q valószínűség számítása (ugyanaz a logika, mint korábban a kiíratásnál)
                q_jump_probability_value = np.nan # Alapértelmezett, ha nem számolható
                if not np.isnan(params['lambda_jump_daily']) and not np.isnan(params['dt_days_resample']) and params['dt_days_resample'] > 0:
                    lambda_dt = params['lambda_jump_daily'] * params['dt_days_resample']
                    if lambda_dt >= 0:
                        q_jump_probability_value = 1 - np.exp(-lambda_dt)

                print(f"  Kinyert napi diffúziós volatilitás (sigma_gbm_daily): {estimated_sigma_gbm_daily:.6f}")
                print(f"  Kinyert napi diffúziós drift (mu_gbm_daily): {estimated_mu_gbm_daily:.6f}")
                print(f"  Kinyert ugrás valószínűsége dt alatt (q_jump): {q_jump_probability_value:.6f}" if not np.isnan(q_jump_probability_value) else "  Kinyert ugrás valószínűsége dt alatt (q_jump): N/A")


                try:
                    with open(output_params_filename, 'w') as f_out:
                        f_out.write(f"CONST_SIGMA={estimated_sigma_gbm_daily:.8f}\n" if not np.isnan(estimated_sigma_gbm_daily) else "CONST_SIGMA=N/A\n")
                        f_out.write(f"CONST_MU={estimated_mu_gbm_daily:.8f}\n" if not np.isnan(estimated_mu_gbm_daily) else "CONST_MU=N/A\n")
                        f_out.write(f"Q_JUMP_PROBABILITY={q_jump_probability_value:.8f}\n" if not np.isnan(q_jump_probability_value) else "Q_JUMP_PROBABILITY=N/A\n") # ÚJ SOR
                    print(f"  Becsült sigma, mu és q értékek kiírva a(z) '{output_params_filename}' fájlba.")
                except IOError as e:
                    print(f"  Hiba a(z) '{output_params_filename}' fájlba írásakor: {e}")
            else:
                print(f"  Nem sikerült paramétereket becsülni a {TARGET_THRESHOLD_FOR_PARAMS:.1f} küszöbértékkel a kiíráshoz.")

        if df_s is None:
            print("\nHiba történt az adatok betöltése/feldolgozása közben. A program leáll.")
        else:
            print(f"\nAdatok feldolgozva/betöltve ({len(df_s)} másodperces adatpont).")
            if len(df_s) < 50:
                 print("Figyelmeztetés: Nagyon kevés adatpont áll rendelkezésre a további elemzéshez.")

            thresholds_to_analyze = [2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0]

            print("\n--- Ugráseloszlások elemzése különböző küszöbértékekkel ---")
            for current_threshold in thresholds_to_analyze:
                print(f"\n--- Elemzés jump_threshold_stds = {current_threshold:.1f} értékkel ---")

                params, identified_jumps = estimate_jump_diffusion_params_daily_from_df(
                                                              df_s,
                                                              jump_threshold_stds=current_threshold)

                if params:
                    # Előre formázott stringek vagy N/A
                    sigma_gbm_daily_str = f"{params['sigma_gbm_daily']:.6f}" if not np.isnan(params['sigma_gbm_daily']) else 'N/A'
                    jump_thresh_abs_str = f"{params['jump_threshold_value_abs']:.6f}" if not np.isnan(params['jump_threshold_value_abs']) else "N/A (nem volt ugrásdetektálás)"
                    lambda_jump_daily_val = params['lambda_jump_daily'] # Érték kinyerése
                    lambda_jump_daily_str = f"{lambda_jump_daily_val:.4f}" if not np.isnan(lambda_jump_daily_val) else 'N/A'
                    mu_jump_str = f"{params['mu_jump']:.6f}" if not np.isnan(params['mu_jump']) else 'N/A'
                    sigma_jump_str = f"{params['sigma_jump']:.6f}" if not np.isnan(params['sigma_jump']) else 'N/A'

                    # q valószínűség számítása
                    q_jump_probability_str = 'N/A'
                    if not np.isnan(lambda_jump_daily_val) and not np.isnan(params['dt_days_resample']) and params['dt_days_resample'] > 0:
                        lambda_dt = lambda_jump_daily_val * params['dt_days_resample']
                        if lambda_dt >= 0: # Biztosítjuk, hogy ne legyen negatív lambda_dt (bár nem kellene)
                            q_jump_probability = 1 - np.exp(-lambda_dt)
                            q_jump_probability_str = f"{q_jump_probability:.6f}"
                        else:
                            q_jump_probability_str = "N/A (negatív lambda_dt)"


                    print(f"  Becsült diffúziós napi volatilitás (sigma_gbm_daily): {sigma_gbm_daily_str}")
                    print(f"  Ugrások azonosítási küszöbértéke (abs. log-hozam): {jump_thresh_abs_str}")
                    print(f"  Azonosított ugrások száma: {params['num_jumps']}")
                    print(f"  Napi ugrási intenzitás (lambda_jump_daily): {lambda_jump_daily_str}")
                    print(f"  Ugrás valószínűsége dt alatt (q = P(Z_t=1)): {q_jump_probability_str} (dt={params['dt_days_resample']:.2e} nap)")
                    print(f"  Ugrásméretek átlaga (mu_jump): {mu_jump_str}")
                    print(f"  Ugrásméretek szórása (sigma_jump): {sigma_jump_str}")

                    # ... (a hisztogram kódja változatlan marad) ...
                    if PLOT_JUMP_HISTOGRAM and identified_jumps is not None and not identified_jumps.empty:
                        plt.figure(figsize=(10, 7))
                        num_bins = min(50, max(10, len(identified_jumps)//5 if len(identified_jumps) > 20 else 15))
                        plt.hist(identified_jumps, bins=num_bins, alpha=0.7, color='skyblue', edgecolor='black', density=True, label=f'Empirikus ugrások ({params["num_jumps"]} db)')

                        est_mu_jump = params['mu_jump']
                        est_sigma_jump = params['sigma_jump']

                        if not np.isnan(est_mu_jump) and not np.isnan(est_sigma_jump) and est_sigma_jump > 0:
                            from scipy.stats import norm as gaussian_dist
                            xmin, xmax = plt.xlim()
                            x_fit = np.linspace(xmin, xmax, 200)
                            pdf_fit = gaussian_dist.pdf(x_fit, est_mu_jump, est_sigma_jump)
                            plt.plot(x_fit, pdf_fit, 'r--', linewidth=2, label=f'Illesztett Gauss\n(μ={est_mu_jump:.4f}, σ={est_sigma_jump:.4f})')

                            plt.axvline(est_mu_jump, color='blue', linestyle='dashed', linewidth=1.5, label=f'Átlag (μ_jump): {est_mu_jump:.4f}')
                            plt.axvline(est_mu_jump + est_sigma_jump, color='green', linestyle='dotted', linewidth=1.5, label=f'μ+σ: {(est_mu_jump + est_sigma_jump):.4f}')
                            plt.axvline(est_mu_jump - est_sigma_jump, color='green', linestyle='dotted', linewidth=1.5, label=f'μ-σ: {(est_mu_jump - est_sigma_jump):.4f}')

                        plt.title(f'Ugrásméretek Eloszlása (Jump Threshold Stds = {current_threshold:.1f})')
                        plt.xlabel('Log-hozam (Ugrásméret)')
                        plt.ylabel('Sűrűség')
                        plt.grid(True, linestyle='--', alpha=0.6)
                        plt.legend()
                        plt.tight_layout()
                        plt.show()
                    elif PLOT_JUMP_HISTOGRAM:
                        print(f"  Nincsenek azonosított ugrások a {current_threshold:.1f} küszöbértékkel a hisztogramhoz.")
                else:
                    print("  Nem sikerült paramétereket becsülni ezzel a küszöbértékkel.")

        print(f"\n===== Paraméterbecslési folyamat befejeződött =====")

    except ImportError:
        print("Hiba: Egy vagy több szükséges könyvtár (pandas, numpy, matplotlib, scipy) nincs telepítve.")
    except Exception as e:
        print(f"Fő programblokkban kritikus hiba történt: {e}")
        import traceback
        traceback.print_exc()

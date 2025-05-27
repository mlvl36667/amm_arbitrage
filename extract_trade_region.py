import os
import re
import glob

def parse_filename_and_get_ratio(filepath):
    """
    Kinyeri a paramétereket a fájlnévből és a ratio értéket a fájl tartalmából.
    Visszaad egy dictionary-t a paraméterekkel és a ratio-val, vagy None-t, ha hiba van
    vagy ha a P értéke 0.05.
    """
    filename = os.path.basename(filepath)
    
    match = re.match(
        r"tails_ratio_"
        r"Q(?P<q_raw>0p0+)_" 
        r"P(?P<p_raw>[0-9p]+)_"
        r"GP(?P<gp_raw>[0-9p]+)_"
        r"GM(?P<gm_raw>[0-9p]+)_"
        r"CS(?P<cs_raw>[0-9p]+)_"
        r"CM(?P<cm_raw>[0-9p]+)\.txt",
        filename
    )

    if not match:
        return None

    params = {}
    raw_params = match.groupdict()

    # P (Probability) feldolgozása és szűrés 0.05-re
    try:
        p_str_numeric = raw_params['p_raw'].replace('p', '.')
        p_numeric = float(p_str_numeric)
        
        # Szűrés: Ha P értéke 0.05, akkor kihagyjuk ezt a fájlt
        # Lebegőpontos számok összehasonlításánál érdemes epsilonnal dolgozni
        epsilon = 1e-9 
        if abs(p_numeric - 0.05) < epsilon:
            # print(f"  Kihagyva (P=0.05): {filename}")
            return None # Kihagyjuk ezt az eredményt

        params['p_raw_str'] = raw_params['p_raw'] # Eredeti string a P-hez

        if p_numeric != 0:
            p_inv_seconds = 1.0 / p_numeric
            if p_inv_seconds >= 60.0:
                p_inv_minutes = p_inv_seconds / 60.0
                if p_inv_minutes > 10000 or (p_inv_minutes < 0.001 and p_inv_minutes !=0):
                    params['p_inv_formatted'] = f"{p_inv_minutes:.1e} min"
                elif p_inv_minutes >= 10:
                    params['p_inv_formatted'] = f"{p_inv_minutes:.1f} min"
                else:
                    params['p_inv_formatted'] = f"{p_inv_minutes:.2f} min"
            else:
                if p_inv_seconds > 10000 or (p_inv_seconds < 0.001 and p_inv_seconds !=0):
                    params['p_inv_formatted'] = f"{p_inv_seconds:.1e} s"
                elif p_inv_seconds >=10:
                    params['p_inv_formatted'] = f"{p_inv_seconds:.1f} s"
                else:
                    params['p_inv_formatted'] = f"{p_inv_seconds:.2f} s"
        else:
            params['p_inv_formatted'] = "inf"
    except ValueError:
        return None

    # GP (Gamma Plus) feldolgozása
    try:
        gp_str_numeric = raw_params['gp_raw'].replace('p', '.')
        gp_numeric = float(gp_str_numeric)
        params['gp_numeric_for_sort'] = gp_numeric # Eredeti numerikus érték a rendezéshez
        params['gp_bps'] = gp_numeric * 10000
    except ValueError:
        return None

    # Tails ratio kiolvasása a fájlból
    try:
        with open(filepath, 'r') as f:
            content = f.read().strip()
            params['tails_ratio_percent'] = float(content)
    except (IOError, ValueError):
        params['tails_ratio_percent'] = "N/A"

    return params


def main():
    search_dir = "." 

    print(f"Fájlok keresése a '{search_dir}' könyvtárban (Q=0 feltétellel, P != 0.05)...")
    print("-" * 70)
    header = f"{'Gamma (GP) [bps]':<16} | {'P nyers':<15} | {'1/P (idő)':<12} | {'Tails Ratio (%)':<15}"
    print(header)
    print("-" * len(header))

    results = []
    processed_files_count = 0
    skipped_p_equals_005_count = 0

    for filepath in glob.glob(os.path.join(search_dir, 'tails_ratio_Q0p0*.txt')):
        data = parse_filename_and_get_ratio(filepath)
        if data:
            results.append(data)
            processed_files_count += 1
        elif os.path.exists(filepath): # Ha a parse None-t adott vissza, de a fájl létezik
            # Ellenőrizzük, hogy a P=0.05 miatt lett-e kihagyva (ez egy egyszerűsített ellenőrzés)
            filename_temp = os.path.basename(filepath)
            match_p_temp = re.search(r"P(?P<p_raw_temp>[0-9p]+)_", filename_temp)
            if match_p_temp:
                p_str_temp = match_p_temp.group('p_raw_temp').replace('p', '.')
                try:
                    if abs(float(p_str_temp) - 0.05) < 1e-9:
                        skipped_p_equals_005_count +=1
                except ValueError:
                    pass # Nem tudtuk ellenőrizni

    # Eredmények rendezése Gamma (numerikus érték) szerint növekvő sorrendbe
    # Majd másodlagosan P nyers string szerint, ha a gamma ugyanaz (opcionális)
    results.sort(key=lambda x: (x.get('gp_numeric_for_sort', float('inf')), x.get('p_raw_str', '')))


    for res in results:
        gp_bps_display = f"{res.get('gp_bps', 0.0):.0f}" 
        p_raw_display = res.get('p_raw_str', 'N/A')
        p_inv_display = res.get('p_inv_formatted', 'N/A')
        
        ratio_display = "N/A"
        if isinstance(res.get('tails_ratio_percent'), float):
            ratio_display = f"{res['tails_ratio_percent']:.5f}"

        print(f"{gp_bps_display:<16} | {p_raw_display:<15} | {p_inv_display:<12} | {ratio_display:<15}")
    
    print("-" * len(header))
    print(f"Feldolgozott és kiírt eredmények száma: {processed_files_count}")
    if skipped_p_equals_005_count > 0:
        print(f"Kihagyott fájlok (P=0.05 miatt): {skipped_p_equals_005_count}")


if __name__ == "__main__":
    main()

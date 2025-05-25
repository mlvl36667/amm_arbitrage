import os
import re

def format_duration_for_output(seconds_val):
    """
    Másodpercben megadott időtartamot formáz olvashatóbbá.
    Ha 60 mp felett van, akkor "X min Y sec" vagy "X min".
    """
    if not isinstance(seconds_val, (int, float)):
        return str(seconds_val) # Pl. "N/A (P=0)"

    if seconds_val <= 0: # Kezeletlen vagy P=0 esetekhez
        return f"{seconds_val:.2f} sec" if isinstance(seconds_val, float) else f"{seconds_val} sec"

    if seconds_val > 60:
        minutes = int(seconds_val // 60)
        remaining_seconds = seconds_val % 60
        
        if abs(remaining_seconds) < 0.001: # Ha a maradék gyakorlatilag nulla
            return f"{minutes} min"
        else:
            # Ha az eredeti másodperc egész volt, a maradék is az
            if isinstance(seconds_val, int) or remaining_seconds.is_integer():
                 return f"{minutes} min {int(remaining_seconds)} sec"
            else:
                 return f"{minutes} min {remaining_seconds:.2f} sec"
    else:
        if isinstance(seconds_val, int) or seconds_val.is_integer():
            return f"{int(seconds_val)} sec"
        else:
            return f"{seconds_val:.2f} sec"

def parse_filename_and_get_value_filtered_by_nit(directory="tails_ratios/", target_nit=16000):
    results_for_output = [] 

    filename_pattern = re.compile(
        r"tails_ratio_"
        r"Q(?P<Q_val>[0-9]+p[0-9]+)_"
        r"P(?P<P_val>[0-9]+p[0-9]+)_"
        r"GP(?P<GP_val>[0-9]+p[0-9]+)_"
        r"GM(?P<GM_val>[0-9]+p[0-9]+)_"
        r"CS(?P<CS_val>[0-9]+p[0-9]+)_"
        r"CM(?P<CM_val>[0-9]+p[0-9]+)_"
        r"NPT(?P<NPT_val>[0-9]+)_"
        r"NIT(?P<NIT_val>[0-9]+)"
        r"\.txt"
    )

    try:
        filenames = os.listdir(directory)
    except FileNotFoundError:
        print(f"# Hiba: A '{directory}' könyvtár nem található.")
        return results_for_output 
    
    filenames.sort()
    processed_any_matching_nit = False

    for filename in filenames:
        match = filename_pattern.match(filename)
        if match:
            params_raw = match.groupdict()
            params_processed = {}
            current_nit_value = -1 

            for key, value_str in params_raw.items():
                param_key_short = key.replace('_val', '')
                if param_key_short in ['Q', 'P', 'GP', 'GM', 'CS', 'CM']:
                    try:
                        params_processed[param_key_short] = float(value_str.replace('p', '.'))
                    except ValueError:
                        params_processed[param_key_short] = None
                elif param_key_short in ['NPT', 'NIT']:
                    try:
                        val_int = int(value_str)
                        params_processed[param_key_short] = val_int
                        if param_key_short == 'NIT':
                            current_nit_value = val_int
                    except ValueError:
                        params_processed[param_key_short] = None
            
            if current_nit_value == target_nit:
                processed_any_matching_nit = True
                filepath = os.path.join(directory, filename)
                try:
                    with open(filepath, 'r') as f:
                        content = f.read().strip()
                        trade_value = float(content)
                    
                    output_data = {"params": {}, "trade_region_value": trade_value}

                    # P reciprok (másodpercben tárolva)
                    p_value = params_processed.get('P')
                    if p_value is not None:
                        if p_value != 0:
                            p_reciprocal = 1 / p_value
                            output_data["params"]["P_reciprocal_seconds"] = p_reciprocal
                        else:
                            output_data["params"]["P_reciprocal_seconds"] = "N/A (P=0)"
                    else:
                         output_data["params"]["P_reciprocal_seconds"] = "N/A (P hiányzik)"
                    
                    # GP bázispontban (feltételezzük GP = GM)
                    gp_value = params_processed.get('GP')
                    if gp_value is not None:
                        output_data["params"]["GP_bps"] = gp_value * 10000 
                    
                    results_for_output.append(output_data)

                except FileNotFoundError:
                    print(f"# Figyelmeztetés: A '{filepath}' fájl nem található.")
                except ValueError:
                    print(f"# Figyelmeztetés: Nem sikerült a '{filepath}' tartalmát float értékké konvertálni: '{content}'")
                except Exception as e:
                    print(f"# Hiba a '{filepath}' fájl feldolgozása közben: {e}")
    
    if not processed_any_matching_nit and os.path.exists(directory):
        print(f"# Nem található 'tails_ratio_*.txt' fájl a '{directory}' könyvtárban, amelynek NIT értéke {target_nit}.")

    return results_for_output

if __name__ == "__main__":
    target_directory = "tails_ratios/"
    nit_filter_value = 16000
    
    extracted_data = parse_filename_and_get_value_filtered_by_nit(target_directory, nit_filter_value)

    if extracted_data:
        for item_index, item in enumerate(extracted_data):
            if item_index > 0: 
                print("-" * 30) 

            params = item['params']
            
            # P reciprokának (várakozási idő) kiírása formázva
            p_reciprocal_sec_val = params.get("P_reciprocal_seconds", "N/A")
            formatted_duration = format_duration_for_output(p_reciprocal_sec_val)
            print(f"  Várakozási idő: {formatted_duration}")
            
            # GP kiírása bázispontban
            if 'GP_bps' in params:
                gp_bps_val = params['GP_bps']
                # Ha a bázispont egész szám, ne írjunk ki felesleges tizedeseket
                if isinstance(gp_bps_val, float) and gp_bps_val.is_integer():
                    print(f"  GP (bps): {int(gp_bps_val)}")
                elif isinstance(gp_bps_val, float):
                     # Általában a bázispontot néhány tizedesjegyig érdemes kiírni, ha nem egész
                    print(f"  GP (bps): {gp_bps_val:.2f}") # Módosítsd a .2f-et, ha más pontosság kell
                else: # int
                    print(f"  GP (bps): {gp_bps_val}")


            # Trade Region Érték kiírása
            print(f"  Trade Region Érték: {item['trade_region_value']:.8f}")

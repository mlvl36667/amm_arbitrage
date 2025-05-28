import re
import pandas as pd
import numpy as np
import os
import glob # Import the glob module for listing files

# List of filenames to process
# Now dynamically collecting all CSV files from the current directory.
# Only considering files where NIT = 1000 and NPT = 801.
all_csv_files = glob.glob('simulation_results/*.csv')
filtered_filenames = []

# Regex to extract NIT and NPT values
nit_regex = re.compile(r"_NIT(\d+)\.csv$")
npt_regex = re.compile(r"_NPT(\d+)_")

for filename in all_csv_files:
    nit_match = nit_regex.search(filename)
    npt_match = npt_regex.search(filename)
    
    if nit_match and npt_match:
        nit_value = int(nit_match.group(1))
        npt_value = int(npt_match.group(1))
        
        if nit_value == 1000 and npt_value == 801:
            filtered_filenames.append(filename)

filenames = filtered_filenames

# Sorting filenames for consistent processing
filenames.sort()

def parse_filename(filename):
    """
    Parses the filename to extract S, M, G, Q, P, XMIN, XMAX, NPT, NIT parameters.
    Handles 'p' as decimal separator and 'neg' for negative numbers.
    """
    
    match = re.search(
        r"S(?P<S>[\d\w.]+)_M(?P<M>[\d\w.]+)_G(?P<G>[\d\w.]+)_QH(?P<QH>[\d\w.]+)_P(?P<P>[\d\w.]+)_XMIN(?P<XMIN>[neg\d\w.]+)_XMAX(?P<XMAX>[\d\w.]+)_NPT(?P<NPT>\d+)_NIT(?P<NIT>\d+)\.csv",
        filename
    )
    if not match:
        raise ValueError(f"The filename does not match the expected pattern: {filename}")

    params = match.groupdict()

    parsed_params = {}
    for key, value in params.items():
        if 'p' in value:
            value = value.replace('p', '.')
        if 'neg' in value:
            value = value.replace('neg', '-')
        
        if key in ['NPT', 'NIT']:
            parsed_params[key] = int(value)
        else:
            parsed_params[key] = float(value)
            
    parsed_params['p_reciprocal'] = 1 / parsed_params['P'] # Store as p_reciprocal for clarity
    parsed_params['Gamma'] = parsed_params['G']

    return parsed_params

def read_f_values(filepath):
    """
    Reads the last column of a CSV file as f(t) values, ensuring they are floating-point numbers.
    Attempts to infer the header and handles non-numeric values by coercing them to NaN.
    """
    df = pd.read_csv(filepath, header='infer')
    
    last_col_name = df.columns[-1]
    
    # Use pd.to_numeric with errors='coerce' to convert non-numeric values to NaN
    f_values = pd.to_numeric(df[last_col_name], errors='coerce').values
    
    # Drop NaN values if you don't want them in calculations, or handle them specifically.
    # For trapz, NaN values will propagate, often resulting in NaN for the integral.
    # If a row is entirely NaN, it might be better to drop it.
    # For now, we'll let NaN propagate and trust trapz handles it gracefully.
    return f_values

def calculate_t_values(xmin, xmax, npt):
    """Generates evenly spaced t values."""
    return np.linspace(xmin, xmax, npt)

def calculate_arb(p_val, gamma_val, t_values, f_values):
    """Calculates the ARB integral."""
    theta = 0.5
    c1 = (theta / (1 - theta))**(1 - theta)
    c2 = ((1 - theta) / theta)**theta

    # First integral part: t >= gamma
    mask1 = t_values >= gamma_val
    t1 = t_values[mask1]
    f1 = f_values[mask1]
    
    integrand1 = c1 * (np.exp(t1 * (1 - theta)) - np.exp(gamma_val * (1 - theta))) + \
                 c2 * (np.exp(gamma_val - t1 * theta) - np.exp(gamma_val * (1 - theta)))
    
    min_len1 = min(len(integrand1), len(f1))
    integral1 = np.trapz(integrand1[:min_len1] * f1[:min_len1], t1[:min_len1])

    # Second integral part: t <= -gamma
    mask2 = t_values <= -gamma_val
    t2 = t_values[mask2]
    f2 = f_values[mask2]
    
    integrand2 = c1 * (np.exp(t2 * (1 - theta)) - np.exp(gamma_val * (theta - 1))) + \
                 c2 * (np.exp(-t2 * theta - gamma_val) - np.exp(gamma_val * (theta - 1)))
    
    min_len2 = min(len(integrand2), len(f2))
    integral2 = np.trapz(integrand2[:min_len2] * f2[:min_len2], t2[:min_len2])

    arb = p_val * (integral1 + integral2)
    return arb

# Process files and collect results
results_data = []
print("\n--- Calculation Results ---")
if not filenames:
    print("No CSV file in the current directory matches the 'NIT = 1000' and 'NPT = 801' conditions.")
else:
    for filename in filenames:
        try:
            params = parse_filename(filename)
            
            f_values = read_f_values(filename)
            t_values = calculate_t_values(params['XMIN'], params['XMAX'], params['NPT'])
            
            if len(f_values) != len(t_values):
                min_len = min(len(f_values), len(t_values))
                f_values = f_values[:min_len]
                t_values = t_values[:min_len]
                print(f"Warning: Mismatch in length between f_values ({len(f_values)}) and t_values ({len(t_values)}) for {filename}. Truncating to {min_len} points.")

            arb_value = calculate_arb(params['p_reciprocal'], params['Gamma'], t_values, f_values) # Use p_reciprocal here
            
            results_data.append({
                'p_reciprocal': round(params['p_reciprocal'], 1), # Round p_reciprocal to 1 decimal place
                'Gamma': params['Gamma'],
                'Calculated ARB': arb_value
            })

        except ValueError as ve:
            print(f"Error processing {filename} (filename parsing or data conversion error): {ve}\n")
        except Exception as e:
            print(f"Error processing {filename}: {e}\n")

# Create a DataFrame from the collected results
if results_data:
    df_results = pd.DataFrame(results_data)
    
    # Define mappings for display names
    gamma_display_map = {
        0.0001: '1 bp',
        0.0005: '5 bp',
        0.001: '10 bp',
        0.003: '30 bp',
        0.01: '100 bp'
    }
    
    p_reciprocal_display_map = {
        600.0: '10 min',
        120.0: '2 min',
        12.0: '12 sec',
        2.0: '2 sec'
    }

    # Define the desired order for p_reciprocal for table display
    # This list explicitly defines the row order for the final table.
    desired_p_reciprocal_order = [600.0, 120.0, 12.0, 2.0]

    # Pivot the table
    pivot_table = df_results.pivot_table(
        index='p_reciprocal',
        columns='Gamma',
        values='Calculated ARB'
    )

    # Reindex the pivot table to ensure the desired row order
    # Any p_reciprocal values not in desired_p_reciprocal_order will be dropped or appear as NaN if reindex has fill_value=None.
    pivot_table = pivot_table.reindex(desired_p_reciprocal_order)

    # Ensure Gamma columns are in the correct order for display
    gamma_sorted_keys = sorted(gamma_display_map.keys())
    pivot_table = pivot_table.reindex(columns=gamma_sorted_keys)

    # Rename columns and index using the display maps
    pivot_table.rename(columns=gamma_display_map, inplace=True)
    pivot_table.rename(index=p_reciprocal_display_map, inplace=True)

    # Format the cell values using .map for deprecation warning fix
    # Use f-string formatting for scientific notation if values are very small
    # Otherwise, format to 7 decimal places
    formatted_table = pivot_table.map(lambda x: f"{x:.1e}" if abs(x) < 1e-4 and x != 0 else f"{x:.7f}")

    print("\n--- Sorted ARB Results Table ---")
    # Print the LaTeX-like table header
    # Using a raw string (r"...") to avoid SyntaxWarning for '\D'
    header_gamma_names = [gamma_display_map[g] for g in gamma_sorted_keys]
    header = r"$\Delta t \setminus \gamma$" + " & " + " & ".join(header_gamma_names) + " \\\\"
    print(header)
    print("\\hline")

    # Print the table rows
    for index, row in formatted_table.iterrows():
        # Handle NaN values in the formatted table if any data was coerced
        row_values_str = [val if pd.notna(val) else 'NaN' for val in row.astype(str)]
        row_str = f"{index} & " + " & ".join(row_values_str) + " \\\\"
        print(row_str)
else:
    print("No results to display in a table.")


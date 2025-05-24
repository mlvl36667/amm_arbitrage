import os
import re
import glob

def parse_filename_and_get_ratio(filepath):
    """
    Extracts parameters from the filename and the ratio value from the file content.
    Returns a dictionary with the parameters and the ratio, or None if there is an error
    or if the P value is 0.05.
    """
    filename = os.path.basename(filepath)

    match = re.match(
        r"tails_ratio_"
        r"Q(?P<q_raw>0p0+)_"  # Matches Q0p0, Q0p00 etc. (Q=0.0)
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

    # Processing P (Probability) and filtering for 0.05
    try:
        p_str_numeric = raw_params['p_raw'].replace('p', '.')
        p_numeric = float(p_str_numeric)

        # Filter: If P value is 0.05, then skip this file
        # When comparing floating-point numbers, it's good to use an epsilon
        epsilon = 1e-9
        if abs(p_numeric - 0.05) < epsilon:
            # print(f"  Skipped (P=0.05): {filename}")
            return None  # Skip this result

        params['p_raw_str'] = raw_params['p_raw']  # Original string for P

        if p_numeric != 0:
            p_inv_seconds = 1.0 / p_numeric
            if p_inv_seconds >= 60.0:
                p_inv_minutes = p_inv_seconds / 60.0
                if p_inv_minutes > 10000 or (p_inv_minutes < 0.001 and p_inv_minutes != 0):
                    params['p_inv_formatted'] = f"{p_inv_minutes:.1e} min"
                elif p_inv_minutes >= 10:
                    params['p_inv_formatted'] = f"{p_inv_minutes:.1f} min"
                else:
                    params['p_inv_formatted'] = f"{p_inv_minutes:.2f} min"
            else:
                if p_inv_seconds > 10000 or (p_inv_seconds < 0.001 and p_inv_seconds != 0):
                    params['p_inv_formatted'] = f"{p_inv_seconds:.1e} s"
                elif p_inv_seconds >= 10:
                    params['p_inv_formatted'] = f"{p_inv_seconds:.1f} s"
                else:
                    params['p_inv_formatted'] = f"{p_inv_seconds:.2f} s"
        else:
            params['p_inv_formatted'] = "inf"
    except ValueError:
        return None

    # Processing GP (Gamma Plus)
    try:
        gp_str_numeric = raw_params['gp_raw'].replace('p', '.')
        gp_numeric = float(gp_str_numeric)
        params['gp_numeric_for_sort'] = gp_numeric  # Original numeric value for sorting
        params['gp_bps'] = gp_numeric * 10000
    except ValueError:
        return None

    # Reading tails ratio from the file
    try:
        with open(filepath, 'r') as f:
            content = f.read().strip()
            params['tails_ratio_percent'] = float(content)
    except (IOError, ValueError):
        params['tails_ratio_percent'] = "N/A" # Indicate if reading failed

    return params


def main():
    search_dir = "./tails_ratios"

    print(f"Searching for files in the '{search_dir}' directory (with Q=0.0 condition, P != 0.05)...")
    print("-" * 70)
    header = f"{'Gamma (GP) [bps]':<16} | {'P raw':<15} | {'1/P (time)':<12} | {'Tails Ratio (%)':<15}"
    print(header)
    print("-" * len(header))

    results = []
    processed_files_count = 0
    skipped_p_equals_005_count = 0

    # Find files matching the pattern for Q=0.0
    for filepath in glob.glob(os.path.join(search_dir, 'tails_ratio_Q0p0*.txt')):
        data = parse_filename_and_get_ratio(filepath)
        if data:
            results.append(data)
            processed_files_count += 1
        elif os.path.exists(filepath):  # If parse returned None, but the file exists
            # Check if it was skipped due to P=0.05 (this is a simplified check for counting)
            filename_temp = os.path.basename(filepath)
            match_p_temp = re.search(r"P(?P<p_raw_temp>[0-9p]+)_", filename_temp)
            if match_p_temp:
                p_str_temp = match_p_temp.group('p_raw_temp').replace('p', '.')
                try:
                    if abs(float(p_str_temp) - 0.05) < 1e-9:
                        skipped_p_equals_005_count += 1
                except ValueError:
                    pass  # Could not verify the P value for skipping count

    # Sort results by Gamma (numeric value) in ascending order
    # Then secondarily by P raw string if gamma is the same (optional, but good for consistency)
    results.sort(key=lambda x: (x.get('gp_numeric_for_sort', float('inf')), x.get('p_raw_str', '')))

    for res in results:
        gp_bps_display = f"{res.get('gp_bps', 0.0):.0f}"
        p_raw_display = res.get('p_raw_str', 'N/A')
        p_inv_display = res.get('p_inv_formatted', 'N/A')

        ratio_value = res.get('tails_ratio_percent')
        ratio_display = "N/A"
        if isinstance(ratio_value, float):
            ratio_display = f"{ratio_value:.5f}"

        print(f"{gp_bps_display:<16} | {p_raw_display:<15} | {p_inv_display:<12} | {ratio_display:<15}")

    print("-" * len(header))
    print(f"Number of processed and printed results: {processed_files_count}")
    if skipped_p_equals_005_count > 0:
        print(f"Skipped files (due to P=0.05): {skipped_p_equals_005_count}")


if __name__ == "__main__":
    main()

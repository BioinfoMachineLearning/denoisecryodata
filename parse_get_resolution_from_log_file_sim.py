"""

author: nabin 
timestamp: Tue March 13 2024 10:09 AM


- rparse output file of phenix.mtriage and extracts relevant informations


"""


import re
import csv
import os
# Define the regular expressions to match the required values
regex_patterns = {
    'd99': r'using map alone \(d99\) +: +(\d+\.\d+) +(\d+\.\d+)',
    'fsc_0.143': r'FSC\(map,model map\)=0\.143 +: +(\d+\.\d+) +(\d+\.\d+)',
    'fsc_0.5': r'FSC\(map,model map\)=0\.5 +: +(\d+\.\d+) +(\d+\.\d+)',
}

# Function to extract values from the log file
def extract_values(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        values = {}
        for key, pattern in regex_patterns.items():
            match = re.search(pattern, content)
            if match:
                masked_value = match.group(1)
                unmasked_value = match.group(2)
                values[key] = {'masked': masked_value, 'unmasked': unmasked_value}
        return values

# Function to save extracted values to CSV
def save_to_csv(data, output_file, target):
    with open(output_file, 'a', newline='') as csvfile:
        fieldnames = ['Target','Simulated Map Alone (Masked)', 'Simulated Map Alone (Unmasked)', 'Simulated FSC_0143 (Masked)', 'Simulated FSC_0143 (Unmasked)','Simulated FSC_05 (Masked)', 'Simulated FSC_05 (Unmasked)']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        # writer.writeheader()
        writer.writerow({'Target' : target,
                        'Simulated Map Alone (Masked)': data['d99'].get('masked', ''),
                        'Simulated Map Alone (Unmasked)': data['d99'].get('unmasked', ''),
                        'Simulated FSC_0143 (Masked)': data.get('fsc_0.143', {}).get('masked', ''),
                        'Simulated FSC_0143 (Unmasked)': data.get('fsc_0.143', {}).get('unmasked', ''),
                        'Simulated FSC_05 (Masked)': data.get('fsc_0.5', {}).get('masked', ''),
                        'Simulated FSC_05 (Unmasked)': data.get('fsc_0.5', {}).get('unmasked', '') })

        # for metric, values in data.items():

            # writer.writerow({'Metric': metric, 'Masked': values['masked'], 'Unmasked': values['unmasked']})

# Main function
def main():
    dataset_dir = "/situs_simulated_pdb"
    output_csv  = "/curr_situs_simulated_fsc_values.csv"
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['Target','Simulated Map Alone (Masked)', 'Simulated Map Alone (Unmasked)', 'Simulated FSC_0143 (Masked)', 'Simulated FSC_0143 (Unmasked)','Simulated FSC_05 (Masked)', 'Simulated FSC_05 (Unmasked)']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
    map_names = os.listdir(dataset_dir)
    map_directories = [item for item in map_names if os.path.isdir(os.path.join(dataset_dir, item))]
    for m in map_directories:
        density_map_path = f"{dataset_dir}/{m}"
        log_file = f"{density_map_path}/sim_mtriage.out"
        if os.path.exists(log_file):
            extracted_values = extract_values(log_file)
            # print(extracted_values)
            save_to_csv(extracted_values, output_csv, target=m)
            # print("Values extracted and saved to", output_csv)
        else:
            print("NO Output File genrated for: ", m)

if __name__ == "__main__":
    main()
    print("DONE!")

"""

author: nabin 
timestamp: Tue March 12 2024 08:17 PM


- runs phenix.mtriage on maps


"""

data_directory = "/situs_simulated_pdb"

import pandas as pd
import os
import subprocess
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: pass in atleast one density map")
        sys.exit(1)
    else:
        values = sys.argv[1]
    
        map_data_directory = f"{data_directory}/{values}"
        density_map = f"{map_data_directory}/{values}_situs_simulated.mrc"
        pdb_file = f"{map_data_directory}/{values}.pdb1"
        logs = f"{map_data_directory}/sim_mtriage.out"
        if os.path.exists(logs):
            os.remove(logs)
        os.chdir(map_data_directory)
        phenix_mtriage_cmd = ['phenix.mtriage', pdb_file, density_map, 'nproc=12']
        print(phenix_mtriage_cmd)    

        try:
            result = subprocess.run(phenix_mtriage_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
            stdout = result.stdout
            stderr = result.stderr
            return_code = result.returncode
            if return_code == 0:
                print("Command executed successfully.", density_map)
                with open(logs, "w") as output_file:
                    output_file.write(stdout)

            else:
                print(f"Command failed with exit code {return_code}.")
                print("Standard Error:")
                print(stderr)
        except subprocess.CalledProcessError as e:
            print(f"Command failed with exit code {e.returncode}.")
            print("Standard Error:")
            print(e.stderr)
        except Exception as e:
            print(f"An error occurred: {str(e)}")




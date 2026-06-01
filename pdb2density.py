import pexpect
import re
import os
import sys

input_file = "density_maps.txt"

def run_simulation(pdb_path, resolution, output_path, voxel_spacing=1):
    pdb2vol_path = "path/to/pdb2vol"

    if not os.path.exists(pdb2vol_path):
        print(f"Error: pdb2vol execution file not found at {pdb2vol_path}")
        return

    if not os.path.exists(pdb_path):
        print(f"PDB file not found: {pdb_path}")
        return

    if os.path.exists(output_path):
        print(f"Simulated map already exists: {output_path}, skipping.")
        return

    data_number = os.path.splitext(os.path.basename(pdb_path))[0]
    print(f"Processing {data_number} with resolution {resolution}...")

    # Regular expressions to match expected prompts
    prompts = [
        re.compile(r"exclude.*?water|ignore.*?water", re.I),           # 0
        re.compile(r"exclude.*?hydrogen|ignore.*?hydrogen", re.I),     # 1
        re.compile(r"exclude.*?codebook|ignore.*?codebook", re.I),     # 2
        re.compile(r"mass-weight", re.I),            # 3
        re.compile(r"B-factor", re.I),               # 4
        re.compile(r"spacing", re.I),                # 5
        re.compile(r"width|resolution", re.I),       # 6
        re.compile(r"Kernel type", re.I),            # 7
        re.compile(r"correct for lattice", re.I),    # 8
        re.compile(r"amplitude", re.I),              # 9
        pexpect.EOF,                                 # 10
        pexpect.TIMEOUT                              # 11
    ]

    # The command is: ./Situs_3.2/bin/pdb2vol input.pdb1 output.mrc
    cmd = f"{pdb2vol_path} {pdb_path} {output_path}"

    try:
        # We start the command and give it a timeout of 300 seconds (5 minutes) per step
        # to allow for heavy processing times if PDB involves many atoms
        child = pexpect.spawn(cmd, encoding='utf-8', timeout=300)

        # To see what is actually happening behind the scenes (useful for debugging problems)
        child.logfile = sys.stdout

        while True:
            idx = child.expect(prompts)

            if idx in [0, 1, 2]:
                # Ignore water, hydrogen, codebook -> Yes (2)
                child.sendline("2")
            elif idx in [3, 4]:
                # Mass-weight, B-factor threshold -> No (1)
                child.sendline("1")
            elif idx == 5:
                # Voxel spacing
                child.sendline(str(voxel_spacing))
            elif idx == 6:
                # Kernel width / resolution
                child.sendline(f"-{resolution}")
            elif idx == 7:
                # Kernel type: Gaussian is usually '1'
                child.sendline("1")
            elif idx == 8:
                # Correct for lattice interpolation smoothing
                child.sendline("1")
            elif idx == 9:
                # Maximum amplitude
                child.sendline("1.0")
            elif idx == 10:
                print(f"\n--- Success: Finished processing {data_number} ---\n")
                break
            elif idx == 11:
                print(f"\n--- Error: Timeout reached processing {data_number} ---\n")
                break

    except Exception as e:
        print(f"Failed to process {data_number}: {e}")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb_file", help="Path to input PDB file")
    parser.add_argument("resolution", help="Resolution value")
    parser.add_argument("output_file", help="Path to output MRC file")
    parser.add_argument("--voxel-spacing", type=float, default=1, help="Voxel spacing (default: 1)")
    args = parser.parse_args()
    run_simulation(args.pdb_file, args.resolution, args.output_file, args.voxel_spacing)
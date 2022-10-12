"""Script to store sting data in a more convenient .mat format."""

from scipy.io import savemat

ANGLES = ["-5", "00", "05", "10", "12", "14", "15", "18", "20"]

PREFIXES = {
    "force": "force_mes"
}

FOLDER = "Force Measurements"
src_fns = [f"./{FOLDER}/{PREFIXES['force']}_a_{angle}.dat" for angle in ANGLES]
dest_fns = [f"./{FOLDER}/{PREFIXES['force']}_a_{angle}.mat" for angle in ANGLES]

KEYS = {
    "Static Pressure": "P",
    "Density": "rho",
    "Fixed Pitot Probe Speed": "v",
    "Normal Force": "Fn",
    "Axial Force": "Fa",
    "Angle of Attack": "a",
}

for i, fn in enumerate(src_fns):
    with open(fn, "r") as f:
        lines = f.readlines()
        data = {}
        for j, line in enumerate(lines):
            key = line.split("=")[0].strip()
            if key in KEYS.keys():
                try:
                    val = line.split("=")[1].strip()
                    val = val.split("\t")[0] # Remove tabs if present
                except Exception:
                    # print(f"Error: Couldn't parse val of line {j}, setting as 'None'")
                    val = None
                data[KEYS[key]] = float(val) # Create dict with symbols
    print(data)
    savemat(dest_fns[i], data)
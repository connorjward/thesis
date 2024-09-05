import pathlib

import matplotlib.pyplot as plt
import pandas as pd
import tomli

from experiments.common.roofline import compute_roofline_data


machine_name = "pyramus"
dir = pathlib.Path(__file__).parent
batchfile = f"{dir}/batch-elitebook-20240905-135923.csv"

with open(f"{dir.parent}/machines/{machine_name}.toml", "rb") as f:
    machine_info = tomli.load(f)

peak_bandwidth = machine_info["peak_bandwidth"]
peak_flops_scalar = machine_info["peak_flops_scalar"]
peak_flops_avx = machine_info["peak_flops_avx"]

data = pd.read_csv(batchfile)
compute_roofline_data(data, peak_bandwidth)

print(data)

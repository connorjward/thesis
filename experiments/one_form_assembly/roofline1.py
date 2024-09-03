import os

import matplotlib.pyplot as plt
import pandas as pd
import tomli

from experiments.common.roofline import parse_csv

dir = os.path.dirname(__file__)
machine_file = f"{dir}/../elitebook.toml"
# batch_file_pyop2 = "???"
batch_file_pyop3 = f"{dir}/batch-elitebook-pyop3-20240903-135749.csv"

with open(machine_file, "rb") as f:
    machine = tomli.load(f)
bandwidth = machine["main_memory"]

# roofline_info_pyop2 = parse_csv(batch_file_pyop2, bandwidth)
roofline_info_pyop3 = parse_csv(batch_file_pyop3, bandwidth)

print(roofline_info_pyop3)

# now plot things
xmin, xmax = 1e-5, 5e-3
ymin, ymax = 5e4,1e8

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xscale("log")
plt.yscale("log")

# peak bandwidth
plt.plot([xmin, xmax], [xmin*bandwidth*1e6, xmax*bandwidth*1e6])

# pyop3, optimal
plt.scatter(
    roofline_info_pyop3["optimal_ai"],
    roofline_info_pyop3["flop_rate"],
)

# pyop3, pessimal
plt.scatter(
    roofline_info_pyop3["pessimal_ai"],
    roofline_info_pyop3["flop_rate"],
)


# TODO: make a pgf
plt.savefig(f"{dir}/roofline1.png")

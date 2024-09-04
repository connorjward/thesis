import os

import matplotlib.pyplot as plt
import pandas as pd
import tomli

from experiments.common.roofline import parse_csv

dir = os.path.dirname(__file__)
machine_file = f"{dir}/../elitebook.toml"
batch_file_pyop2 = f"{dir}/batch-elitebook-pyop2-20240903-173912.csv"
batch_file_pyop3 = f"{dir}/batch-elitebook-pyop3-20240903-155820.csv"

with open(machine_file, "rb") as f:
    machine = tomli.load(f)
bandwidth = machine["main_memory"]

roofline_info_pyop2 = parse_csv(batch_file_pyop2, bandwidth)
roofline_info_pyop3 = parse_csv(batch_file_pyop3, bandwidth)

print(roofline_info_pyop3)

# now plot things
xmin, xmax = 1e-1, 1e2
ymin, ymax = 1e8,1e11

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)
plt.xscale("log")
plt.yscale("log")

# peak bandwidth
plt.plot([xmin, xmax], [xmin*bandwidth*1e6, xmax*bandwidth*1e6])

# peak throughput (from linpack)
plt.plot([xmin, xmax], [45e9, 45e9])

# pyop3, optimal
# plt.scatter(
#     roofline_info_pyop3["optimal_ai"],
#     roofline_info_pyop3["flop_rate"],
# )

# pyop3, pessimal
# plt.scatter(
#     roofline_info_pyop3["pessimal_ai"],
#     roofline_info_pyop3["flop_rate"],
# )
mins = roofline_info_pyop2["pessimal_ai"]
maxs = roofline_info_pyop2["optimal_ai"]
means = (mins+maxs) / 2

plt.errorbar(
    means,
    roofline_info_pyop2["flop_rate"],
    xerr=(means-mins, maxs-means),
    fmt=".",
    capsize=4,
    color="blue",
)

mins = roofline_info_pyop3["pessimal_ai"]
maxs = roofline_info_pyop3["optimal_ai"]
means = (mins+maxs) / 2

plt.errorbar(
    means,
    roofline_info_pyop3["flop_rate"],
    xerr=(means-mins, maxs-means),
    fmt=".",
    capsize=4,
    color="red",
)


# TODO: make a pgf
plt.savefig(f"{dir}/roofline2.png")

import pathlib

import matplotlib.pyplot as plt
import pandas as pd

from experiments.common.roofline import *


BANDWIDTH_COLOR = "red"
THROUGHPUT_COLOR = "blue"


dir = pathlib.Path(__file__).parent

plt.style.use(MPL_PARAMS)
fig = plt.figure(**FIGURE_PARAMS)

ax1 = fig.add_subplot()
ax1.set_xlabel("Core count")
ax1.set_xticks(range(2, 33, 2))

streamdata = pd.read_csv(f"{dir}/streamdata.csv")
ax1.set_ylim(0, 5e10)
ax1.set_ylabel("Memory bandwidth (byte/s)", color=BANDWIDTH_COLOR)
ax1.tick_params(axis="y", labelcolor=BANDWIDTH_COLOR)
ax1.plot(streamdata["nthreads"], streamdata["bandwidth"], marker="o", color=BANDWIDTH_COLOR, linestyle="dashed",markersize=4, linewidth=1)

ax2 = ax1.twinx()
flopdata = pd.read_csv(f"{dir}/flopdata.csv")
ax2.set_ylabel("Arithmetic throughput (FLOP/s)", color=THROUGHPUT_COLOR)
ax2.tick_params(axis="y", labelcolor=THROUGHPUT_COLOR)
ax2.plot(flopdata["nthreads"], flopdata["throughput"], marker="D", color=THROUGHPUT_COLOR, linestyle="dashed",markersize=4, linewidth=1)

plt.savefig(f"{dir}/bandwidth.pdf", dpi=200)

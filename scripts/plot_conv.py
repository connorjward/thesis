import pathlib

import matplotlib.pyplot as plt
import pandas as pd

from experiments.common.roofline import MPL_PARAMS, FIGURE_PARAMS


def h(nref):
    return 1/2**nref


dir = pathlib.Path(__file__).parent

FILENAME = "helmholtz_conv.csv"
raw = pd.read_csv(f"{dir}/{FILENAME}")

plt.style.use(MPL_PARAMS)
fig = plt.figure(**FIGURE_PARAMS)

ax1 = fig.add_subplot()
ax1.set_xlabel("Refinement")
ax1.set_xticks([4, 5, 6, 7])
ax1.set_ylabel("$L_2$ error")
ax1.set_yscale("log")


LINE_STYLE = {"linewidth": 1, "marker": "x"}
CONV_LINE_STYLE = {"linewidth": 1, "linestyle": "dashed"}

# degree 1
filtered = raw.loc[raw["degree"] == 1]
ax1.plot(filtered["refinement"], filtered["norm"], **LINE_STYLE)
ax1.plot([4, 7], [h(4)**1, h(7)**2], **CONV_LINE_STYLE)

# degree 2
# ax1.plot(raw["refinement"], raw["norm"], linewidth=1)
#
# # degree 3
# ax1.plot(raw["refinement"], raw["norm"], linewidth=1)

plt.savefig(f"{dir}/conv.pdf", dpi=200)

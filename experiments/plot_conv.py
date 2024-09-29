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


LINE_STYLE = {"linewidth": 1.5, "marker": "x"}
CONV_LINE_STYLE = {"linewidth": 1, "linestyle": "dashed", "color": "black"}

# convergence lines
ax1.plot([4, 7], [h(4)**1, h(7)**1], **CONV_LINE_STYLE)
ax1.plot([4, 7], [h(4)**2, h(7)**2], **CONV_LINE_STYLE)
ax1.plot([4, 7], [h(4)**3, h(7)**3], **CONV_LINE_STYLE)
ax1.plot([4, 7], [h(4)**4, h(7)**4], **CONV_LINE_STYLE)

for conv in range(1, 5):
    ax1.annotate(
        f"$h^{conv}$", (6.8, h(6.8)**conv),
        (0, 2),
        textcoords="offset points",
        verticalalignment="bottom",
    )

# degree 1
filtered = raw.loc[raw["degree"] == 1]
ax1.plot(filtered["refinement"], filtered["norm"], **LINE_STYLE, label="$p=1$")

# degree 2
filtered = raw.loc[raw["degree"] == 2]
ax1.plot(filtered["refinement"], filtered["norm"], **LINE_STYLE, label="$p=2$")

# degree 3
filtered = raw.loc[raw["degree"] == 3]
ax1.plot(filtered["refinement"], filtered["norm"], **LINE_STYLE, label="$p=3$")

ax1.legend()
plt.savefig(f"{dir}/conv.pdf", dpi=200)

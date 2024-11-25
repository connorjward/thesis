import pathlib

import matplotlib.pyplot as plt

from experiments.common.roofline import plot_roofline


dir = pathlib.Path(__file__).parent
plot_roofline(
    "pyramus",
    f"{dir}/batch-pyramus-pyop2-20241125-140514.csv",
    f"{dir}/batch-pyramus-pyop3-20241125-152110.csv",
    experiment=1,
)
plt.savefig(f"{dir}/roofline1_new.pdf", dpi=200, backend="pgf")

import pathlib

import matplotlib.pyplot as plt

from experiments.common.roofline import plot_roofline


dir = pathlib.Path(__file__).parent
plot_roofline(
    "pyramus",
    f"{dir}/batch-pyramus-pyop2-20241206-114346.csv",
    f"{dir}/batch-pyramus-pyop3-20241206-112037.csv",
    experiment=2,
)
plt.savefig(f"{dir}/roofline2_new.pdf", dpi=200, backend="pgf")

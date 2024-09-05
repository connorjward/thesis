import pathlib

import matplotlib.pyplot as plt

from experiments.common.roofline import plot_roofline


dir = pathlib.Path(__file__).parent
plot_roofline(
    "pyramus",
    f"{dir}/batch-pyramus-pyop2-20240904-172928.csv",
    f"{dir}/batch-pyramus-pyop3-20240904-160249.csv",
    experiment=2,
)
plt.savefig(f"{dir}/roofline2.pdf", dpi=200)

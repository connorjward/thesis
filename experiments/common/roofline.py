import pathlib

import matplotlib.pyplot as plt
import pandas as pd
import tomli


MPL_PARAMS = {
    "font.family": "serif",
    "font.size": 10,
    "text.usetex": True,
}

FIGURE_PARAMS = {
    "figsize": (6, 3.5),
    "layout": "tight",
}

ROOFLINE_LINE_STYLE = {"color": "black", "solid_capstyle": "round"}

SCATTER_STYLE = {"s": 15}

PYOP2_COLOR = "blue"
PYOP2_SCATTER_STYLE = {
    "color": PYOP2_COLOR,
    "marker": "D",  # diamond
    **SCATTER_STYLE,
}
PYOP2_LINE_STYLE = {
    "color": PYOP2_COLOR,
}

PYOP3_COLOR = "red"
PYOP3_SCATTER_STYLE = {
    "color": PYOP3_COLOR,
    "s": 15,
    "marker": "o",
}
PYOP3_LINE_STYLE = {
    "color": PYOP3_COLOR,
}

XMIN = 1e-1
XMAX = 1e2
YMIN = 5e8
YMAX = 1e11

def plot_roofline(machine_name, pyop2_batchfile, pyop3_batchfile, experiment):
    assert experiment in [1, 2]

    path = pathlib.Path(__file__)

    with open(f"{path.parent.parent}/machines/{machine_name}.toml", "rb") as f:
        machine_info = tomli.load(f)

    peak_bandwidth = machine_info["peak_bandwidth"]
    peak_flops_scalar = machine_info["peak_flops_scalar"]
    peak_flops_avx = machine_info["peak_flops_avx"]

    pyop2_data = pd.read_csv(pyop2_batchfile)
    pyop3_data = pd.read_csv(pyop3_batchfile)

    compute_roofline_data(pyop2_data, peak_bandwidth)
    compute_roofline_data(pyop3_data, peak_bandwidth)

    plt.style.use(MPL_PARAMS)
    plt.figure(**FIGURE_PARAMS)

    plt.xlim(XMIN, XMAX)
    plt.ylim(YMIN, YMAX)
    plt.xscale("log")
    plt.yscale("log")

    plt.xlabel("Arithmetic intensity (FLOP/byte)")
    plt.ylabel("Arithmetic throughput (FLOP/s)")

    roofline_corner_scalar = calc_roofline_corner(peak_bandwidth, peak_flops_scalar)
    roofline_corner_avx = calc_roofline_corner(peak_bandwidth, peak_flops_avx)

    # peak bandwidth
    plt.plot(
        [XMIN, roofline_corner_avx],
        [XMIN*peak_bandwidth, roofline_corner_avx*peak_bandwidth],
        **ROOFLINE_LINE_STYLE,
    )
    # peak throughput (scalar)
    plt.plot([roofline_corner_scalar, XMAX], [peak_flops_scalar, peak_flops_scalar], **ROOFLINE_LINE_STYLE)
    # peak throughput (vector)
    plt.plot([roofline_corner_avx, XMAX], [peak_flops_avx, peak_flops_avx], **ROOFLINE_LINE_STYLE)

    plot_roofline_points(pyop2_data, mode="pyop2", experiment=experiment)
    plot_roofline_points(pyop3_data, mode="pyop3", experiment=experiment)


def plot_roofline_points(data, mode, experiment):
    if mode == "pyop2":
        scatter_style = PYOP2_SCATTER_STYLE
        line_style = PYOP2_LINE_STYLE
    else:
        assert mode == "pyop3"
        scatter_style = PYOP3_SCATTER_STYLE
        line_style = PYOP3_LINE_STYLE

    pessimal_ai = data["pessimal_ai"]
    optimal_ai = data["optimal_ai"]
    flop_rate = data["flop_rate"]

    if experiment == 1:
        plt.scatter(pessimal_ai, flop_rate, **scatter_style)
        plt.scatter(optimal_ai, flop_rate, **scatter_style)
        for xmin, xmax, y in zip(pessimal_ai, optimal_ai, flop_rate):
            plt.plot([xmin, xmax], [y, y], **line_style)
    else:
        assert experiment == 2
        plt.scatter((pessimal_ai+optimal_ai)/2, flop_rate, **scatter_style)


def arithmetic_intensity(nflops, memory):
    return nflops / memory


def flop_rate(nflops, time):
    return nflops / time


def peak_flop_rate(arithmetic_intensity, bandwidth):
    """
    arithmetic_intensity: FLOP/byte
    bandwidth: byte/s
    """
    return arithmetic_intensity * bandwidth


def compute_roofline_data(batch_info, peak_bandwidth):
    total_time = batch_info["total_time"]
    ncalls = batch_info["count"]
    average_time = total_time / ncalls
    nflops = batch_info["nflops"]
    pessimal_memory = batch_info["pessimal_memory"]

    batch_info["flop_rate"] = flop_rate(nflops, average_time)

    batch_info["optimal_ai"] = arithmetic_intensity(nflops, batch_info["optimal_memory"])
    batch_info["optimal_peak_flops"] = peak_flop_rate(batch_info["optimal_ai"], peak_bandwidth)

    batch_info["pessimal_ai"] = arithmetic_intensity(nflops, pessimal_memory)
    batch_info["pessimal_peak_flops"] = peak_flop_rate(batch_info["pessimal_ai"], peak_bandwidth)


def calc_roofline_corner(bandwidth, peak_flops):
    # Lines meet when bandwidth * x = peak_flops
    return peak_flops / bandwidth

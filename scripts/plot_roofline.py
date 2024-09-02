#! /bin/env python

import math

import matplotlib.pyplot as plt
import typer


DEFAULT_WIDTH = 345.0 * 0.95 # document width in points
DEFAULT_HEIGHT = DEFAULT_WIDTH * 0.8

PEAK_SCALAR = 2
PEAK_VECTOR = 16
MAIN_MEMORY_BANDWIDTH = 30
CACHE_BANDWIDTH = 100

MIN_ARITHMETIC_INTENSITY = 1/1064
MAX_ARITHMETIC_INTENSITY = 16
MIN_THROUGHPUT = 1/100
MAX_THROUGHPUT = 80


plt.style.use({
    "font.family": "serif",
    "font.size": 10,
    "text.usetex": True,
})


def inches_to_pts(inches):
    return inches * 72


def pts_to_inches(pts):
    return pts / 72


def plot_line(ax, xmin, xmax, ymin, ymax):
    ax.plot(
        [xmin, xmax],
        [ymin, ymax],
        color="black",
        solid_capstyle="round",
    )


def log_diff(low, high):
    return math.log(high) - math.log(low)


def compute_angle(xmin, xmax, ymin, ymax, width, height):
    opp = log_diff(ymin, ymax) / log_diff(MIN_THROUGHPUT, MAX_THROUGHPUT) * height
    adj = log_diff(xmin, xmax) / log_diff(MIN_ARITHMETIC_INTENSITY, MAX_ARITHMETIC_INTENSITY) * width
    return math.degrees(math.atan(opp/adj))


def main(
    outfile: str = typer.Option("", "-o"),
    width: int = DEFAULT_WIDTH,
    height: int = DEFAULT_HEIGHT,
):
    fig, ax = plt.subplots(
        1, 1,
        figsize=(pts_to_inches(width), pts_to_inches(height))
    )

    ax.set_xlabel("Arithmetic intensity (FLOP/byte)")
    ax.set_ylabel("Arithmetic throughput (FLOP/s)")

    ax.set_xscale("log")
    ax.set_yscale("log")

    ax.set_xlim(MIN_ARITHMETIC_INTENSITY, MAX_ARITHMETIC_INTENSITY)
    ax.set_ylim(MIN_THROUGHPUT, MAX_THROUGHPUT)

    # hide tick labels
    ax.set_xticklabels([""]*len(ax.get_xticks()))
    ax.set_yticklabels([""]*len(ax.get_yticks()))

    # cache bandwidth
    xmin = MIN_ARITHMETIC_INTENSITY
    ymin = MIN_ARITHMETIC_INTENSITY * CACHE_BANDWIDTH
    xmax = PEAK_VECTOR / CACHE_BANDWIDTH
    ymax = PEAK_VECTOR
    plot_line(ax, xmin, xmax, ymin, ymax)

    ax.annotate(
        "Peak memory bandwidth",
        (MIN_ARITHMETIC_INTENSITY, ymin),
        (12, 15),
        textcoords="offset points",
        rotation=compute_angle(xmin, xmax, ymin, ymax, width, height),
        rotation_mode="anchor",
        verticalalignment="bottom",
    )

    # peak throughput
    xmin = PEAK_VECTOR / CACHE_BANDWIDTH
    xmax = MAX_ARITHMETIC_INTENSITY
    ymin = ymax = PEAK_VECTOR
    plot_line(ax, xmin, xmax, ymin, ymax)
    ax.annotate(
        "Peak throughput",
        (MAX_ARITHMETIC_INTENSITY, ymax),
        (-90, 2),
        textcoords="offset points",
        verticalalignment="bottom",
    )

    # memory-/compute-bound
    xmin = xmax = PEAK_VECTOR / CACHE_BANDWIDTH
    ymin = 0
    ymax = PEAK_VECTOR
    ax.plot(
        [xmin, xmax],
        [ymin, ymax],
        color="black",
        linestyle="dashed",
        linewidth=1,
    )

    ax.annotate(r"\textit{memory-bound}", (0.01, 0.01), (-20, 30), textcoords="offset points")
    ax.annotate(r"\textit{compute-bound}", (0.01, 0.01), (110, 30), textcoords="offset points")


    if outfile:
        fig.savefig(outfile, backend="pgf")
    else:
        plt.show()


if __name__ == "__main__":
    typer.run(main)

import argparse

import pandas as pd
import tomli

from experiments.common.roofline import *
from experiments.common.roofline import arithmetic_intensity
from experiments.common.roofline import flop_rate
from experiments.common.roofline import peak_flop_rate


parser = argparse.ArgumentParser()
parser.add_argument("--machine", type=str, required=True)
parser.add_argument("--batch", type=str, required=True)
args = parser.parse_args()

with open(args.machine, "rb") as f:
    machine = tomli.load(f)

bandwidth = machine["main_memory"]

batch_info = pd.read_csv(args.batch)
total_time = batch_info["total_time"]
ncalls = batch_info["count"]
average_time = total_time / ncalls
nflops = batch_info["nflops"]
optimal_memory = batch_info["optimal_memory"]
pessimal_memory = batch_info["pessimal_memory"]

ai = arithmetic_intensity(nflops, optimal_memory)

# NOTE: This doesn't include DoF maps etc

flop_rate_ = flop_rate(nflops, average_time)
peak_flops = peak_flop_rate(ai, bandwidth)

aibad = arithmetic_intensity(nflops, pessimal_memory)
peak_flopsbad = peak_flop_rate(aibad, bandwidth)

breakpoint()

print(f"Fraction of peak (optimal):")
print(flop_rate_/peak_flops)

print(f"Fraction of peak (pessimal):")
print(flop_rate_/peak_flopsbad)

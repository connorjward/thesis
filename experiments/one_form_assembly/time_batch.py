import argparse
import itertools
import os

import pandas as pd

import experiments.common as utils
from experiments.one_form_assembly.run import make_function_space


parser = argparse.ArgumentParser()
parser.add_argument("--machine-name", required=True)
parser.add_argument("--mode", required=True)
parser.add_argument("--nfuncs", required=True)
parser.add_argument("--degree", required=True)
parser.add_argument("--ncells", required=True)
args = parser.parse_args()

dir = os.path.dirname(__file__)
flop_data = pd.read_csv(f"{dir}/flops.csv", comment="#")

batch_file = f"{dir}/batch-{args.machine_name}-{args.mode}-{utils.timestamp()}.csv"
with open(batch_file, "w") as f:
    f.write("nfuncs,degree,ncells,name,count,total_time,nflops,optimal_memory,pessimal_memory,ndofs\n")

nfuncss = tuple(map(int, args.nfuncs.split(",")))
degrees = tuple(map(int, args.degree.split(",")))
ncellss = tuple(map(int, args.ncells.split(",")))

for nfuncs, degree, ncells in itertools.product(nfuncss, degrees, ncellss):
    cmd = f"python -m {__package__}.run --nfuncs {nfuncs} --degree {degree} --mode {args.mode} --ncells {ncells}"
    print("Running\n", cmd)
    events = utils.petsclog.profile_script(cmd, "Experiment")
    loop_event, = (ev for ev in events if ev.name == "pyop3_loop")

    fs = make_function_space(ncells, degree)
    optimal_memory = fs.dim() * (nfuncs + 1) * 8
    # if we miss cache every time (ignoring the maps)
    pessimal_memory = (nfuncs + 1) * fs.ufl_domain().num_cells() * fs.finat_element.space_dimension() * 8

    flop_data_single = flop_data[(flop_data["nfuncs"] == nfuncs) & (flop_data["degree"] == degree)]
    nflops = int(flop_data_single["flops_per_cell"].iloc[0])
    nflops *= fs.ufl_domain().num_cells()

    ndofs = fs.dim()

    ncalls = loop_event.count
    # if args.mode == "pyop2":
    #     # PyOP2 records parloops twice as often since they are called twice per invocation
    #     ncalls //= 2

    with open(batch_file, "a") as f:
        f.write(f"{nfuncs},{degree},{ncells},{loop_event.name},{ncalls},{loop_event.total_time},{nflops},{optimal_memory},{pessimal_memory},{ndofs}\n")

print(f"Run complete, timings written to {batch_file}")

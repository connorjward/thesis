import argparse
import subprocess
import os

import pandas as pd
from pyop3.utils import just_one

import experiments.common as utils
from experiments.one_form_assembly.run import make_function_space


parser = argparse.ArgumentParser()
parser.add_argument("--machine-name", required=True)
args = parser.parse_args()

dir = os.path.dirname(__file__)
flop_data = pd.read_csv(f"{dir}/flops.csv", comment="#")

# NOTE: I am not sure what should be varied for this experiment, so just
# run the experiment once for now.

# for now
nfuncs = 3
degree = 1

cmd = f"python -m {__package__}.run --nfuncs {nfuncs} --degree {degree}"
events = utils.petsclog.profile_script(cmd, "Experiment")
loop_event = just_one(ev for ev in events if ev.name == "pyop3_loop")

fs = make_function_space(degree)
optimal_memory = fs.dim() * (nfuncs + 1) * 8
# if we miss cache every time (ignoring the maps)
pessimal_memory = (nfuncs + 1) * fs.ufl_domain().num_cells() * fs.finat_element.space_dimension() * 8

flop_data_single = flop_data[(flop_data["nfuncs"] == nfuncs) & (flop_data["degree"] == degree)]
nflops = int(flop_data_single["flops_per_cell"].iloc[0])

batch_file = f"{dir}/batch-{args.machine_name}-{utils.timestamp()}.csv"

with open(batch_file, "w") as f:
    f.write("name,count,total_time,nflops,optimal_memory,pessimal_memory\n")
    f.write(f"{loop_event.name},{loop_event.count},{loop_event.total_time},{nflops},{optimal_memory},{pessimal_memory}\n")

print(f"Run complete, timings written to {batch_file}")

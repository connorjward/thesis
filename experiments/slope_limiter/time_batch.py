import argparse
import subprocess
import os

from pyop3.utils import just_one

import experiments.common as utils


parser = argparse.ArgumentParser()
parser.add_argument("--machine-name", required=True)
args = parser.parse_args()

# NOTE: I am not sure what should be varied for this experiment, so just
# run the experiment once for now.

cmd = f"python -m {__package__}.run"
events = utils.petsclog.profile_script(cmd, "Experiment")
loop_event = just_one(ev for ev in events if ev.name == "pyop3_loop")

cmd = f"python -m {__package__}.run --memory"
result = subprocess.run(cmd.split(), capture_output=True, text=True)
nflops, optimal_memory, pessimal_memory = map(int, result.stdout.splitlines()[0].split())

dir = os.path.dirname(__file__)
batch_file = f"{dir}/batch-{args.machine_name}-{utils.timestamp()}.csv"

with open(batch_file, "w") as f:
    f.write("name,count,total_time,nflops,optimal_memory,pessimal_memory\n")
    f.write(f"{loop_event.name},{loop_event.count},{loop_event.total_time},{nflops},{optimal_memory},{pessimal_memory}\n")

print(f"Run complete, timings written to {batch_file}")

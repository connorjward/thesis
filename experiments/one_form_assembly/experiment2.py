import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--machine-name", required=True)
parser.add_argument("--mode", required=True)
args = parser.parse_args()

nfuncs = 1
degree = "1,2,3,4,5,6,7"
# ncells = 200
ncells = 100  # for now, use 200 in the real thing
cmd = (
    f"python -m {__package__}.time_batch "
    f"--machine-name {args.machine_name} "
    f"--mode {args.mode} "
    f"--nfuncs {nfuncs} "
    f"--degree {degree} "
    f"--ncells {ncells}"
).split()
subprocess.run(cmd)

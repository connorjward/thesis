import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--machine-name", required=True)
parser.add_argument("--mode", required=True)
args = parser.parse_args()

nfuncs = 1
degree = 1
ncells = "50,100,150,200,250,300"
cmd = (
    f"python -m {__package__}.time_batch "
    f"--machine-name {args.machine_name} "
    f"--mode {args.mode} "
    f"--nfuncs {nfuncs} "
    f"--degree {degree} "
    f"--ncells {ncells}"
).split()
subprocess.run(cmd)

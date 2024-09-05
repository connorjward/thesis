import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--machine-name", required=True)
args = parser.parse_args()

cmd = (
    f"python -m {__package__}.time_batch "
    f"--machine-name {args.machine_name} "
).split()
subprocess.run(cmd)

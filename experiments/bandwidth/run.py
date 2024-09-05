import argparse
import pathlib
import re
import subprocess


MAX_NTHREADS = 32

parser = argparse.ArgumentParser()
parser.add_argument("--mode", required=True)
args = parser.parse_args()

assert args.mode in {"stream", "peakflops"}

data = []
for nthreads in range(1, MAX_NTHREADS+1):
    if args.mode == "stream":
        cmd = f"likwid-bench -t stream -w S0:1GB:{nthreads}"
    else:
        cmd = f"likwid-bench -t peakflops -w S0:1GB:{nthreads}"

    print("Running\n", cmd)
    result = subprocess.run(cmd.split(), capture_output=True, text=True)

    found = False
    for line in result.stdout.splitlines():
        if args.mode == "stream":
            match = re.match(r"MByte/s:\s+(\d+.\d+)", line)
            if match is not None:
                found = True
                bandwidth = float(match.group(1)) * 1e6
                data.append((nthreads, bandwidth))
                break
        else:
            match = re.match(r"MFlops/s:\s+(\d+.\d+)", line)
            if match is not None:
                found = True
                throughput = float(match.group(1)) * 1e6
                data.append((nthreads, throughput))
                break
    assert found

dir = pathlib.Path(__file__).parent

if args.mode == "stream":
    with open(f"{dir}/streamdata.csv", "w") as f:
        f.write("nthreads,bandwidth\n")
        for nthreads, value in data:
            f.write(f"{nthreads},{value}\n")
else:
    with open(f"{dir}/flopdata.csv", "w") as f:
        f.write("nthreads,throughput\n")
        for nthreads, value in data:
            f.write(f"{nthreads},{value}\n")

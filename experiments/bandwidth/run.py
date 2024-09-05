import pathlib
import re
import subprocess


MAX_NTHREADS = 32

data = []
for nthreads in range(1, MAX_NTHREADS+1):
    cmd = f"likwid-bench -t stream -w S0:1GB:{nthreads}"
    print("Running\n", cmd)
    result = subprocess.run(cmd.split(), capture_output=True, text=True)
    found = False
    for line in result.stdout.splitlines():
        match = re.match(r"MByte/s:\s+(\d+.\d+)", line)
        if match is not None:
            found = True
            bandwidth = float(match.group(1)) * 1e6
            data.append((nthreads, bandwidth))
            break
    assert found

dir = pathlib.Path(__file__).parent
with open(f"{dir}/streamdata.csv", "w") as f:
    f.write("nthreads,bandwidth\n")
    for nthreads, bandwidth in data:
        f.write(f"{nthreads},{bandwidth}\n")

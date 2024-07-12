import re
import subprocess


# NFUNCS = [1, 3, 5]
NFUNCS = [1]

for nfuncs in NFUNCS:
    cmd = f"python run.py --nfuncs {nfuncs} -log_view :profile.txt:ascii_flamegraph".split()
    subprocess.run(cmd)

    with open("profile.txt") as f:
        for line in f:
            # match = re.match("\s*;run experiment;\s*;parloop core (\d+)", line)
            match = re.match(r".*;run experiment;.* (\d+)\n", line)
            if match is not None:
                time = int(match.group(1))

# now save nfuncs and time together - repeats are done by run.py
# need to add parloop core event!

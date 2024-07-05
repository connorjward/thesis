import argparse

import numpy as np
from firedrake import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nfuncs", type=int, default=0)
    parser.add_argument("--degree", type=int, default=1)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    mesh = UnitSquareMesh(5, 5)
    V = FunctionSpace(mesh, "P", args.degree)
    v = TestFunction(V)

    if args.nfuncs == 0:
        form = v * dx
    else:
        form = np.prod([Function(V) for _ in range(args.nfuncs)]) * v * dx

    with PETSc.Log.Event("run experiment")
        assemble(form)

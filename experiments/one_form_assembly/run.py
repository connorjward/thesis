import argparse

from firedrake import *


# TODO: Accept likwid option, otherwise do a loop, hot start etc
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nfuncs", type=int, default=0)
    parser.add_argument("--degree", type=int, default=1)
    return parser.parse_known_args()[0]


if __name__ == "__main__":
    args = parse_args()

    mesh = UnitSquareMesh(5, 5)
    V = FunctionSpace(mesh, "P", args.degree)
    v = TestFunction(V)

    if args.nfuncs == 0:
        form = v * dx
    else:
        form = sum(Function(V) for _ in range(args.nfuncs)) * v * dx

    with PETSc.Log.Event("run experiment"):
        assemble(form, form_compiler_parameters={"add_likwid_markers": True})

import argparse

from firedrake import *
from firedrake.assemble import OneFormAssembler


NREPEATS = 3


# TODO: Accept likwid option, otherwise do a loop, hot start etc
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nfuncs", type=int, default=0)
    parser.add_argument("--degree", type=int, default=1)
    parser.add_argument("--debug", action="store_true")
    return parser.parse_known_args()[0]


def make_function_space(degree, *, debug=False):
    if debug:
        mesh = UnitSquareMesh(5, 5)
    else:
        mesh = UnitSquareMesh(50, 50)
    return FunctionSpace(mesh, "P", degree)


if __name__ == "__main__":
    args = parse_args()

    V = make_function_space(args.degree, debug=args.debug)
    v = TestFunction(V)

    if args.nfuncs == 0:
        form = v * dx
    else:
        form = sum(Function(V) for _ in range(args.nfuncs)) * v * dx

    assembler = OneFormAssembler(form, pyop3_compiler_parameters={"add_petsc_event": True})

    # warm start
    assembler.assemble()

    with PETSc.Log.Stage("Experiment"):
        for _ in range(NREPEATS):
            assembler.assemble()

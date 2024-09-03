import argparse
import os

from firedrake import *
from firedrake.assemble import OneFormAssembler


NREPEATS = 3


# TODO: Accept likwid option, otherwise do a loop, hot start etc
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--nfuncs", type=int, default=0)
    parser.add_argument("--degree", type=int, default=1)
    parser.add_argument("--ncells", required=True)
    parser.add_argument("--mode", required=True)
    return parser.parse_known_args()[0]


def make_function_space(ncells, degree):
    mesh = UnitSquareMesh(ncells, ncells)
    return FunctionSpace(mesh, "P", degree)


if __name__ == "__main__":
    args = parse_args()

    V = make_function_space(int(args.ncells), args.degree)
    v = TestFunction(V)

    if args.nfuncs == 0:
        form = v * dx
    else:
        form = sum(Function(V) for _ in range(args.nfuncs)) * v * dx

    if args.mode == "likwid":
        assemble(form, pyop3_compiler_parameters={"add_likwid_markers": True})
    else:
        if args.mode == "pyop3":
            assembler = OneFormAssembler(form, pyop3_compiler_parameters={"add_petsc_event": True})
        else:
            assert args.mode == "pyop2"
            os.environ["PYOP2_ADD_PETSC_EVENT"] = "1"
            assembler = OneFormAssembler(form)

        # warm start
        assembler.assemble()

        with PETSc.Log.Stage("Experiment"):
            for _ in range(NREPEATS):
                assembler.assemble()

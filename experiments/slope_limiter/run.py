# This must be at the top of the file to make sure PETSc logging
# is initialised correctly (else -log_view has no output)
from firedrake import *


import argparse

import loopy as lp
import numpy as np
import pyop3 as op3
from firedrake.petsc import PETSc
from pyop3.ir import LOOPY_LANG_VERSION, LOOPY_TARGET


NREPEATS = 10


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", action="store_true")
    parser.add_argument("--memory", action="store_true")
    parser.add_argument("--degree", type=int, default=1)
    return parser.parse_known_args()[0]


# NOTE: This could be integrated into pyop3 more directly as a scalar function
def make_max_kernel():
    lpy_kernel = lp.make_kernel(
        [],
        "out[0] = out[0] if out[0] > in[0] else in[0]",
        [
            lp.GlobalArg("in", shape=(1,), dtype=ScalarType),
            lp.GlobalArg("out", shape=(1,), dtype=ScalarType),
        ],
        target=LOOPY_TARGET,
        lang_version=LOOPY_LANG_VERSION,
    )
    return op3.Function(lpy_kernel, [op3.READ, op3.RW])


def make_inc_kernel():
    # add 3 because we load both CG and DG store CG
    lpy_kernel = lp.make_kernel(
        [],
        "out[0] = out[0] + 1",
        [
            lp.GlobalArg("out", shape=(1,), dtype=ScalarType),
        ],
        target=LOOPY_TARGET,
        lang_version=LOOPY_LANG_VERSION,
    )
    return op3.Function(lpy_kernel, [op3.INC])


def make_loop_expr(mesh, cg, dg):
    max_ = make_max_kernel()
    return op3.loop(
        v := mesh.vertices.index(),
        op3.loop(
            c := mesh.star(v, k=2).index(),
            max_(dg.dat[c], cg.dat[v]),
        ),
        compiler_parameters={"add_petsc_event": True},
    )


def num_cells_visited(mesh, cg, dg):
    cg = Function(cg.function_space())
    dg = Function(dg.function_space())

    inc = make_inc_kernel()
    op3.do_loop(
        v := mesh.vertices.index(),
        op3.loop(
            mesh.star(v, k=2).index(),
            inc(cg.dat[v]),
        ),
    )
    return int(sum(cg.dat.data_ro))


def pessimal_memory(mesh, cg, dg):
    # We can safely assume that cg is only visited nverts times as it will always
    # remain in cache for the inner loop. Therefore the pessimal memory is
    # 2*nverts (for cg, load+store) + total number of cells visited (including dups.)
    return (2*mesh.num_vertices() + num_cells_visited(mesh, cg, dg)) * 8


def optimal_memory(cg, dg):
    return (cg.function_space().dim() + dg.function_space().dim()) * 8


def flop_count(mesh, cg, dg):
    return num_cells_visited(mesh, cg, dg)


if __name__ == "__main__":
    args = parse_args()

    if args.validate:
        mesh = UnitSquareMesh(1, 1)
    else:
        mesh = UnitSquareMesh(100, 100)

    V_cg = FunctionSpace(mesh, "CG", 1)
    V_dg = FunctionSpace(mesh, "DG", 0)
    cg = Function(V_cg)
    dg = Function(V_dg)

    loop_expr = make_loop_expr(mesh, cg, dg)

    if args.validate:
        # Set the vertex values to the maximum x coordinate of the adjacent cells:
        #
        #    .33 --- .66
        #     |     / |
        #     |   /   |
        #     | /     |
        #    .66 --- .66
        dg.interpolate(mesh.coordinates.sub(0))
        assert np.allclose(sorted(dg.dat.data_ro), [0.33, 0.66], atol=0.01)

        loop_expr()
        assert np.allclose(sorted(cg.dat.data_ro), [0.33, 0.66, 0.66, 0.66], atol=0.01)

        # check pessimal and optimal memory
        assert optimal_memory(cg, dg) == (4 + 2) * 8
        assert pessimal_memory(mesh, cg, dg) == ((2+1) + (2+2) + (2+2) + (2+1)) * 8

        print("Passes checks")

    elif args.memory:
        print(flop_count(mesh, cg, dg), optimal_memory(cg, dg), pessimal_memory(mesh, cg, dg))
    else:
        # warm start
        loop_expr()

        with PETSc.Log.Stage("Experiment"):
            # increasing this will not make things better
            for _ in range(NREPEATS):
                loop_expr()

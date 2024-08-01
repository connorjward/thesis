import argparse

import loopy as lp
import numpy as np
import pyop3 as op3
from firedrake import *
from pyop3.ir import LOOPY_LANG_VERSION, LOOPY_TARGET


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate", action="store_true")
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


def make_loop_expr(mesh, cg, dg):
    max_ = make_max_kernel()
    return op3.loop(
        v := mesh.vertices.index(),
        op3.loop(
            c := mesh.star(v, k=2).index(),
            max_(dg.dat[c], cg.dat[v]),
        ),
    )


if __name__ == "__main__":
    args = parse_args()

    if args.validate:
        mesh = UnitSquareMesh(1, 1)
    else:
        raise NotImplementedError()

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
        print("Passes checks")

    else:
        with PETSc.Log.Event("run experiment"):
            loop_expr()

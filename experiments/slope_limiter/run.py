import argparse

import loopy as lp
import numpy as np
import pyop3 as op3
from firedrake import *


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--degree", type=int, default=1)
    return parser.parse_args()


# NOTE: This could be integrated into pyop3 more directly as a scalar function
def make_max_kernel():
    lpy_kernel = lp.make_kernel(
        "{ [i]: 0 <= i < 1 }",
        "out[0] = max(out[0], in[0])",
        [
            lp.GlobalArg("in", shape=(1,), dtype=ScalarType),
            lp.GlobalArg("out", shape=(1,), dtype=ScalarType),
        ],
        target=lp.CTarget(),
    )
    return op3.Function(lpy_kernel, [op3.READ, op3.RW])


def make_loop_expr():
    mesh = UnitSquareMesh(5, 5)

    V_cg = FunctionSpace(mesh, "CG", 1)
    V_dg = FunctionSpace(mesh, "DG", 0)

    cg = Function(V_cg)
    dg = Function(V_dg)

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
    loop_expr = make_loop_expr()

    with PETSc.Log.Event("run experiment"):
        loop_expr()

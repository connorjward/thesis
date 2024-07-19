import loopy as lp
import pyop3 as op3
from pyrsistent import freeze

from common.utils import artefact_directory, strip_ansi_chars


a = op3.Axis(5)
b = op3.Axis(3)

axes = op3.AxisTree.from_nest({a: b})
datA = op3.HierarchicalArray(axes, dtype=float)
datB = op3.HierarchicalArray(a, dtype=float)

target_axis = a.label
target_component = a.component.label
map_axes = op3.AxisTree.from_nest({a: op3.Axis(2)})
map_dat = op3.HierarchicalArray(map_axes, dtype=int)
mapX = op3.Map({freeze({target_axis: target_component}): [op3.TabulatedMapComponent(target_axis, target_component, map_dat)]})

kernel = op3.Function(
    lp.make_kernel(
        "{ [i]: 0 <= i < 6 }",
        "y[0] = y[0] + x[i]",
        [lp.GlobalArg("x", shape=(6,), dtype=float), lp.GlobalArg("y", shape=(1,), dtype=float)],
        target=op3.ir.LOOPY_TARGET,
        lang_version=op3.ir.LOOPY_LANG_VERSION,
    ),
    [op3.READ, op3.INC],
)

loop = op3.loop(
    p := a.index(),
    kernel(datA[mapX(p)], datB[p]),
)

loopy_kernel = loop.loopy_code.ir

# loopy kernel
with open(f"{artefact_directory()}/codegen_example2_loopy_kernel_default.txt", "w") as f:
    f.write(strip_ansi_chars(str(loopy_kernel)))

# C code
with open(f"{artefact_directory()}/codegen_example2_c_code_default.c", "w") as f:
    f.write(lp.generate_code_v2(loopy_kernel).device_code())

import loopy as lp
import pyop3 as op3
from pyrsistent import freeze

from common.utils import artefact_directory, strip_ansi_chars


axes = op3.AxisTree.from_iterable([5, 3])
dat = op3.HierarchicalArray(axes, dtype=float)

# dat[::2, 1:].assign(666)
loop = op3.loop(op3.AxisTree(op3.Axis(1)).index(), dat[::2, 1:].assign(666, eager=False))

loop()
expected = [0, 666, 666, 0, 0, 0, 0, 666, 666, 0, 0, 0, 0, 666, 666]
assert (dat.data_ro == expected).all()

loopy_kernel = loop.loopy_code.ir

# loopy kernel
with open(f"{artefact_directory()}/codegen_example1_loopy_kernel_default.txt", "w") as f:
    f.write(strip_ansi_chars(str(loopy_kernel)))

# C code
with open(f"{artefact_directory()}/codegen_example1_c_code_default.c", "w") as f:
    f.write(lp.generate_code_v2(loopy_kernel).device_code())

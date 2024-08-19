import loopy as lp
import numpy as np
from pyop2 import op2
from pyop2.codegen.rep2loopy import generate
from pyop2.configuration import target

from common.utils import artefact_directory, strip_ansi_chars


loopy_kernel = lp.make_kernel(
    "{ [i]: 0 <= i < 3 }",  # domains
    "y[i] = 2*x[i]",        # instructions
    [                       # arguments
     lp.GlobalArg("x", shape=(3,),  dtype=float),
     lp.GlobalArg("y", shape=(3,), dtype=float),
    ],
    name="mykernel",
    target=target,
    lang_version=(2018, 2)
)

kernel = op2.Kernel(loopy_kernel, "mykernel")
src_set = op2.Set(10)
target_set = op2.Set(30)

data_set = op2.DataSet(target_set, 2)
dat0 = op2.Dat(data_set)
dat1 = op2.Dat(data_set)

src_to_target_map = op2.Map(src_set, target_set, 3, values=np.zeros((10, 3), dtype=int))

parloop = op2.LegacyParloop(kernel, src_set, dat0(op2.READ, src_to_target_map), dat1(op2.WRITE, src_to_target_map))
global_kernel = generate(parloop.global_kernel.builder)

with open(f"{artefact_directory()}/pyop2_example_loopy_kernel_default.txt", "w") as f:
  f.write(strip_ansi_chars(str(global_kernel)))

with open(f"{artefact_directory()}/pyop2_example_c_code_default.c", "w") as f:
  f.write(lp.generate_code_v2(global_kernel).device_code())

import loopy as lp

from common.utils import artefact_directory


knl = lp.make_kernel(
    "{ [i]: 0 <= i < n }",  # domains
    "y[i] = 2*x[i]",        # instructions
    [                       # arguments
      lp.GlobalArg("x", dtype=float),
      lp.GlobalArg("y", dtype=float),
      lp.ValueArg("n", dtype=int),
    ],
    target=lp.CTarget(),
    lang_version=(2018, 2)
)


with open(f"{artefact_directory()}/loopy_example_c_code_default.c", "w") as f:
    f.write(lp.generate_code_v2(knl).device_code())

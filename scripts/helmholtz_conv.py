from firedrake import *


DEGREES = [1, 2, 3]
NREFS = [4, 5, 6, 7]
FILENAME = "helmholtz_conv.csv"


with open(FILENAME, "w") as f:
    f.write("degree,refinement,norm\n")

for degree in DEGREES:
    for nrefs in NREFS:
        mesh = UnitSquareMesh(2**nrefs, 2**nrefs)
        V = FunctionSpace(mesh, "CG", degree)

        u = TrialFunction(V)
        v = TestFunction(V)

        f = Function(V)
        x, y = SpatialCoordinate(mesh)
        f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))

        a = (inner(grad(u), grad(v)) + inner(u, v)) * dx
        L = inner(f, v) * dx

        u = Function(V)

        solve(a == L, u, solver_parameters={'ksp_type': 'cg', 'pc_type': 'none'})

        f.interpolate(cos(x*pi*2)*cos(y*pi*2))
        norm = errornorm(f, u)

        print(f"degree={degree}, nrefs={nrefs}, norm={norm}")
        with open(FILENAME, "a") as f:
            f.write(f"{degree},{nrefs},{norm}\n")

print("done")

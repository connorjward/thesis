from firedrake import *

# initialise the mesh and function spaces
mesh = UnitSquareMesh(10, 10)
V = VectorFunctionSpace(mesh, "P", 3)
Q = FunctionSpace(mesh, "DP", 2)
W = MixedFunctionSpace([V, Q])

u, p = TrialFunctions(W)
v, q = TestFunctions(W)

# viscosity
nu = Constant(666)

# define lhs and rhs
a = nu * inner(grad(u), grad(v)) * dx - p * div(v) * dx + q * div(u) * dx
L = Form([])  # empty rhs

# construct the boundary condition
g = Constant([666, 666])
bc = DirichletBC(W.sub(0), g, "on_boundary")

# assemble and solve the problem
solution = Function(W)
solve(a == L, solution, bcs=[bc])

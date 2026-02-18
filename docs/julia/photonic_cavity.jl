using Gridap
using GridapGmsh

model = GmshDiscreteModel("mesh.msh")

order = 1

V = TestFESpace(
    model,
    ReferenceFE(nedelec, order),
    #vector_type = Vector{ComplexF64},
    #dirichlet_tags = "PML___None"
)
U = TrialFESpace(V1)

lhs(u, v) = ∇ × u ⋅ ∇ × v
rhs(u, v) = u ⋅ v

assem = Gridap.FESpaces.SparseMatrixAssembler(U, V)
A = assemble_matrix(lhs, assem, U, V)
B = assemble_matrix(rhs, assem, U, V)

vals, vecs = eigs(A, B, sigma = 1 / 1.5)

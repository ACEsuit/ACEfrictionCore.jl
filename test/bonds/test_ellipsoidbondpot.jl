using ACEfrictionCore.ACEbonds, ACEfrictionCore, ACEbase, Test, StaticArrays, LinearAlgebra, JuLIP
using ACEbase.Testing
using ACEfrictionCore.ACEbonds.BondCutoffs: EllipsoidCutoff
using ACEfrictionCore: discrete_jacobi, Rn1pBasis, scal1pbasis, Scal1pBasis,
           evaluate, # evaluate_d, 
           Trig1pBasis, λ, 
           Product1pBasis, Categorical1pBasis, SimpleSparseBasis, 
           SymmetricBasis, Invariant, PIBasis

using ACEfrictionCore.ACEbonds: SymmetricEllipsoidBondBasis
using ACEfrictionCore: polytransform    
##

r0cut = 3.2
rcut = 4.0
zcut = 5.0

cutoff = EllipsoidCutoff(r0cut, rcut, zcut)

maxdeg = 5
maxorder = 3


Bsel = ACEfrictionCore.SparseBasis(; maxorder=maxorder, p = 2, default_maxdeg = maxdeg ) 

# import ACEfrictionCore: filter
# using ACEfrictionCore: Prodb, Onepb,minorder, maxorder

basis = SymmetricEllipsoidBondBasis(Invariant(), Bsel; species=[:Si] );

@show length(basis)
# for s in ACEfrictionCore.get_spec(basis)
#     println(s)
# end

model = ACEfrictionCore.LinearACEModel(basis)
θ = ACEfrictionCore.params(model)
θ = randn(length(θ)) ./ (1:length(θ)).^2
ACEfrictionCore.set_params!(model, θ)

##
using ACEfrictionCore.ACEbonds: ACEBondPotentialBasis, ACEBondPotential

zSi = AtomicNumber(:Si)
_bases = Dict((zSi, zSi) => basis)
_models = Dict((zSi, zSi) => model)
inds = Dict((zSi, zSi) => 1:length(basis))
pot = ACEBondPotential(_models, cutoff)
potbasis = ACEbonds.basis(pot)
# potbasis = ACEBondPotentialBasis(models, inds, cutoff)

@info("Test parameter setter and getter functions")
θ1 = ACEfrictionCore.ACEbonds.params(pot)
ACEfrictionCore.ACEbonds.set_params!(pot,θ1)
println_slim(@test all(θ1 .== ACEfrictionCore.ACEbonds.params(pot)))


at = rattle!(set_pbc!(bulk(:Si, cubic=true) * (2, 2, 1), false), 0.2)


##
@info("Testing energy pot vs energy basis")
println_slim(@test energy(pot, at) ≈ dot(energy(potbasis, at), θ))
# forces(pot, at)

# @info("Testing forces pot vs basis")
# dB = forces(potbasis, at)
# println_slim(@test( sum( θ[n] * dB[n, :] for n = 1:length(θ) ) ≈ forces(pot, at) ))

# @info("Finite-difference forces test")
# println_slim(@test JuLIP.Testing.fdtest(pot, at))

##

# @info("Finite-difference virial test")
# at = set_pbc!(bulk(:Si, cubic=true) * (2, 1, 1), true)
# JuLIP.set_variablecell!(at, true)
# println_slim(@test JuLIP.Testing.fdtest(pot, at))

# at = rattle!(set_pbc!(bulk(:Si, cubic=true) * (2, 1, 1), true), 0.2)
# JuLIP.set_variablecell!(at, true)
# println_slim(@test JuLIP.Testing.fdtest(pot, at))

# @info("Virial pt vs basis")
# at = rattle!(set_pbc!(bulk(:Si, cubic=true) * (2, 2, 2), true), 0.2)
# vB = virial(potbasis, at)
# println_slim(@test( sum( θ[n] * vB[n] for n = 1:length(potbasis) ) ≈ virial(pot, at) ))


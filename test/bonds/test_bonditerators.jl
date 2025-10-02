# very rudimentary test to see check whether BondsIterator and FilteredBondsIterator return the same environments.
using JuLIP
using ACEfrictionCore.ACEbonds
using ACEfrictionCore
using ACEfrictionCore.ACEbonds.BondCutoffs: EllipsoidCutoff, env_transform
using ACEfrictionCore.ACEbonds: bonds
using ACEbase.Testing
using Test
r0cut = 3.2
rcut = 4.0
zcut = 5.0

cutoff = EllipsoidCutoff(r0cut, rcut, zcut)


maxdeg = 5
maxorder = 3


Bsel = ACEfrictionCore.SparseBasis(; maxorder=maxorder, p = 2, default_maxdeg = maxdeg ) 


basis = SymmetricEllipsoidBondBasis(ACEfrictionCore.Invariant(), Bsel; species=[:Si] );



model = ACEfrictionCore.LinearACEModel(basis)
θ = ACEfrictionCore.params(model)
θ = randn(length(θ)) ./ (1:length(θ)).^2
ACEfrictionCore.set_params!(model, θ)

using ACEfrictionCore.ACEbonds: ACEBondPotentialBasis, ACEBondPotential

zSi = AtomicNumber(:Si)
_bases = Dict((zSi, zSi) => basis)
_models = Dict((zSi, zSi) => model)
inds = Dict((zSi, zSi) => 1:length(basis))
pot = ACEBondPotential(_models, cutoff)
potbasis = ACEfrictionCore.ACEbonds.basis(pot)


at = rattle!(set_pbc!(bulk(:Si, cubic=true) * (2, 2, 1), false), 0.2)

E1 = zeros(Float64, length(potbasis))
E1t = zeros(Float64, length(potbasis))

for (i, j, rrij, Js, Rs, Zs) in bonds(at, potbasis)
    global E1, E1t
    local env
    # find the right ace model 
    ace = ACEfrictionCore.ACEbonds._get_model(potbasis, at.Z[i], at.Z[j])
    # transform the euclidean to cylindrical coordinates
    #env = eucl2cyl(rrij, at.Z[i], at.Z[j], Rs, Zs)
    env = env_transform(rrij, at.Z[i], at.Z[j], Rs, Zs, potbasis.cutoff)
    # evaluate 
    ACEfrictionCore.evaluate!(E1t, ace, ACEfrictionCore.ACEConfig(env))
    E1 += E1t 
end

E2 = zeros(Float64, length(potbasis))
E2t = zeros(Float64, length(potbasis))
for (i, j, rrij, Js, Rs, Zs) in bonds(at, potbasis, Array(1:length(at)))
    global E2, E2t
    local env
    # find the right ace model 
    ace = ACEfrictionCore.ACEbonds._get_model(potbasis, at.Z[i], at.Z[j])
    # transform the euclidean to cylindrical coordinates
    #env = eucl2cyl(rrij, at.Z[i], at.Z[j], Rs, Zs)
    env = env_transform(rrij, at.Z[i], at.Z[j], Rs, Zs, potbasis.cutoff)
    # evaluate 
    ACEfrictionCore.evaluate!(E2t, ace, ACEfrictionCore.ACEConfig(env))
    E2 += E2t 
end

println_slim(@test all(E1 .== E2));

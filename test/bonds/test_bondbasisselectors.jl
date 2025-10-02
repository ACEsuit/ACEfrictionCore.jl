using ACEfrictionCore, Test, ACEbase, ACEbase.Testing, StaticArrays
using ACEfrictionCore.Random: rand_rot, rand_refl
using Random: shuffle
##
#=
@info("Rudimentary tests for sparse basis selectors and intersections of such")
r0cut = 2.0
rcut = 1.0
zcut = 2.0
env = ACEfrictionCore.EllipsoidBondEnvelope(r0cut, rcut, zcut;floppy=false, λ= .5)

maxorder = 3
dmaxdeg = 4
Bsel_p2 = ACEfrictionCore.SparseBasis(maxorder; p = 2, default_maxdeg = dmaxdeg) 
Bsel_p1 = ACEfrictionCore.SparseBasis(maxorder; p = 1, default_maxdeg = dmaxdeg) 

Bsel_intesect = ACEfrictionCore.intersect(Bsel_p1,Bsel_p2)
B1p = ACEfrictionCore.Utils.RnYlm_1pbasis(; maxdeg=dmaxdeg)
basis_intersect_inv = ACEfrictionCore.SymmetricBasis(ACEfrictionCore.Invariant(), B1p, Bsel_intesect)
=#

@info("Test bond basis selectors")

r0cut = 2.0
rcut = 1.0
zcut = 2.0
env = ACEfrictionCore.EllipsoidBondEnvelope(r0cut, rcut, zcut;floppy=false, λ= .5)

maxorder = 3
Bsel = ACEfrictionCore.SparseBasis(; maxorder=maxorder, p = 2, default_maxdeg = 4) 

##
@info("Test rotation-equivariance properties")

tol = 10^-14
for property in [ACEfrictionCore.Invariant(), ACEfrictionCore.EuclideanVector(), ACEfrictionCore.EuclideanMatrix()]
    local basis, cfg, B1 
    basis = ACEfrictionCore.Utils.SymmetricBond_basis(property, env, Bsel; )
    @show length(basis)

    rr0 = SVector{3}(rand(Float64,3))
    cfg = [ ACEfrictionCore.State(rr = SVector{3}(rand(Float64,3)), rr0 = rr0,
                    be = rand([:bond,:env])) 
            for _ = 1:10 ] |> ACEConfig
    B1 = ACEfrictionCore.evaluate(basis, cfg)

    println("------------------------------------------------------------")
    @info("Test rotation-equivariance for property $(typeof(property)) with tol = $tol")

    for ntest = 1:30
        Q = rand_refl() * rand_rot()
        Xs2 = ACEfrictionCore.shuffle([ ACEfrictionCore.State(rr = Q * X.rr, rr0 = Q * X.rr0, be = X.be)  for X in cfg.Xs ])
        B2 = ACEfrictionCore.evaluate(basis, ACEConfig(Xs2))
        if property == ACEfrictionCore.Invariant()
            print_tf(@test isapprox(B1, B2, rtol=tol))
        elseif property == ACEfrictionCore.EuclideanVector()
            print_tf(@test isapprox( map(x->Q*x, B1), B2, rtol=tol))
        elseif property == ACEfrictionCore.EuclideanMatrix()
            print_tf(@test isapprox( map(x->Q * x * Q', B1), B2, rtol=tol))
        end
    end
    println()
end


env0 = ACEfrictionCore.EllipsoidBondEnvelope(r0cut, rcut, zcut; p0=1, pr=1, floppy=false, λ= 0.0)

maxorder = 2
Bsel = ACEfrictionCore.SparseBasis(;maxorder=maxorder, p = 2, default_maxdeg = 5) 

RnYlm = ACEfrictionCore.Utils.RnYlm_1pbasis(;   r0 = ACEfrictionCore.cutoff_radialbasis(env), 
                                           rin = 0.0,
                                           trans = PolyTransform(1, ACEfrictionCore.cutoff_radialbasis(env)), 
                                           pcut = 0,
                                           pin = 0, Bsel = Bsel
                                       )

function get_config(X::Vector{SVector{3, Float64}}, k::Int, j::Int)
    n = length(X)
    Js = [i for i = 1:n if i != k]
    Rs = [X[i] - X[k] for i = 1:n if i != k]
    ba = (j=j,r= X[j] - X[k])
    config = [ ACEfrictionCore.State(rr = (j==ba.j ? ba.r :  r-.5 * ba.r), rr0 = ba.r, be = (j==ba.j ? :bond : :env ))  for (j,r) in zip(Js, Rs)] |> ACEConfig
    return config
end

##

n_particle = 11
tol = 10^-14
#env0 = ACEfrictionCore.EllipsoidBondEnvelope(r0cut, rcut, zcut;floppy=false, λ= 0.0)
@info("Test equivariance properties under bond symmetry constraints");
for bs in ["Invariant", "Covariant"]
    for property in [ACEfrictionCore.Invariant(), ACEfrictionCore.EuclideanVector(), ACEfrictionCore.EuclideanMatrix()]
        local basis, X, cfg, B1
        k= rand(1:n_particle)
        j= rand([ i for i =1:n_particle if i != k])
        basis = ACEfrictionCore.Utils.SymmetricBond_basis(property, env0, Bsel; bondsymmetry=bs )
        @show length(basis)

        rr0 = SVector{3}(rand(Float64,3))
        X = [ SVector{3}(rand(Float64,3)) for _ = 1:n_particle ]
        cfg = get_config(X, k, j)
        B1 = ACEfrictionCore.evaluate(basis, cfg)

        println("------------------------------------------------------------")
        @info( string("Test rotation-equivariance for property $(typeof(property)) and bond symmetry of type ", bs))

        for ntest = 1:30
            Q = rand_refl() * rand_rot()
            X2 = [Q * x for x in X]
            cfg2 = get_config(X2, k, j)
            B2 = ACEfrictionCore.evaluate(basis, cfg2)
            if property == ACEfrictionCore.Invariant()
                print_tf(@test isapprox(B1, B2, rtol=tol))
            elseif property == ACEfrictionCore.EuclideanVector()
                print_tf(@test isapprox( map(x->Q*x, B1), B2, rtol=tol))
            elseif property == ACEfrictionCore.EuclideanMatrix()
                print_tf(@test isapprox( map(x->Q * x * Q', B1), B2, rtol=tol))
            end
        end
        println()
    end
end

##
@info("Test bond-symmetry conditions")

tol = 10^-10
for property in [ACEfrictionCore.Invariant(), ACEfrictionCore.EuclideanVector(), ACEfrictionCore.EuclideanMatrix()]

    println("------------------------------------------------------------")
    @info("Test bond-symmetry conditions for property $(typeof(property)) and tolerance tol = $tol")
    println()
    basis_bondinv = ACEfrictionCore.Utils.SymmetricBond_basis(property, env0, Bsel; bondsymmetry="Invariant");
    #@show length(basis_bondinv );

    @info("Test for invariance under bond inversion");
    for ntest = 1:300
        local B1
        rr0 = SVector{3}(rand(Float64,3));
        randX = [ SVector{3}(rand(Float64,3)) for _ =1:10];
        cfg1 = vcat( [ ACEfrictionCore.State(rr = r-.5*rr0, rr0 =  rr0, be = :env) for r in randX ] , [ACEfrictionCore.State(rr =  rr0, rr0 =  rr0, be = :bond)]) |> ACEConfig;
        cfg2 = vcat( [ ACEfrictionCore.State(rr = r-.5*rr0, rr0 = -rr0, be = :env) for r in randX ] , [ACEfrictionCore.State(rr = -rr0, rr0 = -rr0, be = :bond)]) |> ACEConfig;
        B1 = ACEfrictionCore.evaluate(basis_bondinv , cfg1);
        B2 = ACEfrictionCore.evaluate(basis_bondinv, cfg2);
        print_tf(@test all([ isapprox(b1.val,b2.val, rtol = tol) for (b1,b2) in zip(B1,B2)]));
    end
    println()

    
    @info("Test for covariance under bond inversion");
    basis_bondcov = ACEfrictionCore.Utils.SymmetricBond_basis(property, env0, Bsel; bondsymmetry="Covariant");
    #@show length(basis_bondcov );
    for ntest = 1:300
        local B1
        rr0 = SVector{3}(rand(Float64,3));
        randX = [ SVector{3}(rand(Float64,3)) for _ =1:10];
        cfg1 = vcat( [ ACEfrictionCore.State(rr = r-.5*rr0, rr0 =  rr0, be = :env) for r in randX ] , [ACEfrictionCore.State(rr =  rr0, rr0 =  rr0, be = :bond)]) |> ACEConfig;
        cfg2 = vcat( [ ACEfrictionCore.State(rr = r-.5*rr0, rr0 = -rr0, be = :env) for r in randX ] , [ACEfrictionCore.State(rr = -rr0, rr0 = -rr0, be = :bond)]) |> ACEConfig;
        B1 = ACEfrictionCore.evaluate(basis_bondcov, cfg1);
        B2 = ACEfrictionCore.evaluate(basis_bondcov, cfg2);
        print_tf(@test all([ isapprox(b1.val,-b2.val, rtol = tol) for (b1,b2) in zip(B1,B2)]));
    end
    println()
    

end








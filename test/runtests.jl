

using ACEfrictionCore, Test, Printf, LinearAlgebra, StaticArrays, BenchmarkTools

##
@testset "ACEfrictionCore.jl" begin
    # --------------------------------------------
    #   basic polynomial basis building blocks
    @testset "Ylm" begin include("polynomials/test_ylm.jl") end
    @testset "TestWigner" begin include("testing/test_wigner.jl") end
    @testset "Transforms" begin include("transforms/test_transforms.jl") end
    @testset "Lambdas" begin include("transforms/test_lambdas.jl") end
    @testset "OrthogonalPolynomials" begin include("polynomials/test_orthpolys.jl") end

    # --------------------------------------------
    #  states, .. 
    @testset "States" begin include("test_states.jl") end 

    # --------------------------------------------
    # core permutation-invariant functionality
    @testset "1-Particle Basis"  begin include("test_1pbasis.jl") end
    @testset "Categorical1pBasis" begin include("test_discrete.jl") end
    @testset "PIBasis" begin include("test_pibasis.jl") end

    # --------------------------------------------
    #   O(3) equi-variance
    @testset "Clebsch-Gordan" begin include("test_cg.jl") end
    @testset "SymmetricBasis" begin include("test_symmbasis.jl") end
    @testset "EuclideanVector" begin include("test_euclvec.jl") end
    @testset "EuclideanMatrix" begin include("test_EuclideanMatrix.jl") end
    @testset "Multiple SH Bases" begin include("test_multish.jl") end

    # Model tests

    @testset "LinearACEModel"  begin include("test_linearmodel.jl") end
    @testset "MultipleProperties"  begin include("test_multiprop.jl") end
    # MS: removed autodiff tests for now. Auto-diff breaks after Julia 1.7 -> 1.11 upgrade and is not needed for ACEfriction.jl. Too much hassle to maintain here.
    # @testset "AD-LinearACEModel"  begin include("test_admodel.jl") end 

    # Experimental material
    @testset "Sparsification" begin include("test_sparsify.jl") end 
    # @testset "Multipliers" begin include("test_multiplier.jl") end

    @testset "Bonds basics" begin include("./bonds/test_bonds.jl"); end
    @testset "Invariant Cylindrical" begin include("./bonds/test_invcyl.jl"); end
    @testset "Bond Iterators" begin include("./bonds/test_bonditerators.jl"); end
    @testset "Calculator (Cylindrical env.)" begin include("./bonds/test_cylindricalbondpot.jl"); end
    @testset "Calculator (Ellipsoid env.)" begin include("./bonds/test_ellipsoidbondpot.jl"); end
end



    # -----------------------------------------
    #    old tests to be re-introduced - maybe
    # include("test_real.jl")
    # include("test_orth.jl")
    # include("bonds/test_cylindrical.jl")
    # include("bonds/test_fourier.jl")

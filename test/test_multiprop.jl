using LinearAlgebra: length
using ACEfrictionCore, ACEbase, ACEfrictionCore.ACEbase024, Test, ACEfrictionCore.Testing
using ACEfrictionCore: evaluate, SymmetricBasis, PIBasis, O3, rand_vec3, rand_radial
using ACEbase.Testing: println_slim
using StaticArrays
using LinearAlgebra: isdiag


##
    
@info(" Testset for Multiple Properties in a Linear ACEModel")

# construct the 1p-basis
maxdeg = 6
ord = 3
Bsel = SimpleSparseBasis(ord, maxdeg)

B1p = ACEfrictionCore.Utils.RnYlm_1pbasis(; maxdeg=maxdeg)

# generate a configuration
#TODO check if this tests account/test u, and will it work with more things?
nX = 30
Xs = () -> ACEfrictionCore.State(rr = rand_vec3(B1p["Rn"]), u = rand())
cfg = ACEConfig([Xs() for i in 1:nX])

for φ in [ACEfrictionCore.Invariant(),ACEfrictionCore.EuclideanVector(Float64),ACEfrictionCore.EuclideanMatrix(Float64), ACEfrictionCore.SymmetricEuclideanMatrix(Float64)]
    @info("Test for $(typeof(φ))-valued property.")
    basis = SymmetricBasis(φ, B1p, O3(), Bsel)

    ##

    BB = evaluate(basis, cfg)
    # dBB = ACEfrictionCore.evaluate_d(basis, cfg)

    c_m = rand(SVector{3,Float64}, length(BB))

    ##

    singlProp = [ACEfrictionCore.LinearACEModel(basis, rand(length(BB)), evaluator = :standard) for i in 1:length(c_m[1])]
    multiProp = ACEfrictionCore.LinearACEModel(basis, c_m, evaluator = :standard)

    ##

    @info("set_params!")
    c_s = [[c_m[j][i] for j in 1:length(c_m)] for i in 1:length(c_m[1])]

    ACEfrictionCore.set_params!(multiProp,c_m)
    for i in 1:length(c_m[1])
        ACEfrictionCore.set_params!(singlProp[i],c_s[i])
        print_tf(@test(c_s[i] ≈ singlProp[i].c))
    end
    println()

    ##

    # We compare a model with a sinlge property, to the corresponding solution
    # of the multiple property model. We do this by giving the same parameters
    # to both models and then using isapprox 

    @info("evaluate")

    for i in 1:length(c_m[1])
        print_tf(@test(evaluate(singlProp[i],cfg).val ≈ evaluate(multiProp,cfg)[i].val))
    end
    println()

    # ##

    # @info("grad_params")

    # multiGradP = ACEfrictionCore.grad_params(multiProp,cfg)

    # println_slim(@test all(isdiag, multiGradP))

    # for i in 1:length(c_m[1])
    #     singl = ACEfrictionCore.grad_params(singlProp[i], cfg)
    #     mult_i = [ multiGradP[j][i,i] for j = 1:length(c_m)] 
    #     print_tf(@test(singl ≈ mult_i))
    # end
    # println()

    # ##

    # @info("grad_config")

    # mgcfg = ACEfrictionCore.grad_config(multiProp, cfg)

    # for i in 1:length(c_m[1])
    #     singl = ACEfrictionCore.grad_config(singlProp[i],cfg)
    #     mult_i = [ g[i] for g in mgcfg ]
    #     print_tf(@test(singl ≈ mult_i))
    # end
    # println()
end
##

# @info("adjoint_EVAL_D 1 prop")
# #we check by contracting the full matrix that adjoint_eval_config works
# #this assumes that grad_params_config works for one porperty. This is 
# #tested elsewhere. 

# i = 2

# for i in 1:length(c_m[1])

#     #find the jacobian
#     Jac = ACEfrictionCore.grad_params_config(singlProp[i],cfg)

#     #create a random input emulating the pullback input
#     w = rand(SVector{3, Float64}, length(Jac[1,:]))
#     w = [ACEfrictionCore.DState(rr = w[j], u = 0.0) for j in 1:length(w)]

#     #calculate the adjoint and make a zeros to fill
#     grad = ACEfrictionCore.adjoint_EVAL_D(singlProp[i], cfg, w)
#     Jgrad = zeros(length(grad))

#     #contract w into the jacobian to get the solution
#     for j in 1:length(Jac[:,1])
#         Jgrad[j] = sum([ACEfrictionCore.contract(Jac[j,:][k], w[k]) for k in 1:length(w)])
#     end

#     print_tf(@test(grad ≈ Jgrad))
# end
# println()

# ##

# @info("adjoint_EVAL_D >2 prop")

# #now we check that multiple properties work. For this, like before, we simply
# #compare the single property result to each of the multiple properties. 

# function wMaker()
#     wtmp = rand(SVector{3, Float64}, length(cfg))
#     wtmp = [ACEfrictionCore.DState(rr = wtmp[j], u = 0.0) for j in 1:length(wtmp)]
#     return wtmp
# end
# TDX1 = ACEfrictionCore.DState{NamedTuple{(:rr, :u), Tuple{SVector{3, Float64}, Float64}}}
# wo = Matrix{TDX1}(undef, (54,7))
# wt  = [wMaker() for i in 1:length(c_m[1])]

# for i in 1:length(wt)
#     for j in 1:length(wt[i])
#         wo[j,i] = wt[i][j]
#     end
# end

# multiEval = ACEfrictionCore.adjoint_EVAL_D(multiProp, cfg, wo)

# for i in 1:length(c_m[1])
#     singl = ACEfrictionCore.adjoint_EVAL_D(singlProp[i], cfg, wo[:,i])
#     multi = [multiEval[j][i] for j in 1:length(c_m)]
#     print_tf(@test(singl ≈ multi))
# end
# println() 

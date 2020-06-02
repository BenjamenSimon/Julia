
using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames

using Distances
using LinearAlgebra

test_inf_times = [Inf, 3, Inf, 5, 7, 9, Inf, 1, Inf]
test_rem_times = [Inf, 6, Inf, 8, 9.5, 27, Inf, 16, Inf]
test_inf_inds = [false, true, false, true, true, true, false, true, false]

test_xcoords = [1,1,1,2,2,2,3,3,3]
test_ycoords = [1,2,3,1,2,3,1,2,3]
test_xycoords = [test_xcoords  test_ycoords]

test_distmat = pairwise(Euclidean(), test_xycoords, test_xycoords, dims = 1)
test_betamat = betamatform(test_distmat, [0.002, 0.001], 1.5)


plot(test_xcoords, test_ycoords, line = (:scatter))

prod_part(test_inf_times, test_rem_times, test_betamat, test_inf_inds)
#-23.472138032568875

interval_intersect(test_inf_times, test_rem_times, test_inf_inds)


integral_part(test_inf_times, test_rem_times, test_betamat, test_inf_inds)
#0.31900000000000006

log_likelihood(test_inf_times, test_rem_times, test_betamat)
#-23.791138032568874



pop = length(test_inf_times)
waifw = reduce(hcat, map(x -> (test_inf_times .< x .< test_rem_times), test_inf_times))
#Correct


lambdaj = sum((test_betamat[:, test_inf_inds] .* waifw[:, test_inf_inds]), dims = 1)

I0 = findmin(test_inf_times)[2]

I0 = findmin(test_inf_times[test_inf_inds])[2]

sum(log.(lambdaj[Not(I0)]))




interval_intersect(test_inf_times, test_rem_times, test_inf_inds)

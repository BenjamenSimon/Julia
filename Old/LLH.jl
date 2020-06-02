
using Distributions
using Distances
using LinearAlgebra
using InvertedIndices
using Random
using Plots
gr()

inf_times = [1, 3, Inf, 5, 7, 9, Inf, Inf]
rem_times = [16, 6, Inf, 8, 9.5, 27, Inf, Inf]
Bmat1 = fill(0.002, 8, 8)
infected_inds = [true, true, false, true, true, true, false, false]

function prod_part(inf_times, rem_times, Bmat, infected_inds)
    waifw = reduce(hcat, map(x -> (inf_times .< x .< rem_times), inf_times))

    lambdaj = sum((Bmat[:, infected_inds] .* waifw[:, infected_inds]), dims = 1)
    #println("lambdaj = $lambdaj")
    I0 = findmin(inf_times[infected_inds])[2]
    sum(log.(lambdaj[Not(I0)]))
end

function interval_intersect(inf_times, rem_times, infected_inds)
    interval_j = [inf_times rem_times]
    interval_i = interval_j[infected_inds, :]

    int_start = reduce(hcat, map(x -> (min.(x, interval_i[:, 1])), interval_j[:, 1]))
    int_end = reduce(hcat, map(x -> (min.(x, interval_i[:, 2])), interval_j[:, 1]))

    int_end - int_start
end

function integral_part(inf_times, rem_times, Bmat, infected_inds)
    E = interval_intersect(inf_times, rem_times, infected_inds)
    integral = E .* Bmat[infected_inds, :]
    sum(integral)
end

function log_likelihood(inf_times, rem_times, Bmat)
    infected_inds = inf_times .< Inf
    prod = prod_part(inf_times, rem_times, Bmat, infected_inds)
    integral = integral_part(inf_times, rem_times, Bmat, infected_inds)
    prod - integral
end

@time begin
    log_likelihood(inf_times, rem_times, Bmat1)
end

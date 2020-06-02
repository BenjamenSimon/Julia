using Distributions
using Distances
using LinearAlgebra
using InvertedIndices
using Random

#####################
### SOPHISTICATED ###
#####################

### (Sophisticated) Product part ###

function sophisticated_prod_part(inf_times, rem_times, Bmat, infected_inds)

    # Calculate the `who acquired infection from whom` matrix
    # Which essentially tells you for each pair of individuals i and j
    # (row i, column j) can individual j have infected individual i
    # given individual j was infected
    waifw = reduce(hcat, map(x -> (inf_times .< x .< rem_times), inf_times))

    # Calculate the value for the sum of the βij over the sets of individuals
    # who were infectious just before individual j became infected.
    # We also refer to this quantity as λj's for each infectd individual j.
    # We output a vector here, a λj value for each infected individuals.
    λj = sum((Bmat[:, infected_inds] .* waifw[:, infected_inds]), dims = 1)

    # Find the initial infected individual as we do not want to include them.
    # We will need their position in the vector of infected indivudals rather
    # than all the individuals/
    I0 = findmin(inf_times[infected_inds])[2]

    # Take the sum of the logs of the λj values.
    int_part = sum(log.(λj[Not(I0)]))

    return(int_part)
end

#Test function: Example 1
infected_individuals_ex1 = findall(res_ex1[:,2] .< Inf)
sophisticated_prod_part(res_ex1[:,2], res_ex1[:,3], betamat_ex1, infected_individuals_ex1)

#Test function: Example 2
infected_individuals_ex2 = findall(res_ex2[:,2] .< Inf)
sophisticated_prod_part(res_ex2[:,2], res_ex2[:,3], betamat_ex2, infected_individuals_ex2)

#Test function: Example 3
infected_individuals_ex3 = findall(res_ex3[:,2] .< Inf)
sophisticated_prod_part(res_ex3[:,2], res_ex3[:,3], betamat_ex3, infected_individuals_ex3)

#Test function: Example 4
infected_individuals_ex4 = findall(res_ex4[:,2] .< Inf)
sophisticated_prod_part(res_ex4[:,2], res_ex4[:,3], betamat_ex4, infected_individuals_ex4)


### (Sophisticated) Integral part ###

#~~# Interval intersect function #~~#
# This function is used to calculate the minimums of the infection and removal
# times and find the appropriate value of their differences

function sophisticated_interval_intersect(inf_times, rem_times, infected_inds)

    # Create a matrix of the difference between the infection times (in i)
    # and infection times (in j) for i in the vector of infected individuals
    # and j in the vector of all individuals.
    int_start = reduce(hcat, map(x -> (min.(x, inf_times[infected_inds])), inf_times))

    # Create a matrix of the difference between the removal times (in i)
    # and infection times (in j) for i in the vector of infected individuals
    # and j in the vector of all individuals.
    int_end = reduce(hcat, map(x -> (min.(x, rem_times[infected_inds])), inf_times))

    return(int_end - int_start)
end

#~~# Interval intersect function #~~#

function sophisticated_integral_part(inf_times, rem_times, Bmat, infected_inds)

    # Calculate the matrix of the times each individual spends in the
    # infectious period of each other individual while susceptible, using
    # the function we defined before.
    E = sophisticated_interval_intersect(inf_times, rem_times, infected_inds)

    # Calculate the values of the appropriate βij values mutliplied by the
    # intervals (function of the minimums).
    integral = E .* Bmat[infected_inds, :]

    # Take the sum of all the values
    return(sum(integral))
end

#Test function: Example 1
infected_individuals_ex1 = findall(res_ex1[:,2] .< Inf)
sophisticated_integral_part(res_ex1[:,2], res_ex1[:,3], betamat_ex1, infected_individuals_ex1)

#Test function: Example 2
infected_individuals_ex2 = findall(res_ex2[:,2] .< Inf)
sophisticated_integral_part(res_ex2[:,2], res_ex2[:,3], betamat_ex2, infected_individuals_ex2)

#Test function: Example 3
infected_individuals_ex3 = findall(res_ex3[:,2] .< Inf)
sophisticated_integral_part(res_ex3[:,2], res_ex3[:,3], betamat_ex3, infected_individuals_ex3)

#Test function: Example 4
infected_individuals_ex4 = findall(res_ex4[:,2] .< Inf)
sophisticated_integral_part(res_ex4[:,2], res_ex4[:,3], betamat_ex4, infected_individuals_ex4)


### (Sophisticated) Log-likelihood ###

function sophisticated_log_likelihood(inf_times, rem_times, Bmat)

    # Infected individuals
    infected_inds = findall(inf_times .< Inf)

    # Calculate the product part
    prod = sophisticated_prod_part(inf_times, rem_times, Bmat, infected_inds)

    # Calculate the integral part
    integral = sophisticated_integral_part(inf_times, rem_times, Bmat, infected_inds)

    # Calculate the log likelikehood
    return(prod - integral)
end

#Test function: Example 1
sophisticated_log_likelihood(res_ex1[:,2], res_ex1[:,3], betamat_ex1)

#Test function: Example 2
sophisticated_log_likelihood(res_ex2[:,2], res_ex2[:,3], betamat_ex2)

#Test function: Example 3
sophisticated_log_likelihood(res_ex3[:,2], res_ex3[:,3], betamat_ex3)

#Test function: Example 4
sophisticated_log_likelihood(res_ex4[:,2], res_ex4[:,3], betamat_ex4)

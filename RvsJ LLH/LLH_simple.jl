using Distributions
using Distances
using LinearAlgebra
using InvertedIndices
using Random

##############
### SIMPLE ###
##############

### (Simple) Product part ###

function simple_product_part(infection_times, removal_times, beta_matrix)
    #=
        This function will calculate the product part of the equation of interest
        in its logarithmic form. To do this we will need
        * Vector of IDs of infected individuals
        * ID of intial infected
        * For each infected individual, sum up the infectious pressure exerted
          on them just before they became infected
            - Thus we will need a vector for each infected individual j of the
              individuals who were infectious just before j became infected
            - We call this vector a `who acquired infection from whom` vector
              (waifw)
        * We then sum up the logs of all of these infectious pressures (λj's)
          for every individual except the initial infected
            - As there was no infectious pressure prior to the initial infection.
    =#

    # Vector of indexs of infected individuals in infection times vector
    infected_inds = findall(infection_times .< Inf)

    # Vector of infection times for infected individuals only
    infection_times_infeds = infection_times[infected_inds]

    # Index of the initial infected in the vector of infected individuals
    I0 = findmin(infection_times_infeds)[2]

    #Vector of indexs of infected individuals in infection times vector
    # Excluding initial infected
    infected_inds_exc_I0 = infected_inds[Not(I0)]

    ## `Who acquired infection from whom` vectors
    ## Which individuals are infectious just before individual j becomes infected

        # Initialise the vector of vectors, one for each infected
        waifw_vectors = Array{Array{Int64}}(undef, 1, length(infected_inds))

        #=
            For each infected individual j, find the IDs in the infection times
            vector of all the indivudals who have infection times earlier than
            the infection time of individual j, and removal time later than
            the infection time of individual j.
        =#
        index = 1
        for j in infected_inds
            waifw_vectors[index] = findall(infection_times .< infection_times[j] .< removal_times)
            index = index + 1
        end
        # Output: Vector of vectors. One vector for each infected individual containing
        #         a list of the individual ids who were infectious just before they got
        #         infected.

    ## Sum up the infectious pressure (βij) exerted on each infected individual
    ## just before they became infected

        # Intialise the vector to record the λjs
        λj_vec = Array{Float64}(undef, 1, length(infected_inds))

        #=
            For each individual find their λj value using the waifw vectors
            that we have calculated in the previous part.
        =#
        index = 1
        for j in infected_inds
            βij_vec = deepcopy(beta_matrix[j, waifw_vectors[index]])
            λj_vec[index] = sum(βij_vec)
            index = index + 1
        end

    # Take the logs of the λjs
    logλjs = log.(λj_vec)

    # The log of the product term that we are interested in is then given by
    # the sum of the λjs excluding the initial infected.
    product = sum(logλjs[Not(I0)])
end

#Test function: Example 1
simple_product_part(res_ex1[:,2], res_ex1[:,3], betamat_ex1)

#Test function: Example 2
simple_product_part(res_ex2[:,2], res_ex2[:,3], betamat_ex2)

#Test function: Example 3
simple_product_part(res_ex3[:,2], res_ex3[:,3], betamat_ex3)

#Test function: Example 4
simple_product_part(res_ex4[:,2], res_ex4[:,3], betamat_ex4)


### (Simple) Integral part ###

function simple_integral_part(infection_times, removal_times, beta_matrix)

    #=
        This function will calculate the integral part of the equation of interest
        in its logarithmic form, using the double sum of minimums form.
        To do this we will need
        * Vector of IDs of infected individuals
        * For each infected individual i
            - The minimum of their removal times and the infection times of all
              individuals j
            - The minimum of their infection times and the infection times of all
              individuals j
            - The difference of these 2 values for all pairs of i and j
            - ...multiplied by the appropriate infection rate βij
        * The sum of all of these values
    =#

    # Vector of indexs of infected individuals in infection times vector
    infected_inds = findall(infection_times .< Inf)

    # Initialise a vector of vectors for calculating minimums, one for each infected
    min_vectors = Array{Array{Float64}}(undef, 1, length(infected_inds))

    # Calculate the minimum of each of the infected individuals infection
    # and removal times with the infection times of all individuals, respectively.
    # Then take the difference.
    index = 1
    for i in infected_inds
        minRI = min.(removal_times[i], infection_times)
        minII = min.(infection_times[i], infection_times)

        min_vectors[index] = minRI - minII

        index = index + 1
    end

    # Initialise another vector of vectors for calculating product of minimums
    # with the appropriate infection rates, one for each infected
    βtimesmin_vectors = Array{Array{Float64}}(undef, 1, length(infected_inds))

    # For each minimum calculation between i and j, multiply it by the appropriate
    # beta value from the infection rate matrix
    index = 1
    for i in infected_inds
        βtimesmin_vectors[index] = beta_matrix[i,:] .* min_vectors[index]
        index = index + 1
    end

    # Calculate the value of the double sum
    sumsum = sum(sum(βtimesmin_vectors))

    return(sumsum)
end

#Test function: Example 1
simple_integral_part(res_ex1[:,2], res_ex1[:,3], betamat_ex1)

#Test function: Example 2
simple_integral_part(res_ex2[:,2], res_ex2[:,3], betamat_ex2)

#Test function: Example 3
simple_integral_part(res_ex3[:,2], res_ex3[:,3], betamat_ex3)

#Test function: Example 4
simple_integral_part(res_ex4[:,2], res_ex4[:,3], betamat_ex4)

### (Simple) Log-likelihood ###

function simple_log_likelihood(infection_times, removal_times, beta_matrix)

  # Calculate the product part
  prod = simple_product_part(infection_times, removal_times, beta_matrix)

  # Calculate the integral part
  integral = simple_integral_part(infection_times, removal_times, beta_matrix)

  # Calculate the log likelikehood
  return(prod - integral)
end

#Test function: Example 1
simple_log_likelihood(res_ex1[:,2], res_ex1[:,3], betamat_ex1)

#Test function: Example 2
simple_log_likelihood(res_ex2[:,2], res_ex2[:,3], betamat_ex2)

#Test function: Example 3
simple_log_likelihood(res_ex3[:,2], res_ex3[:,3], betamat_ex3)

#Test function: Example 4
simple_log_likelihood(res_ex4[:,2], res_ex4[:,3], betamat_ex4)

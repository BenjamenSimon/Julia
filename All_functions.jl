using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames

using Distances
using LinearAlgebra

## Distance Matrix function ##

function unifdistmat(N, xlim, ylim)
    xcoords = xlim * rand(N)
    ycoords = ylim * rand(N)

    xcoords = reshape(xcoords, :, 1)
    ycoords = reshape(ycoords, :, 1)
    xycoords = [xcoords  ycoords]

    distmat = pairwise(Euclidean(), xycoords, xycoords, dims = 1)

    return [xycoords, distmat]
end

## Beta Matrix function ##

function betamatform(distmat, betas, d)

    betamat = Array{Float64}(undef, size(distmat))

    index = findall(distmat .< d)
    betamat[index] .= betas[1]

    index = findall(distmat .>= d)
    betamat[index] .= betas[2]

    betamat[diagind(betamat)] .= 0

    return betamat
end

## GSE function ##

### Initialise the "Time until contact" matrix ###

function generatetuc(betamat)
    popn = size(betamat, 1)
    W = Array{Float64}(undef, popn, popn)
    uniqueness = popn*popn-popn+1

    while size(unique(W)) != Tuple{Int64}(uniqueness)
        for i = 1:popn
            for j = (1:popn)[Not(i)]
                global W[i,j] = rand(Exponential(1/(betamat[i,j])))
            end
        end
        W[diagind(W)] .= Inf

        println("Initialised")
    end

    return W
end

## Simulaion function ##

function GSEsim(N, betamat, gamma)

    ### Initialise the "Time until contact" matrix ###

    Wmat = generatetuc(betamat)

    ### Initialise the "Infectious period" matrix ###

    Q = rand(Exponential(1/gamma), N)

    Qmat = reshape(repeat(Q, outer = N), N, N)

    ### Initialise the "infection times" vector ###

    inftimesfull = [0 ; fill(Inf, (N-1))]

    ### Initialise functional objects ###

    pop = 1:N
    infect = [1]
    sus = pop[Not(infect)]

    ## The simulation ##

    STOP = 0

    while STOP == 0
        ### Calculate the current submatrices that are relevant ###

        Wcur = Wmat[infect, sus]
        Qcur = Qmat[infect, sus]

        ### Calculate which individuals can become infected ###

        within = Wcur .< Qcur
        whichinf = findall(within .== true)
        whoinf = getindex.(whichinf, [1 2])
        track = [ infect[whoinf[:,1]]  sus[whoinf[:,2]] ]

        ### Calculate the infection times of all potential infectees

        inftimes = inftimesfull[track[:,1]] + Wcur[whichinf]

        ### Recover the index of the newly infected individual ###

        if length(inftimes) != 0
            newinftime, up = findmin(inftimes)
            newinf = track[up, :][2]

        ### Update population ###

            inftimesfull[newinf] = newinftime

            infect = [infect ; newinf]
            sus = pop[Not(infect)]
        else
            STOP = 1
        end
    end

    ## Record the results ##

    remtimesfull = inftimesfull + Q

    res = [1:N inftimesfull remtimesfull Q]

    return(res)
end

## Simulation function (choose initial infected) ##

function GSEsim(N, betamat, gamma, init_inf)

    ### Initialise the "Time until contact" matrix ###

    Wmat = generatetuc(betamat)

    ### Initialise the "Infectious period" matrix ###

    Q = rand(Exponential(1/gamma), N)

    Qmat = reshape(repeat(Q, outer = N), N, N)

    ### Initialise the "infection times" vector ###

    inftimesfull = fill(Inf, (N))
    inftimesfull[init_inf] = 0

    ### Initialise functional objects ###

    pop = 1:N
    infect = [init_inf]
    sus = pop[Not(infect)]

    ## The simulation ##

    STOP = 0

    while STOP == 0
        ### Calculate the current submatrices that are relevant ###

        Wcur = Wmat[infect, sus]
        Qcur = Qmat[infect, sus]

        ### Calculate which individuals can become infected ###

        within = Wcur .< Qcur
        whichinf = findall(within .== true)
        whoinf = getindex.(whichinf, [1 2])
        track = [ infect[whoinf[:,1]]  sus[whoinf[:,2]] ]

        ### Calculate the infection times of all potential infectees

        inftimes = inftimesfull[track[:,1]] + Wcur[whichinf]

        ### Recover the index of the newly infected individual ###

        if length(inftimes) != 0
            newinftime, up = findmin(inftimes)
            newinf = track[up, :][2]

        ### Update population ###

            inftimesfull[newinf] = newinftime

            infect = [infect ; newinf]
            sus = pop[Not(infect)]
        else
            STOP = 1
        end
    end

    ## Record the results ##

    remtimesfull = inftimesfull + Q

    res = [1:N inftimesfull remtimesfull Q]

    return(res)
end

## Seed finder function ##

function seedfinder(N, betas, ds, tinf, Δ, seeds)
    sizes = Array{Any}(undef, seeds)

    for i = 1:seeds
        Random.seed!(i)
        distmat1 = unifdistmat(N, 20, 20)[2]
        betamat1 = betamatform(distmat1, betas, ds)
        results = GSEsim(N, betamat1, 0.15)
        global sizes[i] = length(findall(results[:,2] .< Inf))
    end

    display(histogram(sizes))

    whichseeds = findall(x -> x > (tinf - Δ) && x < (tinf + Δ), sizes)
    return([whichseeds sizes[whichseeds]])
end

## Seed finder function (choose initial infected) ##

function seedfinder(N, betas, ds, tinf, Δ, seeds, init_inf)
    sizes = Array{Any}(undef, seeds)

    for i = 1:seeds
        Random.seed!(i)
        distmat1 = unifdistmat(N, 20, 20)[2]
        betamat1 = betamatform(distmat1, betas, ds)
        results = GSEsim(N, betamat1, 0.15, init_inf)
        global sizes[i] = length(findall(results[:,2] .< Inf))
    end

    display(histogram(sizes))

    whichseeds = findall(x -> x > (tinf - Δ) && x < (tinf + Δ), sizes)
    return([whichseeds sizes[whichseeds]])
end


##~ Log-likelihood function ~##

### Product part function ###

function prod_part(inf_times, rem_times, Bmat, infected_inds)
    pop = length(inf_times)
    waifw = reduce(hcat, map(x -> (inf_times .< x .< rem_times), inf_times))

    lambdaj = sum((Bmat[:, infected_inds] .* waifw[:, infected_inds]), dims = 1)
    I0 = findmin(inf_times[infected_inds])[2]
    sum(log.(lambdaj[Not(I0)]))
end

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

### Interval intersect function ###

function interval_intersect(inf_times, rem_times, infected_inds)
    interval_j = [inf_times rem_times]
    interval_i = interval_j[infected_inds, :]

    int_start = reduce(hcat, map(x -> (min.(x, interval_i[:, 1])), interval_j[:, 1]))
    int_end = reduce(hcat, map(x -> (min.(x, interval_i[:, 2])), interval_j[:, 1]))

    int_end - int_start
end

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

### Intergral part funciton ###

function integral_part(inf_times, rem_times, Bmat, infected_inds)
    E = interval_intersect(inf_times, rem_times, infected_inds)
    integral = E .* Bmat[infected_inds, :]
    sum(integral)
end

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

### Calculate LLH function ###

function log_likelihood(inf_times, rem_times, Bmat)
    infected_inds = inf_times .< Inf
    prod = prod_part(inf_times, rem_times, Bmat, infected_inds)
    integral = integral_part(inf_times, rem_times, Bmat, infected_inds)
    prod - integral
end

function simple_log_likelihood(infection_times, removal_times, beta_matrix)

  # Calculate the product part
  prod = simple_product_part(infection_times, removal_times, beta_matrix)

  # Calculate the integral part
  integral = simple_integral_part(infection_times, removal_times, beta_matrix)

  # Calculate the log likelikehood
  return(prod - integral)
end

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

### Plot infectious periods to validate epidemic ###

function plot_inf_periods(res)

    inf_inds = findall(res[:,2] .< Inf)
    tot = length(inf_inds)

    inf_times_infs = res[:,2][inf_inds]
    rem_times_infs = res[:,3][inf_inds]

    ev_infs = cat(inf_times_infs, rem_times_infs, dims = 2)
    ev_infs_sort = ev_infs[sortperm(ev_infs[:, 1]), :]

    plot(ev_infs_sort[1, :], [1, 1])
    for i in 2:(tot-1)
        plot!(ev_infs_sort[i, :], [i, i])
    end
    display(plot!(ev_infs_sort[tot, :], [tot, tot], legend = false))
end

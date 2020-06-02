
## Packages ##

using Distributions
using Distances
using LinearAlgebra
using InvertedIndices
using Random
using Plots
gr()

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

## Beta Matrix functio ##

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

## Seed finder function

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

@time begin
    seedfinder(100, [0.002 0.001], 10, 25, 1, 1000)
end

#N, betas, ds, tinf, Δ, seeds = 500, [0.0004 0.0002], 12, 125, 10, 250

# 60-80s first run
#35s second run
# different output...


### Example epidemic ###

Random.seed!(68)
exN = 100
exdistmat = unifdistmat(exN, 20, 20)[2]
exbetamat = betamatform(exdistmat, [0.0004 0.0002], 12)
exresults = GSEsim(exN, exbetamat, 0.15)

totinf = Array{Any}(undef, 500)
for i in 1:exN
    global totinf[i] = length(findall(exresults[:, 2] .< exresults[i, 2]))
end







### Bonus ###

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

@time begin
    seedfinder(100, [0.002 0.001], 10, 25, 1, 1000, 99)
end

using Distributions
using Distances
using LinearAlgebra
using InvertedIndices

## Distance Matrix ##

function unifdistmat(N, xlim, ylim)
    xcoords = xlim * rand(N)
    ycoords = ylim * rand(N)

    xcoords = reshape(xcoords, :, 1)
    ycoords = reshape(ycoords, :, 1)
    xycoords = [xcoords  ycoords]

    distmat = pairwise(Euclidean(), xycoords, xycoords, dims = 1)

    return [xycoords, distmat]
end

testxy, testdistmat = unifdistmat(100, 20, 20)

## Beta Matrix ##

function betamatform(distmat, betas, d)

    betamat = Array{Float64}(undef, size(distmat))

    index = findall(distmat .< d)
    betamat[index] .= betas[1]

    index = findall(distmat .>= d)
    betamat[index] .= betas[2]

    betamat[diagind(betamat)] .= 0

    return betamat
end

testbetamat = betamatform(testdistmat, [1 2], [10])



## GSE function ##

### Initialise the "Time until contact" matrix ###

function generatetuc(betamat)
    popn = size(betamat, 1)
    W = Array{Float64}(undef, popn, popn)
    uniqueness = popn*popn-popn+1

    while size(unique(W)) != Tuple{Int64}(uniqueness)
        println("try again")
        for i = 1:popn
            for j = (1:popn)[Not(i)]
                global W[i,j] = rand(Exponential(1/(betamat[i,j])))
            end
        end
        W[diagind(W)] .= Inf

        println("The matrix has ", "$(size(unique(W)))", " unique elements")
        println("It should have ", "$uniqueness", " unique elements")
    end

    return W
end

Wmat = generatetuc(testbetamat)

### Initialise the "Infectious period" matrix ###

Q = rand(Exponential(gamma), N)

Qmat = reshape(repeat(Q, outer = N), N, N)

### Initialise the "infection times" vector ###

inftimesfull = [0 ; fill(Inf, (N-1))]

### Initialise functional objects ###

pop = 1:N
infect = [1, 2]
sus = pop[Not(infect)]

## The simulation ##

### Calculate the current submatrices that are relevant ###

Wcur = Wmat[infect, sus]
Qcur = Qmat[infect, sus]

### Calculate which individuals can become infected ###

within = Wcur .< Qcur
whichinf = findall(within .== true)
whoinf = reshape(reinterpret(Int64, whichinf), length(whichinf), 2)
track = [ infect[whoinf[:,1]]  sus[whoinf[:,2]] ]

### Calculate the infection times of all potential infectees

inftimes = inftimesfull[track[:,1]] + Wcur[whichinf]

### Recover the index of the newly infected individual ###

if length(inftimes != 0)
    newinftime, up = findmin(inftimes)
    newinf = track[up, :][2]

### Update population ###

    inftimesfull[newinf] = newinftime

    infect = [infect ; newinf]
    sus = sus[Not(infect)]

end


#################
####################
#######################
####################
##################

##params

N=100
gamma = 0.15
betamat = betamat1
Random.seed!(2)
## function

Wmat = generatetuc(betamat)

### Initialise the "Infectious period" matrix ###

Q = rand(Exponential(gamma), N)

Qmat = reshape(repeat(Q, outer = N), N, N)

### Initialise the "infection times" vector ###

inftimesfull = [0 ; fill(Inf, (N-1))]

### Initialise functional objects ###

pop = 1:N
infect = [1, 2]
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

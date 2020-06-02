using Distributions
using InvertedIndices
using Random
using Plots
using DataFrames

using Distances
using LinearAlgebra

# Functions needed:
    # betamatform
    # unifdistmat
    # GSEsim (1&2)
    # plot_inf_periods

## Example Epidemic 1 ##
# Simple dummy epidemic, hand made

function Create_epidemic_1()
    inf_times_1 = [0, 3, Inf, 5, 7, 9, Inf, Inf, Inf]
    rem_times_1 = [11, 6, Inf, 8, 10, 12, Inf, Inf, Inf]
    inf_inds_1 = [true, true, false, true, true, true, false, false, false]

    inf_periods_1 = rem_times_1 - inf_times_1

    xcoords_1 = [1,2,3,1,2,3,1,2,3]
    ycoords_1 = [3,3,3,2,2,2,1,1,1]
    xycoords_1 = [xcoords_1  ycoords_1]

    distmat_1 = pairwise(Euclidean(), xycoords_1, xycoords_1, dims = 1)
    betamat_1 = betamatform(distmat_1, [0.002, 0.001], 1.5)

    res_1 = [1:9 inf_times_1 rem_times_1 inf_periods_1]

    return(res_1, xycoords_1, distmat_1, betamat_1)
end

res_ex1, xycords_ex1, distmat_ex1, betamat_ex1 = Create_epidemic_1()

plot_inf_periods(res_ex1)


## Example Epidemic 2 ##
# Different simple dummy epidemic, initial infected not first individual, hand made

function Create_epidemic_2()
    inf_times_2 = [Inf, 3, Inf, 5, 13, 9, Inf, 1, Inf]
    rem_times_2 = [Inf, 16, Inf, 6, 18, 12, Inf, 15, Inf]
    inf_inds_2 = [false, true, false, true, true, true, true, false, false]

    inf_periods_2 = rem_times_2 - inf_times_2

    xcoords_2 = [1,2,3,1,2,3,1,2,3]
    ycoords_2 = [3,3,3,2,2,2,1,1,1]
    xycoords_2 = [xcoords_2  ycoords_2]

    distmat_2 = pairwise(Euclidean(), xycoords_2, xycoords_2, dims = 1)
    betamat_2 = betamatform(distmat_2, [0.002, 0.001], 1.5)

    res_2 = [1:9 inf_times_2 rem_times_2 inf_periods_2]

    return(res_2, xycoords_2, distmat_2, betamat_2)
end

res_ex2, xycords_ex2, distmat_ex2, betamat_ex2 = Create_epidemic_2()

plot_inf_periods(res_ex2)


## Example Epidemic 3 ##
# Real simulated epidemic, initial individual first individual

function Create_epidemic_3()

    Random.seed!(403)

    xycoords_3, distmat_3 = unifdistmat(100, 20, 20)
    betamat_3 = betamatform(distmat_3, [0.002, 0.001], 10)

    res_3 = GSEsim(100, betamat_3, 0.15)

    return(res_3, xycoords_3, distmat_3, betamat_3)
end

res_ex3, xycords_ex3, distmat_ex3, betamat_ex3 = Create_epidemic_3()

plot_inf_periods(res_ex3)


## Example Epidemic 3 ##
# Real simulated epidemic, initial individual 99th individual

function Create_epidemic_4()

    Random.seed!(89)

    xycoords_4, distmat_4 = unifdistmat(100, 20, 20)
    betamat_4 = betamatform(distmat_4, [0.002, 0.001], 10)

    res_4 = GSEsim(100, betamat_4, 0.15, 99)

    return(res_4, xycoords_4, distmat_4, betamat_4)
end

res_ex4, xycords_ex4, distmat_ex4, betamat_ex4 = Create_epidemic_4()

plot_inf_periods(res_ex4)



### EXPORT EPIDEMICS ###
using DelimitedFiles

# Epidemic 1
DelimitedFiles.writedlm("res_ex1.csv", res_ex1, ',')
DelimitedFiles.writedlm("xycords_ex1.csv", xycords_ex1, ',')
DelimitedFiles.writedlm("distmat_ex1.csv", distmat_ex1, ',')
DelimitedFiles.writedlm("betamat_ex1.csv", betamat_ex1, ',')

# Epidemic 2
DelimitedFiles.writedlm("res_ex2.csv", res_ex2, ',')
DelimitedFiles.writedlm("xycords_ex2.csv", xycords_ex2, ',')
DelimitedFiles.writedlm("distmat_ex2.csv", distmat_ex2, ',')
DelimitedFiles.writedlm("betamat_ex2.csv", betamat_ex2, ',')

# Epidemic 3
DelimitedFiles.writedlm("res_ex3.csv", res_ex3, ',')
DelimitedFiles.writedlm("xycords_ex3.csv", xycords_ex3, ',')
DelimitedFiles.writedlm("distmat_ex3.csv", distmat_ex3, ',')
DelimitedFiles.writedlm("betamat_ex3.csv", betamat_ex3, ',')

# Epidemic 4
DelimitedFiles.writedlm("res_ex4.csv", res_ex4, ',')
DelimitedFiles.writedlm("xycords_ex4.csv", xycords_ex4, ',')
DelimitedFiles.writedlm("distmat_ex4.csv", distmat_ex4, ',')
DelimitedFiles.writedlm("betamat_ex4.csv", betamat_ex4, ',')

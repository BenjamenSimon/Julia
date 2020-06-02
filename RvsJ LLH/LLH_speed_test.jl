using Distributions
using Distances
using LinearAlgebra
using InvertedIndices
using Random
using BenchmarkTools



#Test function: Example 1
@benchmark simple_log_likelihood(res_ex1[:,2], res_ex1[:,3], betamat_ex1)

#Test function: Example 2
@benchmark simple_log_likelihood(res_ex2[:,2], res_ex2[:,3], betamat_ex2)

#Test function: Example 3
@benchmark simple_log_likelihood(res_ex3[:,2], res_ex3[:,3], betamat_ex3)

#Test function: Example 4
@benchmark simple_log_likelihood(res_ex4[:,2], res_ex4[:,3], betamat_ex4)


#Test function: Example 1
@benchmark sophisticated_log_likelihood(res_ex1[:,2], res_ex1[:,3], betamat_ex1)

#Test function: Example 2
@benchmark sophisticated_log_likelihood(res_ex2[:,2], res_ex2[:,3], betamat_ex2)

#Test function: Example 3
@benchmark sophisticated_log_likelihood(res_ex3[:,2], res_ex3[:,3], betamat_ex3)

#Test function: Example 4
@benchmark sophisticated_log_likelihood(res_ex4[:,2], res_ex4[:,3], betamat_ex4)

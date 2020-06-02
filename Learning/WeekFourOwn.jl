using Distributions
using StatsBase
using CSV
using DataFrames
using HypothesisTests
using StatsPlots
using GLM


age = rand(18:80, 100) # Uniform Distribuion
wcc = rand(Distributions.Normal(12, 2), 100) #Normal
crp = rand(Distributions.Chisq(4), 100) #Chi-Sq
treatment = rand(["A", "B"], 100) #uniformly weighted
result = rand(["Improved", "Static", "Worse"], 100)


StatsBase.describe(age)

StatsBase.summarystats(age)

data = DataFrame(Age = age, WCC = wcc, CRP = crp, Trt = treatment, Result = result)

size(data)

first(data, 6)

dataA = data[data[:Trt] .== "A", :] # only patient in trt group A
dataB = data[ data[:Trt] .== "B", :] # only trt B

describe(data)

#####

by(data, :Trt, df -> DataFrame(N = size(df, 1)))

by(data, :Trt, size)

by(data, :Trt, df -> mean(df.Age))

by(data, :Trt, df -> describe(df.Age))

@df data density(:Age, group = :Trt,
                  title = "nadhd",
                  xlab = "ndbsa",
                  ylab = "jbaad",
                  legend = :topright)

@df data density(:Age, group = :Result,
                  title = "nadhd",
                  xlab = "ndbsa",
                  ylab = "jbaad",
                  legend = :topright)

@df data density(:Age, group = (:Trt, :Result),
                title = "nadhd",
                xlab = "ndbsa",
                ylab = "jbaad",
                legend = :topright)

@df data boxplot(:Trt, :WCC, lab = "WCC",
                title = "nadhd",
                xlab = "ndbsa",
                ylab = "jbaad",
                legend = :topright)

@df data boxplot(:Result, :WCC, lab = "WCC",
                title = "nadhd",
                xlab = "ndbsa",
                ylab = "jbaad",

@df data corrplot([:Age :WCC :CRP], grid = false)

@df data cornerplot([:Age :WCC :CRP], grid = false, compact = true)

####

# ## Inferential statistics
# We will begin by using Student's _t_ test to compare the mean of a numerical variable between two groups.

# Difference in age between patients in groups A and B
HypothesisTests.EqualVarianceTTest(dataA[:, :Age], dataB[:Age])

# Only the p value for the difference in white cell count between patients in groups A and B
pvalue(EqualVarianceTTest(dataA[:WCC], dataB[:WCC]))

# Difference in c-reactive protein level between patients in groups A and B for unequal variances
UnequalVarianceTTest(dataA[:CRP], dataB[:CRP])

# We can create a variety of linear models using the `GLM.fit()` function.

# Simple model predicting CRP
fit(LinearModel, @formula(CRP ~ 1), data)

# Adding Age as a predictor variable
fit(LinearModel, @formula(CRP ~ Age), data)

# Adding Age and WCC as predictor variables
fit(LinearModel, @formula(CRP ~ Age + WCC), data)

# We can conduct a $\chi^2$ test for independence using the `HypothesisTests.ChisqTest()` function.  First we need to look at the counts.  Below we calculate the number of unique values for the Result variable sample space for patients in groups A and B.

by(dataA, :Result, df -> DataFrame(N = size(df, 1)))

by(dataB, :Result, df -> DataFrame(N = size(df, 1)))

# Enter the data in similar order here
observed = reshape([22, 17, 18, 18, 11, 14], (2, 3))
observed

ChisqTest(observed)

####


#mat = Array{Float64}(undef, 3, 3)
mat = [1. 2. 3.; 4. 5. 6.; 7. 8. 9.]
W = Array{Float64}(undef, 3, 3)

for i = 1:size(mat,1)
    for j = 1:size(mat,2)
        global W[i,j] = rand(Exponential(mat[i,j]))
    end
end

W2 = generatetuc(mat)

mat2 = rand(1:2, 100, 100)

W3 = generatetuc(mat2)

for i = 1:100
    W3 = generatetuc(mat2)
    sW3 = size(unique(W3))
    if  sW3 != Tuple{Int64}(100^2)
        println(i)
        println(sW3)
    end
end

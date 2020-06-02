b = 4
f(x) = b*x
f(6)
#don't need to define b inside the function

function updateSIR(popn)
    sus = popn[1];
    inf = popn[2]
    rem = popn[3]
    newS = sus - lambda*sus*inf*dt
    newI = inf + lambda*sus*inf*dt - gam*inf*dt #gamma is a built in function so don't want to overwrite
    newR = rem + gam*inf*dt
    return [newS newI newR] #no commas to make this a one row of a two dimensional array
end

#set params
dt = 0.5
lambda = 1/200
gam = 1/10

#specify input vector
s, i, r = 1000., 10., 20.
vec = [s i r]
updateSIR(vec)

#Full
lambda = 1/20000   # infection rate parameter (assumes rates are per day)
gam = 1/10       # recovery rate parameter  (ditto)
dt = 0.5         # length of time step in days
tfinal = 610;    # respecting community values: lowercase only in the names
s0 = 10000.0     # initial susceptibles, note that we use the  type Float64 from the start
i0 = 4.          # initial infecteds; set this to 1. to  mimic an epidemic with an index case
r0 = 0.          # not always the case, of course

# initialise the current run
nsteps = round(Int64, tfinal/dt)    # note the use of round() with type Int64 to ensure that nsteps is an integer
resultvals = Array{Float64}(undef, nsteps+1, 3)  # values, rows, columns
timevec = Array{Float64}(undef, nsteps+1)        # values, rows
resultvals[1,:] = [s0, i0, r0]  # ... and assign them to the first row
timevec[1] = 0.                 # also Float64, of course.

# execute the current run
for step  = 1:nsteps
    resultvals[step+1, :] = updateSIR(resultvals[step, :])  # NB! pay careful attention to the rows being used
    timevec[step+1] = timevec[step] + dt
end

using Plots
gr()

plot(timevec, resultvals)

plot(timevec, resultvals,  # we should of course at a minimum provide some labels
title  = "Example of SIR results",
xlabel = "Epidemic day",
ylabel = "Population size",
label  = ["Susceptibles" "Infecteds" "Removeds"]
)

######

tempvar = Array{Any}(undef, 4)
fill!(tempvar, "hello")

#only one active plot can exist

function approxcos(x)
    #initialise the output ... note the use of size() to specify the dimensions of the output vector
    outval = Array{Any}(undef,size(x))

    # now we loop over the input vector, and for each  element calculate and store the approximation
    ii = 0  # this will be the index into the vector
    for aa in x   # this aa is just a number, an element of the vector
        y = 1 - aa^2/2 + aa^4/24 - aa^6/720 + aa^8/(56*720) # the approximation ...
        ii = ii+1            #this sets the index correctly
        outval[ii] = y     # and this stores the approximation in the right place

    end

    return outval
end

x1 = 4*rand(10)  # rand() is one of several random number functions in Julia.
#                  It returns numbers that uniformly fill the interval [0, 1)
#                   .... here we use it get a set of sampling points in the interval [0, 4)

x2 = range(0., stop=4., step=0.01)   # look up range() using "?" ... it's a nice way to get evenly spaced points


y1 = approxcos(x1)
y2 = cos.(x2);

# now the plots
using Plots; gr()    # it is sometimes convenient to cram a line in this way

#first the plot of the approximation points
scatter(x1, y1, legend=:false, title="Illustrating 6th-order approximation to cos")

plot!(x2, y2)   #then add the accurate line with plot!()

using Plots

f(x) = 5 * x^2 + 3 * x - 20
plot(f, -10:10)
plot!(zero, -10, 10)

f(x) = -9*x^2 + 3*x + 6
plot(f, -2, 2)
plot(f, -2:2)
plot(f)

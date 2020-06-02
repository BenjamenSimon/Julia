using DelimitedFiles
# load the package used for importing csvs

wikiEVDraw = DelimitedFiles.readdlm("wikipediaEVDraw.csv", ',')
#getting right quotes is important

using Dates

Dates.DateTime(wikiEVDraw[1,1], "d u y")
# the " " tell Julia what format the date is in
#= d means day
    u means month abbreviated
    y means year
=#

# Frequently a loop will run faster than vectorised code in Julia

for num = 3:7
    println("num is now $num")
end

testvalues = [23, "my name is not a name", 'α']
for x in testvalues
    println("The value of x is now $x")
end
# a loop can iterate over an array


##########

col1 = wikiEVDraw[:, 1] # colon means alll rows

for i = 1:length(col1)
    col1[i] = Dates.DateTime(col1[i], "d u y")
end

Dates.datetime2rata(col1[1])

dayssincemar22(x) = Dates.datetime2rata(x) - Dates.datetime2rata(col1[54])
epidays = Array{Int64}(undef, 54)
for i = 1:length(col1)
    epidays[i] = dayssincemar22(col1[i])
end

wikiEVDraw[:, 1] = epidays
DelimitedFiles.writedlm("wikipediaEVDdatesconverted2.csv", wikiEVDraw, ',')
# note the delimiter, Julias default is a tab, to get .csv must specify comma


mylist = [3, 2, 1]
count = 1
for i = mylist
  mylist[i] = count
  global count=count+1
end

mylist

count=0
for i=1:3
  for j=1:3
    global count=count+1
  end
end
count

summedvals = 3
for k = 1:2:5
    global summedvals = summedvals + k
end
summedvals

for k = 1:2:5
    println(k)
end
# k = (1,3,5)
# 1:5 by 2

function test(x)
    if x==1
        return 1
    end
    if x<=0
        return 0
    end
    return test(x-1)+test(x-2)
end
test(3)

############

using DelimitedFiles
EVDdata = DelimitedFiles.readdlm("wikipediaEVDdatesconverted.csv", ',')
epidays = EVDdata[:,1]
allcases = EVDdata[:, 2]

using Pkg
Pkg.add("Plots")
using Plots
# lots of 'backends' for plots
# used to use pyplot backend but now use gr
gr()

plot(epidays, allcases)

plot(epidays, allcases, linetype = :scatter, marker = :diamond)

plot(epidays, allcases,
title       = "West",
xlabel      = "Days since",
ylabel      = "Total cases",
marker      = (:diamond, 8, "orange"),
line        = (:path, "gray"),
legend      = false,
grid        = false)

savefig("plot1") #default png
savefig("plot2.pdf")
savefig("plot3.png")

##########

EVDdata[end-9:end, :]
# last row - 9 to the last row

a = rand() #a random value
println("a now has vlaue $a")
if a > 0.5
    println("this is quite a large value")
end

for k in 1:8
    b = rand()
    println("b now has the value $b")
    if b>0.5
        println("this is quite a large value")
    end
end

rows, cols = size(EVDdata)
for j = 1:cols
    for i = 1:rows
        if !isdigit(string(EVDdata[i, j])[1]) #! means not
            EVDdata[i,j] = 0
        end
    end
end
# This is a bit of a hacky way to check for non-digits and then
# change them to be 0

#extract the data
epidays = EVDdata[:, 1]
EVDcasesbycountry = EVDdata[:, [4,6,8]]
typeof(EVDcasesbycountry)
EVDcasesbycountry = convert(Array{Float64, 2}, EVDcasesbycountry)
typeof(EVDcasesbycountry)
#new_array = Array{Float64}(array)
#Float64.(array)

#load plots and plot them
using Plots
gr()
plot(epidays, EVDcasesbycountry)
# doesn't like matrix of type any

plot(epidays, EVDcasesbycountry,
marker      = ([:octagon :star7 :square], 4), #no commas
label       = ["Guinea" "liberia" "Sierra Leone"], #no commas
title       = "title",
xlabel      = "x label",
ylabel      = "y label",
line        = (:scatter),
legend      = :topleft
)

savefig("plot4.pdf")

f(x) = 3 * x^2 + 6 * x - 9
for x = -5:5
  println("(",x, ", ", f(x), ")")
end
using Plots
gr() # Activate the GR backend for use with Plots
plot(f, -4, 3) # plot f over [-4,4]
plot!(zero, -4, 3) #plots a horizontal line at zero
# and alters the current plot with !

data = [1.6800483  -1.641695388;
        0.501309281 -0.977697538;
        1.528012113 0.52771122;
        1.70012253 1.711524991;
        1.992493625 1.891000015;
        2.706075824 -0.463427794;
        2.994931927 -0.443566619;
        3.491852811 -1.275179133;
        3.501191722 -0.690499597;
        4.459924502 -5.516130799;
        4.936965851 -6.001703074;
        5.023289852 -8.36416901;
        5.04233698 -7.924477517;
        5.50739285 -10.77482371;
        5.568665171 -10.9171878]

using Plots
gr() # Activate the GR backend for use with Plots
# Use the data array to assign values for x and y here
x,y = data[:,1], data[:,2]
plot(x, y, linetype = :scatter, leg = false)
# scatter(x, y) # this is an alternative method, but does make a legend

data[:,1][end - 3:end]

n = 20
x = sort(rand(20)); y = rand(20)
Plots.scatter(x, y)
plot!(x,y, title = "dnd")
plot!(x,y, leg = false, title = "nf")


1+1

1+1
2+1

println("Hello, world!")

5

5+5

println(5+2+3)

println("5+2+3")

println(5-3.8)
println("5+2+3 ", "equals ", 5+2+3)

println("Hello, world!" * " Oh no!")
println("Hello, world!" * " Oh no!"^5)

5+9 #comment

#?println #help

0.2 + 0.1 - 3 * 6.7 / 4 - 1 - 2 * 3
# Little number means mixing types
# 6.7 is Float64
# and 3 is Int64
# Mixing types majorly slows julia down

println( 3^2^3 )        # exponentiation right-to-left is the Julia convention
println( (3^2)^3 )      # forcing left-to-right using parentheses

2^2^1^1.3^1.5^1.7^20

!!!!true

# First expression
3 & 5 > 0

# Second expression
8 & 5 > 0

1+3<5||2+2<1

#α #\alpha tab

Array{Int64}(undef, 3)

Array{String}(undef, 3)

greeting = "Hello, world!"  # creates a variable called "greeting" whose value is a string

println(greeting)          # println uses the value of greeting when it prints the message~

人 = 20
生 = 11.111
[人, 生]          # another way to make a 1-dimensional array: comma-separated list inside brackets

abstypevariable = Array{Integer}(undef, 2,3) # A two-dimensional array with 2 rows and 3 columns

abstypevariable[2,1]    #NB --- note the brackets, that's how to access elements of an array

'a'==="a"

x = Array{Int64}(undef,11, 12)

typeof(x)

a, b, c = cos(0.2), log(10), abs(-1.22)

a

#= multi
   line
   comment
=#

#= a function with an exclamation at the end is a
   julia convention that says the function will
   alter the input
   such as delete!()
=#


? muladd
#won't work in script
#but does work in console

methods(muladd)

myfunc(firstvar) = 20*firstvar

myfunc(2.5)

addxtoy(x,y) = x + y #don't have to say its going to be a func

addxtoy(2, 2.5) #mixed types, slow, but it works

function nextfunc(a, b, c)
    a*b + c

    #white space, additional comments, all valid
end

nextfunc(7, 5, 3)

nextfunc(7.0, 5, 3)

function showdebugprintln(testvar)
    println("inside the showdebugprint() now")   #this line announces where the report is coming from
    println("The type of testvar is $(typeof(testvar)) and the value of testvar is $testvar")
    #                  and this line reports what value, and hence what type, testvar actually has here
end

# $ is an escape character

a = ['1', 2.]
showdebugprintln(a)
# should be any since mixing characters and float64

mycos(x) = cos(x)

mycos(adj, hyp) = adj/hyp

mycos(12, 13)

methods(mycos)
# tells you the line where it is defined

mycos(theta::Float64) = cos(theta)

#to clear anything have to clear everything
#go to kernel and restart and clear output

#So now can restart
mycos(theta::Float64) = cos(theta)
mycos(hyp, adj) = adj/hyp

mycos(1)
#Doesn't work now because there is no method for Int64

a,b=1,2

a,b=b,a

a

b

function test(input)
  println("$input"^2)
end

test(3)

function test(input)
  println("input"^2)
end

test(3)

add2(x,y) = return x+y
println(add2(5,7))

function add2(x,y) = x+y

function coordinates(x, y=0, z=0)
  println("($x, $y, $z)")
end

methods(coordinates)

coordinates(1,2)

coordinates(1,0,2)

coordinates(1,0)

coordinates(1,0,0)

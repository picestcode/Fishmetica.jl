using Fishmetica
using Base.Test

tic()

println("Test 1")
@time  include("test1.jl")

println("Test 2")
@time include("test2.jl")

println("Test 3")
@time include("test3.jl")

toc()

# write your own tests here
#@test 1 == 2


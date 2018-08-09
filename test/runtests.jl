using Fishmetica
using Base.Test

tic()

println("Test 1")
@time  include("/home/igor/.julia/v0.6/Fishmetica/test/test1.jl")

println("Test 2")
@time include("/home/igor/.julia/v0.6/Fishmetica/test/test2.jl")

println("Test 3")
@time include("/home/igor/.julia/v0.6/Fishmetica/test/test3.jl")

toc()

# write your own tests here
#@test 1 == 2


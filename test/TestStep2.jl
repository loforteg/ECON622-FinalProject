# Perform some tests on Step2

using Test
include("../src/Step2.jl")

@testset "Step2.jl" begin

    # Setting
    eH = 1.0
    eL = 0.1
    E = [eL; eH]
    e_size = size(E)[1]
    # Transition matrix: Pr(j|i) = Π(i,j) (columns sum to 1)
    πHH = 0.925 # e_H | e_H
    πHL = 0.5   # e_H | e_L
    Π = [1 - πHL  πHL
        1 - πHH πHH]
    q = 1.3
    a_lb = -2.0
    a_max = 5
    a_size = 100
    A = range(a_lb, a_max, length = a_size)
    A = [A; ];

    # Test 1: check first dimension of policy function is respected
    pol_func = ones(a_size - 10, e_size)
    @test_throws AssertionError step2(A, E, Π, pol_func)

    # Test 2: check second dimension of policy function is respected
    pol_func = ones(a_size, e_size - 1)
    @test_throws AssertionError step2(A, E, Π, pol_func)

end

# Perform some tests on Step1EndGrid

using Test
include("../src/Step1EndGrid.jl")

@testset "Step1EndGrid.jl" begin

    # Setting
    β = 0.9932
    eH = 1.0
    eL = 0.1
    E = [eL; eH]
    # Transition matrix: Pr(j|i) = Π(i,j) (columns sum to 1)
    πHH = 0.925 # e_H | e_H
    πHL = 0.5   # e_H | e_L
    Π = [1 - πHL  πHL
        1 - πHH πHH]
    q = 1.3
    σ = 1.5
    a_lb = -2.0
    a_max = 5
    a_size = 100
    A = range(a_lb, a_max, length = a_size)
    A = [A; ]
    Aprime = repeat(A, 1, e_size)

    # Test 1: check β makes sense
    β = 1.1
    @test_throws AssertionError EndGridSearch(A, E, Π, β, σ, q)

    # Test 2: check transition matrix
    β = 0.9932
    Π = [0.3 0.5
        0.5 0.5]
    @test_throws AssertionError EndGridSearch(A, E, Π, β, σ, q)

    # Test 3: check that the price makes sense
    q = 0.8
    @test_throws AssertionError EndGridSearch(A, E, Π, β, σ, q)

    # Test 4: check that σ is appropriate
    σ = 0.9
    @test_throws AssertionError EndGridSearch(A, E, Π, β, σ, q)

end

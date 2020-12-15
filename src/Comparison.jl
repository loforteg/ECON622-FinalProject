## Comparison of exogenous grid search and endogenous grid search

using BenchmarkTools

# Set parameters
β = 0.9932
σ = 1.5


# Exogenous shocks
eH = 1.0
eL = 0.1
E = [eL; eH]

# Transition matrix: Pr(j|i) = Π(i,j) (columns sum to 1)
πHH = 0.925 # e_H | e_H
πHL = 0.5   # e_H | e_L
Π = [1 - πHL  πHL
    1 - πHH πHH]

# Assets
a_lb = -2.0
a_max = 5
a_size = 350
A = range(a_lb, a_max, length = a_size)
A = [A; ];

# Guess initial price
q = 1.1


# Check performance simple grid search method
@btime step1(A, E, Π, q, β, σ)

# Check performance endogenous grid search method
@btime EndGridSearch(A, E, Π, β, σ, q)



include("Step1Grid.jl")
include("Step1EndGrid.jl")

## Step3 is the iteration to find the price q

module Step3

using LinearAlgebra, Random, Distributions, Statistics, Plots
using BenchmarkTools, Interpolations

"""
    Relies on step2 and finds the optimal price.
"""

# Define function that computes the market for savings
function mktclearing(A, E, pol_func, F_update)

    """
    Function that computes market clearing conditions.
    Inputs:
        - A: set of assets
        - E: set of shocks
        - pol_func: optimal policy function
        - F_update: stationary distribution function Ψ
    Output:
        - market: demand for assets
    """

    # retrieve dimensions
    a_size = length(A)
    e_size = size(E)[1]

    # initialize market for update
    market = 0.0

    for a = 2:a_size
        for e = 1:e_size

            market += pol_func[a, e] * (F_update[a, e] - F_update[a-1, e])

        end
    end

    return market

end


# Define function for Step3
function step3(A, E, Π, q0, β, σ; maxT = 100, tol_iterations = 0.01,
                tol_market = 0.0025, tol_q = 0.001, weight = 0.5)

    """
    Function that finds the optimal price q.
    Relies on market in the update procedure.
    Inputs:
        - A: set of assets
        - E: set of shocks
        - q0: guess price
        - β: discount factor
        - σ: preference parameter
    Output:
        - q: optimal price
        - policy: optimal policy
        - Ψ: stationary distribution Ψ
    """

    # check q > β
    @assert q0 > β

    q = q0
    solStep1 = step1(A, E, Π, q, β, σ; tol = tol_iterations)
    solStep2 = step2(A, E, Π, solStep1; tol = tol_iterations)

    # Start price iteration
    normdiff = Inf
    t = 1

    while normdiff > tol_q && t <= maxT

        # solve step1
        solStep1 = step1(A, E, Π, q, β, σ; tol = tol_iterations)

        # solve step2
        solStep2 = step2(A, E, Π, solStep1; tol = tol_iterations)

        # compute market
        market = mktclearing(A, E, solStep1, solStep2)


        # update q based on market clearing:
        # increase guess if market > 0 (but just by a little bit)
        if market > 0 && abs(market) > tol_market
            q = q + (1 - weight) * tol_market

        # decrease guess if market < 0 (but just by a little bit)
        elseif market < 0 && abs(market) > tol_market
            q = q - (1 - weight) * tol_market

        # if abs(market) < tol_market, then no need to update q
        else
            q = q
        end

        t = t + 1
        normdiff = abs(q - q0)

        q0 = q

    end

    return (q = q, policy = solStep1, Ψ = solStep2)

end


include("Step2.jl")
include("Step1.jl")
include("Step1EndGrid.jl")

end

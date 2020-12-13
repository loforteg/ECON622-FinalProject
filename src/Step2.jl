## Step2 is iteration to find the stationary distribution

module Step2

using LinearAlgebra, Random, Distributions, Statistics, Plots
using BenchmarkTools, Interpolations

"""
    Relies on step1 (either with simple grid search or endogenous) and finds
    stationary distribution.
"""

# Define function for inverse of policy function a^{-1}
function pol_inv(A, pol_func, aprime, eprime)

    """
    Function that finds the inverse of the optimal policy function.
    Inputs:
        - A: set of assets
        - pol_func: optimal policy function
        - aprime: assets in the next period
        - eprime: exogenous state variable in the next period
    Output: inverse of the policy function
    """

    # if a' < minimum value of policy function, then a must be a_lowerbar
    if minimum(pol_func[:, eprime]) > A[aprime]

        return 1

    # when e'=2, distinguish pre vs post cross 45° line (endogenous a_UpperBar)
    elseif eprime == 2

        return min(findlast(pol_func[:, eprime] .<= A[aprime]),
            findfirst((pol_func[:, eprime] - A) .< 0.01))

    # any other case, find element of A s.t. policy function is closest to a'
    else

        return findlast(pol_func[:, eprime] .<= A[aprime])

    end

end


# Define function for Step2
function step2(A, E, Π, pol_func; tol = 0.01, maxT = 600)

    """
    Function that finds the stationary distribution Ψ.
    Relies on pol_inv in the code.
    Inputs:
        - A: set of assets
        - E: set of exogenous shocks
        - Π: transition matrix
        - pol_func: optimal policy function
    Output:
        - F_update: stationary distribution Ψ
    """

    # retrieve dimensions
    a_size = length(A)
    e_size = size(E)[1]

    # Find endogenous position where a' = a
    aBar = findlast((pol_func[:,2] - A) .< tol)
    cut = max(aBar, 2)           # if aBar = 1 then range will not work

    aux1 = range(zero(eltype(A)), one(eltype(A))/e_size, length = cut)
    aux2 = ones(eltype(A), a_size - cut)/e_size
    aux = [aux1; aux2]           # 350x1
    F0 = repeat(aux, 1, e_size)  # 350x2

    F_update = copy(F0)

    # Iteration of distribution
    normdiff = Inf
    t = 1

    while normdiff > tol && t <= maxT

        for a = 1:a_size
            for e = 1:e_size

                F_update[a,e] = Π[1, e] * F0[pol_inv(A, pol_func, a, 1), 1] +
                          Π[2, e] * F0[pol_inv(A, pol_func, a, 2), 2]

            end
        end

        t += 1
        normdiff = maximum(abs.(F_update - F0))

        F0 = copy(F_update)

    end

    return F_update

end


include("Step1Grid.jl")




end

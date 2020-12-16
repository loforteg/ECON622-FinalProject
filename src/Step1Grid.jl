## Step1Grid is a simple grid search to find optimal policy function

#module Step1Grid

using LinearAlgebra, Random, Distributions, Statistics, Plots
using BenchmarkTools, Interpolations

"""
    Performs a grid search to find the optimal policy function.
    Takes price q as given.
"""

# Define the function for Step1: iteration to find the optimal policy function

function step1(A, E, Π, q, β, σ; v_guess = 0.0, maxT = 600, tol = 0.01)

    """
    Function that finds the optimal policy function for a given price q.
    Inputs:
        - A: set of possible assets (endog. state var.)
        - E: set of exogenous shocks (exog. state var.)
        - Π: transition matrix
        - q: price
        - β: discount factor
        - σ: preferences parameter
    Output:
        - pol_func: optimal policy function
    """

    # Check that parameters and price make sense
    @assert β < 1
    @assert σ > 1
    @assert q > β

    a_size = length(A)
    e_size = size(E)[1]

    # Check that transition matrix is well defined
    for e = 1:e_size
        @assert sum(Π[e,:]) == 1
    end

    v0 = v_guess .* ones(a_size, e_size)
    v_update = zeros(a_size, e_size)
    pol_func = zeros(a_size, e_size)

    # Define the utility function (in the future, put it outside)
    u(c) = c.^(1-σ) ./ (1 - σ)

    # Iteration of value function
    normdiff = Inf
    t = 1

    while normdiff > tol && t <= maxT

        # Do linear interpolation of the value function on knots (A,E)
        v_interpol = LinearInterpolation((A, E), v0)

        # Replace new value function with solution
        for a = 1:a_size
            for e = 1:e_size

                # Generate value function V
                # If c <= 0, set value function to -∞
                V(x) = ((E[e] + A[a] - q*x) <= 0) ? -Inf :
                   u(E[e] + A[a] - q*x) + β*(v_interpol(x, E[1])* Π[e,1] +
                        v_interpol(x, E[2])* Π[e,2])

                # Find the maximum
                aux = findmax(V.(A))      # (maximum, index)
                v_update[a,e] = aux[1]    # actual maximum value
                pol_func[a,e] = A[aux[2]] # policy function a'

            end
        end

        t += 1
        normdiff = maximum(abs.(v_update - v0))

        v0 = copy(v_update)

    end

    return pol_func

end

#export step1

#end

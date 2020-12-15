## Step1EndGrid is an endogenous grid search to find optimal policy function

module Step1EndGrid

using LinearAlgebra, Random, Distributions, Statistics, Plots
using BenchmarkTools, Interpolations, Roots

"""
    Performs an endogenous search to find the optimal policy function.
    Takes price q as given.
"""

function RHSEuler(A, E, Π, q, β, σ, a, aprime, Aprimetol)

    # a is the actual value of assets
    # aprime is the position of the new assets

    """
    This function computes the RHS of the Euler Equation.
    Inputs:
        - A: grid of assets today
        - E: exogenous shocks
        - Π: transtion matrix
        - q: price (guessed one)
        - β: discount factor
        - σ: preference parameter
        - a: actual value of assets today
        - aprime: position of assets tomorrow
        - Aprimetol: provisional grid of future assets
    Output:
        - RHS: right hand side of the Euler Equation, for 2 exogenous state var.
    """

    a_size = length(A)
    e_size = size(E)[1]

    # Define the first derivative of utility function
    # (in the future, make it more general)
    up(c) = c.^(-σ)

    Uprime = zeros(e_size, 1)
    RHS = zeros(e_size, 1)

    for e = 1:e_size

        adoubleprime = LinearInterpolation(A[:, e], Aprimetol[:, e])[aprime]
        Uprime[e] = up(A[aprime, e] + E[e] - q * adoubleprime)
        RHS[e] = up(a + E[e] - A[aprime, e] * q) - β * dot(Π[e, :], Uprime)

    end

end



function initialpolicy(A, E, Π, β, σ, q; maxT = 600, tol = 0.01)

    """
    This function creates the first grid of assets in the following period.
    Inputs:
        - A: grid of assets today
        - E: exogenous shocks
        - Π: transtion matrix
        - β: discount factor
        - σ: preference parameter
        - q: price (guessed one)
    Output:
        - Aprime: grid of assets tomorrow
    """

    a_size = length(A)
    e_size = size(E)[1]

    Aprime = copy(A)
    Aprimetol = copy(A) .+ 2*tol

    # Iteration
    normdiff = Inf
    t = 1

    while normdiff > tol && t <= maxT

        for a = 1:a_size
            for e = 1:e_size

                FOCaprime(aprime) = RHSEuler(A, E, Π, q, β, σ, A[a, e], aprime, Aprimetol)

                if FOCaprime(A[1, e]) > 0.0
                    Aprime[a, e] = A[1, e]
                elseif FOCaprime(A[end, e]) < 0.0
                    Aprime[a, e] = A[end, e]
                else
                    Aprime[a, e] = fzero(FOCaprime, A[1, e], A[end, e])
                end

            end
        end

        t += 1
        normdiff = maximum(abs.(Aprimetol - Aprime))
        Aprimetol = copy(Aprime)

    end

    return Aprime

end


function EndGridSearch(Aprime, E, Π, β, q; maxT = 600, tol = 0.01)

    """
    Function that finds the optimal policy function.
    Inputs:
        - Aprime: grid of assets tomorrow
        - E: exogenous shocks
        - Π: transition matrix
        - β: discount factor
        - q: price (guess)
    Output:
        -
    """

    # Define the inverse of the first derivative of utility function
    # (in the future, make it more general)
    upinv(B) = B.^(-1/σ)

    # Check that discount factor and price make sense
    @assert β < 1
    @assert q > β

    a_size = length(A)
    e_size = size(E)[1]

    # Check that transition matrix is well defined
    for e = 1:e_size
        @assert sum(Π[e,:]) == 1
    end

    C = Aprime .+ E'
    Cupdate = 0.0 .* similar(Aprime)
    Aupdate = similar(Aprime)

    # Iteration
    normdiff = Inf
    t = 1

    while normdiff > tol && t <= maxT

        # Find Ctilda from the Euler Equation
        B = β .* up(C) * Π'
        Ctilda = upinv(B ./ q)

        # Solve for A from the budget constraint
        A = Ctilda + Aprime .* q .- E'

        # Adjust for whether budget constraint binds (only low shock)
        binds = Aprime .< A[:,1]'
        Cupdate[binds] .= (Aprime .+ E' .- Aprime[:,1]' .* q')[binds]

        # Finish this!!
        # Cupdate[.!binds] .= LinearInterpolation()[.!binds]

        t += 1
        normdiff = maximum(abs.(Cupdate - C))



    end

    # Finish this too!
    # Aprime = LinearInterpolation()

    return pol_func = Aprime

end


end

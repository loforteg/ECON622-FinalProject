## Step1EndGrid is an endogenous grid search to find the optimal policy function

#module Step1EndGrid

using LinearAlgebra, Random, Distributions, Statistics, Plots
using BenchmarkTools, Interpolations, Roots

"""
    Performs an endogenous search to find the optimal policy function.
    To do so, 3 functions are nested together:
    1. RHSEuler: computes the RHS of the Euler equation
    2. initialpolicy: for a guess value of consumption, finds a'
    3. EndGridSearch: solves for optimal consumption and a'
    Takes price q as given.
"""

function RHSEuler(A, E, Π, q, β, σ, e, avalue, aprime; tol = 0.01)

    # avalue is the actual value of assets
    # aprime is the actual value of the new assets

    """
    This function computes the RHS of the Euler Equation.
    Inputs:
        - A: assets today
        - E: exogenous shocks
        - Π: transtion matrix
        - q: price (guessed one)
        - β: discount factor
        - σ: preference parameter
        - e: position parameter for the shock (position, not value)
        - avalue: actual value of assets today (value, not position)
        - aprime: actual value of assets tomorrow
        - Aprimetol: provisional grid of future assets
    Output:
        - RHS: right hand side of the Euler Equation, for 2 exogenous state var.
    """

    a_size = length(A)
    e_size = size(E)[1]

    # Create a grid of A multiplied by the length of E
    Agrid = repeat(A, 1, e_size)

    Aprimetol = copy(Agrid) .+ 2*tol

    # Define the first derivative of utility function
    # (in the future, make it more general)
    lowerbound = 1E-7
    up(c::Real) = (max(c, lowerbound)).^(-σ)

    Uprime = zeros(e_size, 1)
    RHS = zeros(e_size, 1)

    aux = LinearInterpolation(Agrid[:, e], Aprimetol[:, e])
    adoubleprime = aux(aprime)
    Uprime = up(aprime .+ E[e] .- q .* adoubleprime)
    temp = Uprime .* ones(e_size, 1)
    temp2 = (Π[e, :]' * temp)[1,1]
    RHS = up(avalue + E[e] - aprime * q) - β * temp2

    return RHS

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

    # Create a grid of A multiplied by the length of E
    Agrid = repeat(A, 1, e_size)
    Aprime = copy(Agrid)
    Aprimetol = copy(Agrid) .+ 2*tol

    # Iteration
    normdiff = Inf
    t = 1

    while normdiff > tol && t <= maxT

        for a = 1:a_size
            for e = 1:e_size

                FOCfunction(aprime) = RHSEuler(A, E, Π, q, β, σ, e, Agrid[a, e], aprime)

                if FOCfunction(Agrid[1, e]) > 0.0
                    Aprime[a, e] = Agrid[1, e]
                elseif FOCfunction(Agrid[end, e]) < 0.0
                    Aprime[a, e] = Agrid[end, e]
                else
                    Aprime[a, e] = fzero(FOCfunction, Agrid[1, e], Agrid[end, e])
                end

            end
        end

        t += 1
        normdiff = maximum(abs.(Aprimetol - Aprime))
        Aprimetol = copy(Aprime)

    end

    return Aprime

end


function EndGridSearch(A, E, Π, β, σ, q; maxT = 600, tol = 0.01)

    """
    Function that finds the optimal policy function.
    Inputs:
        - A: assets today
        - E: exogenous shocks
        - Π: transition matrix
        - β: discount factor
        - σ: preference parameter
        - q: price (guess)
    Output:
        - pol_func: optimal policy function
    """

    # Check that parameters and price make sense
    @assert β < 1
    @assert σ > 1
    @assert q > β

    # Check that transition matrix is well defined
    for e = 1:e_size
        @assert sum(Π[e,:]) == 1
    end

    a_size = length(A)
    e_size = size(E)[1]

    # Define the first derivative of utility function and its inverse
    # (in the future, make it more general)
    lb = 1E-7
    lowerbound = lb .* ones(a_size, e_size)
    up(c) = max.(c, lowerbound).^(-σ)
    upinv(B) = B.^(-1/σ)

    # Compute Aprime: initial grid of next period's assets
    Aprime = initialpolicy(A, E, Π, β, σ, q)

    # Compute initial consumption (use it as guess for iteration)
    C = Aprime .+ E'
    Cupdate = 0.0 .* similar(Aprime)
    Aupdate = similar(Aprime)

    # Iteration to find optimal consumption
    normdiff = Inf
    t = 1

    while normdiff > tol && t <= maxT

        # Find Ctilda from the Euler Equation
        B = β .* up(C) * Π'
        Ctilda = upinv(B ./ q)

        # Solve for A from the budget constraint
        A = Ctilda + Aprime .* q .- E'

        # Adjust for whether budget constraint binds (only low shock)
        binds = Aprime .< A[:,1]
        Cupdate[binds] .= (Aprime .+ E' .- Aprime[:,1] .* q)[binds]
        # Generate auxiliary matrix for when binding
        auxbinding = similar(A)
        for e = 1:e_size
            auxfunc = LinearInterpolation(A[:, e], Ctilda[:, e])
            auxbinding[:, e] = auxfunc.(Aprime[:, e])
        end
        # Fill Cupdate for binding
        Cupdate[.!binds] .= auxbinding[.!binds]

        # values for next iteration
        t += 1
        normdiff = maximum(abs.(Cupdate - C))
        C = copy(Cupdate)

    end

    # Now find the optimal policy function in my original grid
    pol_func = zeros(a_size, e_size)
    for e = 1:e_size
        auxfunc = LinearInterpolation(A[:, e], Aprime[:, e])
        pol_func[:, e] = auxfunc.(Aprime[:, e])
    end

    return pol_func

end






#export RHSEuler, initialpolicy, EndGridSearch

#end

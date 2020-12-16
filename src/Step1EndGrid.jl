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

function RHSEuler(A, E, Π, q, β, σ, a, aprime; tol = 0.01)

    # a is the actual value of assets
    # aprime is the actual of the new assets

    """
    This function computes the RHS of the Euler Equation.
    Inputs:
        - A: assets today
        - E: exogenous shocks
        - Π: transtion matrix
        - q: price (guessed one)
        - β: discount factor
        - σ: preference parameter
        - a: actual value of assets today
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
    #up(a, e, ap, q) = ((a + e - ap*q) <= 0) ? -Inf : (a + e - ap*q).^(-σ)

    Uprime = zeros(e_size, 1)
    RHS = zeros(e_size, 1)

    for e = 1:e_size

        aux = LinearInterpolation(Agrid[:, e], Aprimetol[:, e])
        adoubleprime = aux(aprime)
        #adoubleprime = LinearInterpolation(Agrid[:, e], Aprimetol[:, e])(aprime)
        Uprime[e] = up(aprime .+ E[e] .- q .* adoubleprime)
        RHS[e] = up(a + E[e] - aprime * q) - β * dot(Π[e, :], Uprime)

    end

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

                FOCaprime(aprime) = RHSEuler(A, E, Π, q, β, σ, Agrid[a, e], aprime)

                if FOCaprime(Agrid[1, e])[e] > 0.0
                    Aprime[a, e] = Agrid[1, e]
                elseif FOCaprime(Agrid[end, e])[e] < 0.0
                    Aprime[a, e] = Agrid[end, e]
                else
                    if e == 1
                        auxfunc1(aprime) = RHSEuler(A, E, Π, q, β, σ, Agrid[a, e], aprime)[1]
                        Aprime[a, e] = fzero(auxfunc1, Agrid[175, e])
                    else
                        auxfunc2(aprime) = RHSEuler(A, E, Π, q, β, σ, Agrid[a, e], aprime)[2]
                        Aprime[a, e] = fzero(auxfunc2, Agrid[1, e], Agrid[end, e])
                    end
                end
                # 110-116 may not be working! Check how to adjust it.
                # Maybe use something different from a bisection method.
                #    Aprime[a, e] = fzero(FOCaprime, Agrid[1, e], Agrid[end, e])
                # end

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

    # Define the inverse of the first derivative of utility function
    # (in the future, make it more general)
    upinv(B) = B.^(-1/σ)

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
        binds = Aprime .< A[:,1]'
        Cupdate[binds] .= (Aprime .+ E' .- Aprime[:,1]' .* q')[binds]
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

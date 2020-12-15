# Replicate Huggett (1993)

# --------------------------------------------------------------------------- #
# Set parameters:
β = 0.9932
BoundSet = (-2.0, -4.0, -6.0, -8.0)

# Exogenous state variable grid
eH = 1.0
eL = 0.1
E = [eL; eH]

# Transition matrix: Pr(j|i) = Π(i,j) (columns sum to 1)
πHH = 0.925 # e_H | e_H
πHL = 0.5   # e_H | e_L
Π = [1 - πHL  πHL
    1 - πHH πHH]

# Guess initial price
q0 = 1.1


# --------------------------------------------------------------------------- #
# Solve for σ = 1.5
σ = 1.5

# Start loop over lower bound
for a_lb in BoundSet

    # Generate grid points for endogenous state variable
    a_max = 5
    a_size = 350
    A = range(a_lb, a_max, length = a_size)
    A = [A; ];

    (q, policy, Ψ) = step3(A, E, Π, q0, β, σ)

    println("- Solving the model for σ = $σ and lower bound of a = $a_lb.")
    println("  The optimal price is q = $q.")
    println("  ")

end


# --------------------------------------------------------------------------- #
# Solve for σ = 3
σ = 3

# Start loop over lower bound
for a_lb in BoundSet

    # Generate grid points for endogenous state variable
    a_max = 5
    a_size = 350
    A = range(a_lb, a_max, length = a_size)
    A = [A; ];

    (q, policy, Ψ) = step3(A, E, Π, q0, β, σ)

    println("- Solving the model for σ = $σ and lower bound of a = $a_lb.")
    println("  The optimal price is q = $q.")
    println("  ")

end



# --------------------------------------------------------------------------- #
# Graphs
q0 = 1.1
σ = 1.5
a_lb = -2.0
a_max = 5
a_size = 350
A = range(a_lb, a_max, length = a_size)
A = [A; ];

(q, policy, Ψ) = step3(A, E, Π, q0, β, σ)

# Find endogenous a upper bar (use aUB+20 in graphs to see the intersection)
aUB = findfirst(policy[:,2] - A[:] .<= 0.001)

# Unconstrained plot for Policy Function
p1 = plot(A, [policy A], linewidth = 2, title = "Policy function",
    label = ["eL" "eH" "45°"], legend =:topleft)

# Constrained  plot for Policy Function
p1B = plot(A[1:aUB+20], [policy[1:aUB+20,:] A[1:aUB+20]], linewidth = 2,
    title = "Policy function", label = ["eL" "eH" "45°"], legend =:topleft)

# Unconstrained plot for Stationary Distribution Ψ
p2 = plot(A, Ψ, linewidth = 2, title = "Stationary distribution",
    label = ["eL" "eH"], legend =:topleft)

# Constrained plot for Stationary Distribution Ψ
p2B = plot(A[1:aUB+20], Ψ[1:aUB+20,:], linewidth = 2,
    title = "Stationary distribution", label = ["eL" "eH"], legend =:topleft)

# Find rescaling value
scale = sum(Ψ[end, :])
F = Ψ ./ scale

# Unconstrained plot for Stationary Distribution Ψ
p3 = plot(A, F, linewidth = 2, title = "Stationary distribution",
    label = ["eL" "eH"], legend =:topleft)

# Constrained plot for Stationary Distribution Ψ
p3B = plot(A[1:aUB+20], F[1:aUB+20], linewidth = 2,
    title = "Stationary distribution", label = ["eL" "eH"], legend =:topleft)


include("Step3.jl")
include("Step2.jl")
include("Step1Grid.jl")

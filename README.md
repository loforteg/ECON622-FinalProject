# ECON622-FinalProject
Final Project for the course ECON622. Replication of Huggett (1993) ([paper](https://github.com/loforteg/ECON622-FinalProject/blob/main/literature/Huggett%201993)), with simple grid search and endogenous grid search method.

----
### Table of contents
[Model](#model)

[Computation](#computation)

[Files description](#files-description)

[References](#references)

----
## Model
The model encompasses the following elements:
* a continuum of agents with mass 1;
* in each period, the agent is endowed with ![formula](https://render.githubusercontent.com/render/math?math=e%20\in%20E=(e_h%3Be_l)) according to a Markov Process ![formula](https://render.githubusercontent.com/render/math?math=\pi(e'|e)=Pr(e_{t%2B1}=e'|e_t=e)%3E0) independent across agents;
* the utility function is ![formula](https://render.githubusercontent.com/render/math?math=E%5B\sum_{t=0}^{\infty}\beta^tu(c_t)%5D), with ![formula](https://render.githubusercontent.com/render/math?math=\beta\in(0,1)), ![formula](https://render.githubusercontent.com/render/math?math=u(c)=\frac{c^{1-\sigma}}{-\sigma}), and ![formula](https://render.githubusercontent.com/render/math?math=\sigma%3E1);
* the budget constraint in each period is ![formula](https://render.githubusercontent.com/render/math?math=c%2Ba'q%20\leq%20a%2Be) with ![formula](https://render.githubusercontent.com/render/math?math=a'\geq\underline{a}) and ![formula](https://render.githubusercontent.com/render/math?math=\underline{a}%3C0).

Define the individual state variable as ![formula](https://render.githubusercontent.com/render/math?math=x=(a,e)%20\in%20X) with X being the combination of sets A and E, and set A having ![formula](https://render.githubusercontent.com/render/math?math=\underline{a}) as lower bound.

Let ![formula](https://render.githubusercontent.com/render/math?math=q%3E0) be the constant price for next period's credit.
The decision problem for the individual is

![formula](https://render.githubusercontent.com/render/math?math=v(x%3Bq)=\max_{(c,a')%20\in%20\Gamma(x%3Bq)}u(c)%2B\beta\sum_{e'}v(a',e'%3Bq)\pi(e'|e))

where ![formula](https://render.githubusercontent.com/render/math?math=\Gamma(x%3Bq)=((c,a')%3Ac%2Bqa'%20\leq%20a%2Be%3Bc%20\geq%200%3Ba'%20\geq%20a)).

Let ![formula](https://render.githubusercontent.com/render/math?math=\overline{a}) be an endogenously determined level of credit balance that individual agents never overtake, then ![formula](https://render.githubusercontent.com/render/math?math=S=%5B\underline{a},\overline{a}%5DxE) and ![formula](https://render.githubusercontent.com/render/math?math=\beta_S) is the corresponding Borel sigma-algebra.
For all elements B of the Borel sigma-algebra, ![formula](https://render.githubusercontent.com/render/math?math=\psi(B)) is the mass of agents whose individual state vector lies in B.

Define the transition probabilities P(x,B) as the probability that an agent with state x will have a state vector lying in B in the next period.


A *stationary equilibrium* is defined as c(x), a(x), q, and ![formula](https://render.githubusercontent.com/render/math?math=\psi) such that:
* c(x) and a(x) are optimal decision rules, given q;
* markets clear
* ![formula](https://render.githubusercontent.com/render/math?math=\psi) is stationary

----
## Computation
The computational procedure is broken down in 3 steps:
1. Given q, compute the optimal policy function a(x;q) by iterating
![formula](https://render.githubusercontent.com/render/math?math=(Tv)(x%3Bq)=\max_{(c,a')%20\in%20\Gamma(x%3Bq)}u(c)%2B\beta\sum_{e'}v(a',e'%3Bq)\pi(e'|e))
1. Given the optimal policy function, compute the stationary distribution.
1. Update the price q and repeat the first two steps until convergence is found.

In this repository, I am going to perform the first step in two different ways and compare their performance:
* simple grid search, as in the original Huggett (1993) paper;
* endogenous grid search.

----
## Files description
In the src folder you can find the following files:
* [Step1Grid.jl](https://github.com/loforteg/ECON622-FinalProject/blob/main/src/Step1Grid.jl): this file performs the first step described in the Computation section using a simple grid search method.
Specifically, in each point of the grid A in which consumption is weakly positive, I compute the value function V and then find the value of a' that corresponds to the maximum value function attainable by the individual starting from a given point in the grid.
When consumption is negative, I replace the value function with - ![formula](https://render.githubusercontent.com/render/math?math=\infty).
I then repeat the process until convergence of the maximum value function, that is convergence of the policy function.

* [Step1EndGrid.jl](https://github.com/loforteg/ECON622-FinalProject/blob/main/src/Step1EndGrid.jl): this file performs the first step described in the Computation section using an endogenous grid search method.
The essential idea of the endogenous grid method is that I first need to construct a grid of next period's asset holdings and use that, rather than the grid of current assets holdings, in order to find the optimal policy function.
Specifically, given an initial guess for the policy function, for any pair of a' and e, I construct the right hand side of the Euler Equation and use it to solve for the value of consumption that satisfies the Euler Equation.
Finally, from the budget constraint I derive the optimal choice of assets in the current period.
I use this solution to update the initial guess, as long as I get convergence of the policy function.

For more information on endogenenous grid method, please refer to the following [notes](https://github.com/loforteg/ECON622-FinalProject/blob/main/literature/Notes%20on%20Endogenous%20Grid%20Method.pdf).

* [Step2.jl](https://github.com/loforteg/ECON622-FinalProject/blob/main/src/Step2.jl): this file performs the second step described in the Computation section.
In practice, I first define function `pol_inv` which computes the inverse of the policy function.
I compute it as the element in the grid of A such that `pol_func(A, E)` is the closest to a'.
However, note some adjustments:
  * if a' is smaller than all the values of the policy function, then we set a to ![formula](https://render.githubusercontent.com/render/math?math=\underline{a});
  * if e' is the best shock and a' is bigger than any value of the optimal policy function, then we set a to the endogenous ![formula](https://render.githubusercontent.com/render/math?math=\overline{a}) corresponding to the intersection with the 45° line.

Once I have defined the inverse of the policy function, I generate the function `step2` which finds the stationary distribution ![formula](https://render.githubusercontent.com/render/math?math=\Psi).
As stated above, the code computes the updated distribution using the initial guess distribution, the inverse of the policy function, and the transition probabilities.
The code then iterates over the updated distributions until they converge to the same value.
For the initial guess distribution, note that ![formula](https://render.githubusercontent.com/render/math?math=\Psi) is the cumulative distribution over any (a,e).
Therefore, the initial guess distribution must reach value 1 over all the grid points for a and e.
I am assuming that individuals split equally across the two values of the exogenous state variable and that then they distribute uniformly across the values of a.
In practice, the initial guess cumulative distribution function is given by a range going from 0.0 to 0.5 with as many steps as the number of gridpoints preceding the endogenous ![formula](https://render.githubusercontent.com/render/math?math=\underline{a}) for each value of e.
Summing across the two values of e then gives 1.

* [Step3.jl](https://github.com/loforteg/ECON622-FinalProject/blob/main/src/Step3.jl): this file performs the third step described in the Computation section.
First, define function `mktclearing` that computes the demand for assets using the gridpoints, the policy function, and the stationary distribution as inputs.
Recall that, for market clearing, we only need the demand for assets to be equal to 0.
Then, define the function `step3` that, starting from a guess initial value of prices, solves for the previous two steps and iterates up until the market for assets sums to zero.
I try not to have a "jumpy" update in prices by using the following updating procedure:
  * if there is excess demand, the new price is equal to the old price plus a fraction of the tolerance level to determine that there is no excess demand;
  * if there is excess supply, the new price is equal to the old price diminished by the tolerance level.
  
  
* [Replication.jl](https://github.com/loforteg/ECON622-FinalProject/blob/main/src/Replication.jl): this file replicates the main table and figures reported in Huggett (1993).


----
## References
Huggett, Mark, *The risk-free rate in heterogeneous-agent incomplete-insurance economies*, Journal of economic Dynamics and Control 17.5-6 (1993): 953-969.

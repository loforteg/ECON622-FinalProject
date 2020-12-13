# ECON622-FinalProject
Final Project for the course ECON622. Replication of Huggett (1993), with simple grid search and endogenous grid search method.

----
### Table of contents
[Model](#model)

[Computation](#computation)

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
## References
Huggett, Mark, *The risk-free rate in heterogeneous-agent incomplete-insurance economies*, Journal of economic Dynamics and Control 17.5-6 (1993): 953-969.

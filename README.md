# ECON622-FinalProject
Final Project for the course ECON622. Replication of Huggett (1993), with simple grid search and endogenous grid search method.

----

## Model
The model encompasses the following elements:
* a continuum of agents with mass 1;
* in each period, the agent is endowed with ![formula](https://render.githubusercontent.com/render/math?math=e%20\in%20E=(e_h%3Be_l)) according to a Markov Process ![formula](https://render.githubusercontent.com/render/math?math=\pi(e'|e)=Pr(e_{t%2B1}=e'|e_t=e)%3E0) independent across agents;
* the utility function is ![formula](https://render.githubusercontent.com/render/math?math=E%5B\sum_{t=0}^{\infty}\beta^tu(c_t)%5D), with ![formula](https://render.githubusercontent.com/render/math?math=\beta\in(0,1)), ![formula](https://render.githubusercontent.com/render/math?math=u(c)=\frac{c^{1-\sigma}}{-\sigma}), and ![formula](https://render.githubusercontent.com/render/math?math=\sigma%3E1);
* the budget constraint in each period is ![formula](https://render.githubusercontent.com/render/math?math=c%2Ba'q%20\leq%20a%2Be) with ![formula](https://render.githubusercontent.com/render/math?math=a'\geq\underline{a}) and ![formula](https://render.githubusercontent.com/render/math?math=\underline{a}%3C0).

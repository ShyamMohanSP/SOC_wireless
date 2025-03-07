This repository contains the official implementation of the methods described in our paper:​

**​Optimal Power Procurement for Green Cellular Wireless Networks under Uncertainty and Chance Constraints**
Nadhir Ben Rached, Shyam Mohan Subbiah Pillai, and Raúl Tempone
arXiv:2503.03051​

Abstract

Given the increasing global emphasis on sustainable energy usage and the rising energy demands of cellular wireless networks, this work seeks an optimal short-term, continuous-time power procurement schedule to minimize operating expenditure and the carbon footprint of cellular wireless networks equipped with energy storage capacity and hybrid energy systems comprising uncertain renewable energy sources. Despite the stochastic nature of wireless fading channels, the network operator must ensure a certain quality-of-service (QoS) constraint with high probability. This probabilistic constraint prevents using the dynamic programming principle to solve the stochastic optimal control problem. This work introduces a novel time-continuous Lagrangian relaxation approach tailored for real-time, near-optimal energy procurement in cellular networks, overcoming tractability problems associated with the probabilistic QoS constraint. The numerical solution procedure includes an efficient upwind finite-difference solver for the Hamilton–Jacobi–Bellman equation corresponding to the relaxed problem, and an effective combination of the limited memory bundle method (LMBM) for handling nonsmooth optimization and the stochastic subgradient method (SSM) to navigate the stochasticity of the dual problem. Numerical results, based on the German power system and daily cellular traffic data, demonstrate the computational efficiency of the proposed numerical approach, providing a near-optimal policy in a practical timeframe.

Repository Overview

This repository provides the following components:​

Source Code: MATLAB implementation of the proposed Lagrangian relaxation approach for optimal power procurement.​

Numerical Solver: Upwind finite-difference solver for the Hamilton–Jacobi–Bellman equation.​

Optimization Methods: Implementation of the Limited Memory Bundle Method (LMBM) and Stochastic Subgradient Method (SSM).​

Data: Sample datasets of 2023-24 wind power production and forecasts of 2023-24, and 2024 energy spot-prices from 50Hertz (https://www.50hertz.com).​

Main source file to run is ell_adaptive_alg.m. Follow the instructions in https://napsu.karmitsa.fi/lmbm/ to install the LMBM driver for MATLAB.

If you find this repository useful in your research, please cite our paper:​

@article{BenRached2025Optimal,

  title     = {Optimal Power Procurement for Green Cellular Wireless Networks under Uncertainty and Chance Constraints},
  
  author    = {Ben Rached, Nadhir and Subbiah Pillai, Shyam Mohan and Tempone, Raúl},
  
  journal   = {arXiv preprint arXiv:2503.03051},
  
  year      = {2025},
  
  url       = {https://arxiv.org/abs/2503.03051}
  
}

This work was supported by the KAUST Office of Sponsored Research (OSR) under Award No. URF/1/2584-01-01 and the Alexander von Humboldt Foundation. It was also partially performed as part of the Helmholtz School for Data Science in Life, Earth and Energy (HDS-LEE) and received funding from the Helmholtz Association of German Research Centres.​

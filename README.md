# CBO for Saddle Point Problems (CBO-SP)
Implementation of a consensus-based optimization method for saddle point problems (CBO-SP).

CBO-SP is a novel multi-particle metaheuristic derivative-free optimization method capable of provably finding global Nash equilibria. Following the idea of swarm intelligence, the method employs a group of interacting particles, which perform a minimization over one variable and a maximization over the other.

Version 1.0

Date 23.12.2022

------

## R e f e r e n c e s

### Consensus-Based Optimization for Saddle Point Problems

https://arxiv.org/abs/2212.12334

by

- Hui &nbsp; H u a n g &nbsp; (University of Graz), 
- Jinniao &nbsp; Q i u &nbsp; (University of Calgary),
- Konstantin &nbsp; R i e d l &nbsp; (Technical University of Munich, Munich Center for Machine Learning)

------

## D e s c r i p t i o n

MATLAB implementation of a consensus-based optimization method for saddle point problems.

For the reader's convenience we describe the folder structure in what follows:

BenchmarkFunctions
* objective_function.m: objective function generator
* ObjectiveFunctionPlot.m: plotting routine for objective function

EnergyBasedCBOAnalysis
* analyses: convergence and parameter analyses of CBO-SP
    * CBOSPNumericalExample.m: testing script
* CBOSP: code of CBO-SP optimizer
    * compute_consensus.m: computation of consensus point
    * CBOSP_update: one CBO-SP step
    * CBOSP.m: CBO-SP optimizer
* visualizations: visualization of the CBO-SP dynamics
    * CBOSPIllustrative.m: Illustration of the CBO-SP at work

------

## C i t a t i o n s

```bibtex
@article{CBOSPQiuHuangRiedl,
      title = {Consensus-Based Optimization for Saddle Point Problems},
     author = {Jinniao Qiu and Hui Huang and Konstantin Riedl},
       year = {2022},
    journal = {arXiv preprint arXiv:2212.12334},
}
```

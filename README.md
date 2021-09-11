# MPLS (Minneapolis) Framework for Robust Rotation Averaging

## Introduction and References

MPLS is a powerful alternative to the popular IRLS (Iteratively Reweighted Least Squares) algorithm for solving the general problem of group synchronization (and rotation averaging as a special case).

This repo contains matlab files for implementing the method of the following papers

[Message Passing Least Squares Framework and Its Applications in Rotation Synchronization](https://arxiv.org/pdf/2007.13638.pdf), Yunpeng Shi and Gilad Lerman, ICML 2020.

[Robust Group Synchronization via Cycle-Edge Message Passing](https://arxiv.org/pdf/1912.11347.pdf), Gilad Lerman and Yunpeng Shi, arXiv preprint, 2019.

If you would like to use our code for your paper, please cite the above two works.

## Usage
Download matlab files to the same directory. Checkout and run the following demo code. 
```
Examples/Compare_algorithms.m
```
The following is a sample output (error in degrees) under Nonuniform and Adversarial corruption, which shows great advantage of CEMP and MPLS:

```
Algorithms     MeanError     MedianError
___________    __________    ___________

"Spectral"         13.745        4.4287 
"IRLS-GM"           34.54        1.9641 
"IRLS-L0.5"        5.0656       0.96395 
"CEMP+MST"     2.2115e-06    2.4148e-06 
"CEMP+GCW"       0.010659     0.0034167 
"MPLS"          0.0001131     5.622e-05 

```
## Various Corruption Models
We provide 5 different corruption models. 3 for nonuniform topology and 2 for uniform toplogy (see ``Uniform_Topology.m`` and ``Nonuniform_Topology.m``). Uniform/Nonuniform toplogy refers to whether the corrupted subgraph is Erdos Renyi or not. In other words, the choice of Uniform/Nonuniform toplogy decides how to select edges for corruption. In ``Uniform_Topology.m``, two nodes are connected with probability ``p``. Then edges are independently drawn with probability ``q`` for corruption. In ``Nonuniform_Topology.m``, two nodes are connected with probability ``p``. Then with probability ``p_node_crpt`` a node is selected so that its neighboring edges are candidates for corruption. Next, for each selected node, with probability ``p_edge_crpt`` an edge (among the neighboring edges of the selected node) is corrupted. This is a more malicious scenario where corrupted edges have cluster behavior (so local coruption level can be extremely high). 

One can also optionally add noise to inlier graph and outlier graph through setting ``sigma_in`` and ``sigma_out``.

The argument ``crpt_type`` in the two functions determines how the corrupted relative rotations are generated for those selected edges. In ``Uniform_Topology.m``, there are 2 options of ``crpt_type``: ``uniform`` and ``self-consistent``.
In ``Nonuniform_Topology.m``, there are the following 3 options of ``crpt_type``.

``uniform``: The corrupted relative permutations <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{R_{ij}}"> are i.i.d follows uniform distribution over the space of SO(3).

``self-consistent``: The corrupted <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{R_{ij}}"> are relative rotations of another set of absolute rotations. Namely <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{R_{ij} = R_i^{crpt} R_j^{crpt}'}"> where those absolute rotations are different from the ground truth and are i.i.d drawn from the uniform distribution in the space of SO(3). In this way, the corrupted relative rotations are also cycle-consistent.

``adv``: Extremely malicious corruption that replaces the underlying absolute rotations from ground truth <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{R_i^*}"> to <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{R_i^{crpt}}">. Namely <img src="https://render.githubusercontent.com/render/math?math=\color{red} \mathbf{R_{ij} = R_i^{crpt} R_j^{* }'}">. Additional high noise is added to the corruption, otherwise the recovery of the ground truth can be ill-posed.

The demo code uses the following function for implementing MPLS:
```
MPLS(Ind, RijMat, CEMP_parameters, MPLS_parameters)
```
Each row of ``Ind`` matrix is an edge index (i,j). The edge indices (the rows of Ind) MUST be sorted in ``row-major order``. That is, the edge indices are sorted as (1,2), (1,3), (1,4),..., (2,3), (2,5), (2,8),..., otherwise the code may crash when some edges are not contained in any 3-cycles. Make sure that i<j. If some edges have indices (3,1), then change it to (1,3) and take a transpose to the corresponding Rij.


## Choices of Parameters
Please also see in the above file for different choices of parameters that we recommend. MPLS is not sensitive to those parameters. The general rule for reweighting parameters
```
CEMP_parameters.reweighting
MPLS_parameters.reweighting
```
is the following:

As the iteration # increases, both of them increases (or at least be nondecreasing) and then stop at a constant between 20-50 (cannot be arbitrarily large due to numerical instability). 

The following parameter determines the fraction of cycle-consistency information that MPLS takes account.
```
MPLS_parameters.cycle_info_ratio
```
That is, at each iteration
```
Estimated corruption level of each edge
= (1-MPLS_parameters.cycle_info_ratio) * Residual + MPLS_parameters.cycle_info_ratio * (Cycle inconsistency measure)

```

For sparse graphs (in many real scenarios) we recommend (rely more on residuals)
```
MPLS_parameters.reweighting = CEMP_parameters.reweighting(end);
MPLS_parameters.cycle_info_ratio = 1./((1:MPLS_parameters.max_iter)+1); % in the end ignore 3-cycle information
```

For very dense graphs, we find the following parameters perform better (rely more on 3-cycle information)
```
MPLS_parameters.reweighting = 0.1*1.5.^((1:15)-1);    % more and more aggressive reweighting
MPLS_parameters.cycle_info_ratio = 1-1./((1:MPLS_parameters.max_iter)+1); % in the end focus only on 3-cycle information
```


## Dependencies
The following files are dependencies for running Lie-Algebraic Averaging method that were written by AVISHEK CHATTERJEE (revised and included in this repo). See also [Robust Rotation Averaging](http://www.ee.iisc.ac.in/labs/cvl/papers/robustrelrotavg.pdf) and [Efficient and Robust Large-Scale Rotation Averaging](https://www.cv-foundation.org/openaccess/content_iccv_2013/papers/Chatterjee_Efficient_and_Robust_2013_ICCV_paper.pdf) for details.
```
Weighted_LAA.m
Build_Amatrix.m
R2Q.m
q2R.m
BoxMedianSO3Graph.m
L12.m
RobustMeanSO3Graph.m
```

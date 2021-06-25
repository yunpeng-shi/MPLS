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
Demo_MPLS.m
```
The demo code uses the function
```
MPLS_rotation(Ind, RijMat, CEMP_parameters, MPLS_parameters)
```
Each row of ``Ind`` matrix is an edge index (i,j). The edge indices (the rows of Ind) must be sorted in ``row-major order``. That is, the edge indices are sorted as (1,2), (1,3), (1,4),..., (2,3), (2,5), (2,8),..., otherwise the code may crash when some edges are not contained in any 3-cycles. Make sure that i<j. If some edges have indices (3,1), then change it to (1,3) and take a transpose to the corresponding Rij.

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

In general, one can take this parameter to 0 as iteration # increases. However, for some dense graphs one may let it approach 1 and gradually ignore the residual information. See Demo_MPLS.m for details. 


## Dependencies
The following files are dependencies for running Lie-Algebraic Averaging method that were written by AVISHEK CHATTERJEE (revised and included in this repo). See also [Robust Rotation Averaging](http://www.ee.iisc.ac.in/labs/cvl/papers/robustrelrotavg.pdf) and [Efficient and Robust Large-Scale Rotation Averaging](https://www.cv-foundation.org/openaccess/content_iccv_2013/papers/Chatterjee_Efficient_and_Robust_2013_ICCV_paper.pdf) for details.
```
Weighted_LAA.m
Build_Amatrix.m
R2Q.m
q2R.m
```

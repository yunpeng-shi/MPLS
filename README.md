# MPLS (Minneapolis) Framework for Robust Rotation Averaging

## Intro and References

This repo contains matlab files for implementing the method of the papers

[Message Passing Least Squares Framework and Its Applications in Rotation Synchronization](https://arxiv.org/pdf/2007.13638.pdf), Yunpeng Shi and Gilad Lerman, ICML 2020.
[Robust Group Synchronization via Cycle-Edge Message Passing](https://arxiv.org/pdf/1912.11347.pdf), Gilad Lerman and Yunpeng Shi, arXiv preprint, 2019.

MPLS is a powerful alternative for the popular IRLS (Iteratively Reweighted Least Squares) algorithm for solving the general problem of group synchronization (and rotation averaging as a special case).

## Usage
Download matlab files to the same directory. Checkout and run the following demo code. 
```
Demo_MPLS.m
```
## Choices of Parameters
Please also see in the above file for different choices of parameters that we recommend. MPLS is not sensitive to those parameters. The general rule for reweighting parameters
```
CEMP_parameters.reweighting
MPLS_parameters.reweighting
```
is that:
as the iteration # increases, both of them increases (or at least nondecreasing) and then stops at a fixed number between 20-50.

## Dependencies
The following files are dependencies for running Lie-Algebraic Averaging method that were written by AVISHEK CHATTERJEE (revised and included in this repo). See also [Robust Rotation Averaging](http://www.ee.iisc.ac.in/labs/cvl/papers/robustrelrotavg.pdf) and [Efficient and Robust Large-Scale Rotation Averaging](https://www.cv-foundation.org/openaccess/content_iccv_2013/papers/Chatterjee_Efficient_and_Robust_2013_ICCV_paper.pdf) for details.
```
Weighted_LAA.m
Build_Amatrix.m
R2Q.m
q2R.m
```

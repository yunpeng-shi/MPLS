# MPLS (Minneapolis) Framework for Robust Rotation Averaging

This repo contains matlab files for implementing the method of the paper

[Message Passing Least Squares Framework and Its Applications in Rotation Synchronization, Yunpeng Shi and Gilad Lerman, ICML 2020](https://arxiv.org/pdf/2007.13638.pdf).

MPLS is a powerful alternative for the popular IRLS (Iteratively Reweighted Least Squares) algorithm. 

## Usage
Download matlab files to the same directory. Run the following demo code.
```
Demo_MPLS.m
```
Please also see in above file for different choices of parameters that we recommend.


The following files are dependencies for running Lie-Algebraic Averaging method that were written by AVISHEK CHATTERJEE (Included in this repo). See also [Robust Rotation Averaging](http://www.ee.iisc.ac.in/labs/cvl/papers/robustrelrotavg.pdf) and [https://www.cv-foundation.org/openaccess/content_iccv_2013/papers/Chatterjee_Efficient_and_Robust_2013_ICCV_paper.pdf](Efficient and Robust Large-Scale Rotation Averaging) for details.
```
Weighted_LAA.m;  R2Q.m; q2R.m; Build_Amatrix.m
```

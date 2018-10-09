# Sphere-Sampling-reduction
A lattice sampling algorithm


This is a implementation of a sampling algorithm 
described in the following article:Random Sampling on an High-Dimensional Sphere for Solving SVP.  This process is a sampling reduction algorithm for solving SVP.  The algorithm is based on a new algorithm named Sphere Sampling Reduction. Although this is a algoithm which can be parallel implement. But here we show only a  serial variant. This implementation in C++ requires NTL and is strictly a research-oriented
implementation; The input basis should be given in a .txt and in our test our basis comes from http://www.latticechallenge.org/svp-challenge/ which is a website for SVP challenge and can generate random lattice basis online. The reduction part is comes from NTL's BKZ implementation. And can change its parameters based on requirement. In our algorithm, we made the blocksize of BKZ is 20.


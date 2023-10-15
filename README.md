# Interpolation Filters for MATLAB

This repository includes MATLAB functions to perform interpolation of a discrete signal. Built-in function [`mkpp`](https://www.mathworks.com/help/matlab/ref/mkpp.html) is used to construct a piecewise polynomial. Included methods are:
1. `shlinear.m` — shifted linear interpolation [(Blu et al., 2004)](https://doi.org/10.1109/TIP.2004.826093).
2. `bspline.m` — cubic B-spline interpolation [(Unser, 1999)](https://doi.org/10.1109/79.799930).
3. `cspline.m` — causal cubic spline interpolation in [(Petrinovic, 2008)](https://doi.org/10.1109/TSP.2008.929133) formulation.
4. `cspline2.m` — causal cubic spline interpolation in [(Meinsma et al., 2012)](https://doi.org/10.1109/TSP.2012.2185228) formulation.
5. `moms.m` — cubic MOMS interpolation [(Blu et al., 2004)](https://ieeexplore.ieee.org/document/7079929).

All of the functions operate on the whole data vector. Even though most of the methods are causal, provided functions are designed for evaluation of interpolation filters rather then for real-time implementation.  

Overview of these interpolation kernels is provided in Chapter 2 of [my MSc thesis](https://www.researchgate.net/publication/374194858_Interpolation_Filters_for_Antiderivative_Antialiasing).

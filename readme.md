 

 # SDPDAL beta 1.0
 A MATLAB software for solving Low-rank SDP.

 # Problems and solvers
 - The package contains codes for optimization problems on manifold:   
    ​        $$\min  f(X) + h(X)  s.t. X \in D, A(X) = b,$$ 
    
    where $D = \{X: B(X)=b_0, X\succeq 0\}$, which has the manifold structure by factorizing $X = R^TR$. 
 - This code is currently in the experimental phase and subject to further development and refinement. 
 - The manifold structure of SDPDAL is taken from the ARNT and manopt package.  We thank the authors for kindly sharing their codes. 

Applications have been solved by these solvers:
- SDP relaxation of clustering problems.
  $$\min \;  <-W,X>\;  s.t., Xe = e, Tr(X) = K, X \geq 0, X\succeq 0$$
- Maxcut SDP: $\min  \mathrm{Tr}(CX), s.t., X_{ii}=1, X \succeq 0$
- Low-Rank Nearest Correlation  Estimation: $ \min_{ X \succeq 0} \; \frac{1}{2} \| H \odot (X - C) \|_F^2, \; X_{ii} = 1, \; i = 1, \ldots, n, \; \mathrm{rank}(X) \le p.$
- SDP relaxation of sparse PCA
- Theta problems

## Summary of Datasets 
We list the downloading links to the datasets used in the papers for reproduction:
* Gset, NCM, RCP and Theta: http://faculty.bicmr.pku.edu.cn/~wenzw/code/sdp_data.zip

# References
- [Yifei Wang, Kangkang Deng, Haoyang Liu, Zaiwen Wen, A Decomposition Augmented Lagrangian Method for Low-rank Semidefinite Programming, SIAM Journal on Optimization, Vol. 33, No. 3, 1361-1390, 2023](https://arxiv.org/pdf/2109.11707.pdf)

- [Jiang Hu, Andre Milzarek, Zaiwen Wen, Yaxiang Yuan. Adaptive Quadratically Regularized Newton Method for Riemannian Optimization. SIAM Journal on Matrix Analysis and Applications, Vol. 39, No. 3, pp. 1181–1207](https://epubs.siam.org/doi/10.1137/17M1142478)

- [Nicolas Boumal , Bamdev Mishra, P.-A. Absil and Rodolphe Sepulchre. Manopt, a Matlab Toolbox for Optimization on Manifolds. Journal of Machine Learning Research (2014) 1455-1459](http://jmlr.org/papers/v15/boumal14a.html)




 # The Authors
 We hope that the package is useful for your application.  If you have any bug reports or comments, please feel free to email one of the toolbox authors:
 * Haoyang Liu, liuhaoyang at pku.edu.cn
 * Kangkang Deng, dengkangkang at pku.edu.cn
 * Zaiwen Wen, wenzw at pku.edu.cn

 # Installation 

 `>> cd example` 

 `>> test_ncm`


 # Copyright
-------------------------------------------------------------------------
   Copyright (C) 2022, Haoyang Liu, Kangkang Deng, Zaiwen Wen

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without  even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>


### Optimal Diagonal Preconditioning Toolbox

This repository contains implementations and detailed experiment results for the paper [Scalable Approximate Optimal Diagonal Preconditioning](https://arxiv.org/abs/2312.15594)

**Prerequisites**

The toolbox relies on solving semidefinite programs and needs to invoke optimization solvers through CVX. For finding the optimal digaonal preconditioner, users need to install

- CVX  MATLAB toolbox  at http://cvxr.com/cvx/

For finding the *subspace* optimal preconditioner with semi-infinite programming, users need to additionally install

- Gurobi MATLAB interface at https://docs.gurobi.com/projects/optimizer/en/current/reference/matlab.htmll

**Install**

To install the package, simply clone the repo with

```
git clone https://github.com/Gwzwpxz/opt_subspacepcd
```

and run 

```
optpcd_setup
```

in MATLAB. Upon installation, a toy example will be solved to test the installation.

```
Setting up the optimal diagonal preconditioning repo 
The repo depends on:  

> CVX   matlab toolbox    at http://cvxr.com/cvx/

To run semi-infinite programming-based approach, install 

> Gurobi MATLAB at https://docs.gurobi.com/projects/optimizer/en/current/reference/matlab.html

Running a toy example.

Using the semi-infinite approach for SDP
Condition number of X                            : 1.47e+03 
Condition number of DX for subspace constrained D: 1.47e+03 

Installation completes. Check README.md for usage details. 
```

**Usage**

The optimal preconditioning toolbox provides several utilities that allow users to 

- Compute the optimal Left/Right/Two-sided preconditioner
- Compute the subspace optimal preconditioner using semi-infinite programming

The basic usage  can be demonstrated by the following lines of code

```matlab
% Choose preconditioning type
param.ptype = 'S'; 
% Use Semi-inf solver to solve the problem 
param.method = "SINF";  
% Alternatively, use SDP solver to solve the problem
% param.method = "SDP";  

% Generate preconditioning problem. X is the user data
% Add identity and Jacobian preconditioner
n = size(X, 1);
param.subspace = [ones(n, 1), full(diag(X))];

% Get problem
DSprob = getoptprob(X, param);
% Solve problem
DSopt = optprecond(DSprob);

% Get performance
fprintf("Condition number of X                            : %5.2e \n", cond(full(X)));
fprintf("Condition number of DX for subspace constrained D: %5.2e \n", cond(full(DSopt.pX)));
```

**Contact**

Please contact `gwz@stanford.edu` for questions on the toolbox.

**Cite as**

> @article{gao2023scalable,
>   title={Scalable Approximate Optimal Diagonal Preconditioning},
>   author={Gao, Wenzhi and Qu, Zhaonan and Udell, Madeleine and Ye, Yinyu},
>   journal={arXiv preprint arXiv:2312.15594},
>   year={2023}
> }

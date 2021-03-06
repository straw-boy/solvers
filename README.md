# Google Summer of Code 2020
## solvers
The R package **solvers** is a collection of implementation of various algorithms to solve Sorted L-One Penalised Estimation(SLOPE) problem. This package was intended to be a benchmarking suite of the algorithms, the best of which could then be merged into the [SLOPE] package. Do checkout the articles on how to use this package for comparing the performances and benchmarking. Source code is available [here].

### Installation
You can install the package via devtools.
```sh
> devtools::install_github("straw-boy/solvers")
```

### Targets Achieved
* In addition to existing **FISTA** in SLOPE package, I managed to implement the following algorithms
  * **Alternating Direction Method of Multipliers (ADMM)**. Our implementation further provides the choice of optimization algorithm to solve the subproblem :         **Newton-Raphson**'s algorithm, Broyden–Fletcher–Goldfarb–Shanno (**BFGS**) algorithm and Limited-memory BFGS (**L-BFGS**) algorithm.
  * **Proximal Newton** algorithm. This also supports Quasi flavoured variant which uses FISTA to solve the sub problem.
* Benchmarked the algorithms (vignettes available in the article tab) and discovered usecases where FISTA is outperformed.
* I also discovered few bugs in the existing SLOPE package's code.

### Challenges Faced
- The biggest challenge was dealing with badly conditioned Hessian matrix. This problem first arose in Proximal Newton. I tried out various fixes and have implemented one such fix, but it is not fool-proof. This problem blew up in much severe way when I tried using Matrix Inversion Lemma (in both ADMM and PN) where in one instance, I encountered a temporary matrix whose inverse condition number was exactly 0. This was the major reason ADMM and PN lost out horribly to FISTA when number of features (p) was greater than number of training points (n) and hence benchmarking was done for only n>p cases to save compute time.
- As there are a lot of moving parts and algorithms are primarily numerical, a lot of time went into debugging the implementations.
- In the Proximal Quasi-Newton algorithm, solving the subproblem is highly inefficient right now and I am yet to notice a practical case where PN is outperformed by PQN.
- When large scale problems of poisson regression are solved using second order algorithms, gradients blow up and also more often than not, Hessian turns out to be singular. This renders usage of ADMM(BFGS/L-BFGS) and PN highly unstable. This is the reason I haven't added benchmarks for poisson regression.
- For multinomial family, there is no closed form expression for computing Hessian. I had to used a rather inefficient method there.

### Future Work
This still room for future work in SLOPE+solvers package.
- If one could find a good strategy to tackle singular Hessian, it would solve major problems in both ADMM and PN.
- Due to time constraint, I couldn't merge the solvers' code into SLOPE package. My recommendation would be to definitely consider Proximal Newton for gaussian and binomial likelohood functions when n is significantly greater than p as PN scales really well with n.
- Although PN scaled well with n in multinomial case as well, its performance degrades rapidly as p increases. There are two reasons for it: dimension of Hessian is actually p times number of classes so size of Hessian increases rather faster. Secondly, as mentioned above, the Hessian computation itself is also poor due to absense of closed form expression. So FISTA is undoubtedly the recommended choice here.
- solvers package's default values are not fined tuned for different algorithms right now. One possible line of work could be to find these default values which are in sync with each other.
- If there could be a way to incorporate caching of Gaussian Hessian in the current design, it would be a huge performance boost for gaussian case as Hessian won't change for it.
[SLOPE]: https://jolars.github.io/SLOPE/
[here]: https://github.com/straw-boy/solvers

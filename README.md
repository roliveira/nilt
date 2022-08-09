# Numerical Inverse Laplace Transform Methods

This repository is used to benchmark different algorithms for the numerical Laplace Transform Inversion. This work was partly developed in [Oliveira, R. (2021)](https://doi.org/10.25560/92253).

## Stehfest Algorithm

The algorithm described in [Stehfest, H. (1970)](https://doi.org/10.1145/361953.361969) is implemented in `StehfestAlgorithm` class. To exemplify the results, the entries entry for `f(t) = 1/sqrt(PI*t)` in Table 1 of the paper is obtained by numerically inverting the Laplace Transform of `F(s) = 1/sqrt(s)` and is reproduced below:

| t  | fa(t)          | fn(t)          | err               |
|----|----------------|----------------|-------------------|
| 1  | 0.564189583548 | 0.564186508758 | 5.44992307092e-06 |
| 2  | 0.398942280401 | 0.398939582368 | 6.76296805803e-06 |
| 3  | 0.325735007935 | 0.325734581605 | 1.3088245455e-06  |
| 4  | 0.282094791774 | 0.282093254379 | 5.44992307092e-06 |
| 5  | 0.252313252202 | 0.252312282254 | 3.84421978257e-06 |
| 6  | 0.230329432981 | 0.230329170603 | 1.13914346481e-06 |
| 7  | 0.213243618623 | 0.213243982566 | 1.70670322864e-06 |
| 8  | 0.199471140201 | 0.199469791184 | 6.76296805803e-06 |
| 9  | 0.188063194516 | 0.188063264428 | 3.71747548209e-07 |
| 10 | 0.178412411615 | 0.178411614082 | 4.47016807695e-06 |

## References

- Rodolfo Oliveira. 2021. *Modelling of reactive transport in porous media using continuous time random walks*. PhD Thesis (Mar. 2021). https://doi.org/10.25560/92253
- Harald Stehfest. 1970. *Algorithm 368: Numerical inversion of Laplace transforms [D5]*. Commun. ACM 13, 1 (Jan. 1970), 47–49. https://doi.org/10.1145/361953.361969

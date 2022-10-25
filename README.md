# copulaboost

[![](https://cranlogs.r-pkg.org/badges/grand-total/copulaboost)](https://cran.r-project.org/web/packages/copulaboost/index.html)
[![](https://cranlogs.r-pkg.org/badges/last-month/copulaboost)](https://cran.r-project.org/web/packages/copulaboost/index.html)

Additive copula regression for regression
    problems with binary outcome via gradient boosting 
    [Brant, Hobæk Haff (2022); <arXiv:2208.04669>]. The fitting process
    includes a specialised model selection algorithm for each component, where
    each component is found (by greedy optimisation) among all the D-vines with
    only Gaussian pair-copulas of a fixed dimension, as specified by the user.
    When the variables and structure have been selected, the algorithm then
    re-fits the component where the pair-copula distributions can be different
    from Gaussian, if specified.

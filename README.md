# copulaboost

[![](https://cranlogs.r-pkg.org/badges/grand-total/copulaboost)](https://cran.r-project.org/web/packages/copulaboost/index.html)
[![](https://cranlogs.r-pkg.org/badges/last-month/copulaboost)](https://cran.r-project.org/web/packages/copulaboost/index.html)

This package implements Additive copula regression for regression
    problems with binary outcome via gradient boosting, as detailed in
    [Brant, Hobæk Haff (2022); <arXiv:2208.04669>]. The fitting process
    includes a specialised model selection algorithm for each component, where
    each component is found (by greedy optimisation) among all the D-vines with
    only Gaussian pair-copulas of a fixed dimension, as specified by the user.
    When the variables and structure have been selected, the algorithm then
    re-fits the component where the pair-copula distributions can be different
    from Gaussian, if specified.

In addition to additive copula regression, implemented in the method copulaboost::copulaboost, this package also contains a standalone method that implements copula regression, where $(Y, \mathbf{X})$ can be modeled via a general R-vine (or a D-vine, if specified), and the conditional expectation computed from the vine copula model. In order to do this, i.e., to be able to compute the conditional expectation, the structure of the copula regression model is selected under the condition that $Y$ must be in the top node of the vine.

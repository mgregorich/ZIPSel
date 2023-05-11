<h1 align="center">  The non-negative Garrote for 'too' many zero-inflated predictors </h1>

 
## Background
Variable selection can be problematic in data obtained from mass spectronomy due to the large number of missing values, posing statistical challenges for standard analysis procedures. The available information on a given variable is usually split into two components: i) a binary indicator of the absence/presence of the entity in the sample and ii) its log2-transformed intensity if detected. We propose an extension of the non-negative Garrote (Breiman 1995) for variable selection in the presence of zero-inflation in high-dimensional data. 

## Methods
A simulation study is presented in order to assess the performance of difference variable selection procedures mainly based on penalized likelihood estimation when 'too' many zero-inflated predictors are available for prediction modelling.

The following methods are implemented:

+ the nonnegative garrote (ridge-garrote)
+ ridge-lasso penalization
+ lasso-ridge penalization
+ network-based group penalization
+ multiple univariable modelling ( no penalization)

<br>

**Implementation of the nonnegative Garrote (by Georg Heinze)**

Protogarrote: R software to fit prediction models with or without interactions using high-troughput proteomic biomarkers 

The following functions are available:

* `protogarrote` ... fits  a linear or logistic model with proteomics and clinical predictor variable. May include an interaction of a variable with the proteomics. Uses internal validation to determine the best values of lambda1 and lambda2.
* `coefficients.protogarrote` ... shows the estimated coefficients at the optimal lambda1, lambda2
* `predict.protogarrote` ... predicts outcomes using a `protogarrote` object and new data
* `plot.coefficients.protogarrote` ... shows the estimated coefficients for the proteomics by means of a plot and a table.

<br>

## References
Breiman, L. (1995). "Better subset regression using the nonnegative garrote." Technometrics 37(4): 373-384.

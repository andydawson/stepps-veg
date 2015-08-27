# stepps-veg

This repository contains the code to run the statistical model that estimate spaatial covariance parameters for forest composition at settlement. Covariance estimates are used to constrain vegetation composition estimates in the STEPPS prediction model. 

Covariance estimates are based on the PLS-estimates from the statistical composition model. Here, underlying composition of each taxon is modelled as a Gaussian process, and these processes are linked via a sum-to-one constraint.

Model is fit using Stan 2.6.2. Generated code has been modified to avoid automatic differentiation and use parallel computing via openMP.

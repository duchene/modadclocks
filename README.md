Assessing model adequacy for clock and substitution models in R
================================================================

David Duchene and Sebastian Duchene

david.duchene[at]anu.edu.au

sebastian.ducnene[at]sydney.edu.au

6 April 2015

Summary
-------

This application was designed to assess model adequacy of clock and substitution models used in phylogenetics. It requires the posterior of trees and parameter estimates, and provides the substitution model adequacy as well as the overall and branch-wise clock model adequacy.

The application is still in the process of being built. It was initially designed to test the clock model adequacy method on simulated data, so at present it is only suitable for the output of BEAST2. Analyses with a single gene are currently supported. Clock models supported are the strict clock, the uncorrelated lognormal clock, and the random local clock. Substitution models supported inlcude the JC, HKY, and GTR models, including gamma-distributed rates.
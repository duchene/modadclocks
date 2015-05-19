This repository was made to accompany the article 'Evaluating the adequacy of molecular-clock models using posterior predictive simulations'
============================================================================================================================================

David Duchene and Sebastian Duchene

david.duchene[at]anu.edu.au

sebastian.ducnene[at]sydney.edu.au

19 May 2015

Summary
-------

This repository was designed to assess the adequacy of clock and substitution models used in phylogenetics. It requires the posterior of trees and parameter estimates, and provides the substitution model adequacy, the overall (A index) and branch-wise clock model adequacy, the effect size of the branch wise adequacy, and the uncertainty in posterior predictive simulations.

A flexible application for assessing clock model adequacy is being built. The present repository is specifically to accmpany the article; it is designed to test the clock model adequacy method on simulated data, so it is only suitable for the output of analyses from BEAST2 and a single gene partition. Clock models supported are the strict clock, the uncorrelated lognormal clock, and the random local clock. Substitution models supported inlcude the JC, HKY, and GTR models, including gamma-distributed rates.

The following are a set of examples of how to use the functions and how to present the results.


Example of usage
----------------



Graphical examples of results
-----------------------------




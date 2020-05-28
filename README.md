ImmuneSignatures2
=================

This package provides the code needed to generate the expressionSet objects used in the ImmuneSignatures2 analysis by the Human Immunology Project Consortium by running the vignette "GenerateAnalysisExpressionSets.Rmd" that pre-processes data from the ImmuneSpace platform.

The output of this markdown file are rds objects "IS2_immdata.rds" and then expressionSets in eight variations of age (young / old), cross-study normalization status (normalized / non-normalized) and response status (with response / no response).  The expressionSets with response data includes response calls in the pData slot that have been discretized "low", "moderate", and "high" responder from HAI, NAb, or ELISA datasets (preference given in that order to subjects with multiple response data).

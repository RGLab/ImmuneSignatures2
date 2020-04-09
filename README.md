ImmuneSignatures2
=================

This package provides the code needed to generate the expressionSet objects used in the ImmuneSignatures2 analysis by the Human Immunology Project Consortium by running the vignette "GenerateAnalysisExpressionSets.Rmd" that pre-processes data from the ImmuneSpace platform.

The output of this markdown file are two saved rda objects "IS2_immdata.rda" and "IS2_esets.rda", which when loaded are the objects `immdata` and `esets`.  The `esets` object is a nested list of expressionSets that vary according to three parameters: 1 - with or without response, 2 - cross-study normalized and batch-corrected, and 3 - age (young or older).  To navigate this nested list, you can descend through the three options as shown: `esets$withResponse$norm$young`.  This will give you an expressionSet for the young cohort that has response calls that have been discretized "low", "moderate", and "high" responder from either HAI, NAb, or ELISA datasets.

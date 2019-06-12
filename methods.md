# Methods

For each comparison done the data table was restricted to only the counts from the samples that belong to the specific groups being compared. Then it was filtered to only include CDR’s that have a length of greater than or equal to 8 and less than or equal to 16. Then it was filtered for CDR's that have more than 10 counts in at least `(minGroup - 1)` samples, where minGroup is the size of the smaller group.

The differential analysis was done using R/Bioconductors edgeR package. The standard edgeR normalization function `calcNormFactors` was called. For the differential testing the GLM method was done . Since the data consisted of paired samples from each of several mice a paired test was done with the following design:

```{R}
design ← model.matrix( ~ mouse + group )
```

Where _group_ is the “treatment” grouping variable and _mouse_ is the mouse id. We then processed the data using the GLM version of the dispersion estimators

```{R}
d ← estimateGLMCommonDisp(d,design)
d ← estimateGLMTrendedDisp(d,design)
d ← estimateGLMTagwiseDisp(d,design)
```

and then for differential testing:

```{R}
fit ← glmFit(d,design)
results ← glmLRT(fit,coef=group.cont)
```

where group.cont is the contrast that selects for comparisions between groups.

For the purposes of computing average expression in a natural scale rescale the normalized counts to the original size of the dataset. This is done by setting a scaling factor that is equal to the geometric mean of the total sample counts for the samples. The table of significant CDR’s was filtered to an FDR less than 0.05.


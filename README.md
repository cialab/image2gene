# Deep Learning Predicts Gene Expression as an Intermediate Data Modality to Identify Susceptibility Patterns in Mycobacterium tuberculosis Infected Mouse Models 

## Background
Machine learning has seen sustained and successful application to many diagnostic and prognostic problems in computational histopathology, and yet relatively few efforts have been made to model gene expression from histopathology. This study proposes a methodology which predicts selected gene expression values (microarray) from hematoxylin and eosin whole-slide images as an intermediate data modality to identify fulminant-like tuberculosis ('supersusceptible') in an experimentally infected cohort of diversity outbred mice (n=77).

## Methods
Gradient-boosted trees were utilized as a novel feature selector to identify genes predictive of fulminant-like tuberculosis. A novel attention-based multiple instance learning model for regression was used to predict selected genes' expression from whole-slide images. Gene expression predictions were shown to be sufficiently replicated to identify supersusceptible mice using gradient-boosted trees trained on ground truth gene expression data.

## Findings
Overall, the model was accurate, showing high positive correlations with ground truth gene expression on both cross-validation (n = 77) and external testing sets (n = 33). Furthermore, the sensitivity and specificity of utilizing gene expression predictions to identify supersusceptible mice (n = 77) were 0.88 and 0.95, respectively, and for an external set of mice (n = 33) 0.88 and 0.93, respectively.

## Implications
Our methodology can be used to map histopathology to gene expression with sufficient accuracy to predict a clinical outcome. The proposed methodology exemplifies a computational template for gene expression panels, in which relatively cheap tissue histopathology may be mapped to specific genes' expression to serve as a diagnostic or prognostic tool.


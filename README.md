# Machine Learning Classification Framework of Lipids in Untargeted Metabolomics Using Mass-to-Charge Ratios and Retention Times
Authors: Christelle Colin-Leitzinger, Yonatan Ayalew Mekonnen, Isis Narvaez-Bandera, Vanessa Y. Rubio, Dalia Ercan, 
Eric A. Welsh, Lancia N.F. Darville, Min Liu, Hayley D. Ackerman, Julian Avila-Pacheco, Clary B. Clish, Kevin Hicks, 
John M. Koomen, Nancy Gillis, Brooke L. Fridley, Elsa R. Flores, Oana A. Zeleznik, Paul A. Stewart.

# Objective
To explore the potential of using the mass-to-charge ratio (m/z) and retention time (RT) from LC-MS data as standalone 
predictors for metabolite classification and propose a modeling framework which can be implemented internally on standalone datasets.

# Data availability
The data necessary to reproduce the framework presented in our manuscript are downloadable 
through the Zenodo website https://doi.org/10.5281/zenodo.17237416. 
The code to download them is present in the R script 01.data_preparation.R.

# How to run
To reproduce the modeling framework used in our manuscript, four scripts need to be run.
First, the R script 01.data_preparation.R applies the cleaning, filtering and annotation steps as calculating the
solvent composition to substitute RT.
Then the script 02.modeling.R uses the `tidymodels` R package to split the data, create the 10-fold cross-validation 
and apply a combination of 10 models (naïve Bayes, decision tree and trimmed decision tree, support-vector machine, random forest, 
boosted tree, lasso regression, ridge regression, elastic net, and nearest neighbor) and together with 12 different data preprocessing approaches.
Features were classified as “lipid” or “non-lipid” based on HMDB taxonomy, and the performance of each resulting model was evaluated 
for their accuracy, ROC AUC, and PR AUC for both lipid categories ('lipid”/”non-lipid”). The script 03.figures_and_tables.R generates the figures and tables, while 04.independent_data.R applies the framework to an independent dataset.

# Conclusion
Our result demonstrates that metabolites can be classified as “lipid”, “non-lipid” using only m/z and RT from untargeted LC-MS data, 
without requiring MS2 spectra. Although this study focused on lipid classification, the approach shows potential for broader 
application, which warrants further investigation across diverse compound classes, detection methods, and chromatographic conditions.

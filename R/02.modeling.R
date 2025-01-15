# Import library
library(tidyverse)
library(tidymodels)
library(themis) # for step_smote
library(discrim) # for naive_Bayes
# library(klaR) 


################################################################ Load data----
two_predictor_data <- 
  read_rds(paste0(here::here(), 
                  "/two_predictor_data.rds"))


################################################################ Modeling----
set.seed(1234)

# 1. Splitting the data----
# 3/4 of the data into the training set but split evenly within is_lipids
data_split <- initial_split(two_predictor_data, prop = .75, strata = is_lipids)
# write_rds(data_split, "data_split.rds")
data_split <- 
  read_rds(paste0(here::here(), 
                  "/data_split.rds"))
# Create training and testing data sets:
train_data <- training(data_split)
test_data  <- testing(data_split)


# 2. Prep recipes----
recipe_basic <-
  # 1.model formula
  recipe(is_lipids ~ ., data = train_data)  %>%
  # 2.keep these variables but not use them as either outcomes or predictors
  update_role(hmdb, new_role = "ID")

recipe_smote <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_smote(is_lipids)

recipe_scale <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_scale(all_numeric_predictors())

recipe_norm <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_normalize(all_numeric_predictors())

recipe_corr <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_corr(all_numeric_predictors())

recipe_scale_smote <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_scale(all_numeric_predictors()) %>%
  step_smote(is_lipids)

recipe_norm_smote <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_normalize(all_numeric_predictors()) %>%
  step_smote(is_lipids)

recipe_corr_smote <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_corr(all_numeric_predictors()) %>%
  step_smote(is_lipids)

recipe_scale_corr <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_scale(all_numeric_predictors()) %>%
  step_corr(all_numeric_predictors())

recipe_norm_corr <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_normalize(all_numeric_predictors()) %>%
  step_corr(all_numeric_predictors())

recipe_norm_smote_corr <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_normalize(all_numeric_predictors()) %>%
  step_corr(all_numeric_predictors()) %>%
  step_smote(is_lipids)

recipe_scale_smote_corr <-
  recipe(is_lipids ~ ., data = train_data)  %>%
  update_role(hmdb, new_role = "ID") %>%
  # 3. Set differential steps
  step_scale(all_numeric_predictors()) %>%
  step_corr(all_numeric_predictors()) %>%
  step_smote(is_lipids)

# Generate List of recipes
recipe_list <-
  list(basic = recipe_basic,
       balanced = recipe_smote,
       scaled = recipe_scale, normalized = recipe_norm,
       corr.removed = recipe_corr,
       balanced.scaled = recipe_scale_smote, balanced.normalized = recipe_norm_smote,
       balanced.corr = recipe_corr_smote,
       scaled.corr = recipe_scale_corr, normalized.corr = recipe_norm_corr,
       balanced.normalized.corr = recipe_norm_smote_corr,
       balanced.scaled.corr = recipe_scale_smote_corr)

# 3. Generate spec ----
# tune() is used for model tuning via grid search
nb_spec <-
  naive_Bayes(smoothness = tune(), Laplace = tune()) %>%
  set_engine("klaR") %>%
  set_mode("classification")

dt_spec <- decision_tree(cost_complexity = tune(), tree_depth = tune(), min_n = tune()) %>%
  set_engine("rpart") %>%
  set_mode("classification")

n2_dt_spec <- decision_tree(cost_complexity = tune(), tree_depth = 2, min_n = 2) %>%
  set_engine("rpart") %>%
  set_mode("classification")

svm_spec <-
  svm_rbf(cost = tune(), rbf_sigma = tune(), margin = tune()) %>%
  set_engine("kernlab") %>%
  set_mode("classification")

ranger_spec <- rand_forest(
  mtry = tune(),
  min_n = tune(),
  trees = tune()) %>%
  set_mode("classification") %>%
  set_engine("ranger")

xgboost_spec <-
  boost_tree(trees = tune(), min_n = tune(), tree_depth = tune(), learn_rate = tune(),
             loss_reduction = tune(), sample_size = tune()) %>%
  set_mode("classification") %>%
  set_engine("xgboost")

glmnet_spec <-
  logistic_reg(penalty = tune(), mixture = tune()) %>%
  set_mode("classification") %>%
  set_engine("glmnet")

ridge_spec <-
  logistic_reg(penalty = tune(), mixture = 0) %>%
  set_mode("classification") %>%
  set_engine("glmnet")

lasso_spec <-
  logistic_reg(penalty = tune(), mixture = 1) %>%
  set_mode("classification") %>%
  set_engine("glmnet")

knn_spec <-
  nearest_neighbor(neighbors = tune(), weight_func = tune()) %>%
  set_mode("classification") %>%
  set_engine("kknn")

# 4. Generate list of models----
model_list <-
  list(Random_Forest = ranger_spec, SVM = svm_spec, Naive_Bayes = nb_spec,
       Decision_Tree = dt_spec, n2_Decision_Tree = n2_dt_spec,
       Boosted_Trees = xgboost_spec, KNN = knn_spec, Elasticnet = glmnet_spec,
       Ridge_Regression = ridge_spec, Lasso_Regression = lasso_spec)

model_set <- workflow_set(preproc = recipe_list, models = model_list, cross = T)

# 5. Create cv fold----
set.seed(1234)
mldata_folds <- vfold_cv(train_data, strata = is_lipids)

# 6. Run workflows on the training data set of 1799 using default grid for tuning
# Set the metrics we want the models to be evaluated on
class_metric <- metric_set(accuracy, roc_auc, recall, precision, #pr_auc,
                           sensitivity, specificity)

doParallel::registerDoParallel()
set.seed(1234)
all_workflows <-
  model_set %>% workflow_map(resamples = mldata_folds,
                             metrics = class_metric,
                             verbose = TRUE,
                             control = control_resamples(save_pred = TRUE))


# write_rds(all_workflows, "all_workflows.rds")

# 6. Explore workflows----
all_workflows <- 
  read_rds(paste0(here::here(), 
                  "/all_workflows.rds"))
# Here we can make a plot like figure S1 - code in next R script

# 7. Modeling on training and evaluating on testing with fold----
# First we pick the best spec - here to get best accuracy
set.seed(1234)
final_rf <- all_workflows %>%
  extract_workflow("basic_Random_Forest") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Random_Forest") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_knn <- all_workflows %>%
  extract_workflow("basic_KNN") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_KNN") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_elastinet <- all_workflows %>%
  extract_workflow("basic_Elasticnet") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Elasticnet") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_lasso <- all_workflows %>%
  extract_workflow("basic_Lasso_Regression") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Lasso_Regression") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_ridge <- all_workflows %>%
  extract_workflow("basic_Ridge_Regression") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Ridge_Regression") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_svm <- all_workflows %>%
  extract_workflow("basic_SVM") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_SVM") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_dt <- all_workflows %>%
  extract_workflow("basic_Decision_Tree") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Decision_Tree") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_trimmed_dt <- all_workflows %>%
  extract_workflow("basic_n2_Decision_Tree") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_n2_Decision_Tree") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_boosted_tree <- all_workflows %>%
  extract_workflow("basic_Boosted_Trees") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Boosted_Trees") %>%
                      select_best(metric = "accuracy"))

set.seed(1234)
final_nb <- all_workflows %>%
  extract_workflow("basic_Naive_Bayes") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Naive_Bayes") %>%
                      select_best(metric = "accuracy"))

# Model spec for best accuracy are :
final_rf$fit$actions$model$spec
final_knn$fit$actions$model$spec
final_dt$fit$actions$model$spec
final_trimmed_dt$fit$actions$model$spec
final_elastinet$fit$actions$model$spec
final_lasso$fit$actions$model$spec
final_ridge$fit$actions$model$spec
final_boosted_tree$fit$actions$model$spec
final_svm$fit$actions$model$spec
final_nb$fit$actions$model$spec

# Run models/spec on data fold
set.seed(1234)
rf_results <- final_rf %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
set.seed(1234)
knn_results <- final_knn %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE)
  )
set.seed(1234)
elastinet_results <- final_elastinet %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
set.seed(1234)
lasso_results <- final_lasso %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE)
  )
set.seed(1234)
ridge_results <- final_ridge %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE)
  )
set.seed(1234)
svm_results <- final_svm %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
set.seed(1234)
dt_results <- final_dt %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
set.seed(1234)
trimmed_dt_results <- final_trimmed_dt %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
set.seed(1234)
boosted_tree_results <- final_boosted_tree %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
set.seed(1234)
nb_results <- final_nb %>%
  fit_resamples(
    resamples = mldata_folds,
    metrics = metric_set(roc_auc, accuracy, sensitivity, specificity, precision, recall#, pr_auc
    ),
    control = control_resamples(save_pred = TRUE, verbose = TRUE)
  )
# write_rds(list(rf_results, knn_results, elastinet_results, 
#                lasso_results, ridge_results, svm_results), 
#           "dat.rds")

# a <- read_rds(paste0(here::here(), "/dat.rds"))
# write_rds(rf_results, "rf_results.rds")
# write_rds(knn_results, "knn_results.rds")
# write_rds(elastinet_results, "elastinet_results.rds")
# write_rds(lasso_results, "lasso_results.rds")
# write_rds(ridge_results, "ridge_results.rds")
# write_rds(svm_results, "svm_results.rds")
# write_rds(dt_results, "dt_results.rds")
# write_rds(trimmed_dt_results, "trimmed_dt_results.rds")
# write_rds(boosted_tree_results, "boosted_tree_results.rds")
# write_rds(nb_results, "nb_results.rds")

# Load model results----
rf_results <- read_rds(paste0(here::here(), "/rf_results.rds"))
knn_results <- read_rds(paste0(here::here(), "/knn_results.rds"))
elastinet_results <- read_rds(paste0(here::here(), "/elastinet_results.rds"))
lasso_results <- read_rds(paste0(here::here(), "/lasso_results.rds"))
ridge_results <- read_rds(paste0(here::here(), "/ridge_results.rds"))
svm_results <- read_rds(paste0(here::here(), "/svm_results.rds"))
dt_results <- read_rds(paste0(here::here(), "/dt_results.rds"))
trimmed_dt_results <- read_rds(paste0(here::here(), "/trimmed_dt_results.rds"))
boosted_tree_results <- read_rds(paste0(here::here(), "/boosted_tree_results.rds"))
nb_results <- read_rds(paste0(here::here(), "/nb_results.rds"))


# 8. Fit on the whole training / predict on testing data----
class_metric <- metric_set(accuracy, roc_auc, #recall, precision, #pr_auc, 
                           sensitivity, specificity)
set.seed(1234)
final_rf <- all_workflows %>% 
  extract_workflow("basic_Random_Forest") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_Random_Forest") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_knn <- all_workflows %>% 
  extract_workflow("basic_KNN") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_KNN") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_svm <- all_workflows %>% 
  extract_workflow("basic_SVM") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_SVM") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_eleastinet <- all_workflows %>% 
  extract_workflow("basic_Elasticnet") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_Elasticnet") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_ridge <- all_workflows %>% 
  extract_workflow("basic_Ridge_Regression") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_Ridge_Regression") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_lasso <- all_workflows %>% 
  extract_workflow("basic_Lasso_Regression") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_Lasso_Regression") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_dt <- all_workflows %>% 
  extract_workflow("basic_Decision_Tree") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_Decision_Tree") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_trimmed_dt <- all_workflows %>% 
  extract_workflow("basic_n2_Decision_Tree") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_n2_Decision_Tree") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_boosted_tree <- all_workflows %>% 
  extract_workflow("basic_Boosted_Trees") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_Boosted_Trees") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

set.seed(1234)
final_nb <- all_workflows %>% 
  extract_workflow("basic_Naive_Bayes") %>% 
  finalize_workflow(all_workflows %>% 
                      extract_workflow_set_result("basic_Naive_Bayes") %>% 
                      select_best(metric = "accuracy")) %>% 
  last_fit(data_split,
           metrics = class_metric)

# Now we can explore prediction results on the testing data


################################################################ III. Validation data----
validation3990_data <- 
  read_rds(paste0(here::here(), 
                  "/validation3990_data.rds"))

# Apply our models/spec to the validation data----
final_fitted_RF <- extract_workflow(final_rf)
final_fitted_DT <- extract_workflow(final_dt)
final_fitted_n2DT <- extract_workflow(final_trimmed_dt)
final_fitted_KNN <- extract_workflow(final_knn)

# Get prediction
set.seed(1234)
predicted_validation3990_data_RF <- augment(final_fitted_RF, validation3990_data)
set.seed(1234)
predicted_validation3990_data_DT <- augment(final_fitted_DT, validation3990_data)
set.seed(1234)
predicted_validation3990_data_n2DT <- augment(final_fitted_n2DT, validation3990_data)
set.seed(1234)
predicted_validation3990_data_KNN <- augment(final_fitted_KNN, validation3990_data)



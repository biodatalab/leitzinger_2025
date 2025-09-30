# Import library
library(tidyverse)
library(tidymodels)
library(themis) # for step_smote
library(discrim) # for naive_Bayes
library(klaR)
library(kableExtra)
theme_set(theme_classic())

################################################################ Load data----
validation_ovca <- validation_ovca %>% 
  filter(method == "HILIC-pos")

validation_ovca %>% 
  ggplot(aes(x= row_retention_time, y=row_m_z))+
  geom_point(color= "black")

validation_ovca <- validation_ovca %>% 
  dplyr::select(hmdb, row_m_z, buffer_percent, is_lipids)

################################################################ Modeling----
set.seed(1234)

# 1. Splitting the data----
# 3/4 of the data into the training set but split evenly within is_lipids
data_split <- initial_split(validation_ovca, prop = .75, strata = is_lipids)
# write_rds(data_split, "data_split_ovca.rds")
data_split <- 
  read_rds(paste0(here::here(), 
                  "/data_split_ovca.rds"))
# Create training and testing data sets:
train_data <- training(data_split)
test_data  <- testing(data_split)


# Table training/testing
train_data %>% 
  mutate(data_subset = "OVCA Training set") %>% 
  bind_rows(test_data %>% 
              mutate(data_subset = "OVCA Testing set")) %>% 
  mutate(data_subset = factor(data_subset, 
                              levels = c("OVCA Training set", "OVCA Testing set"))) %>% 
  left_join(., metabolites_classification,
            by = c("hmdb" = "accession")) %>%
  dplyr::select(row_m_z, buffer_percent, is_lipids,
                average_molecular_weight, taxonomy_super_class,
                taxonomy_molecular_framework, taxonomy_sub_class,
                taxonomy_class,
                data_subset) %>% 
  `colnames<-`(str_replace_all(c(colnames(.)), "_", " ")) %>% 
  `colnames<-`(str_to_sentence(colnames(.))) %>% 
  tbl_summary(by = `Data subset`, 
              type = list(`Is lipids` ~ "categorical"),
              sort = list(all_categorical() ~ "frequency"),
              statistic=list(all_continuous() ~ "{median} ({min}, {max})"),
              missing_text = "Unknown") %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_overall() %>% kableExtra::kable()



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


write_rds(all_workflows, "all_workflows_ovca.rds")

# 6. Explore workflows Figure 2----
all_workflows <- 
  read_rds(paste0(here::here(), 
                  "/all_workflows_ovca.rds"))
# Here we can make a plot like figure S1 - code in next R script
collect_metrics(all_workflows) %>%
  separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") %>%
  filter(.metric == "accuracy") %>%
  group_by(wflow_id) %>%
  filter(mean == max(mean)) %>% # equivalent to select_best after tuning
  group_by(model) %>%
  dplyr::select(-.config) %>%
  distinct() %>%
  ungroup() %>%
  mutate(Workflow_Rank =  row_number(-mean),
         .metric = str_to_upper(.metric)) %>%
  
  mutate(Model_Type = factor(Model_Type,
                             levels = c("Decision_Tree", "Random_Forest", 
                                        "KNN", "Elasticnet", "Lasso_Regression", 
                                        "n2_Decision_Tree", 
                                        "Boosted_Trees",
                                        "Ridge_Regression", 
                                        "SVM", "Naive_Bayes"
                             ))) %>%
  filter(!is.na(Model_Type)) %>%
  mutate(Recipe = factor(Recipe,
                         levels = c("basic", "corr.removed", "normalized", "scaled", "balanced",
                                    "normalized.corr", "scaled.corr", "balanced.corr",
                                    "balanced.normalized", "balanced.scaled",
                                    "balanced.normalized.corr", "balanced.scaled.corr"
                         ))) %>% arrange(desc(mean)) %>% 
  distinct(wflow_id, .keep_all = TRUE) %>% 
  # mutate(mn = 0.912) %>% mutate(mx = 0.94) %>%
  ggplot(aes(x=Workflow_Rank, y = mean, shape = Recipe, color = Model_Type)) +
  # geom_rect(aes(xmin = -Inf, xmax = Inf,
  #           ymin = mn, ymax = mx),
  #           fill= "grey90", color= "transparent")+
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean-std_err, ymax = mean+std_err)) +
  # geom_hline(aes(yintercept = 0.926), color= "darkgrey", linetype= 2)+
  # ylim(0.57, 0.95)+
  scale_colour_manual(labels = c("Decision tree", "Random forest", 
                                 "KNN", "Elasticnet", "Lasso regression", 
                                 "Trimmed decision tree", 
                                 "Boosted trees",
                                 "Ridge regression", 
                                 "SVM", "Naive bayes"), 
                      values=c("#8DA0CB", "#FC8D62", "#B2DF8A", "#CAB2D6", "#E5C494",
                               "#66C2A5", "#FFD92F", "#E78AC3", "#A6CEE3", "#B3B3B3"))+
  scale_shape_manual(labels = c("Basic", "Reduced", "Normalized", "Scaled", "Balanced",
                                "Normalized+Reduced", "Scaled+Reduced", "Balanced+Reduced",
                                "Balanced+Normalized", "Balanced+Scaled",
                                "Balanced+Normalized+Reduced", "Balanced+Scaled+Reduced"),
                     values=c(0, 15, 1, 16, 2, 17, 8, 18, 3, 7, 25, 4))+
  labs(#title = "Performance Comparison of Workflow Sets", 
    x = "Workflow Rank",
    y = "Accuracy", color = "Models", shape = "Preprocessors")+
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow = 6, byrow = TRUE))
ggsave("Figure accuracy all models-preprocessors ovca 082425.pdf",
       width = 9,
       height = 9)

collect_metrics(all_workflows) %>%
  separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") %>%
  filter(.metric == "roc_auc") %>% # accuracy   precision      recall     roc_auc sensitivity specificity
  group_by(wflow_id) %>%
  filter(mean == max(mean)) %>% # equivalent to select_best after tuning
  group_by(model) %>%
  dplyr::select(-.config) %>%
  distinct() %>%
  ungroup() %>%
  mutate(Workflow_Rank =  row_number(-mean),
         .metric = str_to_upper(.metric)) %>%
  mutate(Model_Type = factor(Model_Type,
                             levels = c("Decision_Tree", "Random_Forest", 
                                        "KNN", "Elasticnet", "Lasso_Regression", 
                                        "n2_Decision_Tree", 
                                        "Boosted_Trees",
                                        "Ridge_Regression", 
                                        "SVM", "Naive_Bayes"
                             ))) %>%
  filter(!is.na(Model_Type)) %>%
  mutate(Recipe = factor(Recipe,
                         levels = c("basic", "corr.removed", "normalized", "scaled", "balanced",
                                    "normalized.corr", "scaled.corr", "balanced.corr",
                                    "balanced.normalized", "balanced.scaled",
                                    "balanced.normalized.corr", "balanced.scaled.corr"
                         ))) %>% arrange(desc(mean)) %>% 
  ggplot(aes(x=Workflow_Rank, y = mean, shape = Recipe, color = Model_Type)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean-std_err, ymax = mean+std_err)) +
  scale_colour_manual(labels = c("Decision tree", "Random forest", 
                                 "KNN", "Elasticnet", "Lasso regression", 
                                 "Trimmed decision tree", 
                                 "Boosted trees",
                                 "Ridge regression", 
                                 "SVM", "Naive bayes"), 
                      values=c("#8DA0CB", "#FC8D62", "#B2DF8A", "#CAB2D6", "#E5C494",
                               "#66C2A5", "#FFD92F", "#E78AC3", "#A6CEE3", "#B3B3B3"))+
  scale_shape_manual(labels = c("Basic", "Reduced", "Normalized", "Scaled", "Balanced",
                                "Normalized+Reduced", "Scaled+Reduced", "Balanced+Reduced",
                                "Balanced+Normalized", "Balanced+Scaled",
                                "Balanced+Normalized+Reduced", "Balanced+Scaled+Reduced"),
                     values=c(0, 15, 1, 16, 2, 17, 8, 18, 3, 7, 25, 4))+
  # ylim(0.7, 1.0)+
  labs(#title = "Performance Comparison of Workflow Sets", 
    x = "Workflow Rank",
    y = "ROC AUC", color = "Models", shape = "Preprocessors")+
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow = 6, byrow = TRUE))
ggsave("Figure roc_auc all models-preprocessors ovca 082425.pdf",
       width = 9,
       height = 9)


compare_all <- collect_metrics(all_workflows) %>%
  separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") %>%
  group_by(Model_Type, Recipe, .metric) %>% 
  mutate(selected_model = case_when(
    mean == max(mean)            ~ "selected"
  )) %>%
  filter(selected_model == "selected") %>% 
  group_by(Model_Type) %>%
  distinct(pick(!contains(".config")), .keep_all = TRUE) %>%
  ungroup() %>% 
  arrange(Model_Type, Recipe, .metric, mean, std_err) %>% 
  distinct(Model_Type, Recipe, .metric, .keep_all = TRUE) %>% 
  filter(.metric %in% c("roc_auc", "accuracy"))

prauc_data <- collect_predictions(all_workflows, summarize = FALSE) %>%
  separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") %>%
  left_join(compare_all %>%
              dplyr::select(wflow_id, .config, .estimator, Model_Type), .,
            by = c("wflow_id", ".config", "Model_Type"))
prauc_data1 <- prauc_data %>%
  group_by(wflow_id, Recipe, Model_Type, id) %>%
  pr_auc(is_lipids, .pred_No, event_level = "first") %>%
  group_by(wflow_id, Recipe, Model_Type) %>%
  summarize(mean = mean(`.estimate`),
            std_err = sd(`.estimate`)) %>%
  mutate(`.metric` = "pr_auc_No") %>%
  ungroup()
prauc_data2 <- prauc_data %>%
  group_by(wflow_id, Recipe, Model_Type, id) %>%
  pr_auc(is_lipids, .pred_Yes, event_level = "first") %>%
  group_by(wflow_id, Recipe, Model_Type) %>%
  summarize(mean = mean(`.estimate`),
            std_err = sd(`.estimate`)) %>%
  mutate(`.metric` = "pr_auc_Yes") %>%
  ungroup()



prauc_data1 %>%
  mutate(Workflow_Rank =  row_number(-mean),
         .metric = str_to_upper(.metric)) %>%
  
  mutate(Model_Type = factor(Model_Type,
                             levels = c("Decision_Tree", "Random_Forest", 
                                        "KNN", "Elasticnet", "Lasso_Regression", 
                                        "n2_Decision_Tree", 
                                        "Boosted_Trees",
                                        "Ridge_Regression", 
                                        "SVM", "Naive_Bayes"
                             ))) %>%
  filter(!is.na(Model_Type)) %>%
  mutate(Recipe = factor(Recipe,
                         levels = c("basic", "corr.removed", "normalized", "scaled", "balanced",
                                    "normalized.corr", "scaled.corr", "balanced.corr",
                                    "balanced.normalized", "balanced.scaled",
                                    "balanced.normalized.corr", "balanced.scaled.corr"
                         ))) %>% arrange(desc(mean)) %>% 
  ggplot(aes(x=Workflow_Rank, y = mean, shape = Recipe, color = Model_Type)) +
  
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean-std_err, ymax = mean+std_err)) +
  
  scale_colour_manual(labels = c("Decision tree", "Random forest", 
                                 "KNN", "Elasticnet", "Lasso regression", 
                                 "Trimmed decision tree", 
                                 "Boosted trees",
                                 "Ridge regression", 
                                 "SVM", "Naive bayes"), 
                      values=c("#8DA0CB", "#FC8D62", "#B2DF8A", "#CAB2D6", "#E5C494",
                               "#66C2A5", "#FFD92F", "#E78AC3", "#A6CEE3", "#B3B3B3"))+
  scale_shape_manual(labels = c("Basic", "Reduced", "Normalized", "Scaled", "Balanced",
                                "Normalized+Reduced", "Scaled+Reduced", "Balanced+Reduced",
                                "Balanced+Normalized", "Balanced+Scaled",
                                "Balanced+Normalized+Reduced", "Balanced+Scaled+Reduced"),
                     values=c(22, 15, 1, 16, 2, 17, 8, 18, 3, 7, 25, 4))+
  labs(#title = "Performance Comparison of Workflow Sets", 
    x = "",
    y = "PR AUC (non-lipid)", color = "Models", shape = "Preprocessors")+
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow = 6, byrow = TRUE))
ggsave("Figure pr_auc no all models-preprocessors ovca 082425.pdf",
       width = 9,
       height = 9)


prauc_data2 %>%
  mutate(Workflow_Rank =  row_number(-mean),
         .metric = str_to_upper(.metric)) %>%
  
  mutate(Model_Type = factor(Model_Type,
                             levels = c("Decision_Tree", "Random_Forest", 
                                        "KNN", "Elasticnet", "Lasso_Regression", 
                                        "n2_Decision_Tree", 
                                        "Boosted_Trees",
                                        "Ridge_Regression", 
                                        "SVM", "Naive_Bayes"
                             ))) %>%
  filter(!is.na(Model_Type)) %>%
  mutate(Recipe = factor(Recipe,
                         levels = c("basic", "corr.removed", "normalized", "scaled", "balanced",
                                    "normalized.corr", "scaled.corr", "balanced.corr",
                                    "balanced.normalized", "balanced.scaled",
                                    "balanced.normalized.corr", "balanced.scaled.corr"
                         ))) %>% arrange(desc(mean)) %>% 
  ggplot(aes(x=Workflow_Rank, y = mean, shape = Recipe, color = Model_Type)) +
  
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean-std_err, ymax = mean+std_err)) +
  
  scale_colour_manual(labels = c("Decision tree", "Random forest", 
                                 "KNN", "Elasticnet", "Lasso regression", 
                                 "Trimmed decision tree", 
                                 "Boosted trees",
                                 "Ridge regression", 
                                 "SVM", "Naive bayes"), 
                      values=c("#8DA0CB", "#FC8D62", "#B2DF8A", "#CAB2D6", "#E5C494",
                               "#66C2A5", "#FFD92F", "#E78AC3", "#A6CEE3", "#B3B3B3"))+
  scale_shape_manual(labels = c("Basic", "Reduced", "Normalized", "Scaled", "Balanced",
                                "Normalized+Reduced", "Scaled+Reduced", "Balanced+Reduced",
                                "Balanced+Normalized", "Balanced+Scaled",
                                "Balanced+Normalized+Reduced", "Balanced+Scaled+Reduced"),
                     values=c(22, 15, 1, 16, 2, 17, 8, 18, 3, 7, 25, 4))+
  labs(#title = "Performance Comparison of Workflow Sets", 
    x = "",
    y = "PR AUC (lipid)", color = "Models", shape = "Preprocessors")+
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow = 6, byrow = TRUE))
ggsave("Figure pr_auc yes all models-preprocessors ovca 082425.pdf",
       width = 9,
       height = 9)



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

write_rds(rf_results, "rf_results_ovca.rds")
write_rds(knn_results, "knn_results_ovca.rds")
write_rds(elastinet_results, "elastinet_results_ovca.rds")
write_rds(lasso_results, "lasso_results_ovca.rds")
write_rds(ridge_results, "ridge_results_ovca.rds")
write_rds(svm_results, "svm_results_ovca.rds")
write_rds(dt_results, "dt_results_ovca.rds")
write_rds(trimmed_dt_results, "trimmed_dt_results_ovca.rds")
write_rds(boosted_tree_results, "boosted_tree_results_ovca.rds")
write_rds(nb_results, "nb_results_ovca.rds")

# Load model results----
rf_results <- read_rds(paste0(here::here(), "/rf_results_ovca.rds"))
knn_results <- read_rds(paste0(here::here(), "/knn_results_ovca.rds"))
elastinet_results <- read_rds(paste0(here::here(), "/elastinet_results_ovca.rds"))
lasso_results <- read_rds(paste0(here::here(), "/lasso_results_ovca.rds"))
ridge_results <- read_rds(paste0(here::here(), "/ridge_results_ovca.rds"))
svm_results <- read_rds(paste0(here::here(), "/svm_results_ovca.rds"))
dt_results <- read_rds(paste0(here::here(), "/dt_results_ovca.rds"))
trimmed_dt_results <- read_rds(paste0(here::here(), "/trimmed_dt_results_ovca.rds"))
boosted_tree_results <- read_rds(paste0(here::here(), "/boosted_tree_results_ovca.rds"))
nb_results <- read_rds(paste0(here::here(), "/nb_results_ovca.rds"))


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


fig_roc <- final_rf %>%
  collect_predictions() %>%
  mutate(model = "Random forest") %>%
  bind_rows(final_dt %>%
              collect_predictions() %>%
              mutate(model = "Decision tree")) %>%
  bind_rows(., final_trimmed_dt %>%
              collect_predictions() %>%
              mutate(model = "Trimmed decision tree")) %>%
  bind_rows(final_knn %>%
              collect_predictions() %>%
              mutate(model = "KNN")) %>%
  mutate(model = factor(model, levels= c("Decision tree", "Trimmed decision tree", "Random forest","KNN"))) %>% 
  
  group_by(model) %>%
  roc_curve(is_lipids, .pred_No) %>%
  autoplot()+
  scale_color_manual(values = c("#8DA0CB", "#66C2A5", "#FC8D62", "#B2DF8A"), name= NULL)+
  theme_classic()+
  ggtitle("Area under Curve ROC \n(Receiver Operating Characteristics)")

# pr
fig_pr_yes <- final_rf %>%
  collect_predictions() %>%
  mutate(model = "Random forest") %>%
  bind_rows(final_dt %>%
              collect_predictions() %>%
              mutate(model = "Decision tree")) %>%
  bind_rows(., final_trimmed_dt %>%
              collect_predictions() %>%
              mutate(model = "Trimmed decision tree")) %>%
  bind_rows(final_knn %>%
              collect_predictions() %>%
              mutate(model = "KNN")) %>%
  mutate(model = factor(model, levels= c("Decision tree", "Trimmed decision tree", "Random forest","KNN"))) %>% 
  group_by(model) %>%
  pr_curve(is_lipids, .pred_Yes) %>%
  autoplot()+
  scale_color_manual(values = c("#8DA0CB", "#66C2A5", "#FC8D62", "#B2DF8A"), name= NULL)+
  theme_classic()+
  ggtitle("Area under the precision recall curve - \nYes prediction")

fig_pr_no <- final_rf %>%
  collect_predictions() %>%
  mutate(model = "Random forest") %>%
  bind_rows(final_dt %>%
              collect_predictions() %>%
              mutate(model = "Decision tree")) %>%
  bind_rows(., final_trimmed_dt %>%
              collect_predictions() %>%
              mutate(model = "Trimmed decision tree")) %>%
  bind_rows(final_knn %>%
              collect_predictions() %>%
              mutate(model = "KNN")) %>%
  mutate(model = factor(model, levels= c("Decision tree", "Trimmed decision tree", "Random forest","KNN"))) %>% 
  group_by(model) %>%
  pr_curve(is_lipids, .pred_No) %>%
  autoplot()+
  scale_color_manual(values = c("#8DA0CB", "#66C2A5", "#FC8D62", "#B2DF8A"), name= NULL)+
  theme_classic()+
  ggtitle("Area under the precision recall curve - \nNo prediction")

# Figure Scurves testing ----
library(patchwork)
(fig_roc + theme(text = element_text(size = 10),
                 title = element_text(size = 8))) /
  (fig_pr_yes + theme(text = element_text(size = 10),
                      title = element_text(size = 8))) /
  (fig_pr_no + theme(text = element_text(size = 10),
                     title = element_text(size = 8)) )+
  plot_annotation(tag_levels = "A")
ggsave("Figure ROC curves ovca testing.pdf",
       width = 5,
       height = 7)


# Examine over fitting ----
# Figure S5A training and testing -----
testing_performance <- collect_metrics(final_rf) %>% 
  # testing merics
  bind_rows(collect_predictions(final_rf) %>% 
              pr_auc(is_lipids, .pred_No) %>% 
              mutate(.metric = "pr_auc_No")) %>% 
  bind_rows(collect_predictions(final_rf) %>% 
              pr_auc(is_lipids, .pred_Yes) %>% 
              mutate(.metric = "pr_auc_Yes")) %>%
  mutate(model = "Random Forest") %>% 
  
  bind_rows(collect_metrics(final_dt) %>% 
              bind_rows(collect_predictions(final_dt) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_dt) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "Decision Tree")) %>% 
  
  bind_rows(collect_metrics(final_knn) %>% 
              bind_rows(collect_predictions(final_knn) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_knn) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "KNN")) %>% 
  
  bind_rows(collect_metrics(final_eleastinet) %>% 
              bind_rows(collect_predictions(final_eleastinet) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_eleastinet) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "Elasticnet")) %>% 
  
  bind_rows(collect_metrics(final_lasso) %>% 
              bind_rows(collect_predictions(final_lasso) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_lasso) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "Lasso regression")) %>% 
  
  bind_rows(collect_metrics(final_ridge) %>% 
              bind_rows(collect_predictions(final_ridge) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_ridge) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "Ridge regression")) %>% 
  
  bind_rows(collect_metrics(final_trimmed_dt) %>% 
              bind_rows(collect_predictions(final_trimmed_dt) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_trimmed_dt) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "Trimmed Decision Tree")) %>% 
  
  bind_rows(collect_metrics(final_boosted_tree) %>% 
              bind_rows(collect_predictions(final_boosted_tree) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_boosted_tree) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "Boosted Trees")) %>% 
  
  bind_rows(collect_metrics(final_svm) %>% 
              bind_rows(collect_predictions(final_svm) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_svm) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "SVM")) %>% 
  bind_rows(collect_metrics(final_nb) %>% 
              bind_rows(collect_predictions(final_nb) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          mutate(.metric = "pr_auc_No")) %>% 
              bind_rows(collect_predictions(final_nb) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          mutate(.metric = "pr_auc_Yes")) %>%
              mutate(model = "Naive Bayes")) %>% 
  `colnames<-`(c("metric", "estimator", "performance", "config", "model")) %>% 
  mutate(data = "Testing")

training_performance <- collect_metrics(rf_results) %>% 
  `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
  bind_rows(., 
            rf_results %>% # Compare both models
              collect_predictions(summarize = FALSE) %>%
              group_by(id) %>% 
              pr_auc(is_lipids, .pred_No) %>% 
              ungroup() %>% 
              summarize(performance = mean(.estimate),
                        std_err = sd(.estimate)) %>% 
              mutate(metric = "pr_auc_No")) %>% 
  bind_rows(., rf_results %>% # Compare both models
              collect_predictions(summarize = FALSE) %>%
              group_by(id) %>% 
              pr_auc(is_lipids, .pred_Yes) %>% 
              ungroup() %>% 
              summarize(performance = mean(.estimate),
                        std_err = sd(.estimate)) %>% 
              mutate(metric = "pr_auc_Yes")) %>% 
  mutate(model = "Random Forest") %>% 
  
  bind_rows(collect_metrics(dt_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        dt_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., dt_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "Decision Tree")) %>% 
  
  bind_rows(collect_metrics(knn_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        knn_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., knn_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "KNN")) %>% 
  
  bind_rows(collect_metrics(elastinet_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        elastinet_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., elastinet_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "Elasticnet")) %>% 
  
  bind_rows(collect_metrics(lasso_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        lasso_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., lasso_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "Lasso regression")) %>% 
  
  bind_rows(collect_metrics(ridge_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        ridge_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., ridge_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "Ridge regression")) %>% 
  
  bind_rows(collect_metrics(trimmed_dt_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        trimmed_dt_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., trimmed_dt_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "Trimmed Decision Tree")) %>% 
  
  bind_rows(collect_metrics(boosted_tree_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        boosted_tree_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., boosted_tree_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "Boosted Trees")) %>% 
  
  bind_rows(collect_metrics(svm_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        svm_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., svm_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "SVM")) %>% 
  
  bind_rows(collect_metrics(nb_results) %>% 
              `colnames<-`(c("metric", "estimator", "performance", "n", "std_err", "config")) %>% 
              bind_rows(., 
                        nb_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_No) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_No")) %>% 
              bind_rows(., nb_results %>% # Compare both models
                          collect_predictions(summarize = FALSE) %>%
                          group_by(id) %>% 
                          pr_auc(is_lipids, .pred_Yes) %>% 
                          ungroup() %>% 
                          summarize(performance = mean(.estimate),
                                    std_err = sd(.estimate)) %>% 
                          mutate(metric = "pr_auc_Yes")) %>% 
              mutate(model = "Naive Bayes")) %>% 
  
  mutate(data = "Training")

training_performance %>% 
  pivot_wider(id_cols = model,
              names_from = metric, 
              values_from = performance) %>% 
  dplyr::select(-c(precision, recall, sensitivity, specificity)) %>% 
  mutate(rank_prauc_no = rank(-pr_auc_No), .after= pr_auc_No) %>% 
  mutate(rank_prauc_yes = rank(-pr_auc_Yes), .after= pr_auc_Yes) %>%
  mutate(rank_accuracy = rank(-accuracy), .after= accuracy) %>%
  mutate(rank_rocauc = rank(-roc_auc), .after= roc_auc) %>%
  mutate(rank_sum = rowSums(dplyr::select(., starts_with("rank_"))), .after= 1) %>%
  mutate(final_rank = rank(rank_sum), .after= 1) %>%
  arrange(final_rank) %>% 
  kable()

testing_performance %>% 
  mutate(performance = round(performance, 3)) %>% 
  pivot_wider(id_cols = model,
              names_from = metric, 
              values_from = performance) %>% 
  dplyr::select(-c(sensitivity, specificity)) %>%
  mutate(rank_prauc_no = rank(-pr_auc_No), .after= pr_auc_No) %>% 
  mutate(rank_prauc_yes = rank(-pr_auc_Yes), .after= pr_auc_Yes) %>%
  mutate(rank_accuracy = rank(-accuracy), .after= accuracy) %>%
  mutate(rank_rocauc = rank(-roc_auc), .after= roc_auc) %>%
  mutate(rank_sum = rowSums(dplyr::select(., starts_with("rank_"))), .after= 1) %>%
  mutate(final_rank = rank(rank_sum), .after= 1) %>%
  arrange(final_rank) %>% 
  dplyr::select("Models" = model, 
                "Final rank" = final_rank, "Rank sum" = rank_sum, "Accuracy" = accuracy,
                "Accuracy rank" = rank_accuracy, "ROC AUC" = roc_auc, "ROC AUC rank" = rank_rocauc, 
                "PR AUC (No)" = pr_auc_No, "PR AUC (No) rank" = rank_prauc_no, 
                "PR AUC (Yes)" = pr_auc_Yes, "PR AUC (Yes) rank" = rank_prauc_yes) %>% 
  kable()

# Figure S5A + associated Table----
testing_vs_training_performance <- bind_rows(testing_performance, training_performance) %>% 
  mutate(data = factor(data, levels = c("Training", "Testing"))) %>% 
  filter(metric %in% c("accuracy", "roc_auc", "pr_auc_No", "pr_auc_Yes")) %>% 
  mutate(metric = case_when(
    metric == "accuracy"           ~ "Accuracy",
    metric == "roc_auc"            ~ "ROC AUC",
    metric == "pr_auc_No"          ~ "PR AUC (No)",
    metric == "pr_auc_Yes"          ~ "PR AUC (Yes)"
  )) 
testing_vs_training_performance %>%
  mutate(performance = round(performance, 3)) %>% 
  mutate(std_err = round(std_err, 3)) %>% 
  dplyr::select(-c(estimator, config)) %>%
  pivot_wider(names_from = data, 
              values_from = -c(metric, model, data),
              names_glue = "{.value}_{data}") %>% 
  dplyr::select(model, metric, performance_Testing, performance_Training, std_err_Training) %>% 
  kable()
# write_csv(testing_vs_training_performance, "testing_vs_training_performance ovca table.csv")
figs5_panel_A <- testing_vs_training_performance %>% 
  ggplot(aes(x= metric, y=performance, fill = data))+
  geom_bar(stat = "identity",
           position = position_dodge())+
  geom_errorbar(aes(x = metric,
                    ymin = performance - std_err,
                    ymax = performance + std_err),
                width = 0.2, alpha = 0.5,
                position=position_dodge(width=0.9))+
  scale_fill_viridis_d(option = "G", begin = 0.3, end = 0.8, 
                       name = NULL)+
  labs(x = NULL, y = "Performance")+
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))+
  facet_wrap(. ~ model, ncol = 2,
  )
figs5_panel_A

figs5_panel_A+
  theme(legend.position = "bottom")
ggsave("Figure Performance training vs testing ovca2.pdf",
       width = 4,
       height = 8)

# Table ranked testing performance ----
testing_performance <-
  testing_performance %>% 
  dplyr::select(metric, performance, model) %>% 
  pivot_wider(names_from = metric, values_from = performance) %>% 
  mutate(rank_prauc_no = rank(-pr_auc_No), .after= pr_auc_No) %>%
  mutate(rank_prauc_yes = rank(-pr_auc_Yes), .after= pr_auc_Yes) %>%
  mutate(rank_accuracy = rank(-accuracy), .after= accuracy) %>%
  mutate(rank_rocauc = rank(-roc_auc), .after= roc_auc) %>%
  mutate(rank_sum = rowSums(dplyr::select(., starts_with("rank_"))), .after= 1) %>%
  mutate(final_rank = rank(rank_sum), .after= 1) %>%
  dplyr::select(-c(sensitivity, specificity)) %>% 
  arrange(final_rank)
write_csv(testing_performance, "ranked testing_performance ovca.csv")

# 2.Get predictions out on the testing dataset
# Figure S5 panel B Prediction - conf matrix ----
df <- final_rf %>%
  collect_predictions() %>% 
  mutate_at(c("is_lipids"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-lipid"
  )) %>% 
  mutate_at(c(".pred_class"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-\nlipid"
  ))
conf_tab <- table(df$is_lipids, df$.pred_class)
conf_tab <- conf_tab / rowSums(conf_tab)
conf_tab <- as.data.frame(conf_tab, stringsAsFactors = TRUE)
conf_tab$Var2 <- factor(conf_tab$Var2, rev(levels(conf_tab$Var2)))

conf_result_rf_testing <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted", 
       fill = "Overlap")+
  theme_classic()+
  theme(legend.position = "bottom")
conf_result_rf_testing
ggsave("Figure Prediction - conf matrix RF ovca.pdf",
       width = 7,
       height = 6)
final_rf %>%
  collect_predictions() %>%
  conf_mat(is_lipids, .pred_class)

df <- final_boosted_tree %>%
  collect_predictions() %>% 
  mutate_at(c("is_lipids"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-lipid"
  )) %>% 
  mutate_at(c(".pred_class"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-\nlipid"
  ))
conf_tab <- table(df$is_lipids, df$.pred_class)
conf_tab <- conf_tab / rowSums(conf_tab)
conf_tab <- as.data.frame(conf_tab, stringsAsFactors = TRUE)
conf_tab$Var2 <- factor(conf_tab$Var2, rev(levels(conf_tab$Var2)))

conf_result_boosted_tree_testing <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted", 
       fill = "Overlap")+
  theme_classic()+
  theme(legend.position = "bottom")
conf_result_boosted_tree_testing
ggsave("Figure Prediction - conf matrix boosted_tree ovca.pdf",
       width = 7,
       height = 6)

df <- final_nb %>%
  collect_predictions() %>% 
  mutate_at(c("is_lipids"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-lipid"
  )) %>% 
  mutate_at(c(".pred_class"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-\nlipid"
  ))
conf_tab <- table(df$is_lipids, df$.pred_class)
conf_tab <- conf_tab / rowSums(conf_tab)
conf_tab <- as.data.frame(conf_tab, stringsAsFactors = TRUE)
conf_tab$Var2 <- factor(conf_tab$Var2, rev(levels(conf_tab$Var2)))

conf_result_nb_testing <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted", 
       fill = "Overlap")+
  theme_classic()+
  theme(legend.position = "bottom")
conf_result_nb_testing
ggsave("Figure Prediction - conf matrix NB ovca.pdf",
       width = 7,
       height = 6)

df <- final_svm %>%
  collect_predictions() %>% 
  mutate_at(c("is_lipids"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-lipid"
  )) %>% 
  mutate_at(c(".pred_class"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-\nlipid"
  ))
conf_tab <- table(df$is_lipids, df$.pred_class)
conf_tab <- conf_tab / rowSums(conf_tab)
conf_tab <- as.data.frame(conf_tab, stringsAsFactors = TRUE)
conf_tab$Var2 <- factor(conf_tab$Var2, rev(levels(conf_tab$Var2)))

conf_result_smv_testing <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted", 
       fill = "Overlap")+
  theme(legend.position = "bottom")
conf_result_smv_testing
ggsave("Figure Prediction - conf matrix SMV Dec 2024.pdf",
       width = 7,
       height = 6)

legend <- ggpubr::get_legend(conf_result_smv_testing)
# Convert to a ggplot and print
ggpubr::as_ggplot(legend)
ggsave("Figure S5B conf matrix legend ovca.pdf",
       width = 2.4,
       height = 1)

conf_result_rf_testing <- conf_result_rf_testing+
  theme(legend.position = "none")
conf_result_boosted_tree_testing <- conf_result_boosted_tree_testing+
  theme(legend.position = "none")
conf_result_nb_testing <- conf_result_nb_testing+
  theme(legend.position = "none")
conf_result_smv_testing <- conf_result_smv_testing+
  theme(legend.position = "none")


conf_result_rf_testing /
  conf_result_boosted_tree_testing /
  conf_result_nb_testing /
  conf_result_smv_testing +
  plot_layout(axes = "collect")
ggsave("Figure S5B conf matrix ovca.pdf",
       width = 2.2,
       height = 7)




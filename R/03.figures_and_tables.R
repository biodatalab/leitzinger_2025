# Import library
library(tidyverse)
library(tidymodels)
library(gtsummary)
library(kableExtra)
library(klaR)
library(discrim)
theme_set(theme_classic())
theme_gtsummary_compact()

################################################################ Load data----

# Figure S1A----
omics1799_data <- omics1799_data %>% 
  mutate(excluded = case_when(
    row_retention_time > 1 & 
      row_retention_time < 14         ~ "Included",
    TRUE                              ~ "Excluded"
  ))
ini_validation3990_data <- ini_validation3990_data %>% 
  mutate(excluded = case_when(
    row_retention_time > 1 & 
      row_retention_time < 14         ~ "Included",
    TRUE                              ~ "Excluded"
  ))

# figS1_A <-
#   omics1799_data %>% 
#   filter(non_heavy_identified_flag == 0) %>% 
#   ggplot(aes(x= row_retention_time, y=row_m_z, color = excluded))+
#   geom_point(size = 0.5)+
#   scale_color_manual(values = c("red", "grey"))+
#   geom_point(data= omics1799_data %>% filter(non_heavy_identified_flag == 1),
#              aes(x= row_retention_time, y=row_m_z), 
#              color = "black", size = 0.5)+
#   labs(x = "Retention time (mins)", y = "mass/charge number ratio")+
#   xlim(0,16) + ylim(50, 900)+
#   theme(legend.title = element_blank())
# figS1_A
# # ggsave("Figure S1 Panel A training scatter plot.pdf",
# #        width = 7,
# #        height = 6)
# 
# figS1_B <-
#   ini_validation3990_data %>% 
#   filter(non_heavy_identified_flag == 0) %>% 
#   ggplot(aes(x= row_retention_time, y=row_m_z, color = excluded))+
#   geom_point(size = 0.5)+
#   scale_color_manual(values = c("red", "grey"))+
#   geom_point(data= ini_validation3990_data %>% filter(non_heavy_identified_flag == 1),
#              aes(x= row_retention_time, y=row_m_z), 
#              color = "black", size = 0.5)+
#   labs(x = "Retention time (mins)", y = "mass/charge number ratio")+
#   xlim(0,16) + ylim(50, 900)+
#   theme(legend.title = element_blank())
# figS1_B
# # ggsave("Figure S1 Panel B validation scatter plot.pdf",
# #        width = 7,
# #        height = 6)
# figS1_B <-
#   ini_validation_ovca %>% 
#   filter(non_heavy_identified_flag == 0) %>% 
#   ggplot(aes(x= row_retention_time, y=row_m_z, color = excluded))+
#   geom_point(size = 0.5)+
#   scale_color_manual(values = c("red", "grey"))+
#   geom_point(data= ini_validation_ovca %>% filter(non_heavy_identified_flag == 1),
#              aes(x= row_retention_time, y=row_m_z), 
#              color = "black", size = 0.5)+
#   labs(x = "Retention time (mins)", y = "mass/charge number ratio")+
#   xlim(0,16) + ylim(50, 900)+
#   theme(legend.title = element_blank())
# figS1_B
# ggsave("Figure ini_validation_ovca scatter plot.pdf",
#        width = 7,
#        height = 6)

# Fig S1A ----
non_identified <- omics1799_data %>% 
  filter(non_heavy_identified_flag == 0) %>% 
  mutate(data = "Discovery") %>% 
  bind_rows(ini_validation3990_data %>% 
              filter(non_heavy_identified_flag == 0) %>% 
              mutate(data = "Validation"))
identified <- omics1799_data %>% 
  filter(non_heavy_identified_flag == 1) %>% 
  mutate(data = "Discovery") %>% 
  bind_rows(ini_validation3990_data %>% 
              filter(non_heavy_identified_flag == 1) %>% 
              mutate(data = "Validation"))

colors <- c("Identified" = "black", "Excluded" = "red", "Included" = "grey")
figS1_A <- non_identified %>% 
  bind_rows(identified %>% 
              mutate(metalolites = "identified")) %>% 
  ggplot(aes(x= row_retention_time, y=row_m_z, color = excluded))+
  geom_point(size = 0.5)+
  geom_point(data= identified,
             aes(x= row_retention_time, y=row_m_z,
                 color = "Identified"),
             size = 0.5)+
  scale_color_manual(values = colors,
                     breaks=c("Identified","Included","Excluded"))+
  labs(x = "Retention time (mins)", y = "Mass/Charge ratio")+
  xlim(0,16) + ylim(50, 900)+
  facet_wrap(. ~ data)+
  theme(legend.title = element_blank())
figS1_A
ggsave("Figure 1 Panel A scatter plot.pdf",
       width = 11,
       height = 6)


################################################################ Modeling----
data_split <- 
  read_rds(paste0(here::here(), 
                  "/data_split.rds"))
# Create training and testing datasets:
train_data <- training(data_split)
test_data  <- testing(data_split)

# Data subset description

# Table S2----
tbl1 <- train_data %>% 
  mutate(data_subset = "Training set") %>% 
  bind_rows(test_data %>% 
              mutate(data_subset = "Testing set")) %>% 
  mutate(data_subset = factor(data_subset, 
                              levels = c("Training set", "Testing set"))) %>% 
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
              type = `Is lipids` ~ "categorical",
              sort = everything() ~ "frequency",
              statistic=list(all_continuous() ~ "{median} ({min}, {max})"),
              missing_text = "Unknown") %>% 
  bold_labels() %>% add_stat_label() %>% 
  add_overall()

tbl2 <- validation3990_data %>% 
  mutate(data_subset = "Validation set") %>%
  dplyr::select(row_m_z, buffer_percent, is_lipids,
                average_molecular_weight, taxonomy_super_class,
                taxonomy_molecular_framework, taxonomy_sub_class,
                taxonomy_class,
                data_subset) %>% 
  `colnames<-`(str_replace_all(c(colnames(.)), "_", " ")) %>% 
  `colnames<-`(str_to_sentence(colnames(.))) %>% 
  tbl_summary(by = `Data subset`, 
              type = `Is lipids` ~ "categorical",
              sort = everything() ~ "frequency",
              statistic=list(all_continuous() ~ "{median} ({min}, {max})"),
              missing_text = "Unknown") %>% 
  bold_labels() %>% add_stat_label()

tbl_merge(list(tbl1, tbl2), tab_spanner = c("**Training dataset**", 
                                            "**Validation 1 datset**")) %>% 
  kableExtra::kable()


# Figure 1B----
a <- omics1799_data %>% 
  mutate(data = "Initial dataset") %>% 
  bind_rows(full_data %>%
              mutate(data = "Cleaned dataset")) %>% 
  mutate(data = factor(data, 
                       levels = c("Initial dataset",
                                  "Cleaned dataset"))) %>% 
  separate(row_id, into = c("charge", "id")) %>%
  mutate(charge = case_when(
    charge == "neg"            ~ "Negative",
    charge == "pos"            ~ "Positive"
  )) %>% 
  group_by(data) %>% 
  mutate(total = n(), .after = 1) %>% 
  group_by(data, charge) %>% 
  mutate(sum = n(), .after = 1,
         perc = sum / total * 100) %>% 
  distinct(data, charge, .keep_all = TRUE)
b <- ini_validation3990_data %>% 
  mutate(data = "Initial dataset") %>% 
  bind_rows(validation3990_data %>%
              mutate(data = "Cleaned dataset")) %>% 
  mutate(data = factor(data, 
                       levels = c("Initial dataset",
                                  "Cleaned dataset"))) %>% 
  separate(row_id, into = c("charge", "id")) %>%
  mutate(charge = case_when(
    charge == "neg"            ~ "Negative",
    charge == "pos"            ~ "Positive"
  )) %>% 
  group_by(data) %>% 
  mutate(total = n(), .after = 1) %>% 
  group_by(data, charge) %>% 
  mutate(sum = n(), .after = 1,
         perc = sum / total * 100) %>% 
  distinct(data, charge, .keep_all = TRUE)

figS1_B <- a %>% mutate(dat = "Discovery") %>% 
  bind_rows(., b %>% mutate(dat = "Validation 1")) %>% 
  ggplot(aes(x=data, y = perc, fill= charge))+
  labs(x = NULL, y = "Metabolites (%)")+
  scale_fill_viridis_d(option = "H", begin = 0.05, end = 0.75, 
                       direction = -1)+
  geom_bar(stat="identity")+
  theme(#legend.position = "bottom",
    legend.title = element_blank(), 
    axis.text.x = element_text(angle = 45, 
                               hjust = 1),
    axis.title.x = element_text(vjust = 10))+
  facet_wrap(. ~ dat)
figS1_B
ggsave("Figure 1 Panel B metabolite charge content.pdf",
       width = 8,
       height = 6)


# Figure S1C----
validation3990_data <- 
  read_rds(paste0(here::here(), 
                  "/validation3990_data_with_nonidentified_metabolites.rds"))

figS1_C <- full_data %>% 
  dplyr::select(non_heavy_identified_flag) %>% 
  mutate(dat = "Discovery") %>% 
  bind_rows(., validation3990_data %>% 
              dplyr::select(non_heavy_identified_flag) %>% 
              mutate(dat = "Validation")) %>% 
  mutate(non_heavy_identified_flag = case_when(
    non_heavy_identified_flag == 1            ~ "Identified",
    TRUE                                      ~ "Non identified"
  )) %>% 
  mutate(non_heavy_identified_flag = factor(non_heavy_identified_flag, 
                                            levels = c("Non identified",
                                                       "Identified"))) %>% 
  group_by(dat) %>% 
  mutate(total = n(), .after = 1) %>% 
  group_by(dat, non_heavy_identified_flag) %>% 
  mutate(sum = n(), .after = 1,
         perc = sum / total * 100) %>% 
  distinct(dat, non_heavy_identified_flag, .keep_all = TRUE) %>% 
  ggplot(aes(x=dat, y = perc, fill= non_heavy_identified_flag))+
  labs(x = NULL, y = "Metabolites (%)")+
  scale_fill_viridis_d(option = "F", begin = 0.2,
                       end = 0.95, direction = -1,
                       name = NULL)+
  geom_bar(stat="identity")+
  theme(#legend.position = "bottom",
    legend.title = element_blank(), 
    axis.text.x = element_text(angle = 45, 
                               hjust = 1),
    axis.title.x = element_text(vjust = 5))
figS1_C
ggsave("Figure 1 Panel C identified metabolites content.pdf",
       width = 8,
       height = 6)
  
  
a <- full_data %>% 
  dplyr::select(non_heavy_identified_flag) %>% 
  mutate(non_heavy_identified_flag = case_when(
    non_heavy_identified_flag == 1            ~ "identified",
    TRUE                                      ~ "non identified"
  )) %>% 
  mutate(non_heavy_identified_flag = factor(non_heavy_identified_flag, 
                                            levels = c("non identified",
                                                       "identified"))) %>% 
  group_by(non_heavy_identified_flag) %>% 
  mutate(ypos = n()) %>% 
  mutate(ypos = ypos - 0.5*ypos) %>% 
  ungroup() %>% 
  mutate(total = n(), .after = 1) %>% 
  group_by(non_heavy_identified_flag) %>% 
  mutate(sum = n(), .after = 1,
         perc = sum / total * 100)
b <- validation3990_data %>% 
  dplyr::select(non_heavy_identified_flag) %>% 
  mutate(non_heavy_identified_flag = case_when(
    non_heavy_identified_flag == 1            ~ "identified",
    TRUE                                      ~ "non identified"
  )) %>% 
  mutate(non_heavy_identified_flag = factor(non_heavy_identified_flag, 
                                            levels = c("non identified",
                                                       "identified"))) %>% 
  group_by(non_heavy_identified_flag) %>% 
  mutate(ypos = n()) %>% 
  mutate(ypos = ypos - 0.5*ypos) %>% 
  ungroup() %>% 
  mutate(total = n(), .after = 1) %>% 
  group_by(non_heavy_identified_flag) %>% 
  mutate(sum = n(), .after = 1,
         perc = sum / total * 100)

legend <-
  a %>% 
  ggplot(aes(x= "", fill= non_heavy_identified_flag))+
  geom_bar(width=1, color= "white") +
  coord_polar("y", start=0)+
  theme_void()+
  scale_fill_viridis_d(option = "A", begin = 0.3, end = 0.8, 
                       name = NULL)+
  theme(legend.position = "bottom")
legend <- ggpubr::get_legend(legend)
# Convert to a ggplot and print
ggpubr::as_ggplot(legend)
# ggsave("Figure 1C legend.pdf", # Not used
#        width = 6,
#        height = 3)


# Figure 1D----
validation3990_data <- 
  read_rds(paste0(here::here(), 
                  "/validation3990_data_with_nonidentified_metabolites.rds"))
two_predictor_data <- 
  read_rds(paste0(here::here(), 
                  "/two_predictor_data.rds"))

figS1_D <-
  two_predictor_data %>% 
  mutate(dat = "Discovery") %>% 
  bind_rows(., validation3990_data %>%
              mutate(dat = "Validation"))%>%
  mutate(is_lipids = case_when(
    is_lipids == "Yes"      ~ "Lipid",
    is_lipids == "No"       ~ "Non-lipid"
  )) %>% 
  group_by(dat) %>% 
  mutate(total = n(), .after = 1) %>% 
  group_by(dat, is_lipids) %>% 
  mutate(sum = n(), .after = 1,
         perc = sum / total * 100) %>% 
  distinct(dat, is_lipids, .keep_all = TRUE) %>% 
  ggplot(aes(x=is_lipids, y= perc, fill= is_lipids))+
  geom_bar(stat = "identity")+
  scale_fill_viridis_d(option = "H", 
                       begin = 0.2,
                       end = 0.8, #direction = -1,
                       name = NULL)+
  labs(x = NULL, y = "Metabolites (%)")+
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1))+
  facet_wrap(. ~ dat)
figS1_D
ggsave("Figure 1 Panel E Lipids count.pdf",
       width = 8,
       height = 6)


# Figure 1E----
figS1_E <-
  training(data_split) %>% 
  mutate(data = "Training") %>% 
  bind_rows(testing(data_split) %>% 
              mutate(data = "Testing")) %>% 
  mutate(is_lipids = case_when(
    is_lipids == "Yes"      ~ "Lipid",
    is_lipids == "No"       ~ "Non-lipid"
  )) %>% 
  # bind_rows(validation3990_data %>% 
  #             mutate(data = "Validation")) %>% 
  # mutate(data = factor(data, levels = c("Validation", "Training", "Testing"))) %>% 
  ggplot(aes(x = buffer_percent, y = row_m_z, color= data#, shape = data
  ))+
  geom_point(size = 0.5)+
  labs(x = "Gradient Composition", y = "Mass/Charge Ratio")+
  # scale_color_discrete(name = NULL)+
  # scale_shape_manual(values=c(3, 16,16), name = NULL)+
  # scale_color_viridis_d(option = "G", begin = 0.3, end = 0.8, 
  #                       name = NULL)+
  scale_color_manual(values = c(#"grey", 
    "#414081FF", "#60CEACFF"), name = NULL)+
  theme(axis.title.x = element_text(vjust = 8))+
  facet_wrap(. ~ is_lipids, ncol = 2)
figS1_E
# ggsave("Figure 1 Panel E ML data.pdf",
#        width = 7,
#        height = 6)


# 1. Visualize performance comparison of combined workflows - Explore Tuning param
# Figure S2 Accuracy, ROC AUC, PR AUC all models-preprocessors----
all_workflows <- 
  read_rds(paste0(here::here(), 
                  "/all_workflows.rds"))

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
  labs(#title = "Performance Comparison of Workflow Sets", 
    x = "Workflow Rank",
    y = "Accuracy", color = "Models", shape = "Preprocessors")+
  theme(legend.position = "bottom", legend.box = "vertical") +
  guides(fill=guide_legend(nrow = 6, byrow = TRUE))
ggsave("Figure S2 Accuracy all models-preprocessors Dec2024.pdf",
       width = 8,
       height = 7)

collect_metrics(all_workflows) %>%
  separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") %>%
  filter(.metric == "roc_auc") %>% 
  group_by(wflow_id) %>%
  filter(mean == max(mean)) %>% 
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
  ylim(0.7, 1.0)+
  labs(#title = "Performance Comparison of Workflow Sets", 
    x = "Workflow Rank",
    y = "ROC AUC", color = "Models", shape = "Preprocessors")+
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow = 6, byrow = TRUE))
ggsave("Figure roc_auc all models-preprocessors 082425.pdf",
       width = 9,
       height = 9)

# For PR AUC
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
  ylim(0.7, 1)+
  labs(#title = "Performance Comparison of Workflow Sets", 
    x = "",
    y = "PR AUC (non-lipid)", color = "Models", shape = "Preprocessors")+
  theme(legend.position = "bottom", legend.box = "vertical",
        axis.text.x = element_blank()) +
  guides(fill=guide_legend(nrow = 6, byrow = TRUE))
ggsave("Figure pr_auc no all models-preprocessors 082425.pdf",
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
ggsave("Figure pr_auc yes all models-preprocessors 082425.pdf",
       width = 9,
       height = 9)
  
# Not run
# compare_all <- collect_metrics(all_workflows) %>%
#   separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") %>% 
#   filter(Recipe == "basic") %>% 
#   left_join(compare_all %>% 
#               dplyr::select(wflow_id, .config, Model_Type), .,
#             by = c("wflow_id", ".config", "Model_Type")) %>% 
#   bind_rows(prauc_data1, prauc_data2, .)
# rm(prauc_data, prauc_data1, prauc_data2)
# 
# 
# metrics_data <- compare_all %>% 
#   dplyr::select(wflow_id, .metric, mean) %>% 
#   mutate(mean = round(mean, 3)) %>% 
#   # group_by(wflow_id) %>% 
#   pivot_wider(id_cols = wflow_id, 
#               names_from = .metric, 
#               values_from = mean) %>% 
#   mutate(id = c(4, 1, 5, 8, 7, 10, 3, 6, 9, 2)) %>% 
#   arrange(id) %>% 
#   dplyr::select(wflow_id, pr_auc_No, pr_auc_Yes, accuracy, sensitivity,
#                 specificity, precision, recall, roc_auc)
# 
# library(kableExtra)
# library(formattable)
# metrics_data %>% 
#   mutate(pr_auc_No = ifelse(pr_auc_No > 0.93,
#                             cell_spec(pr_auc_No, background = "lightgreen"),
#                             cell_spec(pr_auc_No)
#   ),
#   pr_auc_Yes = ifelse(pr_auc_Yes > 0.61,
#                       cell_spec(pr_auc_Yes, background = "lightgreen"),
#                       cell_spec(pr_auc_Yes)
#   ),
#   accuracy = ifelse(accuracy > 0.87,
#                     cell_spec(accuracy, background = "lightgreen"),
#                     cell_spec(accuracy)
#   ),
#   sensitivity = ifelse(sensitivity > 0.95,
#                        cell_spec(sensitivity, background = "lightgreen"),
#                        cell_spec(sensitivity)
#   ),
#   specificity = ifelse(specificity > 0.6,
#                        cell_spec(specificity, background = "lightgreen"),
#                        cell_spec(specificity)
#   )
#   ) %>% 
#   kableExtra::kable(escape = F) %>%
#   kableExtra::kable_classic(full_width = T, html_font = "Cambria", font_size = 16) %>%
#   kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
# 
# metrics_data <- metrics_data %>% 
#   mutate(rank_prauc_no = rank(-pr_auc_No), .after= pr_auc_No) %>% 
#   mutate(rank_prauc_yes = rank(-pr_auc_Yes), .after= pr_auc_Yes) %>% 
#   mutate(rank_accuracy = rank(-accuracy), .after= accuracy) %>% 
#   mutate(rank_rocauc = rank(-roc_auc), .after= roc_auc) %>% 
#   mutate(rank_sum = rowSums(dplyr::select(., starts_with("rank_"))), .after= 1) %>% 
#   mutate(final_rank = rank(rank_sum), .after= 1) %>% 
#   arrange(final_rank)
# write_rds(metrics_data, "metrics_data.rds")
# write_csv(metrics_data, "metrics_data.csv")
library(kableExtra)
metrics_data <- read_rds(paste0(here::here(), "/metrics_data.rds"))
metrics_data %>% 
  mutate(wflow_id = str_remove(wflow_id, "basic_")) %>% 
  dplyr::select("Models" = wflow_id, 
                "Accuracy" = accuracy,
                "ROC AUC" = roc_auc, 
                "PR AUC (No)" = pr_auc_No, 
                "PR AUC (Yes)" = pr_auc_Yes) %>% 
  kable()

two_predictor_perfmean_RF <- 
  metrics_data %>% 
  filter(wflow_id == "basic_Random_Forest") %>% 
  dplyr::select(accuracy) %>% as.numeric()



# Load model results
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


# Examine performance on the training set----
# accuracy, roc_auc, sensitivity, specificity, precision, recall
trimmed_dt_results %>% # Compare both models
  collect_metrics() %>%
  dplyr::select(".metric", n2_tree_mean = mean) %>%
  full_join(., rf_results %>%
              collect_metrics() %>%
              dplyr::select(".metric", rf_mean = mean), by = ".metric") %>%
  full_join(., dt_results %>% 
              collect_metrics() %>%
              dplyr::select(".metric", tree_mean = mean), by = ".metric") %>% 
  full_join(., knn_results %>%
              collect_metrics() %>%
              dplyr::select(".metric", knn_mean = mean), by = ".metric") %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))

# roc
fig_roc <- trimmed_dt_results %>% # Compare both models
  collect_predictions() %>%
  mutate(model = "Trimmed decision tree") %>%
  bind_rows(rf_results %>%
              collect_predictions() %>%
              mutate(model = "Random forest")) %>%
  bind_rows(., dt_results %>% 
              collect_predictions() %>%
              mutate(model = "Decision tree")) %>%
  bind_rows(knn_results %>%
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
fig_pr_yes <- trimmed_dt_results %>% # Compare both models
  collect_predictions() %>%
  mutate(model = "Trimmed decision tree") %>%
  bind_rows(rf_results %>%
              collect_predictions() %>%
              mutate(model = "Random forest")) %>%
  bind_rows(., dt_results %>% 
              collect_predictions() %>%
              mutate(model = "Decision tree")) %>%
  bind_rows(knn_results %>%
              collect_predictions() %>%
              mutate(model = "KNN")) %>%
  mutate(model = factor(model, levels= c("Decision tree", "Trimmed decision tree", "Random forest","KNN"))) %>% 
  group_by(model) %>%
  pr_curve(is_lipids, .pred_Yes) %>%
  autoplot()+
  scale_color_manual(values = c("#8DA0CB", "#66C2A5", "#FC8D62", "#B2DF8A"), name= NULL)+
  theme_classic()+
  ggtitle("Area under the precision recall curve - \nYes prediction")
trimmed_dt_results %>% # Compare both models
  collect_predictions() %>%
  mutate(model = "n2_tree") %>%
  bind_rows(rf_results %>%
              collect_predictions() %>%
              mutate(model = "rf")) %>%
  bind_rows(., dt_results %>% 
              collect_predictions() %>%
              mutate(model = "tree")) %>%
  bind_rows(knn_results %>%
              collect_predictions() %>%
              mutate(model = "knn")) %>%
  group_by(model) %>%
  pr_auc(is_lipids, .pred_Yes) %>% 
  mutate(.metric = "pr_auc_Yes") %>% 
  arrange(-.estimate) %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))

fig_pr_no <- trimmed_dt_results %>% # Compare both models
  collect_predictions() %>%
  mutate(model = "Trimmed decision tree") %>%
  bind_rows(rf_results %>%
              collect_predictions() %>%
              mutate(model = "Random forest")) %>%
  bind_rows(., dt_results %>% 
              collect_predictions() %>%
              mutate(model = "Decision tree")) %>%
  bind_rows(knn_results %>%
              collect_predictions() %>%
              mutate(model = "knn")) %>%
  mutate(model = factor(model, levels= c("Decision tree", "Trimmed decision tree", "Random forest","KNN"))) %>% 
  group_by(model) %>%
  pr_curve(is_lipids, .pred_No) %>%
  autoplot()+
  scale_color_manual(values = c("#8DA0CB", "#66C2A5", "#FC8D62", "#B2DF8A"), name= NULL)+
  theme_classic()+
  ggtitle("Area under the precision recall curve - \nNo prediction")
trimmed_dt_results %>% # Compare both models
  collect_predictions() %>%
  mutate(model = "n2_tree") %>%
  bind_rows(rf_results %>%
              collect_predictions() %>%
              mutate(model = "rf")) %>%
  bind_rows(., dt_results %>% 
              collect_predictions() %>%
              mutate(model = "tree")) %>%
  bind_rows(knn_results %>%
              collect_predictions() %>%
              mutate(model = "knn")) %>%
  group_by(model) %>%
  pr_auc(is_lipids, .pred_No) %>% 
  mutate(.metric = "pr_auc_No") %>% 
  arrange(-.estimate) %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))


# 5. Fit on the whole training / predict on testing data (1799 data)----
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

# Examine performance
# accuracy, roc_auc, sensitivity, specificity, precision, recall
# final_performance <- final_rf %>% # Compare both models
#   collect_metrics() %>%
#   dplyr::select(".metric", rf = .estimate) %>%
#   full_join(., final_knn %>%
#               collect_metrics() %>%
#               dplyr::select(".metric", knn = .estimate), by = ".metric") %>%
#   full_join(., final_eleastinet %>% 
#               collect_metrics() %>%
#               dplyr::select(".metric", elasticnet = .estimate), by = ".metric") %>% 
#   full_join(., final_lasso %>%
#               collect_metrics() %>%
#               dplyr::select(".metric", lasso = .estimate), by = ".metric") %>% 
#   full_join(., final_ridge %>%
#               collect_metrics() %>%
#               dplyr::select(".metric", ridge = .estimate), by = ".metric") %>% 
#   full_join(., final_svm %>%
#               collect_metrics() %>%
#               dplyr::select(".metric", svm = .estimate), by = ".metric") %>% 
#   mutate(across(where(is.numeric), ~ round(., 3)))
# final_performance
# 
# final_rf %>%
#   collect_metrics() %>%
#   mutate(model = "rf") %>%
#   bind_rows(final_knn %>%
#               collect_metrics() %>%
#               mutate(model = "knn")) %>%
#   bind_rows(., final_eleastinet %>% 
#               collect_metrics() %>%
#               mutate(model = "elasticnet")) %>%
#   bind_rows(final_lasso %>%
#               collect_metrics() %>%
#               mutate(model = "lasso")) %>%
#   bind_rows(final_ridge %>%
#               collect_metrics() %>%
#               mutate(model = "ridge")) %>%
#   bind_rows(final_svm %>%
#               collect_metrics() %>%
#               mutate(model = "svm")) %>%
#   mutate(model = factor(model, levels= c("rf", "knn","elasticnet", "lasso", "ridge", "svm"))) %>%
#   
#   ggplot(aes(x= .metric, y= .estimate, fill = model))+
#   geom_bar(stat = "identity",
#            position = position_dodge())+
#   # geom_errorbar(aes(x = .metric,
#   #                   ymin = mean - std_err,
#   #                   ymax = mean + std_err),
#   #               width = 0.2, alpha = 0.5,
#   #               position=position_dodge(width=0.9))+
#   scale_fill_manual(values = c("#8DA0CB", "#FC8D62", "#B2DF8A", "#CAB2D6", "#E5C494",
#                                "#66C2A5"), 
#                     name= NULL)
# 

# Figure S3
# roc
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
fig_pr_yes <-
  final_rf %>%
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
  geom_hline(yintercept = 0.23,  linetype = "dotted")+
  scale_color_manual(values = c("#8DA0CB", "#66C2A5", "#FC8D62", "#B2DF8A"), name= NULL)+
  theme_classic()+
  ggtitle("Area under the precision recall curve - \n'lipid' prediction")

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
  geom_hline(yintercept = 0.77,  linetype = "dotted")+
  scale_color_manual(values = c("#8DA0CB", "#66C2A5", "#FC8D62", "#B2DF8A"), name= NULL)+
  theme_classic()+
  ggtitle("Area under the precision recall curve - \n'non-lipid' prediction")

# Figure S3 curves testing
library(patchwork)
(fig_roc + theme(text = element_text(size = 10),
                 title = element_text(size = 8))) /
  (fig_pr_yes + theme(text = element_text(size = 10),
                      title = element_text(size = 8))) /
  (fig_pr_no + theme(text = element_text(size = 10),
                     title = element_text(size = 8)) )+
  plot_annotation(tag_levels = "A")
ggsave("Figure S3 curves testing 082525.pdf",
       width = 5,
       height = 7)


# Examine over fitting ----
# Figure 3 training and testing -----
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

# Table 1----
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

# Figure 3 + associated Table S4----
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
# write_csv(testing_vs_training_performance, "testing_vs_training_performance table.csv")
fig1_panelF <- testing_vs_training_performance %>% 
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
fig1_panelF
# ggsave("Figure 3 Performance training vs testing.pdf",
#        width = 5,
#        height = 7)


# FIGURE 1------------------

design <- "
  1155
2255
3455
"
(figS1_A + theme(legend.box.margin=margin(0,0,0,-20)) +
    (figS1_B+ figS1_C + plot_layout(guides = "collect"))+
    figS1_D + figS1_E+ theme(legend.box.margin=margin(0,0,0,-20)))+
  fig1_panelF + theme(legend.box.margin=margin(0,0,0,-10)) +
  plot_annotation(tag_levels = "A")+
  plot_layout(design = design)

ggsave("Figure 1 New landscape Dec.pdf", # pick this one
       width = 12,
       height = 10)

design <- "
##55
1155
  1155
2255
2255
3455
3455
##55
"
figS1_A +
  (figS1_B + figS1_C) + 
  figS1_D + figS1_E +
  fig1_panelF  + plot_layout(guides = "collect")+
  # plot_annotation(tag_levels = "A")+ 
  plot_layout(design = design)

ggsave("Figure 1 New landscape guide collect Dec.pdf",
       width = 11,
       height = 10)

design <- "
1166
2266
3466
5566
"

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
write_csv(testing_performance, "ranked testing_performance.csv")


################################################################ III. Validation data----
validation3990_data <- 
  read_rds(paste0(here::here(), 
                  "/validation3990_data_with_nonidentified_metabolites.rds"))

# Apply our model to the validation data----
final_fitted_RF <- extract_workflow(final_rf)
final_fitted_DT <- extract_workflow(final_dt)
final_fitted_n2DT <- extract_workflow(final_trimmed_dt)
final_fitted_KNN <- extract_workflow(final_knn)

set.seed(1234)
predicted_validation3990_data_RF <- augment(final_fitted_RF, validation3990_data)
set.seed(1234)
predicted_validation3990_data_DT <- augment(final_fitted_DT, validation3990_data)
set.seed(1234)
predicted_validation3990_data_n2DT <- augment(final_fitted_n2DT, validation3990_data)
set.seed(1234)
predicted_validation3990_data_KNN <- augment(final_fitted_KNN, validation3990_data)


# Performance on validation 3990, table 2----
## RF
perf_validation_rf <- predicted_validation3990_data_RF %>% 
  accuracy(truth = is_lipids, estimate = .pred_class) %>% 
  bind_rows(predicted_validation3990_data_RF %>%
              pr_auc(is_lipids, .pred_No) %>%
              mutate(.metric = "pr_auc_No")) %>%
  bind_rows(predicted_validation3990_data_RF %>% 
              pr_auc(is_lipids, .pred_Yes) %>%
              mutate(.metric = "pr_auc_Yes")) %>% 
  bind_rows(predicted_validation3990_data_RF %>% 
              roc_auc(is_lipids, .pred_No)) %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))
perf_validation_dt <- predicted_validation3990_data_DT %>% 
  accuracy(truth = is_lipids, estimate = .pred_class) %>% 
  bind_rows(predicted_validation3990_data_DT %>%
              pr_auc(is_lipids, .pred_No) %>%
              mutate(.metric = "pr_auc_No")) %>%
  bind_rows(predicted_validation3990_data_DT %>% 
              pr_auc(is_lipids, .pred_Yes) %>%
              mutate(.metric = "pr_auc_Yes")) %>% 
  bind_rows(predicted_validation3990_data_DT %>% 
              roc_auc(is_lipids, .pred_No)) %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))
perf_validation_n2dt <- predicted_validation3990_data_n2DT %>% 
  accuracy(truth = is_lipids, estimate = .pred_class) %>% 
  bind_rows(predicted_validation3990_data_n2DT %>%
              pr_auc(is_lipids, .pred_No) %>%
              mutate(.metric = "pr_auc_No")) %>%
  bind_rows(predicted_validation3990_data_n2DT %>% 
              pr_auc(is_lipids, .pred_Yes) %>%
              mutate(.metric = "pr_auc_Yes")) %>% 
  bind_rows(predicted_validation3990_data_n2DT %>% 
              roc_auc(is_lipids, .pred_No)) %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))
perf_validation_knn <- predicted_validation3990_data_KNN %>% 
  accuracy(truth = is_lipids, estimate = .pred_class) %>% 
  bind_rows(predicted_validation3990_data_KNN %>%
              pr_auc(is_lipids, .pred_No) %>%
              mutate(.metric = "pr_auc_No")) %>%
  bind_rows(predicted_validation3990_data_KNN %>% 
              pr_auc(is_lipids, .pred_Yes) %>%
              mutate(.metric = "pr_auc_Yes")) %>% 
  bind_rows(predicted_validation3990_data_KNN %>% 
              roc_auc(is_lipids, .pred_No)) %>% 
  mutate(across(where(is.numeric), ~ round(., 3)))
a <- full_join(perf_validation_rf %>% rename(rf_estimate = .estimate),
               perf_validation_dt %>% rename(dt_estimate = .estimate), 
               by = c(".metric", ".estimator")) %>% 
  full_join(., perf_validation_n2dt %>% rename(n2dt_estimate = .estimate), 
            by = c(".metric", ".estimator")) %>% 
  full_join(., perf_validation_knn %>% rename(knn_estimate = .estimate), 
            by = c(".metric", ".estimator"))
a %>%
  mutate(.metric = case_when(
    .metric == "accuracy"           ~ "Accuracy",
    .metric == "roc_auc"            ~ "ROC AUC",
    .metric == "pr_auc_No"          ~ "PR AUC (No)",
    .metric == "pr_auc_Yes"         ~ "PR AUC (Yes)",
    TRUE                           ~ .metric
  )) %>%
  dplyr::select("Performance Metrics" = .metric,
                "Random Forest" = rf_estimate,
                "Decision Tree" = dt_estimate,
                "Trimmed Decision Tree" = n2dt_estimate,
                "KNN" = knn_estimate) %>% 
  kable()
# write_csv(a, "Performance table on validation set.csv")

# Figure 2 panel D
perf_validation_rf <- perf_validation_rf %>% 
  ggplot(aes(x=.metric, y=.estimate, fill=.metric))+
  geom_bar(stat = "identity")+
  labs(x= "", y= "Performance estimate", fill = NULL)+
  scale_fill_viridis_d(option = "F")+
  coord_flip()+
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        axis.ticks.y = element_blank())
perf_validation_rf
# ggsave("Figure 3 Panel AA Performance validation RF.pdf",
#        width = 2,
#        height = 3)
perf_validation_dt %>% 
  ggplot(aes(x=.metric, y=.estimate, fill=.metric))+
  geom_bar(stat = "identity")+
  labs(x= "", y= "Performance estimate", fill = NULL)+
  scale_fill_viridis_d(option = "F")+
  coord_flip()+
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        axis.ticks.y = element_blank())
perf_validation_dt
# ggsave("Figure 3 Panel AB Performance validation DT.pdf",
#        width = 2,
#        height = 3)
perf_validation_n2dt %>% 
  ggplot(aes(x=.metric, y=.estimate, fill=.metric))+
  geom_bar(stat = "identity")+
  labs(x= "", y= "Performance estimate", fill = NULL)+
  scale_fill_viridis_d(option = "F")+
  coord_flip()+
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        axis.ticks.y = element_blank())
perf_validation_n2dt
# ggsave("Figure 3 Panel AC Performance validation n2DT.pdf",
#        width = 2,
#        height = 3)
perf_validation_knn <-
  perf_validation_knn %>% 
  mutate(.metric = case_when(
    .metric == "accuracy"           ~ "Accuracy",
    .metric == "roc_auc"            ~ "ROC AUC",
    .metric == "pr_auc_No"          ~ "PR AUC (No)",
    .metric == "pr_auc_Yes"         ~ "PR AUC (Yes)",
    TRUE                            ~ .metric
  )) %>% 
  ggplot(aes(x=.metric, y=.estimate, fill=.metric))+
  geom_bar(stat = "identity")+
  labs(x= "", y= "Performance estimate", fill = NULL)+
  scale_fill_viridis_d(option = "F")+
  ylim(0,0.95)+
  coord_flip()+
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 2,byrow=TRUE))
perf_validation_knn
legend <- ggpubr::get_legend(perf_validation_knn)
# Convert to a ggplot and print
ggpubr::as_ggplot(legend)
# ggsave("Figure 3 - perf validation legend.pdf",
#        width = 2,
#        height = 0.5)

perf_validation_knn <- perf_validation_knn+
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank())
perf_validation_knn
# ggsave("Figure 3 Panel AD Performance validation KNN.pdf",
#        width = 2,
#        height = 3)

perf_validation_n2dt <- perf_validation_n2dt %>% 
  theme(axis.title = element_blank())


perf_validation_dt /
  perf_validation_n2dt /
  perf_validation_rf /
  perf_validation_knn +
  plot_layout(#axes = "collect", 
    axis_titles = "collect")
# ggsave("Figure 3 - perf validation.pdf",
#        width = 2.2,
#        height = 7)




# 2.Get predictions out on the testing dataset
# Figure 2 panel A Prediction - conf matrix - RF, DT, n2DT, KNN
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

conf_resuult_rf_testing <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted", 
       fill = "Overlap")+
  theme_classic()+
  theme(legend.position = "bottom")
conf_resuult_rf_testing
ggsave("Figure 2 Panel AA Prediction - conf matrix RF Dec2024.pdf",
       width = 7,
       height = 6)
final_rf %>%
  collect_predictions() %>%
  conf_mat(is_lipids, .pred_class)

df <- final_dt %>%
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

conf_resuult_dt_testing <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted", 
       fill = "Overlap")+
  theme_classic()+
  theme(legend.position = "bottom")
conf_resuult_dt_testing
ggsave("Figure 2 Panel AB Prediction - conf matrix DT Dec2024.pdf",
       width = 7,
       height = 6)

df <- final_trimmed_dt %>%
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

conf_resuult_n2dt_testing <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted", 
       fill = "Overlap")+
  theme_classic()+
  theme(legend.position = "bottom")
conf_resuult_n2dt_testing
ggsave("Figure 2 Panel AC Prediction - conf matrix n2DT Dec2024.pdf",
       width = 7,
       height = 6)

df <- final_knn %>%
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

conf_resuult_knn_testing <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted", 
       fill = "Overlap")+
  theme(legend.position = "bottom")
conf_resuult_knn_testing
ggsave("Figure 2 Panel AD Prediction - conf matrix KNN Dec 2024.pdf",
       width = 7,
       height = 6)

legend <- ggpubr::get_legend(conf_resuult_knn_testing)
# Convert to a ggplot and print
ggpubr::as_ggplot(legend)
# ggsave("Figure 3 - conf matrix legend.pdf",
#        width = 2.4,
#        height = 1)

conf_resuult_rf_testing <- conf_resuult_rf_testing+
  theme(legend.position = "none")
conf_resuult_dt_testing <- conf_resuult_dt_testing+
  theme(legend.position = "none")
conf_resuult_n2dt_testing <- conf_resuult_n2dt_testing+
  theme(legend.position = "none")
conf_resuult_knn_testing <- conf_resuult_knn_testing+
  theme(legend.position = "none")


conf_resuult_dt_testing /
  conf_resuult_n2dt_testing /
  conf_resuult_rf_testing /
  conf_resuult_knn_testing +
  plot_layout(axes = "collect")
ggsave("Figure 2A - conf matrix Dec2024.pdf",
       width = 2.2,
       height = 7)





# Get predictions validation 3990----
# Figure 2C
predicted_validation3990_data_RF <- predicted_validation3990_data_RF %>% 
  mutate_at(c("is_lipids"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-lipid"
  )) %>% 
  mutate_at(c(".pred_class"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-\nlipid"
  ))
conf_tab <- table(predicted_validation3990_data_RF$is_lipids, predicted_validation3990_data_RF$.pred_class)
conf_tab <- conf_tab / rowSums(conf_tab)
conf_tab <- as.data.frame(conf_tab, stringsAsFactors = TRUE)
conf_tab$Var2 <- factor(conf_tab$Var2, rev(levels(conf_tab$Var2)))

conf_resuult_rf_val <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted",
       fill = "Overlap")+
  theme(legend.position = "none")
conf_resuult_rf_val
# ggsave("Figure 2 Panel C Prediction validation - conf matrix RF Dec2024.pdf",
#        width = 7,
#        height = 6)
predicted_validation3990_data_RF %>%
  conf_mat(is_lipids, .pred_class)

predicted_validation3990_data_DT <- predicted_validation3990_data_DT %>% 
  mutate_at(c("is_lipids"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-lipid"
  )) %>% 
  mutate_at(c(".pred_class"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-\nlipid"
  ))
conf_tab <- table(predicted_validation3990_data_DT$is_lipids, predicted_validation3990_data_DT$.pred_class)
conf_tab <- conf_tab / rowSums(conf_tab)
conf_tab <- as.data.frame(conf_tab, stringsAsFactors = TRUE)
conf_tab$Var2 <- factor(conf_tab$Var2, rev(levels(conf_tab$Var2)))

conf_resuult_dt_val <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted",
       fill = "Overlap")+
  theme(legend.position = "none")
conf_resuult_dt_val
# ggsave("Figure 2 Panel C Prediction validation - conf matrix DT dec2024.pdf",
#        width = 7,
#        height = 6)

predicted_validation3990_data_n2DT <- predicted_validation3990_data_n2DT %>% 
  mutate_at(c("is_lipids"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-lipid"
  )) %>% 
  mutate_at(c(".pred_class"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-\nlipid"
  ))
conf_tab <- table(predicted_validation3990_data_n2DT$is_lipids, predicted_validation3990_data_n2DT$.pred_class)
conf_tab <- conf_tab / rowSums(conf_tab)
conf_tab <- as.data.frame(conf_tab, stringsAsFactors = TRUE)
conf_tab$Var2 <- factor(conf_tab$Var2, rev(levels(conf_tab$Var2)))

conf_resuult_n2dt_val <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted",
       fill = "Overlap")+
  theme(legend.position = "none")
conf_resuult_n2dt_val
# ggsave("Figure 2 Panel C Prediction validation - conf matrix n2DT Dec2024.pdf",
#        width = 7,
#        height = 6)

predicted_validation3990_data_KNN <- predicted_validation3990_data_KNN %>% 
  mutate_at(c("is_lipids"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-lipid"
  )) %>% 
  mutate_at(c(".pred_class"), ~ case_when(
    . == "Yes"       ~ "Lipid",
    . == "No"        ~ "Non-\nlipid"
  ))
conf_tab <- table(predicted_validation3990_data_KNN$is_lipids, predicted_validation3990_data_KNN$.pred_class)
conf_tab <- conf_tab / rowSums(conf_tab)
conf_tab <- as.data.frame(conf_tab, stringsAsFactors = TRUE)
conf_tab$Var2 <- factor(conf_tab$Var2, rev(levels(conf_tab$Var2)))

conf_resuult_knn_val <- ggplot(conf_tab, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = scales::percent(Freq)), size = 3) +
  scale_fill_gradient(low = "white", high = "#3575b5")+
  labs(x = "Expected", y = "Predicted",
       fill = "Overlap")+
  theme(legend.position = "none")
conf_resuult_knn_val
# ggsave("Figure 2 Panel C Prediction validation - conf matrix KNN Dec2024.pdf",
#        width = 7,
#        height = 6)

conf_resuult_dt_val /
  conf_resuult_n2dt_val /
  conf_resuult_rf_val /
  conf_resuult_knn_val +
  plot_layout(axes = "collect")
# ggsave("Figure 2C - conf matrix validation.pdf",
#        width = 2.2,
#        height = 7)


# Look at variable importance Figure 2B----
library(vip) 
# RF
all_workflows %>%
  extract_workflow("basic_Random_Forest") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Random_Forest") %>%
                      select_best(metric = "accuracy"))
# mtry = 1
# trees = 768
# min_n = 8
set.seed(123)
last_rf_mod <- 
  rand_forest(mtry = 1, min_n = 8, trees = 768) %>% 
  set_engine("ranger", importance = "impurity") %>% 
  set_mode("classification")

# the last workflow
last_rf_workflow <- 
  all_workflows %>%
  extract_workflow("basic_Random_Forest") %>%
  finalize_workflow(all_workflows %>%
                      extract_workflow_set_result("basic_Random_Forest") %>%
                      select_best(metric = "accuracy")) %>% 
  update_model(last_rf_mod)

# the last fit
set.seed(345)
last_rf_fit <- 
  last_rf_workflow %>% 
  last_fit(data_split)

last_rf_fit %>% 
  extract_fit_parsnip() %>% 
  vip(aesthetics = list(alpha = 0.5))+
  geom_text(aes(label = c("       Gradient \nComposition", 
                          "Mass /\ncharge")),
            nudge_y = c(-15, -10)
  )+
  scale_x_discrete(labels = c("buffer_percent" = "Gradient Composition", "row_m_z" = "Mass /charge"))+
  theme(axis.text.y = element_blank())
ggsave("Figure 2B VI RF Dec2024.pdf",
       width = 2,
       height = 2)

# DT
final_dt <- final_dt1

attributes(final_dt[[6]][[1]][["fit"]][["fit"]][["fit"]])[["ylevels"]][2] <- "L"
attributes(final_dt[[6]][[1]][["fit"]][["fit"]][["fit"]])[["ylevels"]][1] <- "NL"

dt_fit <-
  final_dt %>%
  extract_fit_parsnip()
dt_fit[["fit"]][["frame"]][["var"]]
dt_fit[["fit"]][["frame"]][["var"]][1] <- "GC"
dt_fit[["fit"]][["frame"]][["var"]][2] <- "GC"
dt_fit[["fit"]][["frame"]][["var"]][4] <- "M/Z"
#Generate Decision Tree Plot Using rpart.plot package
rpart.plot::rpart.plot(dt_fit$fit)
rpart.plot(dt_fit$fit, tweak = 1.6, type = 1)
# ggsave("Figure 2B DT plotw.pdf",
#        width = 2,
#        height = 2)

# DT
final_trimmed_dt <- final_trimmed_dt1
attributes(final_trimmed_dt[[6]][[1]][["fit"]][["fit"]][["fit"]])[["ylevels"]][2] <- "L"
attributes(final_trimmed_dt[[6]][[1]][["fit"]][["fit"]][["fit"]])[["ylevels"]][1] <- "NL"
dt_fit <-
  final_trimmed_dt %>%
  extract_fit_parsnip()
dt_fit[["fit"]][["frame"]][["var"]]
dt_fit[["fit"]][["frame"]][["var"]][1] <- "GC"
dt_fit[["fit"]][["frame"]][["var"]][3] <- "M/Z"
rpart.plot::rpart.plot(dt_fit$fit, tweak = 1.5, type = 1)
# ggsave("Figure 2B n2DT plot Dec2024.pdf",
#        width = 3,
#        height = 3)


# What is the m/z-rentention time pattern of predicted data over real data?---- 
# Figure S4
legend <- final_rf %>%
  collect_predictions() %>% dplyr::select(-is_lipids) %>% 
  bind_cols(testing(data_split)) %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point()+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  # scale_color_viridis_d(option = "A", 
  #                       begin = 0.2, end = 0.7,
  #                       name = NULL)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "bottom")
legend <- ggpubr::get_legend(legend)
# Convert to a ggplot and print
ggpubr::as_ggplot(legend)
ggsave("Figure S4 legend.pdf",
       width = 3,
       height = 3)

fig_S3_d3 <- final_rf %>%
  collect_predictions() %>% dplyr::select(-is_lipids) %>% 
  bind_cols(testing(data_split)) %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  mutate(is_lipids = ifelse(is_lipids == "Yes", "Lipid", "Non-lipid"), is_lipids = factor(is_lipids, levels = c("Non-lipid", "Lipid"))) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point(size = 0.5)+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  ylim(23, 70)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "none")
fig_S3_d3
ggsave("Figure S4 training prediction scatter plot RF.pdf",
       width = 3,
       height = 3)

# final_rf %>%
#   collect_predictions() %>% dplyr::select(-is_lipids) %>% 
#   bind_cols(testing(data_split)) %>% 
#   mutate(concordant = case_when(
#     is_lipids == .pred_class      ~ "Correctly classified",
#     is_lipids != .pred_class      ~ "Incorrectly classified",
#   )) %>% 
#   dplyr::select(hmdb, row_m_z, buffer_percent, is_lipids, .pred_class, concordant) %>% 
#   arrange(desc(concordant))

fig_S3_d1 <- final_dt %>%
  collect_predictions() %>% dplyr::select(-is_lipids) %>% 
  bind_cols(testing(data_split)) %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  mutate(is_lipids = ifelse(is_lipids == "Yes", "Lipid", "Non-lipid"), is_lipids = factor(is_lipids, levels = c("Non-lipid", "Lipid"))) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point(size = 0.5)+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  ylim(23, 70)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "none")
fig_S3_d1
ggsave("Figure S4 training prediction scatter plot DT.pdf",
       width = 3,
       height = 3)

fig_S3_d2 <- final_trimmed_dt %>%
  collect_predictions() %>% dplyr::select(-is_lipids) %>% 
  bind_cols(testing(data_split)) %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  mutate(is_lipids = ifelse(is_lipids == "Yes", "Lipid", "Non-lipid"), is_lipids = factor(is_lipids, levels = c("Non-lipid", "Lipid"))) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point(size = 0.5)+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  ylim(23, 70)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "none")
fig_S3_d2
ggsave("Figure S4 training prediction scatter plot n2DT.pdf",
       width = 3,
       height = 3)

fig_S3_d4 <- final_knn %>%
  collect_predictions() %>% dplyr::select(-is_lipids) %>% 
  bind_cols(testing(data_split)) %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  mutate(is_lipids = ifelse(is_lipids == "Yes", "Lipid", "Non-lipid"), is_lipids = factor(is_lipids, levels = c("Non-lipid", "Lipid"))) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point(size = 0.5)+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  ylim(23, 70)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "none")
fig_S3_d4
ggsave("Figure S4 training prediction scatter plot KNN.pdf",
       width = 3,
       height = 3)

# 3990 data----
set.seed(1234)
predicted_validation3990_data_RF <- augment(final_fitted_RF, validation3990_data)
set.seed(1234)
predicted_validation3990_data_DT <- augment(final_fitted_DT, validation3990_data)
set.seed(1234)
predicted_validation3990_data_n2DT <- augment(final_fitted_n2DT, validation3990_data)
set.seed(1234)
predicted_validation3990_data_KNN <- augment(final_fitted_KNN, validation3990_data)

fig_S3_r3 <- predicted_validation3990_data_RF %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  mutate(is_lipids = ifelse(is_lipids == "Yes", "Lipid", "Non-lipid"), is_lipids = factor(is_lipids, levels = c("Non-lipid", "Lipid"))) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point(size = 0.5)+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  ylim(23, 70)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "none")
fig_S3_r3
ggsave("Figure S4 validation prediction scatter plot RF.pdf",
       width = 3,
       height = 3)

fig_S3_r1 <- predicted_validation3990_data_DT %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  mutate(is_lipids = ifelse(is_lipids == "Yes", "Lipid", "Non-lipid"), is_lipids = factor(is_lipids, levels = c("Non-lipid", "Lipid"))) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point(size = 0.5)+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  ylim(23, 70)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "none")
fig_S3_r1
ggsave("Figure S4 validation prediction scatter plot DT.pdf",
       width = 3,
       height = 3)

fig_S3_r2 <- predicted_validation3990_data_n2DT %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  mutate(is_lipids = ifelse(is_lipids == "Yes", "Lipid", "Non-lipid"), is_lipids = factor(is_lipids, levels = c("Non-lipid", "Lipid"))) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point(size = 0.5)+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  ylim(23, 70)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "none")
fig_S3_r2
ggsave("Figure S4 validation prediction scatter plot n2DT.pdf",
       width = 3,
       height = 3)

fig_S3_r4 <- predicted_validation3990_data_KNN %>% 
  mutate(concordant = case_when(
    is_lipids == .pred_class      ~ "Correctly classified",
    is_lipids != .pred_class      ~ "Incorrectly classified",
  )) %>% 
  mutate(is_lipids = ifelse(is_lipids == "Yes", "Lipid", "Non-lipid"), is_lipids = factor(is_lipids, levels = c("Non-lipid", "Lipid"))) %>% 
  ggplot(aes(x=row_m_z, buffer_percent, color= concordant))+
  geom_point(size = 0.5)+
  scale_color_manual(name = NULL, values= c("#3399FF", "tomato"))+
  labs(x = "Mass/Charge", y = "Gradient Composition")+
  ylim(23, 70)+
  facet_wrap(. ~ is_lipids, ncol = 2)+
  theme(legend.position = "none")
fig_S3_r4
ggsave("Figure S4 validation prediction scatter plot KNN.pdf",
       width = 3,
       height = 3)


design <- "
  12
34
56
78
"

fig_S3_d1 + fig_S3_r1 + fig_S3_d2 + fig_S3_r2+
  fig_S3_d3 + fig_S3_r3 + fig_S3_d4 + fig_S3_r4 +
  plot_layout(axes = "collect_y",
              axis_titles = "collect")+ 
  plot_layout(design = design)

ggsave("Figure S4 All panels Dec2024.pdf",
       width = 7,
       height = 9)



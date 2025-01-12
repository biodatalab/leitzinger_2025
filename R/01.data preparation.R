# Import library
library(tidyverse)
library(gtsummary)
theme_set(theme_classic())
theme_gtsummary_compact()


################################################################ I. Load data----
# classification data
metabolites_classification <- 
  read_delim(paste0(here::here(), "/hmdb_keep_v4_python_blessed.txt.zip")) %>% 
  janitor::clean_names()

# metabolites expression data
omics1799_data <- 
  read_delim(
    paste0(
      here::here(), 
      "/flores_1799_tap73_metabolomics_reanalyze_2023-08-17/flores_1799_tap73_metabolomics_reanalysis_2023-08-17_iron_log2_merged.txt")) %>% 
  janitor::clean_names()
#
validation3990_data <- 
  read_delim(
    paste0(
      here::here(), 
      "/flores_3990_metabolomics_tissue_2023-08-17/flores_3990_metabolomics_tissue_iron_log2_merged.txt")) %>% 
  janitor::clean_names()


################################################################ II. Data exploration----
omics1799_data %>% 
  mutate(row_retention_time = round(row_retention_time, 1)) %>% 
  ggplot(aes(x= row_retention_time))+
  geom_bar(stat = "count")+
  geom_vline(aes(xintercept= 1), color= "blue")+
  geom_vline(aes(xintercept= 14), color= "blue")+
  xlim(0, 15)+
  theme_classic()

omics1799_data %>% 
  ggplot(aes(x= row_retention_time, y=row_m_z
  ))+
  geom_point(size= 1)+
  geom_rect(aes(xmin = 0, xmax = 1+0.3, ymin = min(row_m_z)-10, ymax = max(row_m_z) +2), 
            color= "red", fill="transparent") +
  geom_rect(aes(xmin = 14-0.2, xmax = max(row_retention_time)+0.3, ymin =  min(row_m_z)-10, ymax = max(row_m_z) +2),
            color= "red", fill="transparent")

omics1799_data %>%
  separate(row_id, into = c("charge", "id")) %>%
  ggplot(aes(x=charge, fill= charge))+
  geom_bar()+
  ggtitle("In overall metabolites")

omics1799_data %>%
  filter(non_heavy_identified_flag == 1) %>% 
  separate(row_id, into = c("charge", "id")) %>%
  ggplot(aes(x=charge, fill= charge))+
  geom_bar()+
  ggtitle("In non_heavy_identified_flag == 1")


################################################################ III. Clean data----
clean_metabolites <- omics1799_data %>% 
  # cleaning retention time to eliminate false detection during set up, washing and equilibrating the HPLC
  filter(row_retention_time > 1 & row_retention_time < 14)

validation3990_data <- validation3990_data %>% 
  filter(row_retention_time > 1 & row_retention_time < 14)

# Plot resulting data
clean_metabolites %>% 
  mutate(row_retention_time = round(row_retention_time, 1)) %>% 
  ggplot(aes(x= row_retention_time))+
  geom_bar(stat = "count")+
  geom_vline(aes(xintercept= 1), color= "blue")+
  geom_vline(aes(xintercept= 14), color= "blue")+
  xlim(0, 15)+
  theme_classic()

clean_metabolites %>% 
  ggplot(aes(x= row_retention_time, y=row_m_z
  ))+
  geom_point(size= 1)+
  xlim(0, 15)

# Minimize overlap in identified metabolites between both polarities
# Priority is given to metabolite found in the negative ion mode due to the larger in dataset
clean_metabolites1 <- clean_metabolites %>%
  filter(non_heavy_identified_flag == 1) %>% 
  # remove same metabolite if present as neg and pos
  separate(row_id, into = c("charge", "id"), remove = FALSE) %>%
  # give priority to negative metabolite
  arrange(id, charge) %>%
  distinct(id, .keep_all = TRUE)

# select 1 in non identified metabolites
clean_metabolites2 <- clean_metabolites %>%
  filter(non_heavy_identified_flag == 0) %>% 
  # remove same metabolite if present as neg and pos
  separate(row_id, into = c("charge", "id"), remove = FALSE) %>%
  arrange(id, charge) %>%
  distinct(id, .keep_all = TRUE)

# Bind
clean_metabolites <- bind_rows(clean_metabolites1, clean_metabolites2) %>% 
  distinct(id, .keep_all = TRUE)

# Plot resulting data
clean_metabolites %>%
  # filter metabolites in positive or negative detection
  separate(row_id, into = c("charge", "id")) %>%
  filter(duplicated(id)) %>%
  full_join(., omics1799_data %>%
              separate(row_id, into = c("charge", "id")) %>%
              filter(!duplicated(id)),
            by="id") %>%
  dplyr::select(charge.x, id, charge.y) %>%
  filter(is.na(charge.x) | is.na(charge.y)) %>%
  pivot_longer(cols = c(charge.x, charge.y)) %>%
  filter(!is.na(value)) %>%
  ggplot(aes(x=value, fill= value))+
  geom_bar()+
  ggtitle("In non_heavy_identified_flag == 1")

# Do the same for the validation data
clean_metabolites1 <- validation3990_data %>%
  filter(non_heavy_identified_flag == 1) %>% 
  # remove same metabolite if present as neg and pos
  separate(row_id, into = c("charge", "id"), remove = FALSE) %>%
  arrange(id, charge) %>%
  distinct(id, .keep_all = TRUE)

clean_metabolites2 <- validation3990_data %>%
  filter(non_heavy_identified_flag == 0) %>% 
  # remove same metabolite if present as neg and pos
  separate(row_id, into = c("charge", "id"), remove = FALSE) %>%
  arrange(id, charge) %>%
  distinct(id, .keep_all = TRUE)

# Bind
validation3990_data <- bind_rows(clean_metabolites1, clean_metabolites2) %>% 
  distinct(id, .keep_all = TRUE)

# Clean
rm(clean_metabolites1, clean_metabolites2)

# Add taxonomy to our data and create the dichotomous lipids variable Yes/No aka "Lipid/Not Lipid" in paper
full_data <- clean_metabolites %>%
  # Join with known metabolites
  mutate(hmdb = str_remove(hmdb, "^\\|")) %>% 
  separate_wider_delim(cols = hmdb, delim = "|",
                       names = c("hmdb"), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  left_join(., metabolites_classification,
            by = c("hmdb" = "accession"))

# Do the same for validation dataset
validation3990_data <- validation3990_data %>%
  mutate(hmdb = str_remove(hmdb, "^\\|")) %>% 
  separate_wider_delim(cols = hmdb, delim = "|",
                       names = c("hmdb"), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  left_join(., metabolites_classification,
            by = c("hmdb" = "accession"))

full_data <- full_data %>%
  mutate(is_lipids = case_when(
    row_id == "neg_00049"                                   ~ "No",
    str_detect(taxonomy_super_class, "Lipids")              ~ "Yes",
    TRUE                                                    ~ "No"
  )) %>%
  mutate_if(is.character, factor)

validation3990_data <- validation3990_data %>%
  filter(non_heavy_identified_flag == 1) %>% # Can already do this step here
  mutate(is_lipids = case_when(
    str_detect(taxonomy_super_class, "Lipids")              ~ "Yes",
    TRUE                                                    ~ "No"
  )) %>%
  mutate_if(is.character, factor)

# Calculate buffer concentration (aka gradient concentration) based on retention time
full_data <- full_data %>%
  mutate(buffer_percent = case_when(
    row_retention_time >= 0 &
      row_retention_time <= 13            ~  4.615 * row_retention_time + 20
  ))

validation3990_data <- validation3990_data %>%
  mutate(buffer_percent = case_when(
    row_retention_time >= 0 &
      row_retention_time <= 13            ~  4.615 * row_retention_time + 20
  ))
# write_rds(validation3990_data, "validation3990_data_with_nonidentified_metabolites.rds")
# write_rds(validation3990_data, "validation3990_data.rds")

################################################################ IV. Explore clean data----
full_data %>% 
  filter(non_heavy_identified_flag == 0) %>% 
  ggplot(aes(x= row_retention_time, y=row_m_z))+
  geom_point(color= "yellow")+
  geom_point(data= omics1799_data %>% filter(non_heavy_identified_flag == 1),
             aes(x= row_retention_time, y=row_m_z), color= "red")+
  ggtitle("non_heavy_identified_flag == 1 are in red")

tbl <- full_data %>% 
  dplyr::select(non_heavy_identified_flag) %>% 
  mutate(non_heavy_identified_flag = case_when(
    non_heavy_identified_flag == 1            ~ "identified",
    TRUE                                      ~ "non identified"
  )) %>% 
  tbl_summary()
tbl

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
  ungroup()

a %>% 
  ggplot(aes(x= "", fill= non_heavy_identified_flag))+
  geom_bar(width=1, color= "white") +
  coord_polar("y", start=0)+
  theme_void()+
  
  geom_text(data = a %>%
              distinct(non_heavy_identified_flag, .keep_all = TRUE),
            mapping = aes(x=1.6, y= ypos, label = non_heavy_identified_flag
            ))

# Create Final minimal data including 2 predictors (m/z and gradient concentration) and the outcome (is_lipids)
two_predictor_data <- full_data %>%
  # Filter identified metabolites
  filter(non_heavy_identified_flag == 1) %>%
  dplyr::select(hmdb,
                row_m_z, buffer_percent,
                is_lipids)

# Save minimal data
write_rds(two_predictor_data, "two_predictor_data.rds")


# Last check

two_predictor_data %>%
  ggplot(aes(x=is_lipids, fill= is_lipids))+
  geom_bar()

two_predictor_data %>% 
  dplyr::select(is_lipids) %>% 
  tbl_summary(type = list(is_lipids ~ "categorical"))

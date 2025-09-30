# Import library
install.packages("tidyverse")
install.packages("gtsummary")
install.packages("janitor")

library(tidyverse)
library(gtsummary)
theme_set(theme_classic())
theme_gtsummary_compact()


################################################################ I. Load data----
# classification data
url1 <- "https://zenodo.org/records/17236618/files/hmdb_keep_v4_python_blessed.txt?download=1"
metabolites_classification <- 
  read.delim(url1, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
  janitor::clean_names()

# Discovery dataset
url2 <- "https://zenodo.org/records/17236618/files/omics1799_data.txt?download=1"
omics1799_data <- 
  read.delim(url2, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
  janitor::clean_names()

# Validation dataset
url3 <- "https://zenodo.org/records/17236618/files/validation3990_data.txt?download=1"
validation3990_data <- 
  read.delim(url3, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
  janitor::clean_names()

# Independent dataset
url4 <- "https://zenodo.org/records/17236618/files/ovca_metabolomics.tsv?download=1"
validation_ovca <- 
  read.delim(url4, header = TRUE, sep = "\t", stringsAsFactors = FALSE) %>% 
  janitor::clean_names() %>% 
  rename(row_retention_time = rt, row_m_z = m_z, 
         hmdb = hmdb_id) %>% 
  mutate(hmdb1 = str_extract(hmdb, "(\\d+)")) %>% 
  mutate(hmdb = paste0("HMDB00", hmdb1))

rm(url1, url2, url3, url4)


################################################################ II. Data exploration----
# Quick check of the first few rows
head(metabolites_classification)
head(omics1799_data)
head(validation3990_data)
head(validation_ovca)

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
# Exclude RT <1/>14 min
# cleaning retention time to eliminate false detection during set up, washing and equilibrating the HPLC
clean_metabolites <- omics1799_data %>% 
  filter(row_retention_time > 1 & row_retention_time < 14)

validation3990_data <- validation3990_data %>% 
  filter(row_retention_time > 1 & row_retention_time < 14)

validation_ovca <- validation_ovca %>% 
  filter(row_retention_time > 1 & row_retention_time < 14)

# Quick look at the resulting data
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
# Priority is given to metabolite found in the negative ion mode due to the larger dataset
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

# Do the same for the validation data - NOT RUN - all metabolites are detected in the positive polarity
# clean_metabolites1 <- validation_ovca %>% 
#   # filter(non_heavy_identified_flag == 1) %>% # Not needed for this type of data
#   # remove same metabolite if present as neg and pos
#   separate(metabolite_id, into = c("charge", "id"), remove = FALSE) %>%
#   arrange(id, charge) %>%
#   distinct(id, .keep_all = TRUE)

# clean_metabolites2 <- validation3990_data %>%
#   filter(non_heavy_identified_flag == 0) %>% 
#   # remove same metabolite if present as neg and pos
#   separate(row_id, into = c("charge", "id"), remove = FALSE) %>%
#   arrange(id, charge) %>%
#   distinct(id, .keep_all = TRUE)

# Bind
# validation_ovca <- bind_rows(clean_metabolites1, clean_metabolites2) %>% 
#   distinct(id, .keep_all = TRUE)

# Clean
rm(clean_metabolites1, clean_metabolites2)

# Add annotation (taxonomy) to our data for creating the dichotomous lipids variable Yes/No aka "Lipid/Not Lipid" in paper
full_data <- clean_metabolites %>%
  mutate(hmdb = str_remove(hmdb, "^\\|")) %>% 
  # data "hdmi" ids is a string with multiple hdmi ids separated with a "|" character
  # need to separate them
  separate_wider_delim(cols = hmdb, delim = "|",
                       names = c("hmdb"), 
                       too_few = "align_start", too_many = "drop", 
                       cols_remove = TRUE) %>% 
  # Annotate with known metabolites HMDB/KEGG/PubChem
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

# Do the same for validation dataset
validation_ovca <- validation_ovca %>%
  # data "hdmi" ids are unique so no need to separate
  # mutate(hmdb = str_remove(hmdb, "^\\|")) %>% 
  # separate_wider_delim(cols = hmdb, delim = "|",
  #                      names = c("hmdb"), 
  #                      too_few = "align_start", too_many = "drop", 
  #                      cols_remove = TRUE) %>% 
  # Annotate with known metabolites HMDB/KEGG/PubChem
  inner_join(., metabolites_classification,
            by = c("hmdb" = "accession"))

# Create the dichotomous lipids variable Yes/No aka "Lipid/Not Lipid" in paper
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

validation_ovca <- validation_ovca %>%
  filter(!is.na(taxonomy_super_class)) %>% # Can already do this step here
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
write_rds(validation3990_data, "validation3990_data_with_nonidentified_metabolites.rds")

validation_ovca <- validation_ovca %>%
  mutate(buffer_percent = case_when(
    row_retention_time >= 0 &
      row_retention_time <= 13            ~  4.615 * row_retention_time + 20
  ))
write_rds(validation_ovca, "clean_validation_ovca.rds")

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

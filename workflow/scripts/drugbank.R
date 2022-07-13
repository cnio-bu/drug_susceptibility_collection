log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

library("tidyverse")

### SNAKEMAKE I/O ###
drugbank <- snakemake@input[[1]]
out <- snakemake@output[[1]]
blacklist <- snakemake@params[["blacklist"]]

### CODE ###
# Read DrugBank object and blacklist
load(drugbank)
blacklist <- read.csv(blacklist, header = TRUE, sep = "\t")

# DrugBank names
dnames <- drugBank$general_information %>% 
  rename(preferred.drug.bank = name, 
         DrugBank.ID = primary_key) %>%
         select(preferred.drug.bank, DrugBank.ID) %>%
         unique

# Synonyms
synonyms <- drugBank$synonyms %>% 
  separate_rows(language, sep = "/") %>% 
  filter(language == "english") %>% 
  rename(drug = synonym, DrugBank.ID = "drugbank-id") %>% 
  select(drug, DrugBank.ID) %>% 
  unique

# Brands
brands <- drugBank$international_brands %>% 
  rename(drug = brand, DrugBank.ID = "drugbank-id") %>% 
  select(-company) %>% 
  unique

# Merge dataframes
drugbank <- dnames %>%
  merge(synonyms, all = TRUE) %>% 
  merge(brands, all = TRUE)

# Collapsed names
drugbank <- drugbank %>% 
  group_by(DrugBank.ID) %>%
  mutate(drug = toupper(drug),
         preferred.drug.bank = toupper(unique(na.omit(preferred.drug.bank))),
         collapsed.name = str_replace_all(drug, "[^[:alnum:]]", ""),
         collapsed.name = 
           case_when(collapsed.name == drug | 
                       collapsed.name == preferred.drug.bank ~ as.character(NA),
                     TRUE ~ collapsed.name)) %>% 
  unique

# Remove incorrect synonyms
drugbank <- anti_join(drugbank, blacklist)

# Save output
write.table(drugbank, file = out, sep = "\t", quote = FALSE, row.names = FALSE)
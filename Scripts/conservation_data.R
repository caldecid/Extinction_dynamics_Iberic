# IUCN data ---------------------------------------------------------------

######libraries
library(tidyverse)
library(readxl)



####Iberic Peninsula data
excel_sheets("Data/Raw/Lista_Peninsula_Andalucia_Oriental_Final.xls")

peninsula <- read_excel("Data/Raw/Lista_Peninsula_Andalucia_Oriental_Final.xls",
                        sheet = "Lista_Peninsula_Final.txt" )


##ordering the conservation status
peninsula$Status <- factor(peninsula$Status, labels = c("DD", "NE",
                                                        "LC", "NT",
                                                        "VU", "EN",
                                                        "CR", "RE",
                                                        "EW", "EX"),
                           ordered = TRUE)

peninsula <- peninsula %>% mutate(threat_category = case_when(
  Status %in% c("DD", "NE",
                "LC", "NT") ~ "not_threatened",
  Status %in% c("VU", "EN", "CR", "EX",
                "RE", "EW") ~ "threatened"))


########Genus arrangements############
peninsula <- peninsula %>%
  mutate(species = Species) %>%  # Keep the original column
  separate(col = Species, into = c("Genus", "Epithet"), sep = " ",
           extra = "merge")


###saving
write_csv(peninsula, file = "Data/Processed/peninsula_iucn.csv")


#######Andalucia data

andalucia <- read_excel("Data/Raw/Lista_Peninsula_Andalucia_Oriental_Final.xls",
                        sheet = "FLORANDOR")

andalucia$Status <- factor(andalucia$Status, labels = c("DD",
                                                        "LC", "NT",
                                                        "VU", "EN",
                                                        "CR", "EX"),
                           ordered = TRUE)

andalucia <- andalucia %>% mutate(threat_category = case_when(
  Status %in% c("DD", 
                "LC", "NT") ~ "not_threatened",
  Status %in% c("VU", "EN", "CR", "EX") ~ "threatened"))

##Genus arrangements
andalucia <- andalucia  %>%
  mutate(species = Species) %>%  # Keep the original column
  separate(col = Species, into = c("Genus", "Epithet"), sep = " ",
           extra = "merge")

##saving
write_csv(andalucia, file = "Data/Processed/andalucia_iucn.csv")

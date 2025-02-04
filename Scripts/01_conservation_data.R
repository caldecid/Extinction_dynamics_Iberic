# IUCN data

## Libraries
library(tidyverse)
library(readxl)

## Iberic Peninsula data
path <- "Data/Raw/Lista_Peninsula_Andalucia_Oriental_Final.xls"
excel_sheets(path)
peninsula <- read_excel(path,
                        sheet = "Lista_Peninsula_Final.txt" )


### Ordering the conservation status
not_threatened <- c("DD", "NE", "LC", "NT")
threatened <- c("VU", "EN", "CR", "EX", "RE", "EW")


peninsula <- peninsula %>%
  mutate(
    Status = factor(
      Status,
      labels = c(not_threatened, threatened),
      ordered = TRUE
    ),
    threat_category = case_when(
      Status %in% not_threatened ~ "not_threatened",
      Status %in% threatened ~ "threatened"
    )
  )


### Genus arrangements
peninsula <- peninsula %>%
  mutate(species = Species) %>%  # Keep the original column
  separate(col = Species, into = c("Genus", "Epithet"), sep = " ",
           extra = "merge")


### Save
write_csv(peninsula, file = "Data/Processed/peninsula_iucn.csv")


## Andalucia data

andalucia <- read_excel(path,
                        sheet = "FLORANDOR")


andalucia <- andalucia %>%
  mutate(
    Status = factor(
      Status,
      labels = c(not_threatened[-2], threatened[-(5:6)]),
      ordered = TRUE
    ),
    threat_category = case_when(
      Status %in% not_threatened ~ "not_threatened",
      Status %in% threatened ~ "threatened"
    )
  )

### Genus arrangements
andalucia <- andalucia  %>%
  mutate(species = Species) %>%  # Keep the original column
  separate(
    col = Species,
    into = c("Genus", "Epithet"),
    sep = " ",
    extra = "merge"
  )

##saving
write_csv(andalucia, file = "Data/Processed/andalucia_iucn.csv")

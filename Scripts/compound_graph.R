library(stringi)
library(tidyverse)
library(tidyr)
library(xlsx)
library(chemodiv)
library(ggplot2)
library(ggplus)
library(ggforce)
library(httr)
library(jsonlite)
library(webchem)
# devtools::install_github("ricardo-bion/ggradar", 
#                          dependencies = TRUE)
library("ggradar")
library("fmsb")
library(scales)
source("Scripts/NP_classifier_Seb.R")

#READING EVERYTHING IN
source("Scripts/Reading_in_data.R")

MS_corrected <- read.xlsx("Output/Key_compounds.xlsx", sheetIndex = 1) %>%
  rename(Compound.Name = query) %>%
  mutate(Compound.Name = tolower(Compound.Name))
EAD_corrected <- read.xlsx("Output/EAD_Key_compounds.xlsx",sheetName = "corrected") %>%
  mutate(Compound.Name = tolower(Compound.Name))

All_corrected_old <- MS_corrected %>%
  full_join(EAD_corrected) %>% 
  mutate(Compound.Name = ifelse(is.na(corrected_name),
                                Compound.Name,
                                corrected_name)) %>%
  distinct(Compound.Name, .keep_all = T)


CIDs <- get_cid(All_corrected_old$Compound.Name, match = "first")



CIDs <- CIDs %>%
  filter(!is.na(cid))

Com_props <- pc_prop(CIDs$cid)


Grouped_compounds <- Npclassifier_API(Com_props) %>%
  rename(cid = CID) %>%
  mutate(cid = as.character(cid))

Grouped_compounds <- Grouped_compounds %>%
  left_join(CIDs) %>%
  rename(Compound.Name = query) %>%
  mutate(Compound.Name = tolower(Compound.Name)) %>%
  unique()

All_MS <- All_MS %>%
  mutate(Compound.Name = tolower(Compound.Name)) %>%
  left_join(MS_corrected) %>%
  mutate(Compound.Name = ifelse(is.na(corrected_name),
                                Compound.Name,
                                corrected_name)) %>%
  unique() %>%
  mutate(Compound.Name = tolower(Compound.Name)) %>% 
  left_join(Grouped_compounds[7,42:48])
# 
# MS_compounds <- Grouped_compounds %>% 
#   full_join(All_MS, by = "Compound.Name")

MS_compounds <- read.xlsx("Data/MS_compounds_properties.xlsx", sheetIndex = 1) %>% 
  mutate(Superclass_results2 = as.character(Superclass_results2))
  
 
All_metcalfa_data_norm <- All_metcalfa_data_norm %>%
    # mutate(Compound.Name = tolower(Compound.Name)) %>%
    left_join(EAD_corrected, by = "Compound.Name") %>%
    mutate(Compound.Name = ifelse(is.na(corrected_name),
                                  Compound.Name,
                                  corrected_name)) %>% 
  left_join(Grouped_compounds)
# EAD_compounds <- All_metcalfa_data_norm %>%
#   left_join(Grouped_compounds, by = "Compound.Name")

EAD_compounds <- read.xlsx("Data/EAD_compounds_properties.xlsx", sheetIndex = 1) %>% 
  mutate(Superclass_results2 = as.character(Superclass_results2),
         Comments = as.character(Comments))

All_metcalfa_data_norm <- All_metcalfa_data_norm %>% 
  left_join(EAD_compounds)

max_EAD_responses <- max(All_metcalfa_data_norm$Normalized_response)
max_MS_area <- max(All_MS$Area.from.MS)

filtered_EAD <- All_metcalfa_data_norm %>%
  filter(Sex != "nimph") %>%
  group_by(Compound.Name, Sample_name) %>%
  
  summarise(Value = mean(Normalized_response))
max_EAD_bla <- max(filtered_EAD$Value)

compounds_order <- Grouped_compounds %>% 
  group_by(Compound.Name, Superclass_results1) %>%
  summarise() %>% 
  arrange(desc(Superclass_results1)) 
  # filter(is.na(Superclass_results1))
         # !Compound.Name %in% c("benzaldehyde","benzyl alcohol", "o-methylanisole","p-benzoquinone", "2,6-di-tert-butyl-","p-xylene"))
compounds_order_vector <- compounds_order$Compound.Name

EAD_radar<- All_metcalfa_data_norm %>%
  # mutate(Compound.Name = ifelse(is.na(corrected_name), 
  #                               Compound.Name,
  #                               corrected_name)) %>% 
  # mutate(Compound.Name = paste(Superclass_results1, Compound.Name, sep = "_")) %>%
  group_by(Compound.Name, Sample_name) %>% 
  # filter(Sample_name == "TAGETES") %>% 
  summarise(Value = mean(Normalized_response)) %>% 
  mutate(r_radarplot = Value/max_EAD_bla) %>% 
  mutate(group = "EAD") %>% 
  dplyr::select(- Value) %>% 
  rename(Value = r_radarplot) %>% 
  mutate(Compound.Name = tolower(Compound.Name)) 

filtered_MS <- All_MS  
  # filter(`Sample descriptor` == "TAGETES")
max_MS <- max(filtered_MS$Area.from.MS)

MS_radar <- All_MS %>% 
  # filter(`Sample descriptor` == "TAGETES") %>%
  mutate(r_radarplot = (Area.from.MS * 1)/max(Area.from.MS)) %>% 
  mutate(Value = r_radarplot) %>% 
  dplyr::select(Compound.Name, Value, 'Sample descriptor') %>% 
  mutate(group = "MS")  %>% 
  rename(Sample_name = 'Sample descriptor')

FID_ng <- readxl::read_excel("Data/FID_NG.xls", sheet = "conc", col_types = c("text", "skip", "text", "skip", "skip", "skip", "skip", "skip", "numeric", "skip", "numeric", "skip", "numeric"), n_max = 30)  
FID_ng2 <- FID_ng %>% 
  gather(-c(Compound.Name,`functional groups`), key = Sample_name, value = ng)

Max_NG <- max(FID_ng2$ng)
FID_radar <- FID_ng2 %>% 
  mutate(Value = (ng * 1)/max(ng),
         group = "FID") %>% 
  dplyr::select(Compound.Name, Value, Sample_name, group)
  
  

Compounds_EAD_ <- EAD_radar %>% 
  group_by(Compound.Name) %>%
  summarise() %>% 
  left_join(Grouped_compounds) %>% 
  arrange(desc(Superclass_results1)) 

Compounds_EAD_ <- read.xlsx("Output/Compounds_for_radar_order.xlsx", sheetIndex = 1)

Compounds_EAD_ <- Compounds_EAD_ %>% 
  arrange(desc(Superclass_results2))
 
Compounds_EAD_vector <- Compounds_EAD_$Compound.Name
Superclass_for_coloring <- Compounds_EAD_$Superclass_results2
factor(Superclass_for_coloring)

Radar_comps_all <- FID_radar %>% 
  bind_rows(EAD_radar) %>% 
  mutate(Compound.Name = tolower(Compound.Name)) %>% 
  filter(group == "EAD") %>% 
  group_by(Compound.Name, group) %>% 
  summarise(Value = mean(Value)) %>% 
  ungroup() %>% 
  mutate(Sample_name = "All") %>% 
  spread(Compound.Name, Value, fill = 0) %>% 
  dplyr::select(Sample_name, group, Compounds_EAD_vector)
  
  

Radar_data <- FID_radar %>% 
  bind_rows(EAD_radar) %>% 
  mutate(Compound.Name = tolower(Compound.Name)) %>%
  # filter(group == "EAD") %>% 
  # group_by(Compound.Name) %>% summarise
  # group_by(group, Sample_name) %>%
  # summarise()
  # mutate(across(where(is.numeric), scales::rescale(.,to = c(0,1)))) %>%
  # ungroup() %>%
  # arrange(match(Compound.Name, compounds_order_vector)) %>% 
  # pivot_wider(names_from = Compound.Name, 
  #             values_from = Value, values_fill = 0) %>% 
  # group_by(group, Sample_name) %>%
  # mutate(Compound.Name = factor(Compound.Name, levels = compounds_order_vector)) %>% 
  # ungroup() %>% 
  filter(!is.na(Compound.Name)) %>% 
  spread(Compound.Name, Value, fill = 0) %>% 
  dplyr::select(Sample_name, group, Compounds_EAD_vector)
  # select_if(~ !any(is.na(.))) 
  

RADAR_EAD_TAG <- Radar_data[5,]
RADAR_MS_TAG <- Radar_data[6,]
RADAR_EAD_AIL <- Radar_data[1,]
RADAR_MS_AIL <- Radar_data[2,]
RADAR_EAD_ARIST <- Radar_data[3,]
RADAR_MS_ARIST <- Radar_data[4,]

RADAR_EAD_ALL <- Radar_data[c(5,1,3),] %>% 
  gather(-c(Sample_name,group), key = "Compound", value = "mV") %>% 
  group_by(Compound) %>% 
  summarise(meanMV = mean(mV)) %>% 
  mutate(Sample_name = "All") %>% 
  spread(key = Compound, value = meanMV) %>% 
  dplyr::select(Sample_name, Compounds_EAD_vector)


circle_coords <- function(r, n_axis = ncol(RADAR_EAD_TAG) - 2){
  fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
  x <- r*cos(fi)
  y <- r*sin(fi)
  
  tibble(x, y, r)
}

central_distance <- 0.1

step_1 <- map_df(seq(0, 1, 0.25) + central_distance, circle_coords) %>%
  ggplot(aes(x, y)) +
  geom_polygon(data = circle_coords(1 + central_distance), 
               alpha = 1, fill = 'transparent') +
  geom_path(aes(group = r), lty = 2, alpha = 0.5) +
  theme_void()

 axis_coords <- function(n_axis){
  fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2
  x1 <- central_distance*cos(fi)
  y1 <- central_distance*sin(fi)
  x2 <- (1 + central_distance)*cos(fi)
  y2 <- (1 + central_distance)*sin(fi)
  
  tibble(x = c(x1, x2), y = c(y1, y2), id = rep(1:n_axis, 2))
}

step_2 <- step_1 + geom_line(data = axis_coords(ncol(RADAR_EAD_TAG) - 2), 
                             aes(x, y, group = id), alpha = 0.3) +
                        geom_line(data = axis_coords(ncol(RADAR_MS_TAG) - 2), 
                                  aes(x, y, group = id), alpha = 0.3)

text_data <- RADAR_EAD_TAG %>%
  dplyr::select(-group, -Sample_name) %>%
  map_df(~ min(.) + (max(.) - min(.)) * seq(0, 1, 0.25)) %>%
  mutate(r = seq(0, 1, 0.25)) %>%
  pivot_longer(-r, names_to = "parameter", values_to = "value")

text_coords <- function(r, n_axis = ncol(RADAR_EAD_TAG) - 2){
  fi <- seq(0, (1 - 1/n_axis)*2*pi, (1/n_axis)*2*pi) + pi/2 + 0.01*2*pi/r
  x <- r*cos(fi)
  y <- r*sin(fi)
  
  tibble(x, y, r = r - central_distance)
}

labels_data <- map_df(seq(0, 1, 0.25) + central_distance, text_coords) %>%
  bind_cols(text_data %>% dplyr::select(-r))
labels_axis_EAD <- map_df(seq(0, 1, 0.25) + central_distance, text_coords) %>%
  bind_cols(text_data %>% dplyr::select(-r)) %>% filter(parameter == "alpha-copaene")
labels_axis_MS <- map_df(seq(0, 1, 0.25) + central_distance, text_coords) %>%
  bind_cols(text_data %>% dplyr::select(-r)) %>% filter(parameter == "1-hexanol")

step_3 <- step_2 + 
  geom_text(data = labels_axis_EAD, aes(x, 
                                        y, 
                                        label = c(round(x = c(0, (0.25*max_EAD_bla), (.5*max_EAD_bla), (.75*max_EAD_bla), max_EAD_bla), digits = 2)),
                                                  # round(x = c(0, (0.25*bla), (.05*bla), (.75*bla), bla), digits = 2)),
                                                  # round(x = c(0, (0.25*bla), (.05*bla), (.75*bla), bla), digits = 2)),
                                        
                                        ),colour = "#b10026") +
  geom_text(data = labels_axis_MS, aes( x, 
                                        y, 
                                        label = c(round(x = c(0, (0.25*Max_NG), (.5*Max_NG), (.75*Max_NG),Max_NG), digits = 0)),
                                         
                                        ),colour = "#fd8d3c") +
  annotate("text", x=.5, y=-1.5, label= "FID (ng/ul)", colour = "#fd8d3c") +
  annotate("text", x=-.5, y=-1.5, label= "EAD (mV)", colour = "#b10026") +
  geom_text(data = text_coords(1 + central_distance + 0.2), 
            aes(x, y, color = Superclass_for_coloring), 
            label = labels_data$parameter[1:(ncol(RADAR_EAD_TAG)-2)]) +
  scale_colour_manual(values=c("#005a32","#41ab5d","#3690c0","#fb8072","#034e7b", "#67001f", "#000000", "#8c96c6"))

rescaled_coords <- function(r, n_axis){
  fi <- seq(0, 2*pi, (1/n_axis)*2*pi) + pi/2
  tibble(r, fi) %>% mutate(x = r*cos(fi), y = r*sin(fi)) %>% dplyr::select(-fi)
}

rescaled_data_EAD_TAG <- RADAR_EAD_TAG %>% 
  # mutate(across(-group, rescale)) %>%
  mutate(copy = pull(., 3)) %>% 
  pivot_longer(-c(group, Sample_name), names_to = "parameter", values_to = "value") %>%
  # group_by(group) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(RADAR_EAD_TAG) - 2)) %>%
  unnest(cols = c(coords))

rescaled_data_MS_TAG<- RADAR_MS_TAG %>% 
  # mutate(across(-group, rescale)) %>%
  mutate(copy = pull(., 3)) %>% 
  pivot_longer(-c(group, Sample_name), names_to = "parameter", values_to = "value") %>%
  group_by(group) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(RADAR_MS_TAG) - 2)) %>%
  unnest(cols = c(coords))

rescaled_data_EAD_AIL<- RADAR_EAD_AIL %>% 
  # mutate(across(-group, rescale)) %>%
  mutate(copy = pull(., 3)) %>% 
  pivot_longer(-c(group, Sample_name), names_to = "parameter", values_to = "value") %>%
  group_by(group) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(RADAR_EAD_TAG) - 2)) %>%
  unnest(cols = c(coords))

rescaled_data_MS_AIL<- RADAR_MS_AIL %>% 
  # mutate(across(-group, rescale)) %>%
  mutate(copy = pull(., 3)) %>% 
  pivot_longer(-c(group, Sample_name), names_to = "parameter", values_to = "value") %>%
  group_by(group) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(RADAR_MS_TAG) - 2)) %>%
  unnest(cols = c(coords))

rescaled_data_EAD_ARIST<- RADAR_EAD_ARIST %>% 
  # mutate(across(-group, rescale)) %>%
  mutate(copy = pull(., 3)) %>% 
  pivot_longer(-c(group, Sample_name), names_to = "parameter", values_to = "value") %>%
  group_by(group) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(RADAR_EAD_TAG) - 2)) %>%
  unnest(cols = c(coords))

rescaled_data_MS_ARIST<- RADAR_MS_ARIST %>% 
  # mutate(across(-group, rescale)) %>%
  mutate(copy = pull(., 3)) %>% 
  pivot_longer(-c(group, Sample_name), names_to = "parameter", values_to = "value") %>%
  group_by(group) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(RADAR_MS_TAG) - 2)) %>%
  unnest(cols = c(coords))

rescaled_data_EAD_ALL<- RADAR_EAD_ALL %>% 
  # mutate(across(-group, rescale)) %>%
  mutate(copy = pull(., 2)) %>% 
  pivot_longer(-Sample_name, names_to = "parameter", values_to = "value") %>%
  group_by(Sample_name) %>%
  mutate(coords = rescaled_coords(value + central_distance, ncol(RADAR_EAD_ALL) - 1)) %>%
  unnest(cols = c(coords))



Comps_all_step3 <- step_2 + 
  geom_text(data = labels_axis_EAD, aes(x, 
                                        y, 
                                        label = c(round(x = c(0, (0.25*max_EAD_bla), (.5*max_EAD_bla), (.75*max_EAD_bla), max_EAD_bla), digits = 2)),
                                        # round(x = c(0, (0.25*bla), (.05*bla), (.75*bla), bla), digits = 2)),
                                        # round(x = c(0, (0.25*bla), (.05*bla), (.75*bla), bla), digits = 2)),
                                        
  ),colour = "#252525",
  hjust = -1) +
  annotate("text", x=0, y=-1.5, label= "EAD (mV)", colour = "#252525") +
  geom_text(data = text_coords(1 + central_distance + 0.2), 
            aes(x, y, color = Superclass_for_coloring), 
            label = labels_data$parameter[1:(ncol(RADAR_EAD_TAG)-2)]) +
  scale_colour_manual(values=c("#034e7b",
                               "#005a32",
                               "#737373",
                               "#ef3b2c", 
                               "#41ab5d",
                               "#99000d", 
                               "#4a1486",
                               "#fc9272"))

Comps_all_radar <- Comps_all_step3 +
  geom_point(data = rescaled_data_EAD_ALL, 
             aes(x, y,), 
             size = 1,
             colour = "#252525") +
  geom_path(data = rescaled_data_EAD_ALL, 
            aes(x, y,), 
            size = 1,
            color = "#252525")

ggsave(filename = "Output/Compounds_average_chainlenght.pdf", width = 250, height = 200, units = "mm")


Antenna_radar <- step_3 +
  geom_point(data = rescaled_data_EAD_TAG, 
             aes(x, y, group = group, col = group,), 
             size = 1,
             colour = "orange") +
  geom_path(data = rescaled_data_EAD_TAG, 
            aes(x, y, group = group, col = group,), 
            size = 1,
            color = "orange") +
  geom_point(data = rescaled_data_EAD_AIL, 
             aes(x, y, group = group, col = group,), 
             size = 1,
             colour = "green") +
  geom_path(data = rescaled_data_EAD_AIL, 
            aes(x, y, group = group, col = group,), 
            size = 1,
            color = "green") +
  geom_point(data = rescaled_data_EAD_ARIST, 
             aes(x, y, group = group, col = group,), 
             size = 1,
             colour = "blue") +
  geom_path(data = rescaled_data_EAD_ARIST, 
            aes(x, y, group = group, col = group,), 
            size = 1,
            color = "blue")


Mass_spec_radar <- step_3 +
  geom_point(data = rescaled_data_MS_TAG, 
             aes(x, y, group = group, col = group), 
             size = 1,
             color = "orange") +
  geom_path(data = rescaled_data_MS_TAG, 
            aes(x, y, group = group, col = group), 
            size = 1,
            color = "orange") +
  geom_point(data = rescaled_data_MS_AIL, 
             aes(x, y, group = group, col = group), 
             size = 1,
             color = "green") +
  geom_path(data = rescaled_data_MS_AIL, 
            aes(x, y, group = group, col = group), 
            size = 1,
            color = "green") +
  geom_point(data = rescaled_data_MS_ARIST, 
             aes(x, y, group = group, col = group), 
             size = 1,
             color = "blue") +
  geom_path(data = rescaled_data_MS_ARIST, 
            aes(x, y, group = group, col = group), 
            size = 1,
            color = "blue")
  



Tagetes_radar <- step_3 + 
  geom_point(data = rescaled_data_EAD_TAG, 
             aes(x, y, group = group, col = group,), 
             size = 1,
             colour = "#b10026") +
  geom_path(data = rescaled_data_EAD_TAG, 
            aes(x, y, group = group, col = group,), 
            size = 1,
            color = "#b10026") +
  geom_point(data = rescaled_data_MS_TAG, 
             aes(x, y, group = group, col = group), 
             size = 1,
             color = "#fd8d3c") +
  geom_path(data = rescaled_data_MS_TAG, 
            aes(x, y, group = group, col = group), 
            size = 1,
            color = "#fd8d3c") + 
  labs(title = "Tagetes\n") +
  theme(
    plot.title = element_text(size = 20, face = "bold", color = "darkgreen")
    )
  
Tagetes_radar %>% ggsave(filename = "Output/Tagetes_radar_chainlenght.pdf", width = 250, height = 200, units = "mm")

Ailanthus_radar <- step_3 + 
  geom_point(data = rescaled_data_EAD_AIL, 
             aes(x, y, group = group, col = group,), 
             size = 1,
             colour = "#b10026") +
  geom_path(data = rescaled_data_EAD_AIL, 
            aes(x, y, group = group, col = group,), 
            size = 1,
            colour = "#b10026") +
  geom_point(data = rescaled_data_MS_AIL, 
             aes(x, y, group = group, col = group), 
             size = 1,
             color = "#fd8d3c") +
  geom_path(data = rescaled_data_MS_AIL, 
            aes(x, y, group = group, col = group), 
            size = 1,
            color = "#fd8d3c") +
  labs(title = "Ailanthus\n") +
  theme(
    plot.title = element_text(size = 20, face = "bold", color = "darkgreen")
  )

Ailanthus_radar %>% ggsave(filename = "Output/Ailanthus_radar_chainlenght.pdf", width = 250, height = 200, units = "mm")

Aristolochia_radar <- step_3 + 
  geom_point(data = rescaled_data_EAD_ARIST, 
             aes(x, y, group = group, col = group,), 
             size = 1,
             colour = "#b10026") +
  geom_path(data = rescaled_data_EAD_ARIST, 
            aes(x, y, group = group, col = group,), 
            size = 1,
            colour = "#b10026") +
  geom_point(data = rescaled_data_MS_ARIST, 
             aes(x, y, group = group, col = group), 
             size = 1,
             color = "#fd8d3c") +
  geom_path(data = rescaled_data_MS_ARIST, 
            aes(x, y, group = group, col = group), 
            size = 1,
            color = "#fd8d3c")+
  labs(title = "Aristolochia\n") +
  theme(
    plot.title = element_text(size = 20, face = "bold", color = "darkgreen")
  )

Aristolochia_radar %>% 
  ggsave(filename = "Output/Aristolochia_radar_chainlenght.pdf", width = 250, height = 200, units = "mm")

step_5 <- Tagetes_radar + 
  labs(col = "tagetes") +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12))





Sum_Aristolochia <- MS_Aristolochia %>% 
  group_by(`Sample descriptor`) %>% 
  summarise(SUM_area = sum(Area.from.MS)) %>% 
  ungroup() 
MS_Aristolochia1_ <- MS_Aristolochia %>% 
  left_join(Sum_Aristolochia) %>%
  mutate(Compound.Name = tolower(Compound.Name)) %>%
  left_join(MS_corrected) %>%
  mutate(Compound.Name = ifelse(is.na(corrected_name),
                                Compound.Name,
                                corrected_name)) %>%
  mutate(percent = Area.from.MS*100/SUM_area) %>% 
  left_join(Grouped_compounds[42:48]) %>% 
  group_by(Superclass_results1) %>% 
  summarise(n = n())
Sum_Ailanthus <- MS_Ailanthus %>% 
  group_by(`Sample descriptor`) %>% 
  summarise(SUM_area = sum(Area.from.MS)) %>% 
  ungroup() 
MS_Ailanthus1_ <- MS_Ailanthus %>% 
  left_join(Sum_Ailanthus) %>% 
  mutate(Compound.Name = tolower(Compound.Name)) %>%
  left_join(MS_corrected) %>%
  mutate(Compound.Name = ifelse(is.na(corrected_name),
                                Compound.Name,
                                corrected_name)) %>%
  mutate(percent = Area.from.MS*100/SUM_area) %>% 
  left_join(Grouped_compounds[42:48]) %>% 
  group_by(Superclass_results1) %>% 
  summarise(n = n())

Sum_Tagetes <- MS_Tagetes %>% 
  group_by(`Sample descriptor`) %>% 
  summarise(SUM_area = sum(Area.from.MS)) %>% 
  ungroup() 
MS_Tagetes1_ <- MS_Tagetes %>% 
  left_join(Sum_Tagetes) %>% 
  mutate(Compound.Name = tolower(Compound.Name)) %>%
  left_join(MS_corrected) %>%
  mutate(Compound.Name = ifelse(is.na(corrected_name),
                                Compound.Name,
                                corrected_name)) %>%
  mutate(percent = Area.from.MS*100/SUM_area) %>% 
  left_join(Grouped_compounds[42:48]) %>% 
  group_by(Superclass_results1) %>% 
  summarise(n = n())

Supp_Table_1 <- All_MS %>% 
  left_join(Grouped_compounds[42:48]) %>%
  filter(!is.na(`Published.kovats(closest)`)) %>% 
  dplyr::select(Compound.Name, 'Sample descriptor', Area.from.MS, Superclass_results1, Np_pathway1, 'Published.kovats(closest)') %>% 
  pivot_wider( names_from = 'Sample descriptor', values_from = Area.from.MS) %>% 
  left_join(FID_ng, by = "Compound.Name") %>% 
  filter(!is.na(Superclass_results1))

write.xlsx(Supp_Table_1, "Output/Supp_table1.xlsx")  














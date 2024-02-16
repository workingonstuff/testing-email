# AUTHOR: Jenny Kay with support from Anaga Dinesh
# PURPOSE: Adapt Kristin Knox's NHANES code for mapping BCRC list to NHANES urine and blood analytes
# STARTED: 2023-08-28
# written in version: R version 4.2.2 (2022-10-31 ucrt)

# Use Zach Stanfield's metabolite mapping 

library(tidyverse)
library(readxl)
library(devtools)
library(openxlsx) #just for writing excel file output at end
install_github("silentspringinstitute/RNHANES") # Download the Git version of RNHANES

# Set the working directory
workingdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workingdir)

options(stringsAsFactors = FALSE)
options(scipen = 999)

#### Read in files for mapping BCRCs to NHANES analytes ####
urinemap <- read_xlsx("inputs/NHANEScodes_urine.xlsx", sheet=3)%>%
  select(-Added) %>% 
  rename(CASRN = CAS, CASRN.1 = CAS.1) %>% 
  #metsulfuron-methyl was added twice: Victoria used current CASRN 74223-64-6, Woody put deprecated 82197-07-7, delete deprecated entry
  filter(CASRN != "82197-07-7")

bloodmap <- read_xlsx("inputs/NHANEScodes_BloodSerum.xlsx") %>% 
  #blood mapping file does not include ethylene oxide, but it's a BCRC, so add that too
  #   but only 2017-2018 cycle because I'm only going to use up to that cycle
  rbind(c("DTXSID0020600", "Ethylene oxide", "75-21-8", "Red blood cells", "Ethylene oxide", 
          "LBXEOA", "ETHOX_J", "Ethylene Oxide", "2017-2018", "pmol/g_Hb"))

bcrc <- read_csv("inputs/BCRelList.csv") %>%
  rename(BCRC_name = preferred_name) %>%
  select(CASRN, DTXSID, BCRC_name)

#### Identify BCRC urine analytes as parents or metabolites ####
# code adapted from Stanfield: it classifies every chemical in the mapping as exclusively parent, exclusively metabolite,
# or chemical directly measured in urine (could have both parents and metabolites)

# parentchem takes the parent chemicals and their IDs
parentchem <- urinemap %>%
  select(Name, CASRN, DTXSID) %>%
  mutate(type = "Parent") %>%
  distinct()

#chemanalyte takes the chemicals measured in urine
chemanalyte <- urinemap %>%
  select(Name.1, CASRN.1, DTXSID.1) %>%
  rename(Name=Name.1, CASRN=CASRN.1, DTXSID=DTXSID.1) %>%
  mutate(type2 = "Metabolite") %>%
  distinct()

#list of unique parents and measured chems
chemlist <- full_join(parentchem, chemanalyte, by = c("Name", "CASRN", "DTXSID")) %>%
  mutate(type = case_when(!is.na(type) & !is.na(type2) ~ "Directly measured",
                          !is.na(type) ~ type, 
                          !is.na(type2) ~ type2)) %>%
  select(-type2)
  
bcrc_NHANESurine_chems <- inner_join(bcrc, chemlist,  by="CASRN")

# there are 11 chemicals with type = "Metabolite" - in Stanfield terminology that means
# these are metabolites; but, they themselves are on the BCRC list, so these I think we just match directly
bcrc_urine_metabolites <- bcrc_NHANESurine_chems %>%
  filter(type=="Metabolite") %>%
  rename(DTXSID=DTXSID.x) %>%
  mutate(DTXSID_analyte = DTXSID.y, CASRN_analyte = Name) %>%
  select(CASRN, DTXSID, BCRC_name, CASRN_analyte, DTXSID_analyte, type) #for these cases, DTXSID=DTSXID_metabolite
#KK note: I checked the BCRC_name columns; they are all on the BCRC list as preferred names; so this sub-list is all set to go


# there are 21 chemicals that have type = "Directly measured" - these could have both parents and metabolites - 
# but for us the relevant info is that these are in NHANES so I think we can treat them like the "Metabolite"
# Examples: atrazine is measured directly in NHANES in urine, but there is also an atrazine metabolite, Atrazine mercapturate,
# that is also included in urine; But ethyl paraben and genistein are only measured directly
# So maybe I should just check all of these
#check from here-------------------------------------------------------------------
directly_measured_cas <- bcrc_NHANESurine_chems %>%
  filter(type=="Directly measured")# %>%
  select(CASRN)

# see if there are BCRCs that are metabolites measured in NHANES, but the corresponding 
#    parents are not on the BCRC list
directly_measured_check <- urinemap %>%
  mutate(parent_flag = case_when(CASRN %in% directly_measured_cas$CASRN ~ 1, 
                                 TRUE ~ 0),
         metabolite_flag = case_when(CASRN.1 %in% directly_measured_cas$CASRN ~ 1,
                                     TRUE ~ 0)) %>%
  filter(parent_flag==1 | metabolite_flag==1) %>%
  left_join(bcrc, by="CASRN") %>% 
  select(-DTXSID.y) %>% 
  rename(DTXSID = DTXSID.x)
# the cases in above that have NA for BCRC_name are because the parent is not on the BCRC list but the metabolite is
# there's only one BCRC chemical on this list that this is true for: NHANES analyte pentachlorophenol


directly_measured_chems <- directly_measured_check %>%
  filter(is.na(BCRC_name) == FALSE) %>% 
  #parent_flag==1 and metabolite_flag==1 - then Name is both the bcrc and the nhanes analyte
  #parent_flag==1 and metabolite_flag==0 - then Name is the bcrc and Name.1 is the nhanes
  #parent_flag==0 and metabolite_flag==1 - then Name.1 is the bcrc and Name.1 is the nhanes
  mutate(NHANES_analyte=case_when(parent_flag==1 & metabolite_flag==1 ~ Name,
                                  parent_flag==1 & metabolite_flag==0 ~ Name.1,
                                  parent_flag==0 & metabolite_flag==1 ~ Name.1),
         DTXSID_analyte=case_when(parent_flag==1 & metabolite_flag==1 ~ DTXSID,
                                  parent_flag==1 & metabolite_flag==0 ~ DTXSID.1,
                                  parent_flag==0 & metabolite_flag==1 ~ DTXSID.1),
         CASRN_analyte=case_when(parent_flag==1 & metabolite_flag==1 ~ CASRN,
                               parent_flag==1 & metabolite_flag==0 ~ CASRN.1,
                               parent_flag==0 & metabolite_flag==1 ~ CASRN.1)) %>%
  mutate(type="Directly measured") %>%
  select(CASRN, DTXSID, BCRC_name, NHANES_analyte, CASRN_analyte, DTXSID_analyte, type)

#sort by BCRC_name to look at this; it's a small number of cases
# 1. Atrazine - directly measured, but also has metabolite Atrazine mercapturate
# 2. 7,4'-Dihydroxyisoflavone (aka Daidzein) - directly measured, 
#       but also 2 metabolites O-Desmethylangolensin and Equol
# 3. Pentachlorophenol - on BCRC list, directly measured, also 2 metabolites: 
#       2,4,6-Trichlorophenol and 2,4,5-Trichlorophenol in NHANES, and also a metabolite 
#       of hexachlorobenzene and lindane (not on bcrc list), and also beta-Hexachlorocyclohexane 
#       which is on the BCRC list



# there are 54 chemicals that have type="Parent" - for these, I need to match them back up with their
# NHANES metabolite(s); there might be more than one
parent_chemical_cas <- bcrc_NHANESurine_chems %>%
  filter(type=="Parent") %>%
  select(CASRN)

parent_chemicals <- urinemap %>%
  filter(CASRN %in% parent_chemical_cas$CASRN) %>% 
  # In this list, the parents are all BCRCs, they themselves are not in NHANES, and
  #    they are matched to their metabolites in NHANES (in some cases, they have 
  #    more than 1 metabolite; also in some cases one metabolite is linked to multiple parents)
  rename(NHANES_analyte=Name.1,
         DTXSID_analyte=DTXSID.1,
         CASRN_analyte=CASRN.1) %>%
  mutate(type="Parent") %>%
  select(CASRN, DTXSID, Name, NHANES_analyte, CASRN_analyte, DTXSID_analyte, type) %>% 
  # merge in the preferred names from the BCRC list
  left_join(bcrc, by="CASRN") %>% #checked, none of the BCRC names are missing
  select(-DTXSID.y) %>% 
  rename(DTXSID = DTXSID.x) %>% 
  select(CASRN, DTXSID, BCRC_name, NHANES_analyte, CASRN_analyte, DTXSID_analyte, type)

#unique(parent_chemicals$BCRC_name) #checked; there are 54

# So now we need to put the three pieces back together, so that we have one list that lists all of the 
# analytes we are pulling from NHANES, linked to their corresponding BCRC chemicals
full_chemlist <- rbind(bcrc_urine_metabolites, directly_measured_chems, parent_chemicals) %>%
  select(-type) %>%
  distinct() %>% 
  # if we stopped at previous line and ran unique(full_chemlist$BCRC_name), would find 86 unique BCRCs
  # if we ran unique(full_chemlist$NHANES_analyte), would get 74 unique NHANES metabolites
  #make a couple of indicators just so I can keep things straight
  mutate(parent_is_bcrc = case_when(CASRN %in% bcrc$CASRN ~ 1, TRUE ~ 0),
         analyte_is_BCRC = case_when(CASRN_analyte %in% bcrc$CASRN ~ 1, TRUE ~ 0)) 



#### match urine chems to NHANES data files ####
# Read in NHANES variable codes and data file names
urinecodes <- read_xlsx("inputs/NHANEScodes_urine.xlsx", sheet=1) %>%
  select(-Added.By, -Media, -MaxAge) %>% 
  #fix some incorrect CASRN#s KK identified (below)
  rename(CASRN = CAS) %>% 
  mutate(CASRN = case_when(Name == "Enterolactone" ~ "78473-71-9",
                         Name == "Triasulfuron" ~ "82097-50-5", 
                         TRUE ~ CASRN))%>% 
  mutate(NHANESfile = ifelse(NHANESfile == "l26PP_B", "L26PP_B", NHANESfile))  

# 4 that did not match:
#Enterolactone; 78473-71-9; DTXSID0048183 
#Triasulfuron; 82097-50-5; DTXSID0024345
#NNAL (4-(methylnitrosamino)-1-(3-pyridyl)-1-butanol); 76014-81-8; DTXSID8020880; not present in codes file
#Diphenyl phosphate; 838-85-7; DTXSID1048207; not present in codes file


#I think I'm missing some things... ethylene oxide, 

# sheet 2 has the names of the relevant MEC weight variables 
urineweights <-read_xlsx("inputs/NHANEScodes_urine.xlsx", sheet=2) %>% 
  rename(NHANESfile = file, weights_column = wtvariable) %>% 
  mutate(NHANESfile = gsub(".XPT", "", NHANESfile),
         NHANESfile = ifelse(NHANESfile == "l26PP_B", "L26PP_B", NHANESfile)) %>% 
  select(NHANESfile, weights_column)


# merge bcrcs with data file names
bcrc_urinefiles <- full_chemlist %>%
  rename(BCRC_CASRN = CASRN, CASRN=CASRN_analyte) %>% 
  select(BCRC_CASRN, BCRC_name, NHANES_analyte, CASRN:analyte_is_BCRC) %>%
  left_join(urinecodes, by="CASRN") %>% 
  left_join(urineweights, by = "NHANESfile") %>% 
  select(-c(DTXSID, Name)) %>% 
  rename(CASRN_analyte = CASRN, cycle = recent_sample) %>% 
  #filter out NNAL and diphenyl phosphate since they're not in codes file
  filter(!is.na(NHANEScode)) 


#most recent data for each analyte
bcrc_mostrecentdata <- bcrc_urinefiles %>% 
  # make a variable for when data was collected so I can pick the most recent for each
  mutate(cyclenumber = case_when(cycle == "99-00" ~ 1,
                                 cycle == "01-02" ~ 2, 
                                 cycle == "03-04" ~ 3,
                                 cycle == "05-06" ~ 4, 
                                 cycle == "07-08" ~ 5,
                                 cycle == "09-10" ~ 6, 
                                 cycle == "11-12" ~ 7,
                                 cycle == "13-14" ~ 8, 
                                 cycle == "15-16" ~ 9),
         NHANESfile = ifelse(NHANESfile =="l26PP_B", "L26PP_B", NHANESfile)) %>% 
  group_by(NHANES_analyte) %>% 
  filter(cyclenumber == max(cyclenumber)) %>% 
  ungroup() %>% 
  select(-cyclenumber)
#yay - we have most recent cycle for each analyte matched to the corresponding BCRC
#note that some analytes map to multiple BCRCs, and I'm sure many analytes map to other parents

#save this output for later
write_csv(bcrc_mostrecentdata, "outputs/mostrecenturineBCRCs_test.csv")



#work back from here and maybe smooth out the OG steps to input


#^ BCRCS in NHANES and most recent cycles ------------------------------------------------------------------------------------

### Gather relevant urine data for BCRCs ####

bcrc_mostrecentdata<-read.csv("./outputs/mostrecenturineBCRCs.csv")

NHANESanalytes <- bcrc_mostrecentdata %>% 
  #comment column basically just takes the analyte code (URX___) and
  #   replaces X with D and adds LC (URD___LC)
  #   except for PHIP and Trp-P-2 in 13-14 cycle 
  mutate(comment_column = case_when(NHANEScode == "URXPHIP" ~ "URDPHPLC",
                                    NHANEScode == "URXTRP2" ~ "URDTP2LC",
                                    TRUE ~ paste0("URD", substring(NHANEScode, first = 4), "LC"))) %>% 
  #rename columns so they match what the nhanes_quantile function uses for inputs
  rename(file_name = NHANESfile, column = NHANEScode) %>% 
  select(NHANES_analyte, file_name, column, comment_column, weights_column) %>% 
  unique() 


# find all the unique NHANES file names and their corresponding cycles
urinefilestoload <- bcrc_mostrecentdata %>% 
  select(NHANESfile, cycle) %>% 
  unique()

# 19 files to load, let's do in order of cycle year
nhanes_load_data("EPHPP_G", "2011-2012")
## 15-16 personal care/consumer products
EPHPP_I <- nhanes_load_data("EPHPP_I", "2015-2016", demographics = TRUE)

## 15-16 phthalates and plasticizers
PHTHTE_I <- nhanes_load_data("PHTHTE_I", "2015-2016", demographics = TRUE)

## 15-16 PAHs
PAH_I <- nhanes_load_data("PAH_I", "2015-2016", demographics = TRUE)

## 15-16 VOCs
UVOC_I <- nhanes_load_data("UVOC_I", "2015-2016", demographics = TRUE)

## 13-14 Pyrethroids, Herbicides, & Organophosphorus Metabolites 
UPHOPM_H <- nhanes_load_data("UPHOPM_H", "2013-2014", demographics = TRUE)

## 13-14 heterocyclic amines
HCAA_H <- nhanes_load_data("HCAA_H", "2013-2014", demographics = TRUE)

# 11-12 organophsphate insecticides
OPD_G <- nhanes_load_data("OPD_G", "2011-2012", demographics = TRUE)

# 11-12 phthalates - most recent cycle for DMP (metab MMP)
PHTHTE_G <- nhanes_load_data("PHTHTE_G", "2011-2012", demographics = TRUE)

# 09-10 pesticides - environmental
PP_F <- nhanes_load_data("PP_F", "2009-2010", demographics = TRUE)

# 09-10 environmental phenols
EPH_F <- nhanes_load_data("EPH_F", "2009-2010", demographics = TRUE)

# 09-10 phytoestrogens
PHYTO_F <- nhanes_load_data("PHYTO_F", "2009-2010", demographics = TRUE)

# 09-10 Pyrethroids, Herbicides, & Organophosphorus Metabolites 
UPHOPM_F <- nhanes_load_data("UPHOPM_F", "2009-2010", demographics = TRUE)

# 09-10 phthalates 
PHTHTE_F <- nhanes_load_data("PHTHTE_F", "2009-2010", demographics = TRUE)

# 07-08 "current use" pesticides
UPP_E <- nhanes_load_data("UPP_E", "2007-2008", demographics = TRUE)

# 07-08 atrazine and metabolites
UAM_E <- nhanes_load_data("UAM_E", "2007-2008", demographics = TRUE)

# 07-08 pesticides - carbamates & organophosphorous metabs
CARB_E <- nhanes_load_data("CARB_E", "2007-2008", demographics = TRUE)

# 03-04 "current use" pesticides
L26UPP_C <- nhanes_load_data("L26UPP_C", "2003-2004", demographics = TRUE)

# 01-02 "current use" pesticides
# note that unlike others, 99-00 and 01-02 pesticide samples have special 
#    2- and 4-year weights depending on whether the cycles are grouped. 
#    I'm only pulling from one cycle at a time, so using the 2-yr weights
L26PP_B <- nhanes_load_data("L26PP_B", "2001-2002", demographics = TRUE)

# 99-00 "current use" pesticides
LAB26PP <- nhanes_load_data("LAB26PP", "1999-2000", demographics = TRUE)

#probs make this a function-- two lists

### Compute quantiles and detection frequencies with pmap ###

#first create the input lists for calculating quantiles and detection frequencies
#"get" creates a vector of the df column without adding quotes around entries
#    (needed for nhanes_quantile function to run properly)
nhanes_data <- as.list(map(NHANESanalytes$file_name, get))
column <- as.list(NHANESanalytes$column)
comments <- as.list(NHANESanalytes$comment_column)
weights <- as.list(NHANESanalytes$weights_column)


#now calc median and 95th %ile
map_quantiles_list <- pmap(list(nhanes_data, column, comments, weights),
                             .f = nhanes_quantile, quantiles = c(0.5, 0.95))

mapped_quantiles <- do.call(rbind, map_quantiles_list) %>% 
  pivot_wider(names_from = quantile, values_from = c(value, below_lod)) %>%
  rename(NHANESfile = file_name, NHANEScode = column, most_recent_cycle = cycle, 
         median = `value_50%`, percentile_95th = `value_95%`) %>% 
  mutate(median = ifelse(`below_lod_50%` == TRUE, "belowLOD", median),
         percentile_95th = ifelse(`below_lod_95%` == TRUE, "belowLOD", percentile_95th),
         media = "Urine") %>% 
  full_join(bcrc_mostrecentdata, by = c("NHANESfile", "NHANEScode")) %>% 
  select(BCRC_CASRN:CASRN_analyte, analyte_is_BCRC, most_recent_cycle,
         median, percentile_95th, units, media)

#now calc detection frequencies
map_detectfreq_list <- pmap(list(nhanes_data, column, comments, weights),
                           .f = nhanes_detection_frequency)

mapped_detectfreqs <- do.call(rbind, map_detectfreq_list) %>% 
  rename(NHANESfile = file_name, NHANEScode = column, most_recent_cycle = cycle, 
         detectionfrequency = value) %>% 
  mutate(detectionfrequency = round(detectionfrequency, digits = 5)) %>% 
  full_join(bcrc_mostrecentdata, by = c("NHANESfile", "NHANEScode")) %>% 
  select(BCRC_CASRN:CASRN_analyte, analyte_is_BCRC, most_recent_cycle,
         detectionfrequency)


# YAY collect quantiles and detection frequencies - done!
BCRC_NHANES_urine <- full_join(mapped_quantiles, mapped_detectfreqs)


# write_csv(BCRC_NHANES_urine, "outputs/BCRC_NHANES_urinedata.csv")



#### BCRCs measured in blood/serum ####
#Mapping from Z. Stanfield from EPA indicates that blood only includes BCRC parents (no metabolites) 

#find the BCRCs measured in blood 
bcrc_blood_files <- left_join(bloodmap, bcrc, by = c("CASRN", "DTXSID")) %>% 
  filter(!is.na(BCRC_name)) %>% 
  # make a variable for when data was collected so I can pick the most recent for each
  mutate(cyclenumber = case_when(cohort == "1999-2000" ~ 1,
                                 cohort == "2001-2002" ~ 2, 
                                 cohort == "2003-2004" ~ 3,
                                 cohort == "2005-2006" ~ 4, 
                                 cohort == "2007-2008" ~ 5,
                                 cohort == "2009-2010" ~ 6, 
                                 cohort == "2011-2012" ~ 7,
                                 cohort == "2013-2014" ~ 8, 
                                 cohort == "2015-2016" ~ 9,
                                 cohort == "2017-2018" ~ 10,
                                 TRUE ~ 0)) %>% 
  group_by(BCRC_name) %>% 
  #filter out everything except most recent measurements
  filter(cyclenumber == max(cyclenumber)) %>% 
  ungroup() %>% 
  select(CASRN, BCRC_name, cohort, Media, units, data_file_name, variable_name) %>% 
  #rename things so they're what the nhanes_quantile function likes to see
  rename(file_name = data_file_name, column = variable_name) %>% 
  #add comment and weights codes
  mutate(comment_column = ifelse(file_name == "ALD_H", "LBD7ALC",
                                 paste0("LBD", substring(column, first = 4), "LC")),
         weights_column = case_when(file_name == "ALD_H" ~ "WTALD2YR", 
                                    file_name == "ETHOX_J" ~ "WTSA2YR",
                                    TRUE ~ "WTSVOC2Y"))


unique(bcrc_blood_files$file_name)
#6 files to load - "VOCWB_J"  "VOCWB_G"  "VOCWB_E"  "VOCWB_I"  "ALD_H"  "VOCMWB_G"
# all VOCs except one aldehyde


#in order of cycle
VOCWB_J <- nhanes_load_data("VOCWB_J", "2017-2018", demographics = TRUE)
#eth ox
ETHOX_J <- nhanes_load_data("ETHOX_J", "2017-2018", demographics = TRUE)

VOCWB_I <- nhanes_load_data("VOCWB_I", "2015-2016", demographics = TRUE)

ALD_H <- nhanes_load_data("ALD_H", "2013-2014", demographics = TRUE)

VOCWB_G <- nhanes_load_data("VOCWB_G", "2011-2012", demographics = TRUE)

VOCMWB_G <- nhanes_load_data("VOCMWB_G", "2011-2012", demographics = TRUE)

VOCWB_E <- nhanes_load_data("VOCWB_E", "2007-2008", demographics = TRUE)



### Compute quantiles and detection frequencies with pmap ###
blood_data <- as.list(map(bcrc_blood_files$file_name, get))
blood_column <- as.list(bcrc_blood_files$column)
blood_comments <- as.list(bcrc_blood_files$comment_column)
blood_weights <- as.list(bcrc_blood_files$weights_column)


## Median and 95th %ile
blood_quantiles_list <- pmap(list(blood_data, blood_column, blood_comments, blood_weights),
                             .f = nhanes_quantile, quantiles = c(0.5, 0.95))

blood_quantiles <- do.call(rbind, blood_quantiles_list) %>% 
  pivot_wider(names_from = quantile, values_from = c(value, below_lod)) %>%
  rename(most_recent_cycle = cycle, median = `value_50%`, percentile_95th = `value_95%`) %>% 
  mutate(median = ifelse(`below_lod_50%` == TRUE, "belowLOD", median),
         percentile_95th = ifelse(`below_lod_95%` == TRUE, "belowLOD", percentile_95th)) %>% 
  full_join(bcrc_blood_files) %>% 
  select(CASRN, BCRC_name, most_recent_cycle, Media, median, percentile_95th, units)


## Detection frequencies
blood_detectfreq_list <- pmap(list(blood_data, blood_column, blood_comments, blood_weights),
                              .f = nhanes_detection_frequency)

blood_detectfreqs <- do.call(rbind, blood_detectfreq_list) %>% 
  rename(most_recent_cycle = cycle, detectionfrequency = value) %>% 
  mutate(detectionfrequency = round(detectionfrequency, digits = 5)) %>% 
  full_join(bcrc_blood_files, by = c("file_name", "column", "comment_column", "weights_column")) %>% 
  select(CASRN, detectionfrequency)


# YAY collect quantiles and detection frequencies - done!
BCRC_NHANES_blood <- full_join(blood_quantiles, blood_detectfreqs) %>% 
  rename(BCRC_CASRN = CASRN, media = Media) %>% 
  select(BCRC_CASRN:detectionfrequency)
  


# write_csv(BCRC_NHANES_blood, "outputs/BCRC_NHANES_blood.csv")


#### combine urine and blood data and make output files ####
BCRC_NHANES_all <- full_join(BCRC_NHANES_urine, BCRC_NHANES_blood) %>% 
  mutate(NHANES_analyte = ifelse(is.na(NHANES_analyte), BCRC_name, NHANES_analyte),
         CASRN_analyte = ifelse(is.na(CASRN_analyte), BCRC_CASRN, CASRN_analyte),
         analyte_is_BCRC = ifelse(is.na(analyte_is_BCRC), "1", analyte_is_BCRC)) %>% 
  select(BCRC_CASRN:most_recent_cycle, media, detectionfrequency, median:units) %>% 
  arrange(BCRC_CASRN)


# write_csv(BCRC_NHANES_all, "outputs/BCRC_NHANES_analysis.csv")
# 
# write.xlsx(BCRC_NHANES_all, file = "outputs/BCRC_NHANES_analysis.xlsx")


# add in predicted median intake (upper 95% credible interval) from Ring 2019
BCrel_intakes <- read.csv("outputs/BCrel_Effects_and_Sources.csv") %>% 
  select(CASRN, MammaryTumorEvidence, EDC, Genotoxicity, U95intake_mg.kg.d) %>% 
  rename(BCRC_CASRN = CASRN, BCRC_U95intake_mg.kg.d = U95intake_mg.kg.d)


# create output of BCRC levels in people including NHANES measurements and predicted intakes 
BCRC_levels <- left_join(BCRC_NHANES_all, BCrel_intakes, by = "BCRC_CASRN") %>% 
  select(BCRC_CASRN, BCRC_name, MammaryTumorEvidence:BCRC_U95intake_mg.kg.d, NHANES_analyte:units)

write.xlsx(BCRC_levels, file = "outputs/BCRC_LevelsInPeople.xlsx")

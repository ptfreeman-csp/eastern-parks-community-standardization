###Pre-analysis standardization of taxa sampling to equal sampling depth by rarefaction###
setwd(dir = "Google Drive/My Drive/00_CSP_Projects/PaleoTransformation/00_PaleotransformationCSP/Paleotransformation_Update/")
library(tidyverse)
library(readr)
library(RRatepol)
library(vegan)

### Sites in long format 
sites <- read_csv("/Volumes/GoogleDrive/.shortcut-targets-by-id/1N-M82M43r2si2QN8hfEkTluNAAA2ThOH/2019.NC CASC Paleotransformation/RoC_Analysis/00_Data/00_CommunityData_inputs/eco_17_21_harmonized_filtered_pollen-20211201.csv")

### Split into list ###
sites.list <- split( sites , f = sites$dataset.id)

### For loop for standardizing all community datasets to a minimum sample count ###
for (i in 1:length(sites.list)) {
  
  ## Original long-format community data
  community.long <- sites.list[[i]] 
  
  ##Pull site name
  site.name <- unique(community.long$site.name)
  site.id <- unique(community.long$dataset.id)
  
  ##Determine the minimum sample count by summarizing counts at each depth
  min.sample.count <- community.long %>%
    group_by(depth) %>% 
    dplyr::summarise(Total.Count = sum(count, na.rm=T)) %>%
    ungroup() %>%
    dplyr::summarise(min = min(Total.Count))
  
  ## Which depths have the low counts 
  problem.rows <- community.long %>%
    group_by(depth) %>%
    dplyr::summarise(Total.Count = sum(count, na.rm=T)) %>%
    ungroup() %>%
    dplyr::filter(Total.Count < 95) %>% 
    pull(., depth)
  
  problem.counts <- community.long %>%
    group_by(depth) %>%
    dplyr::summarise(Total.Count = sum(count, na.rm=F)) %>%
    ungroup() %>%
    dplyr::filter(Total.Count < 95) %>% 
    pull(., Total.Count)
  
  ##Write out a warning message 
  if(min.sample.count < 95) 
    
    cat(paste0(site.name, site.id, " has ",problem.counts," grains at depth =", problem.rows, ".\n"), file=paste0("00_Data/01_SiteData_NonSpatial/StandardizedCommunityData/Standardized_RemovedStrata/", site.name, site.id, "_removedstrata.txt"))
  
  community.long.clean <- community.long %>% 
    dplyr::filter(!depth %in% problem.rows)
  
  ## The section below is used to write out files with site metadata for cluster analysis ### 
#  site_meta <- community.long.clean[,c(1:10)] %>% 
#    group_by(depth) %>% 
#    dplyr::slice(., 1) %>%
#    ungroup() %>%
#    dplyr::mutate(X = row_number(), .before = .id)
  
  if(nrow(community.long.clean) != nrow(community.long))
    cat(paste0(site.name, site.id, " problem rows removed. Let's roll yeehaw!", ".\n"))
  
  min.sample.count.clean <- community.long.clean %>%
    group_by(depth) %>% 
    dplyr::summarise(Total.Count = sum(count, na.rm=T)) %>%
    ungroup() %>%
    dplyr::summarise(min = min(Total.Count))
  
  ##Skip to next if the error occurs
  #if(min.sample.count.clean < 95) next
  
  community.wide <- community.long.clean %>% 
    group_by(depth, harmon_cat) %>%
    dplyr::mutate(count = ceiling(count)) %>%
    slice(., 1) %>% 
    pivot_wider(., names_from = "harmon_cat", 
                values_from = "count") %>%
    ungroup() %>% 
    dplyr::select(c(11:max(length(names(.))))) %>%
    dplyr::mutate(sample.id = row_number(), .before = 1)
  
  #### Extract ages ###
  ages <- community.long.clean %>% 
    group_by(depth) %>% 
    dplyr::select(age) %>% 
    dplyr::slice(., 1)
  
  ### Change path as necessary 
  age.path <- "/Volumes/GoogleDrive/My Drive/00_CSP_Projects/PaleoTransformation/00_PaleotransformationCSP/Paleotransformation_Update/00_Data/01_SiteData_NonSpatial/StandardizedCommunity_dates/"
  
  ### Write out to csv
  write.csv(ages, paste0(age.path,site.name,"_",site.id, "_standardized_dates.csv"), row.names=F)
  
  ### Rarefy community to 95 individuals per row (sample)
  
  community.prep <- community.wide %>% 
    dplyr::select(-sample.id) #[,-max(length(names(community.wide)))] 
 
   ### replace all NA with 0
  community.prep[is.na(community.prep)] = 0
  
  ## Perform rarefaction 
  community.rarefied <- (as.data.frame(rrarefy(community.prep, sample=95))) %>%
    dplyr::mutate(sample.id = community.wide$sample.id, .before=1)
  
  community.rarefied.summary <- community.rarefied %>%
    dplyr::select(-sample.id) %>%
    dplyr::rowwise() %>%
    dplyr::summarise(total = sum(c_across(where(is.numeric))))
  
  rarefied.problem.rows <- community.rarefied.summary %>%
    dplyr::filter(total < 95) %>% 
    dplyr::mutate(row = row_number()) %>%
    pull(., row)
  
  if(length(rarefied.problem.rows) == 0) 
    cat(paste0(site.name, site.id," has no problem rows in rarefied data",".\n"))
  
  if(length(rarefied.problem.rows) > 0 ) 
    cat(paste0(site.name, site.id," has an rarefaction issue due to ",  rarefied.problem.rows, ".\n"))
  
  #### Change path as necessary 
  
  path = "/Volumes/GoogleDrive/My Drive/00_CSP_Projects/PaleoTransformation/00_PaleotransformationCSP/Paleotransformation_Update/00_Data/01_SiteData_NonSpatial/StandardizedCommunityData/Standardized_for_Ratepol/"
  
  #### Append the metadata to match structure needed for clustering analysis 
  
  #community.rarefied <- bind_cols(site_meta, community.rarefied)
  
  ### Write out to csv
  write.csv(community.rarefied, paste0(path,site.name,"_",site.id, "_standardized.csv"), row.names=F)
  
}

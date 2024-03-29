#### This script is designed to prepare inputs for site-level rate-of-change analyses with the R-Ratepol package. It was originally performed by Patrick Freeman on a local machine within its own R project linked to a personal GitHub repo for easier file management and version control. 
#### It does the following things: 
#### 1. It ingests the harmonized pollen data tables prepared by Amanda Kissel that include information on the pollen samples at each depth in the soil core.
#### 2. It does an assessment of whether any rows have fewer than 100 total pollen grains, flags those rows and removes them. 
#### 3. It then performs a rarefaction procedure to resample all rows to equal depth with replacement. 
#### 4. Once that is done, it also gets the final set of standardized ages for the remaining depths (if any were removed above) and writes these to a prepared file in the format necessary for R-Ratepol.
#### 5. Finally, if any rows were removed in the filtering step (Step 2), it pulls in the site-level age uncertainty draws from the Bacon models run by Shelley Crausbay from their directory and filters them to remove those rows from the sampling dataframe and writes out a new file for use in R-Ratepol. To accomplish this, Patrick Freeman downloaded this directory of Bacon age-depth model outputs from the Google Drive folder provided by Shelley Crausbay to his local machine and this directory path should be updated if any updates to those models are completed. 

########## Script written by Patrick Freeman (patrick[at]csp-inc.org) ##########

setwd("/Users/patrickfreeman-csp/Documents/GitHub/eastern-parks-community-standardization")

### Load libraries
library(tidyverse)
library(vegan)

sites.list <- list.files("/Users/patrickfreeman-csp/Documents/GitHub/eastern-parks-community-standardization/data/inputs/site_pollen_updated_ages_harmonized", pattern=".csv", full.names=T)

### Read sites in
sites.read.list <- list()

for(i in 1:length(sites.list)){
  
  site <- read_csv(sites.list[[i]]) %>% 
    dplyr::mutate_at("dataset.id", as.character)
  
  sites.read.list[[i]] <- site
}

### Bind all together and pivot into long format 
sites.read.df <- bind_rows(sites.read.list)

### Split pollen dataframe by site
pollen_split <- split(sites.read.df, f=sites.read.df$dataset.id)

for(i in 1:length(pollen_split)){
  
  site.df <- pollen_split[[i]] %>% 
    ungroup() 
  dataset.id <- unique(site.df$dataset.id)
  site.name <- unique(site.df$site.name)
  
  cat(paste0(site.name, " is starting.", "\n"))
  
  ### pivot longer to get strata/taxa level counts 
  community.long <- site.df %>% 
    group_by(site.name, dataset.id, depth) %>%
    pivot_longer(., 
                 cols=10:tail(names(.), 1), 
                 names_to="harmon_taxa",
                 values_to="count")
  
  ##Determine the minimum sample count by summarizing counts at each depth
  min.sample.count <- community.long %>%
    group_by(depth) %>% 
    dplyr::summarise(Total.Count = sum(count, na.rm=T)) %>%
    ungroup() %>%
    dplyr::summarise(min = min(Total.Count)) %>%
    dplyr::pull(min)
  
  ## Which depths have the low counts 
  problem.rows <- community.long %>%
    group_by(depth) %>%
    dplyr::summarise(Total.Count = sum(count, na.rm=T)) %>%
    ungroup() %>%
    dplyr::filter(Total.Count < 100) %>% 
    pull(., depth)
  
  ## Pull the problem counts 
  problem.counts <- community.long %>%
    group_by(depth) %>%
    dplyr::summarise(Total.Count = sum(count, na.rm=T)) %>%
    ungroup() %>%
    dplyr::filter(Total.Count < 100) %>% 
    pull(., Total.Count)
  
  ##Write out a warning message 
  if(min.sample.count < 100) 
    
    cat(paste0(site.name," ", dataset.id, " has ",problem.counts," grains at depth = ", problem.rows, ".\n"))
  
  if(min.sample.count > 100)
    
    cat(paste0(site.name," ", dataset.id, " has >100 grains at all depths", ".\n"))
  
  ### prep community df
  community.df.prep <- site.df %>% 
    ungroup() %>%
    dplyr::filter(!depth %in% problem.rows) %>%
    replace(is.na(.), 0) %>%
    dplyr::mutate(sample.id = row_number()) %>% 
    dplyr::relocate(sample.id, .before=everything())
    
  ## Set up for rarefaction
  community.prep <- community.df.prep %>% 
    dplyr::select(-c(1:10)) %>% 
    dplyr::mutate_if(is.numeric, round) ## Round fractional individuals up to make integer counts
  
  ### Set seed to make sure sampling is the same
  set.seed(1234)
  ## Perform rarefaction
  community.rarefied <-
    (as.data.frame(rrarefy(community.prep, sample = 100))) %>%
    dplyr::mutate(sample.id = row_number(), .before = 1)

  ## Ensure that all rows now have 100
  community.rarefied.summary <- community.rarefied %>%
    dplyr::select(-sample.id) %>%
    dplyr::rowwise() %>%
    dplyr::summarise(total = sum(c_across(where(is.numeric))))
  
  ## Flag if problems
  rarefied.problem.rows <- community.rarefied.summary %>%
    dplyr::filter(total < 100) %>% 
    dplyr::mutate(row = row_number()) %>%
    pull(., row)
  
  if(length(rarefied.problem.rows) == 0) 
    cat(paste0(site.name, " ", dataset.id," has no problem rows in rarefied data",".\n"))
  
  if(length(rarefied.problem.rows) > 0 ) 
    cat(paste0(site.name, " ", dataset.id," has an rarefaction issue due to ",  rarefied.problem.rows, ".\n"))
  
  community.path = "/data/outputs/standardized_community/"
  
  ### Write out to csv for rarefied community data 
  write.csv(community.rarefied, paste0("/Users/patrickfreeman-csp/Documents/GitHub/eastern-parks-community-standardization",community.path,site.name,"_",dataset.id, "_standardized.csv"), row.names=F)

  if(min.sample.count >= 100) {
    
    #### Creating standardized age/depth and Bacon sampling draws files to feed to R-Ratepol
    
    ### Adjust the site name to match the way it's captured in the Bacon output folders
    site.name.bacon <- str_replace_all(site.name, " ", "_")
    
    ### list the files in the Bacon run directory
    bacon.list <- list.files("/Users/patrickfreeman-csp/Downloads/Bacon_runs 2", full.names=T)
    
    ### Convert to dataframe
    bacon.df <- data.frame(bacon_folder = bacon.list)
    
    ### Get path to site's folder
    target.bacon.folder.path <- bacon.df %>%
      dplyr::filter(str_detect(bacon_folder, dataset.id)) %>%
      dplyr::filter(str_detect(bacon_folder, site.name.bacon)) %>%
      pull(bacon_folder)
    
    ### Print as gut check 
    print(target.bacon.folder.path)
    
    ### Get list of files within site folder for both ages and sampled draws csvs 
    target.folder.file.list <- list.files(target.bacon.folder.path, full.names=T, pattern=".csv")
    target.folder.file.ages <- list.files(target.bacon.folder.path, full.names=T, pattern=".txt")
    
    ### Get path to ages file and pull file and convert to dataframe
    bacon.ages.path <- data.frame(file=target.folder.file.ages) %>% 
      dplyr::filter(str_detect(file, "_ages.txt")) %>% 
      dplyr::pull(file)
    
    ### Select depth and median age and filter to exclude the rows that have been pulled 
    bacon.ages.df <- read_delim(bacon.ages.path) %>% 
      dplyr::select(depth, median) %>% 
      dplyr::rename(age = median) %>% 
      #dplyr::filter(!depth %in% problem.rows) %>% 
      dplyr::select(-depth) %>% 
      dplyr::mutate(sample.id = row_number()) %>%
      dplyr::relocate(sample.id, .before=everything())
    
    age.path = "/data/outputs/standardized_ages/"
    
    ### Write out to csv for age data for standardized dataset with rows <100 grains removed
    write.csv(bacon.ages.df, paste0("/Users/patrickfreeman-csp/Documents/GitHub/eastern-parks-community-standardization",age.path,site.name,"_",dataset.id, "_ages_standardized.csv"), row.names=F)
    
    ### Convert Bacon samples file to dataframe
    target.folder.df <- data.frame(file=target.folder.file.list)
    
    ### Get path to sampled 1000 draws csv
    bacon.sampled.1000.path <- target.folder.df %>%
      dplyr::filter(str_detect(file, "_1000.csv")) %>%
      dplyr::pull(file)
    
    ### Read in the csv
    bacon.sampled.file <- read_csv(bacon.sampled.1000.path)
    
    ### Append the depth from the original dataframe to enable filtering
    bacon.sampled.file2 <- bacon.sampled.file %>%
      dplyr::mutate(depth = site.df$depth) %>%
      dplyr::relocate(depth, .before=everything())
    
    ### Filter out rows that had less than 100 grains 
    bacon.sampled.file.filtered <- bacon.sampled.file2 %>%
      dplyr::filter(!depth %in% problem.rows) %>% 
      dplyr::select(-depth)
    
    bacon.samples.path = "/data/outputs/standardized_bacon_samples"
    
    ### Write out filtered dataframe as csv
    write_csv(bacon.sampled.file.filtered, paste0("/Users/patrickfreeman-csp/Documents/GitHub/eastern-parks-community-standardization", bacon.samples.path,"/",site.name.bacon,"_",dataset.id,"_1000_filtered.csv"))
    
  }
  
  #### Adjusting the sampled bacon rows if needed
  if(min.sample.count < 100) {
    cat(paste0(site.name," ", dataset.id," starting the alignment of Bacon age/depth samples", ".\n"))
    
    ### Adjust the site name to match the way it's captured in the Bacon output folders
    site.name.bacon <- str_replace_all(site.name, " ", "_")
    
    ### list the files in the Bacon run directory
    bacon.list <- list.files("/Users/patrickfreeman-csp/Downloads/Bacon_runs 2", full.names=T)
    
    ### Convert to dataframe
    bacon.df <- data.frame(bacon_folder = bacon.list)
    
    ### Get path to site's folder
    target.bacon.folder.path <- bacon.df %>%
      dplyr::filter(str_detect(bacon_folder, dataset.id)) %>%
      dplyr::filter(str_detect(bacon_folder, site.name.bacon)) %>%
      pull(bacon_folder)
    
    ### Print as gut check 
    print(target.bacon.folder.path)
    
    ### Get list of files within site folder for both ages and sampled draws csvs 
    target.folder.file.list <- list.files(target.bacon.folder.path, full.names=T, pattern=".csv")
    target.folder.file.ages <- list.files(target.bacon.folder.path, full.names=T, pattern=".txt")
    
    ### Get path to ages file and pull file and convert to dataframe
    bacon.ages.path <- data.frame(file=target.folder.file.ages) %>% 
      dplyr::filter(str_detect(file, "_ages.txt")) %>% 
      dplyr::pull(file)
    
    ### Select depth and median age and filter to exclude the rows that have been pulled 
    bacon.ages.df <- read_delim(bacon.ages.path) %>% 
      dplyr::select(depth, median) %>% 
      dplyr::rename(age = median) %>% 
      dplyr::filter(!depth %in% problem.rows) %>% 
      dplyr::select(-depth) %>% 
      dplyr::mutate(sample.id = row_number()) %>%
      dplyr::relocate(sample.id, .before=everything())
    
    age.path = "/data/outputs/standardized_ages/"
    
    ### Write out to csv for age data for standardized dataset with rows <100 grains removed
    write.csv(bacon.ages.df, paste0("/Users/patrickfreeman-csp/Documents/GitHub/eastern-parks-community-standardization",age.path,site.name,"_",dataset.id, "_ages_standardized.csv"), row.names=F)
    
    ### Convert Bacon samples file to dataframe
    target.folder.df <- data.frame(file=target.folder.file.list)
    
    ### Get path to sampled 1000 draws csv
    bacon.sampled.1000.path <- target.folder.df %>%
      dplyr::filter(str_detect(file, "_1000.csv")) %>%
      dplyr::pull(file)
    
    ### Read in the csv
    bacon.sampled.file <- read_csv(bacon.sampled.1000.path)
    
    ### Append the depth from the original dataframe to enable filtering
    bacon.sampled.file2 <- bacon.sampled.file %>%
      dplyr::mutate(depth = site.df$depth) %>%
      dplyr::relocate(depth, .before=everything())
    
    ### Filter out rows that had less than 100 grains 
    bacon.sampled.file.filtered <- bacon.sampled.file2 %>%
      dplyr::filter(!depth %in% problem.rows) %>% 
      dplyr::select(-depth)
    
    bacon.samples.path = "/data/outputs/standardized_bacon_samples"
    
    ### Write out filtered dataframe as csv
    write_csv(bacon.sampled.file.filtered, paste0("/Users/patrickfreeman-csp/Documents/GitHub/eastern-parks-community-standardization", bacon.samples.path,"/",site.name.bacon,"_",dataset.id,"_1000_filtered.csv"))
  }
    
    
}



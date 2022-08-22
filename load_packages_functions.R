# 20201211, Dany Chauvin

# INSTALLING AND LOADING NECESSARY PACKAGES 

# Following packages are necessary to import data from deepMoma and analyze these data.

# WARNING: libstdc++.so.6 library, with proper GLIBCXX versions is necessary to use some packages imported here. 
# Be sure that GCCcore/8.3.0 is loaded/installed.
# To do so while using Rstudio on the scicore server "service06", do the following, add the following to your ~/.bashrc file:
# `if [[ "$HOSTNAME" = *service06* ]]; then
#     ml GCCcore/8.3.0
#  fi'

# Currently, packages are installed in my HOME folder, from centOS7 binaries at packagemanager.rstudio.com. I am not yet using renv, planned in the future.
# On Rstudio, set R version to 4.0.3.
# If this is the first time you use these codes, uncomment to following lines to install the corresponding packages.

# Uncomment below when using for the first time.
#options(repos = c(REPO_NAME = "https://packagemanager.rstudio.com/all/__linux__/centos7/latest"))
#install.packages("tidyverse")
#install.packages(c("here","cowplot"))
#install.packages("devtools")
#remotes::install_github(c('hadley/multidplyr'))
#install.packages("ccaPP")
#install.packages("~/R/vngMoM.tar", repos = NULL)
#remotes::install_github(c('julou/ggCustomTJ'))
#renv::init()
#install.packages("reticulate")
#install.packages("ggcorrplot")
#install.packages("lemon")
#install.packages("LSD")
#install.packages("parallel")
#install.packages("ggpubr")
# Comment above if the packages are already installed.

# Necessary libraries
library(tidyverse)
library(RcppArmadillo)
library(tools)
library(here)
library(cowplot)
library(devtools)
library(multidplyr)
library(renv)
library(svglite)
#Sys.setenv(RETICULATE_PYTHON = '/scicore/home/nimwegen/rocasu25/.local/share/r-miniconda/envs/r-reticulate/bin/python')
library(reticulate)
library(ggcorrplot)
library(lemon)
library(parallel)
library(broom)
library(stats)
library(ggpubr)
library(comprehenr)

read_Biotek_Synergy_stacker <- function(.path,.ch) {
  #.path
  #.path <- "./data/220818_20220816_acetate0.2.txt"
  .lines <- readLines(.path)
  .l_idx <- stringr::str_detect(.lines, "Time") %>% which %>% (function(.x) .x-1)
  #.ch <- "GFP"
  .data <- lapply(.l_idx, function(.i) {
    if(grepl(.ch,.lines[.i])){
      # extract timelapse step and timestamp
      .m <- stringr::str_match(.lines[.i+1], "Time (\\d+) \\((.*)\\)")
      .time <- .m[3]
      .step <- .m[2] # add time step (because different channels are recorded at different times)
      # prepend to rows and concatenate
      paste(.ch, .time, .step, .lines[(.i+3):(.i+10)], sep=',') %>% 
        paste(collapse='\n')}
    else{return(NA)}
  })
  .data <- .data[!is.na(.data)] %>% #get rid of rows belonging to other channels
    paste(collapse='\n') %>% 
    utils::read.table(text=., sep=",", header=FALSE, stringsAsFactors=FALSE) %>% 
    stats::setNames(c("channel", "time", "step", "row", 1:12, 'last_col')) %>% 
    dplyr::select(-last_col) %>% 
    # convert time to float
    tidyr::extract(time, c('hours', 'min', 'sec'), '(\\d+):(\\d+):(\\d+)') %>%
    dplyr::mutate(time=as.numeric(hours)*3600 + as.numeric(min)*60 + as.numeric(sec)) %>%
    dplyr::select(-hours, -min, -sec) %>% 
    # reshape wide to long
    tidyr::gather(col, value, dplyr::matches("\\d+")) %>% 
    dplyr::mutate(well=paste0(row, col), col=as.numeric(col), channel=.ch) %>% 
    filter(!is.na(value))
  return(as.tibble(.data))
}

read_plate_layout <- function(.path) {
  #Extract the folder name where the spec data are stored: this is the expId. Must be of the sort: date_description. Ex: 20191115_test
  #  read_spec_kinetic() is written such as to call read.table only once 
  # (much faster than calling it for each timepoint and channel)
  .lines <- readLines(.path)
  #print(.lines)
  .l_idx <- stringr::str_detect(.lines, "Type:") %>% which %>% (function(.x) .x)
  
  noFirst <- function(.l){
    new_l <- str_split(.l,",")[[1]]
    new_l <- new_l[2:length(new_l)]
    return(new_l)}
  
  return_plate <- function(.index){
    .col_title <- str_match(.lines[[.index]],"Type: ([a-zA-Z0-9]{1,}),,,,,,,,,,,,$")[[2]]
    .data <- .lines[c((.index+2):(.index+9))]
    .values <- lapply(.data,noFirst) %>% unlist()
    .rows <- rep(LETTERS[c(1:8)],each=12)
    .cols <- rep(c(1:12),times=8)
    new_df <- data_frame(row=.rows,column=.cols,values=.values,type=rep(c(.col_title),each=96))
    return(new_df)}
  
  new_df <- .l_idx %>% lapply(return_plate) %>% 
    bind_rows() %>% 
    pivot_wider(id_cols=c(row,column),names_from=type,values_from=values)
  
  return(new_df)}



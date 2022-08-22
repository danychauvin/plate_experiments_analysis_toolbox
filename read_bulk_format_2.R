# Import necessary packages and functions to read plate reader data
# Dany Chauvin (dany.chauvin@unibas.ch) 20220409

myconditions <- readr::read_csv(path_to_data_list,
                                col_types = cols(
                                  date=col_character(),
                                  description=col_character(),
                                  data_path=col_character(),
                                  plate_map_path=col_character()))

# Import the data
mydata <- myconditions %>%
  group_by(date,description,plate) %>% #date here is as per datalist, this is not the date of the measurement, description and plate should refer, together, to a unique experiment
    do((function(.df){
    plate_ <- as.character(unique(.df$plate))
    date_ <- as.character(unique(.df$date))
    description_ <- as.character(unique(.df$description))
    path_ <- as.character(unique(.df$data_path))
    path_<-"./data/220818_20220816_acetate0.2.txt"
    new_df <- lapply(channels, function(c.){read_Biotek_Synergy_stacker(path_,c.)}) %>% 
      bind_rows() %>% 
      mutate(date=date_,
             description=description_,
             plate=plate_)
    return(new_df)})(.)) %>% 
  select(date,description,plate,row,col,value,channel,time,step) %>% 
  ungroup() %>% 
  pivot_wider(id_cols=c("date","description","plate","col","row","step"),names_from=c("channel"),values_from=c("value","time")) %>% 
  rename(od="value_OD600:600",
         time="time_OD600:600",
         column=col) %>% 
  mutate(time_min=time/60)

# Import plate layout
mylayouts <- myconditions %>%
  group_by(date,description,plate) %>% 
  do((function(.df){
    plate_ <- as.character(unique(.df$plate))
    date_ <- as.character(unique(.df$date))
    description_ <- as.character(unique(.df$description))
    path_ <- as.character(unique(.df$plate_map_path))
    new_df <- read_plate_layout(path_) %>% 
      bind_rows() %>% 
      mutate(date=date_,
             description=description_,
             plate=plate_)
    return(new_df)})(.)) %>% 
  ungroup()
  #rename(exp_start_date=date)

# Completed data

mydata <- mydata %>% 
  left_join(mylayouts,by=c("date","plate","description","row","column")) %>% 
  mutate(well=paste(row,column,sep=""))


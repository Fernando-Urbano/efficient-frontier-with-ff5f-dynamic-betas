
library(tidyverse)
library(lubridate)
library(readr)
library(readxl)
library(rmgarch)
library(quantmod)
library(fGarch)
library(cowplot)
library(tseries)

####################################### T-TEST #######################################

#### Returns ####
osr_6m = tibble(Model = as.character(), `Average Return` = as.numeric(), `In-Sample` = as.character())
osr_12m = tibble(Model = as.character(), `Average Return` = as.numeric(), `In-Sample` = as.character())
osr_24m = tibble(Model = as.character(), `Average Return` = as.numeric(), `In-Sample` = as.character())

for (x in c("1999-2000", "2001-2002", "2003-2004",
            "2005-2006", "2009-2010", "2011-2012",
            "2013-2014", "2015-2016", "2017-2018")){
  
  setwd(paste0("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/",x,"/Efficient Frontier"))
  
  osr_6m = osr_6m %>%
    rbind(read_excel(paste0('Six Month Average Return Difference Statistics (',x,').xlsx')) %>% 
            mutate(`In-Sample` = x))
  
  osr_12m = osr_12m %>%
    rbind(read_excel(paste0("One Year Average Return Difference Statistics (",x,").xlsx")) %>% 
                       mutate(`In-Sample` = x))
  
  osr_24m = osr_24m %>%
    rbind(read_excel(paste0("Two Years Average Return Difference Statistics (",x,").xlsx")) %>% 
                       mutate(`In-Sample` = x))
}


statistics = function (x){
  return(osr_6m_statistics = x %>% 
           purrr::set_names("model", "value", "time") %>% 
           mutate(value = case_when(value < 0.05 ~ 1,
                                    value >= 0.05 ~ 0)) %>% 
           group_by(model) %>% 
           summarise(value = sum(value)))
}

osr_6m_statistics = statistics(osr_6m) %>% 
  mutate(model = model %>% str_remove(" P-Value")) %>% 
  rename("6M Returns" = value)

osr_12m_statistics = statistics(osr_12m) %>% 
  mutate(model = model %>% str_remove(" P-Value")) %>% 
  rename("1Y Returns" = value)

osr_24m_statistics = statistics(osr_24m) %>% 
  mutate(model = model %>% str_remove(" P-Value")) %>% 
  rename("2Y Returns" = value)

#### Sharpe ####
oss_6m = tibble(Model = as.character(), Sharpe = as.numeric(), `In-Sample` = as.character())
oss_12m = tibble(Model = as.character(), Sharpe = as.numeric(), `In-Sample` = as.character())
oss_24m = tibble(Model = as.character(), Sharpe = as.numeric(), `In-Sample` = as.character())

for (x in c("1999-2000", "2001-2002", "2003-2004",
            "2005-2006", "2009-2010", "2011-2012",
            "2013-2014", "2015-2016", "2017-2018")){
  
  setwd(paste0("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/",x,"/Efficient Frontier"))
  
  oss_6m = oss_6m %>%
    rbind(read_excel(paste0('Six Month Sharpe Difference Statistics (',x,').xlsx')) %>% 
            mutate(`In-Sample` = x))
  
  oss_12m = oss_12m %>%
    rbind(read_excel(paste0("One Year Sharpe Difference Statistics (",x,").xlsx")) %>% 
            mutate(`In-Sample` = x))
  
  oss_24m = oss_24m %>%
    rbind(read_excel(paste0("Two Years Sharpe Difference Statistics (",x,").xlsx")) %>% 
            mutate(`In-Sample` = x))
}


statistics = function (x){
  return(oss_6m_statistics = x %>% 
           purrr::set_names("model", "value", "time") %>% 
           mutate(value = case_when(value < 0.05 ~ 1,
                                    value >= 0.05 ~ 0)) %>% 
           group_by(model) %>% 
           summarise(value = sum(value)))
}

oss_6m_statistics = statistics(oss_6m) %>% 
  mutate(model = model %>% str_remove(" P-Value")) %>% 
  rename("6M Sharpe" = value)

oss_12m_statistics = statistics(oss_12m) %>% 
  mutate(model = model %>% str_remove(" P-Value")) %>% 
  rename("1Y Sharpe" = value)

oss_24m_statistics = statistics(oss_24m) %>% 
  mutate(model = model %>% str_remove(" P-Value")) %>% 
  rename("2Y Sharpe" = value)

#### Statistics ####
statistics = plyr::join_all(list(
  osr_6m_statistics, osr_12m_statistics, osr_24m_statistics,
  oss_6m_statistics, oss_12m_statistics, oss_24m_statistics),
  by = 'model', type = 'full') %>% 
  select(model, contains("6M"), contains("1Y"), contains("2Y"))

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores")
statistics %>% writexl::write_xlsx("Difference Statistics.xlsx")

####################################### AVERAGE #######################################


ost_6m = tibble(Model = as.character(),
                `Annualized Return` = as.numeric(),
                `Annualized Vol` = as.numeric(),
                `Annualized Sharpe` = as.numeric(),
                `In-Sample` = as.character())

ost_12m = ost_6m
ost_24m = ost_6m


for (x in c("1999-2000", "2001-2002", "2003-2004",
            "2005-2006", "2009-2010", "2011-2012",
            "2013-2014", "2015-2016", "2017-2018")){
  
  setwd(paste0("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/",x,"/Efficient Frontier"))
  
  ost_6m = ost_6m %>%
    rbind(read_excel(paste0('Six Month Out-of-Sample (',x,').xlsx')) %>% 
            mutate(`In-Sample` = x))
  
  ost_12m = ost_12m %>%
    rbind(read_excel(paste0('One Year Out-of-Sample (',x,').xlsx')) %>% 
            mutate(`In-Sample` = x))
  
  ost_24m = ost_24m %>%
    rbind(read_excel(paste0('Two Years Out-of-Sample (',x,').xlsx')) %>% 
            mutate(`In-Sample` = x))
}

statistics = function (x){
  x %>% 
    group_by(Model) %>% 
    summarise(`Annualized Return` = mean(`Annualized Return`),
              `Annualized Vol` = mean(`Annualized Vol`),
              `Annualized Sharpe` = mean(`Annualized Sharpe`)) %>% 
    arrange(-`Annualized Return`) %>% 
    mutate_at(c('Annualized Return', 'Annualized Vol'), ~ scales::percent(
      ., accuracy = 0.1, big.mark = ',', decimal.mark = '.'))
}  

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores")
statistics(ost_6m) %>% 
  rename("6M Out-of-Sample" = Model) %>% 
  writexl::write_xlsx("6M Statistics.xlsx")
statistics(ost_12m) %>% 
  rename("1Y Out-of-Sample" = Model) %>% 
  writexl::write_xlsx("1Y Statistics.xlsx")
statistics(ost_24m) %>% 
  rename("2Y Out-of-Sample" = Model) %>% 
  writexl::write_xlsx("2Y Statistics.xlsx")






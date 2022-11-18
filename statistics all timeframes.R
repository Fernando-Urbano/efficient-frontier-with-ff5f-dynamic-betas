library(tidyverse)
library(lubridate)
library(readr)
library(readxl)
library(quantmod)
library(fGarch)
library(cowplot)
library(tseries)

df = c()

for (x in c("1999-2000", "2001-2002", "2003-2004",
            "2005-2006", "2009-2010", "2011-2012",
            "2013-2014", "2015-2016", "2017-2018")){
  
  setwd(paste0("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/",x))

# GARCH
  df = df %>% rbind(read_excel(paste0("Statistics - GARCH Dist GED (",x,").xlsx")) %>% 
                      mutate(model = 'GARCH GED') %>% 
                      mutate(timeframe = x))
  df = df %>% rbind(read_excel(paste0("Statistics - GARCH Dist Skew GED (",x,").xlsx")) %>% 
                      mutate(model = 'GARCH Skew GED') %>% 
                      mutate(timeframe = x))
  df = df %>% rbind(read_excel(paste0("Statistics - GARCH Dist Gaussian (",x,").xlsx")) %>% 
                      mutate(model = 'GARCH Gaussian') %>% 
                      mutate(timeframe = x))
  df = df %>% rbind(read_excel(paste0("Statistics - GARCH Dist Skew Gaussian (",x,").xlsx")) %>% 
                      mutate(model = 'GARCH Skew Gaussian') %>% 
                      mutate(timeframe = x))
  
  # APARCH
  df = df %>% rbind(read_excel(paste0("Statistics - APARCH Dist GED (",x,").xlsx")) %>% 
                      mutate(model = 'APARCH GED') %>% 
                      mutate(timeframe = x))
  df = df %>% rbind(read_excel(paste0("Statistics - APARCH Dist Skew GED (",x,").xlsx")) %>% 
                      mutate(model = 'APARCH Skew GED') %>% 
                      mutate(timeframe = x))
  df = df %>% rbind(read_excel(paste0("Statistics - APARCH Dist Gaussian (",x,").xlsx")) %>% 
                      mutate(model = 'APARCH Gaussian') %>% 
                      mutate(timeframe = x))
  df = df %>% rbind(read_excel(paste0("Statistics - APARCH Dist Skew Gaussian (",x,").xlsx")) %>% 
                      mutate(model = 'APARCH Skew Gaussian') %>% 
                      mutate(timeframe = x))

}

models_names = c("GARCH Gaussian", "GARCH Skew Gaussian", "GARCH GED", "GARCH Skew GED",
                 "APARCH Gaussian", "APARCH Skew Gaussian", "APARCH GED", "APARCH Skew GED",
                 "Actual Returns")

df %>% 
  dplyr::filter(id == "MSE") %>% 
  group_by(model) %>% 
  summarise(value = mean(value)) %>% 
  arrange(value) %>% 
  mutate(rank = '#' %>% paste0(seq(1,8))) %>% 
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(value, model)) +
  geom_bar(fill='black', stat = 'identity') +
  coord_cartesian(xlim = c(5.8e-5,6.1e-5)) +
  labs(title = 'Average MSE by Model',
       subtitle = 'All timeframes considered',
       x = '',
       y = '') +
  geom_text(aes(x=value-3.5e-7, y=model, label=scales::scientific(value, digits = 3)),
            color='white') +
  geom_text(aes(x=value+1.5e-7, y=model, label=rank),
            color='black') +
  theme_minimal()
setwd("C:/Users/Dinho Urbano/Desktop/Script R/DCC Distributions/timeframes")
ggsave("Average MSE.png", device = 'png', units = 'in',
       width = 7, height = 5)
  

df %>% 
  dplyr::filter(id == "MAE") %>% 
  group_by(model) %>% 
  summarise(value = mean(value)) %>% 
  arrange(value) %>% 
  mutate(rank = '#' %>% paste0(seq(1,8))) %>% 
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(value, model)) +
  geom_bar(fill='black', stat = 'identity') +
  coord_cartesian(xlim = c(0.00515,0.00528)) +
  labs(title = 'Average MAE by Model',
       subtitle = 'All timeframes considered',
       x = '',
       y = '') +
  geom_text(aes(x=value-0.000015, y=model, label=scales::scientific(value, digits = 3)),
            color='white') +
  geom_text(aes(x=value+0.000008, y=model, label=rank),
            color='black') +
  theme_minimal()
setwd("C:/Users/Dinho Urbano/Desktop/Script R/DCC Distributions/timeframes")
ggsave("Average MAE.png", device = 'png', units = 'in',
       width = 7, height = 5)

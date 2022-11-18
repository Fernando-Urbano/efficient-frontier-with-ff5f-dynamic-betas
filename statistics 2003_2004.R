library(tidyverse)
library(lubridate)
library(readr)
library(readxl)
library(quantmod)
library(fGarch)
library(cowplot)
library(tseries)

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/2003-2004")
# GARCH
garch_ged = read_excel("Statistics - GARCH Dist GED (2003-2004).xlsx") %>% 
  mutate(model = 'GARCH GED')
garch_sged = read_excel("Statistics - GARCH Dist Skew GED (2003-2004).xlsx") %>% 
  mutate(model = 'GARCH Skew GED')
garch_norm = read_excel("Statistics - GARCH Dist Gaussian (2003-2004).xlsx") %>% 
  mutate(model = 'GARCH Gaussian')
garch_snorm = read_excel("Statistics - GARCH Dist Skew Gaussian (2003-2004).xlsx") %>% 
  mutate(model = 'GARCH Skew Gaussian')

# APARCH
aparch_ged = read_excel("Statistics - APARCH Dist GED (2003-2004).xlsx") %>% 
  mutate(model = 'APARCH GED')
aparch_sged = read_excel("Statistics - APARCH Dist Skew GED (2003-2004).xlsx") %>% 
  mutate(model = 'APARCH Skew GED')
aparch_norm = read_excel("Statistics - APARCH Dist Gaussian (2003-2004).xlsx") %>% 
  mutate(model = 'APARCH Gaussian')
aparch_snorm = read_excel("Statistics - APARCH Dist Skew Gaussian (2003-2004).xlsx") %>% 
  mutate(model = 'APARCH Skew Gaussian')

# JUNÇÃO
statistics = rbind(garch_norm, garch_snorm,
                   garch_ged, garch_sged,
                   aparch_norm, aparch_snorm,
                   aparch_ged, aparch_sged) %>% 
  mutate(id = case_when(
    id == 'Shapiro Wild - P-Valor' ~ 'Shapiro Wild - P-Value',
    id == 'Jarque Bera - P-Valor' ~ 'Jarque Bera - P-Value',
    id == 'Ljung Box 90 Erro - P-Valor' ~ 'Ljung Box 90 Residuals - P-Value',
    id == 'Ljung Box 30 Erro - P-Valor' ~ 'Ljung Box 30 Residuals - P-Value',
    id == 'Ljung Box 90 Erro ao Quadrado - P-Valor' ~ 'Ljung Box 90 Squared Residuals - P-Value',
    id == 'Ljung Box 30 Erro ao Quadrado - P-Valor' ~ 'Ljung Box 30 Squared Residuals - P-Value',
    id == 'Assimetria' ~ 'Skewness',
    id == 'Curtose' ~ 'Kurtosis',
    TRUE ~ id))

# INDUSTRY STATISTICS
industry_names = statistics %>% distinct(name) %>% unlist()


for (i in 1:length(industry_names)){
  industry_statistics_raw = statistics %>% 
    dplyr::filter(name == industry_names[i]) %>%
    purrr::set_names("name", industry_names[i], "value", "Model") %>% 
    spread(industry_names[i], value) %>% 
    select(-name)
  
  industry_statistics = industry_statistics_raw %>% 
    purrr::set_names(industry_names[i],
                     names(industry_statistics_raw[2:ncol(industry_statistics_raw)]))
  
  setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/2003-2004/Model Statistics")
  writexl::write_xlsx(industry_statistics, path = paste(industry_names[i], 'Statistics (2003-2004).xlsx'))
  
}

models_names = c("GARCH Gaussian", "GARCH Skew Gaussian",
                 "GARCH GED", "GARCH Skew GED",
                 "APARCH Gaussian", "APARCH Skew Gaussian",
                 "APARCH GED", "APARCH Skew GED")

# MSE
average = statistics %>% 
  dplyr::filter(id == 'MSE') %>% 
  group_by(model) %>% 
  summarise(value = mean(value)) %>% 
  mutate(name = 'Average') %>% 
  mutate(id = 'MSE') %>% 
  dplyr::select(name, id, value, model)

statistics %>% 
  dplyr::filter(id == 'MSE') %>% 
  rbind(average) %>% 
  mutate(color_id = case_when(
    name %in% c("Chems", "Cnstr",
                "FabPr", "Food", "Mines",
                "Oil", "Rtail", "Trans") ~ "Red",
    name %in% c("Average") ~ "Black",
    TRUE ~ "Blue")) %>%  
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(value, model)) +
  geom_point(aes(color=color_id), size=4) +
  geom_line(aes(group=name, color=color_id)) +
  scale_color_manual(values = c('darkred', 'gray', 'black')) +
  facet_wrap(~ name, ncol=6, scales = "free_x") +
  theme_minimal() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  scale_x_continuous(labels = scales::scientific_format(digits = 3)) +
  labs(title = 'MSE (2003-2004)',
       x = '',
       y = 'Model')

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/2003-2004/Model Statistics")
ggsave("MSE (2003-2004).png", device = 'png', units = 'in', width = 17, height = 13)

# MAE
average = statistics %>% 
  dplyr::filter(id == 'MAE') %>% 
  group_by(model) %>% 
  summarise(value = mean(value)) %>% 
  mutate(name = 'Average') %>% 
  mutate(id = 'MAE') %>% 
  dplyr::select(name, id, value, model)

statistics %>% 
  dplyr::filter(id == 'MAE') %>% 
  rbind(average) %>% 
  mutate(color_id = case_when(
    name %in% c("Chems", "Cnstr",
                "FabPr", "Food", "Mines",
                "Oil", "Rtail", "Trans") ~ "Red",
    name %in% c("Average") ~ "Black",
    TRUE ~ "Blue")) %>% 
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(value, model)) +
  geom_point(aes(color=color_id), size=4) +
  geom_line(aes(group=name, color=color_id)) +
  scale_color_manual(values = c('darkred', 'gray', 'black')) +
  facet_wrap(~ name, ncol=6, scales = "free_x") +
  theme_minimal() +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  scale_x_continuous(labels = scales::scientific_format(digits = 3)) +
  labs(title = 'MAE - Dynamic Betas Adjusted Excess of Return (2003-2004)',
       x = '',
       y = 'Model')

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/2003-2004/Model Statistics")
ggsave("MAE (2003-2004).png", device = 'png', units = 'in', width = 17, height = 13)

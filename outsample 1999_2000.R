# OUT-OF-SAMPLE - 1999-2000 IN-SAMPLE
library(tidyverse)
library(lubridate)
library(readr)
library(readxl)
library(rmgarch)
library(quantmod)
library(fGarch)
library(cowplot)
library(tseries)
library(gtools)


setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores")

# FATORES
setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores")
factors_df <- read_csv("Fama French 5 Factors.csv") %>% 
  mutate(date = date %>% ymd()) %>% 
  mutate_if(is.numeric, ~ . / 100) %>% 
  rename(MKT = `Mkt-RF`)

# INDÚSTRIAS INSAMPLE
insample_returns <- read_csv("C:/Users/Dinho Urbano/Desktop/TCC/Dados Financeiros/17 Industries.csv") %>% 
  mutate(date = date %>% ymd()) %>%
  mutate_if(is.numeric, ~ case_when(. == -99.99 ~ NA %>% as.numeric(),
                                    TRUE ~ .)) %>% 
  mutate_if(is.numeric, ~ . / 100) %>% 
  dplyr::filter(date >= "1999-01-01" & date <= "2000-12-31") %>% 
  merge(factors_df %>% dplyr::select(date, RF), by='date') %>% 
  mutate_if(is.numeric, ~ (1 + .) / (1 + RF) - 1) %>% 
  dplyr::select(-RF) %>% 
  gather(id, value, -date) %>% 
  spread(id, value) %>% 
  dplyr::select(-Cnstr)

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000")

# GARCH
garch_ged = read_excel("Returns Hat - GARCH Dist GED (1999-2000).xlsx") %>% 
  mutate(model = 'GARCH GED') %>% 
  dplyr::filter(id != "Cnstr")
garch_sged = read_excel("Returns Hat - GARCH Dist Skew GED (1999-2000).xlsx") %>% 
  mutate(model = 'GARCH Skew GED') %>% 
  dplyr::filter(id != "Cnstr")
garch_norm = read_excel("Returns Hat - GARCH Dist Gaussian (1999-2000).xlsx") %>% 
  mutate(model = 'GARCH Gaussian') %>% 
  dplyr::filter(id != "Cnstr")
garch_snorm = read_excel("Returns Hat - GARCH Dist Skew Gaussian (1999-2000).xlsx") %>% 
  mutate(model = 'GARCH Skew Gaussian') %>% 
  dplyr::filter(id != "Cnstr")

# APARCH
aparch_ged = read_excel("Returns Hat - APARCH Dist GED (1999-2000).xlsx") %>% 
  mutate(model = 'APARCH GED') %>% 
  dplyr::filter(id != "Cnstr")
aparch_sged = read_excel("Returns Hat - APARCH Dist Skew GED (1999-2000).xlsx") %>% 
  mutate(model = 'APARCH Skew GED') %>% 
  dplyr::filter(id != "Cnstr")
aparch_norm = read_excel("Returns Hat - APARCH Dist Gaussian (1999-2000).xlsx") %>% 
  mutate(model = 'APARCH Gaussian') %>% 
  dplyr::filter(id != "Cnstr")
aparch_snorm = read_excel("Returns Hat - APARCH Dist Skew Gaussian (1999-2000).xlsx") %>% 
  mutate(model = 'APARCH Skew Gaussian') %>% 
  dplyr::filter(id != "Cnstr")

# MODELO DOS RETORNOS OCORRIDOS
actual_returns = insample_returns %>% 
  gather(id, stock_rt_hat, -date) %>% 
  mutate(model = "Actual Returns") %>% 
  select(date, stock_rt_hat, id, model)

sample_numbers = combinations(n = 16, r = 10, repeats.allowed = FALSE, v = seq(1, 17, 1)) %>% 
  data.frame() %>% 
  purrr::set_names("stock_1", "stock_2", "stock_3", "stock_4", "stock_5",
                   "stock_6", "stock_7", "stock_8", "stock_9", "stock_10")

# FUNÇÃO PARA DEFINIR PESOS E CALCULAR SHARPE OUT-OF-SAMPLE
outsample_statistics = function(dist_df, stock_n){
  
  stock_n = stock_n %>% unlist()
  
  insample_matrix.list = dist_df %>% 
    spread(id, stock_rt_hat) %>% 
    mutate(date = date %>% ymd()) %>% 
    select(-model) %>% 
    column_to_rownames("date")
  
  insample_matrix = insample_matrix.list %>% 
    select(stock_n) %>% 
    purrr::set_names(
      colnames(insample_matrix.list)[stock_n[1]],
      colnames(insample_matrix.list)[stock_n[2]],
      colnames(insample_matrix.list)[stock_n[3]],
      colnames(insample_matrix.list)[stock_n[4]],
      colnames(insample_matrix.list)[stock_n[5]],
      colnames(insample_matrix.list)[stock_n[6]],
      colnames(insample_matrix.list)[stock_n[7]],
      colnames(insample_matrix.list)[stock_n[8]],
      colnames(insample_matrix.list)[stock_n[9]],
      colnames(insample_matrix.list)[stock_n[10]]
    )
  
  weights.list = as.matrix(inv(cov(insample_matrix))) %*% as.matrix(rep(1, 10))
  weights = weights.list / sum(weights.list)
  
  outsample_matrix = (returns %>% 
                        column_to_rownames('date') %>% 
                        select(insample_matrix %>% colnames())) %>% as.matrix()
  
  outsample_portfolio = outsample_matrix %*% weights
  rt_outsample_portfolio = ((outsample_portfolio %>%
                               data.frame() %>% rownames_to_column('date') %>%
                               mutate(date = date %>% ymd()) %>% 
                               arrange(date) %>% 
                               purrr::set_names('date', 'value') %>% 
                               mutate(value = cumprod(1 + value) - 1) %>% 
                               arrange(date %>% desc()))[1,'value'] + 1) ^ (1 / length(outsample_portfolio)) - 1
  vol_outsample_portfolio = sd(outsample_portfolio)
  sharpe_outsample_portfolio = rt_outsample_portfolio / vol_outsample_portfolio
  
  return(data.frame(
    rt = rt_outsample_portfolio,
    vol = vol_outsample_portfolio,
    sharpe = sharpe_outsample_portfolio) %>% 
      mutate(
        model = dist_df$model[1],
        stock_1 =  (insample_matrix %>% colnames)[1],
        stock_2 =  (insample_matrix %>% colnames)[2],
        stock_3 =  (insample_matrix %>% colnames)[3],
        stock_4 =  (insample_matrix %>% colnames)[4],
        stock_5 =  (insample_matrix %>% colnames)[5],
        stock_6 =  (insample_matrix %>% colnames)[6],
        stock_7 =  (insample_matrix %>% colnames)[7],
        stock_8 =  (insample_matrix %>% colnames)[8],
        stock_9 =  (insample_matrix %>% colnames)[9],
        stock_10 =  (insample_matrix %>% colnames)[10]
      )
  )
  
}

##### INDÚSTRIAS - 6 MESES #####
returns <- read_csv("C:/Users/Dinho Urbano/Desktop/TCC/Dados Financeiros/17 Industries.csv") %>% 
  mutate(date = date %>% ymd()) %>%
  mutate_if(is.numeric, ~ case_when(. == -99.99 ~ NA %>% as.numeric(),
                                    TRUE ~ .)) %>% 
  mutate_if(is.numeric, ~ . / 100) %>% 
  dplyr::filter(date >= "2001-01-01" & date <= "2001-06-30") %>% 
  merge(factors_df %>% dplyr::select(date, RF), by='date') %>% 
  mutate_if(is.numeric, ~ (1 + .) / (1 + RF) - 1) %>% 
  dplyr::select(-RF) %>% 
  gather(id, value, -date) %>% 
  spread(id, value) %>% 
  dplyr::select(-Cnstr)


# UTILIZAÇÃO DA FUNÇÃO NOS PORTFÓLIOS

for (i in 1:nrow(sample_numbers)){
  
  new_outsample_statistics = rbind(
    outsample_statistics(garch_norm, unlist(sample_numbers[i,])),
    outsample_statistics(garch_snorm, unlist(sample_numbers[i,])),
    outsample_statistics(garch_ged, unlist(sample_numbers[i,])),
    outsample_statistics(garch_sged, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_norm, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_snorm, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_ged, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_sged, unlist(sample_numbers[i,])),
    outsample_statistics(actual_returns, unlist(sample_numbers[i,]))
  )
  
  if (i == 1){
    all_outsample_statistics = new_outsample_statistics
  } else {
    all_outsample_statistics = all_outsample_statistics %>% rbind(new_outsample_statistics)
  }
  
  setTxtProgressBar(txtProgressBar(
    min = 1, max = nrow(sample_numbers), initial = 0, style = 3), i)
}

models_names = c("GARCH Gaussian", "GARCH Skew Gaussian", "GARCH GED", "GARCH Skew GED",
                 "APARCH Gaussian", "APARCH Skew Gaussian", "APARCH GED", "APARCH Skew GED",
                 "Actual Returns")

all_outsample_statistics %>% 
  select(sharpe, model) %>%
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(x=sharpe*(252^(1/2)), y=model, fill=sharpe)) +
  geom_boxplot(color = 'black') +
  # geom_jitter(color="gray", size=0.4, alpha=0.3) +
  theme_minimal() +
  labs(title = 'Sharpe Out-of-Sample',
       subtitle = 'In-Sample: 1999-01-01 to 2000-12-31 \nOut-of-Sample: 2001-01-01 to 2001-06-30',
       x = 'Annualized Sharpe',
       y = '') +
  theme(legend.position = 'none') +
  scale_x_continuous(labels = scales::number_format(big.mark = ',',
                                                    decimal.mark = '.'))
setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
ggsave("Six Month Out-of-Sample Sharpe (1999-2000).png", device = 'png', units = 'in',
       width = 10, height = 7)

all_outsample_statistics %>% 
  select(rt, model) %>%
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(x=rt*(252), y=model, fill=rt)) +
  geom_boxplot(color = 'black') +
  # geom_jitter(color="gray", size=0.4, alpha=0.3) +
  theme_minimal() +
  labs(title = 'Returns Out-of-Sample',
       subtitle = 'In-Sample: 1999-01-01 to 2000-12-31 \nOut-of-Sample: 2001-01-01 to 2001-06-30',
       x = 'Annualized Returns',
       y = '') +
  theme(legend.position = 'none') +
  scale_x_continuous(labels = scales::number_format(big.mark = ',',
                                                    decimal.mark = '.'))
setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
ggsave("Six Month Out-of-Sample Returns (1999-2000).png", device = 'png', units = 'in',
       width = 10, height = 7)

all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(rt = mean(rt)*(252)) %>% 
  mutate(vol = mean(vol)*(252^(1/2))) %>% 
  mutate(sharpe = mean(sharpe)*(252^(1/2))) %>% 
  distinct(model, .keep_all = TRUE) %>% 
  select(model, rt, vol, sharpe) %>% 
  purrr::set_names("Model", "Annualized Return",
                   "Annualized Vol", "Annualized Sharpe") -> model_six_month


setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
model_six_month %>% writexl::write_xlsx("Six Month Out-of-Sample (1999-2000).xlsx")

# SHARPE DIFFERENCES
all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(n = seq(1, nrow(all_outsample_statistics)/9)) %>% 
  dplyr::select(n, sharpe, model) %>% 
  spread(model, sharpe) -> sharpe_outsample

(pvalue_sharpe_statistics = tibble(
  "GARCH Gaussian P-Value" = (t.test(x = sharpe_outsample$`GARCH Gaussian`,
                                     y = sharpe_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "GARCH Skew Gaussian P-Value" = (t.test(x = sharpe_outsample$`GARCH Skew Gaussian`,
                                          y = sharpe_outsample$`Actual Returns`,
                                          alternative = "greater",
                                          conf.level = 0.95))$p.value,
  "GARCH GED P-Value" = (t.test(x = sharpe_outsample$`GARCH GED`,
                                y = sharpe_outsample$`Actual Returns`,
                                alternative = "greater",
                                conf.level = 0.95))$p.value,
  "GARCH Skew GED P-Value" = (t.test(x = sharpe_outsample$`GARCH Skew GED`,
                                     y = sharpe_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "APARCH Gaussian P-Value" = (t.test(x = sharpe_outsample$`APARCH Gaussian`,
                                      y = sharpe_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
  "APARCH Skew Gaussian P-Value" = (t.test(x = sharpe_outsample$`APARCH Skew Gaussian`,
                                           y = sharpe_outsample$`Actual Returns`,
                                           alternative = "greater",
                                           conf.level = 0.95))$p.value,
  "APARCH GED P-Value" = (t.test(x = sharpe_outsample$`APARCH GED`,
                                 y = sharpe_outsample$`Actual Returns`,
                                 alternative = "greater",
                                 conf.level = 0.95))$p.value,
  "APARCH Skew GED P-Value" = (t.test(x = sharpe_outsample$`APARCH Skew GED`,
                                      y = sharpe_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
) %>% 
    gather(Model, Sharpe))

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
pvalue_sharpe_statistics %>%
  writexl::write_xlsx("Six Month Sharpe Difference Statistics (1999-2000).xlsx")

# MEAN RETURN DIFFERENCES
all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(n = seq(1, nrow(all_outsample_statistics)/9)) %>% 
  dplyr::select(n, rt, model) %>% 
  spread(model, rt) -> rt_outsample

(pvalue_rt_statistics = tibble(
  "GARCH Gaussian P-Value" = (t.test(x = rt_outsample$`GARCH Gaussian`,
                                     y = rt_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "GARCH Skew Gaussian P-Value" = (t.test(x = rt_outsample$`GARCH Skew Gaussian`,
                                          y = rt_outsample$`Actual Returns`,
                                          alternative = "greater",
                                          conf.level = 0.95))$p.value,
  "GARCH GED P-Value" = (t.test(x = rt_outsample$`GARCH GED`,
                                y = rt_outsample$`Actual Returns`,
                                alternative = "greater",
                                conf.level = 0.95))$p.value,
  "GARCH Skew GED P-Value" = (t.test(x = rt_outsample$`GARCH Skew GED`,
                                     y = rt_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "APARCH Gaussian P-Value" = (t.test(x = rt_outsample$`APARCH Gaussian`,
                                      y = rt_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
  "APARCH Skew Gaussian P-Value" = (t.test(x = rt_outsample$`APARCH Skew Gaussian`,
                                           y = rt_outsample$`Actual Returns`,
                                           alternative = "greater",
                                           conf.level = 0.95))$p.value,
  "APARCH GED P-Value" = (t.test(x = rt_outsample$`APARCH GED`,
                                 y = rt_outsample$`Actual Returns`,
                                 alternative = "greater",
                                 conf.level = 0.95))$p.value,
  "APARCH Skew GED P-Value" = (t.test(x = rt_outsample$`APARCH Skew GED`,
                                      y = rt_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
) %>% 
    gather(Model, `Average Return`))

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
pvalue_rt_statistics %>%
  writexl::write_xlsx("Six Month Average Return Difference Statistics (1999-2000).xlsx")

##### INDÚSTRIAS - 12 MESES #####
returns <- read_csv("C:/Users/Dinho Urbano/Desktop/TCC/Dados Financeiros/17 Industries.csv") %>% 
  mutate(date = date %>% ymd()) %>%
  mutate_if(is.numeric, ~ case_when(. == -99.99 ~ NA %>% as.numeric(),
                                    TRUE ~ .)) %>% 
  mutate_if(is.numeric, ~ . / 100) %>% 
  dplyr::filter(date >= "2001-01-01" & date <= "2001-12-31") %>% 
  merge(factors_df %>% dplyr::select(date, RF), by='date') %>% 
  mutate_if(is.numeric, ~ (1 + .) / (1 + RF) - 1) %>% 
  dplyr::select(-RF) %>% 
  gather(id, value, -date) %>% 
  spread(id, value) %>% 
  dplyr::select(-Cnstr)


# UTILIZAÇÃO DA FUNÇÃO NOS PORTFÓLIOS

for (i in 1:nrow(sample_numbers)){
  
  new_outsample_statistics = rbind(
    outsample_statistics(garch_norm, unlist(sample_numbers[i,])),
    outsample_statistics(garch_snorm, unlist(sample_numbers[i,])),
    outsample_statistics(garch_ged, unlist(sample_numbers[i,])),
    outsample_statistics(garch_sged, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_norm, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_snorm, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_ged, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_sged, unlist(sample_numbers[i,])),
    outsample_statistics(actual_returns, unlist(sample_numbers[i,]))
  )
  
  if (i == 1){
    all_outsample_statistics = new_outsample_statistics
  } else {
    all_outsample_statistics = all_outsample_statistics %>% rbind(new_outsample_statistics)
  }
  
  setTxtProgressBar(txtProgressBar(
    min = 1, max = nrow(sample_numbers), initial = 0, style = 3), i)
}

models_names = c("GARCH Gaussian", "GARCH Skew Gaussian", "GARCH GED", "GARCH Skew GED",
                 "APARCH Gaussian", "APARCH Skew Gaussian", "APARCH GED", "APARCH Skew GED",
                 "Actual Returns")

all_outsample_statistics %>% 
  select(sharpe, model) %>%
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(x=sharpe*(252^(1/2)), y=model, fill=sharpe)) +
  geom_boxplot(color = 'black') +
  # geom_jitter(color="gray", size=0.4, alpha=0.05) +
  theme_minimal() +
  labs(title = 'Sharpe Out-of-Sample',
       subtitle = 'In-Sample: 1999-01-01 to 2000-12-31 \nOut-of-Sample: 2001-01-01 to 2001-12-31',
       x = 'Annualized Sharpe',
       y = '') +
  theme(legend.position = 'none') +
  scale_x_continuous(labels = scales::number_format(big.mark = ',',
                                                    decimal.mark = '.'))
setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
ggsave("One Year Out-of-Sample Sharpe (1999-2000).png", device = 'png', units = 'in',
       width = 10, height = 7)

all_outsample_statistics %>% 
  select(rt, model) %>%
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(x=rt*(252), y=model, fill=rt)) +
  geom_boxplot(color = 'black') +
  # geom_jitter(color="gray", size=0.4, alpha=0.3) +
  theme_minimal() +
  labs(title = 'Returns Out-of-Sample',
       subtitle = 'In-Sample: 1999-01-01 to 2000-12-31 \nOut-of-Sample: 2001-01-01 to 2001-12-31',
       x = 'Annualized Returns',
       y = '') +
  theme(legend.position = 'none') +
  scale_x_continuous(labels = scales::number_format(big.mark = ',',
                                                    decimal.mark = '.'))
setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
ggsave("One Year Out-of-Sample Returns (1999-2000).png", device = 'png', units = 'in',
       width = 10, height = 7)

all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(rt = mean(rt)*(252)) %>% 
  mutate(vol = mean(vol)*(252^(1/2))) %>% 
  mutate(sharpe = mean(sharpe)*(252^(1/2))) %>% 
  distinct(model, .keep_all = TRUE) %>% 
  select(model, rt, vol, sharpe) %>% 
  purrr::set_names("Model", "Annualized Return",
                   "Annualized Vol", "Annualized Sharpe") -> model_one_year

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
model_one_year %>% writexl::write_xlsx("One Year Out-of-Sample (1999-2000).xlsx")

# SHARPE DIFFERENCES
all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(n = seq(1, nrow(all_outsample_statistics)/9)) %>% 
  dplyr::select(n, sharpe, model) %>% 
  spread(model, sharpe) -> sharpe_outsample

(pvalue_sharpe_statistics = tibble(
  "GARCH Gaussian P-Value" = (t.test(x = sharpe_outsample$`GARCH Gaussian`,
                                     y = sharpe_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "GARCH Skew Gaussian P-Value" = (t.test(x = sharpe_outsample$`GARCH Skew Gaussian`,
                                          y = sharpe_outsample$`Actual Returns`,
                                          alternative = "greater",
                                          conf.level = 0.95))$p.value,
  "GARCH GED P-Value" = (t.test(x = sharpe_outsample$`GARCH GED`,
                                y = sharpe_outsample$`Actual Returns`,
                                alternative = "greater",
                                conf.level = 0.95))$p.value,
  "GARCH Skew GED P-Value" = (t.test(x = sharpe_outsample$`GARCH Skew GED`,
                                     y = sharpe_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "APARCH Gaussian P-Value" = (t.test(x = sharpe_outsample$`APARCH Gaussian`,
                                      y = sharpe_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
  "APARCH Skew Gaussian P-Value" = (t.test(x = sharpe_outsample$`APARCH Skew Gaussian`,
                                           y = sharpe_outsample$`Actual Returns`,
                                           alternative = "greater",
                                           conf.level = 0.95))$p.value,
  "APARCH GED P-Value" = (t.test(x = sharpe_outsample$`APARCH GED`,
                                 y = sharpe_outsample$`Actual Returns`,
                                 alternative = "greater",
                                 conf.level = 0.95))$p.value,
  "APARCH Skew GED P-Value" = (t.test(x = sharpe_outsample$`APARCH Skew GED`,
                                      y = sharpe_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
) %>% 
    gather(Model, Sharpe))

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
pvalue_sharpe_statistics %>% writexl::write_xlsx("One Year Sharpe Difference Statistics (1999-2000).xlsx")

# MEAN RETURN DIFFERENCES
all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(n = seq(1, nrow(all_outsample_statistics)/9)) %>% 
  dplyr::select(n, rt, model) %>% 
  spread(model, rt) -> rt_outsample

(pvalue_rt_statistics = tibble(
  "GARCH Gaussian P-Value" = (t.test(x = rt_outsample$`GARCH Gaussian`,
                                     y = rt_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "GARCH Skew Gaussian P-Value" = (t.test(x = rt_outsample$`GARCH Skew Gaussian`,
                                          y = rt_outsample$`Actual Returns`,
                                          alternative = "greater",
                                          conf.level = 0.95))$p.value,
  "GARCH GED P-Value" = (t.test(x = rt_outsample$`GARCH GED`,
                                y = rt_outsample$`Actual Returns`,
                                alternative = "greater",
                                conf.level = 0.95))$p.value,
  "GARCH Skew GED P-Value" = (t.test(x = rt_outsample$`GARCH Skew GED`,
                                     y = rt_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "APARCH Gaussian P-Value" = (t.test(x = rt_outsample$`APARCH Gaussian`,
                                      y = rt_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
  "APARCH Skew Gaussian P-Value" = (t.test(x = rt_outsample$`APARCH Skew Gaussian`,
                                           y = rt_outsample$`Actual Returns`,
                                           alternative = "greater",
                                           conf.level = 0.95))$p.value,
  "APARCH GED P-Value" = (t.test(x = rt_outsample$`APARCH GED`,
                                 y = rt_outsample$`Actual Returns`,
                                 alternative = "greater",
                                 conf.level = 0.95))$p.value,
  "APARCH Skew GED P-Value" = (t.test(x = rt_outsample$`APARCH Skew GED`,
                                      y = rt_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
) %>% 
    gather(Model, `Average Return`))

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
pvalue_rt_statistics %>%
  writexl::write_xlsx("One Year Average Return Difference Statistics (1999-2000).xlsx")

##### INDÚSTRIAS - 24 MESES #####
returns <- read_csv("C:/Users/Dinho Urbano/Desktop/TCC/Dados Financeiros/17 Industries.csv") %>% 
  mutate(date = date %>% ymd()) %>%
  mutate_if(is.numeric, ~ case_when(. == -99.99 ~ NA %>% as.numeric(),
                                    TRUE ~ .)) %>% 
  mutate_if(is.numeric, ~ . / 100) %>% 
  dplyr::filter(date >= "2001-01-01" & date <= "2002-12-31") %>% 
  merge(factors_df %>% dplyr::select(date, RF), by='date') %>% 
  mutate_if(is.numeric, ~ (1 + .) / (1 + RF) - 1) %>% 
  dplyr::select(-RF) %>% 
  gather(id, value, -date) %>% 
  spread(id, value) %>% 
  dplyr::select(-Cnstr)


# UTILIZAÇÃO DA FUNÇÃO NOS PORTFÓLIOS

for (i in 1:nrow(sample_numbers)){
  
  new_outsample_statistics = rbind(
    outsample_statistics(garch_norm, unlist(sample_numbers[i,])),
    outsample_statistics(garch_snorm, unlist(sample_numbers[i,])),
    outsample_statistics(garch_ged, unlist(sample_numbers[i,])),
    outsample_statistics(garch_sged, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_norm, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_snorm, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_ged, unlist(sample_numbers[i,])),
    outsample_statistics(aparch_sged, unlist(sample_numbers[i,])),
    outsample_statistics(actual_returns, unlist(sample_numbers[i,]))
  )
  
  if (i == 1){
    all_outsample_statistics = new_outsample_statistics
  } else {
    all_outsample_statistics = all_outsample_statistics %>% rbind(new_outsample_statistics)
  }
  
  setTxtProgressBar(txtProgressBar(
    min = 1, max = nrow(sample_numbers), initial = 0, style = 3), i)
}

models_names = c("GARCH Gaussian", "GARCH Skew Gaussian", "GARCH GED", "GARCH Skew GED",
                 "APARCH Gaussian", "APARCH Skew Gaussian", "APARCH GED", "APARCH Skew GED",
                 "Actual Returns")

all_outsample_statistics %>% 
  select(sharpe, model) %>%
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(x=sharpe*(252^(1/2)), y=model, fill=sharpe)) +
  geom_boxplot(color = 'black') +
  # geom_jitter(color="gray", size=0.4, alpha=0.3) +
  theme_minimal() +
  labs(title = 'Sharpe Out-of-Sample',
       subtitle = 'In-Sample: 1999-01-01 to 2000-12-31 \nOut-of-Sample: 2001-01-01 to 2002-12-31',
       x = 'Annualized Sharpe',
       y = '') +
  theme(legend.position = 'none') +
  scale_x_continuous(labels = scales::number_format(big.mark = ',',
                                                    decimal.mark = '.'))
setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
ggsave("Two Years Out-of-Sample Sharpe (1999-2000).png", device = 'png', units = 'in',
       width = 10, height = 7)

all_outsample_statistics %>% 
  select(rt, model) %>%
  mutate(model = factor(model, levels = rev(models_names))) %>% 
  ggplot(aes(x=rt*(252), y=model, fill=rt)) +
  geom_boxplot(color = 'black') +
  # geom_jitter(color="gray", size=0.4, alpha=0.3) +
  theme_minimal() +
  labs(title = 'Returns Out-of-Sample',
       subtitle = 'In-Sample: 1999-01-01 to 2000-12-31 \nOut-of-Sample: 2001-01-01 to 2002-12-31',
       x = 'Annualized Returns',
       y = '') +
  theme(legend.position = 'none') +
  scale_x_continuous(labels = scales::number_format(big.mark = ',',
                                                    decimal.mark = '.'))
setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
ggsave("Two Years Out-of-Sample Returns (1999-2000).png", device = 'png', units = 'in',
       width = 10, height = 7)

all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(rt = mean(rt)*(252)) %>% 
  mutate(vol = mean(vol)*(252^(1/2))) %>% 
  mutate(sharpe = mean(sharpe)*(252^(1/2))) %>% 
  distinct(model, .keep_all = TRUE) %>% 
  select(model, rt, vol, sharpe) %>% 
  purrr::set_names("Model", "Annualized Return",
                   "Annualized Vol", "Annualized Sharpe") -> model_two_years

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
model_two_years %>% writexl::write_xlsx("Two Years Out-of-Sample (1999-2000).xlsx")

# SHARPE DIFFERENCES
all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(n = seq(1, nrow(all_outsample_statistics)/9)) %>% 
  dplyr::select(n, sharpe, model) %>% 
  spread(model, sharpe) -> sharpe_outsample

(pvalue_sharpe_statistics = tibble(
  "GARCH Gaussian P-Value" = (t.test(x = sharpe_outsample$`GARCH Gaussian`,
                                     y = sharpe_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "GARCH Skew Gaussian P-Value" = (t.test(x = sharpe_outsample$`GARCH Skew Gaussian`,
                                          y = sharpe_outsample$`Actual Returns`,
                                          alternative = "greater",
                                          conf.level = 0.95))$p.value,
  "GARCH GED P-Value" = (t.test(x = sharpe_outsample$`GARCH GED`,
                                y = sharpe_outsample$`Actual Returns`,
                                alternative = "greater",
                                conf.level = 0.95))$p.value,
  "GARCH Skew GED P-Value" = (t.test(x = sharpe_outsample$`GARCH Skew GED`,
                                     y = sharpe_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "APARCH Gaussian P-Value" = (t.test(x = sharpe_outsample$`APARCH Gaussian`,
                                      y = sharpe_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
  "APARCH Skew Gaussian P-Value" = (t.test(x = sharpe_outsample$`APARCH Skew Gaussian`,
                                           y = sharpe_outsample$`Actual Returns`,
                                           alternative = "greater",
                                           conf.level = 0.95))$p.value,
  "APARCH GED P-Value" = (t.test(x = sharpe_outsample$`APARCH GED`,
                                 y = sharpe_outsample$`Actual Returns`,
                                 alternative = "greater",
                                 conf.level = 0.95))$p.value,
  "APARCH Skew GED P-Value" = (t.test(x = sharpe_outsample$`APARCH Skew GED`,
                                      y = sharpe_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
) %>% 
    gather(Model, Sharpe))

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
pvalue_sharpe_statistics %>% writexl::write_xlsx("Two Years Sharpe Difference Statistics (1999-2000).xlsx")

# MEAN RETURN DIFFERENCES
all_outsample_statistics %>% 
  group_by(model) %>% 
  mutate(n = seq(1, nrow(all_outsample_statistics)/9)) %>% 
  dplyr::select(n, rt, model) %>% 
  spread(model, rt) -> rt_outsample

(pvalue_rt_statistics = tibble(
  "GARCH Gaussian P-Value" = (t.test(x = rt_outsample$`GARCH Gaussian`,
                                     y = rt_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "GARCH Skew Gaussian P-Value" = (t.test(x = rt_outsample$`GARCH Skew Gaussian`,
                                          y = rt_outsample$`Actual Returns`,
                                          alternative = "greater",
                                          conf.level = 0.95))$p.value,
  "GARCH GED P-Value" = (t.test(x = rt_outsample$`GARCH GED`,
                                y = rt_outsample$`Actual Returns`,
                                alternative = "greater",
                                conf.level = 0.95))$p.value,
  "GARCH Skew GED P-Value" = (t.test(x = rt_outsample$`GARCH Skew GED`,
                                     y = rt_outsample$`Actual Returns`,
                                     alternative = "greater",
                                     conf.level = 0.95))$p.value,
  "APARCH Gaussian P-Value" = (t.test(x = rt_outsample$`APARCH Gaussian`,
                                      y = rt_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
  "APARCH Skew Gaussian P-Value" = (t.test(x = rt_outsample$`APARCH Skew Gaussian`,
                                           y = rt_outsample$`Actual Returns`,
                                           alternative = "greater",
                                           conf.level = 0.95))$p.value,
  "APARCH GED P-Value" = (t.test(x = rt_outsample$`APARCH GED`,
                                 y = rt_outsample$`Actual Returns`,
                                 alternative = "greater",
                                 conf.level = 0.95))$p.value,
  "APARCH Skew GED P-Value" = (t.test(x = rt_outsample$`APARCH Skew GED`,
                                      y = rt_outsample$`Actual Returns`,
                                      alternative = "greater",
                                      conf.level = 0.95))$p.value,
) %>% 
    gather(Model, `Average Return`))

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/1999-2000/Efficient Frontier")
pvalue_rt_statistics %>%
  writexl::write_xlsx("Two Years Average Return Difference Statistics (1999-2000).xlsx")


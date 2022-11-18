# 2007-2008
library(tidyverse)
library(lubridate)
library(readr)
library(readxl)
library(rmgarch)
library(quantmod)
library(fGarch)
library(cowplot)
library(tseries)

# FATORES
setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores")
factors_df <- read_csv("Fama French 5 Factors.csv") %>% 
  mutate(date = date %>% ymd()) %>% 
  mutate_if(is.numeric, ~ . / 100) %>% 
  rename(MKT = `Mkt-RF`)

# INDÚSTRIAS
returns <- read_csv("C:/Users/Dinho Urbano/Desktop/TCC/Dados Financeiros/17 Industries.csv") %>% 
  mutate(date = date %>% ymd()) %>%
  mutate_if(is.numeric, ~ case_when(. == -99.99 ~ NA %>% as.numeric(),
                                    TRUE ~ .)) %>% 
  mutate_if(is.numeric, ~ . / 100) %>% 
  dplyr::filter(date >= "2007-01-01" & date <= "2008-12-31") %>% 
  merge(factors_df %>% dplyr::select(date, RF), by='date') %>% 
  mutate_if(is.numeric, ~ (1 + .) / (1 + RF) - 1) %>% 
  dplyr::select(-RF) %>% 
  gather(name, value, -date)

factors <- factors_df %>% dplyr::select(-RF) # o RF é retirado, pois não é mais necessário

setwd("C:/Users/Dinho Urbano/Desktop/TCC/Fatores/2007-2008")

#### GARCH NORM ####
returns_hat_df = NULL
statistics = NULL

for (k in 2:18){
  
  returns_df <- returns %>% 
    spread(name, value)
  
  returns_df %>%
    tibble() %>% 
    dplyr::select("date", k) %>%
    purrr::set_names('date', 'stock') %>% 
    inner_join(factors, by='date') %>% 
    na.omit() %>% 
    data.frame() %>% 
    column_to_rownames('date')-> stock_plus_factors
  
  # DCC GARCH Fit
  xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                     variance.model = list(garchOrder = c(1,1), model = 'fGARCH', submodel = 'GARCH'),
                     distribution.model = 'norm')
  uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
  spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
  cl = makePSOCKcluster(ncol(stock_plus_factors))
  multf = multifit(uspec, stock_plus_factors, cluster = cl)
  
  fit_yx = NULL
  n_try = 0
  
  dcc_function = function (){
    dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
    return(dcc_model) 
  }
  
  possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
  
  while (is.null(fit_yx)){
    fit_yx = possibly_dcc_function()
    n_try = n_try + 1
    print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
    if (n_try >= 20){
      break
    }
  }
  
  if (is.null(fit_yx)){
    print(paste0(colnames(returns_df)[k], " did not work."))
  } else {
    
    for (i in 1:nrow(stock_plus_factors)){
      fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
    }
    
    # DEFINIÇÃO DOS BETAS
    for (i in 1:nrow(stock_plus_factors)){
      new_date = (stock_plus_factors %>% rownames())[i]
      
      x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
      y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
      parameters <- solve(x, y)
      y.bar <- mean(stock_plus_factors[1] %>% unlist())
      x.bar <- colMeans(stock_plus_factors[2:6])
      intercept <- y.bar - x.bar %*% parameters
      
      new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
        purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
      rownames(new_betas) <- new_date
      
      if (i == 1){
        all_betas <- new_betas
      } else {
        all_betas <- rbind(all_betas, new_betas)
      }
      setTxtProgressBar(txtProgressBar(
        min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
    }
    
    # RETORNOS AJUSTADOS E TESTES
    test <- (stock_plus_factors %>% 
               purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
      merge(all_betas %>%
              purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
            by = 0) %>%
      rename(date = Row.names) %>% 
      mutate(stock_rt_hat = intercept_bt +
               MKT_rt * MKT_bt +
               SMB_rt * SMB_bt +
               HML_rt * HML_bt +
               RMW_rt * RMW_bt +
               CMA_rt * CMA_bt) %>%
      mutate(error = stock_rt_hat - stock_rt) %>% 
      mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
      mutate(mse = (stock_rt_hat - stock_rt)^2)
    
    if (is.null(returns_hat_df)){
      returns_hat_df <- test %>% 
        dplyr::select(date, stock_rt_hat) %>% 
        mutate(id = colnames(returns_df)[k])
    } else {
      returns_hat_df <- returns_hat_df %>% 
        rbind(test %>% 
                dplyr::select(date, stock_rt_hat) %>% 
                mutate(id = colnames(returns_df)[k]))
    }
    
    test$mae %>% mean() -> mae
    test$mse %>% mean() -> mse
    (shapiro.test(test$error))$p.value -> shapiro_test
    (jarque.bera.test(test$error))$p.value -> jarque_bera_test
    (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
    (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
    (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
    (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
    skewness(test$error)[1] -> skewness
    kurtosis(test$error)[1] -> kurtosis
    
    new_statistics <- tibble(name = colnames(returns_df)[k],
                             MSE = mse, MAE = mae,
                             `Shapiro Wild - P-Valor` = shapiro_test,
                             `Jarque Bera - P-Valor` = jarque_bera_test,
                             `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
                             `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
                             `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
                             `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
                             `Assimetria` = skewness,
                             `Curtose` = kurtosis) %>% 
      gather(id, value, -name)
    
    if (is.null(statistics)){
      statistics <- new_statistics
    } else {
      statistics <- rbind(statistics, new_statistics)
    }
    
    # GRÁFICOS
    # BETAS
    graphic_betas <- all_betas %>%
      rownames_to_column('date') %>% 
      mutate(date = date %>% ymd()) %>% 
      rename(Intercept = intercept) %>% 
      gather(id, value, -date) %>%
      mutate(type = 'normal') %>% 
      ggplot(aes(x=date, y=value)) +
      geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
      geom_line(aes(group=id), color="black", size=0.85) +
      theme_minimal() +
      facet_wrap(~ id) +
      labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - GARCH Gaussian Distribution (2007-2008)'),
           x = '',
           y = '') +
      theme(legend.title = element_blank(),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',',
        decimal.mark = '.')) +
      scale_x_date(labels = scales::date_format("%Y-%m"))
    
    # ERRO
    graphic_erro <- test %>% 
      dplyr::select(date, error) %>%
      purrr::set_names('date', 'error') %>%
      mutate(date = date %>% ymd()) %>% 
      ggplot(aes(x=date, y=error)) +
      geom_line(color = 'black') +
      labs(title = 'Model Residuals',
           x = '',
           y = '') +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',', decimal.mark = '.')) +
      theme_minimal()
    
    # ACF E PACF
    conf.level <- 0.95
    ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
    
    bacf <- acf(test$error)
    df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Residuals')
    
    bpacf <- pacf(test$error)
    df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Residuals')
    
    bacf <- acf(test$error^2)
    df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Squared Residuals')
    
    bpacf <- pacf(test$error^2)
    df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Squared Residuals')
    
    acf_and_pacf <- plyr::join_all(list(
      df_acf, df_pacf, df_acf_sq, df_pacf_sq),
      by = 'lag', type = 'full') %>%
      tidyr::gather(id, value, -lag)
    
    graphic_acf_and_pacf <- acf_and_pacf %>% 
      ggplot(aes(x=lag, y=value)) +
      geom_bar(fill='black', stat = 'identity') +
      facet_wrap(~ id, nrow=2) +
      geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
      geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
                                                        big.mark = ',',
                                                        decimal.mark = '.')) +
      labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
      theme_minimal()
    
    # JUNÇÃO DOS GRÁFICOS
    plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
              nrow = 3, rel_heights = c(2, 1, 2))
    ggsave(paste0(colnames(returns_df)[k], " - GARCH Distribution Gaussian (2007-2008).png"), device = 'png',
           width = 11, height = 16.5)
    
    setTxtProgressBar(txtProgressBar(
      min = 2, max = 18, initial = 0, style = 3), k)
  }
}

returns_hat_df %>% writexl::write_xlsx("Returns Hat - GARCH Dist Gaussian (2007-2008).xlsx")
statistics %>% writexl::write_xlsx("Statistics - GARCH Dist Gaussian (2007-2008).xlsx")

#### GARCH SNORM ####
returns_hat_df = NULL
statistics = NULL

for (k in 2:18){
  
  returns_df <- returns %>% 
    spread(name, value)
  
  returns_df %>%
    tibble() %>% 
    dplyr::select("date", k) %>%
    purrr::set_names('date', 'stock') %>% 
    inner_join(factors, by='date') %>% 
    na.omit() %>% 
    data.frame() %>% 
    column_to_rownames('date')-> stock_plus_factors
  
  # DCC GARCH Fit
  xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                     variance.model = list(garchOrder = c(1,1), model = 'fGARCH', submodel = 'GARCH'),
                     distribution.model = 'snorm')
  uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
  spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
  cl = makePSOCKcluster(ncol(stock_plus_factors))
  multf = multifit(uspec, stock_plus_factors, cluster = cl)
  
  fit_yx = NULL
  n_try = 0
  
  dcc_function = function (){
    dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
    return(dcc_model) 
  }
  
  possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
  
  while (is.null(fit_yx)){
    fit_yx = possibly_dcc_function()
    n_try = n_try + 1
    print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
    if (n_try >= 20){
      break
    }
  }
  
  if (is.null(fit_yx)){
    print(paste0(colnames(returns_df)[k], " did not work."))
  } else {
    
    for (i in 1:nrow(stock_plus_factors)){
      fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
    }
    
    # DEFINIÇÃO DOS BETAS
    for (i in 1:nrow(stock_plus_factors)){
      new_date = (stock_plus_factors %>% rownames())[i]
      
      x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
      y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
      parameters <- solve(x, y)
      y.bar <- mean(stock_plus_factors[1] %>% unlist())
      x.bar <- colMeans(stock_plus_factors[2:6])
      intercept <- y.bar - x.bar %*% parameters
      
      new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
        purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
      rownames(new_betas) <- new_date
      
      if (i == 1){
        all_betas <- new_betas
      } else {
        all_betas <- rbind(all_betas, new_betas)
      }
      setTxtProgressBar(txtProgressBar(
        min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
    }
    
    # RETORNOS AJUSTADOS E TESTES
    test <- (stock_plus_factors %>% 
               purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
      merge(all_betas %>%
              purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
            by = 0) %>%
      rename(date = Row.names) %>% 
      mutate(stock_rt_hat = intercept_bt +
               MKT_rt * MKT_bt +
               SMB_rt * SMB_bt +
               HML_rt * HML_bt +
               RMW_rt * RMW_bt +
               CMA_rt * CMA_bt) %>%
      mutate(error = stock_rt_hat - stock_rt) %>% 
      mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
      mutate(mse = (stock_rt_hat - stock_rt)^2)
    
    if (is.null(returns_hat_df)){
      returns_hat_df <- test %>% 
        dplyr::select(date, stock_rt_hat) %>% 
        mutate(id = colnames(returns_df)[k])
    } else {
      returns_hat_df <- returns_hat_df %>% 
        rbind(test %>% 
                dplyr::select(date, stock_rt_hat) %>% 
                mutate(id = colnames(returns_df)[k]))
    }
    
    test$mae %>% mean() -> mae
    test$mse %>% mean() -> mse
    (shapiro.test(test$error))$p.value -> shapiro_test
    (jarque.bera.test(test$error))$p.value -> jarque_bera_test
    (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
    (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
    (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
    (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
    skewness(test$error)[1] -> skewness
    kurtosis(test$error)[1] -> kurtosis
    
    new_statistics <- tibble(name = colnames(returns_df)[k],
                             MSE = mse, MAE = mae,
                             `Shapiro Wild - P-Valor` = shapiro_test,
                             `Jarque Bera - P-Valor` = jarque_bera_test,
                             `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
                             `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
                             `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
                             `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
                             `Assimetria` = skewness,
                             `Curtose` = kurtosis) %>% 
      gather(id, value, -name)
    
    if (is.null(statistics)){
      statistics <- new_statistics
    } else {
      statistics <- rbind(statistics, new_statistics)
    }
    
    # GRÁFICOS
    # BETAS
    graphic_betas <- all_betas %>%
      rownames_to_column('date') %>% 
      mutate(date = date %>% ymd()) %>% 
      rename(Intercept = intercept) %>% 
      gather(id, value, -date) %>%
      mutate(type = 'normal') %>% 
      ggplot(aes(x=date, y=value)) +
      geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
      geom_line(aes(group=id), color="black", size=0.85) +
      theme_minimal() +
      facet_wrap(~ id) +
      labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - GARCH Skew Gaussian Distribution (2007-2008)'),
           x = '',
           y = '') +
      theme(legend.title = element_blank(),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',',
        decimal.mark = '.')) +
      scale_x_date(labels = scales::date_format("%Y-%m"))
    
    # ERRO
    graphic_erro <- test %>% 
      dplyr::select(date, error) %>%
      purrr::set_names('date', 'error') %>%
      mutate(date = date %>% ymd()) %>% 
      ggplot(aes(x=date, y=error)) +
      geom_line(color = 'black') +
      labs(title = 'Model Residuals',
           x = '',
           y = '') +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',', decimal.mark = '.')) +
      theme_minimal()
    
    # ACF E PACF
    conf.level <- 0.95
    ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
    
    bacf <- acf(test$error)
    df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Residuals')
    
    bpacf <- pacf(test$error)
    df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Residuals')
    
    bacf <- acf(test$error^2)
    df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Squared Residuals')
    
    bpacf <- pacf(test$error^2)
    df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Squared Residuals')
    
    acf_and_pacf <- plyr::join_all(list(
      df_acf, df_pacf, df_acf_sq, df_pacf_sq),
      by = 'lag', type = 'full') %>%
      tidyr::gather(id, value, -lag)
    
    graphic_acf_and_pacf <- acf_and_pacf %>% 
      ggplot(aes(x=lag, y=value)) +
      geom_bar(fill='black', stat = 'identity') +
      facet_wrap(~ id, nrow=2) +
      geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
      geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
                                                        big.mark = ',',
                                                        decimal.mark = '.')) +
      labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
      theme_minimal()
    
    # JUNÇÃO DOS GRÁFICOS
    plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
              nrow = 3, rel_heights = c(2, 1, 2))
    ggsave(paste0(colnames(returns_df)[k], " - GARCH Distribution Skew Gaussian (2007-2008).png"), device = 'png',
           width = 11, height = 16.5)
    
    setTxtProgressBar(txtProgressBar(
      min = 2, max = 18, initial = 0, style = 3), k)
  }
}

returns_hat_df %>% writexl::write_xlsx("Returns Hat - GARCH Dist Skew Gaussian (2007-2008).xlsx")
statistics %>% writexl::write_xlsx("Statistics - GARCH Dist Skew Gaussian (2007-2008).xlsx")

# #### GARCH STD ####
# returns_hat_df = NULL
# statistics = NULL

# for (k in 2:18){
#   
#   returns_df <- returns %>% 
#     spread(name, value)
#   
#   returns_df %>%
#     tibble() %>% 
#     dplyr::select("date", k) %>%
#     purrr::set_names('date', 'stock') %>% 
#     inner_join(factors, by='date') %>% 
#     na.omit() %>% 
#     data.frame() %>% 
#     column_to_rownames('date')-> stock_plus_factors
#   
#   # DCC GARCH Fit
#   xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
#                      variance.model = list(garchOrder = c(1,1), model = 'fGARCH', submodel = 'GARCH'),
#                      distribution.model = 'std')
#   uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
#   spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
#   cl = makePSOCKcluster(ncol(stock_plus_factors))
#   multf = multifit(uspec, stock_plus_factors, cluster = cl)
#   
#   fit_yx = NULL
#   n_try = 0
#   
#   dcc_function = function (){
#     dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
#     return(dcc_model) 
#   }
#   
#   possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
#   
#   while (is.null(fit_yx)){
#     fit_yx = possibly_dcc_function()
#     n_try = n_try + 1
#     print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
#     if (n_try >= 20){
#       break
#     }
#   }
#   
#   if (is.null(fit_yx)){
#     print(paste0(colnames(returns_df)[k], " did not work."))
#   } else {
#     
#     for (i in 1:nrow(stock_plus_factors)){
#       fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
#     }
#     
#     # DEFINIÇÃO DOS BETAS
#     for (i in 1:nrow(stock_plus_factors)){
#       new_date = (stock_plus_factors %>% rownames())[i]
#       
#       x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
#       y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
#       parameters <- solve(x, y)
#       y.bar <- mean(stock_plus_factors[1] %>% unlist())
#       x.bar <- colMeans(stock_plus_factors[2:6])
#       intercept <- y.bar - x.bar %*% parameters
#       
#       new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
#         purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
#       rownames(new_betas) <- new_date
#       
#       if (i == 1){
#         all_betas <- new_betas
#       } else {
#         all_betas <- rbind(all_betas, new_betas)
#       }
#       setTxtProgressBar(txtProgressBar(
#         min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
#     }
#     
#     # RETORNOS AJUSTADOS E TESTES
#     test <- (stock_plus_factors %>% 
#                purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
#       merge(all_betas %>%
#               purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
#             by = 0) %>%
#       rename(date = Row.names) %>% 
#       mutate(stock_rt_hat = intercept_bt +
#                MKT_rt * MKT_bt +
#                SMB_rt * SMB_bt +
#                HML_rt * HML_bt +
#                RMW_rt * RMW_bt +
#                CMA_rt * CMA_bt) %>%
#       mutate(error = stock_rt_hat - stock_rt) %>% 
#       mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
#       mutate(mse = (stock_rt_hat - stock_rt)^2)
#     
#     if (is.null(returns_hat_df)){
#       returns_hat_df <- test %>% 
#         dplyr::select(date, stock_rt_hat) %>% 
#         mutate(id = colnames(returns_df)[k])
#     } else {
#       returns_hat_df <- returns_hat_df %>% 
#         rbind(test %>% 
#                 dplyr::select(date, stock_rt_hat) %>% 
#                 mutate(id = colnames(returns_df)[k]))
#     }
#     
#     test$mae %>% mean() -> mae
#     test$mse %>% mean() -> mse
#     (shapiro.test(test$error))$p.value -> shapiro_test
#     (jarque.bera.test(test$error))$p.value -> jarque_bera_test
#     (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
#     (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
#     (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
#     (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
#     skewness(test$error)[1] -> skewness
#     kurtosis(test$error)[1] -> kurtosis
#     
#     new_statistics <- tibble(name = colnames(returns_df)[k],
#                              MSE = mse, MAE = mae,
#                              `Shapiro Wild - P-Valor` = shapiro_test,
#                              `Jarque Bera - P-Valor` = jarque_bera_test,
#                              `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
#                              `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
#                              `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
#                              `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
#                              `Assimetria` = skewness,
#                              `Curtose` = kurtosis) %>% 
#       gather(id, value, -name)
#     
#     if (is.null(statistics)){
#       statistics <- new_statistics
#     } else {
#       statistics <- rbind(statistics, new_statistics)
#     }
#     
#     # GRÁFICOS
#     # BETAS
#     graphic_betas <- all_betas %>%
#       rownames_to_column('date') %>% 
#       mutate(date = date %>% ymd()) %>% 
#       rename(Intercept = intercept) %>% 
#       gather(id, value, -date) %>%
#       mutate(type = 'normal') %>% 
#       ggplot(aes(x=date, y=value)) +
#       geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
#       geom_line(aes(group=id), color="black", size=0.85) +
#       theme_minimal() +
#       facet_wrap(~ id) +
#       labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - GARCH T-Student Distribution (2007-2008)'),
#            x = '',
#            y = '') +
#       theme(legend.title = element_blank(),
#             legend.position = 'right',
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       scale_y_continuous(labels = scales::number_format(
#         accuracy = 0.01,
#         big.mark = ',',
#         decimal.mark = '.')) +
#       scale_x_date(labels = scales::date_format("%Y-%m"))
#     
#     # ERRO
#     graphic_erro <- test %>% 
#       dplyr::select(date, error) %>%
#       purrr::set_names('date', 'error') %>%
#       mutate(date = date %>% ymd()) %>% 
#       ggplot(aes(x=date, y=error)) +
#       geom_line(color = 'black') +
#       labs(title = 'Model Residuals',
#            x = '',
#            y = '') +
#       scale_y_continuous(labels = scales::number_format(
#         accuracy = 0.01,
#         big.mark = ',', decimal.mark = '.')) +
#       theme_minimal()
#     
#     # ACF E PACF
#     conf.level <- 0.95
#     ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
#     
#     bacf <- acf(test$error)
#     df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
#       purrr::set_names('lag', 'ACF Residuals')
#     
#     bpacf <- pacf(test$error)
#     df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
#       purrr::set_names('lag', 'PACF Residuals')
#     
#     bacf <- acf(test$error^2)
#     df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
#       purrr::set_names('lag', 'ACF Squared Residuals')
#     
#     bpacf <- pacf(test$error^2)
#     df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
#       purrr::set_names('lag', 'PACF Squared Residuals')
#     
#     acf_and_pacf <- plyr::join_all(list(
#       df_acf, df_pacf, df_acf_sq, df_pacf_sq),
#       by = 'lag', type = 'full') %>%
#       tidyr::gather(id, value, -lag)
#     
#     graphic_acf_and_pacf <- acf_and_pacf %>% 
#       ggplot(aes(x=lag, y=value)) +
#       geom_bar(fill='black', stat = 'identity') +
#       facet_wrap(~ id, nrow=2) +
#       geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
#       geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
#       scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
#                                                         big.mark = ',',
#                                                         decimal.mark = '.')) +
#       labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
#       theme_minimal()
#     
#     # JUNÇÃO DOS GRÁFICOS
#     plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
#               nrow = 3, rel_heights = c(2, 1, 2))
#     ggsave(paste0(colnames(returns_df)[k], " - GARCH Distribution T-Student (2007-2008).png"), device = 'png',
#            width = 11, height = 16.5)
#     
#     setTxtProgressBar(txtProgressBar(
#       min = 2, max = 18, initial = 0, style = 3), k)
#   }
# }
# 
# returns_hat_df %>% writexl::write_xlsx("Returns Hat - GARCH Dist T-Student (2007-2008).xlsx")
# statistics %>% writexl::write_xlsx("Statistics - GARCH Dist T-Student (2007-2008).xlsx")
# 
# #### GARCH SSTD ####
# returns_hat_df = NULL
# statistics = NULL

# for (k in 2:18){
#   
#   returns_df <- returns %>% 
#     spread(name, value)
#   
#   returns_df %>%
#     tibble() %>% 
#     dplyr::select("date", k) %>%
#     purrr::set_names('date', 'stock') %>% 
#     inner_join(factors, by='date') %>% 
#     na.omit() %>% 
#     data.frame() %>% 
#     column_to_rownames('date')-> stock_plus_factors
#   
#   # DCC GARCH Fit
#   xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
#                      variance.model = list(garchOrder = c(1,1), model = 'fGARCH', submodel = 'GARCH'),
#                      distribution.model = 'sstd')
#   uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
#   spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
#   cl = makePSOCKcluster(ncol(stock_plus_factors))
#   multf = multifit(uspec, stock_plus_factors, cluster = cl)
#   
#   fit_yx = NULL
#   n_try = 0
#   
#   dcc_function = function (){
#     dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
#     return(dcc_model) 
#   }
#   
#   possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
#   
#   while (is.null(fit_yx)){
#     fit_yx = possibly_dcc_function()
#     n_try = n_try + 1
#     print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
#     if (n_try >= 20){
#       break
#     }
#   }
#   
#   if (is.null(fit_yx)){
#     print(paste0(colnames(returns_df)[k], " did not work."))
#   } else {
#     
#     for (i in 1:nrow(stock_plus_factors)){
#       fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
#     }
#     
#     # DEFINIÇÃO DOS BETAS
#     for (i in 1:nrow(stock_plus_factors)){
#       new_date = (stock_plus_factors %>% rownames())[i]
#       
#       x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
#       y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
#       parameters <- solve(x, y)
#       y.bar <- mean(stock_plus_factors[1] %>% unlist())
#       x.bar <- colMeans(stock_plus_factors[2:6])
#       intercept <- y.bar - x.bar %*% parameters
#       
#       new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
#         purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
#       rownames(new_betas) <- new_date
#       
#       if (i == 1){
#         all_betas <- new_betas
#       } else {
#         all_betas <- rbind(all_betas, new_betas)
#       }
#       setTxtProgressBar(txtProgressBar(
#         min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
#     }
#     
#     # RETORNOS AJUSTADOS E TESTES
#     test <- (stock_plus_factors %>% 
#                purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
#       merge(all_betas %>%
#               purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
#             by = 0) %>%
#       rename(date = Row.names) %>% 
#       mutate(stock_rt_hat = intercept_bt +
#                MKT_rt * MKT_bt +
#                SMB_rt * SMB_bt +
#                HML_rt * HML_bt +
#                RMW_rt * RMW_bt +
#                CMA_rt * CMA_bt) %>%
#       mutate(error = stock_rt_hat - stock_rt) %>% 
#       mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
#       mutate(mse = (stock_rt_hat - stock_rt)^2)
#     
#     if (is.null(returns_hat_df)){
#       returns_hat_df <- test %>% 
#         dplyr::select(date, stock_rt_hat) %>% 
#         mutate(id = colnames(returns_df)[k])
#     } else {
#       returns_hat_df <- returns_hat_df %>% 
#         rbind(test %>% 
#                 dplyr::select(date, stock_rt_hat) %>% 
#                 mutate(id = colnames(returns_df)[k]))
#     }
#     
#     test$mae %>% mean() -> mae
#     test$mse %>% mean() -> mse
#     (shapiro.test(test$error))$p.value -> shapiro_test
#     (jarque.bera.test(test$error))$p.value -> jarque_bera_test
#     (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
#     (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
#     (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
#     (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
#     skewness(test$error)[1] -> skewness
#     kurtosis(test$error)[1] -> kurtosis
#     
#     new_statistics <- tibble(name = colnames(returns_df)[k],
#                              MSE = mse, MAE = mae,
#                              `Shapiro Wild - P-Valor` = shapiro_test,
#                              `Jarque Bera - P-Valor` = jarque_bera_test,
#                              `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
#                              `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
#                              `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
#                              `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
#                              `Assimetria` = skewness,
#                              `Curtose` = kurtosis) %>% 
#       gather(id, value, -name)
#     
#     if (is.null(statistics)){
#       statistics <- new_statistics
#     } else {
#       statistics <- rbind(statistics, new_statistics)
#     }
#     
#     # GRÁFICOS
#     # BETAS
#     graphic_betas <- all_betas %>%
#       rownames_to_column('date') %>% 
#       mutate(date = date %>% ymd()) %>% 
#       rename(Intercept = intercept) %>% 
#       gather(id, value, -date) %>%
#       mutate(type = 'normal') %>% 
#       ggplot(aes(x=date, y=value)) +
#       geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
#       geom_line(aes(group=id), color="black", size=0.85) +
#       theme_minimal() +
#       facet_wrap(~ id) +
#       labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - GARCH Skew T-Student Distribution (2007-2008)'),
#            x = '',
#            y = '') +
#       theme(legend.title = element_blank(),
#             legend.position = 'right',
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       scale_y_continuous(labels = scales::number_format(
#         accuracy = 0.01,
#         big.mark = ',',
#         decimal.mark = '.')) +
#       scale_x_date(labels = scales::date_format("%Y-%m"))
#     
#     # ERRO
#     graphic_erro <- test %>% 
#       dplyr::select(date, error) %>%
#       purrr::set_names('date', 'error') %>%
#       mutate(date = date %>% ymd()) %>% 
#       ggplot(aes(x=date, y=error)) +
#       geom_line(color = 'black') +
#       labs(title = 'Model Residuals',
#            x = '',
#            y = '') +
#       scale_y_continuous(labels = scales::number_format(
#         accuracy = 0.01,
#         big.mark = ',', decimal.mark = '.')) +
#       theme_minimal()
#     
#     # ACF E PACF
#     conf.level <- 0.95
#     ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
#     
#     bacf <- acf(test$error)
#     df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
#       purrr::set_names('lag', 'ACF Residuals')
#     
#     bpacf <- pacf(test$error)
#     df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
#       purrr::set_names('lag', 'PACF Residuals')
#     
#     bacf <- acf(test$error^2)
#     df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
#       purrr::set_names('lag', 'ACF Squared Residuals')
#     
#     bpacf <- pacf(test$error^2)
#     df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
#       purrr::set_names('lag', 'PACF Squared Residuals')
#     
#     acf_and_pacf <- plyr::join_all(list(
#       df_acf, df_pacf, df_acf_sq, df_pacf_sq),
#       by = 'lag', type = 'full') %>%
#       tidyr::gather(id, value, -lag)
#     
#     graphic_acf_and_pacf <- acf_and_pacf %>% 
#       ggplot(aes(x=lag, y=value)) +
#       geom_bar(fill='black', stat = 'identity') +
#       facet_wrap(~ id, nrow=2) +
#       geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
#       geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
#       scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
#                                                         big.mark = ',',
#                                                         decimal.mark = '.')) +
#       labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
#       theme_minimal()
#     
#     # JUNÇÃO DOS GRÁFICOS
#     plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
#               nrow = 3, rel_heights = c(2, 1, 2))
#     ggsave(paste0(colnames(returns_df)[k], " - GARCH Distribution Skew T-Student (2007-2008).png"), device = 'png',
#            width = 11, height = 16.5)
#     
#     setTxtProgressBar(txtProgressBar(
#       min = 2, max = 18, initial = 0, style = 3), k)
#   }
# }
# 
# returns_hat_df %>% writexl::write_xlsx("Returns Hat - GARCH Dist Skew T-Student (2007-2008).xlsx")
# statistics %>% writexl::write_xlsx("Statistics - GARCH Dist Skew T-Student (2007-2008).xlsx")

#### GARCH GED ####
returns_hat_df = NULL
statistics = NULL

for (k in 2:18){
  
  returns_df <- returns %>% 
    spread(name, value)
  
  returns_df %>%
    tibble() %>% 
    dplyr::select("date", k) %>%
    purrr::set_names('date', 'stock') %>% 
    inner_join(factors, by='date') %>% 
    na.omit() %>% 
    data.frame() %>% 
    column_to_rownames('date')-> stock_plus_factors
  
  # DCC GARCH Fit
  xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                     variance.model = list(garchOrder = c(1,1), model = 'fGARCH', submodel = 'GARCH'),
                     distribution.model = 'ged')
  uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
  spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
  cl = makePSOCKcluster(ncol(stock_plus_factors))
  multf = multifit(uspec, stock_plus_factors, cluster = cl)
  
  fit_yx = NULL
  n_try = 0
  
  dcc_function = function (){
    dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
    return(dcc_model) 
  }
  
  possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
  
  while (is.null(fit_yx)){
    fit_yx = possibly_dcc_function()
    n_try = n_try + 1
    print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
    if (n_try >= 20){
      break
    }
  }
  
  if (is.null(fit_yx)){
    print(paste0(colnames(returns_df)[k], " did not work."))
  } else {
    
    for (i in 1:nrow(stock_plus_factors)){
      fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
    }
    
    # DEFINIÇÃO DOS BETAS
    for (i in 1:nrow(stock_plus_factors)){
      new_date = (stock_plus_factors %>% rownames())[i]
      
      x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
      y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
      parameters <- solve(x, y)
      y.bar <- mean(stock_plus_factors[1] %>% unlist())
      x.bar <- colMeans(stock_plus_factors[2:6])
      intercept <- y.bar - x.bar %*% parameters
      
      new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
        purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
      rownames(new_betas) <- new_date
      
      if (i == 1){
        all_betas <- new_betas
      } else {
        all_betas <- rbind(all_betas, new_betas)
      }
      setTxtProgressBar(txtProgressBar(
        min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
    }
    
    # RETORNOS AJUSTADOS E TESTES
    test <- (stock_plus_factors %>% 
               purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
      merge(all_betas %>%
              purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
            by = 0) %>%
      rename(date = Row.names) %>% 
      mutate(stock_rt_hat = intercept_bt +
               MKT_rt * MKT_bt +
               SMB_rt * SMB_bt +
               HML_rt * HML_bt +
               RMW_rt * RMW_bt +
               CMA_rt * CMA_bt) %>%
      mutate(error = stock_rt_hat - stock_rt) %>% 
      mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
      mutate(mse = (stock_rt_hat - stock_rt)^2)
    
    if (is.null(returns_hat_df)){
      returns_hat_df <- test %>% 
        dplyr::select(date, stock_rt_hat) %>% 
        mutate(id = colnames(returns_df)[k])
    } else {
      returns_hat_df <- returns_hat_df %>% 
        rbind(test %>% 
                dplyr::select(date, stock_rt_hat) %>% 
                mutate(id = colnames(returns_df)[k]))
    }
    
    test$mae %>% mean() -> mae
    test$mse %>% mean() -> mse
    (shapiro.test(test$error))$p.value -> shapiro_test
    (jarque.bera.test(test$error))$p.value -> jarque_bera_test
    (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
    (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
    (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
    (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
    skewness(test$error)[1] -> skewness
    kurtosis(test$error)[1] -> kurtosis
    
    new_statistics <- tibble(name = colnames(returns_df)[k],
                             MSE = mse, MAE = mae,
                             `Shapiro Wild - P-Valor` = shapiro_test,
                             `Jarque Bera - P-Valor` = jarque_bera_test,
                             `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
                             `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
                             `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
                             `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
                             `Assimetria` = skewness,
                             `Curtose` = kurtosis) %>% 
      gather(id, value, -name)
    
    if (is.null(statistics)){
      statistics <- new_statistics
    } else {
      statistics <- rbind(statistics, new_statistics)
    }
    
    # GRÁFICOS
    # BETAS
    graphic_betas <- all_betas %>%
      rownames_to_column('date') %>% 
      mutate(date = date %>% ymd()) %>% 
      rename(Intercept = intercept) %>% 
      gather(id, value, -date) %>%
      mutate(type = 'normal') %>% 
      ggplot(aes(x=date, y=value)) +
      geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
      geom_line(aes(group=id), color="black", size=0.85) +
      theme_minimal() +
      facet_wrap(~ id) +
      labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - GARCH GED Distribution (2007-2008)'),
           x = '',
           y = '') +
      theme(legend.title = element_blank(),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',',
        decimal.mark = '.')) +
      scale_x_date(labels = scales::date_format("%Y-%m"))
    
    # ERRO
    graphic_erro <- test %>% 
      dplyr::select(date, error) %>%
      purrr::set_names('date', 'error') %>%
      mutate(date = date %>% ymd()) %>% 
      ggplot(aes(x=date, y=error)) +
      geom_line(color = 'black') +
      labs(title = 'Model Residuals',
           x = '',
           y = '') +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',', decimal.mark = '.')) +
      theme_minimal()
    
    # ACF E PACF
    conf.level <- 0.95
    ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
    
    bacf <- acf(test$error)
    df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Residuals')
    
    bpacf <- pacf(test$error)
    df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Residuals')
    
    bacf <- acf(test$error^2)
    df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Squared Residuals')
    
    bpacf <- pacf(test$error^2)
    df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Squared Residuals')
    
    acf_and_pacf <- plyr::join_all(list(
      df_acf, df_pacf, df_acf_sq, df_pacf_sq),
      by = 'lag', type = 'full') %>%
      tidyr::gather(id, value, -lag)
    
    graphic_acf_and_pacf <- acf_and_pacf %>% 
      ggplot(aes(x=lag, y=value)) +
      geom_bar(fill='black', stat = 'identity') +
      facet_wrap(~ id, nrow=2) +
      geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
      geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
                                                        big.mark = ',',
                                                        decimal.mark = '.')) +
      labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
      theme_minimal()
    
    # JUNÇÃO DOS GRÁFICOS
    plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
              nrow = 3, rel_heights = c(2, 1, 2))
    ggsave(paste0(colnames(returns_df)[k], " - GARCH Distribution GED (2007-2008).png"), device = 'png',
           width = 11, height = 16.5)
    
    setTxtProgressBar(txtProgressBar(
      min = 2, max = 18, initial = 0, style = 3), k)
  }
}

returns_hat_df %>% writexl::write_xlsx("Returns Hat - GARCH Dist GED (2007-2008).xlsx")
statistics %>% writexl::write_xlsx("Statistics - GARCH Dist GED (2007-2008).xlsx")

#### GARCH SGED ####
returns_hat_df = NULL
statistics = NULL

for (k in 2:18){
  
  returns_df <- returns %>% 
    spread(name, value)
  
  returns_df %>%
    tibble() %>% 
    dplyr::select("date", k) %>%
    purrr::set_names('date', 'stock') %>% 
    inner_join(factors, by='date') %>% 
    na.omit() %>% 
    data.frame() %>% 
    column_to_rownames('date')-> stock_plus_factors
  
  # DCC GARCH Fit
  xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                     variance.model = list(garchOrder = c(1,1), model = 'fGARCH', submodel = 'GARCH'),
                     distribution.model = 'sged')
  uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
  spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
  cl = makePSOCKcluster(ncol(stock_plus_factors))
  multf = multifit(uspec, stock_plus_factors, cluster = cl)
  
  fit_yx = NULL
  n_try = 0
  
  dcc_function = function (){
    dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
    return(dcc_model) 
  }
  
  possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
  
  while (is.null(fit_yx)){
    fit_yx = possibly_dcc_function()
    n_try = n_try + 1
    print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
    if (n_try >= 20){
      break
    }
  }
  
  if (is.null(fit_yx)){
    print(paste0(colnames(returns_df)[k], " did not work."))
  } else {
    
    for (i in 1:nrow(stock_plus_factors)){
      fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
    }
    
    # DEFINIÇÃO DOS BETAS
    for (i in 1:nrow(stock_plus_factors)){
      new_date = (stock_plus_factors %>% rownames())[i]
      
      x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
      y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
      parameters <- solve(x, y)
      y.bar <- mean(stock_plus_factors[1] %>% unlist())
      x.bar <- colMeans(stock_plus_factors[2:6])
      intercept <- y.bar - x.bar %*% parameters
      
      new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
        purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
      rownames(new_betas) <- new_date
      
      if (i == 1){
        all_betas <- new_betas
      } else {
        all_betas <- rbind(all_betas, new_betas)
      }
      setTxtProgressBar(txtProgressBar(
        min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
    }
    
    # RETORNOS AJUSTADOS E TESTES
    test <- (stock_plus_factors %>% 
               purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
      merge(all_betas %>%
              purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
            by = 0) %>%
      rename(date = Row.names) %>% 
      mutate(stock_rt_hat = intercept_bt +
               MKT_rt * MKT_bt +
               SMB_rt * SMB_bt +
               HML_rt * HML_bt +
               RMW_rt * RMW_bt +
               CMA_rt * CMA_bt) %>%
      mutate(error = stock_rt_hat - stock_rt) %>% 
      mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
      mutate(mse = (stock_rt_hat - stock_rt)^2)
    
    if (is.null(returns_hat_df)){
      returns_hat_df <- test %>% 
        dplyr::select(date, stock_rt_hat) %>% 
        mutate(id = colnames(returns_df)[k])
    } else {
      returns_hat_df <- returns_hat_df %>% 
        rbind(test %>% 
                dplyr::select(date, stock_rt_hat) %>% 
                mutate(id = colnames(returns_df)[k]))
    }
    
    test$mae %>% mean() -> mae
    test$mse %>% mean() -> mse
    (shapiro.test(test$error))$p.value -> shapiro_test
    (jarque.bera.test(test$error))$p.value -> jarque_bera_test
    (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
    (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
    (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
    (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
    skewness(test$error)[1] -> skewness
    kurtosis(test$error)[1] -> kurtosis
    
    new_statistics <- tibble(name = colnames(returns_df)[k],
                             MSE = mse, MAE = mae,
                             `Shapiro Wild - P-Valor` = shapiro_test,
                             `Jarque Bera - P-Valor` = jarque_bera_test,
                             `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
                             `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
                             `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
                             `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
                             `Assimetria` = skewness,
                             `Curtose` = kurtosis) %>% 
      gather(id, value, -name)
    
    if (is.null(statistics)){
      statistics <- new_statistics
    } else {
      statistics <- rbind(statistics, new_statistics)
    }
    
    # GRÁFICOS
    # BETAS
    graphic_betas <- all_betas %>%
      rownames_to_column('date') %>% 
      mutate(date = date %>% ymd()) %>% 
      rename(Intercept = intercept) %>% 
      gather(id, value, -date) %>%
      mutate(type = 'normal') %>% 
      ggplot(aes(x=date, y=value)) +
      geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
      geom_line(aes(group=id), color="black", size=0.85) +
      theme_minimal() +
      facet_wrap(~ id) +
      labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - GARCH Skew GED Distribution (2007-2008)'),
           x = '',
           y = '') +
      theme(legend.title = element_blank(),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',',
        decimal.mark = '.')) +
      scale_x_date(labels = scales::date_format("%Y-%m"))
    
    # ERRO
    graphic_erro <- test %>% 
      dplyr::select(date, error) %>%
      purrr::set_names('date', 'error') %>%
      mutate(date = date %>% ymd()) %>% 
      ggplot(aes(x=date, y=error)) +
      geom_line(color = 'black') +
      labs(title = 'Model Residuals',
           x = '',
           y = '') +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',', decimal.mark = '.')) +
      theme_minimal()
    
    # ACF E PACF
    conf.level <- 0.95
    ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
    
    bacf <- acf(test$error)
    df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Residuals')
    
    bpacf <- pacf(test$error)
    df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Residuals')
    
    bacf <- acf(test$error^2)
    df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Squared Residuals')
    
    bpacf <- pacf(test$error^2)
    df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Squared Residuals')
    
    acf_and_pacf <- plyr::join_all(list(
      df_acf, df_pacf, df_acf_sq, df_pacf_sq),
      by = 'lag', type = 'full') %>%
      tidyr::gather(id, value, -lag)
    
    graphic_acf_and_pacf <- acf_and_pacf %>% 
      ggplot(aes(x=lag, y=value)) +
      geom_bar(fill='black', stat = 'identity') +
      facet_wrap(~ id, nrow=2) +
      geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
      geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
                                                        big.mark = ',',
                                                        decimal.mark = '.')) +
      labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
      theme_minimal()
    
    # JUNÇÃO DOS GRÁFICOS
    plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
              nrow = 3, rel_heights = c(2, 1, 2))
    ggsave(paste0(colnames(returns_df)[k], " - GARCH Distribution Skew GED (2007-2008).png"), device = 'png',
           width = 11, height = 16.5)
    
    setTxtProgressBar(txtProgressBar(
      min = 2, max = 18, initial = 0, style = 3), k)
  }
}

returns_hat_df %>% writexl::write_xlsx("Returns Hat - GARCH Dist Skew GED (2007-2008).xlsx")
statistics %>% writexl::write_xlsx("Statistics - GARCH Dist Skew GED (2007-2008).xlsx")

#### APARCH NORM ####
returns_hat_df = NULL
statistics = NULL

for (k in 2:18){
  
  returns_df <- returns %>% 
    spread(name, value)
  
  returns_df %>%
    tibble() %>% 
    dplyr::select("date", k) %>%
    purrr::set_names('date', 'stock') %>% 
    inner_join(factors, by='date') %>% 
    na.omit() %>% 
    data.frame() %>% 
    column_to_rownames('date')-> stock_plus_factors
  
  # DCC GARCH Fit
  xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                     variance.model = list(garchOrder = c(1,1), model = 'apARCH'),
                     distribution.model = 'norm')
  uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
  spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
  cl = makePSOCKcluster(ncol(stock_plus_factors))
  multf = multifit(uspec, stock_plus_factors, cluster = cl)
  
  fit_yx = NULL
  n_try = 0
  
  dcc_function = function (){
    dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
    return(dcc_model) 
  }
  
  possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
  
  while (is.null(fit_yx)){
    fit_yx = possibly_dcc_function()
    n_try = n_try + 1
    print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
    if (n_try >= 20){
      break
    }
  }
  
  if (is.null(fit_yx)){
    print(paste0(colnames(returns_df)[k], " did not work."))
  } else {
    
    for (i in 1:nrow(stock_plus_factors)){
      fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
    }
    
    # DEFINIÇÃO DOS BETAS
    for (i in 1:nrow(stock_plus_factors)){
      new_date = (stock_plus_factors %>% rownames())[i]
      
      x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
      y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
      parameters <- solve(x, y)
      y.bar <- mean(stock_plus_factors[1] %>% unlist())
      x.bar <- colMeans(stock_plus_factors[2:6])
      intercept <- y.bar - x.bar %*% parameters
      
      new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
        purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
      rownames(new_betas) <- new_date
      
      if (i == 1){
        all_betas <- new_betas
      } else {
        all_betas <- rbind(all_betas, new_betas)
      }
      setTxtProgressBar(txtProgressBar(
        min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
    }
    
    # RETORNOS AJUSTADOS E TESTES
    test <- (stock_plus_factors %>% 
               purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
      merge(all_betas %>%
              purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
            by = 0) %>%
      rename(date = Row.names) %>% 
      mutate(stock_rt_hat = intercept_bt +
               MKT_rt * MKT_bt +
               SMB_rt * SMB_bt +
               HML_rt * HML_bt +
               RMW_rt * RMW_bt +
               CMA_rt * CMA_bt) %>%
      mutate(error = stock_rt_hat - stock_rt) %>% 
      mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
      mutate(mse = (stock_rt_hat - stock_rt)^2)
    
    if (is.null(returns_hat_df)){
      returns_hat_df <- test %>% 
        dplyr::select(date, stock_rt_hat) %>% 
        mutate(id = colnames(returns_df)[k])
    } else {
      returns_hat_df <- returns_hat_df %>% 
        rbind(test %>% 
                dplyr::select(date, stock_rt_hat) %>% 
                mutate(id = colnames(returns_df)[k]))
    }
    
    test$mae %>% mean() -> mae
    test$mse %>% mean() -> mse
    (shapiro.test(test$error))$p.value -> shapiro_test
    (jarque.bera.test(test$error))$p.value -> jarque_bera_test
    (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
    (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
    (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
    (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
    skewness(test$error)[1] -> skewness
    kurtosis(test$error)[1] -> kurtosis
    
    new_statistics <- tibble(name = colnames(returns_df)[k],
                             MSE = mse, MAE = mae,
                             `Shapiro Wild - P-Valor` = shapiro_test,
                             `Jarque Bera - P-Valor` = jarque_bera_test,
                             `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
                             `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
                             `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
                             `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
                             `Assimetria` = skewness,
                             `Curtose` = kurtosis) %>% 
      gather(id, value, -name)
    
    if (is.null(statistics)){
      statistics <- new_statistics
    } else {
      statistics <- rbind(statistics, new_statistics)
    }
    
    # GRÁFICOS
    # BETAS
    graphic_betas <- all_betas %>%
      rownames_to_column('date') %>% 
      mutate(date = date %>% ymd()) %>% 
      rename(Intercept = intercept) %>% 
      gather(id, value, -date) %>%
      mutate(type = 'normal') %>% 
      ggplot(aes(x=date, y=value)) +
      geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
      geom_line(aes(group=id), color="black", size=0.85) +
      theme_minimal() +
      facet_wrap(~ id) +
      labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - APARCH Gaussian Distribution (2007-2008)'),
           x = '',
           y = '') +
      theme(legend.title = element_blank(),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',',
        decimal.mark = '.')) +
      scale_x_date(labels = scales::date_format("%Y-%m"))
    
    # ERRO
    graphic_erro <- test %>% 
      dplyr::select(date, error) %>%
      purrr::set_names('date', 'error') %>%
      mutate(date = date %>% ymd()) %>% 
      ggplot(aes(x=date, y=error)) +
      geom_line(color = 'black') +
      labs(title = 'Model Residuals',
           x = '',
           y = '') +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',', decimal.mark = '.')) +
      theme_minimal()
    
    # ACF E PACF
    conf.level <- 0.95
    ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
    
    bacf <- acf(test$error)
    df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Residuals')
    
    bpacf <- pacf(test$error)
    df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Residuals')
    
    bacf <- acf(test$error^2)
    df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Squared Residuals')
    
    bpacf <- pacf(test$error^2)
    df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Squared Residuals')
    
    acf_and_pacf <- plyr::join_all(list(
      df_acf, df_pacf, df_acf_sq, df_pacf_sq),
      by = 'lag', type = 'full') %>%
      tidyr::gather(id, value, -lag)
    
    graphic_acf_and_pacf <- acf_and_pacf %>% 
      ggplot(aes(x=lag, y=value)) +
      geom_bar(fill='black', stat = 'identity') +
      facet_wrap(~ id, nrow=2) +
      geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
      geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
                                                        big.mark = ',',
                                                        decimal.mark = '.')) +
      labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
      theme_minimal()
    
    # JUNÇÃO DOS GRÁFICOS
    plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
              nrow = 3, rel_heights = c(2, 1, 2))
    ggsave(paste0(colnames(returns_df)[k], " - APARCH Distribution Gaussian (2007-2008).png"), device = 'png',
           width = 11, height = 16.5)
    
    setTxtProgressBar(txtProgressBar(
      min = 2, max = 18, initial = 0, style = 3), k)
  }
}

returns_hat_df %>% writexl::write_xlsx("Returns Hat - APARCH Dist Gaussian (2007-2008).xlsx")
statistics %>% writexl::write_xlsx("Statistics - APARCH Dist Gaussian (2007-2008).xlsx")


#### APARCH SNORM ####
returns_hat_df = NULL
statistics = NULL

for (k in 2:18){
  
  returns_df <- returns %>% 
    spread(name, value)
  
  returns_df %>%
    tibble() %>% 
    dplyr::select("date", k) %>%
    purrr::set_names('date', 'stock') %>% 
    inner_join(factors, by='date') %>% 
    na.omit() %>% 
    data.frame() %>% 
    column_to_rownames('date')-> stock_plus_factors
  
  # DCC GARCH Fit
  xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                     variance.model = list(garchOrder = c(1,1), model = 'apARCH'),
                     distribution.model = 'snorm')
  uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
  spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
  cl = makePSOCKcluster(ncol(stock_plus_factors))
  multf = multifit(uspec, stock_plus_factors, cluster = cl)
  
  fit_yx = NULL
  n_try = 0
  
  dcc_function = function (){
    dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
    return(dcc_model) 
  }
  
  possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
  
  while (is.null(fit_yx)){
    fit_yx = possibly_dcc_function()
    n_try = n_try + 1
    print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
    if (n_try >= 20){
      break
    }
  }
  
  if (is.null(fit_yx)){
    print(paste0(colnames(returns_df)[k], " did not work."))
  } else {
    
    for (i in 1:nrow(stock_plus_factors)){
      fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
    }
    
    # DEFINIÇÃO DOS BETAS
    for (i in 1:nrow(stock_plus_factors)){
      new_date = (stock_plus_factors %>% rownames())[i]
      
      x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
      y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
      parameters <- solve(x, y)
      y.bar <- mean(stock_plus_factors[1] %>% unlist())
      x.bar <- colMeans(stock_plus_factors[2:6])
      intercept <- y.bar - x.bar %*% parameters
      
      new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
        purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
      rownames(new_betas) <- new_date
      
      if (i == 1){
        all_betas <- new_betas
      } else {
        all_betas <- rbind(all_betas, new_betas)
      }
      setTxtProgressBar(txtProgressBar(
        min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
    }
    
    # RETORNOS AJUSTADOS E TESTES
    test <- (stock_plus_factors %>% 
               purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
      merge(all_betas %>%
              purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
            by = 0) %>%
      rename(date = Row.names) %>% 
      mutate(stock_rt_hat = intercept_bt +
               MKT_rt * MKT_bt +
               SMB_rt * SMB_bt +
               HML_rt * HML_bt +
               RMW_rt * RMW_bt +
               CMA_rt * CMA_bt) %>%
      mutate(error = stock_rt_hat - stock_rt) %>% 
      mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
      mutate(mse = (stock_rt_hat - stock_rt)^2)
    
    if (is.null(returns_hat_df)){
      returns_hat_df <- test %>% 
        dplyr::select(date, stock_rt_hat) %>% 
        mutate(id = colnames(returns_df)[k])
    } else {
      returns_hat_df <- returns_hat_df %>% 
        rbind(test %>% 
                dplyr::select(date, stock_rt_hat) %>% 
                mutate(id = colnames(returns_df)[k]))
    }
    
    test$mae %>% mean() -> mae
    test$mse %>% mean() -> mse
    (shapiro.test(test$error))$p.value -> shapiro_test
    (jarque.bera.test(test$error))$p.value -> jarque_bera_test
    (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
    (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
    (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
    (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
    skewness(test$error)[1] -> skewness
    kurtosis(test$error)[1] -> kurtosis
    
    new_statistics <- tibble(name = colnames(returns_df)[k],
                             MSE = mse, MAE = mae,
                             `Shapiro Wild - P-Valor` = shapiro_test,
                             `Jarque Bera - P-Valor` = jarque_bera_test,
                             `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
                             `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
                             `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
                             `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
                             `Assimetria` = skewness,
                             `Curtose` = kurtosis) %>% 
      gather(id, value, -name)
    
    if (is.null(statistics)){
      statistics <- new_statistics
    } else {
      statistics <- rbind(statistics, new_statistics)
    }
    
    # GRÁFICOS
    # BETAS
    graphic_betas <- all_betas %>%
      rownames_to_column('date') %>% 
      mutate(date = date %>% ymd()) %>% 
      rename(Intercept = intercept) %>% 
      gather(id, value, -date) %>%
      mutate(type = 'normal') %>% 
      ggplot(aes(x=date, y=value)) +
      geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
      geom_line(aes(group=id), color="black", size=0.85) +
      theme_minimal() +
      facet_wrap(~ id) +
      labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - APARCH Skew Gaussian Distribution (2007-2008)'),
           x = '',
           y = '') +
      theme(legend.title = element_blank(),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',',
        decimal.mark = '.')) +
      scale_x_date(labels = scales::date_format("%Y-%m"))
    
    # ERRO
    graphic_erro <- test %>% 
      dplyr::select(date, error) %>%
      purrr::set_names('date', 'error') %>%
      mutate(date = date %>% ymd()) %>% 
      ggplot(aes(x=date, y=error)) +
      geom_line(color = 'black') +
      labs(title = 'Model Residuals',
           x = '',
           y = '') +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',', decimal.mark = '.')) +
      theme_minimal()
    
    # ACF E PACF
    conf.level <- 0.95
    ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
    
    bacf <- acf(test$error)
    df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Residuals')
    
    bpacf <- pacf(test$error)
    df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Residuals')
    
    bacf <- acf(test$error^2)
    df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Squared Residuals')
    
    bpacf <- pacf(test$error^2)
    df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Squared Residuals')
    
    acf_and_pacf <- plyr::join_all(list(
      df_acf, df_pacf, df_acf_sq, df_pacf_sq),
      by = 'lag', type = 'full') %>%
      tidyr::gather(id, value, -lag)
    
    graphic_acf_and_pacf <- acf_and_pacf %>% 
      ggplot(aes(x=lag, y=value)) +
      geom_bar(fill='black', stat = 'identity') +
      facet_wrap(~ id, nrow=2) +
      geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
      geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
                                                        big.mark = ',',
                                                        decimal.mark = '.')) +
      labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
      theme_minimal()
    
    # JUNÇÃO DOS GRÁFICOS
    plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
              nrow = 3, rel_heights = c(2, 1, 2))
    ggsave(paste0(colnames(returns_df)[k], " - APARCH Distribution Skew Gaussian (2007-2008).png"), device = 'png',
           width = 11, height = 16.5)
    
    setTxtProgressBar(txtProgressBar(
      min = 2, max = 18, initial = 0, style = 3), k)
  }
}

returns_hat_df %>% writexl::write_xlsx("Returns Hat - APARCH Dist Skew Gaussian (2007-2008).xlsx")
statistics %>% writexl::write_xlsx("Statistics - APARCH Dist Skew Gaussian (2007-2008).xlsx")

# #### APARCH STD ####
# returns_hat_df = NULL
# statistics = NULL

# for (k in 2:18){
#   
#   returns_df <- returns %>% 
#     spread(name, value)
#   
#   returns_df %>%
#     tibble() %>% 
#     dplyr::select("date", k) %>%
#     purrr::set_names('date', 'stock') %>% 
#     inner_join(factors, by='date') %>% 
#     na.omit() %>% 
#     data.frame() %>% 
#     column_to_rownames('date')-> stock_plus_factors
#   
#   # DCC GARCH Fit
#   xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
#                      variance.model = list(garchOrder = c(1,1), model = 'apARCH'),
#                      distribution.model = 'std')
#   uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
#   spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
#   cl = makePSOCKcluster(ncol(stock_plus_factors))
#   multf = multifit(uspec, stock_plus_factors, cluster = cl)
#   
#   fit_yx = NULL
#   n_try = 0
#   
#   dcc_function = function (){
#     dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
#     return(dcc_model) 
#   }
#   
#   possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
#   
#   while (is.null(fit_yx)){
#     fit_yx = possibly_dcc_function()
#     n_try = n_try + 1
#     print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
#     if (n_try >= 20){
#       break
#     }
#   }
#   
#   if (is.null(fit_yx)){
#     print(paste0(colnames(returns_df)[k], " did not work."))
#   } else {
#     
#     for (i in 1:nrow(stock_plus_factors)){
#       fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
#     }
#     
#     # DEFINIÇÃO DOS BETAS
#     for (i in 1:nrow(stock_plus_factors)){
#       new_date = (stock_plus_factors %>% rownames())[i]
#       
#       x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
#       y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
#       parameters <- solve(x, y)
#       y.bar <- mean(stock_plus_factors[1] %>% unlist())
#       x.bar <- colMeans(stock_plus_factors[2:6])
#       intercept <- y.bar - x.bar %*% parameters
#       
#       new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
#         purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
#       rownames(new_betas) <- new_date
#       
#       if (i == 1){
#         all_betas <- new_betas
#       } else {
#         all_betas <- rbind(all_betas, new_betas)
#       }
#       setTxtProgressBar(txtProgressBar(
#         min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
#     }
#     
#     # RETORNOS AJUSTADOS E TESTES
#     test <- (stock_plus_factors %>% 
#                purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
#       merge(all_betas %>%
#               purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
#             by = 0) %>%
#       rename(date = Row.names) %>% 
#       mutate(stock_rt_hat = intercept_bt +
#                MKT_rt * MKT_bt +
#                SMB_rt * SMB_bt +
#                HML_rt * HML_bt +
#                RMW_rt * RMW_bt +
#                CMA_rt * CMA_bt) %>%
#       mutate(error = stock_rt_hat - stock_rt) %>% 
#       mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
#       mutate(mse = (stock_rt_hat - stock_rt)^2)
#     
#     if (is.null(returns_hat_df)){
#       returns_hat_df <- test %>% 
#         dplyr::select(date, stock_rt_hat) %>% 
#         mutate(id = colnames(returns_df)[k])
#     } else {
#       returns_hat_df <- returns_hat_df %>% 
#         rbind(test %>% 
#                 dplyr::select(date, stock_rt_hat) %>% 
#                 mutate(id = colnames(returns_df)[k]))
#     }
#     
#     test$mae %>% mean() -> mae
#     test$mse %>% mean() -> mse
#     (shapiro.test(test$error))$p.value -> shapiro_test
#     (jarque.bera.test(test$error))$p.value -> jarque_bera_test
#     (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
#     (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
#     (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
#     (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
#     skewness(test$error)[1] -> skewness
#     kurtosis(test$error)[1] -> kurtosis
#     
#     new_statistics <- tibble(name = colnames(returns_df)[k],
#                              MSE = mse, MAE = mae,
#                              `Shapiro Wild - P-Valor` = shapiro_test,
#                              `Jarque Bera - P-Valor` = jarque_bera_test,
#                              `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
#                              `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
#                              `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
#                              `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
#                              `Assimetria` = skewness,
#                              `Curtose` = kurtosis) %>% 
#       gather(id, value, -name)
#     
#     if (is.null(statistics)){
#       statistics <- new_statistics
#     } else {
#       statistics <- rbind(statistics, new_statistics)
#     }
#     
#     # GRÁFICOS
#     # BETAS
#     graphic_betas <- all_betas %>%
#       rownames_to_column('date') %>% 
#       mutate(date = date %>% ymd()) %>% 
#       rename(Intercept = intercept) %>% 
#       gather(id, value, -date) %>%
#       mutate(type = 'normal') %>% 
#       ggplot(aes(x=date, y=value)) +
#       geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
#       geom_line(aes(group=id), color="black", size=0.85) +
#       theme_minimal() +
#       facet_wrap(~ id) +
#       labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - APARCH T-Student Distribution (2007-2008)'),
#            x = '',
#            y = '') +
#       theme(legend.title = element_blank(),
#             legend.position = 'right',
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       scale_y_continuous(labels = scales::number_format(
#         accuracy = 0.01,
#         big.mark = ',',
#         decimal.mark = '.')) +
#       scale_x_date(labels = scales::date_format("%Y-%m"))
#     
#     # ERRO
#     graphic_erro <- test %>% 
#       dplyr::select(date, error) %>%
#       purrr::set_names('date', 'error') %>%
#       mutate(date = date %>% ymd()) %>% 
#       ggplot(aes(x=date, y=error)) +
#       geom_line(color = 'black') +
#       labs(title = 'Model Residuals',
#            x = '',
#            y = '') +
#       scale_y_continuous(labels = scales::number_format(
#         accuracy = 0.01,
#         big.mark = ',', decimal.mark = '.')) +
#       theme_minimal()
#     
#     # ACF E PACF
#     conf.level <- 0.95
#     ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
#     
#     bacf <- acf(test$error)
#     df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
#       purrr::set_names('lag', 'ACF Residuals')
#     
#     bpacf <- pacf(test$error)
#     df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
#       purrr::set_names('lag', 'PACF Residuals')
#     
#     bacf <- acf(test$error^2)
#     df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
#       purrr::set_names('lag', 'ACF Squared Residuals')
#     
#     bpacf <- pacf(test$error^2)
#     df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
#       purrr::set_names('lag', 'PACF Squared Residuals')
#     
#     acf_and_pacf <- plyr::join_all(list(
#       df_acf, df_pacf, df_acf_sq, df_pacf_sq),
#       by = 'lag', type = 'full') %>%
#       tidyr::gather(id, value, -lag)
#     
#     graphic_acf_and_pacf <- acf_and_pacf %>% 
#       ggplot(aes(x=lag, y=value)) +
#       geom_bar(fill='black', stat = 'identity') +
#       facet_wrap(~ id, nrow=2) +
#       geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
#       geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
#       scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
#                                                         big.mark = ',',
#                                                         decimal.mark = '.')) +
#       labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
#       theme_minimal()
#     
#     # JUNÇÃO DOS GRÁFICOS
#     plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
#               nrow = 3, rel_heights = c(2, 1, 2))
#     ggsave(paste0(colnames(returns_df)[k], " - APARCH Distribution T-Student (2007-2008).png"), device = 'png',
#            width = 11, height = 16.5)
#     
#     setTxtProgressBar(txtProgressBar(
#       min = 2, max = 18, initial = 0, style = 3), k)
#   }
# }
# 
# returns_hat_df %>% writexl::write_xlsx("Returns Hat - APARCH Dist T-Student (2007-2008).xlsx")
# statistics %>% writexl::write_xlsx("Statistics - APARCH Dist T-Student (2007-2008).xlsx")
# 
# #### APARCH SSTD ####
# returns_hat_df = NULL
# statistics = NULL

# for (k in 2:18){
#   
#   returns_df <- returns %>% 
#     spread(name, value)
#   
#   returns_df %>%
#     tibble() %>% 
#     dplyr::select("date", k) %>%
#     purrr::set_names('date', 'stock') %>% 
#     inner_join(factors, by='date') %>% 
#     na.omit() %>% 
#     data.frame() %>% 
#     column_to_rownames('date')-> stock_plus_factors
#   
#   # DCC GARCH Fit
#   xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
#                      variance.model = list(garchOrder = c(1,1), model = 'apARCH'),
#                      distribution.model = 'sstd')
#   uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
#   spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
#   cl = makePSOCKcluster(ncol(stock_plus_factors))
#   multf = multifit(uspec, stock_plus_factors, cluster = cl)
#   
#   fit_yx = NULL
#   n_try = 0
#   
#   dcc_function = function (){
#     dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
#     return(dcc_model) 
#   }
#   
#   possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
#   
#   while (is.null(fit_yx)){
#     fit_yx = possibly_dcc_function()
#     n_try = n_try + 1
#     print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
#     if (n_try >= 20){
#       break
#     }
#   }
#   
#   if (is.null(fit_yx)){
#     print(paste0(colnames(returns_df)[k], " did not work."))
#   } else {
#     
#     for (i in 1:nrow(stock_plus_factors)){
#       fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
#     }
#     
#     # DEFINIÇÃO DOS BETAS
#     for (i in 1:nrow(stock_plus_factors)){
#       new_date = (stock_plus_factors %>% rownames())[i]
#       
#       x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
#       y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
#       parameters <- solve(x, y)
#       y.bar <- mean(stock_plus_factors[1] %>% unlist())
#       x.bar <- colMeans(stock_plus_factors[2:6])
#       intercept <- y.bar - x.bar %*% parameters
#       
#       new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
#         purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
#       rownames(new_betas) <- new_date
#       
#       if (i == 1){
#         all_betas <- new_betas
#       } else {
#         all_betas <- rbind(all_betas, new_betas)
#       }
#       setTxtProgressBar(txtProgressBar(
#         min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
#     }
#     
#     # RETORNOS AJUSTADOS E TESTES
#     test <- (stock_plus_factors %>% 
#                purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
#       merge(all_betas %>%
#               purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
#             by = 0) %>%
#       rename(date = Row.names) %>% 
#       mutate(stock_rt_hat = intercept_bt +
#                MKT_rt * MKT_bt +
#                SMB_rt * SMB_bt +
#                HML_rt * HML_bt +
#                RMW_rt * RMW_bt +
#                CMA_rt * CMA_bt) %>%
#       mutate(error = stock_rt_hat - stock_rt) %>% 
#       mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
#       mutate(mse = (stock_rt_hat - stock_rt)^2)
#     
#     if (is.null(returns_hat_df)){
#       returns_hat_df <- test %>% 
#         dplyr::select(date, stock_rt_hat) %>% 
#         mutate(id = colnames(returns_df)[k])
#     } else {
#       returns_hat_df <- returns_hat_df %>% 
#         rbind(test %>% 
#                 dplyr::select(date, stock_rt_hat) %>% 
#                 mutate(id = colnames(returns_df)[k]))
#     }
#     
#     test$mae %>% mean() -> mae
#     test$mse %>% mean() -> mse
#     (shapiro.test(test$error))$p.value -> shapiro_test
#     (jarque.bera.test(test$error))$p.value -> jarque_bera_test
#     (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
#     (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
#     (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
#     (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
#     skewness(test$error)[1] -> skewness
#     kurtosis(test$error)[1] -> kurtosis
#     
#     new_statistics <- tibble(name = colnames(returns_df)[k],
#                              MSE = mse, MAE = mae,
#                              `Shapiro Wild - P-Valor` = shapiro_test,
#                              `Jarque Bera - P-Valor` = jarque_bera_test,
#                              `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
#                              `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
#                              `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
#                              `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
#                              `Assimetria` = skewness,
#                              `Curtose` = kurtosis) %>% 
#       gather(id, value, -name)
#     
#     if (is.null(statistics)){
#       statistics <- new_statistics
#     } else {
#       statistics <- rbind(statistics, new_statistics)
#     }
#     
#     # GRÁFICOS
#     # BETAS
#     graphic_betas <- all_betas %>%
#       rownames_to_column('date') %>% 
#       mutate(date = date %>% ymd()) %>% 
#       rename(Intercept = intercept) %>% 
#       gather(id, value, -date) %>%
#       mutate(type = 'normal') %>% 
#       ggplot(aes(x=date, y=value)) +
#       geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
#       geom_line(aes(group=id), color="black", size=0.85) +
#       theme_minimal() +
#       facet_wrap(~ id) +
#       labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - APARCH Skew T-Student Distribution (2007-2008)'),
#            x = '',
#            y = '') +
#       theme(legend.title = element_blank(),
#             legend.position = 'right',
#             axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
#       scale_y_continuous(labels = scales::number_format(
#         accuracy = 0.01,
#         big.mark = ',',
#         decimal.mark = '.')) +
#       scale_x_date(labels = scales::date_format("%Y-%m"))
#     
#     # ERRO
#     graphic_erro <- test %>% 
#       dplyr::select(date, error) %>%
#       purrr::set_names('date', 'error') %>%
#       mutate(date = date %>% ymd()) %>% 
#       ggplot(aes(x=date, y=error)) +
#       geom_line(color = 'black') +
#       labs(title = 'Model Residuals',
#            x = '',
#            y = '') +
#       scale_y_continuous(labels = scales::number_format(
#         accuracy = 0.01,
#         big.mark = ',', decimal.mark = '.')) +
#       theme_minimal()
#     
#     # ACF E PACF
#     conf.level <- 0.95
#     ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
#     
#     bacf <- acf(test$error)
#     df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
#       purrr::set_names('lag', 'ACF Residuals')
#     
#     bpacf <- pacf(test$error)
#     df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
#       purrr::set_names('lag', 'PACF Residuals')
#     
#     bacf <- acf(test$error^2)
#     df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
#       purrr::set_names('lag', 'ACF Squared Residuals')
#     
#     bpacf <- pacf(test$error^2)
#     df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
#       purrr::set_names('lag', 'PACF Squared Residuals')
#     
#     acf_and_pacf <- plyr::join_all(list(
#       df_acf, df_pacf, df_acf_sq, df_pacf_sq),
#       by = 'lag', type = 'full') %>%
#       tidyr::gather(id, value, -lag)
#     
#     graphic_acf_and_pacf <- acf_and_pacf %>% 
#       ggplot(aes(x=lag, y=value)) +
#       geom_bar(fill='black', stat = 'identity') +
#       facet_wrap(~ id, nrow=2) +
#       geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
#       geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
#       scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
#                                                         big.mark = ',',
#                                                         decimal.mark = '.')) +
#       labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
#       theme_minimal()
#     
#     # JUNÇÃO DOS GRÁFICOS
#     plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
#               nrow = 3, rel_heights = c(2, 1, 2))
#     ggsave(paste0(colnames(returns_df)[k], " - APARCH Distribution Skew T-Student (2007-2008).png"), device = 'png',
#            width = 11, height = 16.5)
#     
#     setTxtProgressBar(txtProgressBar(
#       min = 2, max = 18, initial = 0, style = 3), k)
#   }
# }
# 
# returns_hat_df %>% writexl::write_xlsx("Returns Hat - APARCH Dist Skew T-Student (2007-2008).xlsx")
# statistics %>% writexl::write_xlsx("Statistics - APARCH Dist Skew T-Student (2007-2008).xlsx")

#### APARCH GED ####
returns_hat_df = NULL
statistics = NULL

for (k in 2:18){
  
  returns_df <- returns %>% 
    spread(name, value)
  
  returns_df %>%
    tibble() %>% 
    dplyr::select("date", k) %>%
    purrr::set_names('date', 'stock') %>% 
    inner_join(factors, by='date') %>% 
    na.omit() %>% 
    data.frame() %>% 
    column_to_rownames('date')-> stock_plus_factors
  
  # DCC GARCH Fit
  xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                     variance.model = list(garchOrder = c(1,1), model = 'apARCH'),
                     distribution.model = 'ged')
  uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
  spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
  cl = makePSOCKcluster(ncol(stock_plus_factors))
  multf = multifit(uspec, stock_plus_factors, cluster = cl)
  
  fit_yx = NULL
  n_try = 0
  
  dcc_function = function (){
    dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
    return(dcc_model) 
  }
  
  possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
  
  while (is.null(fit_yx)){
    fit_yx = possibly_dcc_function()
    n_try = n_try + 1
    print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
    if (n_try >= 20){
      break
    }
  }
  
  if (is.null(fit_yx)){
    print(paste0(colnames(returns_df)[k], " did not work."))
  } else {
    
    for (i in 1:nrow(stock_plus_factors)){
      fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
    }
    
    # DEFINIÇÃO DOS BETAS
    for (i in 1:nrow(stock_plus_factors)){
      new_date = (stock_plus_factors %>% rownames())[i]
      
      x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
      y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
      parameters <- solve(x, y)
      y.bar <- mean(stock_plus_factors[1] %>% unlist())
      x.bar <- colMeans(stock_plus_factors[2:6])
      intercept <- y.bar - x.bar %*% parameters
      
      new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
        purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
      rownames(new_betas) <- new_date
      
      if (i == 1){
        all_betas <- new_betas
      } else {
        all_betas <- rbind(all_betas, new_betas)
      }
      setTxtProgressBar(txtProgressBar(
        min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
    }
    
    # RETORNOS AJUSTADOS E TESTES
    test <- (stock_plus_factors %>% 
               purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
      merge(all_betas %>%
              purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
            by = 0) %>%
      rename(date = Row.names) %>% 
      mutate(stock_rt_hat = intercept_bt +
               MKT_rt * MKT_bt +
               SMB_rt * SMB_bt +
               HML_rt * HML_bt +
               RMW_rt * RMW_bt +
               CMA_rt * CMA_bt) %>%
      mutate(error = stock_rt_hat - stock_rt) %>% 
      mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
      mutate(mse = (stock_rt_hat - stock_rt)^2)
    
    if (is.null(returns_hat_df)){
      returns_hat_df <- test %>% 
        dplyr::select(date, stock_rt_hat) %>% 
        mutate(id = colnames(returns_df)[k])
    } else {
      returns_hat_df <- returns_hat_df %>% 
        rbind(test %>% 
                dplyr::select(date, stock_rt_hat) %>% 
                mutate(id = colnames(returns_df)[k]))
    }
    
    test$mae %>% mean() -> mae
    test$mse %>% mean() -> mse
    (shapiro.test(test$error))$p.value -> shapiro_test
    (jarque.bera.test(test$error))$p.value -> jarque_bera_test
    (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
    (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
    (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
    (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
    skewness(test$error)[1] -> skewness
    kurtosis(test$error)[1] -> kurtosis
    
    new_statistics <- tibble(name = colnames(returns_df)[k],
                             MSE = mse, MAE = mae,
                             `Shapiro Wild - P-Valor` = shapiro_test,
                             `Jarque Bera - P-Valor` = jarque_bera_test,
                             `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
                             `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
                             `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
                             `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
                             `Assimetria` = skewness,
                             `Curtose` = kurtosis) %>% 
      gather(id, value, -name)
    
    if (is.null(statistics)){
      statistics <- new_statistics
    } else {
      statistics <- rbind(statistics, new_statistics)
    }
    
    # GRÁFICOS
    # BETAS
    graphic_betas <- all_betas %>%
      rownames_to_column('date') %>% 
      mutate(date = date %>% ymd()) %>% 
      rename(Intercept = intercept) %>% 
      gather(id, value, -date) %>%
      mutate(type = 'normal') %>% 
      ggplot(aes(x=date, y=value)) +
      geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
      geom_line(aes(group=id), color="black", size=0.85) +
      theme_minimal() +
      facet_wrap(~ id) +
      labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - APARCH GED Distribution (2007-2008)'),
           x = '',
           y = '') +
      theme(legend.title = element_blank(),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',',
        decimal.mark = '.')) +
      scale_x_date(labels = scales::date_format("%Y-%m"))
    
    # ERRO
    graphic_erro <- test %>% 
      dplyr::select(date, error) %>%
      purrr::set_names('date', 'error') %>%
      mutate(date = date %>% ymd()) %>% 
      ggplot(aes(x=date, y=error)) +
      geom_line(color = 'black') +
      labs(title = 'Model Residuals',
           x = '',
           y = '') +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',', decimal.mark = '.')) +
      theme_minimal()
    
    # ACF E PACF
    conf.level <- 0.95
    ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
    
    bacf <- acf(test$error)
    df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Residuals')
    
    bpacf <- pacf(test$error)
    df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Residuals')
    
    bacf <- acf(test$error^2)
    df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Squared Residuals')
    
    bpacf <- pacf(test$error^2)
    df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Squared Residuals')
    
    acf_and_pacf <- plyr::join_all(list(
      df_acf, df_pacf, df_acf_sq, df_pacf_sq),
      by = 'lag', type = 'full') %>%
      tidyr::gather(id, value, -lag)
    
    graphic_acf_and_pacf <- acf_and_pacf %>% 
      ggplot(aes(x=lag, y=value)) +
      geom_bar(fill='black', stat = 'identity') +
      facet_wrap(~ id, nrow=2) +
      geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
      geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
                                                        big.mark = ',',
                                                        decimal.mark = '.')) +
      labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
      theme_minimal()
    
    # JUNÇÃO DOS GRÁFICOS
    plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
              nrow = 3, rel_heights = c(2, 1, 2))
    ggsave(paste0(colnames(returns_df)[k], " - APARCH Distribution GED (2007-2008).png"), device = 'png',
           width = 11, height = 16.5)
    
    setTxtProgressBar(txtProgressBar(
      min = 2, max = 18, initial = 0, style = 3), k)
  }
}

returns_hat_df %>% writexl::write_xlsx("Returns Hat - APARCH Dist GED (2007-2008).xlsx")
statistics %>% writexl::write_xlsx("Statistics - APARCH Dist GED (2007-2008).xlsx")

#### APARCH SGED ####
returns_hat_df = NULL
statistics = NULL

for (k in 2:18){
  
  returns_df <- returns %>% 
    spread(name, value)
  
  returns_df %>%
    tibble() %>% 
    dplyr::select("date", k) %>%
    purrr::set_names('date', 'stock') %>% 
    inner_join(factors, by='date') %>% 
    na.omit() %>% 
    data.frame() %>% 
    column_to_rownames('date')-> stock_plus_factors
  
  # DCC GARCH Fit
  xspec = ugarchspec(mean.model = list(armaOrder = c(0, 0)),
                     variance.model = list(garchOrder = c(1,1), model = 'apARCH'),
                     distribution.model = 'sged')
  uspec = multispec(replicate(ncol(stock_plus_factors), xspec))
  spec1 = dccspec(uspec = uspec, dccOrder = c(1, 1), distribution = 'mvnorm')
  cl = makePSOCKcluster(ncol(stock_plus_factors))
  multf = multifit(uspec, stock_plus_factors, cluster = cl)
  
  fit_yx = NULL
  n_try = 0
  
  dcc_function = function (){
    dcc_model = dccfit(spec1, data = stock_plus_factors, fit.control = list(eval.se = TRUE), fit = multf, cluster = cl)
    return(dcc_model) 
  }
  
  possibly_dcc_function = purrr::possibly(dcc_function, otherwise = NULL)
  
  while (is.null(fit_yx)){
    fit_yx = possibly_dcc_function()
    n_try = n_try + 1
    print(paste0(colnames(returns_df)[k], ' - Número de tentativas: ', n_try))
    if (n_try >= 20){
      break
    }
  }
  
  if (is.null(fit_yx)){
    print(paste0(colnames(returns_df)[k], " did not work."))
  } else {
    
    for (i in 1:nrow(stock_plus_factors)){
      fit_yx@mfit[["Q"]][[i]] = (fit_yx@mfit["H"] %>% data.frame())[(1+(i-1)*6):(6+(i-1)*6)] %>% as.matrix()
    }
    
    # DEFINIÇÃO DOS BETAS
    for (i in 1:nrow(stock_plus_factors)){
      new_date = (stock_plus_factors %>% rownames())[i]
      
      x <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[2:6, 2:6]
      y <- (fit_yx@mfit[["Q"]][[i]] %>% as.matrix())[1, 2:6]
      parameters <- solve(x, y)
      y.bar <- mean(stock_plus_factors[1] %>% unlist())
      x.bar <- colMeans(stock_plus_factors[2:6])
      intercept <- y.bar - x.bar %*% parameters
      
      new_betas <- data.frame(tibble(append(intercept, parameters)) %>% t()) %>% 
        purrr::set_names('intercept', (stock_plus_factors %>% dplyr::select(-stock)) %>% colnames())
      rownames(new_betas) <- new_date
      
      if (i == 1){
        all_betas <- new_betas
      } else {
        all_betas <- rbind(all_betas, new_betas)
      }
      setTxtProgressBar(txtProgressBar(
        min = 1, max = nrow(stock_plus_factors), initial = 0, style = 3), i)
    }
    
    # RETORNOS AJUSTADOS E TESTES
    test <- (stock_plus_factors %>% 
               purrr::set_names("stock_rt", 'MKT_rt', 'SMB_rt', 'HML_rt', 'RMW_rt', 'CMA_rt')) %>% 
      merge(all_betas %>%
              purrr::set_names("intercept_bt", 'MKT_bt', 'SMB_bt', 'HML_bt', 'RMW_bt', 'CMA_bt'),
            by = 0) %>%
      rename(date = Row.names) %>% 
      mutate(stock_rt_hat = intercept_bt +
               MKT_rt * MKT_bt +
               SMB_rt * SMB_bt +
               HML_rt * HML_bt +
               RMW_rt * RMW_bt +
               CMA_rt * CMA_bt) %>%
      mutate(error = stock_rt_hat - stock_rt) %>% 
      mutate(mae = abs(stock_rt_hat - stock_rt)) %>% 
      mutate(mse = (stock_rt_hat - stock_rt)^2)
    
    if (is.null(returns_hat_df)){
      returns_hat_df <- test %>% 
        dplyr::select(date, stock_rt_hat) %>% 
        mutate(id = colnames(returns_df)[k])
    } else {
      returns_hat_df <- returns_hat_df %>% 
        rbind(test %>% 
                dplyr::select(date, stock_rt_hat) %>% 
                mutate(id = colnames(returns_df)[k]))
    }
    
    test$mae %>% mean() -> mae
    test$mse %>% mean() -> mse
    (shapiro.test(test$error))$p.value -> shapiro_test
    (jarque.bera.test(test$error))$p.value -> jarque_bera_test
    (Box.test(test$error, lag=30, type='Ljung'))$p.value -> ljung_box_30_test
    (Box.test(test$error, lag=90, type='Ljung'))$p.value -> ljung_box_90_test
    (Box.test(test$error^2, lag=30, type='Ljung'))$p.value -> ljung_box_30_test_squared
    (Box.test(test$error^2, lag=90, type='Ljung'))$p.value -> ljung_box_90_test_squared
    skewness(test$error)[1] -> skewness
    kurtosis(test$error)[1] -> kurtosis
    
    new_statistics <- tibble(name = colnames(returns_df)[k],
                             MSE = mse, MAE = mae,
                             `Shapiro Wild - P-Valor` = shapiro_test,
                             `Jarque Bera - P-Valor` = jarque_bera_test,
                             `Ljung Box 30 Erro - P-Valor` = ljung_box_30_test,
                             `Ljung Box 90 Erro - P-Valor` = ljung_box_90_test,
                             `Ljung Box 30 Erro ao Quadrado - P-Valor` = ljung_box_30_test_squared,
                             `Ljung Box 90 Erro ao Quadrado - P-Valor` = ljung_box_90_test_squared,
                             `Assimetria` = skewness,
                             `Curtose` = kurtosis) %>% 
      gather(id, value, -name)
    
    if (is.null(statistics)){
      statistics <- new_statistics
    } else {
      statistics <- rbind(statistics, new_statistics)
    }
    
    # GRÁFICOS
    # BETAS
    graphic_betas <- all_betas %>%
      rownames_to_column('date') %>% 
      mutate(date = date %>% ymd()) %>% 
      rename(Intercept = intercept) %>% 
      gather(id, value, -date) %>%
      mutate(type = 'normal') %>% 
      ggplot(aes(x=date, y=value)) +
      geom_line(aes(y=0), color='gray25', linetype = 'dashed') +
      geom_line(aes(group=id), color="black", size=0.85) +
      theme_minimal() +
      facet_wrap(~ id) +
      labs(title = colnames(returns_df)[k] %>% paste('- Conditional Dynamic Betas - APARCH Skew GED Distribution (2007-2008)'),
           x = '',
           y = '') +
      theme(legend.title = element_blank(),
            legend.position = 'right',
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',',
        decimal.mark = '.')) +
      scale_x_date(labels = scales::date_format("%Y-%m"))
    
    # ERRO
    graphic_erro <- test %>% 
      dplyr::select(date, error) %>%
      purrr::set_names('date', 'error') %>%
      mutate(date = date %>% ymd()) %>% 
      ggplot(aes(x=date, y=error)) +
      geom_line(color = 'black') +
      labs(title = 'Model Residuals',
           x = '',
           y = '') +
      scale_y_continuous(labels = scales::number_format(
        accuracy = 0.01,
        big.mark = ',', decimal.mark = '.')) +
      theme_minimal()
    
    # ACF E PACF
    conf.level <- 0.95
    ciline <- qnorm((1 - conf.level)/2)/sqrt(length(test$error))
    
    bacf <- acf(test$error)
    df_acf <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Residuals')
    
    bpacf <- pacf(test$error)
    df_pacf <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Residuals')
    
    bacf <- acf(test$error^2)
    df_acf_sq <- with(bacf, data.frame(lag, acf))[-1,] %>% 
      purrr::set_names('lag', 'ACF Squared Residuals')
    
    bpacf <- pacf(test$error^2)
    df_pacf_sq <- with(bpacf, data.frame(lag, acf)) %>% 
      purrr::set_names('lag', 'PACF Squared Residuals')
    
    acf_and_pacf <- plyr::join_all(list(
      df_acf, df_pacf, df_acf_sq, df_pacf_sq),
      by = 'lag', type = 'full') %>%
      tidyr::gather(id, value, -lag)
    
    graphic_acf_and_pacf <- acf_and_pacf %>% 
      ggplot(aes(x=lag, y=value)) +
      geom_bar(fill='black', stat = 'identity') +
      facet_wrap(~ id, nrow=2) +
      geom_line(aes(y=-ciline), color = 'gray25', linetype = 'dashed') +
      geom_line(aes(y=ciline), color = 'gray25', linetype = 'dashed') +
      scale_y_continuous(labels = scales::number_format(accuracy = 0.01,
                                                        big.mark = ',',
                                                        decimal.mark = '.')) +
      labs(x = 'Lag', y = '', caption = 'Confidence Interval: 95%') +
      theme_minimal()
    
    # JUNÇÃO DOS GRÁFICOS
    plot_grid(graphic_betas, graphic_erro, graphic_acf_and_pacf,
              nrow = 3, rel_heights = c(2, 1, 2))
    ggsave(paste0(colnames(returns_df)[k], " - APARCH Distribution Skew GED (2007-2008).png"), device = 'png',
           width = 11, height = 16.5)
    
    setTxtProgressBar(txtProgressBar(
      min = 2, max = 18, initial = 0, style = 3), k)
  }
}

returns_hat_df %>% writexl::write_xlsx("Returns Hat - APARCH Dist Skew GED (2007-2008).xlsx")
statistics %>% writexl::write_xlsx("Statistics - APARCH Dist Skew GED (2007-2008).xlsx")



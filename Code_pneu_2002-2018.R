# 0. loading libraries ----------------------------------------------------
library(tidyverse)
library(fpp3)
library(fastDummies)
library(forecast)
library(corrr)
library(dlnm)
library(splines)
library(gridExtra)
library(skimr)
library(latex2exp)
library(Cairo)
library(showtext)
font_add_google("Nunito", "nu")
showtext_auto(enable = TRUE)
Sys.setlocale("LC_TIME", "English")
ggplot2::theme_set(theme_classic())

# 1. loading data ------------------------------------
## (1) sourcing air pollution and meteorlogical factors, holidays
source("../Source_MFAP_holiday.r")

## (2) adjust holiday effect
holidays %>% 
    select(holiday) %>% 
    unique
# consider all same type of holidays without ChuSeok, Korean New-Year
holidays2 <- holidays %>% 
    mutate(holiday = recode(holiday, 추석 = "Chuseok", 설날 = "Korean New Year", .default = "holiday"))

## (3) pneumothroax
group <- read_csv("./data/2005-2018_pneu_group.csv") %>% 
    tsibble(key = Menstruating) %>% 
    fill_gaps(N = 0)
sexage <- read_csv("./data/2005-2018_pneu_sexage.csv") %>% 
    tsibble(key = c(SEX_TYPE, AGE_GROUP)) %>% 
    fill_gaps(N = 0)

# 2. Doing EDA ------------------------------------------------------------------
group %>% 
    ggplot(aes(x = Date, y = N)) +
    geom_line() +
    labs(
        y = "Daily spontaneous pneumothorax cases in South Korea",
        x = "Date"
    ) +
    facet_wrap(~ Menstruating, ncol = 2, scales = "free")

## boxplot of airpollution
daily_air %>% 
    # select(-PM25) %>% 
    pivot_longer(SO2:PM25, names_to = "Airpollution") %>% 
    ggplot(aes(x = Airpollution, y = value)) +
    geom_boxplot() +
    stat_summary(fun = median, geom = "text", col = "tomato",     # Add text to plot
                 vjust = -3, aes(label = paste("Median = ", round(..y.., digits = 3)))) +
    stat_summary(fun = sd, geom = "text", col = "tomato",     # Add text to plot
                 vjust = 1, aes(label = paste("Sd = ", round(..y.., digits = 3)))) +
    facet_wrap(~ Airpollution, nrow = 2, scales = "free") +
    theme(legend.position = "none",
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
# ggsave("./plot/boxplot_airpollution.jpeg", device = "jpeg", dpi = 100, type = "cairo")
summary(daily_air)
## Join and check correlations

group_all <- group %>% 
    as_tibble() %>% 
    group_by(Date) %>% 
    summarize(N = sum(N)) %>% 
    left_join(daily_air)
group_all %>% 
    pivot_longer(cols = SO2:Humid, names_to = "MF_AP") %>% 
    ggplot(aes(x = value, y = N)) +
    geom_point() +
    geom_smooth(method = lm, formula = y ~ ns(x)) +
    labs(
        y = "Daily spontaneous pneumothorax cases in South Korea",
        x = "Meteorological factors and air pollution"
    ) +
    facet_wrap(~ MF_AP, ncol = 2, scales = "free")
group_all %>% 
    select(-Date) %>%
    correlate(method = "spearman") %>% 
    shave() %>% 
    rplot(x, shape = 20, colors = c("red", "green"), legend = TRUE,
          print_cor = TRUE) +
    theme(
        axis.text.x = element_text(angle = 45, hjust=1)
    )


# 3. Joining all data and preparing for fitting the DLNMs --------------------
## for joining all data
season_pneu <- msts(group %>% 
                        filter(Menstruating == "Female") %>%
                        pull(N), 
                    seasonal.periods = c(7, 365.25), 
                    start = c(2005, 1)) %>% 
    fourier(K = c(1, 3)) %>% 
    as_tibble %>% 
    mutate(Date = seq(ymd("2005-01-01"), ymd("2018-12-31"), by = "day")) %>% 
    select(Date, everything())
group2 <- group %>% 
    as_tibble() %>% 
    rename(Sex = Menstruating) %>% 
    mutate(Sex = ifelse(Sex != "Male", "Female", "Male")) %>% 
    group_by(Date, Sex) %>% 
    summarize(N = sum(N)) %>% 
    left_join(daily_air, by = "Date") %>% # 대기오염원, 기상요인
    left_join(holidays2, by = "Date") %>%  # 공휴일
    left_join(season_pneu, by = "Date") %>%  # 계절성
    mutate(
        days = wday(Date, label = TRUE) # 요일 효과
    ) %>% 
    ungroup() %>% 
    nest(data = !Sex) %>% 
    drop_na() %>% 
    mutate(data2 = map(data, ~rowid_to_column(.x, var = "Time"))) # 추세

## for Tuning maximum lag, df of var dimension, df of lag dimension
fqaic <- function(model) { # Q-AIC FUNCTION
    loglik <- sum(dpois(model$y,model$fitted.values,log=TRUE))
    phi <- summary(model)$dispersion
    qaic <- -2*loglik + 2*summary(model)$df[3]*phi
    return(qaic)
}
dlnm_model <- function(.data, .pollution, .index){
    cb <- crossbasis(.pollution, 
                     lag = grid[[.index, ".lag"]],
                     argvar = list(fun = "ns", df = grid[[.index, "var_df"]]),
                     arglag = list(fun = "ns", df = grid[[.index, "lag_df"]]))
    m <- glm(N ~ cb + Temp + Temp_range + Prec + Humid + Time + holiday +
                 `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
             family = quasipoisson(), data = .data)
    fqaic(m)
    
}
grid <- expand.grid(.lag = 7:21, var_df = 2:5, lag_df = 2:5) %>% 
    as_tibble()

# 4-1. Fitting DLNMs and Visualization: PM10 ----------------------------------------------------------------
# Incidences of spontaneous pnx between Male and Female
tune_pm10_f <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                  ~dlnm_model(.data = group2$data2[[1]], 
                              .pollution = group2$data2[[1]]$PM10, 
                              .index = .x))
    ) %>% 
    arrange(qaic)
tune_pm10_m <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                   ~dlnm_model(.data = group2$data2[[2]], 
                               .pollution = group2$data2[[2]]$PM10, 
                               .index = .x))
    ) %>% 
    arrange(qaic)

### fitting the best DLNM models based on QAIC
cb_pm10_f <- crossbasis(group2$data2[[1]]$PM10, 
                   lag = tune_pm10_f[[1, ".lag"]],
                   argvar = list(fun = "ns", df = tune_pm10_f[[1, "var_df"]]),
                   arglag = list(fun = "ns", df = tune_pm10_f[[1, "lag_df"]]))
cb_pm10_m <- crossbasis(group2$data2[[2]]$PM10, 
                   lag = tune_pm10_m[[1, ".lag"]],
                   argvar = list(fun = "ns", df = tune_pm10_m[[1, "var_df"]]),
                   arglag = list(fun = "ns", df = tune_pm10_m[[1, "lag_df"]]))
mod_pm10_f <- glm(N ~ cb_pm10_f + Temp + Temp_range + Prec + Humid + Time + holiday +
                      `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                  family = quasipoisson(), data = group2$data2[[1]])
mod_pm10_m <- glm(N ~ cb_pm10_m + Temp + Temp_range + Prec + Humid + Time + holiday +
                      `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                  family = quasipoisson(), data = group2$data2[[2]])
pred_pm10_f <- crosspred(cb_pm10_f, mod_pm10_f, 
                         cen = 43.5, by = 10, at = 0:100)
pred_pm10_m <- crosspred(cb_pm10_m, mod_pm10_m, 
                         cen = 43.5, by = 10, at = 0:100)

### Visualization - 3D surface plot
jpeg("./plot/3D_pneu_PM10_female.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_pm10_f, xlab = "PM10", zlab = "RR", main = "Female", 
     theta = 210, phi = 30, lphi = 30,
     border = "gray40", ticktype = "detailed", font.main = 1, family = "nu")
dev.off()
jpeg("./plot/3D_pneu_PM10_male.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_pm10_m, xlab = "PM10", zlab = "RR", main = "Male",
     theta = 210, phi = 30, lphi = 30,
     border = "gray40", ticktype = "detailed", font.main = 1, family = "nu")
dev.off()

### Visualization - Contour plot
jpeg("./plot/Contour_pneu_PM10_female.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_pm10_f, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Female", xlab = TeX("$PM_{10} \\, (\\mu g /m^3)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()
jpeg("./plot/Contour_pneu_PM10_male.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_pm10_m, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Male", xlab = TeX("$PM_{10} \\, (\\mu g /m^3)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()

### Visualization - Overall cumulative association
jpeg("./plot/Overall_pneu_PM10_bysex.jpeg", type = "cairo", res = 100, width = 750, height = 800)
par(mfrow = c(2, 1))
plot(pred_pm10_f, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$PM_{10} \\, (\\mu g /m^3)$"), main = "Female", font.main = 1, family = "nu")
rug(group2$data2[[1]]$PM10)
plot(pred_pm10_m, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$PM_{10} \\, (\\mu g /m^3)$"), main = "Male", font.main = 1, family = "nu")
rug(group2$data2[[2]]$PM10)
par(mfrow = c(1, 1))
dev.off()

# 4-2. Fitting DLNMs and Visualization: SO2 ----------------------------------------------------------------
tune_so2_f <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[1]], 
                                   .pollution = group2$data2[[1]]$SO2, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)
tune_so2_m <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[2]], 
                                   .pollution = group2$data2[[2]]$SO2, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)

### fitting the best DLNM models based on QAIC
cb_so2_f <- crossbasis(group2$data2[[1]]$SO2, 
                        lag = tune_so2_f[[1, ".lag"]],
                        argvar = list(fun = "ns", df = tune_so2_f[[1, "var_df"]]),
                        arglag = list(fun = "ns", df = tune_so2_f[[1, "lag_df"]]))
cb_so2_m <- crossbasis(group2$data2[[2]]$SO2, 
                        lag = tune_so2_m[[1, ".lag"]],
                        argvar = list(fun = "ns", df = tune_so2_m[[1, "var_df"]]),
                        arglag = list(fun = "ns", df = tune_so2_m[[1, "lag_df"]]))
mod_so2_f <- glm(N ~ cb_so2_f + Temp + Temp_range + Prec + Humid + Time + holiday +
                      `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                  family = quasipoisson(), data = group2$data2[[1]])
mod_so2_m <- glm(N ~ cb_so2_m + Temp + Temp_range + Prec + Humid + Time + holiday +
                      `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                  family = quasipoisson(), data = group2$data2[[2]])
summary(mod_so2_f)
summary(mod_so2_m)
pred_so2_f <- crosspred(cb_so2_f, mod_so2_f, 
                         cen = 4.3, by = 1, at = 0:10)
pred_so2_m <- crosspred(cb_so2_m, mod_so2_m, 
                        cen = 4.3, by = 1, at = 0:10)

### Visualization - 3D surface plot
jpeg("./plot/3D_pneu_SO2_female.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_so2_f, xlab = "SO2 (ppb)", zlab = "RR", main = "Female", 
     theta = 210, phi = 30, lphi = 30, ticktype = "detailed",
     border = "gray40", font.main = 1, family = "nu")
dev.off()
jpeg("./plot/3D_pneu_SO2_male.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_so2_m, xlab = "SO2 (ppb)", zlab = "RR", main = "Male",
     theta = 210, phi = 30, lphi = 30,
     border = "gray40", ticktype = "detailed", font.main = 1, family = "nu")
dev.off()

### Visualization - Contour plot
jpeg("./plot/Contour_pneu_SO2_female.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_so2_f, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Female", xlab = TeX("$SO_{2} \\, (ppb)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()
jpeg("./plot/Contour_pneu_SO2_male.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_so2_m, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Male", xlab = TeX("$SO_{2} \\, (ppb)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()

### Visualization - Overall cumulative association
jpeg("./plot/Overall_pneu_SO2_bysex.jpeg", type = "cairo", res = 100, width = 750, height = 800)
par(mfrow = c(2, 1))
plot(pred_so2_f, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$SO_{2} \\, (ppb)$"), main = "Female", font.main = 1, family = "nu")
rug(group2$data2[[1]]$SO2)
plot(pred_so2_m, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$SO_{2} \\, (ppb)$"), main = "Male", font.main = 1, family = "nu")
rug(group2$data2[[2]]$SO2)
par(mfrow = c(1, 1))
dev.off()

# 4-3. Fitting DLNMs and Visualization: CO ---------------------------------------------------------------------
tune_co_f <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[1]], 
                                   .pollution = group2$data2[[1]]$CO, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)
tune_co_m <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[2]], 
                                   .pollution = group2$data2[[2]]$CO, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)

### fitting the best DLNM models based on QAIC
cb_co_f <- crossbasis(group2$data2[[1]]$CO, 
                       lag = tune_co_f[[1, ".lag"]],
                       argvar = list(fun = "ns", df = tune_co_f[[1, "var_df"]]),
                       arglag = list(fun = "ns", df = tune_co_f[[1, "lag_df"]]))
cb_co_m <- crossbasis(group2$data2[[2]]$CO, 
                       lag = tune_co_m[[1, ".lag"]],
                       argvar = list(fun = "ns", df = tune_co_m[[1, "var_df"]]),
                       arglag = list(fun = "ns", df = tune_co_m[[1, "lag_df"]]))
mod_co_f <- glm(N ~ cb_co_f + Temp + Temp_range + Prec + Humid + Time + holiday +
                     `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                 family = quasipoisson(), data = group2$data2[[1]])
mod_co_m <- glm(N ~ cb_co_m + Temp + Temp_range + Prec + Humid + Time + holiday +
                     `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                 family = quasipoisson(), data = group2$data2[[2]])
summary(mod_co_f)
summary(mod_co_m)
pred_co_f <- crosspred(cb_co_f, mod_co_f, 
                        cen = 467.3, by = 10, at = 400:600)
pred_co_m <- crosspred(cb_co_m, mod_co_m, 
                        cen = 467.3, by = 10, at = 400:600)

### Visualization - 3D surface plot
jpeg("./plot/3D_pneu_CO_female.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_co_f, xlab = "CO (ppb)", zlab = "RR", main = "Female", 
     theta = 210, phi = 30, lphi = 30, ticktype = "detailed", shade = 0.5,
     border = NA, font.main = 1, family = "nu")
dev.off()
jpeg("./plot/3D_pneu_CO_male.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_co_m, xlab = "CO (ppb)", zlab = "RR", main = "Male",
     theta = 210, phi = 30, lphi = 30,
     border = NA, ticktype = "detailed", font.main = 1, family = "nu")
dev.off()

### Visualization - Contour plot
jpeg("./plot/Contour_pneu_CO_female.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_co_f, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Female", xlab = TeX("$CO \\, (ppb)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()
jpeg("./plot/Contour_pneu_CO_male.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_co_m, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Male", xlab = TeX("$CO \\, (ppb)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()

### Visualization - Overall cumulative association
jpeg("./plot/Overall_pneu_CO_bysex.jpeg", type = "cairo", res = 100, width = 750, height = 800)
par(mfrow = c(2, 1))
plot(pred_co_f, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$CO \\, (ppb)$"), main = "Female", font.main = 1, family = "nu")
rug(group2$data2[[1]]$CO)
plot(pred_co_m, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$CO \\, (ppb)$"), main = "Male", font.main = 1, family = "nu")
rug(group2$data2[[2]]$CO)
par(mfrow = c(1, 1))
dev.off()



# 4-4. Fitting DLNMs and Visualization: NO2  --------------------------------------------------------------------
tune_no2_f <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[1]], 
                                   .pollution = group2$data2[[1]]$NO2, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)
tune_no2_m <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[2]], 
                                   .pollution = group2$data2[[2]]$NO2, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)

### fitting the best DLNM models based on QAIC
cb_no2_f <- crossbasis(group2$data2[[1]]$NO2, 
                       lag = tune_no2_f[[1, ".lag"]],
                       argvar = list(fun = "ns", df = tune_no2_f[[1, "var_df"]]),
                       arglag = list(fun = "ns", df = tune_no2_f[[1, "lag_df"]]))
cb_no2_m <- crossbasis(group2$data2[[2]]$NO2, 
                       lag = tune_no2_m[[1, ".lag"]],
                       argvar = list(fun = "ns", df = tune_no2_m[[1, "var_df"]]),
                       arglag = list(fun = "ns", df = tune_no2_m[[1, "lag_df"]]))
mod_no2_f <- glm(N ~ cb_no2_f + Temp + Temp_range + Prec + Humid + Time + holiday +
                     `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                 family = quasipoisson(), data = group2$data2[[1]])
mod_no2_m <- glm(N ~ cb_no2_m + Temp + Temp_range + Prec + Humid + Time + holiday +
                     `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                 family = quasipoisson(), data = group2$data2[[2]])
summary(mod_no2_f)
summary(mod_no2_m)
pred_no2_f <- crosspred(cb_no2_f, mod_no2_f, 
                        cen = 19.1, by = 1, at = 0:40)
pred_no2_m <- crosspred(cb_no2_m, mod_no2_m, 
                        cen = 19.1, by = 1, at = 0:40)

### Visualization - 3D surface plot
jpeg("./plot/3D_pneu_NO2_female.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_no2_f, xlab = "NO2 (ppb)", zlab = "RR", main = "Female", 
     theta = 210, phi = 30, lphi = 30, ticktype = "detailed",
     border = "gray40", font.main = 1, family = "nu")
dev.off()
jpeg("./plot/3D_pneu_NO2_male.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_no2_m, xlab = "NO2 (ppb)", zlab = "RR", main = "Male",
     theta = 210, phi = 30, lphi = 30,
     border = "gray40", ticktype = "detailed", font.main = 1, family = "nu")
dev.off()

### Visualization - Contour plot
jpeg("./plot/Contour_pneu_NO2_female.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_no2_f, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Female", xlab = TeX("$NO_{2} \\, (ppb)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()
jpeg("./plot/Contour_pneu_NO2_male.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_no2_m, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Male", xlab = TeX("$NO_{2} \\, (ppb)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()

### Visualization - Overall cumulative association
jpeg("./plot/Overall_pneu_NO2_bysex.jpeg", type = "cairo", res = 100, width = 750, height = 800)
par(mfrow = c(2, 1))
plot(pred_no2_f, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$NO_{2} \\, (ppb)$"), main = "Female", font.main = 1, family = "nu")
rug(group2$data2[[1]]$NO2)
plot(pred_no2_m, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$NO_{2} \\, (ppb)$"), main = "Male", font.main = 1, family = "nu")
rug(group2$data2[[2]]$NO2)
par(mfrow = c(1, 1))
dev.off()


# 4-5. Fitting DLNMs and Visualization: O3 --------------------------------------------------------------------
tune_o3_f <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[1]], 
                                   .pollution = group2$data2[[1]]$O3, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)
tune_o3_m <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[2]], 
                                   .pollution = group2$data2[[2]]$O3, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)

### fitting the best DLNM models based on QAIC
cb_o3_f <- crossbasis(group2$data2[[1]]$O3, 
                       lag = tune_o3_f[[1, ".lag"]],
                       argvar = list(fun = "ns", df = tune_o3_f[[1, "var_df"]]),
                       arglag = list(fun = "ns", df = tune_o3_f[[1, "lag_df"]]))
cb_o3_m <- crossbasis(group2$data2[[2]]$O3, 
                       lag = tune_o3_m[[1, ".lag"]],
                       argvar = list(fun = "ns", df = tune_o3_m[[1, "var_df"]]),
                       arglag = list(fun = "ns", df = tune_o3_m[[1, "lag_df"]]))
mod_o3_f <- glm(N ~ cb_o3_f + Temp + Temp_range + Prec + Humid + Time + holiday +
                     `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                 family = quasipoisson(), data = group2$data2[[1]])
mod_o3_m <- glm(N ~ cb_o3_m + Temp + Temp_range + Prec + Humid + Time + holiday +
                     `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                 family = quasipoisson(), data = group2$data2[[2]])
summary(mod_o3_f)
summary(mod_o3_m)
pred_o3_f <- crosspred(cb_o3_f, mod_o3_f, 
                        cen = 25.4, by = 1, at = 0:40)
pred_o3_m <- crosspred(cb_o3_m, mod_o3_m, 
                        cen = 25.4, by = 1, at = 0:40)

### Visualization - 3D surface plot
jpeg("./plot/3D_pneu_O3_female.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_o3_f, xlab = "O3 (ppb)", zlab = "RR", main = "Female", 
     theta = 210, phi = 30, lphi = 30, ticktype = "detailed",
     border = "gray40", font.main = 1, family = "nu")
dev.off()
jpeg("./plot/3D_pneu_O3_male.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_o3_m, xlab = "O3 (ppb)", zlab = "RR", main = "Male",
     theta = 210, phi = 30, lphi = 30,
     border = "gray40", ticktype = "detailed", font.main = 1, family = "nu")
dev.off()

### Visualization - Contour plot
jpeg("./plot/Contour_pneu_O3_female.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_o3_f, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Female", xlab = TeX("$O_{3} \\, (ppb)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()
jpeg("./plot/Contour_pneu_O3_male.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_o3_m, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Male", xlab = TeX("$O_{3} \\, (ppb)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()

### Visualization - Overall cumulative association
jpeg("./plot/Overall_pneu_O3_bysex.jpeg", type = "cairo", res = 100, width = 750, height = 800)
par(mfrow = c(2, 1))
plot(pred_o3_f, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$O_{3} \\, (ppb)$"), main = "Female", font.main = 1, family = "nu")
rug(group2$data2[[1]]$O3)
plot(pred_o3_m, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$O_{3} \\, (ppb)$"), main = "Male", font.main = 1, family = "nu")
rug(group2$data2[[2]]$O3)
par(mfrow = c(1, 1))
dev.off()

# 4-6. Fitting DLNMs and Visualization: PM25  --------------------------------------------------------------------
tune_pm25_f <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[1]], 
                                   .pollution = group2$data2[[1]]$PM25, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)
tune_pm25_m <- grid %>% 
    mutate(
        qaic = map_dbl(seq(nrow(grid)), 
                       ~dlnm_model(.data = group2$data2[[2]], 
                                   .pollution = group2$data2[[2]]$PM25, 
                                   .index = .x))
    ) %>% 
    arrange(qaic)

### fitting the best DLNM models based on QAIC
cb_pm25_f <- crossbasis(group2$data2[[1]]$PM25, 
                      lag = tune_pm25_f[[1, ".lag"]],
                      argvar = list(fun = "ns", df = tune_pm25_f[[1, "var_df"]]),
                      arglag = list(fun = "ns", df = tune_pm25_f[[1, "lag_df"]]))
cb_pm25_m <- crossbasis(group2$data2[[2]]$PM25, 
                      lag = tune_pm25_m[[1, ".lag"]],
                      argvar = list(fun = "ns", df = tune_pm25_m[[1, "var_df"]]),
                      arglag = list(fun = "ns", df = tune_pm25_m[[1, "lag_df"]]))
mod_pm25_f <- glm(N ~ cb_pm25_f + Temp + Temp_range + Prec + Humid + Time + holiday +
                    `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                family = quasipoisson(), data = group2$data2[[1]])
mod_pm25_m <- glm(N ~ cb_pm25_m + Temp + Temp_range + Prec + Humid + Time + holiday +
                    `S1-7` + `C1-7` + `S1-365` + `C1-365` + `S2-365` + `C2-365` + `S3-365` + `C3-365`,
                family = quasipoisson(), data = group2$data2[[2]])
summary(mod_pm25_f)
summary(mod_pm25_m)
pred_pm25_f <- crosspred(cb_pm25_f, mod_pm25_f, 
                       cen = 21.5, by = 1, at = 15:50)
pred_pm25_m <- crosspred(cb_pm25_m, mod_pm25_m, 
                       cen = 21.5, by = 1, at = 15:50)

### Visualization - 3D surface plot
jpeg("./plot/3D_pneu_PM25_female.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_pm25_f, xlab = "PM25", zlab = "RR", main = "Female", 
     theta = 210, phi = 30, lphi = 30, ticktype = "detailed",
     border = "gray40", font.main = 1, family = "nu")
dev.off()
jpeg("./plot/3D_pneu_PM25_male.jpeg", type = "cairo", res = 100, width = 700, height = 750)
plot(pred_pm25_m, xlab = "PM25", zlab = "RR", main = "Male",
     theta = 210, phi = 30, lphi = 30,
     border = "gray40", ticktype = "detailed", font.main = 1, family = "nu")
dev.off()

### Visualization - Contour plot
jpeg("./plot/Contour_pneu_PM25_female.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_pm25_f, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Female", xlab = TeX("$PM_{25} \\, (\\mu g /m^3)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()
jpeg("./plot/Contour_pneu_PM25_male.jpeg", type = "cairo", res = 100, width = 800, height = 500)
plot(pred_pm25_m, "contour", key.title = title("RR", font.main = 1),
     plot.title = title("Male", xlab = TeX("$PM_{25} \\, (\\mu g /m^3)$"), ylab = "Lag", 
                        font.main = 1, family = "nu"))
dev.off()

### Visualization - Overall cumulative association
jpeg("./plot/Overall_pneu_PM25_bysex.jpeg", type = "cairo", res = 100, width = 750, height = 800)
par(mfrow = c(2, 1))
plot(pred_pm25_f, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$PM_{25} \\, (\\mu g /m^3)$"), main = "Female", font.main = 1, family = "nu")
rug(group2$data2[[1]]$PM25)
plot(pred_pm25_m, col = "tomato", lwd = 2, "overall", ci.arg = list(density = 15, lwd = 2),
     ylab = "RR", xlab = TeX("$PM_{25} \\, (\\mu g /m^3)$"), main = "Male", font.main = 1, family = "nu")
rug(group2$data2[[2]]$PM25)
par(mfrow = c(1, 1))
dev.off()


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
Sys.setlocale("LC_TIME", "English")
ggplot2::theme_set(theme_classic())

# 1. loading data ------------------------------------
## (1) air pollution, meteorological factors
daily_air <- read_rds("../../Data/airpollution and meteorological factors/airpollution and meteorlogical factors.rds") %>% 
    rename(Date = date,
           SO2 = SO2_Avg,
           CO = CO_Avg,
           O3 = O3_Avg,
           NO2 = NO2_Avg,
           PM10 = PM10_Avg,
           PM25 = PM25_Avg,
           `Relative humidity` = Humidity_Avg,
           Temperature = Temperature_Avg,
           `Temperature range` = Temperature_range)

## (2) holidays
holidays <- read_csv("../../Data/holidays in South Korea/korean_holiday_2004-2022_remove duplicated.csv") %>% 
    rename(holiday = dateName,
           Date = date) %>% 
    select(Date, holiday)
holidays %>% 
    select(holiday) %>% 
    unique
# consider all same type of holidays without ChuSeok, Korean New-Year
holidays2 <- holidays %>% 
    mutate(holiday = recode(holiday, 추석 = "Chuseok", 설날 = "Korean New Year", .default = "holiday"))

## (3) pneumothroax
group <- read_csv("./data/2005-2018_repneu_group.csv") %>% 
    filter(year(Date) > 2014) %>% 
    tsibble(key = Menstruating) %>% 
    fill_gaps(N = 0)
sexage <- read_csv("./data/2005-2018_repneu_sexage.csv") %>% 
    tsibble(key = c(SEX_TYPE, AGE_GROUP)) %>% 
    fill_gaps(N = 0)

# 2. Doing EDA ------------------------------------------------------------------
re_group %>%
    ggplot(aes(x = Date, y = N)) +
    geom_line()+
    labs(
        y = "Daily spontaneous pneumothorax relapse cases in South Korea",
        x = "Date"
    ) +
    facet_wrap(~ Menstruating, ncol = 2, scales = "free")
# re_group %>% # episode 처리 필요 -> episode 처리 필요없이 30일이아니라 31일 또는 32일 lag로 잡자: 완료
#     filter(Menstruating == "Female", N > 3) %>% 
#     View

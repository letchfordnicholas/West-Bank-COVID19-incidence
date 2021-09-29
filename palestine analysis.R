rm(list=ls())

# libraries

library(tidyverse)
library(readxl)
library(EpiEstim)
library(incidence)
library(ggplot2)
library(dplyr)
library(magrittr)
library(tidyr)
library(readr)
library(zoo)
library(tseries)
library(lubridate)
library(ggpubr)

setwd("~/Desktop/work_related_stuff/WHO/Palestine/data")

# ageg groups
ages <- c(y0_9="0 - 9",
          y10_19="10 - 19",
          y20_29="20 - 29",
          y30_39="30 - 39",
          y40_49="40 - 49",
          y50_59="50 - 59",
          y60_69="60 - 69",
          y70_79="70 - 79",
          y80="80 and over")

# load WB data
case_data <- read_csv("WB Cases Age Group Jan 31.csv") %>%
  gather(key=age_group,value=incidence,-date) %>%
  mutate(date = as.Date(date,"%d/%m/%Y"),
         age_group = recode_factor(age_group,!!!ages)) %>% 
  group_by(age_group) %>%
  arrange(-desc(date)) 


# --------------------------------------------------------------------------------------------------------
# Rt calculation

# merge Israel and WB datasets
case_data_sum <- case_data %>%
  ungroup() %>%
  group_by(date) %>%
  summarise(incidence = sum(incidence)) %>%
  mutate(smooth_incidence = rollapply(incidence,7,mean,align='right',fill=NA)) %>%
  filter(!is.na(incidence)) 

cum_cas<-case_data_sum %>% 
  dplyr::group_by(date) %>%
  arrange(date) %>% 
  plyr:: mutate(cum_incidence = cumsum(incidence))

# where is 12 cumulative cases reached? 12 cases was noted as a requirement for EpiEstim
s <- which.max(cum_cas$cum_incidence > 12)
cum_cas <- cum_cas[s:nrow(cum_cas),]

# add index
cum_cas <- cum_cas %>%
  rowid_to_column() %>%
  rename(index = "rowid")

# 2 week Rt using epiEstim package
T0 <- nrow(cum_cas)
t_start <- seq(2, T0-13) # starting at 2 as conditional on the past observations
t_end <- t_start + 13

res_parametric_si <- estimate_R(cum_cas$incidence,
                                method = "parametric_si",
                                config = make_config(list(
                                  mean_si = 5.12,
                                  std_si = 4.28,
                                  t_start = t_start,
                                  t_end = t_end
                                )))$R[,c("t_end","Mean(R)","Quantile.0.025(R)","Quantile.0.975(R)")]

names(res_parametric_si)<-c("index","mean_R","lower_R","upper_R")

# merge with incidence data for having one data set for both Rt and EpiCurve side x side
epidata <- merge(cum_cas, res_parametric_si, by = "index", all = TRUE, sort = FALSE) %>%
  select(date, incidence, mean_R, lower_R, upper_R)


# -----------------------------------------------------------------------------------------------------------------------
# Figure 2: Reported incidence disaggregated by 10 year age groups

case_data %>%
  mutate(smoothed_incidence = rollapply(incidence,7,mean,align='right',fill=NA)) %>%
  ggplot(aes(x=date,y=smoothed_incidence,fill=age_group)) +
    geom_bar(position = "stack",stat="identity",width=1) +
    ylab("Daily Cases") +
    theme(axis.text.x = element_text(size=12),
          axis.title.x = element_text(size=16),
          axis.text.y = element_text(size=12),
          axis.title.y = element_text(size=16),
          panel.background = element_rect(fill = "white",
                                          colour = "black",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major = element_line(size = 0.5, linetype = 'solid',
                                            colour = "lightgrey")) +
  guides(fill=guide_legend(title="Age Group"))


# -------------------------------------------------------------------------------------------------------------------------------------
# Test Positivity Ratio, Rt, incidence, mortality and PHSM data formating

# load test positivity ratio data 
positivity_test_data <- read_csv("TPR WB T0131.csv") %>%
  rename(date = "Date",
         positivity = "Positivity") %>%
  mutate(date = as.Date(date,"%d/%m/%Y"),
         positivity = positivity * 100,
         smoothed_positivity = rollapply(positivity,7,mean,align='right',fill=NA),
         smoothed_positivity = ifelse(date == as.Date("2020-03-06") | date == as.Date("2020-03-05"),5.428,smoothed_positivity)) 

# load mortality data
deaths_wb <- read_csv("Deaths WB Jan 31.csv") %>%
  rename(date = "Date",
         deaths = "Deaths") %>%
  mutate(date = as.Date(date,"%d/%m/%Y"),
         smooth_deaths = rollapply(deaths,7,mean,align='right',fill=NA))

# format Rt data
epidata <- epidata %>%
  rename(mean_r = "mean_R",
         lower_r = "lower_R",
         upper_r = "upper_R") %>%
  mutate(date = as.Date(date,"%d/%m/%Y"),
         smooth_upper_r = rollapply(upper_r,7,mean,align='right',fill=NA),
         smooth_lower_r = rollapply(lower_r,7,mean,align='right',fill=NA),
         smooth_mean_r = rollapply(mean_r,7,mean,align='right',fill=NA)) %>%
  select(-incidence)

# merge all data sources
merged_data <- case_data_sum %>%
  full_join(deaths_wb,by="date") %>%
  full_join(positivity_test_data,by="date") %>%
  full_join(epidata,by="date") %>%
  arrange(-desc(date)) %>%
  filter(date >= "2020-03-05")

# Figure 3: Epidemiological indicators and PHSM measures as applied in the West Bank: Rt and mortality 
coef_1 <- 5
bottom_plot <- ggplot(merged_data,aes(x=date)) +
  geom_line(aes(y=smooth_mean_r),colour="forestgreen",size=2) +
  geom_ribbon(aes(ymin=smooth_lower_r,ymax=smooth_upper_r),alpha=0.5,fill="forestgreen") +
  geom_line(aes(y=smooth_deaths/coef_1),colour="red",size=2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Rt",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ .*coef_1, name="Deaths")) +
  theme(axis.text.y.left = element_text(colour="forestgreen",size=16),
        axis.text.y.right = element_text(colour="red",size=16),
        axis.title.y.left = element_text(colour = "forestgreen",size=20),
        axis.title.y.right = element_text(colour = "red",size=20),
        axis.title.x = element_text(size=20),
        axis.text.x = element_text(size=16),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid"),
        panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                        colour = "lightgrey")) +
  scale_x_date(date_labels = "%b %Y")

print(bottom_plot)

# Figure 3: Epidemiological indicators and PHSM measures as applied in the West Bank: Incidence and TPR
coef_2 <- 0.02
middle_plot <- ggplot(merged_data,aes(x=date)) +
  geom_line(aes(y=smooth_incidence),colour="blue",size=2) +
  geom_line(aes(y=smoothed_positivity/coef_2),colour="orange",size=2) + 
  scale_y_continuous(
    # Features of the first axis
    name = "Daily Cases",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~ .*coef_2, name="Test Positivity Rate (%)")) +
  theme(axis.text.y.left = element_text(colour="blue",size=16),
        axis.text.y.right = element_text(colour="orange",size=16),
        axis.title.y.right = element_text(colour = "orange",size=20),
        axis.title.y.left = element_text(colour = "blue",size=20)) +
  xlab("") + theme(axis.text.x = element_blank(),
                   axis.ticks = element_blank(),
                   panel.background = element_rect(fill = "white",
                                                   colour = "black",
                                                   size = 0.5, linetype = "solid"),
                   panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                                     colour = "lightgrey"))

print(middle_plot)

# creating data for PHSM figure
response_data <- data.frame(date = seq(from = as.Date("2020-03-01"),to = as.Date("2021-01-31"),by = "day"),
                                curfew = c(NA),
                                school_closure = c(NA),
                                int_travel_ban = c(NA),
                                mass_gather = c(NA),
                                work_from_home = c(NA),
                                business_closure = c(NA)) %>%
  mutate(curfew = ifelse(date >= as.Date("2020-3-22") & date <= as.Date("2020-5-25"),1,
                         ifelse(date >= as.Date("2020-5-26") & date <= as.Date("2020-9-1"), 0.5,
                                ifelse(date >= as.Date("2020-12-1"),0.5,curfew))),
         school_closure = ifelse(date >= as.Date("2020-3-5") & date <= as.Date("2020-6-1"),1,
                                 ifelse(date >= as.Date("2020-8-1") & date <= as.Date("2020-11-30"), 0.5,
                                        ifelse(date >= as.Date("2021-1-12"),0.5,school_closure))),
         int_travel_ban = ifelse(date >= as.Date("2020-3-9") & date <= as.Date("2020-6-30"),1,
                                 ifelse(date >= as.Date("2020-7-1"),0.5,int_travel_ban)),
         mass_gather = ifelse(date >= as.Date("2020-3-9") & date <= as.Date("2020-6-30"),1,
                                 ifelse(date >= as.Date("2020-7-1") & date <= as.Date("2020-8-31"), 0.5,
                                        ifelse(date >= as.Date("2020-9-1"),1,mass_gather))),
         work_from_home = ifelse(date >= as.Date("2020-3-22") & date <= as.Date("2020-5-25"),1,
                                 ifelse(date >= as.Date("2020-5-26") & date <= as.Date("2020-6-30"), 0.5,
                                        ifelse(date >= as.Date("2020-9-1"),1,work_from_home))),
         business_closure = ifelse(date >= as.Date("2020-3-22") & date <= as.Date("2020-5-25"),1,
                                   ifelse(date >= as.Date("2020-6-01") & date <= as.Date("2020-9-1"), 0.5,
                                          ifelse(date >= as.Date("2020-12-15"),0.5,business_closure)))) %>%
  gather(key = phsm,value = value,-date) %>%
  mutate(phsm = recode(phsm,
                       school_closure = "School Closure",
                       work_from_home = "Work From Home",
                       mass_gather = "Restrict Mass Gatherings",
                       curfew = "Curfew",
                       business_closure = "Business Closure",
                       int_travel_ban = "International Travel Ban"),
         phsm = factor(phsm,levels = c("Curfew",
                                       "School Closure",
                                       "Work From Home",
                                       "Business Closure",
                                       "Restrict Mass Gatherings",
                                       "International Travel Ban"))) %>%
#  filter(!is.na(value))
         mutate(value = ifelse(is.na(value),0,value),
                line_width = ifelse(value>=0.5,1,0))

# Figure 3: Epidemiological indicators and PHSM measures as applied in the West Bank: PHSM data
top_plot <- response_data %>%
#  select(date,school_closing,work_closing) %>%
  ggplot(aes(x=date,y=phsm,colour=phsm,alpha=value)) +
    geom_line(size=6) +
    theme(legend.position = "blank") +
    xlab("") + ylab("PHSM") + theme(axis.text.x = element_blank(),
                     axis.ticks = element_blank(),
                     axis.text.y = element_text(size=16),
                     axis.title.y = element_text(size=20),
                     panel.background = element_rect(fill = "white",
                                                     colour = "black",
                                                     size = 0.5, linetype = "solid"),
                     panel.grid.major.x = element_line(size = 0.5, linetype = 'solid',
                                                       colour = "lightgrey"))
  
print(top_plot)

# combining all figures
multi_plot <- ggarrange(top_plot,middle_plot,bottom_plot,nrow=3,ncol=1,align="v")


# ------------------------------------------------------------------------------------------------------------------------
# Cross-correlation comparing Israel and the West Bank incidence data:

# load Israel incidence data
israel_wb_comparison <- read_csv("Mirror Mirror.csv") %>%
  select(date,Israel,Palestine) %>%
  rename(israel = "Israel",
         palestine = "Palestine") %>%
  mutate(date = as.Date(date,"%d/%m/%Y")) %>%
  arrange(-desc(date))

# cross-correlation function
ccf_results <- ccf(israel_wb_comparison$palestine,israel_wb_comparison$israel,lag=90,pl=FALSE)

ccf_df <- data.frame(lag = ccf_results$lag,
                     ccf = ccf_results$acf) %>%
  mutate(smooth_ccf = rollapply(ccf,5,mean,align='right',fill=NA))

# Figure 4: Correlation between Israeli and WB incidence data, against the temporal displacement between them.
ggplot(ccf_df,aes(x=lag,y=smooth_ccf)) +
  geom_line() +
  xlab("Displacement") + ylab("Correlation") +
  theme(axis.text.x = element_text(size=12),
        axis.title.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=16),
        panel.background = element_rect(fill = "white",
                                        colour = "black",
                                        size = 0.5, linetype = "solid")) +
  scale_x_continuous(limits = c(-90, 90),
                     breaks = c(-90, -60, -30, 0, 30, 60, 90),
                     label = c(-90, -60, -30, 0, 30, 60, 90)) +
  scale_y_continuous(limits = c(0, 0.9),
                     breaks = c(0, .3, .6, .9),
                     label = c(0, .3, .6, .9))



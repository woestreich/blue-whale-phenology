#####################################
# song_ts.R
# 
# Identifies periods of behavioral transition (foraging->migration) from CInight:CIday.
# Plots CI (song) time series (2015-2020) at daily resolution.
# Plots CInight:CIday time series (2015-2020) at daily resolution, highlighting periods of behavioral transition.
# Generates Figure 2.
#
# Will Oestreich
# Last update: May 1, 2021
#####################################

## clear variables and load packages
rm(list = ls())
library(tidyverse)
library(zoo)
library(patchwork)
library(lubridate)
library(padr)

## set smoothing windowsize and sample window size for t tests
windowsize <- 15
windowsize_ttest <- 30

## load in daily CI values and dates
ci_daily_seq <- read.csv("outputs/ci_daily_seq.csv",header = FALSE) 
ci_daily_seq <- t(ci_daily_seq)
ci_daily_seq <- data.frame(ci_daily_seq)%>%
  mutate(coln = seq(1,2173, 1))
colnames(ci_daily_seq) <- c("n","dd","d","all","ratio","coln")
a <- read.csv("outputs/ci_daily_seq_dates.csv", header = FALSE) %>%
  mutate(coln = seq(1, 2173, 1))
ci_daily_seq <- ci_daily_seq %>%
  left_join(a) %>%
  mutate(date = V1) %>%
  select(-V1, -coln)
ci_daily_seq$date <- as.Date(ci_daily_seq$date, format = "%d-%b-%Y")

## fill in missing dates w/ Nans
ci_daily_seq <- pad(ci_daily_seq)
ci_daily_seq$n[is.nan(ci_daily_seq$n)] <- NA
ci_daily_seq$d[is.nan(ci_daily_seq$d)] <- NA
ci_daily_seq$dd[is.nan(ci_daily_seq$dd)] <- NA
ci_daily_seq$all[is.nan(ci_daily_seq$all)] <- NA
ci_daily_seq$ratio[is.nan(ci_daily_seq$ratio)] <- NA

## run t tests to identify significant drops in CInight:CIday (behavioral transition)
dci_sig1 <- data.frame(as.numeric(rep(1,2223)))
colnames(dci_sig1) <- c("pv")
dci_sig2 <- data.frame(as.numeric(rep(1,2223)))
colnames(dci_sig2) <- c("pv")
for (d in (1):(2223-windowsize_ttest)) {
  if (d < 209+windowsize_ttest) {
    dci_sig1$pv[d] <- NA
    dci_sig2$pv[d] <- NA 
  }
  else {
    i1 <- d-windowsize_ttest
    i2 <- d-1
    j1 <- d+1
    j2 <- d+windowsize_ttest
    s1 <- ci_daily_seq$ratio[i1:i2]
    s2 <- ci_daily_seq$ratio[j1:j2]
    res1 <- t.test(s1,s2,alternative = "greater", var.equal = FALSE, na.omit =TRUE)
    res2 <- t.test(s1,s2,alternative = "less", var.equal = FALSE, na.omit =TRUE)
    dci_sig1$pv[d] <- res1$p.value
    dci_sig2$pv[d] <- res2$p.value
  }
}
ci_daily_seq$pv1 <- dci_sig1$pv
ci_daily_seq$pv2 <- dci_sig2$pv

## running mean
ci_daily_seq$ratio_rm <- rollapply(ci_daily_seq$ratio,windowsize,mean,fill=NA,na.rm = TRUE)
ci_daily_seq$ratio_rm[1:209] <- NA # correct for running mean creating "artificial" ratio values before the first date of recording (July 29, 2015)
ci_daily_seq$all_rm <- rollapply(ci_daily_seq$all,windowsize,mean,fill=NA,na.rm = TRUE)
ci_daily_seq$all_rm[1:209] <- NA # correct for running mean creating "artificial" values before the first date of recording (July 29, 2015)
ci_daily_seq$ratio_rm[ci_daily_seq$all_rm < 1.01] <- NA
ci_daily_seq$n_rm <- rollapply(ci_daily_seq$n,windowsize,mean,fill=NA,na.rm = TRUE)
ci_daily_seq$n_rm[1:209] <- NA # correct for running mean creating "artificial" values before the first date of recording (July 29, 2015)
ci_daily_seq$d_rm <- rollapply(ci_daily_seq$d,windowsize,mean,fill=NA,na.rm = TRUE)
ci_daily_seq$d_rm[1:209] <- NA # correct for running mean creating "artificial" values before the first date of recording (July 29, 2015)

## only consider CInight:CIday for days with blue whale song presence (CI >= 1.01)
ci_daily_seq$pv1[ci_daily_seq$all_rm < 1.01] <- NA
ci_daily_seq$pv2[ci_daily_seq$all_rm < 1.01] <- NA

## find dates of significant CInight:CIday decrease (restricts search until after mid-August, to avoid issues with start of recording and 2015, and to avoid consideration of first song presence for each year (arrival))
increase<- ci_daily_seq %>%
  filter(pv2<0.05)
decrease<- ci_daily_seq %>%
  filter(pv1<0.05)
song <- ci_daily_seq %>%
  filter(all_rm>1.01)
write.csv(song,"outputs/song.csv", row.names = FALSE) #save dates w/ song for later use 
decrease15 <- decrease %>%
  filter(date<as.Date("2016-03-01"))
decrease16 <- decrease %>%
  filter(date<as.Date("2017-03-01") & date>as.Date("2016-08-15"))
decrease17 <- decrease %>%
  filter(date<as.Date("2018-03-01") & date>as.Date("2017-08-15"))
decrease18 <- decrease %>%
  filter(date<as.Date("2019-03-01") & date>as.Date("2018-08-15"))
decrease19 <- decrease %>%
  filter(date<as.Date("2020-03-01") & date>as.Date("2019-08-15"))
decrease20 <- decrease %>%
  filter(date<as.Date("2021-03-01") & date>as.Date("2020-08-15"))

## save decrease as csv
write.csv(decrease,"outputs/decrease.csv", row.names = FALSE)

## generate Figure 2
tiff("outputs/Fig2.tiff",units="in", width=10,height=4,res=300)
# panel A
pa <- ci_daily_seq %>%
  select(date,n_rm,d_rm) %>% 
  pivot_longer(-date) %>%
  ggplot(aes(x=date,y = value)) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2015-07-30"), xmax = as.Date("2015-12-17"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2016-08-19"), xmax = as.Date("2017-01-26"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2017-07-08"), xmax = as.Date("2018-01-13"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2018-07-10"), xmax = as.Date("2019-01-29"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2019-07-26"), xmax = as.Date("2020-01-27"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2020-07-14"), xmax = as.Date("2021-01-24"),
           ymin = -Inf, ymax = Inf) +
  geom_line(aes(linetype = name)) +
  xlab("") +
  scale_x_date(date_breaks = "1 months",date_labels = "%b",date_minor_breaks = "1 month",expand = c(0, 0)) +
  theme_bw() + 
  theme(axis.text.x=element_text(size = 7, angle=60, hjust=1)) +
  theme(legend.position = "none") +
  geom_vline(xintercept=as.Date("2016-01-01")) +
  geom_vline(xintercept=as.Date("2017-01-01")) +
  geom_vline(xintercept=as.Date("2018-01-01")) +
  geom_vline(xintercept=as.Date("2019-01-01")) +
  geom_vline(xintercept=as.Date("2020-01-01")) +
  geom_vline(xintercept=as.Date("2021-01-01")) +
  ylim(0.98, 1.22) +
  ylab(expression(atop("Blue whale song ", paste("call index (CI)")))) +
  scale_linetype_manual(values = c(1, 3), labels = c("day", "night"), name = "") +
  annotate("text", label = "2015", x = as.Date("2015-05-25"), y = 1.2) +
  annotate("text", label = "2016", x = as.Date("2016-05-25"), y = 1.2) +
  annotate("text", label = "2017", x = as.Date("2017-05-25"), y = 1.2) +
  annotate("text", label = "2018", x = as.Date("2018-05-25"), y = 1.2) +
  annotate("text", label = "2019", x = as.Date("2019-05-25"), y = 1.2) +
  annotate("text", label = "2020", x = as.Date("2020-05-25"), y = 1.2) +
  annotate("text", label = "A", x = as.Date("2015-02-01"), y = 1.2, fontface = 2) +
  geom_segment(aes(x=as.Date("2015-01-15"), y=1.125, xend=as.Date("2015-04-15"), yend=1.125),color="black") +
  geom_segment(aes(x=as.Date("2015-01-15"), y=1.075, xend=as.Date("2015-04-15"), yend=1.075),color="black",linetype="dotted") +
  annotate("text", label = "Day", x = as.Date("2015-05-01"), y = 1.125, hjust=0) +
  annotate("text", label = "Night", x = as.Date("2015-05-01"), y = 1.075, hjust=0) 

# panel b
pb <- ci_daily_seq %>%
  select(date,ratio_rm,ratio) %>% 
  ggplot(aes(x=date,y=ratio_rm)) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2015-07-30"), xmax = as.Date("2015-12-17"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2016-08-19"), xmax = as.Date("2017-01-26"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2017-07-08"), xmax = as.Date("2018-01-13"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2018-07-10"), xmax = as.Date("2019-01-29"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2019-07-26"), xmax = as.Date("2020-01-27"),
           ymin = -Inf, ymax = Inf) +
  annotate("rect", fill = "gray", alpha = 0.4, 
           xmin = as.Date("2020-07-14"), xmax = as.Date("2021-01-24"),
           ymin = -Inf, ymax = Inf) +
  geom_line() + 
  geom_point(data = decrease15,aes(x=date,y=ratio_rm),color="red",size=1) +
  geom_point(data = decrease16,aes(x=date,y=ratio_rm),color="red",size=1) +
  geom_point(data = decrease17,aes(x=date,y=ratio_rm),color="red",size=1) +
  geom_point(data = decrease18,aes(x=date,y=ratio_rm),color="red",size=1) +
  geom_point(data = decrease19,aes(x=date,y=ratio_rm),color="red",size=1) +
  geom_point(data = decrease20,aes(x=date,y=ratio_rm),color="red",size=1) +
  xlab("") +
  scale_x_date(date_breaks = "1 months",date_labels = "%b",date_minor_breaks = "1 month",expand = c(0, 0)) +
  theme_bw() + 
  theme(axis.text.x=element_text(size = 7, angle=60, hjust=1)) +
  geom_vline(xintercept=as.Date("2016-01-01")) +
  geom_vline(xintercept=as.Date("2017-01-01")) +
  geom_vline(xintercept=as.Date("2018-01-01")) +
  geom_vline(xintercept=as.Date("2019-01-01")) +
  geom_vline(xintercept=as.Date("2020-01-01")) +
  geom_vline(xintercept=as.Date("2021-01-01")) +
  ylim(0.8, 2.1) +
  labs(y=expression(paste('CI'[night],':CI'[day])), x = "") +
  annotate("text", label = "2015", x = as.Date("2015-05-25"), y = 2) +
  annotate("text", label = "2016", x = as.Date("2016-05-25"), y = 2) +
  annotate("text", label = "2017", x = as.Date("2017-05-25"), y = 2) +
  annotate("text", label = "2018", x = as.Date("2018-05-25"), y = 2) +
  annotate("text", label = "2019", x = as.Date("2019-05-25"), y = 2) +
  annotate("text", label = "2020", x = as.Date("2020-05-25"), y = 2) +
  geom_point(aes(x=as.Date("2015-01-17"), y=1.7),color="red") +
  annotate("text", label = "Behavioral", x = as.Date("2015-02-01"), y = 1.7, hjust=0) +
  annotate("text", label = "transition", x = as.Date("2015-02-01"), y = 1.5, hjust=0) +
  annotate("text", label = "B", x = as.Date("2015-02-01"), y = 2, fontface = 2) 


# put them together
pa/pb
dev.off()
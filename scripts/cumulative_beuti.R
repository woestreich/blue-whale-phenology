#####################################
# cumulative_beuti.R
# 
# Calculates cumulative sums of BEUTI (Jacox et al., 2018) at 37N over study period (2015-2020). 
# Calculates phenological metrics of upwelling (Bograd et al., 2009) for each of these years.
# Plots and overlays climatological mean and extremes of BEUTI at 37N (1988-2020).
# Also generates a theoretical cumulative BEUTI curve for labeling and explaining phenology metrics.
# Generates Figure 3.
#
# Will Oestreich
# Last update: May 1, 2021
#####################################

## clear variables and load packages
rm(list = ls())
library(tidyverse)
library(patchwork)
library(RColorBrewer) 
library(zoo)
library(lubridate)
library(directlabels)

## set running mean window size
windowsize = 10

## load in daily beuti
beuti_daily <- read.csv("data/BEUTI_daily.csv",header = TRUE) 
beuti_daily$date <- as.Date(with(beuti_daily, paste(year, month, day,sep="-")), "%Y-%m-%d")

## generate idealized beuti curve for explanation of phenology metrics (Figure 3A)
A<-0
K<-2200
B<-0.03
v=1.5
Q=110
curve<-data.frame(matrix(ncol = 2, nrow = 366))
colnames(curve) <- c("yday","bt")
curve$doy<-as.numeric(1:366)
curve$bt<-0
for (i in 1:366) {
  x <- curve$doy[i]
  curve$bt[i] <- A + (K-A)/((1 + Q*exp(-B*x))^(1/v)) - 88
}

## calculate cumulative sum of beuti for each year of time series (1988-2020)
beuti_daily$yday <- yday(beuti_daily$date)
beuti_daily$csum <- 0
i1 <- 1
for (y in 1988:2020) {
  b <- beuti_daily %>% filter(year == y)
  bcsum <- cumsum(b$X37N)
  i2 <- i1 + length(bcsum) - 1
  beuti_daily$csum[i1:i2] <- bcsum
  i1 <- i2 + 1
}

## calculate climatological mean, 5th percentile, and 95th percentile for cumulative beuti curves 
beuti_clim <- data.frame(matrix(ncol = 4, nrow = 366))
colnames(beuti_clim) <- c("yday","csummean","csum5pctl","csum95pctl")
for (i in 1:366) {
  b <- beuti_daily %>% filter(yday == i)
  beuti_clim$yday[i] <- i
  beuti_clim$csummean[i] <- mean(b$csum,na.rm = TRUE)
  beuti_clim$csum5pctl[i] <- quantile(b$csum,.05,na.rm = TRUE)
  beuti_clim$csum95pctl[i] <- quantile(b$csum,.95,na.rm = TRUE)
}
# running mean to smooth
beuti_clim$csummean <- rollapply(beuti_clim$csummean,windowsize,mean,fill=NA,na.rm = TRUE)
beuti_clim$csum5pctl <- rollapply(beuti_clim$csum5pctl,windowsize,mean,fill=NA,na.rm = TRUE)
beuti_clim$csum95pctl <- rollapply(beuti_clim$csum95pctl,windowsize,mean,fill=NA,na.rm = TRUE)

## differences between smoothed daily values, needed for calculating STI, MAX, and END later
beuti_clim$diff[2:366] <- diff(beuti_clim$csummean)
beuti_clim$diff5[2:366] <- diff(beuti_clim$csum5pctl)
beuti_clim$diff95[2:366] <- diff(beuti_clim$csum95pctl)

## "clim" dataframe for storing STI, MAX, and END of the climatoligical mean BEUTI curve, used later in plotting Figure 3B
clim <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(clim) <- c("sti","sti_y","end","end_y","max","max_y")
clim$sti <- which.min(beuti_clim$csummean)
clim$sti_y <- beuti_clim$csummean[clim$sti]
clim$end <- which.max(beuti_clim$csummean)
clim$end_y <- beuti_clim$csummean[clim$end]
clim$max <- which.max(beuti_clim$diff)
clim$max_y <- beuti_clim$csummean[clim$max]
# set climatology "year" = 0 for coloring on plot
beuti_clim$year <- 0

## now calculate and combine climatoloigcal curves for 2015-2020 at 37 N, will be plotted as colored lines on Figure 3B
# 2015
beuti2015 <- filter(beuti_daily, year==2015)
beuti2015 <- select(beuti2015,c("year","date","X37N"))
beuti2015$csum <- cumsum(beuti2015$X37N)
beuti2015$csumrm <- rollapply(beuti2015$csum,windowsize,mean,fill=NA,na.rm = TRUE)
beuti2015$yday <- yday(beuti2015$date)
# 2016
beuti2016 <- filter(beuti_daily, year==2016)
beuti2016 <- select(beuti2016,c("year","date","X37N"))
beuti2016$csum <- cumsum(beuti2016$X37N)
beuti2016$csumrm <- rollapply(beuti2016$csum,windowsize,mean,fill=NA,na.rm = TRUE)
beuti2016$yday <- yday(beuti2016$date)
# 2017
beuti2017 <- filter(beuti_daily, year==2017)
beuti2017 <- select(beuti2017,c("year","date","X37N"))
beuti2017$csum <- cumsum(beuti2017$X37N)
beuti2017$csumrm <- rollapply(beuti2017$csum,windowsize,mean,fill=NA,na.rm = TRUE)
beuti2017$yday <- yday(beuti2017$date)
# 2018
beuti2018 <- filter(beuti_daily, year==2018)
beuti2018 <- select(beuti2018,c("year","date","X37N"))
beuti2018$csum <- cumsum(beuti2018$X37N)
beuti2018$csumrm <- rollapply(beuti2018$csum,windowsize,mean,fill=NA,na.rm = TRUE)
beuti2018$yday <- yday(beuti2018$date)
# 2019
beuti2019 <- filter(beuti_daily, year==2019)
beuti2019 <- select(beuti2019,c("year","date","X37N"))
beuti2019$csum <- cumsum(beuti2019$X37N)
beuti2019$csumrm <- rollapply(beuti2019$csum,windowsize,mean,fill=NA,na.rm = TRUE)
beuti2019$yday <- yday(beuti2019$date)
# 2020
beuti2020 <- filter(beuti_daily, year==2020)
beuti2020 <- select(beuti2020,c("year","date","X37N"))
beuti2020$csum <- cumsum(beuti2020$X37N)
beuti2020$csumrm <- rollapply(beuti2020$csum,windowsize,mean,fill=NA,na.rm = TRUE)
beuti2020$yday <- yday(beuti2020$date)
# combine
beuti_csum <- do.call("rbind", list(beuti2015,beuti2016,beuti2017,beuti2018,beuti2019,beuti2020))

## calculate phenological metrics (STI, MAX, END) for each year 2015-2020
# initialize data frames for cumulutive sum at each BEUTI phenological index (useful for accumulation plots)
sti <- data.frame(matrix(ncol = 4, nrow = 6))
colnames(sti) <- c("year","yday","date","csumrm")
maxi <- data.frame(matrix(ncol = 4, nrow = 6))
colnames(maxi) <- c("year","yday","date","csumrm")
endi <- data.frame(matrix(ncol = 4, nrow = 6))
colnames(endi) <- c("year","yday","date","csumrm")
# calculate phenological metrics for each year
for (i in 1:6) {
  yr1 <- i + 2014
  yr2 <- i + 2015
  
  beuti <- beuti_daily[which(beuti_daily$year >= yr1 & beuti_daily$year < yr2),]
  beuti$mbcsum <- cumsum(beuti$X37N)  #pull out Monterey Bay (37N)
  beuti$mbcsumrm <- rollapply(beuti$mbcsum,windowsize,mean,fill=NA,na.rm = TRUE) #running mean of cumulative sum
  beuti$mbrm <- rollapply(beuti$X37N,windowsize,mean,fill=NA,na.rm = TRUE) #running mean of daily values
  
  # STI, MAXI, and ENDI values for overlay on cumulative sum plots
  sti$year[i] <- yr1
  sti$yday[i] <- which.min(beuti$mbcsumrm)
  sti$date[i] <- as.Date(sti$yday[i], origin = make_date(yr1,1,1))
  sti$csumrm[i] <- beuti$mbcsumrm[sti$yday[i]]
  
  maxi$year[i] <- yr1
  maxi$yday[i] <- which.max(beuti$mbrm) 
  maxi$date[i] <- as.Date(maxi$yday[i], origin = make_date(yr1,1,1))
  maxi$csumrm[i] <- beuti$mbcsumrm[maxi$yday[i]]
  
  endi$year[i] <- yr1
  endi$yday[i] <- which.max(beuti$mbcsumrm)
  endi$date[i] <- as.Date(endi$yday[i], origin = make_date(yr1,1,1))
  endi$csumrm[i] <- beuti$mbcsumrm[endi$yday[i]]
}

## Generate and save figure
tiff("outputs/Figure3.tiff",units="in",width=10,height=4, res=400)
my_palette <- brewer.pal(name="Dark2",n=6)
my_palette <- palette(c(my_palette,"black"))
# Panel A
pa<-ggplot(curve,aes(doy,bt)) +
  geom_line(size=1) +
  geom_point(aes(x=1,y=0),size=4,shape=21,fill="white") +
  geom_point(aes(x=150,y=bt[150]),size=4,shape=24,fill="white") +
  geom_point(aes(x=366,y=bt[366]),size=4,shape=23,fill="white") +
  ylab("BEUTI cumulative sum\n(mmol/m/s)") +
  xlab("Day of year") +  
  theme_classic() +
  theme(legend.position = "none") +
  annotate("text", label = expression(TUMI == sum(BEUTI(t), i=STI, END)), x = 60, y = 2800) +
  annotate("text", label = expression(LUMI == sum(BEUTI(t), i=MAX, END)), x = 60, y = 2200) +
  annotate("text", label = "STI", x = 1, y = curve$bt[1]+200) +
  annotate("text", label = "MAX", x = 125, y = curve$bt[150]+100) +
  annotate("text", label = "END", x = 366, y = curve$bt[366]+200) +
  annotate("text", label = "A", x = 10, y = 3300, fontface = 2) +
  geom_segment(aes(x=150, y=-150, xend=150, yend=curve$bt[150]-20),size=1,color="black",linetype="dashed") +
  geom_segment(aes(x=366, y=-150, xend=366, yend=curve$bt[366]-20),size=1,color="black",linetype="dashed") +
  ylim(-150,3300) 
# Panel B
pb<-ggplot(beuti_clim, aes(yday,csummean)) +
  geom_ribbon(aes(ymin=csum5pctl,ymax=csum95pctl),fill = "grey70") +
  geom_line(color="black",size=1,linetype="dashed") +
  geom_point(data = clim,aes(x=sti,y=sti_y,color="black",size=3,shape=21)) +
  geom_point(data = clim,aes(x=max,y=max_y,color="black",size=3,shape=24)) +
  geom_point(data = clim,aes(x=end,y=end_y,color="black",size=3,shape=23)) +
  geom_line(data=beuti_csum, aes(x=yday,y=csumrm, color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(data = sti,aes(x=yday,y=csumrm,fill=as.factor(year),size=3,shape=21)) +
  geom_point(data = maxi,aes(x=yday,y=csumrm,fill=as.factor(year),size=3,shape=24)) +
  geom_point(data = endi,aes(x=yday,y=csumrm,fill=as.factor(year),size=3,shape=23)) +
  scale_fill_manual(values=my_palette) +
  scale_shape_identity() +
  annotate("text", label = "2015", x = 378, y = 1340, color = my_palette[1], size=3) +
  annotate("text", label = "2016", x = 378, y = 1600, color = my_palette[2], size=3) +
  annotate("text", label = "2017", x = 378, y = 1920, color = my_palette[3], size=3) +
  annotate("text", label = "2018", x = 378, y = 2190, color = my_palette[4], size=3) +
  annotate("text", label = "2019", x = 309, y = 1605, color = my_palette[5], size=3) +
  annotate("text", label = "2020", x = 378, y = 2525, color = my_palette[6], size=3) +
  ylab("BEUTI cumulative sum\n(mmol/m/s)") +
  xlab("Day of year") +
  theme_classic() +
  theme(legend.position = "none")  +
  geom_point(aes(x=20,y=2600),shape=21,color="black",fill="white",size=3) +
  geom_point(aes(x=20,y=2400),shape=24,color="black",fill="white",size=3) +
  geom_point(aes(x=20,y=2200),shape=23,color="black",fill="white",size=3) +
  geom_segment(aes(x=0, y=3000, xend=40, yend=3000),size=1,color="black",linetype="dashed") +
  geom_rect(aes(xmin=0,xmax=40,ymin=2700,ymax=2900),color="grey70",fill="grey70") + 
  annotate("text", label = "BEUTI[STI]", x = 74, y = 2590, parse = TRUE) +
  annotate("text", label = "BEUTI[MAX]", x = 77, y = 2390, parse = TRUE) +
  annotate("text", label = "BEUTI[END]", x = 77, y = 2190, parse = TRUE) +
  annotate("text", label = paste("Climatological","mean"), x = 112, y = 2990) +
  annotate("text", label = paste("Climatological","5th-95th","pctl"), x = 135, y = 2790) +
  annotate("text", label = "B", x = 10, y = 3300, fontface = 2)  +
  ylim(-150,3300) 
# Combine panels
pa+pb
dev.off()
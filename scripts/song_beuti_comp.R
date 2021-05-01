#####################################
# song_bueti_comp.R
# 
# Compares annual song presence and timing of foraging->migration behavioral transition to phenological metrics of upwelling.
# Generates Figure 4 and Supplemental Figure 3
#
# Will Oestreich
# Last update: May 1, 2021
#####################################

## clear variables and load packages
rm(list = ls())
library(tidyverse)
library(gridExtra)
library(zoo)
library(patchwork)
library(RColorBrewer) 
library(lubridate)
library(directlabels)

## load in daily beuti
beuti_daily <- read.csv("data/BEUTI_daily.csv",header = TRUE) 
beuti_daily$date <- as.Date(with(beuti_daily, paste(year, month, day,sep="-")), "%Y-%m-%d")

## load in acoustic song and migration metrics
migr <- read.csv("outputs/migr.csv",header = TRUE)

## set running mean window size
windowsize = 10

# initialize data frames for phenological metrics
d1 <- data.frame(matrix(ncol = 40, nrow = 6))
colnames(d1) <- c("year","sti","stirm","stirmdate","endi","endirm","endirmdate","maxi","maxirm","maxirmdate","tumi","tumi2","lusi","lusirm",
                  "sti_s","stirm_s","endi_s","endirm_s","maxi_s","maxirm_s","tumi_s","tumi2_s","lusi_s","lusirm_s",
                  "sti_n","stirm_n","endi_n","endirm_n","maxi_n","maxirm_n","tumi_n","tumi2_n","lusi_n","lusirm_n",
                  "migr_start","migr_end","migr_mean","song_start","song_peak","song_end")

# initialize data frames for cumulutive sum at each BEUTI phenological index (useful for accumulation plots)
sti <- data.frame(matrix(ncol = 8, nrow = 6))
colnames(sti) <- c("year","yday","date","csumrm","yday_s","csumrm_s","yday_n","csumrm_n")
maxi <- data.frame(matrix(ncol = 8, nrow = 6))
colnames(maxi) <- c("year","yday","date","csumrm","yday_s","csumrm_s","yday_n","csumrm_n")
endi <- data.frame(matrix(ncol = 8, nrow = 6))
colnames(endi) <- c("year","yday","date","csumrm","yday_s","csumrm_s","yday_n","csumrm_n")

## calculate phenological metrics for each year
for (i in 1:6) {
  yr1 <- i + 2014
  yr2 <- i + 2015
  
  # Monterey Bay BEUTI (37 N); also avg BEUTI for regions immediately north (38-39 N) and south (35-36 N)
  beuti <- beuti_daily[which(beuti_daily$year >= yr1 & beuti_daily$year < yr2),]
  beuti$mbcsum <- cumsum(beuti$X37N)
  beuti$mbcsumrm <- rollapply(beuti$mbcsum,windowsize,mean,fill=NA,na.rm = TRUE)
  beuti$mbrm <- rollapply(beuti$X37N,windowsize,mean,fill=NA,na.rm = TRUE)
  beuti$ncsum <- cumsum((beuti$X38N+beuti$X39N)/2)
  beuti$ncsumrm <- rollapply(beuti$ncsum,windowsize,mean,fill=NA,na.rm = TRUE)
  beuti$nrm <- rollapply((beuti$X38N+beuti$X39N)/2,windowsize,mean,fill=NA,na.rm = TRUE)
  beuti$scsum <- cumsum((beuti$X35N+beuti$X36N)/2)
  beuti$scsumrm <- rollapply(beuti$scsum,windowsize,mean,fill=NA,na.rm = TRUE)
  beuti$srm <- rollapply((beuti$X35N+beuti$X36N)/2,windowsize,mean,fill=NA,na.rm = TRUE)
  
  # Phenology metrics data frame (d1)
  d1$year[i] <- yr1
  
  d1$sti[i] <- which.min(beuti$mbcsum)  #MB metrics
  d1$stirm[i] <- which.min(beuti$mbcsumrm)
  d1$endi[i] <- which.max(beuti$mbcsum)
  d1$endirm[i] <- which.max(beuti$mbcsumrm)
  d1$maxi[i] <- which.max(beuti$X37N)
  d1$maxirm[i] <- which.max(beuti$mbrm)
  d1$tumi[i] <- sum(beuti$X37N[d1$stirm[i]:d1$endirm[i]])
  d1$tumi2[i] <- sum(beuti$X37N[d1$maxirm[i]:d1$endirm[i]])
  
  d1$sti_n[i] <- which.min(beuti$ncsum)  #N (38-39) metrics
  d1$stirm_n[i] <- which.min(beuti$ncsumrm)
  d1$endi_n[i] <- which.max(beuti$ncsum)
  d1$endirm_n[i] <- which.max(beuti$ncsumrm)
  d1$maxi_n[i] <- which.max((beuti$X38N+beuti$X39N))
  d1$maxirm_n[i] <- which.max(beuti$nrm)
  d1$tumi_n[i] <- sum(beuti$X37N[d1$stirm_n[i]:d1$endirm_n[i]])
  d1$tumi2_n[i] <- sum(beuti$X37N[d1$maxirm_n[i]:d1$endirm_n[i]])
  
  d1$sti_s[i] <- which.min(beuti$scsum)  #S (35-36) metrics
  d1$stirm_s[i] <- which.min(beuti$scsumrm)
  d1$endi_s[i] <- which.max(beuti$scsum)
  d1$endirm_s[i] <- which.max(beuti$scsumrm)
  d1$maxi_s[i] <- which.max((beuti$X35N+beuti$X36N))
  d1$maxirm_s[i] <- which.max(beuti$srm)
  d1$tumi_s[i] <- sum(beuti$X37N[d1$stirm_s[i]:d1$endirm_s[i]])
  d1$tumi2_s[i] <- sum(beuti$X37N[d1$maxirm_s[i]:d1$endirm_s[i]])
  
  d1$migr_start[i] <- migr$migr_start[i]  #Acoustic metrics
  d1$migr_start_date[i] <- as.Date(d1$migr_start[i], origin = make_date(yr1,1,1))
  d1$migr_end[i] <- migr$migr_end[i]
  d1$migr_end_date[i] <- as.Date(d1$migr_end[i], origin = make_date(yr1,1,1))
  d1$migr_mean[i] <- mean(c(d1$migr_start[i],d1$migr_end[i]))
  d1$migr_mean_date[i] <- as.Date(d1$migr_mean[i], origin = make_date(yr1,1,1))
  d1$song_start[i] <- migr$song_start[i]
  d1$song_start_date[i] <- as.Date(d1$song_start[i], origin = make_date(yr1,1,1))
  d1$song_peak[i] <- migr$song_peak[i]
  d1$song_peak_date[i] <- as.Date(d1$song_peak[i], origin = make_date(yr1,1,1))
  d1$song_end[i] <- migr$song_end[i]
  d1$song_end_date[i] <- as.Date(d1$song_end[i], origin = make_date(yr1,1,1))
  
  # STI, MAXI, and ENDI values for overlay on cumulative sum plots
  sti$year[i] <- yr1
  sti$yday[i] <- d1$stirm[i]  #MB values
  sti$date[i] <- as.Date(sti$yday[i], origin = make_date(yr1,1,1))
  sti$csumrm[i] <- beuti$mbcsumrm[d1$stirm[i]]
  sti$yday_n[i] <- d1$stirm_n[i]  #N values
  sti$csumrm_n[i] <- beuti$ncsumrm[d1$stirm_n[i]]
  sti$yday_s[i] <- d1$stirm_s[i]  #S values
  sti$csumrm_s[i] <- beuti$scsumrm[d1$stirm_s[i]]
  
  maxi$year[i] <- yr1
  maxi$yday[i] <- d1$maxirm[i]  #MB values
  maxi$date[i] <- as.Date(maxi$yday[i], origin = make_date(yr1,1,1))
  maxi$csumrm[i] <- beuti$mbcsumrm[d1$maxirm[i]]
  maxi$yday_n[i] <- d1$maxirm_n[i]  #N values
  maxi$csumrm_n[i] <- beuti$ncsumrm[d1$maxirm_n[i]]
  maxi$yday_s[i] <- d1$maxirm_s[i]  #S values
  maxi$csumrm_s[i] <- beuti$scsumrm[d1$maxirm_s[i]]
  
  endi$year[i] <- yr1
  endi$yday[i] <- d1$endirm[i] #MB values
  endi$date[i] <- as.Date(endi$yday[i], origin = make_date(yr1,1,1))
  endi$csumrm[i] <- beuti$mbcsumrm[d1$endirm[i]]
  endi$yday_n[i] <- d1$endirm_n[i]  #N values
  endi$csumrm_n[i] <- beuti$ncsumrm[d1$endirm_n[i]]
  endi$yday_s[i] <- d1$endirm_s[i]  #S values
  endi$csumrm_s[i] <- beuti$scsumrm[d1$endirm_s[i]]
}

# calculate LUSI (length of upwelling season) from STI and END
d1$lusirm <- d1$endirm - d1$stirm

## load in song presence dates and decrease dates (significant CInight:CIday decrease which defines behavioral transition)
song <- read.csv("outputs/song.csv",header = TRUE) 
decrease <- read.csv("outputs/decrease.csv",header = TRUE) 

## add day of year to decrease and song data frames
decrease$doy <- as.numeric(strftime(decrease$date, format = "%j"))
song$doy <- as.numeric(strftime(song$date, format = "%j"))

## add year to decrease and song data frames
decrease$year <- format(as.Date(decrease$date, format="%Y-%m-%d"),"%Y")
song$year <- format(as.Date(song$date, format="%Y-%m-%d"),"%Y")

## accounts for wraparound at end of year
decrease$doy[length(decrease$doy)] <- 367 
decrease$year[length(decrease$year)] <- 2020
song$month <- format(as.Date(song$date, format="%Y-%m-%d"),"%m")
for (i in 1:length(song$doy)) {
  if (as.numeric(song$month[i])<4) {
    if (as.numeric(song$year[i]==2017)) { #2017 to catch leap year (not 2016 as this is year + 1)
      song$doy[i] <- song$doy[i] + 366
    } else {
      song$doy[i] <- song$doy[i] + 365
    }
    yearn <- as.numeric(song$year[i])-1
    song$year[i]<-format(as.Date(paste(yearn, 1, 1, sep = "-")),"%Y")
  }
}

## remove July 23, 2015 from dataframe (before start of recording, artifact of the rolling mean applied earlier)
song<-song[!(song$date<as.Date("2015-07-30")),]

## add middle date of behavioral transition to decrease data frame
d15 <- decrease %>% filter(year==2015)
decrease$doycentr[decrease$year==2015]<-(max(d15$doy) + min(d15$doy))/2
d16 <- decrease %>% filter(year==2016)
decrease$doycentr[decrease$year==2016]<-(max(d16$doy) + min(d16$doy))/2
d17 <- decrease %>% filter(year==2017)
decrease$doycentr[decrease$year==2017]<-(max(d17$doy) + min(d17$doy))/2
d18 <- decrease %>% filter(year==2018)
decrease$doycentr[decrease$year==2018]<-(max(d18$doy) + min(d18$doy))/2
d19 <- decrease %>% filter(year==2019)
decrease$doycentr[decrease$year==2019]<-(max(d19$doy) + min(d19$doy))/2
d20 <- decrease %>% filter(year==2020)
decrease$doycentr[decrease$year==2020]<-(max(d20$doy) + min(d20$doy))/2

## add phenology metrics to decrease and song data frames
phenology<-merge(d1,decrease,by="year",all=TRUE)
songplus<-merge(d1,song,by="year",all=TRUE)

## Generate panels for Figure 4 and S3
my_palette <- brewer.pal(name="Dark2",n=6)
# extract legend
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# Max index (running mean version)
p1a <- ggplot(phenology, aes(x=maxirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=maxirm,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[MAX]*' (day of year)')) +
  ylab("Behavioral transition\n(day of year)") + 
  xlim(110,190) +
  annotate("text", label = "C", x = 111, y = 362, fontface = 2) +
  theme_classic() + theme(legend.title = element_blank()) + theme(legend.direction = "horizontal")
mylegend<-g_legend(p1a)

p1b <- ggplot(songplus, aes(x=maxirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=maxirm,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[MAX]*' (day of year)')) +
  ylab("Song presence\n(day of year)") + 
  xlim(110,190) +
  annotate("text", label = "A", x = 111, y = 390, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

# Spring transition index
p2a <- ggplot(phenology, aes(x=stirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=stirm,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[STI]*' (day of year)')) +
  ylab("Behavioral transition\n(day of year)") +  
  xlim(0,60) +
  annotate("text", label = "D", x = 1, y = 362, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

p2b <- ggplot(songplus, aes(x=stirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=stirm,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[STI]*' (day of year)')) +
  ylab("Song presence\n(day of year)") + 
  xlim(0,60) +
  annotate("text", label = "A", x = 1, y = 390, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

# End index
p3a <- ggplot(phenology, aes(x=endirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=endirm,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[END]*' (day of year)')) +
  ylab("") + 
  xlim(320,380) +
  annotate("text", label = "E", x = 321, y = 362, fontface = 2) +
  theme_classic()  + theme(legend.title = element_blank()) + theme(legend.direction = "horizontal")

p3b <- ggplot(songplus, aes(x=endirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=endirm,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) + 
  xlab(expression('BEUTI'[END]*' (day of year)')) +
  ylab("") + 
  xlim(320,380) +
  annotate("text", label = "B", x = 321, y = 390, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

# Total upwelling magnitude index
p4a <- ggplot(phenology, aes(x=tumi,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=tumi,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[TUMI]*' (mmol/s/m)')) +
  ylab("") + 
  xlim(1100,2500) +
  annotate("text", label = "D", x = 1125, y = 362, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

p4b <- ggplot(songplus, aes(x=tumi,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=tumi,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[TUMI]*' (mmol/s/m)')) +
  ylab("") + 
  xlim(1100,2500) +
  annotate("text", label = "B", x = 1125, y = 390, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

# Post-max total upwelling magnitude index (tumi2)
p5a <- ggplot(phenology, aes(x=tumi2,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=tumi2,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[LUMI]*' (mmol/s/m)')) +
  ylab("") + 
  xlim(500,1600) +
  annotate("text", label = "F", x = 525, y = 362, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

p5b <- ggplot(songplus, aes(x=tumi2,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=tumi2,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  xlab(expression('BEUTI'[LUMI]*' (mmol/s/m)')) +
  ylab("") + 
  xlim(500,1600) +
  annotate("text", label = "C", x = 525, y = 390, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

## 4-panel version of Figure 4, showing acoustic comparisons to MAX and TUMI
tiff("outputs/Fig4.tiff",units="in", width=7,height=6,res=300)
grid.arrange(p1b+theme(legend.position = 'none'),p4b,p1a+
               theme(legend.position = c(.65,.2),
                     legend.background = element_blank(),
                     legend.box.background = element_rect(color = "black")) +
               guides(color = guide_legend(nrow = 3)),
             p4a,nrow=2,ncol=2,
             heights=c(2.8,2.8),
             widths =c(2.95,2.8))
dev.off()

## 6-panel version of Figure S3, showing acoustic comparisons to STI, END, and LUMI
tiff("outputs/FigS3.tiff",units="in", width=10.5,height=6,res=300)
grid.arrange(p2b+theme(legend.position = 'none'),
             p3b,
             p5b,
             p2a,
             p3a+theme(legend.position = c(.3,.2),
                       legend.background = element_blank(),
                       legend.box.background = element_rect(color = "black")) +
               guides(color = guide_legend(nrow = 3)),
             p5a,
             nrow=2,ncol=3,
             heights=c(2.8,2.8),
             widths =c(2.95,2.8,2.8))
dev.off()

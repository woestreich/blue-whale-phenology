#####################################
# song_bueti_comp.R
# 
# Compares annual song presence and timing of foraging->migration behavioral transition to phenological metrics of upwelling.
# Generates Figure 4 and Supplemental Figures 3 & 4
#
# Will Oestreich
# Last update: June 21, 2021
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

## Statistics for Figure 4
# Panel A (transition vs. MAX)
x1 <- c(128,118,123,160,154,185) # MAX index
y1 <- c(248,256,231,322,296,340) # start of transition
y2 <- c(267.5,277.5,294,335.5,325.5,353.5) # center of transition
y3 <- c(287,299,357,349,355,367) # end of transition
# start of transition lm
mstart_a <- lm(y1~x1)
start_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_a_pts) <- c("x","y")
start_a_pts$x <- c(min(x1),max(x1))
start_a_pts$y <- c(mstart_a$coefficients[1] + mstart_a$coefficients[2]*min(x1),
                 mstart_a$coefficients[1] + mstart_a$coefficients[2]*max(x1))
# center of transition lm
mcenter_a <- lm(y2~x1)
center_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(center_a_pts) <- c("x","y")
center_a_pts$x <- c(min(x1),max(x1))
center_a_pts$y <- c(mcenter_a$coefficients[1] + mcenter_a$coefficients[2]*min(x1),
                   mcenter_a$coefficients[1] + mcenter_a$coefficients[2]*max(x1))
#end of transition lm
mend_a <- lm(y3~x1)
end_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_a_pts) <- c("x","y")
end_a_pts$x <- c(min(x1),max(x1))
end_a_pts$y <- c(mend_a$coefficients[1] + mend_a$coefficients[2]*min(x1),
                    mend_a$coefficients[1] + mend_a$coefficients[2]*max(x1))

# Panel B (transition vs. TUMI)
x2 <- c(1220.241,1538.166,1910.764,2055.806,1645.540,2400.422) #TUMI
#start of transition lm
mstart_b <- lm(y1~x2)
start_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_b_pts) <- c("x","y")
start_b_pts$x <- c(min(x2),max(x2))
start_b_pts$y <- c(mstart_b$coefficients[1] + mstart_b$coefficients[2]*min(x2),
                   mstart_b$coefficients[1] + mstart_b$coefficients[2]*max(x2))
#center of transition lm
mcenter_b <- lm(y2~x2)
center_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(center_b_pts) <- c("x","y")
center_b_pts$x <- c(min(x2),max(x2))
center_b_pts$y <- c(mcenter_b$coefficients[1] + mcenter_b$coefficients[2]*min(x2),
                   mcenter_b$coefficients[1] + mcenter_b$coefficients[2]*max(x2))
#end of transition lm
mend_b <- lm(y3~x2)
end_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_b_pts) <- c("x","y")
end_b_pts$x <- c(min(x2),max(x2))
end_b_pts$y <- c(mend_b$coefficients[1] + mend_b$coefficients[2]*min(x2),
                    mend_b$coefficients[1] + mend_b$coefficients[2]*max(x2))

## Generate Figure 4 
my_palette <- brewer.pal(name="Dark2",n=6)
# extract legend
# https://github.com/hadley/ggplot2/wiki/Share-a-legend-between-two-ggplot2-graphs
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# Max index (running mean version)
pa <- ggplot(phenology, aes(x=maxirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=maxirm,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_a_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=center_a_pts, aes(x=x,y=y)) +
  geom_line(data=end_a_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 188, y = 347, label = paste('R^2 ==',round(summary(mstart_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 207.5, y = 347, label = "*", size=3, parse = FALSE, hjust = 0) +
  annotate("text", x = 188, y = 360, label = paste('R^2 ==',round(summary(mcenter_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 207.5, y = 360, label = "*", size=3, parse = FALSE, hjust = 0) +
  annotate("text", x = 188, y = 373, label = paste("R^2 ==",round(summary(mend_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  xlab(expression('BEUTI'[MAX]*' (yearday)')) +
  ylab("Behavioral transition\n(yearday)") + 
  xlim(110,210) +
  ylim(230,380) +
  annotate("text", label = "A", x = 111, y = 370, fontface = 2) +
  theme_classic() + theme(legend.title = element_blank()) + theme(legend.direction = "horizontal")
mylegend<-g_legend(p1a)

# Total upwelling magnitude index
pb <- ggplot(phenology, aes(x=tumi,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=tumi,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_b_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=center_b_pts, aes(x=x,y=y)) +
  geom_line(data=end_b_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 2450, y = 326, label = paste('R^2 ==',round(summary(mstart_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 2450, y = 353, label = paste('R^2 ==',round(summary(mcenter_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 2780, y = 353, label = "*", size=3, parse = FALSE, hjust = 0) +
  annotate("text", x = 2450, y = 380, label = paste("R^2 ==",round(summary(mend_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 2780, y = 380, label = "*", size=3, parse = FALSE, hjust = 0) +
  xlab(expression('BEUTI'[TUMI]*' (mmol/m/s)')) +
  ylab("") + 
  xlim(1100,2800) +
  ylim(230,380) +
  annotate("text", label = "B", x = 1125, y = 370, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

## 2-panel version of Figure 4, showing behavioral transition in comparison to MAX and TUMI
tiff("outputs/Fig4.tiff",units="in", width=7,height=3,res=300)
grid.arrange(pa+theme(legend.position = c(.65,.2),
                       legend.background = element_blank(),
                       legend.box.background = element_rect(color = "black")) +
               guides(color = guide_legend(nrow = 3)),
             pb,nrow=1,ncol=2,
             heights=c(2.8),
             widths =c(2.95,2.8))
dev.off()


## Supplemental figure 3
# S3 stats
x1 <- c(128,118,123,160,154,185) # MAX index
x2 <- c(1220.241,1538.166,1910.764,2055.806,1645.540,2400.422) #TUMI
y1 <- c(211,230,190,189,205,197) # start of song
y2 <- c(268,302,307,293,315,310) # peak of song
y3 <- c(252,389,380,392,391,386) # end of song

# start of song vs. max lm
mstart_a <- lm(y1~x1)
start_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_a_pts) <- c("x","y")
start_a_pts$x <- c(min(x1),max(x1))
start_a_pts$y <- c(mstart_a$coefficients[1] + mstart_a$coefficients[2]*min(x1),
                   mstart_a$coefficients[1] + mstart_a$coefficients[2]*max(x1))
# peak of song vs. max lm
mpeak_a <- lm(y2~x1)
peak_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(peak_a_pts) <- c("x","y")
peak_a_pts$x <- c(min(x1),max(x1))
peak_a_pts$y <- c(mpeak_a$coefficients[1] + mpeak_a$coefficients[2]*min(x1),
                  mpeak_a$coefficients[1] + mpeak_a$coefficients[2]*max(x1))
# end of song vs. max lm
mend_a <- lm(y3~x1)
end_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_a_pts) <- c("x","y")
end_a_pts$x <- c(min(x1),max(x1))
end_a_pts$y <- c(mend_a$coefficients[1] + mend_a$coefficients[2]*min(x1),
                 mend_a$coefficients[1] + mend_a$coefficients[2]*max(x1))
# start of song vs. tumi lm
mstart_b <- lm(y1~x2)
start_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_b_pts) <- c("x","y")
start_b_pts$x <- c(min(x2),max(x2))
start_b_pts$y <- c(mstart_b$coefficients[1] + mstart_b$coefficients[2]*min(x2),
                   mstart_b$coefficients[1] + mstart_b$coefficients[2]*max(x2))
# peak of song vs. tumi lm
mpeak_b <- lm(y2~x2)
peak_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(peak_b_pts) <- c("x","y")
peak_b_pts$x <- c(min(x2),max(x2))
peak_b_pts$y <- c(mpeak_b$coefficients[1] + mpeak_b$coefficients[2]*min(x2),
                  mpeak_b$coefficients[1] + mpeak_b$coefficients[2]*max(x2))
# end of song vs. tumi lm
mend_b <- lm(y3~x2)
end_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_b_pts) <- c("x","y")
end_b_pts$x <- c(min(x2),max(x2))
end_b_pts$y <- c(mend_b$coefficients[1] + mend_b$coefficients[2]*min(x2),
                 mend_b$coefficients[1] + mend_b$coefficients[2]*max(x2))

# song vs. max
s3a <- ggplot(songplus, aes(x=maxirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=maxirm,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_a_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=peak_a_pts, aes(x=x,y=y)) +
  geom_line(data=end_a_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 188, y = 194, label = paste('R^2 ==',round(summary(mstart_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 188, y = 310, label = paste('R^2 ==',round(summary(mpeak_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 188, y = 395, label = paste("R^2 ==",round(summary(mend_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  xlab(expression('BEUTI'[MAX]*' (yearday)')) +
  ylab("Song presence\n(yearday)") + 
  xlim(85,205) +
  ylim(185,420) +
  annotate("text", label = "A", x = 87, y = 405, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

# song vs. tumi
s3b <- ggplot(songplus, aes(x=tumi,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=tumi,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_b_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=peak_b_pts, aes(x=x,y=y)) +
  geom_line(data=end_b_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 2450, y = 194, label = paste('R^2 ==',round(summary(mstart_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 2450, y = 315, label = paste('R^2 ==',round(summary(mpeak_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 2450, y = 417, label = paste("R^2 ==",round(summary(mend_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  xlab(expression('BEUTI'[TUMI]*' (mmol/m/s)')) +
  ylab("") + 
  xlim(1100,2650) +
  ylim(185,420) +
  annotate("text", label = "B", x = 1125, y = 405, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

## 2-panel version of Figure S3, showing song presence in comparison to MAX and TUMI
tiff("outputs/FigS3.tiff",units="in", width=8,height=3,res=300)
grid.arrange(s3a+theme(legend.position = c(.13,.45),
                       legend.background = element_blank(),
                       legend.title = element_blank(),
                       legend.box.background = element_rect(color = "black")) +
               guides(color = guide_legend(nrow = 6)),
             s3b,
             nrow=1,ncol=2,
             heights=c(2.8),
             widths =c(2.95,2.8))
dev.off()

## Figure S4
# S4 stats
x1 <- c(43,26,42,5,20,5) # STI
x2 <- c(360,361,360,360,329,361) #END
x3 <- c(600,1188.1,1485.8,996.8,1021.6,1018.5) #LUMI
y1 <- c(211,230,190,189,205,197) # start of song
y2 <- c(268,302,307,293,315,310) # peak of song
y3 <- c(252,389,380,392,391,386) # end of song
y4 <- c(248,256,231,322,296,340) # start of transition
y5 <- c(267.5,277.5,294,335.5,325.5,353.5) # center of transition
y6 <- c(287,299,357,349,355,367) # end of transition
# start of song vs. sti lm
mstart_a <- lm(y1~x1)
start_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_a_pts) <- c("x","y")
start_a_pts$x <- c(min(x1),max(x1))
start_a_pts$y <- c(mstart_a$coefficients[1] + mstart_a$coefficients[2]*min(x1),
                   mstart_a$coefficients[1] + mstart_a$coefficients[2]*max(x1))
# peak of song vs. sti lm
mpeak_a <- lm(y2~x1)
peak_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(peak_a_pts) <- c("x","y")
peak_a_pts$x <- c(min(x1),max(x1))
peak_a_pts$y <- c(mpeak_a$coefficients[1] + mpeak_a$coefficients[2]*min(x1),
                   mpeak_a$coefficients[1] + mpeak_a$coefficients[2]*max(x1))
# end of song vs. sti lm
mend_a <- lm(y3~x1)
end_a_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_a_pts) <- c("x","y")
end_a_pts$x <- c(min(x1),max(x1))
end_a_pts$y <- c(mend_a$coefficients[1] + mend_a$coefficients[2]*min(x1),
                 mend_a$coefficients[1] + mend_a$coefficients[2]*max(x1))
# start of song vs. end lm
mstart_b <- lm(y1~x2)
start_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_b_pts) <- c("x","y")
start_b_pts$x <- c(min(x2),max(x2))
start_b_pts$y <- c(mstart_b$coefficients[1] + mstart_b$coefficients[2]*min(x2),
                   mstart_b$coefficients[1] + mstart_b$coefficients[2]*max(x2))
# peak of song vs. end lm
mpeak_b <- lm(y2~x2)
peak_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(peak_b_pts) <- c("x","y")
peak_b_pts$x <- c(min(x2),max(x2))
peak_b_pts$y <- c(mpeak_b$coefficients[1] + mpeak_b$coefficients[2]*min(x2),
                  mpeak_b$coefficients[1] + mpeak_b$coefficients[2]*max(x2))
# end of song vs. end lm
mend_b <- lm(y3~x2)
end_b_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_b_pts) <- c("x","y")
end_b_pts$x <- c(min(x2),max(x2))
end_b_pts$y <- c(mend_b$coefficients[1] + mend_b$coefficients[2]*min(x2),
                 mend_b$coefficients[1] + mend_b$coefficients[2]*max(x2))
# start of song vs. lumi lm
mstart_c <- lm(y1~x3)
start_c_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_c_pts) <- c("x","y")
start_c_pts$x <- c(min(x3),max(x3))
start_c_pts$y <- c(mstart_c$coefficients[1] + mstart_c$coefficients[2]*min(x3),
                   mstart_c$coefficients[1] + mstart_c$coefficients[2]*max(x3))
# peak of song vs. lumi lm
mpeak_c <- lm(y2~x3)
peak_c_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(peak_c_pts) <- c("x","y")
peak_c_pts$x <- c(min(x3),max(x3))
peak_c_pts$y <- c(mpeak_c$coefficients[1] + mpeak_c$coefficients[2]*min(x3),
                  mpeak_c$coefficients[1] + mpeak_c$coefficients[2]*max(x3))
# end of song vs. lumi lm
mend_c <- lm(y3~x3)
end_c_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_c_pts) <- c("x","y")
end_c_pts$x <- c(min(x3),max(x3))
end_c_pts$y <- c(mend_c$coefficients[1] + mend_c$coefficients[2]*min(x3),
                 mend_c$coefficients[1] + mend_c$coefficients[2]*max(x3))
# start of transition vs. sti lm
mstart_d <- lm(y4~x1)
start_d_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_d_pts) <- c("x","y")
start_d_pts$x <- c(min(x1),max(x1))
start_d_pts$y <- c(mstart_d$coefficients[1] + mstart_d$coefficients[2]*min(x1),
                   mstart_d$coefficients[1] + mstart_d$coefficients[2]*max(x1))
# center of transition vs. sti lm
mcenter_d <- lm(y5~x1)
center_d_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(center_d_pts) <- c("x","y")
center_d_pts$x <- c(min(x1),max(x1))
center_d_pts$y <- c(mcenter_d$coefficients[1] + mcenter_d$coefficients[2]*min(x1),
                    mcenter_d$coefficients[1] + mcenter_d$coefficients[2]*max(x1))
# end of transition vs. sti lm
mend_d <- lm(y6~x1)
end_d_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_d_pts) <- c("x","y")
end_d_pts$x <- c(min(x1),max(x1))
end_d_pts$y <- c(mend_d$coefficients[1] + mend_d$coefficients[2]*min(x1),
                 mend_d$coefficients[1] + mend_d$coefficients[2]*max(x1))
# start of transition vs. end lm
mstart_e <- lm(y4~x2)
start_e_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_e_pts) <- c("x","y")
start_e_pts$x <- c(min(x2),max(x2))
start_e_pts$y <- c(mstart_e$coefficients[1] + mstart_e$coefficients[2]*min(x2),
                   mstart_e$coefficients[1] + mstart_e$coefficients[2]*max(x2))
# center of transition vs. end lm
mcenter_e <- lm(y5~x2)
center_e_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(center_e_pts) <- c("x","y")
center_e_pts$x <- c(min(x2),max(x2))
center_e_pts$y <- c(mcenter_e$coefficients[1] + mcenter_e$coefficients[2]*min(x2),
                    mcenter_e$coefficients[1] + mcenter_e$coefficients[2]*max(x2))
# end of transition vs. end lm
mend_e <- lm(y6~x2)
end_e_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_e_pts) <- c("x","y")
end_e_pts$x <- c(min(x2),max(x2))
end_e_pts$y <- c(mend_e$coefficients[1] + mend_e$coefficients[2]*min(x2),
                 mend_e$coefficients[1] + mend_e$coefficients[2]*max(x2))
# start of transition vs. lumi lm
mstart_f <- lm(y4~x3)
start_f_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(start_f_pts) <- c("x","y")
start_f_pts$x <- c(min(x3),max(x3))
start_f_pts$y <- c(mstart_f$coefficients[1] + mstart_f$coefficients[2]*min(x3),
                   mstart_f$coefficients[1] + mstart_f$coefficients[2]*max(x3))
# center of transition vs. lumi lm
mcenter_f <- lm(y5~x3)
center_f_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(center_f_pts) <- c("x","y")
center_f_pts$x <- c(min(x3),max(x3))
center_f_pts$y <- c(mcenter_f$coefficients[1] + mcenter_f$coefficients[2]*min(x3),
                    mcenter_f$coefficients[1] + mcenter_f$coefficients[2]*max(x3))
# end of transition vs. lumi lm
mend_f <- lm(y6~x3)
end_f_pts <- data.frame(matrix(ncol = 2, nrow = 2))
colnames(end_f_pts) <- c("x","y")
end_f_pts$x <- c(min(x3),max(x3))
end_f_pts$y <- c(mend_f$coefficients[1] + mend_f$coefficients[2]*min(x3),
                 mend_f$coefficients[1] + mend_f$coefficients[2]*max(x3))

s4a <- ggplot(songplus, aes(x=stirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=stirm,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_a_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=peak_a_pts, aes(x=x,y=y)) +
  geom_line(data=end_a_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 45, y = 210, label = paste('R^2 ==',round(summary(mstart_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 45, y = 293, label = paste('R^2 ==',round(summary(mpeak_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 45, y = 328, label = paste("R^2 ==",round(summary(mend_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  xlab(expression('BEUTI'[STI]*' (yearday)')) +
  ylab("Song presence\n(yearday)") + 
  xlim(0,60) +
  annotate("text", label = "A", x = 1, y = 390, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

s4b <- ggplot(songplus, aes(x=endirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=endirm,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) + 
  geom_line(data=start_b_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=peak_b_pts, aes(x=x,y=y)) +
  geom_line(data=end_b_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 364, y = 205, label = paste('R^2 ==',round(summary(mstart_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 364, y = 298, label = paste('R^2 ==',round(summary(mpeak_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 364, y = 360, label = paste("R^2 ==",round(summary(mend_b)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  xlab(expression('BEUTI'[END]*' (yearday)')) +
  ylab("") + 
  xlim(320,380) +
  annotate("text", label = "B", x = 321, y = 390, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

s4c <- ggplot(songplus, aes(x=tumi2,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=tumi2,y=song_peak,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_c_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=peak_c_pts, aes(x=x,y=y)) +
  geom_line(data=end_c_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 1520, y = 202, label = paste('R^2 ==',round(summary(mstart_c)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 1520, y = 318, label = paste('R^2 ==',round(summary(mpeak_c)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 1520, y = 425, label = paste("R^2 ==",round(summary(mend_c)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  xlab(expression('BEUTI'[LUMI]*' (mmol/m/s)')) +
  ylab("") + 
  xlim(500,1750) +
  annotate("text", label = "C", x = 525, y = 440, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

s4d <- ggplot(phenology, aes(x=stirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=stirm,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_d_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=center_d_pts, aes(x=x,y=y)) +
  geom_line(data=end_d_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 45, y = 236, label = paste('R^2 ==',round(summary(mstart_d)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 45, y = 275, label = paste('R^2 ==',round(summary(mcenter_d)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 45, y = 315, label = paste("R^2 ==",round(summary(mend_d)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 56.5, y = 236, label = "*", size=3, parse = FALSE, hjust = 0) +
  annotate("text", x = 56.5, y = 275, label = "*", size=3, parse = FALSE, hjust = 0) +
  xlab(expression('BEUTI'[STI]*' (yearday)')) +
  ylab("Behavioral transition\n(yearday)") +  
  xlim(0,60) +
  ylim(230,375) +
  annotate("text", label = "D", x = 1, y = 370, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

s4e <- ggplot(phenology, aes(x=endirm,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=endirm,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_e_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=center_e_pts, aes(x=x,y=y)) +
  geom_line(data=end_e_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 364, y = 280, label = paste('R^2 ==',round(summary(mstart_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 364, y = 307, label = paste('R^2 ==',round(summary(mpeak_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 364, y = 335, label = paste("R^2 ==",round(summary(mend_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  xlab(expression('BEUTI'[END]*' (yearday)')) +
  ylab("") + 
  xlim(320,380) +
  ylim(230,375) +
  annotate("text", label = "E", x = 321, y = 370, fontface = 2) +
  theme_classic()  + theme(legend.title = element_blank()) + theme(legend.direction = "horizontal")

s4f <- ggplot(phenology, aes(x=tumi2,y=doy)) +
  geom_point(aes(color=as.factor(year)),size=1) +
  scale_color_manual(values=my_palette) +
  geom_point(aes(x=tumi2,y=doycentr,fill=as.factor(year)),size=4,shape=22) +
  scale_fill_manual(values=my_palette) +
  geom_line(data=start_f_pts, aes(x=x,y=y), linetype = "dashed") +
  geom_line(data=center_f_pts, aes(x=x,y=y)) +
  geom_line(data=end_f_pts, aes(x=x,y=y), linetype = "dashed") +
  annotate("text", x = 1520, y = 268, label = paste('R^2 ==',round(summary(mstart_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 1520, y = 314, label = paste('R^2 ==',round(summary(mpeak_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  annotate("text", x = 1520, y = 361, label = paste("R^2 ==",round(summary(mend_a)$r.squared, 2)), size=3, parse = TRUE, hjust = 0) +
  xlab(expression('BEUTI'[LUMI]*' (mmol/m/s)')) +
  ylab("") + 
  xlim(500,1750) +
  ylim(230,375) +
  annotate("text", label = "F", x = 525, y = 370, fontface = 2) +
  theme_classic() + theme(legend.position = 'none')

## 6-panel version of Figure S4, showing acoustic comparisons to STI, END, and LUMI
tiff("outputs/FigS4.tiff",units="in", width=10.5,height=6,res=300)
grid.arrange(s4a+theme(legend.position = 'none'),
             s4b,
             s4c,
             s4d,
             s4e+theme(legend.position = c(.3,.2),
                       legend.background = element_blank(),
                       legend.box.background = element_rect(color = "black")) +
               guides(color = guide_legend(nrow = 3)),
             s4f,
             nrow=2,ncol=3,
             heights=c(2.8,2.8),
             widths =c(2.95,2.8,2.8))
dev.off()




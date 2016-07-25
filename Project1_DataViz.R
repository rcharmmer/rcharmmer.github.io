library(dplyr)   ## load dplyr package for data preprocessing
## read the biomarker data set
df1 <- read.csv("C:/Users/ricky/Dropbox/interviews/gilead/Biomarker Data.csv")
## read the averse effect time data set
df2 <- read.csv("C:/Users/ricky/Dropbox/interviews/gilead/AE times.csv")
## join the two data sets
df3 <- left_join(df1,df2,by=c("SUBJECT"="ID"))
## filter those records to be deleted based on adverse effect data set
df4 <- filter(df3,WEEK>EVENT.TIME)
## generate the final data frame to be used to identify useful markers
df5 <- anti_join(df3,df4,by=names(df3))

## add a new column with base marker.value(marker.value at week 0) to the data frame
library(data.table)                                        ## lode data.table package 
key <- paste(df5$SUBJECT,df5$MARKER.NAME,df5$WEEK,sep="_") ## build a unique key for each record with SUBJECT,MARKER.NAME and WEEK info. 
dt  <- data.table(cbind(df5,key))                          ## combine the key to the data frame
setkey(dt,key)                                             ## set the key as index to each record
new_key <- paste(df5$SUBJECT,df5$MARKER.NAME,0,sep="_")    ## add a new key with SUBJECT and MARKER.NAME info, with WEEK=0.

base <- c()                                ## initiate a new column which has MARKER.VALUE when WEEK=0

##look up the table and get MARKER.VALUE indexed with the new key and return base MARKER.VALUE(WEEK=0) for all records 
for (i in 1:nrow(dt)){
  base[i] <- dt[new_key[i]]$MARKER.VALUE  
}

## check the distribution of base value
summary(base)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.00   10.00   22.00   26.03   38.00  103.00

## build the metric FOLD
## since base has zero values, 1 is added to both the numerator and denominator of FOLD=(MARKER.VALUE/BASE), 
## the FOLD is calculated as below:
FOLD <- (df5$MARKER.VALUE+1)/(base+1)

df <- mutate(df5,FOLD)  ## add FOLD to the data frame

## group the data by MARKER.NAME and WEEK, cacluate the mean and sd for MARKER.VALUE, and attach those two values to table 
df_wide<- summarise(group_by(df,MARKER.NAME,WEEK),AVG.MARKER.VALUE=mean(MARKER.VALUE),MARKER.VALUE.SD=sd(MARKER.VALUE),AVG.FOLD=mean(FOLD),FOLD.SD=sd(FOLD))


library("ggplot2") ## load ggplot for data visualization

## draw the "marker.value vs week" line chart
value.chart <- ggplot(data=df_wide, aes(x=WEEK, y=AVG.MARKER.VALUE, group = MARKER.NAME, colour = MARKER.NAME)) +
  geom_line() + 
  labs(x="Week",y="Marker Value")+
  geom_point(size=6, shape=20, fill="white") + 
  geom_text(aes(label=ifelse(df_wide$WEEK==7,as.character(df_wide$MARKER.NAME),"")),hjust=2, vjust=0,size=5) +
  ggtitle("Change of Marker Value during the Treatment")+
  theme(plot.title = element_text(size=32,face="bold",color="deepskyblue2"))+
  theme(legend.text = element_text(size=20))+
  theme(axis.title = element_text(size=28,face="bold",color="deepskyblue2"))+
  theme(axis.text.x= element_text(size=24))+
  theme(axis.text.y= element_text(size=24))+
  theme(legend.title=element_blank())+
  geom_vline(xintercept = 0.5,linetype="longdash")+
  annotate("text", x = 0, y = 90, label = "Control",size=12,color="gray17")+
  annotate("text", x = 2, y = 90, label = "Treatment",size=12,color="orangered1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


value.chart  ## show the chart

## draw the "fold vs week" line chart
fold.chart <- ggplot(data=df_wide, aes(x=WEEK, y=AVG.FOLD, group = MARKER.NAME, colour = MARKER.NAME)) +
  geom_line() + 
  labs(x="Week",y="Fold")+
  geom_point(size=6, shape=20, fill="white") + 
  geom_text(aes(label=ifelse(df_wide$WEEK==7,as.character(df_wide$MARKER.NAME),"")),hjust=2, vjust=0,size=5) +
  ggtitle("Marker Fold Change during the Treatment")+
  theme(plot.title = element_text(size=32,face="bold",color="deepskyblue2"))+
  theme(legend.text = element_text(size=20))+
  theme(axis.title = element_text(size=28,face="bold",color="deepskyblue2"))+
  theme(axis.text.x= element_text(size=24))+
  theme(axis.text.y= element_text(size=24))+
  theme(legend.title=element_blank())+
  geom_vline(xintercept = 0.5,linetype="longdash")+
  annotate("text", x = 0, y = 3.5, label = "Control",size=12,color="gray17")+
  annotate("text", x = 2, y = 3.5, label = "Treatment",size=12,color="orangered1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  

fold.chart  ## show the chart

## make a short list of markers, where M1,M2,M3,M4 are selected 
df_wide_short <- filter(df_wide,MARKER.NAME %in% c("M1","M2","M3","M4"))

## prepare graphs of "MARKER.VALUE vs WEEK" with error bars
## calculate error bars of MARKER.VALUE
limits.value <- aes(ymax = df_wide_short$AVG.MARKER.VALUE + df_wide_short$MARKER.VALUE.SD, ymin = df_wide_short$AVG.MARKER.VALUE - df_wide_short$MARKER.VALUE.SD,
                   colour = df_wide_short$MARKER.NAME)
## draw the "MARKER.VALUE vs WEEK" graph with error bars
value.chart.short <- ggplot(data=df_wide_short, aes(x=WEEK, y=AVG.MARKER.VALUE, group = MARKER.NAME, colour = MARKER.NAME)) +
  geom_line() + 
  labs(x="Week",y="Marker Value")+
  geom_point(size=6, shape=20, fill="white") + 
  facet_wrap(~ MARKER.NAME) +
  geom_errorbar(limits.value, width=0.2)+
  guides(colour=FALSE)+
  ggtitle("Change of Marker Value during the Treatment")+
  theme(plot.title = element_text(size=32,face="bold",color="deepskyblue2"))+
  theme(legend.text = element_text(size=20))+
  theme(axis.title = element_text(size=28,face="bold",color="deepskyblue2"))+
  theme(axis.text.x= element_text(size=24))+
  theme(axis.text.y= element_text(size=24))+
  theme(strip.text = element_text(size=25))+
  theme(legend.title=element_blank())+
  geom_vline(xintercept = 0.5,linetype="longdash")+
  annotate("text", x = 0, y = 150, label = "Control",size=7,color="gray17")+
  annotate("text", x = 2, y = 150, label = "Treatment",size=7,color="orangered1")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

value.chart.short


## prepare graphs of "FOLD vs WEEK" with error bars
## calculate error bars of FOLD
limits.fold <- aes(ymax = df_wide_short$AVG.FOLD + df_wide_short$FOLD.SD, ymin = df_wide_short$AVG.FOLD - df_wide_short$FOLD.SD,
                   colour = df_wide_short$MARKER.NAME)

## draw the "FOLD vs WEEK" graph
fold.chart.short <- ggplot(data=df_wide_short, aes(x=WEEK, y=AVG.FOLD, group = MARKER.NAME, colour = MARKER.NAME)) +
  geom_line() + 
  labs(x="Week",y="Fold")+
  geom_point(size=6, shape=20, fill="white") + 
  facet_wrap(~ MARKER.NAME) +
  geom_errorbar(limits.fold, width=0.2)+
  guides(colour=FALSE)+
  ggtitle("Fold Change during the Treatment")+
  theme(plot.title = element_text(size=32,face="bold",color="deepskyblue2"))+
  theme(legend.text = element_text(size=20))+
  theme(axis.title = element_text(size=28,face="bold",color="deepskyblue2"))+
  theme(axis.text.x= element_text(size=24))+
  theme(axis.text.y= element_text(size=24))+
  theme(strip.text = element_text(size=25))+
  theme(legend.title=element_blank())+
  geom_vline(xintercept = 0.5,linetype="longdash")+
  annotate("text", x = 0, y = 5.5, label = "Control",size=7,color="gray17")+
  annotate("text", x = 2, y = 5.5, label = "Treatment",size=7,color="orangered1")+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

fold.chart.short


## compare the fold change of the short listed markers between subjects having adverse effect 
## and not having adverse effect
df_ae <- filter(df,SUBJECT%in%c(1:12)) ## get SUBJECT(ID=1:12) which has AE 
df_ne <- filter(df,SUBJECT%in%c(13:39)) ## get SUBHECT(ID=13:39) which has no AE

Event <- as.factor(c(rep(1,160),rep(0,180)))  ## label SUBJECT with regard to AE(EFFECT=1 if AE, EFFEC=0 if no AE)
## get records with AE and compute the mean and std.error of fold
ae <- summarize(group_by(df_ae,WEEK,MARKER.NAME),avg.fd=mean(FOLD),se=sd(FOLD)/sqrt(length(FOLD)))
## get records without AE and compute the mean and std.error of fold
ne <- summarize(group_by(df_ne,WEEK,MARKER.NAME),avg.fd=mean(FOLD),se=sd(FOLD)/sqrt(length(FOLD)))

dff <- rbind(ae,ne) ## bind two tables above
dfff <- cbind(dff,Event)  ## bind the EFFECT label to table
dffff <- filter(dfff,MARKER.NAME%in%c("M1","M2","M3","M4")) ## only choose MARKER.NAME= M1,M2,M3,M4.
## draw the "fold vs week" line chart between the two group (Subject with AE vs Subject without AE) with std. error
short.list <- ggplot(data=dffff, aes(x=WEEK, y=avg.fd, group = Event, colour = Event)) +
  geom_line() + 
  labs(x="Week",y="Fold")+
  geom_point(size=6, shape=20, fill="white") + 
  facet_wrap(~ MARKER.NAME) +
  geom_errorbar(aes(ymin=avg.fd-se, ymax=avg.fd+se), width=.1) +
  ggtitle("Comparing Fold Change between Group w/ and w/o Adverse Event")+
  theme(plot.title = element_text(size=32,face="bold",color="deepskyblue2"))+
  theme(legend.text = element_text(size=20))+
  theme(axis.title = element_text(size=28,face="bold",color="deepskyblue2"))+
  theme(axis.text.x= element_text(size=24))+
  theme(axis.text.y= element_text(size=24))+
  theme(strip.text = element_text(size=25))+
  theme(legend.title=element_text(size=20))+
  geom_vline(xintercept = 0.5,linetype="longdash")+
  annotate("text", x = 0, y = 3.5, label = "Control",size=7,color="gray17")+
  annotate("text", x = 2, y = 3.5, label = "Treatment",size=7,color="orangered1")+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

short.list

## box plot to compare fold change one week before adverse effect
df.bp <- df %>% filter(MARKER.NAME %in% c("M1","M2","M3","M4"))  ## get rows with MARKER.NAME as M1,M2,M3 and M4
df.bp <- df.bp %>% mutate(LABEL=ifelse(is.na(EVENT.TIME),0,1))   ## attach adverse event label as a new column
df.bp.w356 <- df.bp %>% filter(WEEK%in%c(3,5,6))   ## get rows with WEEK=3
df.bp.w356$WEEK <- as.factor(df.bp.w356$WEEK)
df.bp.w356$LABEL <- as.factor(df.bp.w356$LABEL)


w356.bp <- ggplot(df.bp.w356,aes(x=WEEK ,y=FOLD,color=LABEL,fill=LABEL)) +
  facet_wrap(~MARKER.NAME,scales="free") + 
  theme_bw()+
  geom_point(position=position_jitterdodge(dodge.width=0.9)) +
  geom_boxplot(fill="white",outlier.colour = NA, position = position_dodge(width=0.9))+
  ggtitle("Marker Fold One Week Before the Event(Group w/ Event vs. Group w/o Event)")+
  theme(plot.title = element_text(size=32,face="bold",color="deepskyblue2"))+
  theme(legend.text = element_text(size=20))+
  theme(axis.title = element_text(size=28,face="bold",color="deepskyblue2"))+
  theme(axis.text.x= element_text(size=24))+
  theme(axis.text.y= element_text(size=24))+
  theme(strip.text = element_text(size=25))+
  theme(legend.title=element_text(size=20))+
  labs(x="Week",y="Fold")+
  scale_fill_discrete(name = "Event")+
  guides(color=guide_legend(title="Event"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

w356.bp


## two sample t-test between two groups at week 6
m3.fold <- df.bp.w356 %>% filter(MARKER.NAME=="M3"&WEEK==6) %>% select(FOLD,LABEL)
m3.fold.0 <- m3.fold[which(m3.fold$LABEL==0),]$FOLD
m3.fold.1 <- m3.fold[which(m3.fold$LABEL==1),]$FOLD

t.test(m3.fold.0,m3.fold.1)






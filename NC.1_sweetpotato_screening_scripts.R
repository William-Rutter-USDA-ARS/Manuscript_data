###Notes### 
#Scripts used for analysis of Phase2 screening of sweetpotato germplasm using M.e. isolate NC.1
#Supplemental data for manuscript entitled "Identification of Sweetpotato Germplasm Resistant to Pathotypically Distinct Isolates of Meloidogyne enterolobii from the Carolinas" Rutter et.al. Plant Disease 2020

###Required Libraries###
library(readxl)
library(ggplot2)

###Load data###
Combined_data <- read_excel("Menterolobii_NC.1_combinded_sweetpotato_screens.xlsx")
#formatting
Combined_data$Screen<-as.factor(Combined_data$Screen)


###Data analysis###

## susceptable control data alone
sus<- c("COVINGTON", "BEAUREGARD")
sus_ind<-which(Combined_data$Entry %in% sus)
sus_data<-Combined_data[sus_ind,]

#generate dot plot to view log egg per gram root between screens for susceptable controls
pl.sus.leggrtwt<-ggplot(sus_data, aes(x=reorder(Entry, legg.rtwt), y=legg.rtwt, color=Screen, fill = Screen)) + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.0) + theme_bw(base_size = 16) + theme(axis.text.x=element_text(angle=90,hjust=1))
pl.sus.leggrtwt

#generate dot plot to view Galling between screens for susceptable controls
Pl.sus.Gall<-ggplot(sus_data, aes(x=reorder(Entry, Gall), y=Gall, color=Screen, fill = Screen)) + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.0) + theme_bw(base_size = 16) + theme(axis.text.x=element_text(angle=90,hjust=1))
Pl.sus.Gall

#Checking to see if Beauregaurd was significantly different between screens 1,2,&3
beau.data<-sus_data[sus_data$Accession == "BEAUREGARD",]
#egg data anova
anova.beau.log.egg.rtwt<- aov(legg.rtwt ~ Screen, data = beau.data) 
summary(anova.beau.log.egg.rtwt)
a.tuk<-TukeyHSD(anova.beau.log.egg.rtwt)
log.beau.legg.rtwt.diffs<-as.data.frame(a.tuk$Screen) 
log.beau.legg.rtwt.diffs #not significantly for egg data

#galling data anova
anova.beau.gall<- aov(Gall ~ Screen, data = beau.data) 
summary(anova.beau.gall)
b.tuk<-TukeyHSD(anova.beau.gall)
log.beau.gall.diffs<-as.data.frame(b.tuk$Screen) 
log.beau.gall.diffs #not significantly for galling data

#Checking to see if Covington was sig different between screens
cov.data<-sus_data[sus_data$Accession == "COVINGTON",]
#egg data anova
anova.cov.log.egg.rtwt<- aov(legg.rtwt ~ Screen, data = cov.data) 
summary(anova.cov.log.egg.rtwt)
c.tuk<-TukeyHSD(anova.cov.log.egg.rtwt)
log.cov.legg.rtwt.diffs<-as.data.frame(c.tuk$Screen)
log.cov.legg.rtwt.diffs  #not significant for egg data

#gall data anova 
anova.cov.gall<- aov(Gall ~ Screen, data = cov.data) 
summary(anova.cov.gall)
d.tuk<-TukeyHSD(anova.cov.gall)
log.cov.gall.diffs<-as.data.frame(d.tuk$Screen) #not significant for egg data
log.cov.gall.diffs

#checking to see if Covington and Beauregard were significantly different in terms of eggs or galling scores accross three screens
#sus legg.rt.wt combined
anova.sus.legg.rtwt<- aov(legg.rtwt ~ Entry, data = sus_data)
summary(anova.sus.legg.rtwt)
e.tuk<-TukeyHSD(anova.sus.legg.rtwt)
sus.legg.rtwt.diffs<-as.data.frame(e.tuk$Entry) 
sus.legg.rtwt.diffs #adjusted p < 0.001

#sus galling combined 
anova.sus.gall<- aov(Gall ~ Entry, data = sus_data)
summary(anova.sus.gall)
f.tuk<-TukeyHSD(anova.sus.gall)
sus.Gall.diffs<-as.data.frame(f.tuk$Entry) 
sus.Gall.diffs #adjusted p < 0.05

##analysis accross all sweetpotato lines

#generate a dot plot to view Galling between lines accross all screens
Pl.all.Gall<-ggplot(Combined_data, aes(x=reorder(Entry, Gall), y=Gall)) + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.0) + theme_bw(base_size = 16) + theme(axis.text.x=element_text(angle=90,hjust=1))
Pl.all.Gall

#generate dot plot to view log egg per gram root between screens for susceptable controls
pl.all.leggrtwt<-ggplot(Combined_data, aes(x=reorder(Entry, legg.rtwt), y=legg.rtwt)) + geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.0) + theme_bw(base_size = 16) + theme(axis.text.x=element_text(angle=90,hjust=1))
pl.all.leggrtwt

#Galling
anova.Galling<- aov(Gall ~ Entry, data = Combined_data) 
summary(anova.Galling)
g.tuk<-TukeyHSD(anova.Galling)
Galling.diffs<-as.data.frame(g.tuk$Entry) 
Galling.diffs

#Eggs/gram root
anova.legg.rtwt<- aov(legg.rtwt ~ Entry, data = Combined_data) 
summary(anova.legg.rtwt)
h.tuk<-TukeyHSD(anova.legg.rtwt)
legg.rtwt.diffs<-as.data.frame(h.tuk$Entry) 
legg.rtwt.diffs

#aggregating the mean values for Table 2 isolated NC.1
#eggs per gram of root
m.egg.rtwt<-aggregate(egg.rtwt~Entry, FUN=mean, data = Combined_data)
sd.egg.rtwt<-aggregate(egg.rtwt~Entry, FUN=sd, data = Combined_data)
summary(m.egg.rtwt$egg.rtwt)
df.egg.rtwt<-merge(m.egg.rtwt, sd.egg.rtwt, by="Entry")
df.names<-c("Entry", "mean_eggs/gram_root", "sd_eggs/gram_root" )
colnames(df.egg.rtwt)<-df.names

# percent galling
m.Gall<-aggregate(Gall~Entry, FUN=mean, data = Combined_data)
sd.Gall<-aggregate(Gall~Entry, FUN=sd, data = Combined_data)
summary(m.Gall$Gall)
df.Gall<-merge(m.Gall, sd.Gall, by="Entry")
df.names<-c("Entry", "mean_percent_galling", "sd_percent_galling" )
colnames(df.Gall)<-df.names













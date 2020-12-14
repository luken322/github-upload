library(vegan)
library(ggplot2)
library(reshape)
library(tidyr)
library(plyr)
library(dplyr)
library(tibble)
library(ggplot2)
?reshape
setwd("C://Users/Heather//Documents//MS Research")
zoops <- read.csv("ZoopBiom_zoop_Fisher.csv")

min.year <- 1996
months.to.consider <- c(5, 6, 7,8)

zoop <- rename(zoops,
                Adult.Copepods=TOTAL.ADULT.COPEPODS,
                Copepods=TOTAL.COPEPODS,
                Crustaceans=TOTAL.CRUSTACEANS,
                Cladocerans=TOTAL.CLADOCERANS,
               Chydoridae=Chydorus)

zoops <- zoop %>%
  filter(Year>=min.year,
         Month %in% months.to.consider) %>% 
  select(-X0, -Sample..,-Season, -Date,-Month,-Nauplii,-Adult.Copepods,-Copepods,-Crustaceans,-Cladocerans,-TOTAL.ZOOPLANKTON) %>%
  pivot_longer(cols=-c(Year), names_to="Measure") %>%
  group_by(Year,Measure) %>%
  summarize(mean=mean(value)) %>%
  pivot_wider(id_cols=c(Year), names_from=Measure, values_from=mean) 

zoops <- column_to_rownames(zoops, 'Year')

zoops.matrix <- data.matrix(zoops)

yrs <- c(1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 
         2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 
         2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)

#Default dissimilarity metric is Bray-Curtis, excludes joint absences 
#Joint absences= any species/taxa that is @ 0 in both groups in pairwise comparison will not be included in calc.
#ie: two communities are no more similar when both do not have a rare species
#Exclusion of joint absences tends to inflate B-diversity (variation) when there are many pairwise 0's
#This is not our case because we only have one taxa that has a lot of pairwise 0's (Scapholeberis)

#Test different metrics?
#rankindex(yrs, wisconsin(sqrt(zoops.matrix)))

#Since our values are large, NMDS performs both square root transformation and Wisconsin tr. (output range: 0-1)
#Set minimum runs with try and max runs with trymax for best solution w/low stress, previous.best to set run 0 (standard)
#To get the same convergence, same ordination points, set seed here (but first explore different starting points), use set seed 343, try=100,trymax=500
set.seed(343)
zoops.nmds <- metaMDS(zoops.matrix,k=2,try=100,trymax =500,previous.best = TRUE)
??metaMDS
summary(zoops.nmds)
names(zoops.nmds)
zoops.nmds$stress
zoops.nmds$distance
zoops.nmds$points
zoops.nmds$data
#Dissimilarity values
View(zoops.nmds$diss)
zoops.nmds$iidx

ordisurf(zoops.nmds, yrs, main = "annual", col="black")

#look at a Shepard plot, showing scatter around the regression between the interpoint distances in the final configuration vs. original dissimilarities
stressplot(zoops.nmds)
#if there is large scatter about line that suggests original dissimilarities are not well preserved in reduced dimensions, but looks good here

plot(zoops.nmds)

ordiplot(zoops.nmds, type="n", main="Annual")
orditorp(zoops.nmds, display = "species", col="red", cex=1.5,air=0.01)
orditorp(zoops.nmds, display="sites", cex=1.25, air=0.01)

#telling R to use 3 groups of 8,8,8 years when drawing similarity polygons (ordihull), can add species and sites, mark years with color 
treat3 <-c(rep("dec1", 8), rep("dec2",8), rep("dec3", 8))
ordiplot(zoops.nmds, type="n", main="Annual")
ordihull(zoops.nmds, groups=treat3, draw="polygon",col="grey",label=F)
orditorp(zoops.nmds, display = "species", col="red", air=0.01,cex=1.5)
orditorp(zoops.nmds, display="sites", col=c(rep("dark green",8),rep("blue",8),rep("orange",8)),cex=1.25, air=0.01)

#Try ordiellipse
treat3 <-c(rep("1996-2003", 8), rep("2004-2011",8), rep("2012-2019", 8))
ordiplot(zoops.nmds, type="n", main="Annual")
ordiellipse(zoops.nmds, groups=treat3,col=c(rep("dark green",1),rep("blue",1),rep("orange",1)),label=F,kind="se",conf = 0.95)
orditorp(zoops.nmds, display = "species", col="red", air=0.01,cex=1.3)
orditorp(zoops.nmds, display="sites", col=c(rep("dark green",8),rep("blue",8),rep("orange",8)),cex=0.95, air=0.01)

#2 groups of 12
treat2 <- c(rep("dec1", 12), rep("dec2",12))
ordiplot(zoops.nmds, type="n", main="Annual")
ordihull(zoops.nmds, groups=treat2, draw="polygon",col="grey",label=F)
orditorp(zoops.nmds, display = "species", col="red", air=0.01,cex=1.5)
orditorp(zoops.nmds, display="sites", col=c(rep("dark green",12),rep("blue",12)),cex=1.25, air=0.01)

############GGplot ordiellipse plot

NMDS = data.frame(NMDS1 = zoops.nmds$points[,1], NMDS2 = zoops.nmds$points[,2],group=treat3)

ord<-ordiellipse(zoops.nmds, treat3, display = "sites",
                 kind = "se", conf = 0.95, label =F) 
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                   veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}
View(NMDS)
species.scores1 <- as.data.frame(scores(zoops.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores1$species <- rownames(species.scores1)
??geom_point
NMDS_Annual <-ggplot(data = NMDS, aes(NMDS1, NMDS2)) +
  geom_point(aes(shape=group, alpha=group),size=1.18, stroke=0.48) +
  scale_alpha_manual(values=c("1996-2003"=1, "2004-2011"=0.6,"2012-2019"=0.6))+
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,linetype=group), size=0.28)+
  geom_text(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),alpha=1.2,size=2.2) +
  scale_shape_manual(values=c("1996-2003"=3, "2004-2011"=16, "2012-2019"=2)) +
  scale_linetype_manual(values=c("1996-2003"="dashed", "2004-2011"="solid", "2012-2019"="dotted"))+
  #scale_colour_manual(values=c("1996-2003 Sum"="red","1996-2003 Spr"="blue","2004-2011 Sum"="red","2004-2011 Spr"="blue","2012-2019 Sum"= "red","2012-2019 Spr"= "blue")) +
  coord_equal() +
  scale_x_continuous(breaks=c(-0.8,-0.4,0,0.4,0.8),limits=c(-0.8,0.9))+
  scale_y_continuous(breaks=c(-0.8,-0.4,0,0.4,0.8),limits=c(-0.8,0.9))+
  theme_bw(11) + 
  theme(axis.text = element_text(size=8, color="black"),  # remove x-axis text
        axis.ticks = element_blank(),  # remove axis ticks
    axis.title = element_text(size=9), 
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.width = unit(2.4, "line"),
    legend.text = element_text(size=7),
    rect = element_rect(colour='black',size=1.2))

setwd("C://Users/Heather//Documents//MS Research//Chapter 1//Rplots")
ggsave(plot = NMDS_Annual, "NMDS_Annual.jpg", device="jpg", h = 4.2/1, w = 4.2/1,dpi=1200)

############Seasonal NMDS################
zoopses <- zoop %>%
  filter(Year>=min.year,
         Month %in% months.to.consider) %>% 
  select(-X0, -Sample.., -Date,-Month,-Nauplii,-Adult.Copepods,-Copepods,-Crustaceans,-Cladocerans,-TOTAL.ZOOPLANKTON) %>%
  pivot_longer(cols=-c(Year,Season), names_to="Measure") %>%
  group_by(Year,Season,Measure) %>%
  summarize(mean=mean(value)) %>%
  pivot_wider(id_cols=c(Year,Season), names_from=Measure, values_from=mean) %>%
  arrange(Season,Year) 

levels(zoopses$Season) <- c("a", "b")

zoopses <- zoopses %>%
  unite(Year_Season,c("Year", "Season"),sep = "")
  

zoopses <- column_to_rownames(zoopses, 'Year_Season')
View(zoopses)
zoopses.matrix <- data.matrix(zoopses)
set.seed(35)
zoopses.nmds <- metaMDS(zoopses.matrix, k=2, try=100,trymax = 500, previous.best = TRUE)

summary(zoopses.nmds)
names(zoopses.nmds)
zoopses.nmds$stress
zoopses.nmds$diss
zoopses.nmds$points
zoopses.nmds$data

yrses <- c(1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 
         2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 
         2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019, 
         1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 
         2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 
         2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019)

ordisurf(zoopses.nmds, yrses, main = "annual", col="black")

#look at a Shepard plot, showing scatter around the regression between the interpoint distances in the final configuration vs. original dissimilarities
stressplot(zoopses.nmds)
#if there is large scatter about line that suggests original dissimilarities are not well preserved in reduced dimensions, but looks good here

plot(zoopses.nmds)

ordiplot(zoopses.nmds, type="n", main="Seasonal")
orditorp(zoopses.nmds, display = "species", col="red", air=0.01, cex=1.2)
orditorp(zoopses.nmds,display="sites", col=c(rep("blue",24),rep("dark green",24)),cex=1, air=0.01)

#telling R to use 4 groups of 8 yrs spring, 12 years spr, same for summer
#when drawing similarity polygons (ordihull), can add species and sites, mark seasons with color 
treat4 <-c(rep("spr1", 8), rep("spr2", 16), rep("sum3", 8), rep("sum4", 16))
ordiplot(zoopses.nmds, type="n", main="Seasonal", xlim=c(-0.5,0.5))
ordihull(zoopses.nmds, groups=treat4, draw="polygon",col=c(rep("lightblue1",2),rep("pale green",2)),label=F)
orditorp(zoopses.nmds, display = "species", col="red",cex=1.2, air=0.01)
orditorp(zoopses.nmds, display="sites", col=c(rep("blue",24),rep("dark green",24)),cex=1, air=0.01)
?ordiplot

treat5 <-c(rep("1996-2003 Spr", 8), rep("2004-2011 Spr", 8), rep("2012-2019 Spr",8),rep("1996-2003 Sum", 8), rep("2004-2011 Sum", 8),rep("2012-2019 Sum",8))
ordiplot(zoopses.nmds, type="n", main="Seasonal", xlim=c(-0.5,0.5))
ordiellipse(zoopses.nmds, groups=treat5, col=c(rep("blue",3),rep("green",3)),label=F,kind="se",conf = 0.95)
orditorp(zoopses.nmds, display = "species", col="red",cex=1.2, air=0.01)
orditorp(zoopses.nmds, display="sites", col=c(rep("blue",24),rep("dark green",24)),cex=1, air=0.01)

#####GGplot ordiellipse plot

NMDS = data.frame(NMDS1 = zoopses.nmds$points[,1], NMDS2 = zoopses.nmds$points[,2],group=treat5)

ord<-ordiellipse(zoopses.nmds, treat5, display = "sites",
                 kind = "se", conf = 0.95, label =F) 
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

df_ell <- data.frame()
for(g in levels(NMDS$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(NMDS[NMDS$group==g,],
                                                  veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                ,group=g))
}

species.scores1 <- as.data.frame(scores(zoopses.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores1$species <- rownames(species.scores1)

NMDS_Seasonal <- ggplot(data = NMDS, aes(NMDS1, NMDS2)) + geom_point(aes(color = group,shape=group),size=1.18, stroke=0.48) +
  geom_path(data=df_ell, aes(x=NMDS1, y=NMDS2,colour=group,linetype=group), size=0.28)+
  geom_text(data=species.scores1,aes(x=NMDS1,y=NMDS2,label=species),alpha=1,size=2.2) +
  scale_shape_manual(values=c("1996-2003 Spr"=3,"1996-2003 Sum"=3, "2004-2011 Spr"=16,"2004-2011 Sum"=16, "2012-2019 Spr"=2,"2012-2019 Sum"=2)) +
  scale_linetype_manual(values=c("1996-2003 Spr"="dashed","1996-2003 Sum"="dashed", "2004-2011 Spr"="solid","2004-2011 Sum"="solid", "2012-2019 Spr"="dotted","2012-2019 Sum"="dotted"))+
  scale_colour_manual(values=c("1996-2003 Sum"="red","1996-2003 Spr"="blue","2004-2011 Sum"="red","2004-2011 Spr"="blue","2012-2019 Sum"= "red","2012-2019 Spr"= "blue")) +
  coord_equal() +
  scale_x_continuous(breaks=c(-0.8,-0.4,0,0.4,0.8),limits=c(-0.8,0.9))+
  scale_y_continuous(breaks=c(-0.8,-0.4,0,0.4,0.8),limits=c(-0.8,0.9))+
  theme_bw(11) + 
  theme(axis.text = element_text(size=8, color="black"),  # remove x-axis text
    axis.ticks = element_blank(),  # remove axis ticks
    axis.title = element_text(size=9), 
    panel.background = element_blank(), 
    panel.grid.major = element_blank(),  #remove major-grid labels
    panel.grid.minor = element_blank(),  #remove minor-grid labels
    plot.background = element_blank(),
    legend.title = element_blank(),
    legend.position = "top",
    legend.key.width = unit(2.4, "line"),
    legend.text = element_text(size=7),
    rect=element_rect(colour='black',size=1.2))

setwd("C://Users/Heather//Documents//MS Research//Chapter 1//Rplots")
ggsave(plot = NMDS_Seasonal, "NMDS_Seasonal.jpg", device="jpg", h = 4.4/1, w = 4.4/1,dpi=1200)

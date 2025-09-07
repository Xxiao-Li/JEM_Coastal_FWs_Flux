rm(list=ls())
library(vegan)
library(ggplot2)
library(ggsci)
library(eoffice)
library(cowplot)
setwd("E:/00-food web functioning and stability/Data_manuscript/PCoA")
PCoAdata <- read.csv("./PCoAdata.csv",header = TRUE,fileEncoding = "UTF-8-BOM")
rownames(PCoAdata)<-paste0("site",1:70)
PCoAdata[is.na(PCoAdata)]<-0
PCoAplot <- PCoAdata[,3:16]
bray_dist <- vegdist(PCoAplot,method = "bray")

library(ape)
df.pcoa <- pcoa(bray_dist,correction = "cailliez")
df.pcoa$vectors
df.pcoa$values
df.plot<-data.frame(df.pcoa$vectors)
head(df.plot)

x_label<-round(df.pcoa$values$Rel_corr_eig[1]*100,2)
y_label<-round(df.pcoa$values$Rel_corr_eig[2]*100,2)
x_label
y_label


df.plot$group<-PCoAdata$Site
df.plot$month <- PCoAdata$Month
#,shape=month
p1 <- ggplot(data=df.plot,aes(x=Axis.1,y=Axis.2))+
  geom_point(aes(fill = group, shape = month),size = 3.5,alpha=0.87)+
  theme_bw()+
  theme(panel.grid = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),
       y=paste0("PCoA2 ",y_label,"%"))+
  stat_ellipse(data=df.plot,
               geom = "polygon",
               aes(fill=group),
               level=0.95,
               alpha=0.2)+
   scale_shape_manual(values = c(21:24))+
   scale_fill_d3()+
   scale_color_manual(values = c("black","black","black"))+
   theme(axis.text = element_text(family="Arial", size=16),
         legend.text = element_text(family="Arial", size=16),
         axis.title = element_text(family="Arial", size=16),
         legend.title = element_text(family="Arial", size=16))

p1

topptx(p1,"E:/00-food web functioning and stability/Data_manuscript/PCoA/p1.pptx",height = 6,width = 6)
  

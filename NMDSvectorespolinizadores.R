############################
#############################
library(vegan)
library(permute)
library(lattice)
##polinizadores composition
pol <- read.table(file.choose(), header = TRUE)
head(pol) # base de datos pol

amb<- read.table(file.choose(), header = TRUE)
head(amb) # vamb


#Analysis of Similarities
pol.dist<-vegdist(pol, "bray") # para datos de abundancia se usa bray, la primer celda de la base de datos debe estar VACIA 
pol.dist
summary(pol.dist)

sitios<- read.table(file.choose(), header = TRUE)
head(sitios)#cargo la info de los sitios. Archivo"Sitios.txt" o"sitios+alt+DAS"
attach(sitios)

dispersion<-betadisper(pol.dist, group=LU)
permutest(dispersion) 
mod.HSD <- TukeyHSD(dispersion) # no da significativa puedo seguir 
plot(mod.HSD)


#La matriz.de.distancias es la que obtienes con vegdist
#El grupo en tu caso sería la condic
data.ano.pol <- anosim(pol.dist, LU)
summary(data.ano.pol)

plot(data.ano.pol)
#R=0.65. P=0.003

#The ANOSIM statistic compares the mean of ranked dissimilarities between groups to 
#the mean of ranked dissimilarities within groups. An R value close to "1.0" suggests
#dissimilarity between groups while an R value close to "0" suggests an even distribution 
#of high and low ranks within and between groups. R values below "0" suggest that dissimilarities are greater within groups
#than between groups. See Clarke and Gorley (2001) for a guide to interpreting ANOSIM R values.

veg.adonis <- adonis2(pol.dist ~ LU,  contr.unordered = "contr.sum", data = pol)
veg.adonis

#Df SumOfSqs     R2      F Pr(>F)    
#LU        3  1.28071 0.6722 5.4685  0.001 ***
#Residual  8  0.62453 0.3278                  
#Total    11  1.90524 1.0000  


###########NMDS
ord <- metaMDS(pol, distance = "bray")
plot(ord, type = "t")

## Fit environmental variables
ef <- envfit(ord, amb)
ef
#valores de significancia los vectores, en este caso son las condiciones amb
#NMDS1    NMDS2     r2 Pr(>r)   
#Arbol    0.84801 -0.52998 0.1621  0.457   
#Abustos  0.92632 -0.37675 0.7031  0.003 **
#Hierbas -0.91609  0.40096 0.4004  0.108   
#CG      -0.92343  0.38377 0.6392  0.012 * 
#SueloD   0.00436 -0.99999 0.2327  0.323   
#rique    0.55941 -0.82889 0.1110  0.606 

ef.plot<- plot(ef, p.max = 0.05) # para los vectores sig

e.plot<- plot(ef) # para todos 

fit<-envfit(ord, amb)
arrow<-data.frame(fit$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p<-filter(arrow, P<=0.05) 

NMDS1 <- ord$points[,1] ##tiro esto por que no me reconoce los NMDS1 Y 2 en ggplot
NMDS2 <- ord$points[,2]



## solo los vectores significativos
library(ggplot2)
p1 <-ggplot(ef.plot, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=LU),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=LU), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()+
  geom_segment(data=arrow.p, aes(x=0, y=0, xend=NMDS1, yend=NMDS2, lty=FG), arrow=arrow(length=unit(.2, "cm")*arrow.p$R)) 

p1 

### todos los vectores
p<-ggplot(e.plot, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=LU),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=LU), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()+
  geom_segment(data=arrow, aes(x=0, y=0, xend=NMDS1, yend=NMDS2, lty=FG), arrow=arrow(length=unit(.2, "cm"))) 

p

##apiformes composition
api <- read.table(file.choose(), header = TRUE)
head(api) # base de datos api


#Analysis of Similarities
api.dist<-vegdist(api, "bray") # para datos de abundancia se usa bray, la primer celda de la base de datos debe estar VACIA 
api.dist
summary(api.dist)

attach(sitios)

dispersion<-betadisper(api.dist, group=LU)
permutest(dispersion) 
mod.HSD <- TukeyHSD(dispersion) # no da significativa puedo seguir 
plot(mod.HSD)


#La matriz.de.distancias es la que obtienes con vegdist
#El grupo en tu caso sería la condic
data.ano.api <- anosim(api.dist, LU)
summary(data.ano.api)

plot(data.ano.api)
#R=0.61. P=0.001

#The ANOSIM statistic compares the mean of ranked dissimilarities between groups to 
#the mean of ranked dissimilarities within groups. An R value close to "1.0" suggests
#dissimilarity between groups while an R value close to "0" suggests an even distribution 
#of high and low ranks within and between groups. R values below "0" suggest that dissimilarities are greater within groups
#than between groups. See Clarke and Gorley (2001) for a guide to interpreting ANOSIM R values.

veg.adonis <- adonis2(api.dist ~ LU,  contr.unordered = "contr.sum", data = api)
veg.adonis




###########NMDS
ord <- metaMDS(api, distance = "bray")
plot(ord, type = "t")

## Fit environmental variables
ef <- envfit(ord, amb)
ef
#valores de significancia los vectores, en este caso son las condiciones amb
#NMDS1    NMDS2     r2 Pr(>r)   
 
#Arbol   -0.996650  0.081835 0.2797  0.201  
#Abustos -0.987280  0.158968 0.6125  0.014 *
#  Hierbas  0.989160 -0.146815 0.2173  0.293  
#CG       0.983840 -0.179073 0.6125  0.019 *


ef.plot<- plot(ef, p.max = 0.05) # para los vectores sig

e.plot<- plot(ef) # para todos 

fit<-envfit(ord, amb)
arrow<-data.frame(fit$vectors$arrows,R = fit$vectors$r, P = fit$vectors$pvals)
arrow$FG <- rownames(arrow)
arrow.p<-filter(arrow, P<=0.05) 

NMDS1 <- ord$points[,1] ##tiro esto por que no me reconoce los NMDS1 Y 2 en ggplot
NMDS2 <- ord$points[,2]



## solo los vectores significativos
library(ggplot2)
p1 <-ggplot(ef.plot, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=LU),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=LU), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()+
  geom_segment(data=arrow.p, aes(x=0, y=0, xend=NMDS1, yend=NMDS2, lty=FG), arrow=arrow(length=unit(.2, "cm")*arrow.p$R)) 

p1 

### todos los vectores
p<-ggplot(e.plot, aes(NMDS1, NMDS2))+
  geom_point(aes(NMDS1, NMDS2, color=LU),position=position_jitter(.1))+##separates overlapping points
  stat_ellipse(aes(fill=LU), alpha=.2,type='t',size =1, geom="polygon")+ ##changes shading on ellipses
  theme_minimal()+
  geom_segment(data=arrow, aes(x=0, y=0, xend=NMDS1, yend=NMDS2, lty=FG), arrow=arrow(length=unit(.2, "cm"))) 

p


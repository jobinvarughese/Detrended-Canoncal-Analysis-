# Detrended-Canonical-Analysis-

####Detrended Canonical Analysis####

#removing rowsum=0


mat<- data[,c(4,6,8,9,10, 11, 12, 14, 15, 25,37,64)]
mat1<-mat[rowSums(mat[,c(-12,-13,-14,-15)])>0,]

length(mat1$Type2)

decor<-decorana(mat1[,c(-11,-12,-13,-14,-15)], iweigh=0, iresc=4, ira=0, mk=26, short=0,
         before=NULL, after=NULL)

fit <- envfit(decor, mat1[,c(-11)], perm = 999)

plot(decor, choices=c(1,2), origin=TRUE,
     display=c("species"),
     cex = .8, cols = c(1,2), type="text",xlim=c(-2,4))
my_colors<-c("#00AFBB",  "purple", "green", "#E7B800", "#FC4E07")
my_pch <- c(15,16,17,18,19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(mat1$Type2)], col = my_colors[factor(mat1$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(mat1$Type2), cex = 1, pch= my_pch[unique(factor(mat1$Type2))],col = my_colors[unique(factor(mat1$Type2))])



####variables differ####

test<- lmer(Moss~Type1.1+(1|Site), data=data)
summary(test)

test<- lmer(data$Slope~Type1.1+(1|Site), data=data)
summary(test)

test<- aov(data$Moss~data$Type1.1)
summary(test)

test<- aov(data$Tot_bas_ar~data$Type1.1)
summary(test)

test<- aov(data$Tot_ct_tree~data$Type1.1)
summary(test)


test<- aov(data$TCI_2x2~data$Type1.1)
summary(test)

test<- aov(data$Fire~data$Type1.1)
summary(test)

test<- aov(data$Long.Deci..E.~data$Type1.1)
summary(test)

test<- aov(data$Ann.Temp~data$Type1.1)
summary(test)

test<- aov(data$Ann.prec~data$Type1.1)
summary(test)

test<- aov(data$Slope~data$Type1.1)
summary(test)


test<- aov(data$Sine.aspect~data$Type1.1)
summary(test)


test<- aov(data$TWI~data$Type1.1)
summary(test)

####DCA####

mat<- data[,c(5,7,9,10,11, 12, 13, 15, 16, 26, 67, 64, 45,46,48,49,50,63,66,43,44, 60, 65)]#with all imp invasive taxa
#mat<- data[,c(2,3,5,7,9,10,11, 12, 13, 15, 16, 26, 67, 64, 45,46,48,49,50,63,66,43,44, 60, 65)]#with all imp invasive taxa


#mat<- data[,c(5,9,10,11, 12,  15, 64, 45,46,48,49,50,63,66)]#with only the few common taxa


mat1<-mat[rowSums(mat[,c(-12,-13, -14,-15,-16,-17,-18, -19,-20,-21, -22, -23)])>0,]
#mat1<-mat[rowSums(mat[,c(-7,-8,-9,-10,-11,-12,-13,-14)])>0,]
#mat1<-mat[rowSums(mat[,c(-1,-2, -14,-15,-16,-17,-18, -19,-20,-21, -22, -23, -24, -25)])>0,]


length(mat1$Type2)
length(mat1$Type2[mat1$Type2=="Acacia"])
length(mat1$Type2[mat1$Type2=="Pine"])
length(mat1$Type2[mat1$Type2=="Acacia-Eucalyptus"])
length(mat1$Type2[mat1$Type2=="Eucalyptus"])
length(mat1$Type2[mat1$Type2=="Mixed"])

mat2<-mat1[mat1$Type2!="Mixed",]
length(mat2$Type2)

write.csv(mat2, "mat2.csv")
NMDS_dat<-mat2[,c(13:19)]

data_trans<- decostand(NMDS_dat, "standardize")
str(data_trans)



library(factoextra)
my_pca <- prcomp(NMDS_dat, scale = TRUE,
                 center = TRUE, retx = T)
summary(my_pca)

fviz_eig(my_pca)

fviz_pca_ind(my_pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
fviz_pca_var(my_pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
groups <- as.factor(mat2$Type2)
fviz_pca_biplot(my_pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = groups , # Individuals color
                geom = c("point"),
                addEllipses = TRUE, # Concentration ellipses
                ellipse.level=0.95,
                legend.title = "Groups",
)


fviz_pca_ind(my_pca,
             col.ind = groups, # color by groups
             palette = c("#00AFBB",  "purple", "green", "#E7B800", "#FC4E07"),
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "Groups",
             geom = c("point"),
             repel = TRUE
)


####DCA####
mat2<-mat1[mat1$Type2!="Mixed",]
length(mat2$Moss)
decor<-decorana(mat2[,c(-12,-13,-14,-15,-16,-17,-18, -19,-20,-21, -22, -23)], iweigh=0, iresc=4, ira=0, mk=50, short=0,
                before=NULL, after=NULL)


#decor<-decorana(mat1[,c(-14:-7)], iweigh=1, iresc=4, ira=0, mk=26, short=0,
 #               before=NULL, after=NULL)
plot(decor, choices=c(1,2), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type="text",xlim=c(-2,4))
fit <- envfit(decor, mat2[,c(20,21,22,23)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07")
my_pch <- c(15,16,17,18)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(mat2$Type2)], col = my_colors[factor(mat2$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(mat2$Type2), cex = 1, pch= my_pch[unique(factor(mat2$Type2))],col = my_colors[unique(factor(mat2$Type2))])

plot(
  decor,
  type = NULL,
  smooth = FALSE,
  span = 0.2,
  style = c("color", "bw"),
  show_ggplot_code = FALSE,
  ...
)




####For moss and tree count
mat<- data[,c(4,6,8,9,10, 11, 12, 14, 15, 25,37,64, 41,59,65)]
mat1<-mat[rowSums(mat[,c(-12,-13,-14,-15)])>0,]
mat2<-mat1[mat1$Type1.1!="Mixed",]
length(mat2$Moss)
NMDS_dat<-mat2[,c(13:15)]

my_pca <- prcomp(NMDS_dat, scale = TRUE,
                 center = TRUE, retx = T)
groups <- as.factor(mat2$Type1.1)
fviz_pca_biplot(my_pca, repel = TRUE,
                col.var = "#2E9FDF", # Variables color
                col.ind = groups , # Individuals color
                geom = c("point"),
                addEllipses = TRUE, # Concentration ellipses
                ellipse.level=0.95,
                legend.title = "Groups",
)

####DCA with %####
data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar2.csv", head = T, row.names = 1)
mat<- data[,c(42,43,44,46:50,51,53,68,79)]
mat1<-data[rowSums(data[,c(42,43,44,46:50,51,53,68,79)])>0,]
decor<-decorana(data[,c(42:79)], iweigh=0, iresc=1, ira=0, mk=50, short=0,
                before=NULL, after=NULL)

plot(decor, choices=c(1,2), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type="text",xlim=c(-2,4))
fit <- envfit(decor, data[,c(82,83,104)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07", "green")
my_pch <- c(15,16,17,18,19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(data$Type2)], col = my_colors[factor(data$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(data$Type2), cex = 1, pch= my_pch[unique(factor(data$Type2))],col = my_colors[unique(factor(data$Type2))])




mat1<-data[rowSums(data[,c(43,47:50,51,53,54,64,68,79)])>0,]

mat2<-mat1[mat1$Type2 !="Mixed",]

decor<-decorana(mat2[,c(43,47:50,51,53,54,64,68,79)], iweigh=0, iresc=4, ira=0, mk=26, short=0,
                before=NULL, after=NULL)

plot(decor, choices=c(1,3), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type="text",xlim=c(-2,4))
fit <- envfit(decor, mat2[,c(82,83,104)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07", "green")
my_pch <- c(15,16,17,18,19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(mat2$Type2)], col = my_colors[factor(mat2$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(mat2$Type2), cex = 1, pch= my_pch[unique(factor(mat2$Type2))],col = my_colors[unique(factor(mat2$Type2))])

library(ggvegan)
autoplot(decor, display = "species", geom = "text")

autoplot(
  decor,
  axes = c(1, 4),
  geom = c("point", "text"),
  layers = c("species", "sites"),
  legend.position = "right",
  title = NULL,
  subtitle = NULL,
  caption = NULL
)


####DCA_complete datset (192)####

data<-read.csv(file ="~/Desktop/iiser/8 sem/Vegetation paper/invasive_noshola_noothers_grtr_thnzerobasalar.csv", head = T, row.names = 1)
decor<-decorana(data[,c(4:38)], iweigh=0, iresc=4, ira=0, mk=50, short=0,
                before=NULL, after=NULL)


#decor<-decorana(mat1[,c(-14:-7)], iweigh=1, iresc=4, ira=0, mk=26, short=0,
#               before=NULL, after=NULL)
plot(decor, choices=c(2,4), origin=TRUE,
     display=c("species"),
     cex = .6, cols = c(1,2), type=c("points"),xlim=c(-2,4))
fit <- envfit(decor, data[,c(100,101,102)], perm = 999)
plot(fit, cex = .6)
my_colors<-c("#00AFBB",  "purple", "#E7B800", "#FC4E07", "green")
my_pch <- c(15,16,17,18, 19)
points(decor, display = c("sites"),
       choices=1:2,cex = 1, pch=my_pch[factor(data$Type2)], col = my_colors[factor(data$Type2)], xlim=c(-2,4),  origin = TRUE) 
legend(1,-1,unique(data$Type2), cex = 1, pch= my_pch[unique(factor(data$Type2))],col = my_colors[unique(factor(data$Type2))])






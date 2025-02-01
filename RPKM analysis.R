library(Seurat)
library(ggpubr)
library(ggsignif)
library(ggplot2)
setwd("D:/703/wwy/soilvirus/virus anysis/碳水化合物酶rpma")
data<- read.csv('碳水化合物酶rpma-1.csv', row.names = 1, check.names = FALSE)
data<-data[,3:9]
data$group<-factor(data$group, c("Alpine grassland","Degraded grassland","Deserted grassland"))
my_comparisons <- list( c("Alpine grassland", "Degraded grassland"), 
                        c("Degraded grassland", "Deserted grassland"),
                        c("Alpine grassland", "Deserted grassland"))
p1 <- ggplot(data, aes(group, GH, fill=group)) +
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons,
              #step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(color='black',size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  xlab("")+ylab("GH")
p1
p2 <- ggplot(data, aes(group, CBM, fill=group)) +
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  xlab("")+ylab("CBM")
p2
p3 <- ggplot(data, aes(group, GT, fill=group)) +
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  xlab("")+ylab("GT ")
p3
p4 <- ggplot(data, aes(group, CE, fill=group)) +
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  xlab("")+ylab("CE")
p4

p5 <- ggplot(data, aes(group, PLs, fill=group)) +
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  xlab("")+ylab("PLs")
p5

p6 <- ggplot(data, aes(group, Phoh, fill=group)) +
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons,
              step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  xlab("")+ylab("Phoh")
p6


setwd("D:/703/wwy/soilvirus/virus anysis/碳水化合物酶rpma")
library(ggpubr)
library(ggsignif)
library(ggplot2)
data<- read.csv('碳水化合物酶rpma.csv', row.names = 1, check.names = FALSE)
data$group<-factor(data$group, c("Alpine","Degraded","Deserted"))
data$cazy<- factor(data$cazy, c("Glycoside Hydrolases","Glycosyl Transferases", "Carbohydrate Binding Modules", 
                                "Polysaccharide Lyases","Carbohydrate Esterases","phoH"))
my_comparisons <- list( c("Alpine", "Degraded"), 
                        c("Degraded", "Deserted"),
                        c("Alpine", "Deserted"))
p1 <- ggplot(data, aes(group, rpkm, fill=group)) +
  geom_boxplot()+
  geom_signif(comparisons = my_comparisons,
              #step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor=element_blank(),legend.title =element_blank() )+
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=12),
        axis.text.y = element_text(color='black',size=12),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  facet_wrap(.~cazy,scales="free")+
  theme(strip.text=element_text(size = 11, color = "black", face = "bold"))+
  xlab("")+ylab("RPKM")
p1
ggsave('病毒AMG-rpkm.pdf', p1, width = 8, height = 6)
ggsave('病毒AMG-rpkm.png', p1, width = 8, height = 6)


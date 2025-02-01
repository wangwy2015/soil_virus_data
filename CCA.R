library(vegan)
library(ggrepel)
library(ggplot2)
library(ggpubr)

windowsFonts(A=windowsFont("Times New Roman"),
             B=windowsFont("Arial"))

sampledata <- read.csv("vOTU.csv",header = TRUE ,row.names=1,sep = ",",stringsAsFactors = TRUE)
env <- read.csv("土壤性质统计.csv",header = TRUE ,row.names=1,sep = ",",stringsAsFactors = TRUE)
env<-env[-c(1:3,5,8,11,12)]
group <- read.csv("group.csv",header = TRUE ,row.names=1,sep = ",",stringsAsFactors = TRUE)
group<-factor(group$sample, levels = c("Alpine", "Degraded","Deserted"))
sampledata <- t(sampledata)
#对OTU数据进行hellinger转化
sampledata <- decostand(sampledata,method = "hellinger")
#group <- as.list(group)
#定义分组的填充颜色
#col <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
#先进行DCA分析
dca <- decorana(veg = sampledata)
dca1 <- max(dca$rproj[,1])
dca2 <- max(dca$rproj[,2])
dca3 <- max(dca$rproj[,3])
dca4 <- max(dca$rproj[,4])
GL <- data.frame(DCA1 = c(dca1), DCA2 = c(dca2), DCA3 = c(dca3), DCA4 = c(dca4))
GL
rownames(GL) <- c("Gradient length")
write.csv(GL, file = "dca.csv")
#再进行CCA分析  
cca <- cca(sampledata, env, scale = TRUE)
ccascore <- scores(cca)
ccascore$sites
cca$CCA$biplot
ccascore$species

write.csv(ccascore$sites, file = "cca.sample.csv")
write.csv(cca$CCA$biplot, file = "cca.env.csv")
write.csv(ccascore$species, file = "cca.species.csv")

CCAE <- as.data.frame(cca$CCA$biplot[,1:2])

CCAS1 <- ccascore$sites[,1]*0.3
CCAS2 <- ccascore$sites[,2]*0.3

plotdata <- data.frame(rownames(ccascore$sites), CCAS1, CCAS2, group)
colnames(plotdata) <- c("sample","CCAS1","CCAS2","group")

cca1 <- round(cca$CCA$eig[1]/sum(cca$CCA$eig)*100,2)
cca2 <- round(cca$CCA$eig[2]/sum(cca$CCA$eig)*100,2)
#绘制CCA图
p <- ggplot(plotdata, aes(CCAS1, CCAS2)) +
  geom_point(aes(fill = group, color =  group),size = 5) + 
  scale_fill_manual(values = group)+
  scale_color_manual(values=c("#668BAB","#41A2A2","#C0A794"))+
  #stat_chull(geom = "polygon", aes(group =  group, color =  group, fill =  group), alpha = 0.1) +
  xlab(paste("CCA1 ( ",cca1,"%"," )", sep = "")) + 
  ylab(paste("CCA2 ( ",cca2,"%"," )", sep = "")) +
  geom_segment(data = CCAE, aes(x = 0, y = 0, xend = CCAE[,1], yend = CCAE[,2]),
               colour = "black", size = 0.8,
               arrow = arrow(angle = 30, length = unit(0.3, "cm"))) +
  geom_text_repel(data = CCAE, segment.colour = "black",
                  aes(x = CCAE[,1], y = CCAE[,2], 
                      label = rownames(CCAE)),size= 5) +
  geom_vline(aes(xintercept = 0), linetype = "dotted") +
  geom_hline(aes(yintercept = 0), linetype = "dotted") +
  theme(panel.background = element_rect(fill = "white", colour = "black"), 
        panel.grid = element_blank(),
        axis.title = element_text(color = "black", size = 18),
        axis.ticks.length = unit(0.2,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.x = element_text(colour = "black", size = 18),
        axis.title.y = element_text(colour="black", size = 18),
        axis.text = element_text(colour = "black", size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 18), legend.key = element_blank(),
        plot.title = element_text(size = 22, colour = "black", 
                                  face = "bold", hjust = 0.5)) +
  theme(text=element_text(size=20))

p
ggsave('CCA分析.pdf', p, width = 8, height = 6)
ggsave('CCA分析.png', p, width = 8, height = 6)

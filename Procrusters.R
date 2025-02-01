
setwd("D:/703/wwy/soilvirus/virus anysis/Procrustes分析")
library(vegan)

##样方-环境属性矩阵
env <- read.csv("土壤性质统计.csv",row = 1,header = T)

#环境变量的 PCA 需要标准化，详情 ?rda
env_pca <- rda(env, scale = TRUE)

##样方-物种丰度矩阵
otu <- read.csv("碳水化合物酶rpma_1.csv",row = 1,header = T)

otu_hel <- decostand(otu, method = 'hellinger')

#对转化后的物种数据执行 PCA，无需标准化，详情 ?rda
otu_pca <- rda(otu_hel, scale = FALSE)

##排序图比较，以 PCA 的 I 型标尺为例
par(mfrow = c(1, 2))
biplot(env_pca, choices = c(1, 2), scaling = 1, 
       main = '环境组成的PCA', col = c('red', 'blue'))
biplot(otu_pca, choices = c(1, 2), scaling = 1, 
       main = '物种组成的PCA', col = c('red', 'blue'))
#Procrustes 分析
#提取两个 PCA 中的样方排序坐标，均以 I 型标尺为例
site_env <- summary(otu_pca, scaling = 1)$site
site_otu <- summary(otu_pca, scaling = 1)$site

#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(X = env_pca, Y = otu_pca, symmetric = TRUE)
summary(proc)

#旋转图
plot(proc, kind = 1, type = 'text')

#一些重要的结果提取
names(proc)

head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置

#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(X = env_pca, Y = otu_pca, permutations = how(nperm = 999))
prot
res<-residuals(prot)
write.csv(res, "病毒AMG与土壤性质参差-1.csv")
#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量
library(ggplot2)

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.csv("group_1.csv",row = 1,header = T)
Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')

#ggplot2 作图
p <- ggplot(Y) +
  geom_point(aes(X1, X2, color = groups), size = 3, shape = 17) +
  geom_point(aes(PC1, PC2, color = groups), size = 3, shape = 16) +
  scale_color_manual(values = c("#668BAD","#41A2A1","#C0A799"), limits = c('Alpine', 'Degraded', 'Deserted')) +
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2),
               color = 'black',alpha = 0.5, size = 0.5) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.5) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.5) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.5) +
  annotate('text', label = sprintf('M^2 == 0.699'),
           x = -0.18, y = 0.32, size = 5, parse = TRUE) +
  annotate('text', label = c('P == 0.004'),
           x = -0.18, y = 0.28, size = 5, parse = TRUE)+
  theme(axis.ticks=element_line(color="black",size=0.5,lineend = 6),
        axis.text.x=element_text(vjust=0.5,size=15,colour = "black"),
        axis.text.y=element_text(vjust=0.5,size=15,colour = "black"),
        axis.title = element_text(size=15,colour = "black"),
        legend.text = element_text(size=15,colour = "black"),
        legend.background = element_rect(fill=rgb(1,1,1, alpha = 0), colour = NA),
        legend.position = c(0.15,0.8))+
  guides(shape=guide_legend(title = ""))+
  guides(fill=guide_legend(title = ""))+
  guides(color=guide_legend(title = ""))

p
#输出图片
ggsave('procrustes_病毒碳水化合物rpkm-土壤性质-1.pdf', p, width = 6, height = 5)
ggsave('procrustes_病毒碳水化合物rpkm-土壤性质-1.png', p, width = 6, height = 5)

library(Seurat)
library(ggpubr)
library(ggsignif)
library(ggplot2)
data<- read.csv('病毒AMG与土壤性质参差-1-整理.csv', row.names = 1, check.names = FALSE)
data$group<-factor(data$group, c("Alpine","Degraded","Deserted"))
my_comparisons <- list( c("Alpine", "Degraded"), 
                        c("Degraded", "Deserted"),
                        c("Alpine", "Deserted"))
p1 <- ggplot(data, aes(group, value, fill=group)) +
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
  xlab("")+ylab("Residuals")
p1





setwd("D:/703/wwy/soilvirus/virus anysis/Procrustes分析")
library(vegan)

##样方-环境属性矩阵
env <- read.csv("碳水化合物酶rpma_1.csv",row = 1,header = T)
#环境变量的 PCA 需要标准化，详情 ?rda
env_pca <- rda(env, scale = TRUE)

##样方-物种丰度矩阵
otu <- read.csv("bacteria community.csv",row = 1,header = T)
otu<-t(otu)
#物种数据 Hellinger 预转化，详情 ?decostand
otu_hel <- decostand(otu, method = 'hellinger')

#对转化后的物种数据执行 PCA，无需标准化，详情 ?rda
otu_pca <- rda(otu_hel, scale = FALSE)

##排序图比较，以 PCA 的 I 型标尺为例
par(mfrow = c(1, 2))
biplot(env_pca, choices = c(1, 2), scaling = 1, 
       main = '环境组成的PCA', col = c('red', 'blue'))
biplot(otu_pca, choices = c(1, 2), scaling = 1, 
       main = '物种组成的PCA', col = c('red', 'blue'))
#Procrustes 分析
#提取两个 PCA 中的样方排序坐标，均以 I 型标尺为例
site_env <- summary(otu_pca, scaling = 1)$site
site_otu <- summary(otu_pca, scaling = 1)$site

#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(X = env_pca, Y = otu_pca, symmetric = TRUE)
summary(proc)

#旋转图
plot(proc, kind = 1, type = 'text')

#一些重要的结果提取
names(proc)

head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置
#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(X = env_pca, Y = otu_pca, permutations = how(nperm = 999))
prot

#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量
res<-residuals(prot)
write.csv(res, "病毒AMG与细菌组成参差_1.csv")

library(ggplot2)

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.csv("group.csv",row = 1,header = T)
Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')

#ggplot2 作图
p <- ggplot(Y) +
  geom_point(aes(X1, X2, color = groups), size = 3, shape = 16) +
  geom_point(aes(PC1, PC2, color = groups), size = 3, shape = 17) +
  scale_color_manual(values = c("#668BAD","#41A2A1","#C0A799"), limits = c('Alpine', 'Degraded', 'Deserted')) +
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2),
               color = 'black', alpha = 0.5, size = 0.5) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.5) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.5) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.5) +
  annotate('text', label = sprintf('M^2 == 0.729'),
           x = -0.18, y = 0.25, size = 5, parse = TRUE) +
  annotate('text', label = c('P == 0.011'),
           x = -0.18, y = 0.22, size = 5, parse = TRUE)+
  theme(axis.ticks=element_line(color="black",size=0.5,lineend = 6),
        axis.text.x=element_text(vjust=0.5,size=15,colour = "black"),
        axis.text.y=element_text(vjust=0.5,size=15,colour = "black"),
        axis.title = element_text(size=15,colour = "black"),
        legend.text = element_text(size=15,colour = "black"),
        legend.background = element_rect(fill=rgb(1,1,1, alpha = 0), colour = NA),
        legend.position = c(0.15,0.8))+
  guides(shape=guide_legend(title = ""))+
  guides(fill=guide_legend(title = ""))+
  guides(color=guide_legend(title = ""))

p

#输出图片
ggsave('procrustes_病毒AMG与细菌组成_1.pdf', p, width = 6, height = 5)
ggsave('procrustes_病毒AMG与细菌组成_1.png', p, width = 6, height = 5)

data<- read.csv('病毒AMG与细菌组成参差_1_整理.csv', row.names = 1, check.names = FALSE)
data$group<-factor(data$group, c("Alpine","Degraded","Deserted"))
my_comparisons <- list( c("Alpine", "Degraded"), 
                        c("Degraded", "Deserted"),
                        c("Alpine", "Deserted"))
p1 <- ggplot(data, aes(group, value, fill=group)) +
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
  xlab("")+ylab("Residuals")
p1



setwd("D:/703/wwy/soilvirus/virus anysis/Procrustes分析")
library(vegan)

##样方-环境属性矩阵
env <- read.csv("土壤性质统计.csv",row = 1,header = T)

#环境变量的 PCA 需要标准化，详情 ?rda
env_pca <- rda(env, scale = TRUE)

##样方-物种丰度矩阵
otu <- read.csv("碳水化合物酶rpma_去掉不显著碳代谢.csv",row = 1,header = T)

otu_hel <- decostand(otu, method = 'hellinger')

#对转化后的物种数据执行 PCA，无需标准化，详情 ?rda
otu_pca <- rda(otu_hel, scale = FALSE)

##排序图比较，以 PCA 的 I 型标尺为例
par(mfrow = c(1, 2))
biplot(env_pca, choices = c(1, 2), scaling = 1, 
       main = '环境组成的PCA', col = c('red', 'blue'))
biplot(otu_pca, choices = c(1, 2), scaling = 1, 
       main = '物种组成的PCA', col = c('red', 'blue'))
#Procrustes 分析
#提取两个 PCA 中的样方排序坐标，均以 I 型标尺为例
site_env <- summary(otu_pca, scaling = 1)$site
site_otu <- summary(otu_pca, scaling = 1)$site

#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(X = env_pca, Y = otu_pca, symmetric = TRUE)
summary(proc)

#旋转图
plot(proc, kind = 1, type = 'text')

#一些重要的结果提取
names(proc)

head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置

#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(X = env_pca, Y = otu_pca, permutations = how(nperm = 999))
prot
res<-residuals(prot)
write.csv(res, "病毒AMG与土壤性质参差-去掉不显著碳代谢.csv")
#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量
library(ggplot2)

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.csv("group_1.csv",row = 1,header = T)
Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')

#ggplot2 作图
p <- ggplot(Y) +
  geom_point(aes(X1, X2, color = groups), size = 3, shape = 17) +
  geom_point(aes(PC1, PC2, color = groups), size = 3, shape = 16) +
  scale_color_manual(values = c("#668BAD","#41A2A1","#C0A799"), limits = c('Alpine', 'Degraded', 'Deserted')) +
  scale_y_continuous(limits = c(-0.3,0.6))+
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2),
               color = 'black',alpha = 0.5, size = 0.5) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.5) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.5) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.5) +
  annotate('text', label = sprintf('M^2 == 0.79'),
           x = -0.18, y = 0.5, size = 5, parse = TRUE) +
  annotate('text', label = c('p == 0.028'),
           x = -0.18, y = 0.35, size = 5, parse = TRUE)+
  theme(axis.ticks=element_line(color="black",size=0.5,lineend = 6),
        axis.text.x=element_text(vjust=0.5,size=15,colour = "black"),
        axis.text.y=element_text(vjust=0.5,size=15,colour = "black"),
        axis.title = element_text(size=15,colour = "black"),
        legend.text = element_text(size=15,colour = "black"),
        legend.background = element_rect(fill=rgb(1,1,1, alpha = 0), colour = NA),
        legend.position = c(0.15,0.8))+
  guides(shape=guide_legend(title = ""))+
  guides(fill=guide_legend(title = ""))+
  guides(color=guide_legend(title = ""))

p
#输出图片
ggsave('procrustes_病毒碳水化合物rpkm-土壤性质-去掉不显著碳代谢.pdf', p, width = 6, height = 5)
ggsave('procrustes_病毒碳水化合物rpkm-土壤性质-去掉不显著碳代谢.png', p, width = 6, height = 5)

library(Seurat)
library(ggpubr)
library(ggsignif)
library(ggplot2)
data<- read.csv('病毒AMG与土壤性质参差-去掉不显著碳代谢_整理.csv', row.names = 1, check.names = FALSE)
data$group<-factor(data$group, c("Alpine","Degraded","Deserted"))
my_comparisons <- list( c("Alpine", "Degraded"), 
                        c("Degraded", "Deserted"),
                        c("Alpine", "Deserted"))
p1 <- ggplot(data, aes(group, value, fill=group)) +
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons,
              #step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(color='black',size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  
  xlab("")+ylab("Residuals")
p1

ggsave('procrustes_病毒碳水化合物rpkm-土壤性质-去掉不显著碳代谢_参差.pdf', p1, width = 4, height = 4)
ggsave('procrustes_病毒碳水化合物rpkm-土壤性质-去掉不显著碳代谢_参差.png', p1, width = 4, height = 4)



setwd("D:/703/wwy/soilvirus/virus anysis/Procrustes分析")
library(vegan)

##样方-环境属性矩阵
env <- read.csv("碳水化合物酶rpma_去掉不显著碳代谢.csv",row = 1,header = T)
#环境变量的 PCA 需要标准化，详情 ?rda
env_pca <- rda(env, scale = TRUE)

##样方-物种丰度矩阵
otu <- read.csv("bacteria community.csv",row = 1,header = T)
otu<-t(otu)
#物种数据 Hellinger 预转化，详情 ?decostand
otu_hel <- decostand(otu, method = 'hellinger')

#对转化后的物种数据执行 PCA，无需标准化，详情 ?rda
otu_pca <- rda(otu_hel, scale = FALSE)

##排序图比较，以 PCA 的 I 型标尺为例
par(mfrow = c(1, 2))
biplot(env_pca, choices = c(1, 2), scaling = 1, 
       main = '环境组成的PCA', col = c('red', 'blue'))
biplot(otu_pca, choices = c(1, 2), scaling = 1, 
       main = '物种组成的PCA', col = c('red', 'blue'))
#Procrustes 分析
#提取两个 PCA 中的样方排序坐标，均以 I 型标尺为例
site_env <- summary(otu_pca, scaling = 1)$site
site_otu <- summary(otu_pca, scaling = 1)$site

#执行 Procrustes 分析，详情 ?procrustes
#以对称分析为例（symmetric = TRUE）
proc <- procrustes(X = env_pca, Y = otu_pca, symmetric = TRUE)
summary(proc)

#旋转图
plot(proc, kind = 1, type = 'text')

#一些重要的结果提取
names(proc)

head(proc$Yrot)  #Procrustes 分析后 Y 的坐标
head(proc$X)  #Procrustes 分析后 X 的坐标
proc$ss  #偏差平方和 M2 统计量
proc$rotation  #通过该值可获得旋转轴的坐标位置
#PROTEST 检验，详情 ?protest
#以 999 次置换为例
#注：protest() 中执行的是对称 Procrustes 分析，X 和 Y 的分配调换不影响 M2 统计量的计算
set.seed(123)
prot <- protest(X = env_pca, Y = otu_pca, permutations = how(nperm = 999))
prot

#重要统计量的提取
names(prot)
prot$signif  #p 值
prot$ss  #偏差平方和 M2 统计量
res<-residuals(prot)
write.csv(res, "病毒AMG与细菌组成参差_去掉不显著碳代谢.csv")

library(ggplot2)

#提取 Procrustes 分析的坐标
Y <- cbind(data.frame(proc$Yrot), data.frame(proc$X))
X <- data.frame(proc$rotation)

#添加分组信息
group <- read.csv("group.csv",row = 1,header = T)
Y$samples <- rownames(Y)
Y <- merge(Y, group, by = 'samples')

#ggplot2 作图
p <- ggplot(Y) +
  geom_point(aes(X1, X2, color = groups), size = 3, shape = 16) +
  geom_point(aes(PC1, PC2, color = groups), size = 3, shape = 17) +
  scale_color_manual(values = c("#668BAD","#41A2A1","#C0A799"), limits = c('Alpine', 'Degraded', 'Deserted')) +
  scale_y_continuous(limits = c(-0.3,0.55))+
  geom_segment(aes(x = X1, y = X2, xend = PC1, yend = PC2),
               color = 'black', alpha = 0.5, size = 0.5) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  labs(x = 'Dimension 1', y = 'Dimension 2', color = '') +
  geom_vline(xintercept = 0, color = 'gray', linetype = 2, size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', linetype = 2, size = 0.5) +
  geom_abline(intercept = 0, slope = X[1,2]/X[1,1], size = 0.5) +
  geom_abline(intercept = 0, slope = X[2,2]/X[2,1], size = 0.5) +
  annotate('text', label = sprintf('M^2 == 0.746'),
           x = -0.18, y = 0.25, size = 5, parse = TRUE) +
  annotate('text', label = c('P == 0.013'),
           x = -0.18, y = 0.22, size = 5, parse = TRUE)+
  theme(axis.ticks=element_line(color="black",size=0.5,lineend = 6),
        axis.text.x=element_text(vjust=0.5,size=15,colour = "black"),
        axis.text.y=element_text(vjust=0.5,size=15,colour = "black"),
        axis.title = element_text(size=15,colour = "black"),
        legend.text = element_text(size=15,colour = "black"),
        legend.background = element_rect(fill=rgb(1,1,1, alpha = 0), colour = NA),
        legend.position = c(0.15,0.8))+
  guides(shape=guide_legend(title = ""))+
  guides(fill=guide_legend(title = ""))+
  guides(color=guide_legend(title = ""))

p

#输出图片
ggsave('procrustes_病毒AMG与细菌组成_去掉不显著碳代谢.pdf', p, width = 6, height = 5)
ggsave('procrustes_病毒AMG与细菌组成_去掉不显著碳代谢.png', p, width = 6, height = 5)

data<- read.csv('病毒AMG与细菌组成参差_去掉不显著碳代谢_整理.csv', row.names = 1, check.names = FALSE)
data$group<-factor(data$group, c("Alpine","Degraded","Deserted"))
my_comparisons <- list( c("Alpine", "Degraded"), 
                        c("Degraded", "Deserted"),
                        c("Alpine", "Deserted"))
p1 <- ggplot(data, aes(group, value, fill=group)) +
  geom_boxplot() +
  geom_signif(comparisons = my_comparisons,
              #step_increase = 0.3,
              map_signif_level = T, #修改参数map_signif_level=TRUE
              test = wilcox.test)+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'),
        legend.key = element_rect(fill = 'transparent')) +
  theme(legend.position = 'none',
        axis.text.x = element_text(color='black',size=15),
        axis.text.y = element_text(color='black',size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15))+
  scale_fill_manual(values=c("#668BAB","#41A2A2","#C0A794","#F3AC66"))+
  
  xlab("")+ylab("Residuals")
p1

ggsave('procrustes_病毒AMG与细菌组成_去掉不显著碳代谢_参差.pdf', p1, width = 4, height = 4)
ggsave('procrustes_病毒AMG与细菌组成_去掉不显著碳代谢_参差.png', p1, width = 4, height = 4)








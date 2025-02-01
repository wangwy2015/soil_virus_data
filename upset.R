
library(AmesHousing)
library(tidyverse)
library(patchwork)  
library(palmerpenguins)
library(naniar)
library(UpSetR)

outPic="upset_2.pdf" 
pdf(file=outPic,onefile = FALSE,width=9,height=6)
votu <- read.csv('confirmed_virus.csv', row.names = 1, check.names = FALSE)
gg_miss_upset(votu,nsets = length(votu),               #展示多少个数据
              nintersects = 50,                       #展示基因集数目
              #order.by = "freq",                      #按照数目排序
              show.numbers = "yes",                   #柱状图上方是否显示数值
              number.angles = 0,                     #字体角度
              point.size = 4,                         #点的大小
              line.size = 1.0,                        #线条粗细
              text.scale = 3,
              color.pal = "black",
              mainbar.y.label = "vOTUs Numbers")
dev.off()
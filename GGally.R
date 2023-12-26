#install.packages("GGally")
library(GGally)
library(ggplot2)
library(ggprism)
library(ggh4x)
setwd("C:/Users/maihuanzhuo/Desktop")
data(tips, package = "reshape")
#可视化数据分布情况
ggbivariate(tips,outcome = "smoker",explanatory = c("day","time","sex","tip"))+
  scale_fill_brewer(type = "qual")

ggbivariate(tips,outcome = "smoker",legend = 3)+
  scale_fill_brewer(type = "qual")+
  theme_prism()+
  labs(color = "smoker")+
  theme(legend.position = "top",legend.title = element_text(size = 14, face = "bold"))+
  theme(strip.background.x = element_blank(),
    strip.background.y = element_rect(fill = "white",
                                      linetype = "blank"),#blank, solid, dashed, dotted, dotdash, longdash, twodash
    strip.text.x = element_blank(),
    strip.text.y = element_text(angle = 0, size = 12))
ggsave(filename = "ggbivariate1.pdf", width = 12, height = 8, dpi = 600)

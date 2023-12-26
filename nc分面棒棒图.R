#复现NC棒棒图
#《Molecular profiling of aromatase inhibitor sensitive and resistant ER+HER2- postmenopausal breast cancers》
#附件有源数据https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10328947/bin/41467_2023_39613_MOESM9_ESM.zip
library(dplyr)
library(ggplot2)
library(magrittr)
library(tidyr)
library(vroom)
library(forcats)

setwd("C:/Users/maihuanzhuo/Desktop/R包")
#这里用到了vroom包的vroom()函数对数据进行读取，与传统的read.table()函数相比
#其主要有两个优势：（1）读取速度更快；（2）自动判断文件分隔符。
data <- vroom(file = 'Source_data_Fig_3a.csv', col_names = TRUE)#用csv读，xlsx报错
head(data)

#按照原文定义顺序
levels <- c(
  'E2F TARGETS', 'G2M CHECKPOINT', 'MITOTIC SPINDLE', 'MYC TARGETS V1', 'MYC TARGETS V2', 'P53 PATHWAY',
  'ESTROGEN RESPONSE EARLY', 'ESTROGEN RESPONSE LATE', 'IL2 STAT5 SIGNALING', 'KRAS SIGNALING DN', 'KRAS SIGNALING UP', 'MTORC1 SIGNALING', 'TNFA SIGNALING VIA NFKB',
  'ALLOGRAFT REJECTION', 'COMPLEMENT', 'IL6 JAK STAT3 SIGNALING', 'INFLAMMATORY RESPONSE', 'INTERFERON ALPHA RESPONSE', 'INTERFERON GAMMA RESPONSE',
  'HYPOXIA',
  'EPITHELIAL MESENCHYMAL TRANSITION',
  'GLYCOLYSIS',
  'APICAL JUNCTION'
)

#吸取颜色
colors <- c(
  rep('#279D77', 6),
  rep('#CF6611', 7),
  rep('#7974A1', 6),
  '#E2348E',
  '#E5B63D',
  '#91793C',
  '#686868'
)

data <- data %>% 
  #将以"NES"开头的列变成长格式，其中列名保存在新的"Compare"列中，列值保存在"NES"列中
  pivot_longer(cols = starts_with('NES'), names_to = 'Compare', values_to = 'NES') %>% 
  pivot_longer(cols = starts_with('GeneRatio'), names_to = 'GeneRatio_Group', values_to = 'Ratio') %>% 
  pivot_longer(cols = starts_with('p.adjust'), names_to = 'p.adjust_Group', values_to = 'p.adjust') %>% 
  #使用gsub()函数将"Description"列中的 "HALLMARK_" 替换为空字符串
  mutate(Description = gsub(pattern = 'HALLMARK_', replacement = '', x = Description)) %>% 
  #将 "Description" 列中的下划线 "_" 替换为空格 " "
  mutate(Description = gsub(pattern = '_', replacement = ' ', x = Description)) %>%
  #保留 "Description" 列中包含在 "levels" 中的值的行，前面level按照原文定义了
  filter(Description %in% levels) %>% 
  #重新排序 "Description" 列的因子水平，以反转它们的顺序
  mutate(Description = fct_relevel(Description, rev(levels))) %>% 
  #将 "Compare" 列中的 "NES " 替换为空字符串
  mutate(Compare = gsub(pattern = 'NES ', replacement = '', x = Compare)) %>% 
  #forcats包对列进行因子化fct_relevel并匹配给定的顺序
  mutate(Compare = fct_relevel(Compare, c('PRs vs GRs', 'PRs ESR1 HIGH vs GRs', 
                                          'PRs ESR1 LOW vs GRs', 'PRs ESR1 LOW vs PRs ESR1 HIGH')))

p1 <- data %>% 
  ggplot(aes(x = NES, y = Description)) +
  geom_vline(xintercept = 0, color = 'grey', linewidth = 1) +#添加竖线
  geom_segment(aes(x = 0, xend = NES, y = Description, yend = Description), 
               color = 'grey', linewidth = 1) +#添加横线，即棒棒糖图
  geom_point(aes(color = -log10(p.adjust), size = Ratio)) +
  scale_color_gradient(low = 'blue', high = 'red', breaks = seq(3, 9, 2), limits = c(1, 9), 
                       name = 'Significance\n(-log10 FDR)') +
  scale_size_continuous(name = 'GeneRatio(%)', range = c(2, 5), limits = c(30, 60)) +
  scale_x_continuous(limits = c(-3, 3), breaks = seq(-3, 3, 1), expand = c(0, 0)) +
  labs(x = 'Normalized Enrichment Score (NES)', y = NULL) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(linewidth = 1.5, color = 'black'),
    panel.spacing.x = unit(0.15, units = 'in'),
    strip.background = element_rect(linewidth = 1.5, color = 'black'),
    axis.ticks = element_line(color = 'black'),
    axis.text.x = element_text(family = 'sans', color = 'black', size = 11),
    axis.text.y = element_text(family = 'sans', size = 11, face = 'bold', color = rev(colors)),#用geom_text代替试试.就不会出现警告
    axis.title.x = element_text(family = 'sans', color = 'black', size = 14),
    legend.title = element_text(family = 'sans', color = 'black', face = 'bold', size = 12),
    legend.text = element_text(family = 'sans', color = 'black', face = 'bold', size = 11)
  ) +
  facet_wrap(~Compare, ncol = 4)
p1

rect.data <- data.frame(
  ymin = c(0, 1.5, 2.5, 3.5, 4.5, 10.5, 17.5),
  ymax = c(1.5, 2.5, 3.5, 4.5, 10.5, 17.5, 23.5),
  colors = letters[1:7]
)
head(rect.data)

p2 <- p1 + 
  geom_rect(data = rect.data, aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = colors), 
            inherit.aes = FALSE, alpha = 0.2, show.legend = FALSE) +
  scale_fill_manual(values = c('#686868', '#91793C', '#E5B63D', '#E2348E', '#7974A1', '#CF6611', '#279D77')) +
  guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))#调整图例顺序
p2
ggsave(filename = "nc分面棒棒图.pdf", width = 14, height = 8, dpi = 600)

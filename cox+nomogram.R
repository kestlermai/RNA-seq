library(rms)
#install.packages("survminer")
library(survival)
library(survminer)
#定义因子
lung$sex <- factor(lung$sex,levels = c(1,2),labels = c("male","female"))
dd <- datadist(lung)#用datadist函数创建datadist对象并把lung数据赋值给data_lung
options(datadist="dd")#options将lung的赋值用到后续模型中

#构建cox比例风险回归模型
#cph(formula, data, subset, weights, na.action, cluster, x = FALSE, y = TRUE, surv = TRUE, time.inc, ...)
#formula：一个公式，指定了生存时间和事件的关系，通常是形如 Surv(time, status) ~ X1 + X2 + ... + Xn 的形式，
#其中 time 是生存时间向量，status 是事件指示符向量，X1, X2, ..., Xn 是预测变量。
#x：一个逻辑值，指定是否返回模型系数。
#y：一个逻辑值，指定是否返回响应变量。
#surv：一个逻辑值，指定是否返回生存数据。

#ph.ecog：ECOG成绩得分（0 =好5 =死）
#ph.karno：医师对Karnofsky成绩评分（差= 0-良好= 100）
#pat.karno：Karnofsky表现评分，按患者评分
cox <- cph(formula=Surv(time,status)~age+sex+ph.ecog+ph.karno,
           data=lung,
           x=T,y=T,surv=T)
#查看cox模型
summary(cox)
#PH检验,基于Schoenfeld残差
res_cox <- cox.zph(cox)
res_cox
plot(res_cox)
ggcoxzph(res_cox)
ggcoxdiagnostics(cox, type = "schoenfeld")
#森林图
fit.cox <- coxph(Surv(time, status) ~ . , data = lung)
ggforest(fit.cox, data = lung,
         main = "Hazard ratio",
         cpositions = c(0.01, 0.15, 0.35), # 更改前三列的相对位置
         fontsize = 0.7,
         refLabel = "reference",
         noDigits = 2)

#若变量P<0.05，说明该变量的风险值具有时间依赖性，应该将该变量进行分层分析
cox_str <- cph(formula=Surv(time,status)~age+sex+ph.ecog+strat(ph.karno),
              data=lung,
              x=T,y=T,surv=T)
  
res_cox_str <- cox.zph(cox_str)
res_cox_str
plot(res_cox_str)
summary(cox_str)#分层变量不参与模型拟合，所以结果不会出现分层变量的系数

# 计算 cumulative incidence 和 competing risk
times <- seq(0, max(lung$time), by = 365.25)
cumulative_incidence <- survfit(cox_str, newdata = lung, 
                                times = times, type = "restricted")
competing_risk <- survfit(cox_str, newdata = lung, 
                          times = seq(0, max(lung$time), by = 365.25), type = "cumulative")

# 将结果合并为一个数据框
cumulative_incidence_df <- data.frame(cumulative_incidence$time, cumulative_incidence$surv, competing_risk$surv)
colnames(cumulative_incidence_df) <- c("Time", "Cumulative_Incidence", "Competing_Risk")

ggplot(cumulative_incidence_df, aes(Time)) + 
  geom_line(aes(y = Cumulative_Incidence, color = "Cumulative Incidence")) +
  geom_line(aes(y = Competing_Risk, color = "Competing Risk")) +
  scale_color_manual(name = "Legend", values = c("blue", "red"), labels = c("Cumulative Incidence", "Competing Risk")) +
  labs(x = "Time", y = "Survival Probability", title = "Cumulative Incidence Plot")

# 计算 cumulative incidence 和 competing risk
cumulative_incidence <- survfit(cox_str, newdata = lung, 
                                type = "cumulative incidence", conf.type = "none")
competing_risk <- survfit(cox_str, newdata = lung, 
                          type = "competing risks", conf.type = "none")

# 绘制 cumulative incidence plot
plot(cuminc(cumulative_incidence, competing_risk), xlab = "Time", ylab = "Cumulative incidence", 
     main = "Cumulative incidence plot", lty = 1:2, col = 1:2, lwd = 2)
legend("topright", legend = c("Cumulative incidence", "Competing risk"), lty = 1:2, col = 1:2, lwd = 2)

#Survival函数是预测模型的生存函数
survival_pred <- Survival(cox_str)
survival_pred_1year <- function(x)survival_pred(365,x)
survival_pred_2year <- function(x)survival_pred(730,x)

nom <- nomogram(cox_str,fun=list(survival_pred_1year,survival_pred_2year),#fun输入自定义函数列表
                fun.at=c(0.05,seq(0.1,0.9,by=0.05),0.95),#at每个函数的轴刻度值
                funlabel = c("1 year survial probability","2 year survial probability"))

plot(nom, 
     lplabel="Linear Predictor",
     xfrac = 0.2, # 左侧标签距离坐标轴的距离
     #varname.label = TRUE, 
     tcl = -0.2, # 刻度长短和方向 
     lmgp = 0.1, # 坐标轴标签距离坐标轴远近
     points.label ='Points', 
     total.points.label = 'Total Points',
     cap.labels = FALSE,
     cex.var = 1, # 左侧标签字体大小
     cex.axis = 1, # 坐标轴字体大小
     col.grid = gray(c(0.8, 0.95))) # 竖线颜色


library(regplot)
#install.packages("regplot")
regplot(cox,
        plots = c("violin", "boxes"), ##连续性变量形状，可选"no plot" "density" "boxes" "ecdf" "bars" "boxplot" "violin" "bean" "spikes"；分类变量的形状，可选"no plot" "boxes" "bars" "spikes"
        observation = lung[1,], #用哪行观测，或者T F
        center = T, # 对齐变量
        subticks = T,
        droplines = T,#是否画竖线
        title = "nomogram",
        points = T, # 截距项显示为0-100
        odds = T, # 是否显示OR值
        showP = T, # 是否显示变量的显著性标记
        rank = "sd", # 根据sd给变量排序
        interval="confidence", # 展示可信区间
        clickable = F) # 是否可以交互

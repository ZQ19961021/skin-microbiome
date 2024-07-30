data1 = read.table("age_data.txt",sep = "\t",check.names = F,
                   stringsAsFactors = F,header = T,row.names = 1)
data1 = data1[!is.na(data1$skin_age),]
data1 = data1[!is.na(data1$age),]

library(ggplot2)

a = lm(data1$skin_age~data1$age)
summary(a)
ggplot(data1,aes(x = data1$age,y = data1$skin_age))+geom_point()+
  geom_smooth(method = "lm")+theme_classic()

data1$point_line = (0.63838*data1$age + 18.99732-data1$skin_age)/sqrt(0.63838^2+1^2)
data1$type = "control"
data1$type[data1$point_line > 3.639] = "Elderly"
data1$type[data1$point_line < -3.639] = "Young"
ggplot(data1,aes(x = data1$age,y = data1$skin_age))+geom_point(aes(color = type))+
  geom_smooth(method = "lm")+theme_classic()

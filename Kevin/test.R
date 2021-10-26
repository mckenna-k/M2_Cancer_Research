library(tidyverse)
library(survival)
data <- read.csv("~/M2SSD/ProjetTutore/M2_Cancer_Research/data/TCGA-CDR_data.csv")
head(data[4])
data = data %>% na_if('#N/A')
data = data %>% mutate_at(c(2,4,13,16,17,22,26:33),as.character)
data = data %>% mutate_at(c(4,13,16,17,22,26:33),as.numeric)

data2 = data %>% mutate('site' = str_sub(data$bcr_patient_barcode,6,7))
table(data2$site)
data_reduc = data[c(1,3:6,9,11,13,16,17,26:33)]


ggplot(data)+
  geom_bar(aes(type))
ggplot(data)+
  geom_bar(aes(race))
ggplot(data)+
  geom_bar(aes(gender))
ggplot(data)+
  geom_bar(aes(age_at_initial_pathologic_diagnosis))

ggplot(data)+
  geom_density(aes(as.numeric(OS.time)))+
  geom_density(aes(as.numeric(DFI.time)), col = 'red')+
  facet_wrap(data$type)
ggplot(data)+
  geom_density(aes(as.numeric(OS.time)))+
  geom_density(aes(as.numeric(DFI.time)), col = 'red')+
  geom_density(aes(as.numeric(PFI.time)), col = 'blue')+
  facet_wrap(data$ajcc_pathologic_tumor_stage, scales = 'free')
ggplot(data)+
  geom_density(aes(as.numeric(OS.time)))+
  facet_wrap(data$gender)
ggplot(data)+
  geom_density(aes(as.numeric(OS.time)))+
  facet_wrap(data$race)
ggplot(data)+
  geom_density(aes(as.numeric(OS.time)))+
  facet_grid(data$race~data$type, scales = 'free')
ggplot(data)+
  geom_density(aes(as.numeric(OS.time)))+
  geom_density(aes(as.numeric(DSS.time)))
  facet_grid(data$gender~data$type, scales = 'free')
# 
# data = data %>% na_if('#N/A')
# data = data %>% mutate_at(c(4,13,16,17,22,26:33),as.character)
# data = data %>% mutate_at(c(4,13,16,17,22,26:33),as.numeric)
# 
# data_reduc = data[c(1,3:6,9,11,13,16,17,26:33)]

types = factor(data_reduc$type)
data_arranged = arrange(data_reduc, OS.time)
data_arranged = data_arranged %>% mutate_at(11,as.numeric)

#Overall
y = Surv(data_arranged$OS.time, event = data_arranged$OS)
test <- survfit(y~1, conf.type = 'plain')
plot(test)
test
#BRCA CANCER
y = Surv(data_arranged$OS.time[data_arranged$type=='BRCA'], event = data_arranged$OS[data_arranged$type=='BRCA'])
test <- survfit(y~1, conf.type = 'plain')
plot(test)
test
# All types
y = Surv(data_arranged$OS.time, event = data_arranged$OS)
test <- survfit(y~data_arranged$type, conf.type = 'plain')
plot(test, col= 1:33)
text(test)
test
# Overall by race
y = Surv(data_arranged$OS.time, event = data_arranged$OS)
test <- survfit(y~data_arranged$race, conf.type = 'plain')
plot(test, col= 1:8)
# legend('topright',legend = unique(factor(data$race)), fill = c(1:8))
test
# Overall by gender
y = Surv(data_arranged$OS.time, event = data_arranged$OS)
test <- survfit(y~data_arranged$gender, conf.type = 'plain')
plot(test, col= unique(factor(data_arranged$gender)))
legend('topright',legend = unique(factor(data$gender)), fill = c(1:2))
test




# Cancers to look at 
# HNSC, CESC, UCEC, LIHC, KIRP, LGG, 
# White people only? - MESO
# High account in black people?  CHOL, ESCA, PRAD, 

tab2 = as.data.frame(tab)






















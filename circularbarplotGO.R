######################################################################
#Author: Natalia Faraj Murad
#This code creates a circular barplot representation of Gene
#Ontology data from topGO enrichment analysis. This code was
#adquired from https://www.r-graph-gallery.com/circular-barplot.html
#and adapted for our data. This analysis is available was published on
#Silva-Brandao et al (2020) Transcriptome differential co-expression 
#reveals distinct molecular response of fall-armyworm strains to DIMBOA.
#DOI: 10.1002/ps.6051 Journal: Pest Management Science
######################################################################

setwd("C:/Users/natalia/Desktop/Karina/up&down")

#Loading packages
library(tidyverse)
library(extrafont)
library(gridExtra)


###############EXP1
exp1 <- read.csv("exp1.csv", h = T, sep = "\t", dec = ".", stringsAsFactors = FALSE)

vec1 <- exp1[,2]
vec2 <- exp1[1:23,4]

individual <- vector(length=54)
individual[1:31] = vec1
individual[32:54]=vec2

group = vector(length = 54)
group[1:31] = "UP"
group[32:54]= "DOWN"

vec3 <- exp1[,1]
vec4 <- exp1[1:23,3] 

value <- vector(length=54)
value[1:31] = vec3
value[32:54]= vec4
value <- abs(value)

data <- data.frame(
  individual,
  group,
  value
)

#Ordering

data = data %>% arrange(group, value)

#Setting the number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
colnames(to_add) = colnames(data)
to_add$group=rep(levels(data$group), each=empty_bar)

data=rbind(data, to_add)
data=data %>% arrange(group)
data$id=seq(1, nrow(data))

#Getting the name and the y position of each label
label_data=data
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     #Substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

#Preparing a data frame for base lines
base_data=data %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

#Preparing a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]
#pdf("exp1.pdf")

#Making the plot
p = ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
    geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
#  scale_fill_manual(values=c("#BFE476FF", "#C0BA99")) +
    scale_fill_manual(values=c("#d8b365", "#5ab4ac")) +
  
#Adding a val=100/75/50/25 lines. Doing it at the beginning to make sur barplots are OVER it.
geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
#Adding text showing the value of each 100/75/50/25 lines
annotate("text", x = rep(max(data$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
ylim(-100,120) +
theme_minimal() +
theme(
  legend.position = "none",
  axis.text = element_blank(),
  axis.title = element_blank(),
  panel.grid = element_blank(),
  plot.margin = unit(rep(-1,4), "cm") 
) +
coord_polar() + 
geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), colour="black", family="Times New Roman", size=3.5, angle= label_data$angle, inherit.aes = FALSE ) +
  
#Adding base line information
geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
annotate("text", x = 15, y = -20, label="DOWN", angle = 270, fontface="bold", family="Times New Roman") +
annotate("text", x = 44, y = -12, label="UP", angle = 90, fontface="bold", family="Times New Roman")
# geom_text(data=base_data, aes(x = title, y = -18, label=group), position = position_stack(vjust = -1), hjust = c(0, 1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
#?geom_text
p
ggsave("exp1.tiff", units="in", width=7, height=7, dpi=600, compression = 'lzw')

dev.off()
###################EXP2
exp2 <- read.csv("exp2.csv", h = T, sep = "\t", dec = ".", stringsAsFactors = FALSE)

vec1 <- exp2[,2]
vec2 <- exp2[1:24,4]

individual <- vector(length=50)
individual[1:26] = vec1
individual[27:50]=vec2

group = vector(length = 50)
group[1:26] = "UP"
group[27:50]= "DOWN"

vec3 <- exp2[,1]
vec4 <- exp2[1:24,3] 
value <- vector(length=50)
value[1:26] = vec3
value[27:50]= vec4

value <- abs(value)

data1 <- data.frame(
  individual,
  group,
  value
)

#Ordering
data1 = data1 %>% arrange(group, value)

#Setting a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(data1$group), ncol(data1)) )
colnames(to_add) = colnames(data1)
to_add$group=rep(levels(data1$group), each=empty_bar)
data1=rbind(data1, to_add)
data1=data1 %>% arrange(group)
data1$id=seq(1, nrow(data1))

#Getting the name and the y position of each label
label_data=data1
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     #Substracting 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

#Preparing a data frame for base lines
base_data=data1 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

#Preparing a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

#Making the plot
p1 = ggplot(data1, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  #  scale_fill_manual(values=c("#BFE476FF", "#C0BA99")) +
  scale_fill_manual(values=c("#d8b365", "#5ab4ac")) +
  
#Adding a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
#Adding text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data1$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), colour="black", family="Times New Roman", size=3.7, angle= label_data$angle, inherit.aes = FALSE ) +
  
#Adding base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  annotate("text", x = 15, y = -18, label="DOWN", angle = 270, fontface="bold", family="Times New Roman") +
  annotate("text", x = 44, y = -12, label="UP", angle = 90, fontface="bold", family="Times New Roman")
# geom_text(data=base_data, aes(x = title, y = -18, label=group), position = position_stack(vjust = -1), hjust = c(0, 1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
#?geom_text
p1

ggsave("exp2.tiff", units="in", width=7, height=7, dpi=600, compression = 'lzw')
#grid.arrange(p1 , p1 , ncol=2)

###################EXP3
exp3 <- read.csv("exp3.csv", h = T, sep = "\t", dec = ".", stringsAsFactors = FALSE)

vec1 <- exp3[,2]
vec2 <- exp3[1:20,4]

individual <- vector(length=59)
individual[1:39] = vec1
individual[40:59]=vec2

group = vector(length = 59)
group[1:39] = "UP"
group[40:59]= "DOWN"

vec3 <- exp3[,1]
vec4 <- exp3[1:20,3] 
value <- vector(length=59)
value[1:39] = vec3
value[40:59]= vec4

value <- abs(value)

data1 <- data.frame(
  individual,
  group,
  value
)


#Ordering
data1 = data1 %>% arrange(group, value)

#Setting a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(data1$group), ncol(data1)) )
colnames(to_add) = colnames(data1)
to_add$group=rep(levels(data1$group), each=empty_bar)
data1=rbind(data1, to_add)
data1=data1 %>% arrange(group)
data1$id=seq(1, nrow(data1))

#Getting the name and the y position of each label
label_data=data1
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     #Subtracting 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

#Preparing a data frame for base lines
base_data=data1 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

#Preparing a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

#Making the plot
p2 = ggplot(data1, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  #  scale_fill_manual(values=c("#BFE476FF", "#C0BA99")) +
  scale_fill_manual(values=c("#d8b365", "#5ab4ac")) +
  
  # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
#Adding text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data1$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), colour="black", family="Times New Roman", size=3.7, angle= label_data$angle, inherit.aes = FALSE ) +
  
 #Adding base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  annotate("text", x = 15, y = -18, label="DOWN", angle = 270, fontface="bold", family="Times New Roman") +
  annotate("text", x = 50, y = -12, label="UP", angle = 90, fontface="bold", family="Times New Roman")
# geom_text(data=base_data, aes(x = title, y = -18, label=group), position = position_stack(vjust = -1), hjust = c(0, 1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
#?geom_text
p2
ggsave("exp3.tiff", units="in", width=8, height=7, dpi=600, compression = 'lzw')
###################EXP4
exp4 <- read.csv("exp4.csv", h = T, sep = "\t", dec = ".", stringsAsFactors = FALSE)

vec1 <- exp4[1:27,2]
vec2 <- exp4[,4]

individual <- vector(length=62)
individual[1:27] = vec1
individual[27:62]=vec2

group = vector(length = 62)
group[1:27] = "UP"
group[28:62]= "DOWN"

vec3 <- exp4[1:27,1]
vec4 <- exp4[,3] 
value <- vector(length=62)
value[1:27] = vec3
value[27:62]= vec4

value <- abs(value)

data1 <- data.frame(
  individual,
  group,
  value
)


#Ordering
data1 = data1 %>% arrange(group, value)

#Setting a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(data1$group), ncol(data1)) )
colnames(to_add) = colnames(data1)
to_add$group=rep(levels(data1$group), each=empty_bar)
data1=rbind(data1, to_add)
data1=data1 %>% arrange(group)
data1$id=seq(1, nrow(data1))

#Getting the name and the y position of each label
label_data=data1
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     #Subtracting 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

#Preparing a data frame for base lines
base_data=data1 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

#Preparing a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

#Making the plot
p3 = ggplot(data1, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  #  scale_fill_manual(values=c("#BFE476FF", "#C0BA99")) +
  scale_fill_manual(values=c("#d8b365", "#5ab4ac")) +
  
#Adding a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
#Adding text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data1$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), colour="black", family="Times New Roman", size=3.7, angle= label_data$angle, inherit.aes = FALSE ) +
  
#Adding base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  annotate("text", x = 17, y = -24, label="DOWN", angle = 270, fontface="bold", family="Times New Roman") +
  annotate("text", x = 50, y = -18, label="UP", angle = 90, fontface="bold", family="Times New Roman")
# geom_text(data=base_data, aes(x = title, y = -18, label=group), position = position_stack(vjust = -1), hjust = c(0, 1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
#?geom_text
p3
ggsave("exp4.tiff", units="in", width=8, height=8, dpi=600, compression = 'lzw')
#grid.arrange(p , p1 , ncol=2)

###################EXP5
exp5 <- read.csv("exp5.csv", h = T, sep = "\t", dec = ".", stringsAsFactors = FALSE)

vec1 <- exp5[,2]
vec2 <- exp5[1:31,4]

individual <- vector(length=63)
individual[1:32] = vec1
individual[33:63]=vec2

group = vector(length = 63)
group[1:32] = "UP"
group[33:63]= "DOWN"

vec3 <- exp5[,1]
vec4 <- exp5[1:31,3] 
value <- vector(length=63)
value[1:32] = vec3
value[33:63]= vec4

value <- abs(value)

data1 <- data.frame(
  individual,
  group,
  value
)


#Ordering
data1 = data1 %>% arrange(group, value)

#Setting a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(data1$group), ncol(data1)) )
colnames(to_add) = colnames(data1)
to_add$group=rep(levels(data1$group), each=empty_bar)
data1=rbind(data1, to_add)
data1=data1 %>% arrange(group)
data1$id=seq(1, nrow(data1))

#Getting the name and the y position of each label
label_data=data1
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     #Subtracting 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

#Preparing a data frame for base lines
base_data=data1 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

#Preparing a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

#Making the plot
p4 = ggplot(data1, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  #  scale_fill_manual(values=c("#BFE476FF", "#C0BA99")) +
  scale_fill_manual(values=c("#d8b365", "#5ab4ac")) +
  
#Adding a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
#Adding text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data1$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), colour="black", family="Times New Roman", size=3.7, angle= label_data$angle, inherit.aes = FALSE ) +
  
#Adding base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  annotate("text", x = 18, y = -18, label="DOWN", angle = 270, fontface="bold", family="Times New Roman") +
  annotate("text", x = 53, y = -12, label="UP", angle = 90, fontface="bold", family="Times New Roman")
# geom_text(data=base_data, aes(x = title, y = -18, label=group), position = position_stack(vjust = -1), hjust = c(0, 1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
#?geom_text
p4
ggsave("exp5.tiff", units="in", width=8, height=8, dpi=600, compression = 'lzw')
###################EXP6
exp6 <- read.csv("exp6.csv", h = T, sep = "\t", dec = ".", stringsAsFactors = FALSE, row.names=NULL)

vec1 <- exp6[1:23,2]
vec2 <- exp6[,4]

individual <- vector(length=61)
individual[1:23] = vec1
individual[24:61]=vec2

group = vector(length = 61)
group[1:23] = "UP"
group[24:61]= "DOWN"

vec3 <- exp6[1:23,1]
vec4 <- exp6[,3] 
value <- vector(length=61)
value[1:23] = vec3
value[24:61]= vec4

value <- as.numeric(value)
value <- abs(value)

data1 <- data.frame(
  individual,
  group,
  value
)

#Ordering
data1 = data1 %>% arrange(group, value)

#Setting a number of 'empty bar' to add at the end of each group
empty_bar=3
to_add = data.frame( matrix(NA, empty_bar*nlevels(data1$group), ncol(data1)) )
colnames(to_add) = colnames(data1)
to_add$group=rep(levels(data1$group), each=empty_bar)
data1=rbind(data1, to_add)
data1=data1 %>% arrange(group)
data1$id=seq(1, nrow(data1))

#Getting the name and the y position of each label
label_data=data1
number_of_bar=nrow(label_data)
angle= 90 - 360 * (label_data$id-0.5) /number_of_bar     #Subtracting 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
label_data$hjust<-ifelse( angle < -90, 1, 0)
label_data$angle<-ifelse(angle < -90, angle+180, angle)

#Preparing a data frame for base lines
base_data=data1 %>% 
  group_by(group) %>% 
  summarize(start=min(id), end=max(id) - empty_bar) %>% 
  rowwise() %>% 
  mutate(title=mean(c(start, end)))

#Preparing a data frame for grid (scales)
grid_data = base_data
grid_data$end = grid_data$end[ c( nrow(grid_data), 1:nrow(grid_data)-1)] + 1
grid_data$start = grid_data$start - 1
grid_data=grid_data[-1,]

#Making the plot
p5 = ggplot(data1, aes(x=as.factor(id), y=value, fill=group)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  #  scale_fill_manual(values=c("#BFE476FF", "#C0BA99")) +
  scale_fill_manual(values=c("#d8b365", "#5ab4ac")) +
  
#Adding a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
  geom_segment(data=grid_data, aes(x = end, y = 80, xend = start, yend = 80), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 60, xend = start, yend = 60), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 40, xend = start, yend = 40), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  geom_segment(data=grid_data, aes(x = end, y = 20, xend = start, yend = 20), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
  
#Adding text showing the value of each 100/75/50/25 lines
  annotate("text", x = rep(max(data1$id),4), y = c(20, 40, 60, 80), label = c("20", "40", "60", "80") , color="grey", size=3 , angle=0, fontface="bold", hjust=1) +
  
  geom_bar(aes(x=as.factor(id), y=value, fill=group), stat="identity", alpha=0.5) +
  ylim(-100,120) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    plot.margin = unit(rep(-1,4), "cm") 
  ) +
  coord_polar() + 
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), colour="black", family="Times New Roman", size=3.7, angle= label_data$angle, inherit.aes = FALSE ) +
  
#Adding base line information
  geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), colour = "black", alpha=0.8, size=0.6 , inherit.aes = FALSE )  +
  annotate("text", x = 18, y = -18, label="DOWN", angle = 270, fontface="bold", family="Times New Roman") +
  annotate("text", x = 50, y = -12, label="UP", angle = 90, fontface="bold", family="Times New Roman")
# geom_text(data=base_data, aes(x = title, y = -18, label=group), position = position_stack(vjust = -1), hjust = c(0, 1), colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
#?geom_text
p5
ggsave("exp6.tiff", units="in", width=8, height=8, dpi=600, compression = 'lzw')


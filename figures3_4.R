
library(ggplot2)
library(patchwork)
library(doubletrouble)
library(ggpubr)

# ============================================================================
#this section creates the figures 4 G and 4 H and calculates their p values 
# ============================================================================


forpaired = read.delim("pairedtest.txt")

aggregate(Normalized_events_MYA ~ Condition, data = forpaired, FUN = mean)
aggregate(Normalized_events_length ~ Condition, data = forpaired, FUN = mean)

duplicated = subset(forpaired, forpaired$Condition == "Yes Event")[,3]
control = subset(forpaired, forpaired$Condition == "No Event")[,3]

length_test <- wilcox.test(duplicated, control, paired = TRUE, conf.int = TRUE)
length_test


duplicated = subset(forpaired, forpaired$Condition == "Yes Event")[,4]
control = subset(forpaired, forpaired$Condition == "No Event")[,4]

mya_test <- wilcox.test(duplicated, control, paired = TRUE, conf.int = TRUE)
mya_test



ggplot(forpaired, aes(x=factor(Condition), y=Normalized_events_length, fill = factor(Condition)))+
  geom_boxplot(notch = FALSE, outlier.shape = 3, coef = 0, linewidth = 1.4)+
  theme_classic()+
  theme(axis.text=element_text(size=20, face="bold"),
        axis.title=element_text(size=20, face="bold"),
        legend.text = element_text(size = 20, face="bold"),
        legend.title = element_blank(),
        legend.position = "none",
        title = element_text(size=20, face="bold"))+
  stat_summary(geom="point", 
               shape=20, size=5, color="black", fill="black",fun = "mean") +
  
  xlab("")+
  ylab("Events/Branch lenght")+
  scale_fill_manual(values= c("#DADAEB", "#6A51A3"))+
  geom_signif(comparisons = list(c("Yes Event", "No Event")),
              map_signif_level = TRUE, size = 1,
              textsize = 10)+
  ylim(0,60)


ggplot(forpaired, aes(x=factor(Condition), y=Normalized_events_MYA, fill = factor(Condition)))+
  geom_boxplot(notch = FALSE, outlier.shape = 3, coef = 0, linewidth = 1.4)+
  theme_classic()+
  theme(axis.text=element_text(size=20, face="bold"),
        axis.title=element_text(size=20, face="bold"),
        legend.text = element_text(size = 20, face="bold"),
        legend.title = element_blank(),
        legend.position = "none",
        title = element_text(size=20))+
  stat_summary(geom="point", 
               shape=20, size=5, color="black", fill="black",fun = "mean") +
  
  xlab("")+
  ylab("Events/MY")+
  scale_fill_manual(values= c("#DADAEB", "#6A51A3"))+
  geom_signif(comparisons = list(c("Yes Event", "No Event")),
              map_signif_level = TRUE,size = 1,
              textsize = 10)+
  ylim(0,200)


# ============================================================================
#this section creates the figures 3 G to J and calculates their p values 
# ============================================================================

nodes = read.delim("nodes.txt")
p7 = ggplot(nodes, aes(x=factor(WGD_node_condition), y=Normalized_events_MYA, fill = factor(WGD_node_condition)))+
  geom_boxplot(notch = FALSE, outlier.shape = 3, coef = 0, linewidth = 1.4)+
  theme_classic()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        title = element_text(size=20,face="bold"))+
  scale_fill_manual(values= c("#DADAEB", "#6A51A3"))+
  stat_summary(geom="point", 
               shape=20, size=5, color="black", fill="black") +
  
  xlab("")+
  ylab("Events/MYA")+
  geom_signif(comparisons = list(c("WGD node", "Other nodes")),
              map_signif_level = TRUE,size = 1,
              textsize = 5)+
  ylim(0,700)

aggregate(Normalized_events_MYA ~ WGD_node_condition, data = nodes, FUN = mean)
aggregate(Normalized_events_length ~ WGD_node_condition, data = nodes, FUN = mean)

wilcox.test(Normalized_events_MYA ~ WGD_node_condition, data = nodes)
wilcox.test(Normalized_events_length ~ WGD_node_condition, data = nodes)

p8 = ggplot(nodes, aes(x=factor(WGD_node_condition), y=Normalized_events_length, fill = factor(WGD_node_condition)))+
  geom_boxplot(notch = FALSE, outlier.shape = 3, coef = 0, linewidth = 1.4)+
  theme_classic()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        title = element_text(size=20,face="bold"))+
  scale_fill_manual(values= c("#DADAEB", "#6A51A3"))+
  stat_summary(geom="point", 
               shape=20, size=5, color="black", fill="black") +
  
  xlab("")+
  ylab("Events/branch")+
  geom_signif(comparisons = list(c("WGD node", "Other nodes")),
              map_signif_level = TRUE,size = 1,
              textsize = 5)+
  ylim(0,70)


p10 = ggplot(nodes, aes(x=factor(WGD_node_condition_2), y=Normalized_events_length, fill = factor(WGD_node_condition_2)))+
  geom_boxplot(notch = FALSE, outlier.shape = 3, coef = 0, linewidth = 1.4)+
  theme_classic()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        title = element_text(size=20,face="bold"))+
  scale_fill_manual(values= c("#DADAEB", "#6A51A3"))+
  stat_summary(geom="point", 
               shape=20, size=5, color="black", fill="black") +
  
  xlab("")+
  ylab("Events/branch")+
  geom_signif(comparisons = list(c("WGD Lagged node", "Other nodes")),
              map_signif_level = TRUE,size = 1,
              textsize = 5)+
  ylim(0,70)




p9 = ggplot(nodes, aes(x=factor(WGD_node_condition_2), y=Normalized_events_MYA, fill = factor(WGD_node_condition_2)))+
  geom_boxplot(notch = FALSE, outlier.shape = 3, coef = 0, linewidth = 1.4)+
  theme_classic()+
  theme(axis.text=element_text(size=20,face="bold"),
        axis.title=element_text(size=20,face="bold"),
        legend.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.position = "none",
        title = element_text(size=20,face="bold"))+
  scale_fill_manual(values= c("#DADAEB", "#6A51A3"))+
  stat_summary(geom="point", 
               shape=20, size=5, color="black", fill="black") +
  
  xlab("")+
  ylab("Events/MYA")+
  geom_signif(comparisons = list(c("WGD Lagged node", "Other nodes")),
              map_signif_level = TRUE, size = 1,
              textsize = 10)+
  ylim(0,700)

aggregate(Normalized_events_MYA ~ WGD_node_condition_2, data = nodes, FUN = mean)
aggregate(Normalized_events_length ~ WGD_node_condition_2, data = nodes, FUN = mean)
wilcox.test(Normalized_events_MYA ~ WGD_node_condition_2, data = nodes)
wilcox.test(Normalized_events_length ~ WGD_node_condition_2, data = nodes)



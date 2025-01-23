#Pp-AIP_metagenome_code
library(ggplot2)
library(RColorBrewer)
library(scales)
library(ggpubr)
library(ggforce)

##Figure S23 (Paenibacillus_polymyxa_abun)
pp  <- read.csv("./00.data/Paenibacillus_polymyxa_abun.csv",header = T,sep=",",row.names = 1)

pp$group= c( rep('blank_control',5), rep('WT',5), rep('dagr',5))
pp$group <- factor(pp$group,levels = c("Blank","WT","dagr"),ordered=TRUE)

ggplot(pp, aes(x=group, y=Abundance,fill=group)) +
  theme_bw()+
  stat_boxplot(mapping=aes(x=group,y=Abundance),
               width=0.15,position=position_dodge(0.8))+ 
  geom_boxplot(aes(x=group,y=Abundance,fill=group),  
               position=position_dodge(0.8),           
               width=0.6,                                 
               outlier.color = "black")+
  stat_compare_means(method = "wilcox.test",comparisons=list(c("blank_control", "WT"),c("blank_control", "dagr"),c("dagr","WT")),label = "p.signif",hide.ns =FALSE)+
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  theme(panel.grid=element_blank(),
        axis.title.y = element_text(color = 'black',size = 16),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = 'black',size = 14),
        axis.text.x = element_text(color = 'black',size = 14,angle = 45,hjust = 1),
        legend.position = 'none',
        legend.text = element_blank())+
  labs(y="Relative abundance of P. polymyxa (%)")

##Figure S24A & B (phylum)
library(reshape2)
td  <- read.csv("./00.data/Phylum_profile_bac.csv",header = T,sep=",",row.names = 1)
td$phylum<-factor(rownames(td), levels = rownames(td))
data <- melt(td, id = 'phylum')
names(data) <- c('Phylum','Sample','Abundance')
group = c( rep('Blank',55), rep('WT',55), rep('dagr',55))
data$group = group
data_p <- data
data_p$Phylum <- factor(data_p$Phylum,levels = c("Pseudomonadota","Actinomycetota","Bacteroidota","Bacillota","Planctomycetota","Verrucomicrobiota","Myxococcota","Cyanobacteriota","Thermodesulfobacteriota","Acidobacteriota","Others"),ordered=TRUE)
data_p$group <- factor(data_p$group,levels = c("Blank","WT","dagr"),ordered=TRUE)

ggplot(data_p, aes(x=Sample,y=Abundance,fill= Phylum)) +
  theme_bw()+
  geom_bar(
    stat = "identity", 
    width = 0.6)+
  facet_row(~group, scales = "free_x", space = "free")+
  scale_fill_manual(
    values = c('Pseudomonadota'="#BEBADA",'Actinomycetota'="#8DD3C7",'Bacteroidota'="#FB8072",'Bacillota'="#B3DE69",'Planctomycetota'="#FFFFB3",'Verrucomicrobiota'="#BC80BD",'Myxococcota'="#FDB462",'Cyanobacteriota'="#80B1D3",'Thermodesulfobacteriota'="#FCCDE5",'Acidobacteriota'="#CCEBC5",'Others'="#D9D9D9"))+
  theme(panel.grid=element_blank(),
        axis.title.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = 'black',size = 11),
        axis.text.x = element_text(color = 'black',size = 11,angle = 45,hjust = 1), #x轴标签偏???45°，并下降0.5
        #panel.grid = element_blank(),
        legend.position = 'right',
        legend.text = element_text(size = 10))+
  labs(y="Relative abundance (%)")

ggplot(data_p, aes(x=group, y=Abundance)) +
  theme_bw()+
  stat_boxplot(mapping=aes(x=group,y=Abundance),
               width=0.15,position=position_dodge(0.8))+  
  geom_boxplot(aes(x=group,y=Abundance,fill=group), 
               position=position_dodge(0.8),       
               width=0.6,                        
               outlier.color = "black")+
  facet_wrap(~Phylum,nrow=2,scale = "free")+
  stat_compare_means(method = "wilcox.test",comparisons=list(c("Blank", "WT"),c("Blank", "dagr"),c("dagr","WT")),label = "p.signif",hide.ns =FALSE)+
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  theme(panel.grid=element_blank(),
        axis.title.y = element_text(face = 'bold',color = 'black',size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(color = 'black',size = 11),
        axis.text.x = element_text(color = 'black',size = 11,angle = 45,hjust = 1),
        legend.position = 'right',
        legend.text = element_text(size = 10))+
  labs(y="Relative abundance (%)")

##Figure S25 (Species_alpha-diversity)
###devtools::install_github("thomazbastiaanssen/Tjazi")
library(Tjazi)
td  <- read.csv("./00.data/Species_profile_abun_3group.csv",header = T,sep=",",row.names = 1)
###Remove features with prevalence < 10% in two steps:
###First, determine how often every feature is absent in a sample
n_zeroes <- rowSums(td == 0)
###Then, remove features that are absent in more than your threshold (90% in this case).
td <- td[n_zeroes <= round(ncol(td) * 0.90),]
###Perform a CLR transformation
td.exp <- clr_c(td)

library(vegan)
Shannon <- diversity(td, index = "shannon", MARGIN = 2, base = exp(1))
Simpson <- diversity(td, index = "simpson", MARGIN = 2, base = exp(1))
invSimpson <- diversity(td, index = "invsimpson", MARGIN = 2, base = exp(1))
Richness <- specnumber(td, MARGIN = 2)#spe.rich =sobs
index <- as.data.frame(cbind(Shannon, Simpson, invSimpson, Richness))
tdf <- t(td)
tdf<-ceiling(as.data.frame(t(td)))

obs_chao_ace <- t(estimateR(tdf))
obs_chao_ace <- obs_chao_ace[rownames(index),]

index$Chao <- obs_chao_ace[,2]
index$Ace <- obs_chao_ace[,4]
index$obs <- obs_chao_ace[,1]

index$Pielou <- Shannon / log(Richness, 2)
index$Goods_coverage <- 1 - colSums(td ==1) / colSums(td)
write.table(cbind(sample=c(rownames(index)),index),'./00.data/diversity.index_abun.csv', row.names = F, sep = ',', quote = F)

data_ggplot<- read.csv("00.data/diversity.index_abun.csv",header = T,sep=",",row.names = 1)
data_ggplot$group= c( rep('Blank',5), rep('WT',5), rep('dagr',5))
data_ggplot$group <- factor(data_ggplot$group,levels = c("Blank","WT","dagr"),ordered=TRUE)

p1<-ggplot(data_ggplot, aes(x=group, y=Simpson,fill=group))+
  geom_boxplot()+
  theme_bw()+
  stat_compare_means(method = "wilcox.test",comparisons=list(c("Blank", "WT"),c("Blank", "dagr"),c("dagr","WT")),label = "p.signif",hide.ns =FALSE)+
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  labs(title="Alpha diversity (Species level)", x=" ", y="Simpson index")+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5), legend.title=element_blank())
p2<-ggplot(data_ggplot, aes(x=group, y=Shannon,fill=group))+
  geom_boxplot()+
  theme_bw()+
  stat_compare_means(method = "wilcox.test",comparisons=list(c("Blank", "WT"),c("Blank", "dagr"),c("dagr","WT")),label = "p.signif",hide.ns =FALSE)+
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  labs(title="Alpha diversity (Species level)", x=" ", y="Shannon index")+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5), legend.title=element_blank())
p3<-ggplot(data_ggplot, aes(x=group, y=Richness,fill=group))+
  geom_boxplot()+
  theme_bw()+
  stat_compare_means(method = "wilcox.test",comparisons=list(c("Blank", "WT"),c("Blank", "dagr"),c("dagr","WT")),label = "p.signif",hide.ns =FALSE)+
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  labs(title="Alpha diversity (Species level)", x=" ", y="Richness index")+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5), legend.title=element_blank())
p4<-ggplot(data_ggplot, aes(x=group, y=Chao,fill=group))+
  geom_boxplot()+
  theme_bw()+
  stat_compare_means(method = "wilcox.test",comparisons=list(c("Blank", "WT"),c("Blank", "dagr"),c("dagr","WT")),label = "p.signif",hide.ns =FALSE)+
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  labs(title="Alpha diversity (Species level)", x=" ", y="Chao index")+
  theme(panel.grid = element_blank(),
        legend.position = 'none',
        plot.title=element_text(hjust=0.5), legend.title=element_blank())
ggarrange(p2,p1,p4,ncol = 3, nrow = 1,common.legend = TRUE,align = "hv",legend = "bottom")

##Figure 5B (Species_beta-diversity)
library(reshape2)
library(ape)
library(vegan)
sp2 = data.frame(td.exp)
sample_id=colnames(sp2)
group=c( rep('Blank',5), rep('WT',5), rep('dagr',5))
data_group=data.frame(sample_id, group)
bray_dist<-vegdist(t(sp2),method = "euclidean")
data_pcoa<-pcoa(bray_dist,correction = "cailliez")
data_plot<-data.frame(data_pcoa$vectors)
###PERMANOVA
dune.div <- adonis2(bray_dist~group, data_group, permutations = 9999, method="euclidean")
dune_adonis <- paste0("adonis R2: ",round(dune.div$R2,2), "; P-value: ", dune.div$`Pr(>F)`)

x_label<-round(data_pcoa$values$Relative_eig[1]*100,2)
y_label<-round(data_pcoa$values$Relative_eig[2]*100,2)

data_plot$group<-data_group$group
pB1<-ggplot(data=data_plot,aes(x=Axis.1,y=Axis.2,color=group,fill=group))+
  geom_point(size=5,)+
  scale_color_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.title.y = element_text(color = 'black',size = 16),
        axis.title.x = element_text(color = 'black',size = 16),
        axis.text.y = element_text(color = 'black',size = 14),
        axis.text.x = element_text(color = 'black',size = 14),
        legend.position = 'none',
        legend.text = element_blank())+
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  labs(x=paste0("PCoA1 ",x_label,"%"),y=paste0("PCoA2 ",y_label,"%"))+
  stat_ellipse(level = 0.95,data=data_plot,geom = "polygon",aes(fill = group),alpha=0.3)+
  xlim(-60,60)+ylim(-65,65)

data_plot$group<-factor(data_plot$group,levels = c("dagr","WT","Blank"),ordered=TRUE)
pB2<-ggplot(data_plot,aes(x=group,y=Axis.1)) +
  stat_boxplot(mapping=aes(x=group,y=Axis.1),
               width=0.15,position=position_dodge(0.8))+ 
  geom_boxplot(aes(x=group,y=Axis.1,fill=group),   
               position=position_dodge(0.8),    
               width=0.6,                          
               outlier.color = "black")+
  stat_compare_means(method = "wilcox.test",comparisons=list(c("Blank", "WT"),c("Blank", "dagr"),c("WT", "dagr")),label = "p.signif",hide.ns =FALSE)+
  coord_flip() +
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+
  ylim(-60,60)

pB3<-ggplot(data_plot,aes(x=rev(group),y=Axis.2)) +
  stat_boxplot(mapping=aes(x=rev(group),y=Axis.2),
               width=0.15,position=position_dodge(0.8))+   
  geom_boxplot(aes(x=rev(group),y=Axis.2,fill=group),   
               position=position_dodge(0.8),          
               width=0.6,                      
               outlier.color = "black")+
  stat_compare_means(method = "wilcox.test",comparisons=list(c("Blank", "WT"),c("Blank", "dagr"),c("WT", "dagr")),label = "p.signif",hide.ns =FALSE)+
  scale_fill_manual(
    values = c('blank_control'="#999999",'WT'="#F9BE80",'dagr'="#6A98CC"))+
  theme_bw()+
  theme(panel.grid=element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y=element_blank(),
        legend.position = "none")+
  ylim(-65,65)

pB4<-ggplot(data_plot, aes(Axis.1, Axis.2)) +
  geom_text(aes(x = -0.5,y = 0.6,label = paste("PERMANOVA:\ndf = ",dune.div$Df[1],  "\nR2 = ",round(dune.div$R2[1],2),  "\np-value = ",dune.div$`Pr(>F)`[1],sep = "")),
            size = 5) +
  theme_bw() +
  xlab("") + ylab("") +
  theme(panel.grid=element_blank(), 
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())
ggarrange(pB2,pB4,pB1,pB3,ncol = 2, nrow = 2, legend = "none",align = "hv",widths = c(3, 1), heights = c(1, 3))

##Figure S26 (Species_sig_ternary)
library(vcd)
data = read.table("00.data/cluster_source_N.csv",head=T,sep=",",row.names = 1,comment.char="$")
color<-color2<-data[,5]
ternaryplot(data[,1:3],cex =0.8,scale = 100,
            col = color,
            #prop_size=6,
            grid_color = "grey50",labels_color = "black",main =" ",labels = "outside")
grid_legend("topright", pch=20, col=c("#BEBADA","#8DD3C7","#FB8072","#B3DE69","#FFFFB3",
                                      "#BC80BD","#FDB462","#80B1D3","#FCCDE5","#CCEBC5",
                                      "#D9D9D9","#FFFFFF"), 
labels=c("Pseudomonadota","Actinomycetota","Bacteroidota","Bacillota",'Planctomycetota',
         "Verrucomicrobiota","Myxococcota","Cyanobacteriota","Thermodesulfobacteriota",
            "Acidobacteriota","Others","nosig"), title = "group")
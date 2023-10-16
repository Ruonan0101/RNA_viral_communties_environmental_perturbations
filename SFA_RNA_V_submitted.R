library(ggplot2)
library(ggsignif)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(scales)
library(ggsankey)
library(dplyr)
library(gridExtra)
library(vegan)
library(ggdendro)
library(ggtext)
library(rpart)
library(partykit)
library(data.tree)
library(ggtree)
library(corrplot)
library(Hmisc)
library(igraph)
library(ggrepel)
library(RColorBrewer)
library(colorRamps)
library(rfPermute)
library(randomForest)
library(A3)
library(gridExtra)
library(reshape2)
library(lavaan)
library(semPlot)
library(pheatmap)


##Fig2##sankey####RNA V overview####
setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/')
d<-read.csv("RNA_overview_forSankey.csv", header = TRUE)
df <- d %>%
  make_long(Novelty,	Viral_taxonomy_phylum, Host_range,	Host_type)
TotalCount = nrow(d)
dagg <- df%>%
  dplyr::group_by(node)%>%
  tally()

dagg <- dagg%>%
  dplyr::group_by(node)%>%
  dplyr::mutate(pct = n/TotalCount)
df2 <- merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)

mycolors<-c(brewer.pal(name="Set1", n = 5),brewer.pal(name="Set2", n = 5), alpha=0.7)

fig2a <- ggplot(df2, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = paste0(node,"\n n=", n, ' (',  round(pct* 100,1), '%)' )))+
  geom_sankey(flow.alpha = 0.5,  color = "gray40", show.legend = TRUE)+
  geom_sankey_label(size = 4, color = "white", fill= "gray40")+  
  theme_bw()+theme_sankey(base_size = 10) + 
  theme(legend.position = "none", axis.text.x = element_blank())+
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank())+ 
  scale_color_manual(values=mycolors)+ labs(fill = 'Nodes')
#+scale_fill_viridis_d(option = "plasma")+ labs(fill = 'Nodes')
fig2a

####pie chart ###
setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/')
data<-read.csv("RNAV_composition_NEW.csv",row.names=1,check.names=FALSE)
data <- data %>% 
  arrange(desc(Group)) %>%
  mutate(prop = Value / sum(data$Value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )%>%
  mutate(ymax = cumsum(prop))%>%
  mutate(ymin = c(0, head(ymax, n=-1)))%>%
  mutate(labelPosition = (ymax + ymin) / 2)%>%
  mutate(label= paste0(Group,'\n', percent(Value/sum(Value))))


fig2b<-ggplot(data, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Group)) +
  geom_rect() +
  geom_label( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  scale_fill_manual(values=c('#9CA700', '#00A7FF', '#7F96FF'))
fig2b

data2<-read.csv("RNAV_family_composition_PieChart_minorGroup.csv",row.names=1,check.names=FALSE)
data2 <- data2 %>% 
  arrange(desc(Group)) %>%
  mutate(prop = Value / sum(data2$Value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )%>%
  mutate(ymax = cumsum(prop))%>%
  mutate(ymin = c(0, head(ymax, n=-1)))%>%
  mutate(labelPosition = (ymax + ymin) / 2)%>%
  mutate(label= paste0(Group,'\n', percent(Value/sum(Value))))


fig2c<-ggplot(data2, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=Group)) +
  geom_rect() +
  geom_label_repel( x=3.5, aes(y=labelPosition, label=label), size=3) +
  scale_fill_brewer(palette=4) +
  coord_polar(theta="y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")+
  scale_fill_manual(values=c('#e5c494','#fc8d62', '#8da0cb','#ffd92f','#e78ac3','#7F96FF'))
fig2c

grid.arrange(fig2a, arrangeGrob(fig2b, fig2c),ncol=2, widths = c(6, 2.7))


#####fig3 PC####
setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/')
pca_5com<-read.csv("cummu_pca_fiveCommunity.csv")
pca_3com<-read.csv("cummu_pca_3Community.csv")

fig3<-ggplot(pca_5com, aes(x=Principle_component, y=Value, color=Group)) + geom_point()+geom_line()+
  labs(x='Principle component',y='Cumulative proportion of variance explained')+
  scale_color_manual(values=c("#008000", '#9CA700', "#0000FF",'#00A7FF','#937cd9'))+
  theme_bw()
fig3

####fig4 
setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/OTU_table/')
abund_table_Total<-read.csv("Total_RNA_V_con_forPCA.csv",row.names=1,check.names=FALSE)
abund_table_Prok<-read.csv("Prok_RNA_V_com_forPCA.csv",row.names=1,check.names=FALSE)
abund_table_Euk<-read.csv("Euk_RNA_V_com_forPCA.csv",row.names=1,check.names=FALSE)


abund_table_Total<-t(abund_table_Total)
abund_table_Prok<-t(abund_table_Prok)
abund_table_Euk<-t(abund_table_Euk)

grouping_info<-data.frame(row.names=rownames(abund_table_Total),t(as.data.frame(strsplit(rownames(abund_table_Total),"_"))))

betad_Total<-vegdist(abund_table_Total,method="bray")
betad_Prok<-vegdist(abund_table_Prok,method="bray")
betad_Euk<-vegdist(abund_table_Euk,method="bray")

sol_Total<-metaMDS(abund_table_Total,distance = "bray", k = 6, trymax = 50)
sol_Euk<-metaMDS(abund_table_Euk,distance = "bray", k = 6, trymax = 50)
sol_Prok<-metaMDS(abund_table_Prok,distance = "bray", k = 6, trymax = 50)

MDS_Total <- cmdscale(vegdist(abund_table_Total, method = "bray"), k = 2, eig = T, add = T )
MDS_Euk <- cmdscale(vegdist(abund_table_Euk, method = "bray"), k = 2, eig = T, add = T )
MDS_Prok <- cmdscale(vegdist(abund_table_Prok, method = "bray"), k = 2, eig = T, add = T )

res_adonis_Total <- adonis2(betad_Total ~ X2, grouping_info)
res_adonis_Euk <- adonis2(betad_Euk ~ X2, grouping_info)
res_adonis_Prok <- adonis2(betad_Prok ~ X2, grouping_info)

res_adonis_Total
res_adonis_Euk
res_adonis_Prok

NMDS_Total=data.frame(x=sol_Total$point[,1],y=sol_Total$point[,2],Root=as.factor(grouping_info[,2]))
NMDS_Euk=data.frame(x=sol_Euk$point[,1],y=sol_Euk$point[,2],Root=as.factor(grouping_info[,2]))
NMDS_Prok=data.frame(x=sol_Prok$point[,1],y=sol_Prok$point[,2],Root=as.factor(grouping_info[,2]))

totalRootE<-ggplot(data=NMDS_Total,aes(x,y))+geom_point(aes(color=Root), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of plant on total RNA viral community, Pr(>F):0.001***")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))
EukRRootE<-ggplot(data=NMDS_Euk,aes(x,y))+geom_point(aes(color=Root), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of plant on eukaryotic RNA viral community, Pr(>F):0.001***")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))
ProkRRootE<-ggplot(data=NMDS_Prok,aes(x,y))+geom_point(aes(color=Root), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of plant on prokaryotic RNA viral community, Pr(>F):0.001***")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))

setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/')
##root
df<-read.csv("JoseAlkar_Bare_T4_0-5_RootE_VC.csv", header = TRUE)
cultivar_root<-ggboxplot(df, x = 'Root', y = 'RNA_Viruses', color = 'Cultivar')+
  facet_grid(Variable~Treatment,scales='free')+
  scale_color_manual(values=c('#fc8d59','#999999','#91bfdb'))+
  #stat_pvalue_manual(stat.test, label = "p.adj.signif",hide.ns = TRUE, tip.length = 0)+
  theme(legend.position = 'right', strip.placement='outside', axis.title = element_blank(),axis.text = element_text(size = 10),
        plot.title = element_text(size=12),strip.text = element_blank())

grid.arrange(totalRootE,EukRRootE,ProkRRootE,cultivar_root, ncol=2)


####fig5####
setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/OTU_table/')
###Depth
abund_table_Total<-read.csv("Total_RNA_V_con_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)
abund_table_Prok<-read.csv("Prok_RNA_V_com_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)
abund_table_Euk<-read.csv("Euk_RNA_V_com_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)


abund_table_Total<-t(abund_table_Total)
abund_table_Prok<-t(abund_table_Prok)
abund_table_Euk<-t(abund_table_Euk)

grouping_info<-data.frame(row.names=rownames(abund_table_Total),t(as.data.frame(strsplit(rownames(abund_table_Total),"_"))))

betad_Total<-vegdist(abund_table_Total,method="bray")
betad_Prok<-vegdist(abund_table_Prok,method="bray")
betad_Euk<-vegdist(abund_table_Euk,method="bray")

sol_Total<-metaMDS(abund_table_Total,distance = "bray", k = 6, trymax = 50)
sol_Euk<-metaMDS(abund_table_Euk,distance = "bray", k = 6, trymax = 50)
sol_Prok<-metaMDS(abund_table_Prok,distance = "bray", k = 6, trymax = 50)

MDS_Total <- cmdscale(vegdist(abund_table_Total, method = "bray"), k = 2, eig = T, add = T )
MDS_Euk <- cmdscale(vegdist(abund_table_Euk, method = "bray"), k = 2, eig = T, add = T )
MDS_Prok <- cmdscale(vegdist(abund_table_Prok, method = "bray"), k = 2, eig = T, add = T )

res_adonis_Total <- adonis2(betad_Total ~ X6, grouping_info)
res_adonis_Euk <- adonis2(betad_Euk ~ X6, grouping_info)
res_adonis_Prok <- adonis2(betad_Prok ~ X6, grouping_info)

res_adonis_Total
res_adonis_Euk
res_adonis_Prok

NMDS_Total=data.frame(x=sol_Total$point[,1],y=sol_Total$point[,2],Depth=as.factor(grouping_info[,6]))
NMDS_Euk=data.frame(x=sol_Euk$point[,1],y=sol_Euk$point[,2],Depth=as.factor(grouping_info[,6]))
NMDS_Prok=data.frame(x=sol_Prok$point[,1],y=sol_Prok$point[,2],Depth=as.factor(grouping_info[,6]))

totalDepthE<-ggplot(data=NMDS_Total,aes(x,y))+geom_point(aes(color=Depth), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of depth on total RNA viral community, Pr(>F):0.005**")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#b37c27','#573603'))
EukRDepthE<-ggplot(data=NMDS_Euk,aes(x,y))+geom_point(aes(color=Depth), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of depth on eukaryotic RNA viral community, Pr(>F):0.047*")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#b37c27','#573603'))
ProkRDepthE<-ggplot(data=NMDS_Prok,aes(x,y))+geom_point(aes(color=Depth), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of depth on prokaryotic RNA viral community, Pr(>F):0.001***")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#b37c27','#573603'))

###Irrigation
abund_table_Total<-read.csv("Total_RNA_V_con_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)
abund_table_Prok<-read.csv("Prok_RNA_V_com_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)
abund_table_Euk<-read.csv("Euk_RNA_V_com_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)


abund_table_Total<-t(abund_table_Total)
abund_table_Prok<-t(abund_table_Prok)
abund_table_Euk<-t(abund_table_Euk)

grouping_info<-data.frame(row.names=rownames(abund_table_Total),t(as.data.frame(strsplit(rownames(abund_table_Total),"_"))))

betad_Total<-vegdist(abund_table_Total,method="bray")
betad_Prok<-vegdist(abund_table_Prok,method="bray")
betad_Euk<-vegdist(abund_table_Euk,method="bray")

sol_Total<-metaMDS(abund_table_Total,distance = "bray", k = 6, trymax = 50)
sol_Euk<-metaMDS(abund_table_Euk,distance = "bray", k = 6, trymax = 50)
sol_Prok<-metaMDS(abund_table_Prok,distance = "bray", k = 6, trymax = 50)

MDS_Total <- cmdscale(vegdist(abund_table_Total, method = "bray"), k = 2, eig = T, add = T )
MDS_Euk <- cmdscale(vegdist(abund_table_Euk, method = "bray"), k = 2, eig = T, add = T )
MDS_Prok <- cmdscale(vegdist(abund_table_Prok, method = "bray"), k = 2, eig = T, add = T )

res_adonis_Total <- adonis2(betad_Total ~ X4, grouping_info)
res_adonis_Euk <- adonis2(betad_Euk ~ X4, grouping_info)
res_adonis_Prok <- adonis2(betad_Prok ~ X4, grouping_info)

res_adonis_Total
res_adonis_Euk
res_adonis_Prok

NMDS_Total=data.frame(x=sol_Total$point[,1],y=sol_Total$point[,2],Irrigation=as.factor(grouping_info[,4]),Depth=as.factor(grouping_info[,6]),Cultivar=as.factor(grouping_info[,3]))
NMDS_Euk=data.frame(x=sol_Euk$point[,1],y=sol_Euk$point[,2],Irrigation=as.factor(grouping_info[,4]))
NMDS_Prok=data.frame(x=sol_Prok$point[,1],y=sol_Prok$point[,2],Irrigation=as.factor(grouping_info[,4]))

totalWaterE<-ggplot(data=NMDS_Total,aes(x,y))+geom_point(aes(color=Irrigation), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of irrigation on total RNA viral community, Pr(>F):0.002**")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#105deb','#7ca1e6'))
EukRWaterE<-ggplot(data=NMDS_Euk,aes(x,y))+geom_point(aes(color=Irrigation), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of irrigation on eukaryotic RNA viral community, Pr(>F):0.008**")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#105deb','#7ca1e6'))
ProkRWaterE<-ggplot(data=NMDS_Prok,aes(x,y))+geom_point(aes(color=Irrigation), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of irrigation on prokaryotic RNA viral community, Pr(>F):0.036*")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#105deb','#7ca1e6'))


####Cultivar
abund_table_Total<-read.csv("Total_RNA_V_con_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)
abund_table_Prok<-read.csv("Prok_RNA_V_com_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)
abund_table_Euk<-read.csv("Euk_RNA_V_com_forPCA_BareRemove.csv",row.names=1,check.names=FALSE)


abund_table_Total<-t(abund_table_Total)
abund_table_Prok<-t(abund_table_Prok)
abund_table_Euk<-t(abund_table_Euk)

grouping_info<-data.frame(row.names=rownames(abund_table_Total),t(as.data.frame(strsplit(rownames(abund_table_Total),"_"))))

betad_Total<-vegdist(abund_table_Total,method="bray")
betad_Prok<-vegdist(abund_table_Prok,method="bray")
betad_Euk<-vegdist(abund_table_Euk,method="bray")

sol_Total<-metaMDS(abund_table_Total,distance = "bray", k = 6, trymax = 50)
sol_Euk<-metaMDS(abund_table_Euk,distance = "bray", k = 6, trymax = 50)
sol_Prok<-metaMDS(abund_table_Prok,distance = "bray", k = 6, trymax = 50)

MDS_Total <- cmdscale(vegdist(abund_table_Total, method = "bray"), k = 2, eig = T, add = T )
MDS_Euk <- cmdscale(vegdist(abund_table_Euk, method = "bray"), k = 2, eig = T, add = T )
MDS_Prok <- cmdscale(vegdist(abund_table_Prok, method = "bray"), k = 2, eig = T, add = T )

res_adonis_Total <- adonis2(betad_Total ~ X3, grouping_info)
res_adonis_Euk <- adonis2(betad_Euk ~ X3, grouping_info)
res_adonis_Prok <- adonis2(betad_Prok ~ X3, grouping_info)

res_adonis_Total
res_adonis_Euk
res_adonis_Prok

NMDS_Total=data.frame(x=sol_Total$point[,1],y=sol_Total$point[,2],Cultivar=as.factor(grouping_info[,3]))
NMDS_Euk=data.frame(x=sol_Euk$point[,1],y=sol_Euk$point[,2],Cultivar=as.factor(grouping_info[,3]))
NMDS_Prok=data.frame(x=sol_Prok$point[,1],y=sol_Prok$point[,2],Cultivar=as.factor(grouping_info[,3]))

totalCultivarE<-ggplot(data=NMDS_Total,aes(x,y))+geom_point(aes(color=Cultivar), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of cultivar on total RNA viral community, Pr(>F):0.033*")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#fc8d59','#91bfdb'))
EukRCultivarE<-ggplot(data=NMDS_Euk,aes(x,y))+geom_point(aes(color=Cultivar), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of cultivar on eukaryotic RNA viral community, Pr(>F):0.044*")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#fc8d59','#91bfdb'))
ProkRCultivarE<-ggplot(data=NMDS_Prok,aes(x,y))+geom_point(aes(color=Cultivar), size=1)+theme_bw()+
  labs(x='NMDS1', y='NMDS2', subtitle="Impact of cultivar on prokaryotic RNA viral community, Pr(>F):0.172")+
  theme(axis.title=element_text(size=8), plot.subtitle=element_text(size=10))+
  scale_color_manual(values=c('#fc8d59','#91bfdb'))

grid.arrange(arrangeGrob(totalDepthE, EukRDepthE, ProkRDepthE), arrangeGrob(totalWaterE, EukRWaterE, ProkRWaterE), arrangeGrob(totalCultivarE, EukRCultivarE, ProkRCultivarE), ncol=3)

####Fig6###
setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/')
parametertable<-read.csv("ENV_meta_4V_Prok_Euk_com_cat_forCorrelation_update.csv",row.names=1,check.names=FALSE)
parametertable<-read.csv("ENV_meta_4V_Prok_Euk_com_cat_forCorrelation_update_reduced.csv",row.names=1,check.names=FALSE)

res <- rcorr(as.matrix(parametertable))
head(res$r)
head(res$P)
###fig6a
corrplot(res$r, type="upper", order = 'alphabet', p.mat = res$P, sig.level = 0.05, insig = "blank",tl.cex=0.7, diag=F,tl.col = "black") 


library(rfPermute)
library(randomForest)
setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/')
parameterall = read.csv('ENV_meta_4V_Prok_Euk_com_cat_forCorrelation_update_reduced.csv', row.names=1)
TotalRNAViralAbundance.rp <- rfPermute(TotalRNAViralAbundance ~ ., data = parameterall, na.action = na.omit, ntree = 100, num.rep = 50)

#plotTrace(TotalRNAViralAbundance.rp)
library(A3)

TotalRNAViralAbundance.rp <- rfPermute(TotalRNAViralAbundance ~ Planted_Bare+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, data = parameterall, na.action = na.omit, ntree = 100, num.rep = 50)
a3(TotalRNAViralAbundance ~ Root+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, parameterall, randomForest, n.folds = 5, p.acc = 0.01)
plotImportance(TotalRNAViralAbundance.rp, scale = TRUE, size=1.5, alpha=0.05, plot.type = 'heatmap', main = 'Total RNA viral abundance\nR^2=41.4%, P<0.01')


TotalRNAViralRichness.rp <- rfPermute(TotalRNAViralRichness ~ Planted_Bare+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, data = parameterall, na.action = na.omit, ntree = 100, num.rep = 50)
a3(TotalRNAViralRichness ~ Root+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, parameterall, randomForest, n.folds = 5, p.acc =0.01)
plotImportance(TotalRNAViralRichness.rp, scale = TRUE, size=1.5, alpha=0.05, plot.type = 'heatmap', main = 'Total RNA viral richness\nR^2=63.8%, P<0.01')

Prok.RNAViralAbundance.rp <- rfPermute(Prok.RNAViralAbundance ~ Planted_Bare+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, data = parameterall, na.action = na.omit, ntree = 100, num.rep = 50)
a3(Prok.RNAViralAbundance ~ Root+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, parameterall, randomForest, n.folds = 5, p.acc = 0.01)
plotImportance(Prok.RNAViralAbundance.rp, scale = TRUE, size=1.5, alpha=0.05, plot.type = 'heatmap', main = 'Prokaryotic RNA viral abundance\nR^2=30.2%, P<0.01')

Prok.RNAViralRichness.rp <- rfPermute(Prok.RNAViralRichness ~ Planted_Bare+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, data = parameterall, na.action = na.omit, ntree = 100, num.rep = 50)
a3(Prok.RNAViralRichness ~ Root+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, parameterall, randomForest, n.folds = 5, p.acc = 0.01)
plotImportance(Prok.RNAViralRichness.rp, scale = TRUE, size=1.5, alpha=0.05, plot.type = 'heatmap', main = 'Prokaryotic RNA viral richness\nR^2=37.1%, P<0.01')

Euk.RNAViralRichness.rp <- rfPermute(Euk.RNAViralRichness ~ Planted_Bare+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, data = parameterall, na.action = na.omit, ntree = 100, num.rep = 50)
a3(Euk.RNAViralRichness ~ Root+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, parameterall, randomForest, n.folds = 5, p.acc = 0.01)
plotImportance(Euk.RNAViralRichness.rp, scale = TRUE, size=1.5, alpha=0.05, plot.type = 'heatmap', main = 'Eukaryotic RNA viral richness\nR^2=66.4%, P<0.01')

Euk.RNAViralAbundance.rp <- rfPermute(Euk.RNAViralAbundance ~ Planted_Bare+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, data = parameterall, na.action = na.omit, ntree = 100, num.rep = 50)
a3(Euk.RNAViralAbundance ~ Root+Cultivar+Depth+Euk.Abundance+Euk.Richness+Prok.Abundance+Prok.Richness+Water+Euk.PC1+Prok.PC1, parameterall, randomForest, n.folds = 5, p.acc = 0.01)
plotImportance(Euk.RNAViralAbundance.rp, scale = TRUE, size=1.5, alpha=0.05, plot.type = 'heatmap', main = 'Eukaryotic RNA viral abundance\nR^2=23.1%, P<0.03')


####fig7
pairtable<-read.csv('sigPair_correlation_RF_reduced_renamed_cleaned_NEW.csv')
library(igraph)
g1 <- graph_from_data_frame(d=pairtable,directed=F)

plot.igraph(g1, edge.width=E(g1)$value*30)

edge_density(g1, loops=F)
deg <- degree(g1, mode="all")
plot(g1, vertex.size=deg*3)
hist(deg, breaks=1:vcount(g1)-1, main="Histogram of node degree")
deg.dist <- degree_distribution(g1, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")

E(g1)$weight <- abs(E(g1)$value)
ceb <- cluster_edge_betweenness(g1, weights = E(g1)$weight) 
dendPlot(ceb, mode="hclust")
plot(ceb, g1) 
class(ceb)
clp <- cluster_label_prop(g1)
plot(clp, g1)
cfg <- cluster_fast_greedy(as.undirected(g1))
###selected###
plot(cfg, as.undirected(g1),vertex.size=9, vertex.label.cex=0.8, layout=layout.fruchterman.reingold(g1), 
     label.degree = pi/2, edge.color = "grey", label.font=2)



setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/')
d<-read.csv("ENV_meta_4V_Prok_Euk_com_cat_forCorrelation_updateSEM_normalized.csv",row.names=1,check.names=FALSE)

mydata<-na.omit(d)
Model3<-'
Euk. abundance+Euk. richness~Depth+Water
Euk.PC1~Cultivar_Bare
Euk. RNA. V.rich.~ Euk. richness+Euk.PC1+Water+Cultivar_Bare
Euk. RNA.V.abund.~Euk. abundance+Water
'
####  chisq:18.178;     df pvalue    gfi    cfi 
#18.178 13.000  0.151  0.856  0.942 
#rmr   srmr  rmsea    nfi 
#0.425  0.116  0.129  0.841 

#chisq=18.17;pvalue=0.15;cfi= 0.94

Model3Fit <- sem(Model3, data=mydata, fixed.x=F, check.gradient = FALSE)
fitMeasures(Model3Fit, c("chisq","df","pvalue","gfi","cfi","rmr","srmr","rmsea",'nfi'))
varTable(Model3Fit)
table<-parameterEstimates(Model3Fit,standardized=TRUE)
table<-table[!table$lhs==table$rhs,]
b<-gettextf('%.2f\np=%.2f', table$std.all, digits=table$pvalue)

semPaths(Model3Fit, what='std',layout='tree3', ask = FALSE,edgeLabels = b,
         edge.label.cex=0.4,sizeMan=10,label.font=15,
         posCol = c("blue", "red"),
         fade=F, nCharNodes=0, intercepts = F,residuals=F, thresholds = 0.05, legend=FALSE)

grps <- list(Abiotic = c("Depth", "Water", 'Cultivar_Bare'), 
             microbe=c('Euk.abundance','Euk.richness','Euk.PC1'),
             virus = c("Euk.RNA.V.rich.",'Euk.RNA.V.abund.')
)
semPaths(Model3Fit, what='std',layout = 'tree3', ask = FALSE,
         edge.label.cex=0.8,sizeMan=9,
         edgeLabels = b,label.font=50,edge.label.position=0.5,
         groups = grps, color = c("#dfc27d", 
                                  '#04631c',
                                  '#64cc9f'),
         posCol = c("blue", "red"),nodeLabels = 1:8,
         fade=F, nCharNodes=0, intercepts = F,residuals=F, thresholds = F, legend=FALSE)

x = c(-1, 0, 1, -0.5, 1, -0.5,0,0.5)
y = c(0, 0, 0, -1, -1, 1,1,1)
ly = matrix(c(x, y), ncol=2)


semPaths(Model3Fit, what='std',layout = ly, ask = FALSE,
         edge.label.cex=0.7,sizeMan=10,
         edgeLabels = b,label.font=60,edge.label.position=0.5,
         groups = grps, color = c("#dfc27d", 
                                  '#04631c',
                                  '#64cc9f'),
         posCol = c("blue", "red"),
         fade=F, nCharNodes=0, intercepts =F,residuals=F, thresholds = F, legend=FALSE)



#####supplementary figures ##########
library(pheatmap)
setwd('/Users/wuru978/Desktop/2022_Manuscripts/22_P2_RNAV_irrigation_maybeDepth/FigureTable_dir/')
data<-read.delim("euk_all_abund_table.txt",header=T, row.names='cluster')
target<-read.delim("target_NEW.txt",header=T,row.names=1)
#target <- target[, -1, drop=F]
target
cal_z_score <- function(x){
  (log(x+1))
}
data_norm <- t(apply(data, 1, cal_z_score))

row.names(target) <- colnames(data_norm)

ann_colors = list(
  Depth = c(X0_5cm ='#b37c27',X15_25cm ='#573603'),
  Irrigation = c(T1='#105deb',T4='#7ca1e6'),
  Cultivar_root= c(Alkar='#fc8d59',Bare='#999999',Jose='#91bfdb'))

colors<- c('#b37c27','#573603', '#105deb','#7ca1e6', '#fc8d59','#999999','#91bfdb')
  
  
######supplementary figure 1#####
###eukaryotic community####
pheatmap(data_norm,
         annotation=target,
         annotation_colors = ann_colors,
         cutree_cols = 6,
         fontsize = 8,cluster_rows=FALSE, cluster_cols=TRUE,show_colnames =F, show_rownames = T)



######supplementary figure 2#####
######host_taxa#####
data<-read.delim("The_7_host_tax_abundance.txt",header=T, row.names='cluster')
data
target<-read.delim("target_NEW.txt",header=T,row.names=1)
#target <- target[, -1, drop=F]
target
cal_z_score <- function(x){
  (log(x+1))
}
data_norm <- t(apply(data, 1, cal_z_score))
row.names(target) <- colnames(data_norm)

pheatmap(data_norm,
         annotation=target,
         annotation_colors = ann_colors,
         cutree_cols = 6,
         fontsize = 8,cluster_rows=FALSE, cluster_cols=TRUE,show_colnames =F, show_rownames = T)

######supplementary figure 3#####
data<-read.delim("euk_v_comp.txt",header=T, row.names='cluster')
data

cal_z_score <- function(x){
  (log(x+1))
}
data_norm <- t(apply(data, 1, cal_z_score))
row.names(target) <- colnames(data_norm)

pheatmap(data_norm,
         annotation=target,
         annotation_colors = ann_colors,
         cutree_cols = 6,
         fontsize = 8,cluster_rows=FALSE, cluster_cols=TRUE,show_colnames =F, show_rownames = T)

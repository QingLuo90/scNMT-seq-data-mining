
setwd("/Users/qinluo/Desktop/seconde_analysis")
table = read.table("Oct4_TSS_5kb_2.txt",header=T)

library(ggplot2)

#how many sites are tested
length(table(table$Pos))
#[1] 1046


cell_loci_num = as.matrix(table(table$cell))
new_order = order(cell_loci_num[,1],decreasing = T)
cell_loci_num = cell_loci_num[new_order,1]
summary(cell_loci_num)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#1.00   34.00   78.00   97.75  144.00  872.00



length(which(cell_loci_num >= 3))
#914
filtered_cells = names(cell_loci_num)[1:914]

cell = match(table$cell,filtered_cells)
cell_filtered = table[!is.na(cell),]
rownames(cell_filtered) = c(1:dim(cell_filtered)[1])




#####
library(reshape2)

###############cell_filtered_2 is used for binning
cell_filtered_2 = cell_filtered[,c(2,5,6)]

cell_filtered_2_1 = cell_filtered_2[which(cell_filtered_2[,2]==1),]
cell_filtered_2_0 = cell_filtered_2[which(cell_filtered_2[,2]==0),]

cell_plot_1 = dcast(cell_filtered_2_1,Pos~cell)
cell_plot_0 = dcast(cell_filtered_2_0,Pos~cell)
rownames( cell_plot_0) = cell_plot_0[,1]
rownames( cell_plot_1) = cell_plot_1[,1]
cell_plot_0 = cell_plot_0[,-1]
cell_plot_1 = cell_plot_1[,-1]

cell_plot_0[!is.na(cell_plot_0)] <- 0
cell_plot_0[is.na(cell_plot_0)] <- 3
cell_plot_1[!is.na(cell_plot_1)] <- 1
cell_plot_1[is.na(cell_plot_1)] <- 3


r = match(rownames(cell_plot_0),rownames(cell_plot_1))
c = match(colnames(cell_plot_0),colnames(cell_plot_1))
cell_plot_1.2 = cell_plot_1[r,]
cell_plot_1.2 = t(t(cell_plot_1.2)[c,])
cell_plot_1.2[is.na(cell_plot_1.2)] <- 3

dimnames(cell_plot_1.2) = dimnames(cell_plot_0)

write.table(cell_plot_1.2,"cell_plot_1.txt",sep="\t",quote=F)
write.table(cell_plot_0,"cell_plot_0.txt",sep="\t",quote=F)

rt1 = read.table("cell_plot_0.txt",header=T)
rt2 = read.table("cell_plot_1.txt",header=T)

cell_plot = rt1 + rt2 - 3

cell_plot = t(cell_plot)

library(pheatmap)

write.table(cell_plot,"cell_plot_N.txt",sep="\t",row.names = T,quote = F)



rt = read.table("cell_plot_N.txt")
meta = read.table("ACC_QC_meta.txt",header=T)
x = match( rownames(rt),meta$sample)
x = sort(x)
meta = meta[x,]
x = match( meta$sample,rownames(rt))
rt = rt[x,]


cell_TSS_5kb_N_meta_order = rt
meta_order = meta
oct4 = read.table("oct4.txt",header=T)
library(limma)
oct4 = avereps(oct4, ID = oct4$Cell)
x = match( meta$sample, oct4[,1])
oct4 = oct4[x,]
x = as.numeric(oct4[,2])
meta_new = cbind(meta,x)
meta_new_2 = na.omit(meta_new)
colnames(meta_new_2)[13] = "Oct4"

x = match( rownames(cell_TSS_5kb_N_meta_order),meta_new_2$sample)
x = sort(x)
meta_new_2 = meta_new_2[x,]
x = match( meta_new_2$sample,rownames(cell_TSS_5kb_N_meta_order))
cell_TSS_5kb_N_meta_order = cell_TSS_5kb_N_meta_order[x,]

write.table(cell_TSS_5kb_N_meta_order ,"cell_TSS_5kb_N_meta_order.txt",sep="\t",quote = F)
write.table(meta_new_2,"TSS_5kb_meta_order.txt",sep="\t",quote = F,row.names = T)


####plot oct4 expression
ggplot(meta_new_2)+geom_boxplot(aes(x=lineage10x_2, y=Oct4,fill=lineage10x_2))+
    theme_bw()+scale_fill_brewer(palette = "Set1")+coord_flip()+facet_grid(~stage)


##########plot heatmap
index = paste(meta_new_2$stage,meta_new_2$lineage10x_2,sep="_")
#table(index)
#index
#E4.5_Epiblast E4.5_Primitive_endoderm           E5.5_Epiblast  E5.5_Visceral_endoderm
#52                      32                      69                      18
#E6.5_Epiblast       E6.5_ExE_ectoderm           E6.5_Mesoderm   E6.5_Primitive_Streak
#90                       7                       5                      27
#E6.5_Visceral_endoderm           E7.5_Ectoderm           E7.5_Endoderm           E7.5_Epiblast
#11                      37                      38                      30
#E7.5_Mesoderm   E7.5_Primitive_Streak
#109                      26


E4.5_Epi = cell_TSS_5kb_N_meta_order [which(index=="E4.5_Epiblast"),]
E4.5_Epi = E4.5_Epi[order(meta_new_2[which(index=="E4.5_Epiblast"),][,13],decreasing = T),]

E5.5_Epi = cell_TSS_5kb_N_meta_order [which(index=="E5.5_Epiblast"),]
E5.5_Epi = E5.5_Epi[order(meta_new_2[which(index=="E5.5_Epiblast"),][,13],decreasing = T),]

E6.5_Epi = cell_TSS_5kb_N_meta_order [which(index=="E6.5_Epiblast"),]
E6.5_Epi = E6.5_Epi[order(meta_new_2[which(index=="E6.5_Epiblast"),][,13],decreasing = T),]

E7.5_Epi = cell_TSS_5kb_N_meta_order [which(index=="E7.5_Epiblast"),]
E7.5_Epi = E7.5_Epi[order(meta_new_2[which(index=="E7.5_Epiblast"),][,13],decreasing = T),]


plot_data = rbind(E4.5_Epi ,E5.5_Epi ,E6.5_Epi ,E7.5_Epi )

annotation_row = data.frame(
  cell = factor(rep(c("E4.5_Epi", "E5.5_Epi", "E6.5_Epi","E7.5_Epi"), c(52,   69 ,    90 , 30 )))
)
rownames(annotation_row)=rownames(plot_data)

pheatmap(plot_data,cluster_cols = F,cluster_rows = F,
         show_rownames = F,annotation_row= annotation_row, color =
           colorRampPalette(  colors = c("#000000" ,"#FFFFFF","#808080","#008080"))  (4),fontsize_col = 5)


###############heatmap via ggplot
####################################################33
plot = cell_filtered_2
plot$cell = as.character(plot$cell)
index = paste(meta_new_2$stage,meta_new_2$lineage10x_2,sep="_")


#####E4.5_Epiblast
cell.use = meta_new_2$sample[which(index=="E4.5_Epiblast")]
cell.use = as.character(cell.use)
E4.5_Epiblast  <- plot[plot$cell %in% cell.use, ]

#####E4.5_Primitive_endoderm
cell.use = meta_new_2$sample[which(index=="E4.5_Primitive_endoderm ")]
cell.use = as.character(cell.use)
E4.5_Primitive_endoderm  <- plot[plot$cell %in% cell.use, ]

#####E5.5_Epiblast
cell.use = meta_new_2$sample[which(index=="E5.5_Epiblast")]
cell.use = as.character(cell.use)
E5.5_Epiblast  <- plot[plot$cell %in% cell.use, ]

#####E6.5_Epiblast
cell.use = meta_new_2$sample[which(index=="E6.5_Epiblast")]
cell.use = as.character(cell.use)
E6.5_Epiblast <- plot[plot$cell %in% cell.use, ]

#####E6.5_Primitive_Streak
cell.use = meta_new_2$sample[which(index=="E6.5_Primitive_Streak")]
cell.use = as.character(cell.use)
E6.5_Primitive_Streak <- plot[plot$cell %in% cell.use, ]

#####E7.5_Epiblast
cell.use = meta_new_2$sample[which(index=="E7.5_Epiblast")]
cell.use = as.character(cell.use)
E7.5_Epiblast <- plot[plot$cell %in% cell.use, ]


#####E7.5_Primitive_Streak
cell.use = meta_new_2$sample[which(index=="E7.5_Primitive_Streak")]
cell.use = as.character(cell.use)
E7.5_Primitive_Streak <- plot[plot$cell %in% cell.use, ]

###################################################################################################
###################################################################################################
#signal binning


Enhancer_filtered = cell_filtered_2
x = match(Enhancer_filtered$cell, meta$sample)
Enhancer_filtered_meta = meta[x,]
Enhancer_filtered_meta = cbind(Enhancer_filtered,Enhancer_filtered_meta)
Enhancer_filtered_meta = na.omit(Enhancer_filtered_meta)
write.table(Enhancer_filtered_meta,"TSS_5kb_filtered_meta.txt",sep="\t",quote = F)

library(windowscanr)
table(Enhancer_filtered_meta$lineage10x_2)
# Ectoderm           Endoderm           Epiblast       ExE_ectoderm
#3806               4461              29922                956
#Mesoderm Primitive_endoderm   Primitive_Streak  Visceral_endoderm
#9768               3019               7112               3326

win_size =200
win_step=200


#######bin signal according to both stage and lineage####################################################
staged_lineage = paste( Enhancer_filtered_meta$stage,Enhancer_filtered_meta$lineage10x_2,sep="_")

#table(staged_lineage)
#E4.5_Epiblast E4.5_Primitive_endoderm           E5.5_Epiblast
#6630                    3019                    7770
#E5.5_Visceral_endoderm           E6.5_Epiblast       E6.5_ExE_ectoderm
#2221                   12498                     956
#E6.5_Mesoderm   E6.5_Primitive_Streak  E6.5_Visceral_endoderm
#763                    4086                    1105
#E7.5_Ectoderm           E7.5_Endoderm           E7.5_Epiblast
#3806                    4461                    3024
#E7.5_Mesoderm   E7.5_Primitive_Streak
#9005                    3026

#############E4.5


start =35501032
de_sig = Enhancer_filtered_meta[which(staged_lineage =="E4.5_Epiblast"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E4.5_Epiblast = rol_pos

de_sig = Enhancer_filtered_meta[which(staged_lineage =="E4.5_Primitive_endoderm"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E4.5_Primitive_endoderm = rol_pos

rol_pos_merge_E4.5 = cbind(rol_pos[,1:5],
                           rol_pos_E4.5_Epiblast$rate_mean,
                           rol_pos_E4.5_Primitive_endoderm$rate_mean
                           
)
colnames(rol_pos_merge_E4.5) = c(colnames(rol_pos)[1:5],
                                 "E4.5_Epiblast",
                                 "E4.5_Primitive_endoderm"
)

#blue;red; green; purple
spline_int1 <- as.data.frame(spline(rol_pos_merge_E4.5$win_mid, rol_pos_merge_E4.5$E4.5_Epiblast))
spline_int2 <- as.data.frame(spline(rol_pos_merge_E4.5$win_mid, rol_pos_merge_E4.5$E4.5_Primitive_endoderm))

ggplot(rol_pos_merge_E4.5,aes(x=win_mid))+
  geom_line(data = spline_int1, aes(x = x, y = y),col="#377EB8")+
  geom_line(data = spline_int2, aes(x = x, y = y),col="#E41A1C")+
  theme_bw()+
  scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), limits=c(0,10000))+
  scale_y_continuous(name="Accessibility",
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))



#############E5.5
de_sig = Enhancer_filtered_meta[which(staged_lineage =="E5.5_Epiblast"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E5.5_Epiblast = rol_pos

de_sig = Enhancer_filtered_meta[which(staged_lineage =="E5.5_Visceral_endoderm"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E5.5_Visceral_endoderm = rol_pos

rol_pos_merge_E5.5 = cbind(rol_pos[,1:5],
                           rol_pos_E5.5_Epiblast$rate_mean,
                           rol_pos_E5.5_Visceral_endoderm$rate_mean
                           
)
colnames(rol_pos_merge_E5.5) = c(colnames(rol_pos)[1:5],
                                 "E5.5_Epiblast",
                                 "E5.5_Visceral_endoderm"
)



spline_int1 <- as.data.frame(spline(rol_pos_merge_E5.5$win_mid, rol_pos_merge_E5.5$E5.5_Epiblast))
spline_int2 <- as.data.frame(spline(rol_pos_merge_E5.5$win_mid, rol_pos_merge_E5.5$E5.5_Visceral_endoderm))


ggplot(rol_pos_merge_E5.5,aes(x=win_mid))+
  geom_line(data = spline_int1, aes(x = x, y = y),col="#377EB8")+
  geom_line(data = spline_int2, aes(x = x, y = y),col="#E41A1C")+
  theme_bw()+
  scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), limits=c(0,10000))+
  scale_y_continuous(name="Accessibility",
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))






#############E6.5
de_sig = Enhancer_filtered_meta[which(staged_lineage =="E6.5_Epiblast"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E6.5_Epiblast= rol_pos

de_sig = Enhancer_filtered_meta[which(staged_lineage =="E6.5_ExE_ectoderm"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E6.5_ExE_ectoderm = rol_pos


de_sig = Enhancer_filtered_meta[which(staged_lineage =="E6.5_Mesoderm"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E6.5_Mesoderm= rol_pos


de_sig = Enhancer_filtered_meta[which(staged_lineage =="E6.5_Primitive_Streak"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E6.5_Primitive_Streak= rol_pos

de_sig = Enhancer_filtered_meta[which(staged_lineage =="E6.5_Visceral_endoderm"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E6.5_Visceral_endoderm= rol_pos



rol_pos_merge_E6.5 = cbind(rol_pos[,1:5],
                           rol_pos_E6.5_Epiblast$rate_mean,
                           rol_pos_E6.5_Primitive_Streak$rate_mean,
                           rol_pos_E6.5_Visceral_endoderm$rate_mean,
                           rol_pos_E6.5_Mesoderm$rate_mean,
                           rol_pos_E6.5_ExE_ectoderm$rate_mean)
colnames(rol_pos_merge_E6.5) = c(colnames(rol_pos)[1:5],
                                 "E6.5_Epiblast",
                                 "E6.5_Primitive_Streak",
                                 
                                 "E6.5_Visceral_endoderm",
                                 "E6.5_Mesoderm",
                                 "E6.5_ExE_ectoderm")


spline_int1 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_Epiblast))
spline_int2 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_Primitive_Streak))
spline_int3 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_Visceral_endoderm))
spline_int4 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_Mesoderm))
spline_int5 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_ExE_ectoderm))

ggplot(rol_pos_merge_E6.5,aes(x=win_mid))+
  geom_line(data = spline_int1, aes(x = x, y = y),col="blue")+
  geom_line(data = spline_int2, aes(x = x, y = y),col="red")+
  geom_line(data = spline_int3, aes(x = x, y = y),col="green")+
  geom_line(data = spline_int4, aes(x = x, y = y),col="orange")+
  geom_line(data = spline_int5, aes(x = x, y = y),col="black")+
  theme_bw()+
  scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), limits=c(0,10000))+
  scale_y_continuous(name="Accessibility",
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))






#############E7.5
de_sig = Enhancer_filtered_meta[which(staged_lineage =="E7.5_Epiblast"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E7.5_Epiblast= rol_pos

de_sig = Enhancer_filtered_meta[which(staged_lineage =="E7.5_Ectoderm"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E7.5_Ectoderm= rol_pos


de_sig = Enhancer_filtered_meta[which(staged_lineage =="E7.5_Endoderm"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E7.5_Endoderm  = rol_pos


de_sig = Enhancer_filtered_meta[which(staged_lineage =="E7.5_Mesoderm"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E7.5_Mesoderm= rol_pos

de_sig = Enhancer_filtered_meta[which(staged_lineage =="E7.5_Primitive_Streak"),1:2]
pos = avereps(de_sig,ID=de_sig$Pos)
pos = matrix(as.numeric(pos), ncol=2,dimnames=dimnames(pos))
pos = as.data.frame(pos)
pos$Pos = pos$Pos - start
pos = na.omit(pos)
rol_pos = winScan(x=pos, position="Pos",values="rate",win_size=win_size,win_step=win_step,funs="mean")
rol_pos_E7.5_Primitive_Streak= rol_pos



rol_pos_merge_E7.5 = cbind(rol_pos[,1:5],
                           rol_pos_E7.5_Epiblast$rate_mean,
                           rol_pos_E7.5_Primitive_Streak$rate_mean,
                           rol_pos_E7.5_Ectoderm$rate_mean,
                           rol_pos_E7.5_Endoderm$rate_mean,
                           rol_pos_E7.5_Mesoderm$rate_mean)
colnames(rol_pos_merge_E7.5) = c(colnames(rol_pos)[1:5],
                                 "E7.5_Epiblast",
                                 "E7.5_Primitive_Streak",
                                 
                                 "E7.5_Ectoderm",
                                 "E7.5_Endoderm",
                                 "E7.5_Mesoderm")


spline_int1 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Epiblast))
spline_int2 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Primitive_Streak))
spline_int3 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Ectoderm))
spline_int4 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Endoderm))
spline_int5 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Mesoderm))

ggplot(rol_pos_merge_E7.5,aes(x=win_mid))+
geom_line(data = spline_int1, aes(x = x, y = y),col="blue")+
 geom_line(data = spline_int2, aes(x = x, y = y),col="red")+
 geom_line(data = spline_int3, aes(x = x, y = y),col="green")+
 geom_line(data = spline_int4, aes(x = x, y = y),col="orange")+
geom_line(data = spline_int5, aes(x = x, y = y),col="black")+
  theme_bw()+
  scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), limits=c(0,10000))+
  scale_y_continuous(name="Accessibility",
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))



spline_int1 <- as.data.frame(spline(rol_pos_merge_E4.5$win_mid, rol_pos_merge_E4.5$E4.5_Epiblast))
spline_int2 <- as.data.frame(spline(rol_pos_merge_E5.5$win_mid, rol_pos_merge_E5.5$E5.5_Epiblast))
spline_int3 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_Epiblast))
spline_int4 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Epiblast))


ggplot(rol_pos_merge_E4.5)+
  geom_line(data = spline_int1, aes(x = x, y = y,col="E4.5 Epiblast"),lwd=.8)+
  geom_line(data = spline_int2, aes(x = x, y = y,col="E5.5 Epiblast"),lwd=.8)+
  geom_line(data = spline_int3, aes(x = x, y = y,col="E6.5 Epiblast"),lwd=.8)+
  geom_line(data = spline_int4, aes(x = x, y = y,col="E7.5 Epiblast"),lwd=.8)+
  theme_bw()+
  scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),limits=c(0,10000))+
  scale_y_continuous(name="Accessibility",
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))+
  scale_color_manual(name = "Colors",
                     values = c("E4.5 Epiblast" = "#E41A1C", "E5.5 Epiblast" = "#4DAF4A","E6.5 Epiblast" = "orange", "E7.5 Epiblast" = "#377EB8"))




fill()spline_int2 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_Primitive_Streak))
spline_int3 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Primitive_Streak))
spline_int4 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Ectoderm))


ggplot(rol_pos_merge_E4.5)+
  geom_line(data = spline_int2, aes(x = x, y = y,col="E6.5 Primitive Streak"),lwd=.8)+
  geom_line(data = spline_int3, aes(x = x, y = y,col="E7.5 Primitive_Streak"),lwd=.8)+
  geom_line(data = spline_int4, aes(x = x, y = y,col="E7.5 Ectoderm"),lwd=.8)+
  theme_bw()+
  scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), limits=c(0,10000))+
  scale_y_continuous(name="Accessibility",
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))+
  scale_color_manual(name = "Colors",
                     values = c("E6.5 Primitive Streak" = "#E41A1C", "E7.5 Primitive_Streak" = "#4DAF4A","E7.5 Ectoderm" = "#377EB8"))


spline_int3 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_ExE_ectoderm))
spline_int1 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_Epiblast))


ggplot(rol_pos_merge_E4.5)+
  geom_line(data = spline_int3, aes(x = x, y = y,col="E6.5 Epiblast"),lwd=.8)+
  geom_line(data = spline_int1, aes(x = x, y = y,col="E6.5 ExE ectoderm"),lwd=.8)+
  
  scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), limits=c(0,10000))+
  scale_y_continuous(name="Accessibility",
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))+
  scale_color_manual(name = "Colors",
                     values = c("E6.5 Epiblast" = "#E41A1C","E6.5 ExE ectoderm" = "#4DAF4A"))






spline_int3 <- as.data.frame(spline(rol_pos_merge_E4.5$win_mid, rol_pos_merge_E4.5$E4.5_Epiblast))
spline_int1 <- as.data.frame(spline(rol_pos_merge_E4.5$win_mid, rol_pos_merge_E4.5$E4.5_Primitive_endoderm))


ggplot(rol_pos_merge_E4.5)+
  geom_line(data = spline_int3, aes(x = x, y = y,col="E4.5 Epiblast"),lwd=.8)+
  geom_line(data = spline_int1, aes(x = x, y = y,col="E4.5 Primitive endoderm"),lwd=.8)+
  scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000), limits=c(0,10000))+
  scale_y_continuous(name="Accessibility",
                     breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))+
  scale_color_manual(name = "Colors",
                     values = c("E4.5 Epiblast" = "#E41A1C","E4.5 Primitive endoderm" = "#377EB8"))




##plot gene isoform annotation
library(gggenes)
library(ggplot2)
genes = read.table("gene_exon.txt",header=T)

p1= ggplot(genes, aes(xmin = start, xmax = end, y = Gene,fill=Gene)) +
    geom_gene_arrow(arrowhead_width=grid::unit(.3, "mm"), arrowhead_height = grid::unit(.7, "mm"),
       arrow_body_height = grid::unit(1, "mm")) +theme_bw()+scale_fill_brewer(palette = "Set3")+
    scale_x_continuous(name="Position", limits=c(35501000,35511000),
          labels=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,
          35508000,35509000,35510000,35511000),
    breaks =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,
          35508000,35509000,35510000,35511000)
)

spline_int1 <- as.data.frame(spline(rol_pos_merge_E4.5$win_mid, rol_pos_merge_E4.5$E4.5_Epiblast))

p2=ggplot(rol_pos_merge_E4.5)+
   geom_line(data = spline_int1, aes(x = x, y = y,col="E4.5 Epiblast"),lwd=.8)+
   theme_bw()+
   scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
       labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
       limits=c(0,10000))+
   scale_y_continuous(name="Accessibility",
        breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))+
       scale_color_manual(name = "Colors",
       values = c("E4.5 Epiblast" = "#E41A1C", "E5.5 Epiblast" = "#4DAF4A","E6.5 Epiblast" = "orange",
       "E7.5 Epiblast" = "#377EB8"))

p3= ggplot(E4.5_Epi)+
   geom_point(aes(x=Pos,y=cell,col=as.factor(rate)),cex=.2,shape=22)+
   theme_classic()+
   theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
   scale_x_continuous(name="Position",
                      breaks=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(35501000,35511000))+
   scale_color_manual(values = c("grey","royalblue"))+labs(col="Methylation State")


p4=ggplot(E4.5_Epi,aes(x=Pos,col=as.factor(rate)))+geom_histogram(bins=100)+theme_classic()+
 scale_x_continuous(name="Position",
     breaks=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
     labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
     limits=c(35501000,35511000))+labs(col="Methylation State")+ylab("Depth")
#library(Seurat)
#CombinePlots(plots=list(p1,p2,p3,p4),ncol=1,rel_heights = c(1,2,2,1))

library(cowplot)
plot_grid (p1,p2,p3,p4,align="v",ncol=1,axis="tblr",rel_heights = c(2,5,1,2))




######################################################
spline_int1 <- as.data.frame(spline(rol_pos_merge_E4.5$win_mid, rol_pos_merge_E4.5$E4.5_Epiblast))
spline_int2 <- as.data.frame(spline(rol_pos_merge_E5.5$win_mid, rol_pos_merge_E5.5$E5.5_Epiblast))
spline_int3 <- as.data.frame(spline(rol_pos_merge_E6.5$win_mid, rol_pos_merge_E6.5$E6.5_Epiblast))
spline_int4 <- as.data.frame(spline(rol_pos_merge_E7.5$win_mid, rol_pos_merge_E7.5$E7.5_Epiblast))


#E5.5######################################################
p5=ggplot(rol_pos_merge_E5.5)+
 geom_line(data = spline_int2, aes(x = x, y = y,col="E5.5 Epiblast"),lwd=.8)+
 theme_bw()+
   scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(0,10000))+
   scale_y_continuous(name="Accessibility",
                      breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))+
   scale_color_manual(name = "Colors",
                      values = c("E4.5 Epiblast" = "#E41A1C", "E5.5 Epiblast" = "#4DAF4A","E6.5 Epiblast" = "orange",
                                 "E7.5 Epiblast" = "#377EB8"))

p6=ggplot(E5.5_Epiblast)+
   geom_point(aes(x=Pos,y=cell,col=as.factor(rate)),cex=.2,shape=22)+
   theme_classic()+
   theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
   scale_x_continuous(name="Position",
                      breaks=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(35501000,35511000))+
   scale_color_manual(values = c("grey","royalblue"))+labs(col="Methylation State")

p7=ggplot(E5.5_Epiblast,aes(x=Pos,col=as.factor(rate)))+geom_histogram(bins=100)+theme_classic()+
   scale_x_continuous(name="Position",
                      breaks=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(35501000,35511000))+labs(col="Methylation State")+ylab("Depth")


#E6.5######################################################
p8=ggplot(rol_pos_merge_E6.5)+
   geom_line(data = spline_int3, aes(x = x, y = y,col="E6.5 Epiblast"),lwd=.8)+
   theme_bw()+
   scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(0,10000))+
   scale_y_continuous(name="Accessibility",
                      breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))+
   scale_color_manual(name = "Colors",
                      values = c("E4.5 Epiblast" = "#E41A1C", "E5.5 Epiblast" = "#4DAF4A","E6.5 Epiblast" = "orange",
                                 "E7.5 Epiblast" = "#377EB8"))

p9=ggplot(E6.5_Epiblast)+
   geom_point(aes(x=Pos,y=cell,col=as.factor(rate)),cex=.2,shape=22)+
   theme_classic()+
   theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
   scale_x_continuous(name="Position",
                      breaks=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(35501000,35511000))+
   scale_color_manual(values = c("grey","royalblue"))+labs(col="Methylation State")

p10=ggplot(E6.5_Epiblast,aes(x=Pos,col=as.factor(rate)))+geom_histogram(bins=100)+theme_classic()+
   scale_x_continuous(name="Position",
                      breaks=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(35501000,35511000))+labs(col="Methylation State")+ylab("Depth")

#E7.5######################################################
p11=ggplot(rol_pos_merge_E7.5)+
   geom_line(data = spline_int4, aes(x = x, y = y,col="E7.5 Epiblast"),lwd=.8)+
   theme_bw()+
   scale_x_continuous(name="Position",breaks=c(0,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(0,10000))+
   scale_y_continuous(name="Accessibility",
                      breaks=c(0,0.2,0.4,0.6,0.8,1.0), limits=c(0,1))+
   scale_color_manual(name = "Colors",
                      values = c("E4.5 Epiblast" = "#E41A1C", "E5.5 Epiblast" = "#4DAF4A","E6.5 Epiblast" = "orange",
                                 "E7.5 Epiblast" = "#377EB8"))

p12=ggplot(E7.5_Epiblast)+
   geom_point(aes(x=Pos,y=cell,col=as.factor(rate)),cex=.2,shape=22)+
   theme_classic()+
   theme(axis.title.x=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())+
   scale_x_continuous(name="Position",
                      breaks=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(35501000,35511000))+
   scale_color_manual(values = c("grey","royalblue"))+labs(col="Methylation State")

p13=ggplot(E7.5_Epiblast,aes(x=Pos,col=as.factor(rate)))+geom_histogram(bins=100)+theme_classic()+
   scale_x_continuous(name="Position",
                      breaks=c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      labels =c(35501000,35502000,35503000,35504000,35505000,35506000,35507000,35508000,35509000,35510000,35511000),
                      limits=c(35501000,35511000))+labs(col="Methylation State")+ylab("Depth")


pdf("E4.5_Epi.pdf",height=9,width = 12)
plot_grid (p1,p2,p3,p4,align="v",ncol=1,axis="tblr",rel_heights = c(1,2,1,2))
dev.off()

pdf("E5.5_Epi.pdf",height=7.5,width = 12)
plot_grid (p5,p6,p7,align="v",ncol=1,axis="tblr",rel_heights = c(2,1,2))
dev.off()


pdf("E6.5_Epi.pdf",height=7.5,width = 12)
plot_grid (p8,p9,p10,align="v",ncol=1,axis="tblr",rel_heights = c(2,1,2))
dev.off()

pdf("E7.5_Epi.pdf",height=7.5,width = 12)
plot_grid (p11,p12,p13,align="v",ncol=1,axis="tblr",rel_heights = c(2,1,2))
dev.off()

plot_grid (p5,p6,p7,align="v",ncol=1,axis="tblr",rel_heights = c(2,1,2))
plot_grid (p8,p9,p10,align="v",ncol=1,axis="tblr",rel_heights = c(2,1,2))
plot_grid (p11,p12,p13,align="v",ncol=1,axis="tblr",rel_heights = c(2,1,2))






setwd("~/becherlab_shares/People/Galli/PROJECTS/MS/DMF/DMF_paper/CellCNN_data/CytoNorm/Normalized/")
#saveRDS(data_all, '062020_data_all_umap.RDS')
data_all <- readRDS('062020_data_all.RDS')
#saveRDS(temp$manual_labels, '062020_manual_labels.RDS')
manual_labels <- readRDS('062020_manual_labels.RDS')
data_all$manual_labels <- manual_labels
dim(data_all)
data_all$cell_id <- 1:21942263

#data_all$manual_labels <- temp$manual_labels
###Subsample per condition for umap
data_umap <- data_all[,c(make.names(intersect(panel_surf$Antigen[panel_surf$uneven =='n'], panel_ics$Antigen[panel_ics$uneven == 'y'])),'pz','manual_labels') ]
data_umap <- data_umap[data_umap[,'manual_labels']!=3 &data_umap[,'manual_labels']!=5 &
                         data_umap[,'manual_labels']!=6 &data_umap[,'manual_labels']!=7 &
                         data_umap[,'manual_labels']!=8, ]



cell_tsne <- matrix(, nrow = 1000*(length(levels(data_umap$pz))), ncol = ncol(data_umap))
set.seed(123)
data_list <- split.data.frame(data_umap, list(data_umap$pz)) 
str(data_list)
sub_data_list <- list()
for(i in 1:max(data_umap$pz)) {
  if (nrow(data_list[[i]]) >= 1000) {
    set.seed(123)
    sub_data_list[[i]] <- data_list[[i]][sample(1:nrow(data_list[[i]]), 1000, replace = F),]
  } else {sub_data_list[[i]] <- data_list[[i]]
  }}
cell_tsne <- rbind.fill(sub_data_list)
remove(data_list)


ex_new <- ggplot(data.frame(cell_tsne[cell_tsne[,'pz'] >94,]), aes(x = CD8, y = CD4)) +
  geom_hex(bins = 100) +
  xlim(-0.2, 1.5) +
  ylim(-0.2, 1.5) +
  coord_fixed(ratio = 1) +
  scale_fill_gradientn(colours = jet.colors(100), trans = "sqrt") +
  theme_facs
ex_new

## Load tsne
tsne <- readRDS('tsne_multipanel.RDS')


# prepare data for UMAP

### Select clustering columns
library(umap)

data_rtsne <- cell_tsne[cell_tsne[,'CD3'] >0  , colnames(data_umap)[-c(18:19)]]
custom.settings = umap.defaults
custom.settings$verbose <- T
custom.settings$n_components <- 2
custom.settings$n_neighbors <- 25
custom.settings$random_state <- 123
custom.settings$n_epochs <- 300
custom.settings$min_dist <- 0.6
data_umap <-  umap(data_rtsne, custom.settings)


# prepare the tSNE data
tsne <- as.data.frame(data_umap$layout)
colnames(tsne) <- c("tSNE1", "tSNE2")
#saveRDS(tsne, 'tsne_multipanel.RDS')



# plot tSNE black
t1 <- ggplot(tsne, aes(x = tSNE1, y = tSNE2)) +
  geom_point(size = 0.5) +
  coord_fixed(ratio = 1) +
  theme_tsne 
t1



# save the plots as png
ggsave(filename = "OUTPUT/tsne_final.png", plot = t1, scale = 0.5)
ggsave(filename = "OUTPUT/gm_tsne_final.pdf", plot = t1, scale = 0.5)


data_umap <- data_all
###Subsample per condition for umap
cell_tsne <- matrix(, nrow = 1000*(length(levels(data_umap$pz))), ncol = ncol(data_umap))
set.seed(123)
data_list <- split.data.frame(data_umap, list(data_umap$pz)) 
str(data_list)
sub_data_list <- list()
for(i in 1:max(data_umap$pz)) {
  if (nrow(data_list[[i]]) >= 1000) {
    set.seed(123)
    sub_data_list[[i]] <- data_list[[i]][sample(1:nrow(data_list[[i]]), 1000, replace = F),]
  } else {sub_data_list[[i]] <- data_list[[i]]
  }}
cell_tsne <- rbind.fill(sub_data_list)
remove(data_list)


# prepare the expression data
data_plot <- cbind(cell_tsne[cell_tsne[,'CD3'] >0  ,], tsne)
head(data_plot)
data.kopf <- merge(data_plot, md, by = 'pz')
head(data.kopf)
saveRDS(data.kopf, 'density.raw.RDS')

data.kopf
### Plot FlowSOM overaly
data_dff <- data.kopf
data_dff[data_dff[,"manual_labels"] == 1,"manual_labels"] <- "CD4" 
data_dff[data_dff[,"manual_labels"] == 2,"manual_labels"] <- "CD8" 
data_dff[data_dff[,"manual_labels"] == 3,"manual_labels"] <- "gdT.cells" 
data_dff[data_dff[,"manual_labels"] == 4,"manual_labels"] <- "B.cell" 
data_dff[data_dff[,"manual_labels"] == 5,"manual_labels"] <- "NK" 
#data_dff[data_dff[,"manual_labels"] == 6,"manual_labels"] <- "Treg" 
data_dff[data_dff[,"manual_labels"] == 6,"manual_labels"] <- "Myeloid" 
data_dff$manual_labels <- factor(data_dff$manual_labels, levels = c("CD4","CD8", "gdT.cells", "B.cell", "NK", 'Treg',  "Myeloid")) 
#data_dff$tpc <- factor(data_dff$diagnosis2, levels = c('HC','CIS','MS','NINDC','INDC')) 

# add Treg in color
names(db_lin)[8] <- 'Treg'
#dev.off()
t3 <- ggplot(data_dff, aes(x = tSNE1, y = tSNE2, color = manual_labels)) +
  #geom_point(size = 0.05) +
  geom_density_2d( size = 0.25, bins = 30) + 
  #geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  #scale_x_discrete(limits = c("CD4","CD8", "TCRgd", "B cell", "NK", "Myeloid"))+
  scale_color_manual(values = db_lin) +
  scale_fill_manual(values = db_lin) +
  #scale_colour_gradientn(colours = jet.colors(100), limits = c(0,6)) +
  theme_facs+
  #facet_grid(rows =  vars(responder), cols = vars(tpc)) +
  # facet_wrap('timepoint') +
  theme_tsne +
  xlim(min(data_dff$tSNE1)-1,max(data_dff$tSNE1)+1 )+
  theme(plot.background = element_rect(color="white"),
        title = element_text(size = rel(1.1)),
        legend.position = "right",
        legend.key.size = unit(5,"point"),
        legend.key.width = unit(5,"point"),
        plot.title = element_text(size = rel(1.1)),
        strip.text = element_text(size = rel(1.1)),
        strip.background = element_blank(),
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 2)))
t3 

# save the plots as png
ggsave(filename = "Memory_cluster_map_density.png", plot = t3, 
       scale = 2.5, width = 4*1.65, height = 4*3.5, units = c("in"))


# plot tSNEs with expression overlayed
plot_cols <- make.names(intersect(panel_surf$Antigen[panel_surf$uneven =='n'], panel_ics$Antigen[panel_ics$uneven == 'y']))
data_melt <- melt(data.kopf[,c(plot_cols, "tSNE1", "tSNE2")], id.vars = c("tSNE1", "tSNE2"))
data_melt$variable <- factor(data_melt$variable, levels = levels(data_melt$variable)[order(levels(data_melt$variable))])
t2 <- ggplot(na.omit(data_melt), aes(x = tSNE1, y = tSNE2, color = value)) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = jet.colors(100), limits = c(-0.2,1.2)) +
  facet_wrap(~ variable, nrow = 3, scales = "fixed") +
  theme_tsne
t2


# save the plots as png
ggsave(filename = "lymph_tsne_expression_final.png", t2, scale = 0.5, 
       width = 6*2, height = 3*2.3, units = c("in"))

################# cellCNN differences
##surf
fm <- ddply(subset(data_all,pz >93), .(pz), function(x)
{df <- data.frame(x) 
(nrow(df[df[,"filter_2_binary"] > 0.25,])/nrow(df[df[,"pz"],]))*100 }
)
colnames(fm)[2] <- "freq"
fm
###ics
fm <- ddply(subset(data_all,pz <=93), .(pz), function(x)
{df <- data.frame(x) 
(nrow(df[df[,"filter_0_continuous"] > 0.25,])/nrow(df[df[,"pz"],]))*100 }
)
colnames(fm)[2] <- "freq"
fm$freq


fm <- merge( md_cn ,fm, by = 'pz')
fm_ord <- fm[ order(fm$tpc) , ]
T2T1 <- (fm_ord$freq[32:62]-fm_ord$freq[1:31])/fm_ord$freq[1:31]
T3T2 <- (fm_ord$freq[63:93]-fm_ord$freq[32:62])/fm_ord$freq[32:62]

# 
# T2T1 <- (fm_ord$freq[32:62]-fm_ord$freq[1:31])
# T3T2 <- (fm_ord$freq[63:93]-fm_ord$freq[1:31])

# 
# T2T1 <- (fm_ord$freq[1:31]-fm_ord$freq[32:62])/fm_ord$freq[1:31]
# T3T2 <- (fm_ord$freq[1:31]-fm_ord$freq[63:93])/fm_ord$freq[1:31]

colnames(fm)
freq <- c(T2T1, T3T2)
df_plot <- cbind(fm_ord[1:62,c(1:3,5,9,10,13)], freq)

# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot)[which(colnames(df_plot)%in% 'freq')])

# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
fmm
order <- as.numeric(rownames(fmm[order(fmm$median, decreasing = T),]))

# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)

# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[order])


# melt for ggplot
df_melt <- melt(df_plot, measure.vars = colnames(df_plot) %in% 'freq')
head(df_melt)


# calculate summary frequencies
fmm <- ddply(df_melt, .(variable), function(x) {
  c(median = median(x$value, na.rm = T),
    sem = b.median(x$value, 1000)
  )}) 
fmm



# reorder according to frequency
fmm <- fmm[order(fmm[,"median"], decreasing = T),]
fmm



# define things for plotting
dodge <- position_dodge(width = 0.9)
limits <- aes(ymax = fmm$median + fmm$sem,
              ymin = fmm$median - fmm$sem)



# reorder
df_melt$variable <- factor(df_melt$variable, levels = levels(df_melt$variable)[fmm$variable])
df_melt_2 <- df_melt
#df_melt_2 <- merge(df_melt, md , by = 'gate_source')

# plot boxplots
df_melt_2$diagnosis.spec <- factor(df_melt_2$diagnosis.spec, levels = levels(df_melt_2$diagnosis.spec)[c(2,5,4,3,6,1)])
symlog_trans <- function(base = 10, thr = 1, scale = 1){
  trans <- function(x)
    ifelse(abs(x) < thr, x, sign(x) * 
             (thr + scale * suppressWarnings(log(sign(x) * x / thr, base))))
  
  inv <- function(x)
    ifelse(abs(x) < thr, x, sign(x) * 
             base^((sign(x) * x - thr) / scale) * thr)
  
  breaks <- function(x){
    sgn <- sign(x[which.max(abs(x))])
    if(all(abs(x) < thr))
      pretty_breaks()(x)
    else if(prod(x) >= 0){
      if(min(abs(x)) < thr)
        sgn * unique(c(pretty_breaks()(c(min(abs(x)), thr)),
                       log_breaks(base)(c(max(abs(x)), thr))))
      else
        sgn * log_breaks(base)(sgn * x)
    } else {
      if(min(abs(x)) < thr)
        unique(c(sgn * log_breaks()(c(max(abs(x)), thr)),
                 pretty_breaks()(c(sgn * thr, x[which.min(abs(x))]))))
      else
        unique(c(-log_breaks(base)(c(thr, -x[1])),
                 pretty_breaks()(c(-thr, thr)),
                 log_breaks(base)(c(thr, x[2]))))
    }
  }
  trans_new(paste("symlog", thr, base, scale, sep = "-"), trans, inv, breaks)
}
b2 <- ggplot(data = droplevels(subset(df_melt_2, exclude_ics != 'y')), aes(x = tpc, y = value, fill = variable,  label= patient)) +
  geom_hline(aes(yintercept=0), colour="#990000", linetype="dashed")+
  geom_point(size = 4, aes(colour = variable)) +
  geom_boxplot(lwd = 0.5, color = "white", outlier.shape = NA, aes(color = variable)) +
  stat_boxplot(geom = 'errorbar', width = 0.5, aes(colour = variable)) +
  #facet_wrap('responder.rad', scales = 'fix', ncol = 7) +
  scale_fill_manual(values = db1[4]) +
  scale_color_manual(values = db1[4]) +
  ylab('fold change')+
  #geom_text(aes(label=patient),hjust=0, vjust=0)+
  #ylim(-1,5) +
  #facet_zoom(ylim = c(-10, 10), zoom.data = ifelse(a <= 10, NA, FALSE))+
  #scale_y_continuous(labels = percent,labels = scientific)+
  scale_y_continuous(trans=symlog_trans(thr = 1.5)) + 
  theme_bar2
b2 


#+ scale_y_log10(limits = c(10e-3,10e+4), breaks=c(10e-4,10e-3,10e-2,0.1,1,10,100,1000,10000))


nrow(fm[complete.cases(fm),])
t <- (fm[complete.cases(fm),])
t <- merge(t, md, by = 'gate_source')
table(t$diagnosis)
##NINDC=29
##RRMS=31

# save
ggsave(filename = "Foldchange_sign1_ics.pdf", plot = b2, width = 1.15, height = 1.8, 
       scale = 2.5, useDingbats = F)

########################
###### STAT block ######
### Summary

# fig <- 'Fig_2B'
# 
# sum <- ddply(df_melt_2, .(variable, diagnosis) , function(x) {
#   c(median = median(x$value, na.rm = T), 
#     sem = b.median(x$value, 1000), 
#     sd = sd(x$value, na.rm = T), 
#     n = length(na.omit(x$value)),
#     max = max(x$value, na.rm = T),
#     min = min(x$value, na.rm = T)
#   )})
# sum
# write.xlsx(sum, 'STATS/Fig_2B_pop_boxplots_sum.xlsx')


#####stats
wres <- ddply(df_melt_2, .(variable), function(x) {
  t <- pairwise.wilcox.test(data = x, x$value, g = x$tpc, p.adjust.method = "BH", paired = F, correct = T, exact = F)
  with(t, data.frame(method, t(as.vector(p.value)))) })

wres
wres <- melt(wres, id.vars = c("variable",'method' ), measure.vars = colnames(wres)[-c(1,2)])

###decide the * pValue system
wres[,4] <- round(wres[,4], 4)
t <- vector()
for(i in 1:nrow(wres)) {
  if (wres$value[i] <= 0.01)  (t[i] <- "***")
  else if (wres$value[i] <= 0.05)  (t[i] <- "**")
  else if (wres$value[i] <= 0.10)  (t[i] <- "*")  
  else (t[i] <- "")}
wres$stars <- t
wres

write.xlsx(wres, 'STATS/Fig_2B_pop_boxplots_stat.xlsx')


addWorksheet(wb, paste(fig,'_sum',sep = ''))
addWorksheet(wb, paste(fig,'_stat',sep = ''))
writeData(wb, sheet = paste(fig,'_sum',sep = ''), x = sum)
writeData(wb, sheet = paste(fig,'_stat',sep = ''), x =wres)

wres <- ddply(df_melt_2, .(variable), function(x) {
  t2 <- coin::wilcox_test(value ~ tpc, data = x, distribution = "exact", conf.int = T)
  p <- pvalue(t2)
  z <- t2@statistic@teststatistic
  r <- as.numeric(z)/sqrt(nrow(x))
  d <- 2*r/sqrt(1-r^2)
  y <- x[!is.na(x$value),]
  cles <- cles.fnc(variable = "value", group = "diagnosis", baseline = "RRMS", 
                   data = y, print = F)
  data.frame(p.value = p, z = -z, r = -r, r2 = r^2, d = -d, cles)
})
wres

write.xlsx(wres, 'STATS/Fig_2B_pop_boxplots_effect_size.xlsx')

addWorksheet(wb, paste(fig,'_eff_size',sep = ''))
writeData(wb, sheet = paste(fig,'_eff_size',sep = ''), x =wres)

##############################################################################

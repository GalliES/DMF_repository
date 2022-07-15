####################
#####data_prep######
####################

library(dplyr)
library(readr)
setwd("/Immunology_shares/Becherlab/People/Galli/PROJECTS/MS/DMF/DMF_paper/CellCNN_data/sigpop-ics")
list.files(path = getwd(), pattern = ".csv$" )

df <- list.files(path=getwd(), full.names = TRUE) %>% 
  lapply(read_csv) 

df1 <- bind_rows(df)
head(df1)

class(df)
str(df)
head(df[1,12])

df1 <- as.data.frame(df1) 
head(df1)

# make a gate.source vector
dim.vec <- as.numeric(lapply(df, nrow))

gate.source <- as.vector(x = NULL)
for(i in 1:length(dim.vec)) {temp.source <- rep(i, dim.vec[i])
gate.source <- c(gate.source, temp.source)}
table(gate.source)
str(gate.source)
pz <- gate.source
data_all_ics <- cbind(df1, pz, gate.source)

data_all_ics <- data_all_ics[!duplicated(data_all_ics[,colnames(data_all_ics)[1:34] ]), ]


setwd("/Immunology_shares/Becherlab/People/Galli/PROJECTS/MS/DMF/DMF_paper/CellCNN_data/sigpop-surf/")
list.files(path = getwd(), pattern = ".csv$" )

df <- list.files(path=getwd(), full.names = TRUE) %>% 
  lapply(read_csv) 

df2 <- bind_rows(df)
head(df2)

class(df)
str(df)
head(df[1,12])

df2 <- as.data.frame(df2) 
head(df1)


# make a gate.source vector
dim.vec <- as.numeric(lapply(df, nrow))

gate.source <- as.vector(x = NULL)
for(i in 1:length(dim.vec)) {temp.source <- rep(i, dim.vec[i])
gate.source <- c(gate.source, temp.source)}
table(gate.source)
str(gate.source)
pz <- gate.source+93
data_all_surf <- cbind(df2,pz, gate.source)
colnames(data_all_surf)
data_all_surf <- data_all_surf[!duplicated(data_all_surf[,colnames(data_all_surf)[1:33] ]), ]
dim(data_all_surf)

data_all <- rbind.fill(data_all_ics, data_all_surf) 
dim(data_all)
head(data_all)
tail(data_all)
#table(data_all$id, data_all$manual_labels)
dim(data_all)
data_all$cell_id <- 1:21942263 





####################
#####data norm######
####################


setwd("/Immunology_shares/Becherlab/People/Galli/PROJECTS/MS/DMF/DMF_paper/")

# read in panel
panel_ics <- read.xlsx("ICS/FILES/panel_ics.xlsx")
# read in panel
panel_surf <- read.xlsx("SURF/FILES/panel_surf.xlsx")


# read meta_cytonorm
setwd("/Immunology_shares/Becherlab/People/Galli/PROJECTS/MS/DMF/DMF_paper/CellCNN_data/")
md_cn <- read.xlsx("05.2020_meta_data_CytoNorm_dmf_combined.xlsx")
dim(md_cn)

setwd("/Immunology_shares/Becherlab/People/Galli/PROJECTS/MS/DMF/DMF_paper/CellCNN_data/CytoNorm/")
data_umap <- data_all[,c(make.names(intersect(panel_surf$Antigen[panel_surf$uneven =='n'], panel_ics$Antigen[panel_ics$uneven == 'y'])),'gate.source','pz', 'cell_id') ]
dim(data_umap)

#BiocManager::install('matrixStats')
#BiocManager::install('latticeExtra')
#saveRDS(data_all,'data_all.RDS')
data_umap[data_umap[,'gate.source'] ==195,]
str(data_umap)
data_umap$cell_id <- as.numeric(data_umap$cell_id)
head(data_umap)


for(i in 1:max(data_umap[,'pz'])) {
  data_x <- data.matrix(data_umap)
  ff <-flowFrame(data_x[data_x[,'pz'] == i, ])
  write.FCS(ff, filename = paste(i,md_cn$diagnosis[md_cn$pz == i],
                                 md_cn$stim[md_cn$pz == i],
                                 "basel.fcs", sep = "_"))
}




library(CytoNorm)
files <- list.files(getwd(), pattern = "fcs$")
batch <- stringr::str_match(files, "_([0-2])*_basel.fcs")[, 2] ### batch stim
type <- as.character(c(rep(2,198)))
type[which(stringr::str_match(files, "HC*")[, 1] == 'HC')] <- 1
type


data <- data.frame(File = files,
                   Path = file.path(getwd(), files),
                   Type = type,
                   Batch = paste(stringr::str_match(files, "_([0-2])*_basel.fcs")[, 2],'basel.fcs',sep = '_') ,
                   stringsAsFactors = FALSE)

# data <- data.frame(File = files,
#                    Path = file.path(getwd(), files),
#                    Type = as.character(c(rep(1,192),2,1,1,2,rep(1,6))),
#                    Batch = paste(type_norm,'basel.fcs',sep = '_') ,
#                    stringsAsFactors = FALSE)
data$Type <- c("1" = "Train", "2" = "Validation")[data$Type]

train_data <- dplyr::filter(data, Type == "Train")
validation_data <- dplyr::filter(data, Type == "Validation")


ff <- flowCore::read.FCS(data$Path[106], truncate_max_range = F)
head(ff)

channels <- flowCore::colnames(ff)[c(1:17)]


library(CytoNorm)
transformList <- flowCore::transformList(channels,
                                         cytofTransform)
transformList.reverse <- flowCore::transformList(channels,
                                                 cytofTransform.reverse)

fsom <- prepareFlowSOM(train_data$Path,
                       channels,
                       nCells = 6000,
                       FlowSOM.params = list(xdim = 10,
                                             ydim = 10,
                                             nClus = 20,
                                             scale = FALSE),
                       transformList = transformList,
                       seed = 1)

cvs <- testCV(fsom, cluster_values = c(5:20)) 

cvs$pctgs$`20`


model <- CytoNorm.train(files = train_data$Path,
                        labels = train_data$Batch,
                        channels = channels,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = 6000, 
                                              xdim = 10,
                                              ydim = 10,
                                              nClus = 12,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)


CytoNorm.normalize(model = model,
                   files = validation_data$Path,
                   labels = validation_data$Batch,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "Normalized",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE)





####################
#####data asse######
####################


setwd("/Immunology_shares/Becherlab/People/Galli/PROJECTS/MS/DMF/DMF_paper/CellCNN_data/CytoNorm/Normalized/")
library(ncdfFlow)
library(gtools)

## read in the files 
files = list.files(path = getwd(), pattern = ".fcs$" )
mixedsort(files)

# myFiles.sorted <- sort(files) # this gives the alphabetic sorting, not what you want
# myFiles.sorted
# # split between the part that comes before the numerics and the "1.img" etc.--adjust appropriately
# split <- strsplit(myFiles.sorted, "Norm_") 
# # strip the "1.img" etc such that only the numeric part is left
# # turn the characters in numeric
# split <- as.numeric(sapply(split, function(x) x <- sub("_basel.fcs", "", x[2])))
# # not you can sort, by using order, that gives the original filenames, ordered on the numeric part of the filename
# myFiles.correct.order <- myFiles.sorted[order(split)]

fs  <- read.ncdfFlowSet(files = mixedsort(files),
                        transformation = F,
                        phenoData = ,
                        truncate_max_range = F)

ff <- flowCore::read.FCS(list.files(path = getwd(), pattern = ".fcs$" )[[1]], truncate_max_range = F)
head(ff)
names(fs)
sampleNames(fs)

# combine data into a matrix
data_norm <- fsApply(fs, exprs)
data_norm <- data.frame(data_norm)
# dim(data_norm)
# head(data_norm)
# tail(data_norm)
# data_norm <- data_norm %>% arrange(cell_id)
# min(data_norm[,'cell_id'])
# max(data_norm[,'cell_id'])

data_norm$cell_id <- data_all$cell_id


# normalize to have everything from 0 to 1
data_norm <- data.matrix(data_norm)
data.trans.new <- data_norm
per.vector <- apply(data.trans.new, 2, function(x) quantile(x, 0.9999, names = F))
per.vector
data.trans.new <- t(t(data.trans.new) / as.numeric(per.vector))
data.trans.new[,c("gate.source",'pz', "cell_id")] <- data_norm[,c("gate.source",'pz', "cell_id")]
head(data.trans.new)
head(data_norm)
max(data_norm[,'CD4'])
max(data.trans.new[,'CD4'])

#data.trans.new[,'IL.9'] <- t(t(data.trans.new[,'IL.9']) / as.numeric(2.2))



# # make example plot
# data.trans.new <- data.frame(data.trans.new)
# f <- subset(data.trans.new, gate.source > 70)
# ex_new <- ggplot(f, aes(x = CD20, y = CD4)) +
#    geom_hex(bins = 100) +
#    xlim(-0.2, 1.5) +
#    ylim(-0.2, 1.5) +
#    coord_fixed(ratio = 1) +
#    scale_fill_gradientn(colours = jet.colors(100), trans = "sqrt") +
#    theme_facs
# ex_new



#### combina
data_base <- data_all
data_norm <- data.trans.new
dim(data_base)
dim(data_norm)

head(data_base[,-c(which(colnames(data_base) %in% colnames(data_norm)[-c(18:19)] ))])
#data_all <- merge(data_base[,-c(which(colnames(data_base) %in% colnames(data_norm)[-c(20)] ))], data_norm, by = 'cell_id')
data_all <- cbind(data_base[,-c(which(colnames(data_base) %in% colnames(data_norm) ))], data_norm)
dim(data_all)
temp <- data_all
data_all[data_all > 1.51] <- 1
data_all[,c(55:57)] <- temp[,c(55:57)]
colnames(data_all)
table(data_all[,'pz'])

####################
##### FlowSOM ######
####################


#saveRDS(data_all, '052020_data_all.RDS')
#data_all <- readRDS('052020_data_all.RDS')
ex_new <- ggplot(data.frame(subset(data_all, pz >=94)), aes(x = CD4, y = FOXP3)) +
  geom_hex(bins = 100) +
  xlim(-0.2, 1.5) +
  ylim(-0.2, 1.5) +
  coord_fixed(ratio = 1) +
  scale_fill_gradientn(colours = jet.colors(100), trans = "sqrt") +
  theme_facs
ex_new

#data_all[data_all[,'pz']>=94, c("CD8","CD4") ] <- data_all[data_all[,'pz']>=94, c("CD8","CD4") ]*1.2
#data_all[data_all[,'pz']>=94, c("FOXP3","CD25") ] <- data_all[data_all[,'pz']>=94, c("FOXP3","CD25") ]/2




####FlowSOM
data <- droplevels(subset(data_all, pz>93))
ex_new <- ggplot(data.frame(subset(data, pz <194)), aes(x = CD4, y = CD8)) +
  geom_hex(bins = 100) +
  xlim(-0.2, 1.5) +
  ylim(-0.2, 1.5) +
  coord_fixed(ratio = 1) +
  scale_fill_gradientn(colours = jet.colors(100), trans = "sqrt") +
  theme_facs
ex_new

data[,c("FOXP3","CD25")] <- data[,c("FOXP3","CD25")]*2
ff_new <- flowFrame(exprs = data.matrix(data), desc = list(FIL = 1))
ff_new


#clustering_cols<- c("CD4",  "CD8", "CD56", "CD14", "CXCR5",  "TCRgd", "CD20")
(clustering_cols <- make.names(as.vector(panel_surf$Antigen[ panel_surf$uneven != 'y' & panel_surf$Antigen %in% colnames(data_all) ])))

# run FlowSOM (with set.seed for reproducibility)
set.seed(123)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = F, scale = F, compensate = F)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clustering_cols)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
labels <- out_fSOM$map$mapping[,1]



# choose the plot_cols for the heatmap
(plot_cols <- make.names(as.vector(panel_surf$Antigen[ panel_surf$uneven != 'y' & panel_surf$Antigen %in% colnames(data_all) ])))
#plot_cols <- c("CD3", "TCRgd", "CD4", "CCR2", "CD8", "CD7", "CD45RA", "CD16", "CD19",  "CD56", "CD14",  "CD20", "CCR7", "GM.CSF", "IFN.g")
(plot_cols <- plot_cols[order(plot_cols)])

table(labels)

# make labels heatmap
heat_mat <- matrix(NA, nrow = 100, ncol = length(plot_cols))
for(i in 1:100) {
  temp_mat <- data[labels == i, plot_cols]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}



# rename
rownames(heat_mat) <- 1:100
colnames(heat_mat) <- plot_cols 



# set the color scheme and plot heatmap
breaks <- seq(0, 1, by = 0.111111111)
white.black <- colorRampPalette(c("white", "black"))(n = 9)



# make a cluster heatmap
heatmap.2(heat_mat, 
          scale = "none",
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks)



# consensus clustering to determine how many clusters are there
set.seed(123)
results <- ConsensusClusterPlus(t(heat_mat), maxK = 30, reps = 1000, pItem = 0.8, 
                                pFeature = 1,  clusterAlg = "hc", verbose = F, 
                                distance = "euclidean", seed = 123, plot = "png", 
                                writeTable = T)

# perform actual metaclustering
choosen_k <- 20
set.seed(123)
out_meta <- metaClustering_consensus(out_fSOM$map$codes, k = choosen_k) 
pop_labels <- out_meta[labels]

(plot_cols <- make.names(as.vector(panel$Antigen[panel$Category == "lineage" & panel$Antigen %in% colnames(data)])))
plot_cols <- make.names(panel$Antigen)
plot_cols <- clustering_cols
(plot_cols <- plot_cols[order(plot_cols)])

# go through all clusters and calculate mean for every channel
set.seed(123)
heat_mat <- matrix(NA, nrow = choosen_k, ncol = length(plot_cols))
for(i in 1:choosen_k) {
  temp_mat <- data[pop_labels == i, plot_cols]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}



# rename
rownames(heat_mat) <- paste("cluster", 1:choosen_k, sep = "")
colnames(heat_mat) <- plot_cols
dim(heat_mat)


#dev.off()
# make a cluster heatmap
heatmap.2(heat_mat, 
          scale = "none",
          #margins = c(5,5),
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks)


table(pop_labels)

CD4 <- c(12,16,20)
CD8 <- c(1,2,9)
TCRgd <- c(7)
B.cells <- c(14,18,19,17)
PC <- c()
NK <- c(3)
Treg <- c(13) #Treg 4
Myeloid <- c(4,5,6,8,10,11,15)
DC <- c()


# give them the right metacluster number
th_labels <- rep(0, length(pop_labels))
th_labels[pop_labels %in% CD4] <- 1
th_labels[pop_labels %in% CD8] <- 2
th_labels[pop_labels %in% TCRgd] <- 3
th_labels[pop_labels %in% B.cells] <- 4
th_labels[pop_labels %in% NK] <- 5
th_labels[pop_labels %in% Treg] <- 6 #Treg
th_labels[pop_labels %in% Myeloid] <- 7
#th_labels[pop_labels %in% cd20t] <- 8

data[,"manual_labels"] <- th_labels
data[,c("FOXP3","CD25")] <- data[,c("FOXP3","CD25")]/2

head(data)
#data[,c("CD8","CD4")] <- data[,c("CD8","CD4")]/1.2
data_surf <- data




####FlowSOM
#metto il peso su CD4 e CD8
data <- droplevels(subset(data_all, pz<= 93))
#data[,c("CD8","CD4")] <- data[,c("CD8","CD4")]*1.2
ff_new <- flowFrame(exprs = data.matrix(data), desc = list(FIL = 1))
ff_new


#clustering_cols<- c("CD4",  "CD8", "CD56", "CD14", "CXCR5",  "TCRgd", "CD20")
(clustering_cols <- make.names(as.vector(panel_ics$Antigen[ panel_ics$uneven != 'n' & panel_ics$Antigen %in% colnames(data_all) ])))

# run FlowSOM (with set.seed for reproducibility)
set.seed(123)
out_fSOM <- FlowSOM::ReadInput(ff_new, transform = F, scale = F, compensate = F)
out_fSOM <- FlowSOM::BuildSOM(out_fSOM, colsToUse = clustering_cols)
out_fSOM <- FlowSOM::BuildMST(out_fSOM)
labels <- out_fSOM$map$mapping[,1]



# choose the plot_cols for the heatmap
(plot_cols <- make.names(as.vector(panel_ics$Antigen[ panel_ics$uneven != 'n' & panel_ics$Antigen %in% colnames(data_all) ])))
#plot_cols <- c("CD3", "TCRgd", "CD4", "CCR2", "CD8", "CD7", "CD45RA", "CD16", "CD19",  "CD56", "CD14",  "CD20", "CCR7", "GM.CSF", "IFN.g")
(plot_cols <- plot_cols[order(plot_cols)])

table(labels)

# make labels heatmap
heat_mat <- matrix(NA, nrow = 100, ncol = length(plot_cols))
for(i in 1:100) {
  temp_mat <- data[labels == i, plot_cols]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}



# rename
rownames(heat_mat) <- 1:100
colnames(heat_mat) <- plot_cols 



# set the color scheme and plot heatmap
breaks <- seq(0, 1, by = 0.111111111)
white.black <- colorRampPalette(c("white", "black"))(n = 9)



# make a cluster heatmap
heatmap.2(heat_mat, 
          scale = "none",
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks)



# consensus clustering to determine how many clusters are there
set.seed(123)
results <- ConsensusClusterPlus(t(heat_mat), maxK = 30, reps = 1000, pItem = 0.8, 
                                pFeature = 1,  clusterAlg = "hc", verbose = F, 
                                distance = "euclidean", seed = 123, plot = "png", 
                                writeTable = T)

# perform actual metaclustering
choosen_k <- 20
set.seed(123)
out_meta <- metaClustering_consensus(out_fSOM$map$codes, k = choosen_k) 
pop_labels <- out_meta[labels]

(plot_cols <- make.names(as.vector(panel$Antigen[panel$Category == "lineage" & panel$Antigen %in% colnames(data)])))
plot_cols <- make.names(panel$Antigen)
plot_cols <- clustering_cols
(plot_cols <- plot_cols[order(plot_cols)])

# go through all clusters and calculate mean for every channel
set.seed(123)
heat_mat <- matrix(NA, nrow = choosen_k, ncol = length(plot_cols))
for(i in 1:choosen_k) {
  temp_mat <- data[pop_labels == i, plot_cols]
  heat_mat[i,] <- as.matrix(apply(temp_mat, 2, function (x) mean(x, na.rm = T)))
}



# rename
rownames(heat_mat) <- paste("cluster", 1:choosen_k, sep = "")
colnames(heat_mat) <- plot_cols
dim(heat_mat)


#dev.off()
# make a cluster heatmap
heatmap.2(heat_mat, 
          scale = "none",
          #margins = c(5,5),
          Colv = F, Rowv = T,
          trace = "none",
          col = white.black,
          breaks = breaks)


table(pop_labels)

CD4 <- c(1,2,3,7,10,12,13,16,19)
CD8 <- c(4,5,6,8,9,11)
TCRgd <- c(14)
B.cells <- c(18)
PC <- c()
NK <- c(15,20)
Treg <- c() #Treg
Myeloid <- c(17)
DC <- c()

# CD4 <- c(1,2,6,8,10,11)
# CD8 <- c(9,3,4,5,7)
# TCRgd <- c(12)
# B.cells <- c(15)
# PC <- c()
# NK <- c(13)
# Treg <- c() #Treg
# Myeloid <- c(14)
# DC <- c()

# give them the right metacluster number
th_labels <- rep(0, length(pop_labels))
th_labels[pop_labels %in% CD4] <- 1
th_labels[pop_labels %in% CD8] <- 2
th_labels[pop_labels %in% TCRgd] <- 3
th_labels[pop_labels %in% B.cells] <- 4
th_labels[pop_labels %in% NK] <- 5
#th_labels[pop_labels %in% Treg] <- 6 #Treg
th_labels[pop_labels %in% Myeloid] <- 7
#th_labels[pop_labels %in% DC] <- 8

data[,"manual_labels"] <- th_labels

head(data)
data_ics <- data

data_all <- rbind(data_ics, data_surf)


####################
#####data_all ######
####################

#saveRDS(data_all, '052020_data_all_umap.RDS')
data_all <- readRDS('052020_data_all_umap.RDS')

dim(data_all)
data_all$cell_id <- 1:21942263

###Subsample per condition for umap
data_umap <- data_all[,c(make.names(intersect(panel_surf$Antigen[panel_surf$uneven =='n'], panel_ics$Antigen[panel_ics$uneven == 'y'])),'pz') ]
# data_umap[,c("CD4","CD8")] <- data_umap[,c("CD4","CD8")]/1.2
# data_umap[,c("CD7","HLADR")] <- data_umap[,c("CD7","HLADR")]*1.2
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

data_rtsne <- cell_tsne[cell_tsne[,'CD3'] >0  , colnames(data_umap)[-18]]
custom.settings = umap.defaults
custom.settings$verbose <- T
custom.settings$n_components <- 2
custom.settings$n_neighbors <- 25
custom.settings$random_state <- 123
custom.settings$n_epochs <- 500
custom.settings$min_dist <- 0.3
data_umap <-  umap(data_rtsne, custom.settings)


# prepare the tSNE data
tsne <- as.data.frame(data_umap$layout)
colnames(tsne) <- c("tSNE1", "tSNE2")
#saveRDS(tsne, 'tsne_multipanel.RDS')



# plot tSNE black
t1 <- ggplot(tsne, aes(x = tSNE1, y = tSNE2)) +
  geom_point(size = 0.5) +
  coord_fixed(ratio = 1) +
  theme_tsne + 
  theme_facs
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


# choose filter
data_cnn_ics <- droplevels(subset(data_plot, filter_0_continuous > 0.5 & tSNE2 <= 1))
data_cnn_surf <- droplevels(subset(data_plot, filter_2_continuous > 0.5))

# plot density tSNE
data_plot$pz <- as.factor(data_plot$pz)
str(data_plot)
t1 <- ggplot(data_plot, aes(x = tSNE1, y = tSNE2 )) +
  geom_density2d(color = "grey", size = 0.2) +
  stat_density_2d( bins=25)+
  #facet_wrap("diagnosis_spec") +
  xlim(min(data_plot$tSNE1)-1.5,max(data_plot$tSNE1)+1.5)+
  ylim(min(data_plot$tSNE2)-1.5,max(data_plot$tSNE2)+1.5)+
  #geom_point(size = 0.5) +
  coord_fixed(ratio = 1) +
  #facet_grid(vars(treatment2), vars(patient)) +
  #facet_wrap('patient') +
  #facet_wrap('timepoint2') +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  theme_tsne + 
  theme_facs
t1


t1 +  geom_point(data = data_cnn_ics, aes(x = tSNE1, y = tSNE2), 
                 size = 0.1, color = db1[4]) +
  geom_point(data = data_cnn_surf, aes(x = tSNE1, y = tSNE2),  
             size = 0.5, color = db1[1])


# plot point tSNE

t1 <- ggplot(data_plot, aes(x = tSNE1, y = tSNE2 )) +
  geom_point(color = "grey", size = 0.05) +
  #stat_density_2d( bins=25)+
  #facet_wrap("diagnosis_spec") +
  xlim(min(data_plot$tSNE1)-1.5,max(data_plot$tSNE1)+1.5)+
  ylim(min(data_plot$tSNE2)-1.5,max(data_plot$tSNE2)+1.5)+
  #geom_point(size = 0.5) +
  coord_fixed(ratio = 1) +
  #facet_grid(vars(treatment2), vars(patient)) +
  #facet_wrap('patient') +
  #facet_wrap('timepoint2') +
  scale_fill_manual(values = db1) +
  scale_color_manual(values = db1) +
  #theme_tsne + 
  theme_facs
t1

t1 +  geom_density2d(data = data_cnn_ics, aes(x = tSNE1, y = tSNE2), 
                     color = db1[4]) +
  geom_density2d(data = data_cnn_surf, aes(x = tSNE1, y = tSNE2), 
                 color = db1[1]) 

# save the plots as png
ggsave(filename = "all_tsne_selected_both.png", scale = 0.5)



# plot tSNEs with expression overlayed
plot_cols <- union(make.names(as.vector(panel_surf$Antigen[ panel_surf$Antigen %in% colnames(data_all) ])),make.names(as.vector(panel_ics$Antigen[ panel_ics$Antigen %in% colnames(data_all) ])))
data_plot[,c("CD4","CD8")] <- data_plot[,c("CD4","CD8")]/1.2
data_melt <- melt(data_plot[,c(plot_cols, "tSNE1", "tSNE2")], id.vars = c("tSNE1", "tSNE2"))
data_melt$variable <- factor(data_melt$variable, levels = levels(data_melt$variable)[order(levels(data_melt$variable))])
t2 <- ggplot(na.omit(data_melt), aes(x = tSNE1, y = tSNE2, color = value)) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = jet.colors(100), limits = c(-0.2,1.2)) +
  facet_wrap(~ variable, ncol = 10, scales = "fixed") +
  theme_tsne
t2


t3 <- ggplot(na.omit(subset(data_melt, variable == "CD14" )), aes(x = tSNE1, y = tSNE2, color = value)) +
  geom_point(size = 0.05) +
  coord_fixed(ratio = 1) +
  scale_colour_gradientn(colours = jet.colors(100), limits = c(-0.2,1.2)) +
  #facet_wrap(~ variable, ncol = 5, scales = "free") +
  theme_tsne
t3


# save the plots as png
ggsave(filename = "tsne_expression_final.png", t2, scale = 2, 
       width = 10*2, height = 5*2.3, units = c("in"))


#plot meta_clusters

# read in metadata
md_cn[,c(1:4,7,8,10)] <- lapply(md_cn[,c(1:4,7,8,10)], factor)
data_dff <- data_plot
data_dff <- merge(data_dff, md_cn, by = 'pz')
data_dff[data_dff[,'tSNE1']> 0 &data_dff[,'tSNE2']> 2&data_dff[,'manual_labels'] == 1,  'manual_labels' ]  <- 7
data_dff[data_dff[,'tSNE2']> 2&data_dff[,'manual_labels'] == 1,  'manual_labels' ]  <- 4
#data_dff[data_dff[,'tSNE1']> 0 &data_dff[,'tSNE2']> 2&data_dff[,'manual_labels'] == 2,  'manual_labels' ]  <- 7

# prepare data for Rtsne
data_dff[data_dff[,"manual_labels"] == 1,"manual_labels"] <- "CD4" 
data_dff[data_dff[,"manual_labels"] == 2,"manual_labels"] <- "CD8" 
data_dff[data_dff[,"manual_labels"] == 3,"manual_labels"] <- "gdT.cells" 
data_dff[data_dff[,"manual_labels"] == 4,"manual_labels"] <- "B.cell" 
data_dff[data_dff[,"manual_labels"] == 5,"manual_labels"] <- "NK" 
data_dff[data_dff[,"manual_labels"] == 6,"manual_labels"] <- "Treg" 
data_dff[data_dff[,"manual_labels"] == 7,"manual_labels"] <- "Myeloid" 
data_dff$manual_labels <- factor(data_dff$manual_labels, levels = c("CD4","CD8", "gdT.cells", "B.cell", "NK", 'Treg',  "Myeloid")) 
#data_dff$diagnosis2 <- factor(data_dff$diagnosis2, levels = c('HC','CIS','MS','NINDC','INDC')) 
#data_dff$pop_labels <- factor(data_dff$pop_labels, levels = c(1:16)) 

data_dff$pz <- as.numeric(data_dff$pz)
#dev.off()
names(db_lin)[8] <- 'Treg'
t3 <- ggplot(subset(data_dff,pz >=1) , aes(x = tSNE1, y = tSNE2, color =  manual_labels)) +
  # geom_density2d(color = "grey", size = 0.5) +
  # stat_density_2d( bins=25)+
  geom_point(size = 0.01) +
  coord_fixed(ratio = 1) +
  #scale_x_discrete(limits = c("CD4","CD8", "TCRgd", "B cell", "NK", "Myeloid"))+
  scale_color_manual(values = c(db_lin)) +
  theme_facs+
  #facet_grid(rows =  vars(responder), cols = vars(tpc)) +
  #facet_wrap('tpc', nrow = 1) +
  theme_tsne +
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
ggsave(filename = "dmf_manual_lab_tpc.png", plot = t3, 
       scale = 1.5, width = 3*1.65, height = 2*3.5, units = c("in"))


















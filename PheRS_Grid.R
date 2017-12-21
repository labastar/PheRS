setwd('~/Dropbox (Personal)/Lisa HPO - PRS/Science Submission/Resubmission 4/PheRS phenotype grid/')
library(pheatmap)

## disease ID for grid
dID_for_grid="OMIM_219700"

## IDs to put in grid
IDs_for_grid=c(11, 234, 813, 213, 455, 891, 44, 233, 377, 610, 987)


PheRS_map=read.table('disease_to_phecodes.txt', sep="\t",header=T,stringsAsFactors = F)
phecode_strings=read.table('phecode_strings.txt', sep="\t",header=T,stringsAsFactors = F)
head(PheRS_map)

icd9s=read.table('icd9s_sample_data.txt',sep="\t",header=T)

### create phecode phenotypes
icd9_to_phecodes=read.table("icd9_to_phecode.txt", sep="\t",header=T,stringsAsFactors = F,colClasses = c('character','character'))
phecodes=merge(icd9s,icd9_to_phecodes,by="icd9")
phecodes=unique(phecodes[c("ID", "phecode")])
phecodes$value = 1
phecodes <- spread(phecodes,  ID,value,fill=0)

## Get phecodes for dID
phecode_list=as.character(PheRS_map[PheRS_map$MIM==dID_for_grid,]$phecodes)
phecode_list=strsplit(phecode_list,",")[[1]]

## Subset phecode matrix by IDs in IDs_for_grid
grid_df=phecodes[, c('phecode',IDs_for_grid)]



## Subset phecode matrix by defining phecodes for dID
grid_df=grid_df[grid_df$phecode %in% phecode_list,]

## Set rownames to phecode strings 
grid_df=merge(grid_df,phecode_strings,by="phecode")
row.names(grid_df)=grid_df$phecode_string
grid_df$phecode=NULL
grid_df$phecode_string=NULL

## Make pheatmap
pheatmap(grid_df,color=c("white","black"),treeheight_row=0,treeheight_col=0,legend=F)

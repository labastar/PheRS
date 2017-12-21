library(tidyr)

## Read in map files
PheRS_map=read.table('disease_to_phecodes.txt', sep="\t",header=T,stringsAsFactors = F)
icd9_to_phecodes=read.table("icd9_to_phecode.txt", sep="\t",header=T,stringsAsFactors = F,colClasses = c('character','character'))

## Read in sample icd9 data for individuals
icd9s=read.table('icd9s_sample_data.txt',sep="\t",header=T)


### map individual icd9s to phecodes
phecodes=merge(icd9s,icd9_to_phecodes,by="icd9")
phecodes=unique(phecodes[c("ID", "phecode")])
phecodes$value = 1
phecodes <- spread(phecodes,  phecode,value,fill=0)


## Calculate weights from icd9 data
pop_count=nrow(phecodes)
weights=data.frame(colSums(phecodes[2:length(phecodes)]))
weights$phecode= rownames(weights)
names(weights)[1]="case_count"
weights$prev = weights$case_count/pop_count
weights$w=log10(pop_count/weights$case_count)
## If you want to use weights generated from the adult VUMC cohort used in the discovery cohort, uncomment the next 2 lines
#weights=read.table('weights_VUMC_discovery.txt', sep="\t", header=T,stringsAsFactors = F)
#rownames(weights)=weights$phecode


## add weights to phecodes
phes=names(phecodes)[-1]
for(i in 1:length(phes)){
  phe=phes[i]
  phecodes[[phe]]=phecodes[[phe]]*weights[weights$phecode==phe,]$w
}



## Create PheRS score
MIMs=unique(PheRS_map$MIM)
PheRS=data.frame(phecodes$ID)
names(PheRS)[1]="ID"
for(i in 1:length(MIMs)){
  MIM=MIMs[i]
  phecode_list=as.character(PheRS_map[PheRS_map$MIM==MIM,]$phecodes)
  phecode_list=strsplit(phecode_list,",")[[1]]
  phecodes[phecode_list]
  PheRS[MIM]=rowSums(phecodes[phecode_list],na.rm=T)
  
}

## Add FID and IID for Plink
PheRS$FID = PheRS$ID
names(PheRS)[1]="IID"

## Reorder columns
PheRS=PheRS[c("FID","IID",MIMs)]

## Write out a file formatted for Plink
write.table(PheRS,file="PheRS_file_sample.txt",quote=F,sep="\t",row.names=F,col.names=T)


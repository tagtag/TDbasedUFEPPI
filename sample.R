require(Matrix)
require(irlba)

#For BIOGRID data set
#human PPI
#PPI for BIOGRID is supposed to be downloaded as suggested in the paper
PPI_human <- read.delim("~/RESEARCH/PPI_tensor_2/BIOGRID/BIOGRID-ORGANISM-Homo_sapiens-4.4.236.tab3.txt.gz", header=FALSE, comment.char="#")
gene_all <- union(PPI_human[,24],PPI_human[,27]) # get union of various gene IDs and uniprotkb
TABLE2 <- sparseMatrix(i=match(PPI_human[,24],gene_all),j=match(PPI_human[,27],gene_all),dims=rep(length(gene_all),2)) #stored as sparse matrix format
TABLE2 <- TABLE2+t(TABLE2)  #PPI is symmetrized
TABLE2[TABLE2==2] <- TABLE2[TABLE2==2]/2  #duplicated PPI is halven 
rownames(TABLE2) <- gene_all  #row names are assigned 
colnames(TABLE2) <- gene_all  #col names are assigned 

#mouse PPI
#PPI for BIOGRID is supposed to be downloaded as suggested in the paper
PPI_mouse <- read.delim("~/RESEARCH/PPI_tensor_2/BIOGRID/BIOGRID-ORGANISM-Mus_musculus-4.4.236.tab3.txt.gz", header=FALSE, comment.char="#")
gene1_all <- union(PPI_mouse[,24],PPI_mouse[,27])  # get union of various gene IDs and uniprotkb
TABLE2p <- sparseMatrix(i=match(PPI_mouse[,24],gene1_all),j=match(PPI_mouse[,27],gene1_all),dims=rep(length(gene1_all),2)) #stored as sparse matrix format
TABLE2p <- TABLE2p+t(TABLE2p)  #PPI is symmetrized
TABLE2p[TABLE2p==2] <- TABLE2p[TABLE2p==2]/2  #duplicated PPI is halven 
rownames(TABLE2p) <- gene1_all  #row names are assigned 
colnames(TABLE2p) <- gene1_all #col names are assigned 


#Unification of human and mouse PPI with orthologs
#Ortholog infoirmation is supposed to be downloaded from uniprot as suggested by the paper
ortholog <- read.delim("~/RESEARCH/PPI_tensor_2/Ortholog/uniprotkb_organism_id_10090_OR_organism_2024_05_17.tsv.gz")
uniprot <- unlist(lapply(strsplit(rownames(TABLE2),"|",fixed=T),"[",1)) #get uniprotkb for human
uniprotp <- unlist(lapply(strsplit(rownames(TABLE2p),"|",fixed=T),"[",1)) #get uniprotkb for mouse

#Mouse and Human PPI are aligned as suggested in the paper
uniprot_common <- intersect(uniprot,uniprotp)
length(uniprot_common)
#[1] 9688 vs mus
index <- match(uniprot_common,uniprot)
index <- c(index,setdiff(seq_len(dim(TABLE2)[1]),index))
index1 <- match(uniprot_common,uniprotp)
index1 <- c(index1,setdiff(seq_len(dim(TABLE2p)[1]),index1))

TABLE2 <- TABLE2[index,index]
uniprot <- uniprot[index]
TABLE2p <- TABLE2p[index1,index1]
uniprotp <- uniprotp[index1]

index12 <- length(uniprot_common) +match(ortholog[match(gsub("HUMAN","MOUSE",ortholog[match(uniprot[-seq_along(uniprot_common)],ortholog[,1]),3]),ortholog[,3]),1],uniprotp[-seq_along(uniprot_common)])
#index12 <- length(uniprot_common) +match(ortholog[match(gsub("HUMAN","RAT",ortholog[match(uniprot[-seq_along(uniprot_common)],ortholog[,1]),3]),ortholog[,3]),1],uniprotp[-seq_along(uniprot_common)])
index11 <- cbind(length(uniprot_common) +seq_along(index12),index12)
index11<-index11[!is.na(index11[,2]),]
index11 <- index11[!duplicated(index11[,2]),]


num <- dim(TABLE2)[1]+dim(TABLE2p)[1]-(length(uniprot_common)+sum(!is.na(index11[,1])))
TABLE3 <-Matrix(0,num,num);Tn <- rep(NA,num)
TABLE3p<-Matrix(0,num,num);Tnp <- rep(NA,num)
inum <- c(seq_along(uniprot_common),index11[,1])
inum1 <- c(seq_along(uniprot_common),index11[,2])
TABLE3[inum,inum] <- TABLE2[inum,inum]
TABLE3p[inum1,inum1] <- TABLE2p[inum1,inum1]
Tn[inum] <- rownames( TABLE2[inum,inum])
Tnp[inum1] <- rownames( TABLE2p[inum1,inum1])
TABLE3[(length(inum)+1):length(uniprot),seq_along(uniprot)] <- TABLE2[-inum,]
TABLE3[seq_along(uniprot),(length(inum)+1):length(uniprot)] <- TABLE2[,-inum]
Tn[(length(inum)+1):length(uniprot)] <- rownames(TABLE2[-inum,])
TABLE3p[(length(uniprot)+1):dim(TABLE3p)[1],
        c(seq_along(inum),(length(uniprot)+1):dim(TABLE3p)[1])] <- TABLE2p[-inum1,]
TABLE3p[ c(seq_along(inum),(length(uniprot)+1):dim(TABLE3p)[1]),
         (length(uniprot)+1):dim(TABLE3p)[1]] <- TABLE2p[,-inum1]
Tnp[(length(uniprot)+1):dim(TABLE3p)[1]] <- rownames(TABLE2p[-inum1,])

#compute HOSVD by applying SVD to unfolded tensor
SVD <- irlba(cbind(TABLE3,TABLE3p),10)
P <- pchisq(scale(SVD$u[,2])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01) #selecting genes by u2i 
FALSE  TRUE 
27611   195 

P <- pchisq(scale(SVD$u[,3])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01) #selecting genes by u3i

FALSE  TRUE 
27331   475 


#For DIP data set
#human PPI
x <- read.csv("Hsapi20170205.txt",sep="\t") #PPI for DIP is supposed to be downloaded as suggested in the paper
TABLE <- table(x[,1],x[,2])
TABLE2 <- data.matrix(TABLE)
gene_all <- union(x[,1],x[,2]) # get union of various gene IDs and uniprotkb
TABLE2 <- sparseMatrix(i=match(x[,1],gene_all),j=match(x[,2],gene_all),dims=rep(length(gene_all),2)) #stored as sparse matrix format
TABLE2 <- TABLE2+t(TABLE2) #PPI is symmetrized
TABLE2[TABLE2==2] <- TABLE2[TABLE2==2]/2 #duplicated PPI is halven 
rownames(TABLE2) <- gene_all #row names are assigned 
colnames(TABLE2) <- gene_all #col names are assigned

#mouse PPI
x1 <- read.csv("Mmusc20170205.txt",sep="\t") #PPI for DIP is supposed to be downloaded as suggested in the paper
#---- symetric ---
gene1_all <- union(x1[,1],x1[,2]) # get union of various gene IDs and uniprotkb
TABLE2p <- sparseMatrix(i=match(x1[,1],gene1_all),j=match(x1[,2],gene1_all),dims=rep(length(gene1_all),2)) #stored as sparse matrix format
TABLE2p <- TABLE2p+t(TABLE2p) #PPI is symmetrized
TABLE2p[TABLE2p==2] <- TABLE2p[TABLE2p==2]/2 #duplicated PPI is halven 
rownames(TABLE2p) <- gene1_all #row names are assigned 
colnames(TABLE2p) <- gene1_all #col names are assigned

#Unification of human and mouse PPI with orthologs
#Ortholog information is supposed to be downloaded from uniprot as suggested by the paper 
ortholog <- read.delim("~/RESEARCH/PPI_tensor_2/Ortholog/uniprotkb_organism_id_10090_OR_organism_2024_05_17.tsv.gz")
uniprot <- unlist(lapply(strsplit(rownames(TABLE2),"|",fixed=T),function(x){x[length(x)]}))
uniprot <-gsub("uniprotkb:","",uniprot) #get uniprotkb for human
uniprotp <- unlist(lapply(strsplit(rownames(TABLE2p),"|",fixed=T),function(x){x[length(x)]}))
uniprotp <-gsub("uniprotkb:","",uniprotp)  #get uniprotkb for mouse

#Mouse and Human PPI are aligned as suggested in the paper
uniprot_common <- intersect(uniprot,uniprotp)
length(uniprot_common)
#[1] 1097
index <- match(uniprot_common,uniprot)
index <- c(index,setdiff(seq_len(dim(TABLE2)[1]),index))
index1 <- match(uniprot_common,uniprotp)
index1 <- c(index1,setdiff(seq_len(dim(TABLE2p)[1]),index1))

TABLE2 <- TABLE2[index,index]
uniprot <- uniprot[index]
TABLE2p <- TABLE2p[index1,index1]
uniprotp <- uniprotp[index1]

index12 <- length(uniprot_common) +match(ortholog[match(gsub("HUMAN","MOUSE",ortholog[match(uniprot[-seq_along(uniprot_common)],ortholog[,1]),3]),ortholog[,3]),1],uniprotp[-seq_along(uniprot_common)])
index11 <- cbind(length(uniprot_common) +seq_along(index12),index12)
index11<-index11[!is.na(index11[,2]),]
index11 <- index11[-369,] #remove duplicated 1190

num <- dim(TABLE2)[1]+dim(TABLE2p)[1]-(length(uniprot_common)+sum(!is.na(index11[,1])))
TABLE3 <-Matrix(0,num,num);Tn <- rep(NA,num)
TABLE3p<-Matrix(0,num,num);Tnp <- rep(NA,num)
inum <- c(seq_along(uniprot_common),index11[,1])
inum1 <- c(seq_along(uniprot_common),index11[,2])
TABLE3[inum,inum] <- TABLE2[inum,inum]
TABLE3p[inum1,inum1] <- TABLE2p[inum1,inum1]
Tn[inum] <- rownames( TABLE2[inum,inum])
Tnp[inum1] <- rownames( TABLE2p[inum1,inum1])
TABLE3[(length(inum)+1):length(uniprot),seq_along(uniprot)] <- TABLE2[-inum,]
TABLE3[seq_along(uniprot),(length(inum)+1):length(uniprot)] <- TABLE2[,-inum]
Tn[(length(inum)+1):length(uniprot)] <- rownames(TABLE2[-inum,])
TABLE3p[(length(uniprot)+1):dim(TABLE3p)[1],
        c(seq_along(inum),(length(uniprot)+1):dim(TABLE3p)[1])] <- TABLE2p[-inum1,]
TABLE3p[ c(seq_along(inum),(length(uniprot)+1):dim(TABLE3p)[1]),
         (length(uniprot)+1):dim(TABLE3p)[1]] <- TABLE2p[,-inum1]
Tnp[(length(uniprot)+1):dim(TABLE3p)[1]] <- rownames(TABLE2p[-inum1,])


#compute HOSVD by applying SVD to unfolded tensor
SVD <- irlba(cbind(TABLE3,TABLE3p),10)
P <- pchisq(scale(SVD$u[,2])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01) #selecting genes by u2i

FALSE  TRUE 
5518   196 

P <- pchisq(scale(SVD$u[,3])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)  #selecting genes by u3i

FALSE  TRUE 
5655    59 

P <- pchisq(scale(SVD$u[,6])^2,1,lower.tail=F)
table(p.adjust(P,"BH")<0.01)  #selecting genes by u6i

FALSE  TRUE 
5665    49 



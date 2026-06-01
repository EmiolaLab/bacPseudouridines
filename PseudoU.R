library(data.table)
library(tidyr)
library(dplyr)
library(getopt)

spec = matrix(c('ufile','u',1,"character",'tfile','t',1,"character",'gfile','g',1,"character",'pval','p',2,"double",'ofile','o',1,"character"),byrow=TRUE,ncol=4)
opt=getopt(spec)
untreated_file <- opt$ufile
treated_file <- opt$tfile
gff_file <- opt$gfile
output_file <- opt$ofile
pval <- opt$pval


## Read the gff file which has strand information

annot <- fread(gff_file,header = F,sep = '\t') 
annot$V9 <- gsub(";.*","", annot$V9)
annot$V9 <- gsub("ID=","", annot$V9)
annot <- annot[is.na(annot$V3) | annot$V3 != "repeat_region", ] # Deleting repeat regions such as crisprs
# Reannote strand ID
annot$V7 <- gsub("\\-","minus", annot$V7)
annot$V7 <- gsub("\\+","plus", annot$V7)
colnames(annot) <- c("chr","Prodigal","Type","Start","End","Blank","Strand","Dash","Gene")

# Analyze BS-treated sample
Treated <- fread(treated_file,header = F,sep = '\t')
colnames(Treated) <- c("position","chr","ref","depth","count","base","positive_strand","negative_strand","percent_bias","vaf","mutation")
    
t <- as.data.frame(Treated %>% group_by(position,chr) %>% summarise(countsum = sum(count)) ) #Aggregate and sum the counts per position. 
Treated <- merge(Treated,t, by=c("position", "chr"))

Treated <- Treated[which(Treated$countsum >= 5),] # Deletion count >= 5
Treated$count <- Treated$countsum
Treated <- unite(Treated, "Genome_pos", c("position","chr"), remove = FALSE) 
    
#####
t <- as.data.frame(Treated %>% group_by(position,chr) %>% summarise(Sum = sum(vaf)) ) #Aggregate and sum the deletion rations per position.
Treated <- merge(Treated,t, by=c("position", "chr"))
Treated <- unite(Treated, "Genome_pos", c("position","chr"), remove = FALSE) # Adding a new column to group Genome and position

# There can be multiple deletion pattern per position. e.g -T, -TTC, -TAA etc. We keep the deletion pattern with the highest del ratio
Treated <- setDT(Treated)[, .SD[which.max(vaf)], by=Genome_pos] 
Treated$vaf <- Treated$Sum
Treated <- Treated[which(Treated$vaf >= 0.02),] # Keeping on positions with >= 2% deletion ratios
Treated <- Treated[,c(1:12)]
Treated <- as.data.frame(Treated)

# Integrating genome annotation information from gff file    
annot_tmp <- annot[annot$chr %in% unique(Treated$chr), ] # subsetting only genomes in the treated samples 

# We retrieve the gene info for each position    
    for(x in Treated$Genome_pos){
      i <- sub(paste0("_", ".*"), "", x)
      i <- as.integer(i)
      z <- sub("^[^_]*_", "", x)
      var <- annot_tmp[with(annot_tmp, Start <= i & End >= i & chr == z),]
      var$Genome_pos <- x
      assign(paste("r_out", x, sep = "."),as.data.frame(var))
    }
    TMP <- as.data.frame(rbindlist(mget(ls(pattern = "^r_out."))))
    rm(list=ls(pattern="^r_out."))
    
Treated <- merge(Treated, TMP, by="Genome_pos")
Treated <- Treated[!duplicated(Treated),]
    
# Get strand info
t1 <- Treated[which(Treated$ref == "T" & Treated$Strand == "plus"),]  
t2 <- Treated[which(Treated$ref == "A" & Treated$Strand == "minus"),] #for genes in reverse strand, putative pseudouridine sites are 'A'
Treated <- rbind(t1,t2)
    
#############
# Analyze Untreated sample
Untreated <- fread(untreated_file, sep = '\t', header = T) # Gene cluster info
colnames(Untreated) <- c("position","chr","ref","depth","count","base","positive_strand","negative_strand","percent_bias","vaf","mutation")
Untreated <- unite(Untreated, "Genome_pos", c("position","chr"), remove = FALSE)
Untreated <- Untreated[Untreated$Genome_pos %in% Treated$Genome_pos, ] # Get positions that are only present in Treated samples
    
# Retrieve data from sense & antisense strand
u1 <- Untreated[which(Untreated$ref == "T"),]
u1 <- u1[grep("^\\(-T|^\\(T", u1$base),] #fetching del and non-del sities. i.e. all T sites
    
u2 <- Untreated[which(Untreated$ref == "A"),]
u2 <- u2[grep("^\\(-A|^\\(A", u2$base),]
    
Untreated <- rbind(u1,u2)

# processing deletion sites in untreated samples    
U_del <- Untreated[grep("del", Untreated$mutation),] 
    
t <- as.data.frame(U_del %>% group_by(position,chr) %>% summarise(countsum = sum(count)) ) #Aggregate and sum the counts per position. 
U_del <- merge(U_del,t, by=c("position", "chr"))
U_del$count <- U_del$countsum
    
t <- as.data.frame(U_del %>% group_by(position,chr) %>% summarise(Sum = sum(vaf)) ) #Aggregate and sum the deletion rations per position.
U_del <- merge(U_del,t, by=c("position", "chr"))
U_del <- setDT(U_del)[, .SD[which.max(vaf)], by=Genome_pos] # keep max in each group
U_del$vaf <- U_del$Sum
U_del <- U_del[,c(1:12)]
    
U_nodel <- Untreated[!Untreated$Genome_pos %in% U_del$Genome_pos, ] #non-del sites

Untreated <- rbind(U_del, U_nodel)
    
DF <- merge(Treated,Untreated,by="Genome_pos")
DF$sample <- treated_file
DF <- DF[!duplicated(DF$Genome_pos), ]
DF$vaf.y[DF$mutation.y != "del"] <- 0 # replace vaf in non-mutation sites with zero

#Deletion ratio in BS-treated samples should be at least 2-fold higher than in untreated samples  
DF$ratio <- DF$vaf.x/DF$vaf.y
DF <- do.call(data.frame,lapply(DF, function(x) replace(x, is.infinite(x),1000))) #replace Inf with arbitrary 1000. This happens when del ratio in untreated sample is zero
DF <- DF[which(DF$ratio >= 2),] 
DF$count.y[DF$count.y == 0] <- DF$depth.y[DF$count.y == 0] 
    
#Fisher test comparing deletion counts and total read depth between BS-treated and untreated samples
    for(x in unique(DF$Genome_pos)){
      Fis <- DF[which(DF$Genome_pos == x),]
      
      if(Fis$mutation.y == "no-mutation")
      {
        Fis$T_nodel <- Fis$depth.x - Fis$count.x
        Fis$U_del <- 0
        Fis <- Fis[,c(6,25,35,36)]
        Fis <- as.data.frame(mapply(c, Fis[,c(1,3)], Fis[,c(4,2)]))
      } else {
        
        Fis$T_nodel <- Fis$depth.x - Fis$count.x
        Fis$U_nodel <- Fis$depth.y - Fis$count.y
        Fis <- Fis[,c(6,26,35,36)]
        Fis <- as.data.frame(mapply(c, Fis[,c(1,3)], Fis[,c(2,4)]))
      }
      colnames(Fis) <- c("deletion","no deletion")
      rownames(Fis) <- c("Treated","Untreated")
      pval <- fisher.test(Fis)$p.value
      assign(paste("r_out", x, sep = "."),as.data.frame(merge(x,pval)))
    }
fisher_result <- as.data.frame(rbindlist(mget(ls(pattern = "^r_out."))))
rm(list=ls(pattern="^r_out."))
colnames(fisher_result) <- c("Genome_pos", "pvalue")
DF <- merge(DF, fisher_result, by="Genome_pos")
    
    # filter by pvalue
DF <- DF[which(DF$pvalue < pval),]
DF$delta_del_ratio <- DF$vaf.x - DF$vaf.y 

cols <- c("sample", "Genome_pos", "position.x", "chr.x", "ref.x","depth.x","count.x", "vaf.x","mutation.x",
          "Type","Start","End","Strand","Gene","depth.y","count.y","vaf.y", "mutation.y", "pvalue","delta_del_ratio")

PsedoUridine_sites <- DF[, cols]
    
colnames(PsedoUridine_sites) <- c("Sample","Position_Genome", "Position", "Genome", "Reference","Depth (BS-treated)",
                                  "Count (BS-treated)","Deletion ratio (BS-treated)", "Mutation (BS-treated)", "Type","Start","End","Strand","Gene",
                                  "Depth (untreated)","Count (untreated)","Deletion ratio (untreated)", "Mutation (Untreated)","Pvalue","Delta deletion ratio")
 

write.table(PsedoUridine_sites, file = output_file, sep = '\t', col.names = T, row.names = F, quote = F)

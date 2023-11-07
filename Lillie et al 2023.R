#### Characterising polymoprhic ERVs in KOALA ####
### 2023/11/07 Mette Lillie mette.lillie@ebc.uu.se

#### SETUP | dir - packaages ####
setwd("/Users/metli609/Documents/PostDoc_UUJern/Koala 2/github/") 
#require(dplyr)
require(data.table)
#require(plyr)
#require(StructuralVariantAnnotation)
#require(VariantAnnotation)
require(GenomicRanges) 
#require(rtracklayer)
require(regioneR)
#require(stringr)
#require(RColorBrewer)
#require(intansv)
require(geodist)
#require(pheatmap)
#

#
#### binnedSum function: ####
binnedSum <- function(bins, numvar, mcolname)
{
  stopifnot(is(bins, "GRanges"))
  stopifnot(is(numvar, "RleList"))
  stopifnot(identical(seqlevels(bins), names(numvar)))
  bins_per_chrom <- split(ranges(bins), seqnames(bins))
  sums_list <- lapply(names(numvar),
                      function(seqname) {
                        views <- Views(numvar[[seqname]],
                                       bins_per_chrom[[seqname]])
                        viewSums(views)
                      })
  new_mcol <- unsplit(sums_list, as.factor(seqnames(bins)))
  mcols(bins)[[mcolname]] <- new_mcol
  bins
}

#### KOALA REFERENCE | fai ####
# read in koala reference fai file:
ref.fai <- fread("./files/GCF_002099425.1_phaCin_unsw_v4.1_genomic.fna.fai")
ref.fai <- ref.fai[,1:2]
names(ref.fai) <- c("contig", "length")
ref.fai # 1907 contigs (1906 autosomal + mt)
## read in reference file for contigs accessions from genbank and refseq 
contig.ref <- fread("./files/koala_reference_contig_reference.csv")
contig.ref


ref.fai.MTS <- ref.fai
ref.fai.MTS$contig <- contig.ref$genbank.accn[match(ref.fai.MTS$contig, contig.ref$refseq.accn)]
ref.fai.MTS$contig[ref.fai.MTS$contig == "na"] <- "MT"

## read in reference file for ervs found by retrotector and their locations on genbank and refseq scaffs
erv.reference <- fread("./files/erv_reference.csv")
head(erv.reference)


#### KOALA REFERENCE | ERV info from retrotector ####
# all the retrotector hits with score over 300
erv.scaff.map.range <- toGRanges(erv.reference[,c("hitrefseq", "scaff.start", "scaff.end", "erv.id")])
erv.scaff.map.range # 991
length(unique(erv.scaff.map.range@seqnames)) # 279 scaffolds with retrotectored ervs



## read in file with details of the retrotector ervs used in the phylo and which clades they cluster into 
erv.clades <- fread("./files/erv_clades.csv")
erv.clades

# nb: only 129; this is only based upon the pCi phylo from Lillie et al 2022 


#### KOALA REFERENCE | ERVs curated ####
# curated list of candidate ervs in the reference 
# from retrotector and blat inputs
# blat inputs is a curated library taking blat hits to three clades LTRs and longer hits to provirus querys for three clades
# NB NW_018344052.1	7608537	7611871	3334	*	pCi954	468	blat	KoRV is a KoRV recomb.. 

reference.erv.candidates.threeclades.final.df <- fread("./files/reference.erv.candidates.threeclades.csv")
reference.erv.candidates.threeclades.final <- toGRanges(reference.erv.candidates.threeclades.final.df[,c("seqnames", "start","end" , "erv" ,     "score"  ,  "method", "clade")])
reference.erv.candidates.threeclades.final
reference.erv.candidates.threeclades.final$clade[reference.erv.candidates.threeclades.final$clade == "phaCinbeta-like"] <- "phaCinBetalike"
reference.erv.candidates.threeclades.final$clade[reference.erv.candidates.threeclades.final$clade == "phaCinbeta"] <- "phaCinBeta"
reference.erv.candidates.threeclades.final$clade[reference.erv.candidates.threeclades.final$clade == "KoRV"] <- "KoRVlike"

# GRanges object with 730 ranges and 4 metadata columns:
# seqnames            ranges strand |         erv     score          method          clade

reference.erv.candidates.threeclades.final.MTS <- reference.erv.candidates.threeclades.final

reference.erv.candidates.threeclades.final.MTS$NW <- reference.erv.candidates.threeclades.final.MTS@seqnames
reference.erv.candidates.threeclades.final.MTS$MTS <- contig.ref$genbank.accn[match(as.character(reference.erv.candidates.threeclades.final.MTS$NW), contig.ref$refseq.accn)]
reference.erv.candidates.threeclades.final.MTS
reference.erv.candidates.threeclades.final.MTS <- toGRanges(reference.erv.candidates.threeclades.final.MTS$MTS, 
                                                            reference.erv.candidates.threeclades.final.MTS@ranges, 
                                                            reference.erv.candidates.threeclades.final.MTS$erv, 
                                                            reference.erv.candidates.threeclades.final.MTS$score, 
                                                            reference.erv.candidates.threeclades.final.MTS$method, 
                                                            reference.erv.candidates.threeclades.final.MTS$clade)
names(mcols(reference.erv.candidates.threeclades.final.MTS)) <- c("width", "number", "erv", "score", "method", "clade")




ref.ervs.bed <- data.frame(chr=reference.erv.candidates.threeclades.final@seqnames, 
                           start=reference.erv.candidates.threeclades.final@ranges@start - 1000, 
                           end=reference.erv.candidates.threeclades.final@ranges@start + reference.erv.candidates.threeclades.final@ranges@width+999)
ref.ervs.bed # 1000 bp up and downstream of locus
#write.table(ref.ervs.bed, file="ref.ervs.bed", quote =F, row.names = F, col.names = F)

ref.ervs.bed$otherchr <- contig.ref$genbank.accn[match(ref.ervs.bed$chr, contig.ref$refseq.accn)]
ref.ervs.bed
#write.table(ref.ervs.bed[,c(4,2,3)], file="ref.ervs.bed", quote =F, row.names = F, col.names = F, sep="\t")


ref.ervs.bed.MST.granges <- toGRanges(ref.ervs.bed[,c(4,2,3)])


#### KOALA REFERENCE | REPEATMASKER  ####
## read in repeatmasker output for the koala reference genome
ref.repeatmasker <- fread("./files/GCF_002099425.1_phaCin_unsw_v4.1_rm_eds.out")
ref.repeatmasker <- ref.repeatmasker[3:nrow(ref.repeatmasker),]
head(ref.repeatmasker)
nrow(ref.repeatmasker) # 7814628
names(ref.repeatmasker) <- c("chr", "start", "end", "rep.class", "rep.fam")
ref.repeatmasker.range <- toGRanges(ref.repeatmasker)
head(ref.repeatmasker.range)
table(ref.repeatmasker.range$rep.fam)
ref.repeatmasker.range.simplerepeats <- ref.repeatmasker.range[ref.repeatmasker.range$rep.fam == "Simple_repeat",]

ref.repeatmasker.range.simplerepeats
ref.repeatmasker.range.simplerepeats$locus <- paste(ref.repeatmasker.range.simplerepeats)



#### RESEQ METADATA | downloaded AWS 2022/11/08 ####
# Sample data from latest release 2022/10/19 downloaded from AWS servers on 2022/11/08 
awgg_metadata <- fread("./files/Koala_Metadata-19-10-2022.csv")
names(awgg_metadata) <- gsub("\ ", "_", names(awgg_metadata))
# where is Featherdale_M_46879?
awgg_metadata[awgg_metadata$AWS_File_Name %like% "Featherdale_F_46879",]
# need to correct
awgg_metadata$AWS_File_Name[awgg_metadata$AWS_File_Name == "Featherdale_F_46879"] <- "Featherdale_M_46879"
# where is Broadwater_M_50382? Is it supposed to be Broad_M_50328??
awgg_metadata$AWS_File_Name[awgg_metadata$AWS_File_Name == "Broadwater_M_50328"] <- "Broadwater_M_50382"

## latitudes and longitudes for locations based on the names of folders
## csv with cols:
## AWS_Folder_Name latitude longitude
awgg_locslatslons <- fread("./files/awggfolders_latslons.csv")
awgg_locslatslons



#### AWGG locations | pairwise distances between regions #### 
awgg_locslatslons_prox <- awgg_locslatslons
names(awgg_locslatslons_prox) <- c("AWS_Folder_Name", "lat", "lon")
awgg_locslatslons_prox <- awgg_locslatslons_prox[awgg_locslatslons_prox$AWS_Folder_Name != "Captive",]
geodist(awgg_locslatslons_prox, measure = "haversine")/1000
# manually made awggnames.order.csv
awgg.the.order.of.things <- fread("./files/awggnames.order.csv", header = F)
names(awgg.the.order.of.things) <- "AWS_Folder_Name"
awgg.the.order.of.things$AWS_Folder_Name



#### RESEQ BAM COV | samtools depth #### 
## read in file of bam sequencing coverage
koala.cov <- fread("./files/awgg_koala_bams_reads.coverage", header = F)
koala.cov
names(koala.cov) <- c("sample", "cov")




#### DELETIONS | load LUMPY results ####
# lumpy deletions called on whole genome
# need to load in lumpy, then subset by overlap with reference ervs 


lumpy.dir <- "./lumpy_vcfs" 
lumpy.vcfs <- list.files(path = lumpy.dir, 
                         pattern = "*.vcf", full.names = T)
length(lumpy.vcfs) # nb only 429 since a PortMacq sample failed running lumpy -_- 

lumpy.summary <- data.frame(samp = gsub("\\..*", "", gsub(".lumpy.vcf", "", basename(lumpy.vcfs))), 
                            lumpy.dels = NA, 
                            lumpy.erv.dels=NA)

for (i in c(1:length(lumpy.vcfs))){
  name <- gsub(".lumpy.vcf", "", basename(lumpy.vcfs[i]))
  loadin <- VariantAnnotation::readVcf(lumpy.vcfs[i])
  loadin.intansv <- readLumpy(lumpy.vcfs[i], regSizeLowerCutoff=100, regSizeUpperCutoff=30000,method="Lumpy")
  loadin.temp <- toGRanges(loadin.intansv$del)
  lumpy.summary$lumpy.dels[i] <- length(loadin.temp)
  rm(loadin.intansv)
  grange.temp.erv <- subsetByOverlaps(loadin.temp, reference.erv.candidates.threeclades.final.MTS)
  lumpy.summary$lumpy.erv.dels[i] <- length(grange.temp.erv)
  if (length(grange.temp.erv) != 0){
    mcols(grange.temp.erv) <- NULL
    grange.temp.erv$method <- "lumpy"
    grange.temp.erv$sample <- name
  }
  assign(paste(name, ".lumpy.erv.ranges", sep=""), grange.temp.erv)
}
# this process takes a long while.. started approx 9am till before 1040.
# will have warnings when seqnames are not exactly the same 


plot(lumpy.summary$lumpy.dels, lumpy.summary$lumpy.erv.dels, pch=16, col="grey")

lumpy.summary[lumpy.summary$lumpy.dels < 2000,]
summary(lumpy.summary) # so between 11 and 43 erv dels per individual.

## lumpy to granges
obj.list = as.list(ls()[sapply(mget(ls(), .GlobalEnv), is.object)])
lumpy.granges.list <- obj.list[obj.list %like% ".lumpy.erv.ranges" ]
lumpy.granges.list

# compile
lumpy.outgrange <- get(lumpy.granges.list[[1]])
for (i in 2:length(lumpy.granges.list)){
  tempgrange <- get(lumpy.granges.list[[i]])
  lumpy.outgrange <- c(lumpy.outgrange, tempgrange)
}
lumpy.outgrange 

lumpy.erv.granges <- reduce(lumpy.outgrange) # 88 ranges 
lumpy.erv.loci.counts <- lumpy.erv.granges
mcols(lumpy.erv.loci.counts) <- NULL

# getting counts - evals <= 1 
for (i in 1:length(unique(lumpy.outgrange$sample))){
  g1 <- lumpy.outgrange[lumpy.outgrange$sample == unique(lumpy.outgrange$sample)[i]]
  lumpy.erv.loci.counts$hits <- countOverlaps(lumpy.erv.loci.counts, g1, type = "any")
  names(mcols(lumpy.erv.loci.counts))[names(mcols(lumpy.erv.loci.counts)) == "hits"] <- unique(lumpy.outgrange$sample)[i]
}
lumpy.erv.loci.counts



## get erv clade into lumpy data 
lumpy.erv.df <- lumpy.erv.loci.counts
mcols(lumpy.erv.df) <- NULL
lumpy.erv.df

lumpy.erv.df$erv<- NA
lumpy.erv.df$range <- NA
lumpy.erv.df$clade <- NA
lumpy.erv.df$method <- NA

for(i in 1:length(lumpy.erv.df)){
  temp <- subsetByOverlaps(reference.erv.candidates.threeclades.final.MTS, lumpy.erv.df[i,])
  if (length(temp) > 0) {
    lumpy.erv.df$erv[i] <- temp$erv
    lumpy.erv.df$size[i] <- temp$width
    lumpy.erv.df$clade[i] <- temp$clade
    lumpy.erv.df$method[i] <- temp$method
    lumpy.erv.df$range[i] <- paste(temp@ranges)
    lumpy.erv.df$intersect[i] <- intersect(reference.erv.candidates.threeclades.final.MTS, lumpy.erv.df[i,])@ranges@width
  }
}
lumpy.erv.df
lumpy.erv.df$size.match <- lumpy.erv.df$size - lumpy.erv.df@ranges@width
lumpy.erv.df$intersect.match <- lumpy.erv.df$intersect - lumpy.erv.df$size
lumpy.erv.df$intersect.match.pct <- lumpy.erv.df$intersect / lumpy.erv.df$size

lumpy.erv.df

# manual visualisation and curation of these results 

# retain those deletions that overlap at least 80% with ERVs in the reference genome 
lumpy.erv.df$locus <- paste(paste(lumpy.erv.df), lumpy.erv.df$clade, sep="_")
lumpy.erv.df <- lumpy.erv.df[lumpy.erv.df$intersect.match.pct > 0.8,]
length(lumpy.erv.df) # 76
lumpy.erv.loci.counts <- lumpy.erv.loci.counts[ranges(lumpy.erv.loci.counts) %in% ranges(lumpy.erv.df)]



#### DELETIONS | load DELLY results ####
# delly has already been subsetted by overlap with bed file: ref.ervs.bed (730 reference ervs with 1000 bp up and down stream regions)

delly.intansv <- readDelly(file="./koalapop_dellygeno_regions.vcf", 
                           regSizeLowerCutoff=100, regSizeUpperCutoff=30000,
                           readsSupport=3, method="Delly")

# also read in delly rowgranges - use chr and start to merge to intansv so we can retain DEL info
delly.vcf <- VariantAnnotation::readVcf("./koalapop_dellygeno_regions.vcf") ## multiple dim: 371 430
delly.rowgranges <- rowRanges(delly.vcf) # 371
delly.GT <- geno(delly.vcf)$GT

DEL.ref <- data.frame(chr=delly.rowgranges@seqnames, start=delly.rowgranges@ranges@start, DEL=names(delly.rowgranges))
delly.intansv.df <- as.data.frame(delly.intansv$del)
names(delly.intansv.df) <- c("chr", "start", "end", "size", "SU")
delly.intansv.df$SU <- as.numeric(gsub("SU=", "", delly.intansv.df$SU))
delly.intansv.df <- merge(delly.intansv.df, DEL.ref, by=c("chr", "start"))
delly.intansv.df <- delly.intansv.df[order(delly.intansv.df$chr),]
delly.intansv.df[duplicated(delly.intansv.df$chrstart),]
delly.intansv.ranges <- toGRanges(delly.intansv.df[,c("chr", "start", "end")])
mcols(delly.intansv.ranges) <- delly.intansv.df[,c("SU","DEL")]
delly.intansv.ranges$locus <- paste(delly.intansv.ranges)
delly.intansv.GT <- delly.GT[rownames(delly.GT) %in% delly.intansv.ranges$DEL,]

dim(delly.intansv.GT) # 66 regions

hist(delly.intansv.ranges@ranges@width, breaks =50)
delly.intansv.df$erv<- NA
delly.intansv.df$erv.size <- NA
delly.intansv.df$range <- NA
delly.intansv.df$clade <- NA
delly.intansv.df$method <- NA

for(i in 1:nrow(delly.intansv.df)){
  temp <- subsetByOverlaps(reference.erv.candidates.threeclades.final.MTS, delly.intansv.ranges[i,])
  if (length(temp) > 0) {
    delly.intansv.df$erv[i] <- temp$erv
    delly.intansv.df$erv.size[i] <- temp$width
    delly.intansv.df$clade[i] <- temp$clade
    delly.intansv.df$method[i] <- temp$method
    delly.intansv.df$range[i] <- paste(temp@ranges)
    delly.intansv.df$intersect[i] <- intersect(reference.erv.candidates.threeclades.final.MTS, delly.intansv.ranges[i,])@ranges@width
  }
  
}

delly.intansv.df$intersect.match <- delly.intansv.df$intersect - delly.intansv.df$size
delly.intansv.df$intersect.match.pct <- delly.intansv.df$intersect / delly.intansv.df$size

delly.intansv.df <- delly.intansv.df[!is.na(delly.intansv.df$clade),]
# retain those deletions that overlap at least 80% with ERVs in the reference genome 
delly.intansv.df <- delly.intansv.df[delly.intansv.df$intersect.match.pct >= 0.8,]
delly.intansv.ranges <- toGRanges(delly.intansv.df[,c("chr", "start", "end")])
mcols(delly.intansv.ranges) <- delly.intansv.df[,c("SU","DEL")]
delly.intansv.ranges$locus <- paste(delly.intansv.ranges)
delly.intansv.GT <- delly.GT[rownames(delly.GT) %in% delly.intansv.ranges$DEL,]

delly.intansv.df$size.match <- delly.intansv.df$erv.size-delly.intansv.df$size
delly.intansv.GT <- delly.intansv.GT[delly.intansv.df$DEL,]

delly.intansv.GT.geno <- delly.intansv.GT
delly.intansv.GT.geno[delly.intansv.GT.geno == "0/0"] <- 0
delly.intansv.GT.geno[delly.intansv.GT.geno == "0/1"] <- 1
delly.intansv.GT.geno[delly.intansv.GT.geno == "1/1"] <- 2
delly.intansv.GT.geno[delly.intansv.GT.geno == "./."] <- NA
delly.intansv.GT.geno.summary <- data.frame(locus=rownames(delly.intansv.GT.geno), 
                                            hom0=NA, 
                                            het=NA, 
                                            hom1=NA, 
                                            NAs=NA)
for (i in 1:nrow(delly.intansv.GT.geno)){
  delly.intansv.GT.geno.summary$hom0[i] <- length(delly.intansv.GT.geno[i,][delly.intansv.GT.geno[i,] == 0])
  delly.intansv.GT.geno.summary$het[i] <- length(delly.intansv.GT.geno[i,][delly.intansv.GT.geno[i,] == 1])
  delly.intansv.GT.geno.summary$hom1[i] <- length(delly.intansv.GT.geno[i,][delly.intansv.GT.geno[i,] == 2])
  delly.intansv.GT.geno.summary$NAs[i] <- length(delly.intansv.GT.geno[i,][is.na(delly.intansv.GT.geno[i,])])
}
delly.intansv.GT.geno.summary$sum <- delly.intansv.GT.geno.summary$het+delly.intansv.GT.geno.summary$hom0+delly.intansv.GT.geno.summary$hom1

delly.intansv.GT.geno.summary$erv <- delly.intansv.df$erv[match(delly.intansv.GT.geno.summary$locus, delly.intansv.df$DEL)]
delly.intansv.GT.geno.summary$clade <- delly.intansv.df$clade[match(delly.intansv.GT.geno.summary$locus, delly.intansv.df$DEL)]
delly.intansv.GT.geno.summary$loc.clade <- paste(delly.intansv.GT.geno.summary$locus, "_", 
                                                 delly.intansv.df$chr[match(delly.intansv.GT.geno.summary$locus, delly.intansv.df$DEL)], ":", 
                                                 delly.intansv.df$start[match(delly.intansv.GT.geno.summary$locus, delly.intansv.df$DEL)], "-", 
                                                 delly.intansv.df$end[match(delly.intansv.GT.geno.summary$locus, delly.intansv.df$DEL)], "_",
                                                 delly.intansv.GT.geno.summary$clade, sep="")




#### DELETIONS | combine DELLY and LUMPY results  ####
### polarise count info and join
# delly gives genotype data:
# delly if  "1/1" does not have ERV 0
#       if  "1/0" has ERV
#       if  "0/0" has ERV
# lumpy is presence absence deletions, so |call-1| will polarise
delly.intansv.GT.ervcalls <- delly.intansv.GT
delly.intansv.GT.ervcalls[delly.intansv.GT.ervcalls == "0/0"] <- 1
delly.intansv.GT.ervcalls[delly.intansv.GT.ervcalls == "0/1"] <- 1
delly.intansv.GT.ervcalls[delly.intansv.GT.ervcalls == "1/1"] <- 0
delly.intansv.GT.ervcalls[delly.intansv.GT.ervcalls == "./."] <- NA
delly.intansv.GT.ervcalls.mat <- matrix(as.numeric(delly.intansv.GT.ervcalls), nrow=nrow(delly.intansv.GT.ervcalls), 
                                        ncol=ncol(delly.intansv.GT.ervcalls))
row.names(delly.intansv.GT.ervcalls.mat) <- delly.intansv.GT.geno.summary$loc.clade
colnames(delly.intansv.GT.ervcalls.mat) <- colnames(delly.intansv.GT.ervcalls)


lumpy.erv.loci.call <- as.data.frame(mcols(lumpy.erv.loci.counts))
rownames(lumpy.erv.loci.call) <- paste(lumpy.erv.loci.counts)
lumpy.erv.loci.call.mat <- as.matrix(lumpy.erv.loci.call)
lumpy.erv.loci.call.mat <- abs(lumpy.erv.loci.call.mat - 1)



### delly / lumpy overlaps
delly.lumpy.overlaps <- subsetByOverlaps(delly.intansv.ranges, lumpy.erv.df)
temp.hits<- findOverlaps(delly.intansv.ranges, lumpy.erv.df)
lumpy.locus <- unique(CharacterList(split(lumpy.erv.df$locus[subjectHits(temp.hits)],
                                          queryHits(temp.hits))))
lumpy.width <- unique(CharacterList(split(lumpy.erv.df@ranges@width[subjectHits(temp.hits)],
                                          queryHits(temp.hits))))
mcols(delly.lumpy.overlaps) <- DataFrame(mcols(delly.lumpy.overlaps), lumpy.locus, lumpy.width)

for(i in 1:length(delly.lumpy.overlaps)){
  delly.lumpy.overlaps$intersect[i] <- intersect(delly.lumpy.overlaps[i,], lumpy.erv.df)@ranges@width
}
delly.lumpy.overlaps$intersect.lumpy.pct <- delly.lumpy.overlaps$intersect/as.numeric(delly.lumpy.overlaps$lumpy.width)
delly.lumpy.overlaps$intersect.delly.pct <- delly.lumpy.overlaps@ranges@width/as.numeric(delly.lumpy.overlaps$lumpy.width)


### delly / no lumpy overlaps 
delly.NOlumpy.overlaps <- delly.intansv.ranges[!delly.intansv.ranges$DEL %in% delly.lumpy.overlaps$DEL]

### lumpy / no delly overlaps 
lumpy.NOdelly.overlaps <- lumpy.erv.df[!lumpy.erv.df$locus %in% unlist(delly.lumpy.overlaps$lumpy.locus),]


## manual visualisation and curation 

loci <- NULL
for (i in 1:length(delly.lumpy.overlaps)){
  loci <- c(loci, unlist(delly.lumpy.overlaps$DEL[i]), as.character(unlist(delly.lumpy.overlaps$lumpy.locus[i])))
}
loci <- c(loci, unlist(delly.NOlumpy.overlaps$DEL), unlist(lumpy.NOdelly.overlaps$locus))
delly.lumpy.overlaps
delly.intansv.GT.ervcalls.mat.cp <- delly.intansv.GT.ervcalls.mat
lumpy.erv.loci.call.mat.cp <- lumpy.erv.loci.call.mat
heatmap(delly.intansv.GT.ervcalls.mat.cp)
heatmap(lumpy.erv.loci.call.mat.cp)






#### RETROSEQ | load vcfs ####
## load retroseq vcfs and convert each into grange object
koala.refgen.retroseq.vcfs <- paste("./retroseq_vcfs/", list.files(path = "./retroseq_vcfs/", 
                                                                   pattern = "*vcf"), sep = "")
for (i in 1:length(koala.refgen.retroseq.vcfs)){
  name <- gsub("_Retroseq", "", gsub("ERVlib.*", "", basename(koala.refgen.retroseq.vcfs[i])))
  temp.vcf <- VariantAnnotation::readVcf(koala.refgen.retroseq.vcfs[i])
  temp.granges <- rowRanges(temp.vcf)
  temp.granges$erv <- lapply(temp.vcf@info$MEINFO, `[[`, 1)
  temp.granges$erv <- gsub("HERV\\-","HERV_",temp.granges$erv)
  temp.granges$GQ <- as.numeric(geno(temp.vcf)$GQ)
  temp.granges$SP <- as.numeric(geno(temp.vcf)$SP)
  temp.granges$CLIP3 <-as.numeric(geno(temp.vcf)$CLIP3)
  temp.granges$CLIP5 <- as.numeric(geno(temp.vcf)$CLIP5)
  temp.granges$FL <- as.numeric(geno(temp.vcf)$FL)
  names(temp.granges) <- NULL
  temp.granges$ALT <- NULL
  temp.granges$FILTER <- NULL
  temp.granges$paramRangeID <- NULL
  temp.granges$sample <- name
  if (sum(as.integer(!temp.granges@seqnames %like% "NW")) != 0) {
    temp.granges$genbank.accn <- temp.granges@seqnames
    temp.granges <- merge(temp.granges, contig.ref[,c("genbank.accn","refseq.accn")], by="genbank.accn")
    temp.granges$seqnames <- as.character(temp.granges$refseq.accn)
    temp.granges <- toGRanges(temp.granges[,c("seqnames","start","end", "REF","QUAL","erv","GQ","SP","CLIP3","CLIP5","FL","sample")])
  } else {
    temp.granges <- toGRanges(temp.granges)
  }
  assign(paste(name, ".retroseq_ref.range", sep = ""), temp.granges)
  rm(temp.vcf)
  rm(temp.granges)
}

#### RETROSEQ | define candidate loci ####
## compile retroseq vcfs | filter | reduce to candidate ERV loci 

# get retroseq ranges
obj.list = as.list(ls()[sapply(mget(ls(), .GlobalEnv), is.object)])
retroseq.granges.list <- obj.list[obj.list %like% ".retroseq_ref.range" ]

# compile
outgrange <- get(retroseq.granges.list[[1]])
for (i in 2:length(retroseq.granges.list)){
  tempgrange <- get(retroseq.granges.list[[i]])
  outgrange <- c(outgrange, tempgrange)
}

outgrange.fl8.clips2 <- outgrange[outgrange$FL == 8 & outgrange$CLIP3 >=2 & outgrange$CLIP5 >=2,]
i=1
temp <- outgrange.fl8.clips2[outgrange.fl8.clips2$sample == unique(outgrange.fl8.clips2$sample)[i]]
bam.cov <- koala.cov$cov[koala.cov$sample == gsub("_pruned.*", "", unique(outgrange.fl8.clips2$sample)[i]) ]
temp <- temp[temp$GQ >= bam.cov/2 & temp$GQ <= 3*bam.cov]
koala.fl8.clips2.adaptiveGQfilter <- temp
for (i in 2:length(unique(outgrange.fl8.clips2$sample))){
  temp <- outgrange.fl8.clips2[outgrange.fl8.clips2$sample == unique(outgrange.fl8.clips2$sample)[i]]
  bam.cov <- koala.cov$cov[koala.cov$sample ==  gsub("_pruned.*", "", unique(outgrange.fl8.clips2$sample)[i]) ]
  temp <- temp[temp$GQ >= bam.cov/2 & temp$GQ <= 3*bam.cov]
  koala.fl8.clips2.adaptiveGQfilter <- c(koala.fl8.clips2.adaptiveGQfilter, temp)
}
koala.fl8.clips2.adaptiveGQfilter 
# 833774 ranges ~ loop takes 15 mins

koala.fl8.clips2.adaptiveGQfilter.pad100red <- GenomicRanges::reduce(resize(koala.fl8.clips2.adaptiveGQfilter, 
                                                                            width = width(koala.fl8.clips2.adaptiveGQfilter)+(100), fix = "center"))
koala.fl8.clips2.adaptiveGQfilter.pad100red 
# 66336 ranges



#### RETROSEQ | assign ERV to candidate loci ##### granges list
# prep out dataframe:
koala.fl8.clips2.adaptiveGQfilter.pad100red$locus <- paste(koala.fl8.clips2.adaptiveGQfilter.pad100red)
koala.fl8.clips2.adaptiveGQfilter.pad100red.out <- data.frame(locus = paste(koala.fl8.clips2.adaptiveGQfilter.pad100red))

for (i in 1:length(retroseq.granges.list)){
  sample.grange <- get(retroseq.granges.list[[i]])
  name <- unique(sample.grange$sample)
  bam.cov <- koala.cov$cov[koala.cov$sample == gsub("_pruned.*", "",  name) ]
  sample.grange <- sample.grange[sample.grange$GQ >= bam.cov/2 & sample.grange$GQ <= 3*bam.cov & sample.grange$FL == 8 & sample.grange$CLIP3 >=2 & sample.grange$CLIP5 >=2,]
  temp.overlap.loci <- subsetByOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red, sample.grange)
  temp.overlap.samp <- subsetByOverlaps(sample.grange, temp.overlap.loci)
  temp.overlap.samp$locus <- paste(temp.overlap.samp)
  temp.overlap.samp <- temp.overlap.samp[!duplicated(temp.overlap.samp$locus),]
  # removing duplicated regions and choosing those with the highest filter values
  temp.overlap.samp.mcols <- as.data.frame(mcols(temp.overlap.samp)) %>%
    dplyr::arrange(locus, desc(FL)) %>% 
    group_by(locus) %>%
    top_n(2, abs(FL)) %>%
    top_n(1, abs(GQ)) %>%
    as.data.frame()
  temp.overlap.samp.mcols <- temp.overlap.samp.mcols[!duplicated(temp.overlap.samp.mcols$locus),]
  temp.overlap.samp.mcols$locus.erv <- paste(temp.overlap.samp.mcols$locus, temp.overlap.samp.mcols$erv, sep="_")
  temp.overlap.samp$locus.erv <- paste(temp.overlap.samp$locus, temp.overlap.samp$erv, sep="_")
  temp.overlap.samp <- temp.overlap.samp[temp.overlap.samp$locus.erv %in% temp.overlap.samp.mcols$locus.erv,]
  
  temp.overlap.loci <- subsetByOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red, temp.overlap.samp)
  temp.hits<- findOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red, temp.overlap.samp)
  erv.col <- unique(CharacterList(split(temp.overlap.samp$erv[subjectHits(temp.hits)],
                                        queryHits(temp.hits))))
  site.col <- unique(CharacterList(split(temp.overlap.samp@ranges@start[subjectHits(temp.hits)],
                                         queryHits(temp.hits))))
  mcols(temp.overlap.loci) <- DataFrame(mcols(temp.overlap.loci), erv.col, site.col)
  names(mcols(temp.overlap.loci)) <- c("locus" , paste(name,".erv", sep = ""), paste(name,".site", sep = ""))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out <- merge(koala.fl8.clips2.adaptiveGQfilter.pad100red.out, 
                                                           mcols(temp.overlap.loci), by="locus", all.x=TRUE)
}


# pull out matrix of ervs
koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out[, grep("erv", names(koala.fl8.clips2.adaptiveGQfilter.pad100red.out))]
rownames(koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout) <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out$locus
koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout.mat <-  as.matrix(koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout)

# pull out matrix of sites
koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out[, grep("site", names(koala.fl8.clips2.adaptiveGQfilter.pad100red.out))]
rownames(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout) <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out$locus
koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat <-  as.matrix(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout)


## erv assignment to loci:
koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp <- data.frame(locus = rownames(koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout), 
                                                                           erv.assigned = NA, erv.assigned.freq = NA, 
                                                                           actual.stop = NA, actual.start= NA)
for (i in 1:nrow(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp)){
  erv.list <- unlist(base::strsplit(as.character(unlist(t(koala.fl8.clips2.adaptiveGQfilter.pad100red.ervout.mat[i,1:ncol(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat)]))),"\\-"))
  erv.list[erv.list %like% "KoRV"] <- "KoRVlike"
  erv.table <- as.data.frame(table(unlist(erv.list), exclude = c("hybrid", "unknown", NA)))
  clade.table <- merge(data.frame(erv.id = erv.table$Var1, Freq=erv.table$Freq), erv.clades[,c("erv.id", "clade")], by="erv.id", all.x=T)
  clade.table$clade[clade.table$erv.id %like% "KoRV"] <- "KoRVlike"
  clade.table$clade[is.na(clade.table$clade)] <- as.character(clade.table$erv.id[is.na(clade.table$clade)])
  clade.table <- aggregate(clade.table$Freq ~ clade.table$clade, clade.table, sum)
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned[i] <- paste(erv.table$Var1[erv.table$Freq == max(erv.table$Freq)], collapse="|")
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned.freq[i] <- max(erv.table$Freq) / sum(erv.table$Freq)
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$actual.stop[i] <- max(na.omit(as.numeric(unlist(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat[i,1:ncol(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat)]))))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$actual.start[i] <- min(na.omit(as.numeric(unlist(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat[i,1:ncol(koala.fl8.clips2.adaptiveGQfilter.pad100red.sitesout.mat)]))))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.ervs[i] <- paste(unlist(erv.table$Var1), collapse="|")
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.erv.counts[i] <- paste(unlist(erv.table$Freq), collapse="|")
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.clades[i] <- paste(unlist(clade.table$`clade.table$clade`), collapse="|")
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.clades.count[i] <- paste(unlist(clade.table$`clade.table$Freq`), collapse="|")
  ervs <- unlist(base::strsplit(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned[i], split = "\\|"))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned.clades[i] <- ifelse(length(unique(erv.clades$clade[erv.clades$erv.id %in% ervs])) == 0, 
                                                                                                koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned[i], 
                                                                                                paste(unique(erv.clades$clade[erv.clades$erv.id %in% ervs]), collapse="|"))
  unique.ervs <- unlist(base::strsplit(koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.ervs[i], split = "\\|"))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned.ed[i] <- ifelse(length(unique(erv.clades$clade[erv.clades$erv.id %in% unique.ervs])) == 1 & length(ervs) != 1, 
                                                                                            ervs[1], koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$erv.assigned[i])
}


#### RETROSEQ | get counts for candidate loci #### 

koala.fl8.clips2.adaptiveGQfilter.pad100red
# ran loop for different values for accepting calls:
runs <- data.frame(FL = c(3, 3), 
                   GQ = c(2, 4))

for (n in 1:nrow(runs)) {
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts <- koala.fl8.clips2.adaptiveGQfilter.pad100red
  for (i in 1:length(retroseq.granges.list)){
    g1 <- get(retroseq.granges.list[[i]])
    g1 <- g1[g1$GQ >= runs$GQ[n] & g1$FL >= runs$FL[n],]
    sample <- unique(g1$sample)
    koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts$hits <- countOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts, g1, type = "any")
    koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts$hits[koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts$hits >= 1] <- 1
    names(mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts))[names(mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts)) == "hits"] <- sample
  }
  assign(paste("koala.granges.fl8.clips2.pad100red.counts.gq", runs$GQ[n], "fl" , runs$FL[n], sep = ""), 
         koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts)
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.df <- as.data.frame(mcols(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus <- data.frame(locus = paste(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts), 
                                                                                counts = rowSums(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.df[,2:ncol(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.df)]))
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$chr <- gsub("\\:.*", "", koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$locus)
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$start <- ranges(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts)@start
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$length <- ranges(koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts)@width
  koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$end <- koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$start + koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus$length -1
  assign(paste("koala.granges.fl8.clips2.pad100red.counts.gq", runs$GQ[n], "fl" , runs$FL[n],  ".perlocus", sep = ""), 
         koala.fl8.clips2.adaptiveGQfilter.pad100red.out.counts.perlocus)
}



#### RETROSEQ | filtering loci: ####
# RETROSEQ | filter for three clades that were purely called for that clade. ####
## filter for those loci that were called for KoRVlike, phaCinBeta or phaCinBetalike
koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered <- koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp[koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.clades == "KoRVlike" |
                                                                                                                       koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.clades == "phaCinBeta" |
                                                                                                                       koala.fl8.clips2.adaptiveGQfilter.pad100red.singleervhits.cp$unique.clades == "phaCinBetalike", ]


koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered.gr <- toGRanges(data.frame(chr=gsub( ":.*", "",koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered$locus), 
                                                                                start=koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered$actual.start,
                                                                                end=koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered$actual.stop, 
                                                                                erv=koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered$erv.assigned, 
                                                                                clade=koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered$unique.clades,
                                                                                locus=koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered$locus))

nrow(koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered) # 12982
table(koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered$erv.assigned.clades)
#       KoRVlike     phaCinBeta phaCinBetalike 
#       9313           3165            504 


## RETROSEQ | filter out those loci that overlap with reference ERVs ####
# using margins 250 bp either way from ERV call. 
reference.erv.candidates.threeclades.final.wmargin <- resize(reference.erv.candidates.threeclades.final, 
                                                             width = width(reference.erv.candidates.threeclades.final)+(500), fix = "center")

# find overlaps: 
retroseqloci.refervswmargin.overlap <- subsetByOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered.gr, reference.erv.candidates.threeclades.final.wmargin)
temp.hits<- findOverlaps(koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered.gr, reference.erv.candidates.threeclades.final.wmargin)
erv.ref <- unique(CharacterList(split(reference.erv.candidates.threeclades.final.wmargin$erv[subjectHits(temp.hits)],
                                      queryHits(temp.hits))))
erv.ref.method <- unique(CharacterList(split(reference.erv.candidates.threeclades.final.wmargin$method[subjectHits(temp.hits)],
                                             queryHits(temp.hits))))
erv.ref.clade <- unique(CharacterList(split(reference.erv.candidates.threeclades.final.wmargin$clade[subjectHits(temp.hits)],
                                            queryHits(temp.hits))))

mcols(retroseqloci.refervswmargin.overlap) <- DataFrame(mcols(retroseqloci.refervswmargin.overlap), erv.ref, erv.ref.method, erv.ref.clade)
specialcase.retroseqloci.refervs.overlaps <- retroseqloci.refervswmargin.overlap[unlist(retroseqloci.refervswmargin.overlap$clade != retroseqloci.refervswmargin.overlap$erv.ref.clade)] # 119 overlaps with same clade
# see also section "filters working: true secondary insertions" into phacinbetalike ERVs # 
retroseqloci.refervswmargin.overlap.toremove <- retroseqloci.refervswmargin.overlap[unlist(retroseqloci.refervswmargin.overlap$clade == retroseqloci.refervswmargin.overlap$erv.ref.clade)] # 119 overlaps with same clade

## filter them out
koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered2 <- koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered[!koala.fl8.clips2.adaptiveGQfilter.pad100red.filtered$locus %in% retroseqloci.refervswmargin.overlap.toremove$locus,]


#### JOINED CALL #### 
# manual curation to remove loci that were 0 counts across samples
    # ensure all datasets from retroseq, lumpy and delly have the same names, in the same order 
    # remove delly loci that overlap with lumpy
    # retroseq loci are retricted to those called for only ONE clade; add back in the loci that have two closely located integrations that have been visually confirmed
#result is poly.ervs.joined.mat.fil2
poly.ervs.joined.mat.fil2 <- fread("koala_polyERVs.csv") # 12990 loci used in analysis



### 2023/11/07 Mette Lillie mette.lillie@ebc.uu.se









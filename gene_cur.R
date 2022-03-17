library(magrittr)
library(tidyverse)
library(GenomicRanges)
gene_count <- read.delim("/wrk/yuanzhen/expressRNA/gene_count.txt")
gene_count$id <- rownames(gene_count)
gene_count$dmso <- (gene_count$dmso1 + gene_count$dmso2) / 2
gene_count$iaa <- (gene_count$iaa01 + gene_count$iaa02) / 2
gene_count$foldg <- (log2(gene_count$iaa + 1)) - (log2(gene_count$dmso + 1))
gene_fg <- filter(gene_count,foldg>1 & foldg<3)
gffara <- GenomicFeatures::makeTxDbFromGFF("/wrk/data/genome/yz_genome_data/aragenome/Athaliana.gff3")
gene_ <- GenomicFeatures::genes(gffara)
gene_1kb <- GenomicRanges::GRanges(seqnames=gene_@seqnames,
                                   ranges = IRanges::IRanges(start=gene_@ranges@start-1000,
                                                             end = gene_@ranges@start+1000))
mcols(gene_1kb)$id <- gene_$gene_id
gene_1kb=gene_1kb[IRanges::start(ranges(gene_1kb))>0]

seqlevels(gene_1kb) %<>% str_replace("Chr","chr")
seqlevels(gene_1kb,pruning.mode="coarse")=seqlevels(gene_1kb)[seqlevels(gene_1kb) %>% str_detect("[0-9]$")]

genome="/wrk/data/genome/yz_genome_data/aragenome/TAIR10_chr_all.fas"
genom_ara <- Biostrings::readDNAStringSet(genome)
names(genom_ara) %<>% str_replace("Chr","chr")
seqs <- genom_ara[gene_1kb]
names(seqs) <- gene_1kb$id

###the data from http://cisbp.ccbr.utoronto.ca/bulk_archive.php
tff_info <- read.delim("/wrk/yuanzhen/motifdb/cisBP/TF_Information_all_motifs_plus.txt")
tff_gene <- data.frame(id=tff_info$DBID,name=tff_info$TF_Name) %>% unique()


getkmer <- function(seqfile="/wrk/wenchenjin/work/Chenjin1_1__ATI_ara.cell_1/pre_processed/cleaned/A9_9_IAA__leaf__conc_0.1nM__time_2h__30n.gz",n=1000,kmer_num=9L){
  seqs=fjComm::getSeq_fqfachrFile(seqfile)
  kmer_cnt=fjComm::kmerCntBit(strings =seqs, k = kmer_num, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 10)
  seqs_sf <- stringi::stri_rand_shuffle(seqs)
  kmer_cnt_sf=fjComm::kmerCntBit(strings =seqs_sf, k = kmer_num, diffLen = T, collapse = T, asDf = T, all_possible_k = T, pseudo = 10) 
  kmer_total <- dplyr::left_join(kmer_cnt,kmer_cnt_sf,by="kmer") %>% filter(!counts.x < counts.y & !counts.x == counts.y)
  kmer_total$kmer_nor <- -log2((kmer_total$counts.x - kmer_total$counts.y) / sum(str_count(seqs)))
  kmer_top1000 <- kmer_total %>% dplyr::arrange(kmer_nor) %>% dplyr::slice(1:n)
  Biostrings::DNAStringSet(kmer_top1000$kmer)
}

kmer <- getkmer()
#BiocParallel::bplapply(seq_along(seqs),function(x){
# seq_n <- names(seqs)[[x]]
# test <- Biostrings::countPDict(kmer,seqs[[x]])
# c(seq_n,sum(test))
#},BPPARAM = BiocParallel::MulticoreParam(workers=30,progressbar=TRUE))
data_gen <- data.frame(matrix(NA,length(seqs),2))
for(i in seq_along(seqs)){
  cat(i,"\n")
  seq_n <- names(seqs)[[i]]
  test <- Biostrings::countPDict(kmer,seqs[[i]]) %>% sum()
  data_gen[i,1] <- seq_n
  data_gen[i,2] <- test %>% as.numeric()
}
colnames(data_gen) <- c("id","count")
ex_count <- left_join(data_gen,gene_fg,by="id") %>% na.omit() %>% filter(dmso > 0 & iaa > 0)
#ex_tf_count <- left_join(tff_gene,ex_count,by="id")

test=data.frame(count=ex_count$count,dmso=(ex_count$dmso1+ex_count$dmso2)/2)
test_o <- test[order(-test$dmso),]
test_o <- test_o[-c(1:5),] %>% filter(dmso >500) %>% filter(dmso<600000)
p <- ggplot(ex_count,aes(x=log(count),y=foldg))+geom_point(fill=ex_count$id)+geom_smooth(method = "lm",linetype=3,se=T,colour="black",formula = y ~ x)
plotly::ggplotly(p)






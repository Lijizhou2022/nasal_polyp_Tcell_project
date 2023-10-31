# for bulk-TCR-seq
## Step1 extract and quantify TCR clones for each sample
sample_igblast_file<-as.data.frame(read.table(file = sample_igblast_path,
                                                 sep = "\t",header = T))
#split sample_igblast_file into sample_TRA_igblast_file and sample_TRB_igblast_file according to C gene

#clone quantification
sample_TRA_igblast_file$germline_cdr3<-as.character(sample_TRA_igblast_file$germline_cdr3)

sample_TRA_igblast_file$germline_cdr3<-as.factor(sample_TRA_igblast_file$germline_cdr3)
sample_TRA_germline_cdr3<-levels(sample_TRA_igblast_file$germline_cdr3)
sample_TRA_germline_cdr3<-data.frame(sample_TRA_germline_cdr3)
sample_TRA_igblast_file$germline_cdr3<-as.character(sample_TRA_igblast_file$germline_cdr3)
sample_TRA_germline_cdr3$sample_TRA_germline_cdr3<-as.character(sample_TRA_germline_cdr3$sample_TRA_germline_cdr3)
sample_TRA_germline_cdr3$umi_count<-0
length<-length(sample_TRA_igblast_file$umi_count)
for (i in 1:length(sample_TRA_igblast_file$sequence_id)){ 
  #j=sample_TRA_germline_cdr3$umi_count[which(sample_TRA_germline_cdr3$TRA_germline_cdr3==sample_TRA_igblast_file$germline_cdr3[i])]
  sample_TRA_germline_cdr3$umi_count[which(sample_TRA_germline_cdr3$sample_TRA_germline_cdr3==sample_TRA_igblast_file$germline_cdr3[i])]<-1+sample_TRA_germline_cdr3$umi_count[which(sample_TRA_germline_cdr3$sample_TRA_germline_cdr3==sample_TRA_igblast_file$germline_cdr3[i])]
  #print(TRA_igblast_file$umi_count[i])
}
row.names(sample_TRA_germline_cdr3)<-sample_TRA_germline_cdr3$sample_TRA_germline_cdr3

sample_TRB_igblast_file$germline_cdr3<-as.character(sample_TRB_igblast_file$germline_cdr3)

sample_TRB_igblast_file$germline_cdr3<-as.factor(sample_TRB_igblast_file$germline_cdr3)
sample_TRB_germline_cdr3<-levels(sample_TRB_igblast_file$germline_cdr3)
sample_TRB_germline_cdr3<-data.frame(sample_TRB_germline_cdr3)
sample_TRB_igblast_file$germline_cdr3<-as.character(sample_TRB_igblast_file$germline_cdr3)
sample_TRB_germline_cdr3$sample_TRB_germline_cdr3<-as.character(sample_TRB_germline_cdr3$sample_TRB_germline_cdr3)
sample_TRB_germline_cdr3$umi_count<-0
length<-length(sample_TRB_igblast_file$umi_count)
for (i in 1:length(sample_TRB_igblast_file$sequence_id)){ 
  #j=sample_TRB_germline_cdr3$umi_count[which(sample_TRB_germline_cdr3$TRB_germline_cdr3==sample_TRB_igblast_file$germline_cdr3[i])]
  sample_TRB_germline_cdr3$umi_count[which(sample_TRB_germline_cdr3$sample_TRB_germline_cdr3==sample_TRB_igblast_file$germline_cdr3[i])]<-1+sample_TRB_germline_cdr3$umi_count[which(sample_TRB_germline_cdr3$sample_TRB_germline_cdr3==sample_TRB_igblast_file$germline_cdr3[i])]
  #print(TRB_igblast_file$umi_count[i])
}
row.names(sample_TRB_germline_cdr3)<-sample_TRB_germline_cdr3$sample_TRB_germline_cdr3


## Step2 build the whole clone matrix
Args= list.files("dir of germline-cdr3 files")
for (i in 1:length(Args)) {
  if(i==1){
    whole<-read.table(Args[i],header = 1,row.names = 1)
    next
  }else{
    add<-read.table(Args[i],header = 1,row.names = 1)
    whole<-merge(x = whole,y = add,
                 by.x = 'row.names',by.y = 'row.names',all = T)
    row.names(whole)<-whole[,1]
    whole<-whole[,-1]
  }
}
## The "whole" matrix is the clone count matrix of all clones across all samples. It can be used to downstream analysis such as common clone analysis and  circos graph.  


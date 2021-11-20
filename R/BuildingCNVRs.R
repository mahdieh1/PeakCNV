#' Runs Building CNVRs step.
BuildingCNVRs<- function()
{
  
  map <- bedtools_intersect("-a case.bed -b control.bed ")
  intersect <-eval(map) 
  write.table(intersect,file="intersect.bed", quote=F, sep="\t", row.names=F, col.names=F)
  mapping <-read.table('intersect.bed')
  mapping <-mapping[,c(1,2,3)]
  write.table(mapping,file="mapping.bed", quote=F, sep="\t", row.names=F, col.names=F)
  
  
  count.case <-bedtools_intersect("-a mapping.bed -b case.bed -c")
  case.eval<-eval(count.case) 
  write.table(case.eval,file="overlap_Bin_case.bed", quote=F, sep="\t", row.names=F, col.names=F)
  case<-read.table('overlap_Bin_case.bed')
  case<-case[,c(1,2,3,6)]
  
  count.control<-bedtools_intersect("-a mapping.bed -b control.bed -c")
  control.eval<-eval(count.control) 
  write.table(control.eval,file="overlap_Bin_Control.bed", quote=F, sep="\t", row.names=F, col.names=F)
  control<-read.table('overlap_Bin_Control.bed')
  control<-control[,c(1,2,3,6)]
  
  new.df<-data.frame();
  new.df<-case;
  new.df$V7 <- control$V6;
 
  C<-read.table('case.bed')
  #No.case=length(unique(C$V4));
  No.case=nrow(C);
  C<-read.table('control.bed')
  #No.control=length(unique(C$V4));
  No.control=nrow(C);
  new.df[,6]=No.case-new.df[,4];
  new.df[,7]=No.control-new.df[,5];
  
  new.df<-new.df[!duplicated(new.df[ , c("V1", "V2","V3")]), ]
  for(i in 1:nrow(new.df))
  {
    a <- rbind( c(new.df[i,4],new.df[i,6]), c(new.df[i,5],new.df[i,7]) );
    fisher.result <- fisher.test(x=a,  alternative = "greater")
    print(i)
    print(fisher.result$p.value)
    new.df[i,8]<-fisher.result$p.value;
  }
  
  
  #Apply filter threshold
  cat("Enter thrshold for p-value : \n")
  threshold <- as.integer(readline(prompt = ""))
  #threshold=0.05;
  CNVRs <- new.df[new.df[,8]< threshold,]
  CNVRs<-CNVRs[,c(1,2,3)]
  CNVRs<-CNVRs[!duplicated(CNVRs), ]
  CNVRs<-CNVRs %>% filter(CNVRs$V3 > CNVRs$V2)

  write.table(CNVRs,file="CNVRs.bed", quote=F, sep="\t", row.names=F, col.names=F)
  #Bins<-bedtools_makewindows("-b CNVRs.bed -w 1") 
  #Bins.eval<-eval(Bins) 
  #write.table(Bins.eval,file="Step1.bed", quote=F, sep="\t", row.names=F, col.names=F)
  #CNVR.1<-read.table('Step1.bed')
  #CNVR.1<-CNVR.1[,c(3,4,5)]
  #write.table(CNVR.1,file="Step1.bed", quote=F, sep="\t", row.names=F, col.names=F)
  #overlapping
  count.case1<-bedtools_intersect("-a CNVRs.bed -b case.bed -c")
  case.eval1<-eval(count.case1) 
  write.table(case.eval1,file="overlap_Bin_Case1.bed", quote=F, sep="\t", row.names=F, col.names=F)
  case1<-read.table('overlap_Bin_Case1.bed')
  case1<-case1[,c(1,2,3,6)]
  case1$V2=case1$V2-1;
  write.table(case1,file="overlap_Bin_Case.bed", quote=F, sep="\t", row.names=F, col.names=F)
  
  #generate new files
  intersect_all <- bedtools_intersect("-a CNVRs.bed -b case.bed -wa -wb")
  intersect.eval<-eval(intersect_all) 
  write.table(intersect.eval,file="intersect.bed", quote=F, sep="\t", row.names=F, col.names=F)
  intersect1<-read.table('intersect.bed')
  intersect1<-intersect1[,c(1,2,3,6,7,8,11)]
  intersect1$V2=intersect1$V2-1;
  write.table(intersect1,file="intersect.bed", quote=F, sep="\t", row.names=F, col.names=F)
}

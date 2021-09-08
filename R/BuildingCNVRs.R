#' Runs Building CNVRs step.
BuildingCNVRs<- function()
{
  # if (!requireNamespace("BiocManager", quietly = TRUE))
  #   install.packages("BiocManager")
  # BiocManager::install("HelloRanges")
  #
  # library(HelloRanges)

  #setwd(system.file("unitTests", "data", "makewindows", package="HelloRanges"))

  #Window <- bedtools_makewindows("-g test.genome -w 20")
  #bin<-eval(Window)
  #write.table(bin, file="window.bed", quote=F, sep="\t", row.names=F, col.names=F)
  #df<-read.table('window.bed')
  #df1<-df[,c(3:5)]
  #write.table(df1,file="window.bed", quote=F, sep="\t", row.names=F, col.names=F)


  count.case<-bedtools_intersect("-a window.bed -b case.bed -c")
  case<-eval(count.case)
  write.table(case,file="overlap_Bin_Case.bed", quote=F, sep="\t", row.names=F, col.names=F)
  case<-read.table('overlap_Bin_Case.bed')
  case<-case[,c(1,2,3,6)]

  count.control<-bedtools_intersect("-a window.bed -b control.bed -c")
  control<-eval(count.control)
  write.table(control,file="overlap_Bin_Control.bed", quote=F, sep="\t", row.names=F, col.names=F)
  control<-read.table('overlap_Bin_Control.bed')
  control<-control[,c(1,2,3,6)]


  new.df<-data.frame();
  new.df<-case;
  new.df$V7 <- case$V6;
  #check
  C<-read.table('case.bed')
  #enter from input
  No.case=length(unique(C$V4));
  C<-read.table('control.bed')
  No.control=length(unique(C$V4));

  new.df[,6]=No.case-new.df[,4];
  new.df[,7]=No.control-new.df[,5];

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
  CNVRs <- new.df[new.df[,8]< threshold,]
  CNVRs<-CNVRs[,c(1,2,3)]
  write.table(CNVRs,file="CNVRs.bed", quote=F, sep="\t", row.names=F, col.names=F)

  #overlapping
  count.case<-bedtools_intersect("-a CNVRs.bed -b case.bed -c")
  case<-eval(count.case)
  write.table(case,file="overlap_Bin_Case.bed", quote=F, sep="\t", row.names=F, col.names=F)
  case<-read.table('overlap_Bin_Case.bed')
  case<-case[,c(1,2,3,6)]
  case$V2=case$V2-2;
  write.table(case,file="overlap_Bin_Case.bed", quote=F, sep="\t", row.names=F, col.names=F)

  #generate new files
  intersect_all <- bedtools_intersect("-a CNVRs.bed -b case.bed -wa -wb")
  intersect<-eval(intersect_all)
  write.table(intersect,file="intersect.bed", quote=F, sep="\t", row.names=F, col.names=F)
  intersect1<-read.table('intersect.bed')
  intersect1<-intersect1[,c(1,2,3,6,7,8,11)]
  intersect1$V2=intersect1$V2-1;
  write.table(intersect1,file="intersect.bed", quote=F, sep="\t", row.names=F, col.names=F)
}

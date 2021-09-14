#' Runs Clustering step.
Clustering<-function()
{
  Chr.cluster<-data.frame();
  win<-read.table("overlap_Bin_Case.bed")#,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
 # win$V2=win$V2-2;#remove it puts in building cnvrs
  #win<-win[,c(1,2,3,6)]
  f<-unique(win$V1)
  
  for(k in f)
  
  { 
    total<-data.frame();
    df<-data.frame();
    row<-data.frame();
    col<-data.frame();
    win_sel <- win[win[,1]==k,];
    col2 <- win_sel$V2;
    row<-as.data.frame(t(col2));
    col<-as.data.frame(win_sel$V3);
    for(i in 1:nrow(col))
    {
      for(j in 1:ncol(row))
      {
        df[i,j]<- abs(col[i,1]-row[1,j])
      }
    }
    
    subtract<-read.table(file=paste(k, "subtract.txt", sep=""),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
    total<-cbind(df,subtract)
    if(nrow(total)<5){
      rate=1;
    }
    else
    {
      rate=5;
    }
    iris_matrix <- as.matrix(total)
    kNNdistplot(iris_matrix, rate)
    cat("Enter value for eps : \n")
    threshold <- as.integer(readline(prompt = ""))
    abline(h=threshold, col="red")
    
    db = dbscan(iris_matrix, threshold, rate)
    c<-db$cluster 
    c<-data.frame(win_sel,c) 
    write.table(c, file=paste(k, "DBSCAN.txt", sep=""), sep="\t", col.names = F, row.names = F)
    Chr.cluster<-rbind(Chr.cluster,c)
  }
  write.table(Chr.cluster, file="clustering.txt", sep="\t", col.names = F, row.names = F)
}

#' Generate Input matrix for Clustering step.

InputMatrix<-function()
{
  df <- read.table("intersect.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  No <- read.table("overlap_Bin_Case.bed")#,header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  Number<-unique(df$V1)
  #No$V2=No$V2-2;#remove it puts in building cnvrs
  #No<-No[,c(1,2,3,6)]#remove it puts in building cnvrs

  for(k in Number)
  {
    print( paste("chromosome: ", k)) 
    row<-data.frame();
    col<-data.frame();
    f<-data.frame();
    subtract<-data.frame();
    df1=df[df[,1]==k,]
    cnv1<-No[,1:3]
    cnv1<-cnv1[cnv1[,1]==k,]
    No1=No[No[,1]==k,]
    print('phase1')
    for (j in (1:nrow(cnv1)))
    {
      col[j,1]<-cnv1[j,2];
      row[1,j]<-cnv1[j,2];
    }
    print('phase2')
    for(i in 1:nrow(col))
    {
      f[i,i]<-No1[i,4];
      subtract[i,i]<-0;
    }
    print('phase3')
    for(i in(1:nrow(col)))
    {
      for(j in (1:ncol(row)))
      {
        if(i!=j)

        {
          c<-df1[df1$V2==col[i,1],];
          #print('length c');
          d<-df1[df1$V2==row[1,j],];
          #print('length d');
          total <- rbind(c,d);
          f[i,j]<-sum(duplicated(total$V7));
          print(f[i,j])
          #print('ok')
          subtract[i,j]<-abs(f[i,i]-f[i,j]);
          print(subtract[i,j])
          print(i);
          print(j);
          print(sum(!duplicated(total$V7)));
        }
      }
    }

    print('writing file')
    write.table(subtract, file=paste(k, "subtract.txt", sep=""), sep="\t", col.names = F, row.names = F)

  }
}


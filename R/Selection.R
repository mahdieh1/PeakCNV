#' Select the best CNVRs within cluster .
Selection<-function()
{
  data<-read.table('clustering.txt')
type<-unique(data$V1)
CNVRs.final<-data.frame()
for(chr in type)
{
  #chr=22
  df1<-data[data$V1==chr,]
  df<-read.table(file=paste(chr, "subtract.txt", sep=""),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
  finalscore<-data.frame();
  No_clus<-unique(df1$V5);
  for (i in No_clus) 
  {
    #i=1
    print(i)
    df2<-cbind(df,df1$V5)
    f<-as.integer(rownames(df2[df2$`df1$V5`==i,]))
    #if(is.null(f)!=1)
    d<-df2[f,f];
    print(nrow(d))
    if(length(d)==1){
      #if(nrow(d)==NULL){
      s<-data.frame();
      #s<-df2[df2$`df1$V5`==i,]
      #single<-cbind(single,s)
      s[1,1]<-0;
      finalscore<-rbind(finalscore,s)
    }
    
    else{
      print('cluster is not empty')
      uniuque<-data.frame();
      for (k in 1:nrow(d)) {
        sum<-0;
        for (j in 1:nrow(d)) {
          if(k!=j){
            f3 <-d[k,];
            f4<- d[j,];
            f3 <-as.numeric(f3);
            f4 <-as.numeric(f4);
            #res <- cor(f3, f4,method = "kendall")#,use = "complete.obs");
            res <- cor.test(f3, f4,method = "kendall")
            print('Correlation')
            print(res$p.value)
            sum<-sum+res$estimate;
            uniuque[k,1]<-sum;
          }
        }
        
      }
      #uniuque <-as.numeric(uniuque);
      score<-data.frame();
      #finalscore<-data.frame();
      #score = rep(1, nrow(uniuque))
      for (n in 1:nrow(uniuque)) 
      {
        f<-(sum(d[n,])-(uniuque[n,1]));
        score[n,1]<-f
      }
      
      finalscore<-rbind(finalscore,score)
    }  
    #   else  {
    #   s<-data.frame();
    #   #s<-df2[df2$`df1$V5`==i,]
    #   #single<-cbind(single,s)
    #   s[1,1]<-40
    #   finalscore<-rbind(finalscore,s)
    # 
    # }
  }
  colnames(finalscore) <- c("score")
  df1<- cbind(finalscore,df1)
  #df1<-df1[,c(8,9,10,11,12,13)]
  final.merge<-data.frame();
  #final.merge1<-dat.frame();
  best_region<-data.frame();
  for (m in No_clus) {
    #m=0
    best_region<-data.frame();
    cluster=df1[df1$V5==m,]
    if(nrow(cluster)>1){
      print('sorting')
      sorted_cluster <- cluster[order(-cluster$score),];
      best_region=head(sorted_cluster,n=1)
      #best_region=sorted_cluster[1,];
      for( h in 2:nrow(sorted_cluster)){
        if ((best_region[1,4]-sorted_cluster[h,3] ) < 1000)
        {
          print('merging');
          #print(sorted_cluster[k,2]- best_region[1,3])
          #print(sorted_cluster[k,2])
          best_region$V3= sorted_cluster[h,4];
          best_region$V2=best_region$V2;
          #best_region$score
          print('start')
          print(best_region$V2)
          print('end')
          print(best_region$V3)
          #final.merge<-best_region;
          
        }
        else
        {
          print('no merge')
        }
        #final.merge<-rbind(best_region,final.merge);
      }
      
    } 
    #else{
    #  print('check')
    #  final.merge<-rbind(cluster,final.merge)
    #}
    #final.merge<-rbind(final.merge,final.merge1);
    #rbind(gain.final,final.merge)
    #write.table(final, file=paste(i, "loss.txt", sep=""), sep="\t", col.names = F, row.names = F)
    final.merge<-rbind(best_region,final.merge)
  }
  write.table(final.merge, file=paste(chr, "selection.txt", sep="") , sep="\t", col.names = F, row.names = F)
  CNVRs.final<-rbind(CNVRs.final,final.merge)
}
#type<-unique(data$V1)
#CNVR.final1<-data.frame();
#for (chromosome in type)
#{
#  print(i)
 # selection<-read.table(file=paste(i, "selection", sep=""),header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
 # CNVR.final1<-rbind(CNVR.final1,selection)
  
#}
write.table(CNVRs.final, file= "FinalCNVRs.txt" , sep="\t", col.names = F, row.names = F)
}

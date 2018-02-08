#name_func is used to find the candidates (to provide the expression data)
name_func<-function(name_base, mirna_base, mrna_base){
  x<-match(name_base[,2],mirna_base[1,])
  y<-match(name_base[,3],mrna_base[1,])
  return(list(x,y))
}
#mirna_matrix is used to get the original mirna expression data
mirna_matrix<-function(name_base, mirna_base, mirna_name, mirna){  
  mirna_name1<-rep(0,nrow(mirna_base))
  for(i in 1:nrow(mirna_base)){
    mirna_name1[i]<-unlist(strsplit(mirna_base[i,1],"\\|"))[1]
  }
  mirna_use<-mirna_base[mirna_name1%in%mirna,mirna_name]
  mirna_name1<-mirna_name1[mirna_name1%in%mirna]
  return(list(mirna_use, mirna_name1))
}
#mrna_matrix is used to get the original mrna expression data
mrna_matrix<-function(name_base, mrna_base, mrna_name, mrna){
  mrna_exist<-rep(0,nrow(mrna_base))
  mrna_name_all<-matrix(data=0,ncol=2,nrow=nrow(mrna_base),byrow=TRUE)
  mrna_name1<-rep(0,nrow(mrna_base))
  
  for(i in 1:nrow(mrna_base)){
    mrna_name1[i]<-unlist(strsplit(mrna_base[i,1],"\\|"))[1]
    mrna_name_all[i,1]<-unlist(strsplit(mrna_base[i,1],"\\|"))[1]
    mrna_name_all[i,2]<-unlist(strsplit(mrna_base[i,1],"\\|"))[2]
    if(mrna_name1[i]%in%mrna==1){
      mrna_exist[i]=i
    }
  }
  mrna_use<-mrna_base[mrna_exist,mrna_name]
  mrna_name1<-mrna_name1[mrna_exist]
  mrna_name_sp<-mrna_name_all[mrna_exist,]
  mrna_fullname<-mrna_base[mrna_exist,1]
  return(list(mrna_use, mrna_name1, mrna_name_sp, mrna_fullname))
}
#when the number of NA is over 0.5, the mirna or mrna is deleted from the list, mirna_mrna_data is used to get the filtered data
mirna_mrna_data_unselected<-function(mirna_use, mrna_use, mirna_name, mrna_name, mrna_name_sp, mrna_fullname){
  mirna_use[is.na(mirna_use)]<-0
  mrna_use[is.na(mrna_use)]<-0
  mirna_use[mirna_use<0]<-0
  mrna_use[mrna_use<0]<-0
  mirna_sgn<-seq(1,nrow(mirna_use),1)
  mrna_sgn<-seq(1,nrow(mrna_use),1)
  for (i in 1:nrow(mirna_use)){
    if(sum(mirna_use[i,]==0)>0.9999*ncol(mirna_use)){
      mirna_sgn[i]=0
    }
  }
  for (i in 1:nrow(mrna_use)){
    if(sum(mrna_use[i,]==0)>0.99999*ncol(mrna_use)){
      mrna_sgn[i]=0
    }
  }  
  mirna_use1<-mirna_use[mirna_sgn,]
  mrna_use1<-mrna_use[mrna_sgn,]
  mirna_name1<-mirna_name[mirna_sgn]
  mrna_name1<-mrna_name[mrna_sgn]
  mrna_name_sp1<-mrna_name_sp[mrna_sgn,]
  mrna_fullname1<-mrna_fullname[mrna_sgn]
  rownames(mirna_use1)<-mirna_name1
  rownames(mrna_use1)<-mrna_fullname1
  return(list(mirna_use1, mrna_use1, mirna_name1, mrna_name1, mrna_name_sp1, mrna_fullname1))
}
#mirna_mrna_loc is used to get the wMRE matrix that fits the expression data
mirna_mrna_loc<-function(mirna_name, mrna_name, mirna, mrna, wMRE, mrna_fullname){
  mirna_loc<-match(mirna_name,mirna)
  mrna_loc<-match(mrna_name,mrna)
  wMRE<-t(wMRE)
  wMRE_use<-wMRE[mirna_loc, mrna_loc]
  rownames(wMRE_use)<-mirna_name
  colnames(wMRE_use)<-mrna_fullname
  return(list(mirna_loc, mrna_loc, wMRE_use))
}
#sum_wepro is used to calculate the denominator of the wepro method
sum_wepro<-function(mirna,mrna,wMRE){
  res<-rep(0,ncol(mirna))
  for(i in 1:ncol(mirna)){
    mirna_use1<-as.numeric(mirna[,i])
    mrna_use1<-as.numeric(mrna[,i])
    temp=mirna_use1%*%t(mrna_use1)
    temp2<-wMRE*temp
    res[i]=sum(temp2)
  }
  return(res)
}
#wepro_rankresult is used to calculate the result of wepro model, output[1] is rank matrix
wepro_rankresult<-function(mirna_use1, mrna_use1, wMRE, sumup, ID, thres_num){
  numm<-floor(length(wMRE[wMRE!=0]))/100
  mirna_use<-matrix(as.numeric(mirna_use1),nrow=nrow(mirna_use1))
  mrna_use<-matrix(as.numeric(mrna_use1),nrow=nrow(mrna_use1))
  #pro_matrix_temp saves the probability matrix of each single sample
  pro_matrix_temp<-matrix(data=0,ncol=nrow(mrna_use),nrow=nrow(mirna_use),byrow=TRUE)
  #pro_matrix_sumup saves the probability matrix of the whole cancer
  #pro_matrix_sumup<-matrix(data=0,ncol=nrow(mrna_use),nrow=nrow(mirna_use),byrow=TRUE)
  #pro_matrix_borda is the total rank of the pairs
  pro_matrix_borda<-rep(0,nrow(mrna_use)*nrow(mirna_use))
  uni<-c()
  sum<-nrow(mrna_use)*nrow(mirna_use)
  inter<-seq(1, sum, 1)
  tot<-0
  sample_num<-matrix(data=0,ncol=ncol(mirna_use),nrow=1000,byrow=TRUE)
  #sort_one_rank saves the position of pairs that are ordered by prob
  sort_one_rank<-matrix(data=0,ncol=ncol(mirna_use),nrow=thres_num,byrow=TRUE)
  #sort_one_prob is the probability number that are ordered by prob
  sort_one_prob<-matrix(data=0,ncol=ncol(mirna_use),nrow=thres_num,byrow=TRUE)
  sort_one_rank_m<-matrix(data=0,ncol=ncol(mirna_use),nrow=numm,byrow=TRUE)
  #sort_one_prob is the probability number that are ordered by prob
  sort_one_prob_m<-matrix(data=0,ncol=ncol(mirna_use),nrow=numm,byrow=TRUE)
  sigg<-rep(0,ncol(mirna_use))
  for(m in 1:ncol(mirna_use)){
    pro_matrix_temp<-(mirna_use[,m]%*%t(mrna_use[,m]))*wMRE/sumup[m]
    pro_matrix_borda <- pro_matrix_borda + rank(-pro_matrix_temp)
    #pro_matrix_sumup <- pro_matrix_sumup + pro_matrix_temp
    sort_one_rank[,m]<-order(pro_matrix_temp,decreasing=TRUE)[1:thres_num]
    sort_one_prob[,m]<-sort(pro_matrix_temp,decreasing=TRUE)[1:thres_num]
    sort_one_rank_m[,m]<-order(pro_matrix_temp,decreasing=TRUE)[1:numm]
    sort_one_prob_m[,m]<-sort(pro_matrix_temp,decreasing=TRUE)[1:numm]
    sigg[m]<-floor(length(pro_matrix_temp[pro_matrix_temp!=0])/100)
    sample_num[,m]<-order(pro_matrix_temp,decreasing=TRUE)[1:1000]
    thres_num1<-floor(length(which(pro_matrix_temp!=0))/100)
    tot<-tot+thres_num1
    temp1<-order(pro_matrix_temp,decreasing=TRUE)[1:thres_num1]
    uni<-union(temp1, uni)
    inter<-intersect(temp1, inter)
  }
  
  thres_num2<-floor(tot/ncol(mirna_use))
  sort_sumup_borda<-order(pro_matrix_borda,decreasing=FALSE)[1:thres_num2]
  
  colnames(sort_one_rank)<-ID
  colnames(sort_one_prob)<-ID
  colnames(sample_num)<-ID
  #sort_sumup_oneave saves the position of pairs in the sumup matrix that are ordered by prob
  sort_sumup_borda<-order(pro_matrix_borda,decreasing=FALSE)[1:thres_num]
  sort_sumup_borda_cor<-order(pro_matrix_borda,decreasing=FALSE)[1:thres_num2]
  #sort_sumup_oneave<-order(pro_matrix_sumup,decreasing=TRUE)[1:thres_num]
  return(list(sort_one_rank, sort_one_prob, sort_sumup_borda, sort_sumup_borda_cor, uni, inter, sample_num, sort_one_rank_m, sort_one_prob_m, sigg))
  #return(list(pro_matrix_sumup, sort_one_sumup, sort_sumup_oneave, prob, sort_one_prob, CPTS_full, APTS_full))
}
#this function is used to get the detailed result
detail_output_five<-function(sort_one_rank, sort_one_prob, ID, mirna_name, mrna_name){
  k0=5*ncol(sort_one_rank)
  m<-matrix(data=0,ncol=k0,nrow=5000,byrow=TRUE)
  temp<-c()
  for(j in 1:ncol(sort_one_rank)){
    for(i in 1:5000){
      k1<- 5*j-4
      k2<- 5*j-3
      k3<- 5*j-2
      k4<- 5*j-1
      k5<- 5*j
      temp1<-ceiling(sort_one_rank[i,j]/length(mirna_name))
      temp2<-sort_one_rank[i,j]-(ceiling(sort_one_rank[i,j]/length(mirna_name))-1)*length(mirna_name)
      m[i, k1]<-i
      m[i, k2]<-mrna_name[temp1,1]
      m[i, k3]<-mrna_name[temp1,2]
      m[i, k4]<-mirna_name[temp2]
      m[i, k5]<-sort_one_prob[i,j]
    }
    temp<-c(temp, ID[j], ID[j], ID[j], ID[j], ID[j])
  }
  colnames(m)<-temp
  return(m)
}

mrna<-as.matrix(read.table("mrna_list.txt",head=TRUE,sep="\t"))
mirna<-as.matrix(read.table("mirna_list.txt",head=TRUE,sep="\t"))
wMRE_matrix_cwcs1<-as.matrix(read.table("wMRE_all.txt",head=TRUE,sep="\t"))
wMRE_matrix_cwcs<-abs(wMRE_matrix_cwcs1)

rm(wMRE_matrix_cwcs1)
gc()

name_cancer<-as.matrix(read.table("TCGA_ESCC_1to1.txt",head=TRUE,sep="\t"))
mirna_cancer<-as.matrix(read.table("ESCA.miRseq_mature_RPM_log2.txt",head=FALSE,sep="\t"))
mrna_cancer<-as.matrix(read.table("ESCA.uncv2.mRNAseq_RSEM_normalized_log2.txt",head=FALSE,sep="\t"))

cutoff<-0.8
x<-name_func(name_cancer, mirna_cancer, mrna_cancer)
z1<-mirna_matrix(name_cancer, mirna_cancer, x[[1]], mirna)
z2<-mrna_matrix(name_cancer, mrna_cancer, x[[2]], mrna)
z3<-mirna_mrna_data_unselected(z1[[1]],z2[[1]],z1[[2]],z2[[2]],z2[[3]],z2[[4]])
z7<-mirna_mrna_loc(z3[[3]], z3[[4]], mirna, mrna, wMRE_matrix_cwcs, z3[[6]])
z9<-sum_wepro(z3[[1]],z3[[2]],z7[[3]])

thres_num<-5000
z16<-wepro_rankresult(z3[[1]], z3[[2]], z7[[3]], z9, name_cancer[,1], thres_num)
z65<-detail_output_five(z16[[1]], z16[[2]], name_cancer[,1], z3[[3]], z3[[5]])
write.table(z65,file="top 5000 detail information of miRTIGO.txt",quote=FALSE,sep="\t",row.names = FALSE)

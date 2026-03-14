setwd('F:/Author_Template__2_main/plos/data')
getwd()
library(psych)
library(igraph)

########Input
mul_edge_list <- read.table('F:/Author_Template__2_main/plos/data/edge_file.txt')
mul_gene_list <- read.table('F:/Author_Template__2_main/plos/data/point_file.txt')
tf_gene <- read.table('F:/Author_Template__2_main/plos/data/TF_trrust_rawdata.human.txt')
MIR <- read.table("F:/Author_Template__2_main/plos/data/mirna_MIR.txt")
PATHWAY <- read.csv("F:/Author_Template__2_main/plos/data/kegg_pathway.csv",F)


####brca
G_mut <- read.table('F:/Author_Template__2_main/plos/data/mut_brca_18846.txt',sep='\t',TRUE)
def_gene <- read.csv("brca_degs_05.csv")
mrna_exp <- read.table('F:/Author_Template__2_main/plos/data/Exp_brca.txt',sep='\t',TRUE)
dif_module <- read.csv("F:/Author_Template__2_main/plos/data/brca_nodewl60_dif_module10_1408_0001.csv")

####luad
G_mut <- read.table('F:/Author_Template__2_main/plos/data/mut_luad_18497.txt',sep='\t',TRUE)
def_gene <- read.csv("luad_degs_05.csv")
mrna_exp <- read.table('F:/Author_Template__2_main/plos/data/Exp_luad.txt',sep='\t',TRUE)
dif_module <- read.csv("F:/Author_Template__2_main/plos/data/luad_nodewl60_dif_module10_1010_0001.csv")

#prad

G_mut <- read.table('F:/Author_Template__2_main/plos/data/mut_prad_12415.txt',sep='\t',TRUE)
def_gene <- read.csv("prad_degs_025.csv")
mrna_exp <- read.table('F:/Author_Template__2_main/plos/data/Exp_prad.txt',sep='\t',TRUE)
dif_module <- read.csv("F:/Author_Template__2_main/plos/data/prad_nodewl60_dif_module10_899_0001.csv")


Construct_PPINET <- function(mul_edge_list,mul_gene_list,tf_gene){
  names(mul_gene_list) <- c('n','genes')  
  mul_genenames <- mul_gene_list$genes    
  AdjMatrix <- matrix(data=0,nrow=length(mul_genenames),ncol=length(mul_genenames)) 
  rownames(AdjMatrix) <- mul_genenames   
  colnames(AdjMatrix) <- mul_genenames
  for(i in 1:nrow(mul_edge_list))
  {
    tem_list <- mul_edge_list[i,]  
    x<-as.numeric(tem_list[1])   
    y<- as.numeric(tem_list[2])
    AdjMatrix[x,y] = 1
    AdjMatrix[y,x] = 1 
  }
  dele_row <- which(tf_gene[,3]=='Unknown')  
  tf_known_r <- nrow(tf_gene) - length(dele_row)
  tf_known <- matrix(0,tf_known_r,3)   
  j <- 1
  for (i in 1:nrow(tf_gene)) {
    if(tf_gene[i,3]=='Unknown')next()    
    tf_known[j,1] <- tf_gene[i,1]
    tf_known[j,2] <- tf_gene[i,2]
    tf_known[j,3] <- tf_gene[i,3]
    j <- j + 1
  }
  
  count <- 0     
  for (i in 1:nrow(tf_known)) {
    gene1 <- tf_known[i,1]
    gene2 <- tf_known[i,2]
    index1 <- which(mul_gene_list[,2]==gene1)
    index2 <- which(mul_gene_list[,2]==gene2)
    AdjMatrix[index1,index2] <- 1
    AdjMatrix[index2,index1] <- 1
    count <- count + 1
  }
  for (i in 1:nrow(AdjMatrix)) {
    AdjMatrix[i,i] <- 0
  }
  return(AdjMatrix)
}

detect_communities <- function(adj_matrix, method = "louvain") {
  # 输入验证
  stopifnot(
    is.matrix(adj_matrix),
    nrow(adj_matrix) == ncol(adj_matrix),
    all(adj_matrix %in% c(0,1))
  )
  
  # 创建图对象
  g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  
  # 选择社区检测算法
  communities <- switch(
    method,
    "louvain" = cluster_louvain(g),
    "walktrap" = cluster_walktrap(g),
    "infomap" = cluster_infomap(g),
    "labelprop" = cluster_label_prop(g),
    stop("Unknown method. Supported: louvain/walktrap/infomap/labelprop")
  )
  
  # 构建结果矩阵
  result <- cbind(
    Node = 1:nrow(adj_matrix),
    Community = communities$membership
  )
  
  return(result)
}

IDGDD <- function(AdjMatrix,G_mut,MIR,dif_module,beta){
  G_mut<-G_mut[,-1]
  gene_mut1 <- matrix(data=0,nrow=ncol(G_mut),ncol=2)
  gene_mut1[,1] <- colnames(G_mut)
  mul_genenames <- rownames(AdjMatrix)
  gene_first <- intersect(mul_genenames,gene_mut1[,1])
  randomwalk <- matrix(0,length(gene_first),length(gene_first))
  colnames(randomwalk) <- gene_first
  rownames(randomwalk) <- gene_first
  i1 <- nrow(randomwalk)-1
  for (i in 1:i1) {
    rindex1 <- which(mul_genenames==gene_first[i])
    rneib <- which(AdjMatrix[rindex1,]!=0)
    if(length(rneib)==0)next
    for (j in 1:length(rneib)) {
      rindex2 <- rneib[j]
      cindex2 <- which(gene_first==mul_genenames[rindex2])
      if(length(rindex2)==0)next
      randomwalk[i,cindex2] <- 1
      randomwalk[cindex2,i] <- 1
    }
  }
  tmirna <- unique(MIR[,2])
  gene_mi_matrix <- matrix(0,length(gene_first),length(tmirna))
  colnames(gene_mi_matrix) <- tmirna
  rownames(gene_mi_matrix) <- gene_first
  
  for (i in 1:nrow(gene_mi_matrix)) {
    rindex1 <- which(MIR[,4]==gene_first[i])
    if(length(rindex1)==0)next
    mir1 <- MIR[rindex1,2]
    mir1 <- unique(mir1)
    for (j in 1:ncol(gene_mi_matrix)) {
      rindex2 <- which(mir1==tmirna[j])
      if(length(rindex2)!=0){
        gene_mi_matrix[i,j] <-1
      }
      
    }
  }
  gg_mi_matrix <- matrix(0,length(gene_first),length(gene_first))
  colnames(gg_mi_matrix) <- gene_first
  rownames(gg_mi_matrix) <- gene_first
  i1 <- nrow(gg_mi_matrix)-1
  for (i in 1:i1) {
    rindex1 <- which(gene_mi_matrix[i,]!=0)
    if(length(rindex1)==0)next
    j1 <- i+1
    for (j in j1:nrow(gg_mi_matrix)) {
      #print(j)
      rindex2 <- which(gene_mi_matrix[j,]!=0)
      if(length(rindex2)==0)next
      in_mi <- intersect(rindex1,rindex2)
      uni_mi <- union(rindex1,rindex2)
      gg_mi_matrix[i,j] <- length(in_mi)/length(uni_mi)
      gg_mi_matrix[j,i] <- length(in_mi)/length(uni_mi)
    }
  }
  rweight <- randomwalk
  randomwalk <- rweight + gg_mi_matrix
  randomwalk <- randomwalk/2
  CCN <- randomwalk
  M_matrix_weight <- CCN
  gene_score <- matrix(0,nrow = length(gene_first),ncol = 6)
  rownames(gene_score)<- gene_first
  for (i in 1:nrow(CCN)) {
    gene_score[i,2] <-sum(as.numeric(CCN[i,]))
  }
  for(i in 1:nrow(CCN))
  {
    sum1<- as.numeric(gene_score[i,2])
    if(sum1==0)next()
    neib_list1 <- which(CCN[i,]!=0)
    for(j in 1:length(neib_list1))
    {
      gene_index1 <- neib_list1[j]
      sum2 <- as.numeric(gene_score[gene_index1,2])
      M_matrix_weight[i,gene_index1] <- abs(CCN[i,gene_index1])/sqrt(sum1*sum2)
    }
  }
  M_matrix1 <- matrix(data=0,nrow=nrow(CCN),ncol=1)
  for(i in 1:length(gene_first))
  {
    index_find <- which(colnames(G_mut)==gene_first[i])
    gene_score[i,3] <- sum(G_mut[,index_find])/nrow(G_mut)
  }
  for(i in 1:nrow(CCN))
  {
    
    M_matrix1[i,1] <- as.numeric(gene_score[i,3])
  }
  Exp_initial <- M_matrix1
  rownames(Exp_initial) <- gene_score[,1]
  
  count2=1
  Exp_tem <- (1-beta)*M_matrix_weight%*%Exp_initial + beta*Exp_initial
  Exp_final <- (1-beta)*M_matrix_weight%*%Exp_tem + beta*Exp_initial
  Exp_text <- abs(Exp_final - Exp_tem)
  result_index <- which((Exp_text)>10^-6)
  while(length(result_index)!=0){
    Exp_tem <- Exp_final
    Exp_final <- (1-beta)*M_matrix_weight%*%Exp_tem + beta*Exp_initial
    Exp_text <- abs(Exp_final - Exp_tem)
    result_index <- which((Exp_text)>10^-6)
    count2 <- count2 + 1
    
  }
  for (i in 1:nrow(gene_score)) {
    
    gene_score[i,4] <- Exp_final[i,1]
  }
  
  PATHWAY <- PATHWAY[,-1]
  PATHWAY <- PATHWAY[,-1]
  
  gene_pathway <- matrix(0, length(gene_first), nrow(PATHWAY))
  rownames(gene_pathway) <- gene_first
  colnames(gene_pathway) <- 1:nrow(PATHWAY)
  
  gene_pathway <- t(apply(PATHWAY, 1, function(pathway_genes) {
    as.integer(gene_first %in% pathway_genes)
  }))
  gene_pathway <- t(gene_pathway)  
  cs <- colSums(gene_pathway)
  re_col <- which(cs>1)
  gene_pathway <- gene_pathway[,re_col]
  rownames(gene_pathway) <- gene_first
  retain_gene <- which(rowSums(gene_pathway)==0)
  gene_pathway <- as.matrix(gene_pathway)
  
  part1 <- matrix(0,length(gene_first),ncol(gene_pathway))
  rownames(part1) <- gene_first
  colnames(part1) <- 1:ncol(gene_pathway)
  
  for (i in 1:length(gene_first)) {
    e1 <- which(gene_pathway[i,]!=0)
    if(length(e1)==0)next
    
    part1[i,e1] <- 1
    s1 <- sum(part1[i,])
    for (j in 1:ncol(part1)) {
      part1[i,j] <- part1[i,j]/s1
    }
  }
  
  part2 <- matrix(0,length(gene_first),ncol(gene_pathway))
  rownames(part2) <- gene_first
  colnames(part2) <- 1:ncol(gene_pathway)
  for (i in 1:ncol(part2)) {
    s1 <- sum(gene_pathway[,i])
    if(s1==0)next
    for (j in 1:nrow(part2)) {
      part2[j,i] <- gene_pathway[j,i]/s1
    }
  }
  
  tP <- part1%*%t(part2)
  reta_gene <- which(rowSums(tP)==0)
  
  #beta <- 0.4
  count2=1
  Exp_initial <- gene_score[,3]
  Exp_tem <- (1-beta)*tP%*%Exp_initial + beta*Exp_initial
  Exp_final <- (1-beta)*tP%*%Exp_tem + beta*Exp_initial
  Exp_text <- abs(Exp_final - Exp_tem)
  result_index <- which((Exp_text)>10^-6)
  while(length(result_index)!=0){
    Exp_tem <- Exp_final
    Exp_final <- (1-beta)*tP%*%Exp_tem + beta*Exp_initial
    Exp_text <- abs(Exp_final - Exp_tem)
    result_index <- which((Exp_text)>10^-6)
    count2 <- count2 + 1
    
  }
  for (i in 1:nrow(gene_score)) {
    gene_score[i,5] <- Exp_final[i,1]
  }
  
  for (i in 1:length(retain_gene)) {
    j <- retain_gene[i]
    gene_score[j,5] <- gene_score[j,4] 
  }
  for (i in 1:nrow(gene_score)) {
    
    gene_score[i,5] <-  (gene_score[i,4] + gene_score[i,5])/2
  }
  
  CCN <- randomwalk
  for (i in 1:nrow(CCN)) {
    del <- which(CCN[i,]!=0)
    CCN[i,del] <- 1
  }
  
  neighbors <- apply(CCN, 1, function(x) which(x == 1))
  
  M4 <- matrix(0, nrow = length(gene_first), ncol = length(gene_first))
  rownames(M4) <- gene_first
  colnames(M4) <- gene_first
  
  for (i in seq_len(nrow(CCN))) {
    gene1_neib <- neighbors[[i]]
    if (length(gene1_neib) == 0) next
    for (gene2 in gene1_neib) {
      common_neighbors <- intersect(neighbors[[i]], neighbors[[gene2]])
      ke <- length(common_neighbors)
      
      if (ke < 2) {
        M4[i, gene2] <- 0
        next
      }
      sub_matrix <- CCN[common_neighbors, common_neighbors]
      M4[i, gene2] <- sum(sub_matrix) / (ke^2)
    }
  }
  
  for (i in 1:nrow(M4)) {
    gene_score[i,6] <- sum(as.numeric(M4[i,]))
    
  }
  
  rownames(mrna_exp) <- mrna_exp[,1]
  mrna_exp <- mrna_exp[,-1]
  for (v1 in 1:ncol(mrna_exp)) {
    aa <- substr(colnames(mrna_exp)[v1],14,15)
    bb <- as.numeric(aa)
    if(bb<11)break
  }
  v2 <- v1 - 1
  v3 <- ncol(mrna_exp)
  
  gene_exp1 <- matrix(data=0,nrow=nrow(mrna_exp),ncol=2)
  rownames(gene_exp1) <- rownames(mrna_exp)
  gene_exp1[,1] <- 1:nrow(mrna_exp)
  for(i in 1:nrow(mrna_exp))
  {
    m1 <- mean(as.numeric(mrna_exp[i,1:v2]))
    m2 <- mean(as.numeric(mrna_exp[i,v1:v3]))
    var1 <- sd(as.numeric(mrna_exp[i,1:v2]))
    var2 <- sd(as.numeric(mrna_exp[i,v1:v3]))
    var_sum <- var1^2+var2^2
    sum1 <- (m1-m2)^2/(4*var_sum) 
    sum2 <-  0.5*log(var_sum/(2*var1*var2))
    gene_exp1[i,2] <- sum1+sum2
  }
  dele_index <- which(gene_exp1[,2]=='Inf')
  gene_exp1[dele_index,2] <- 0
  
  gene_exp1 <- gene_exp1[order(gene_exp1[,2],decreasing = T),]
  
  dif_gene <- def_gene[,1]
  gene_second <- intersect(gene_first,dif_gene)
  dif_gene <- gene_second
  def_exp <- gene_exp1[match(dif_gene,rownames(gene_exp1)),]
  
  exp_score <- matrix(0,length(gene_first),3)
  rownames(exp_score) <- gene_first
  colnames(exp_score) <- c("mrna","mirna","score")
  for (i in 1:nrow(exp_score)) {
    iscore <- 0
    
    r1 <- which(rownames(CCN)==gene_first[i])
    g_neib <- which(CCN[r1,]!=0)
    if(length(g_neib)==0){
      exp_score[i,1] <- iscore
      next
    }
    for (j in 1:length(g_neib)) {
      g2 <- g_neib[j]
      rindex2 <- which(rownames(def_exp)==gene_first[g2])
      
      if(length(rindex2)!=0){
        iscore <- iscore + as.numeric(def_exp[rindex2,2])
      }
    }
    exp_score[i,1] <- iscore
  }
  
  exp_score[,3] <- exp_score[,1]  
  gene_score[,1] <- 1:nrow(gene_score)
  gene_score3 <- cbind(gene_score,exp_score)
  
  max_v <- 9000
  
  rownames(dif_module) <- dif_module[,1]
  dif_module <- dif_module[,-1]
  dif_module <- as.matrix(dif_module)
  for (i in 1:nrow(dif_module)) {
    dif_module[i,i] <- 0
    del_col <- which(dif_module[i,]<max_v)
    dif_module[i,del_col] <- 0
    del_coln <- which(dif_module[i,]>=max_v)
    dif_module[i,del_coln] <- 1
    
  }
  
  comm_module <- detect_communities(dif_module, method = "louvain")
  comm_module[,1] <- dif_gene
  com_freq <- comm_module[,2]
  freq_table <- table(com_freq)
  
  
  elements_to_keep <- names(freq_table)[freq_table >= 3]
  
  
  filtered_data <- unique(com_freq[com_freq %in% as.numeric(elements_to_keep)])
  
  num_c <- length(filtered_data)
  exp_score <- matrix(0,length(gene_first),3)
  rownames(exp_score) <- gene_first
  colnames(exp_score) <- c("mrna","mirna","score")
  for (i in 1:nrow(exp_score)) {
    iscore <- 0
    mv <- 0
    ml <- 0
    cv <- 0
    same_ml <- 0
    rindex1 <- which(rownames(def_exp)==gene_first[i])
    if(length(rindex1)!=0){
      iscore <- as.numeric(def_exp[rindex1,2])
    }
    r1 <- which(rownames(CCN)==gene_first[i])
    g_neib <- which(CCN[r1,]!=0)
    if(length(g_neib)==0){
      exp_score[i,1] <- iscore
      next
    }
    same_com <- matrix(0,1,num_c)
    for (k in 1:num_c) {
      k_n <- filtered_data[k]
      list1 <- comm_module[which(comm_module[,2]==k_n),1]
      
      list2 <- intersect(gene_first[g_neib],list1)
      if(length(list2)==mv && length(list2)>0){
        same_ml <- same_ml + 1
        same_com[1,same_ml] <- k_n
      }
      if(length(list2)>mv){
        mv <- length(list2)
        ml <- k_n
        same_com[1,] <- 0
        same_com[1,1] <- k_n
        same_ml <- 1
      }
      
    }
    houxuan_com <- setdiff(unique(same_com),"0")
    
    if(length(houxuan_com)!=0){
      same_com <- matrix(0,1,length(houxuan_com))
      for (cj in 1:length(houxuan_com)) {
        cv <- 0
        
        list1 <- comm_module[which(comm_module[,2]==houxuan_com[cj]),1]
        list2 <- intersect(gene_first[g_neib],list1)
        
        for (j in 1:length(list2)) {
          g2 <- list2[j]
          rindex2 <- which(rownames(def_exp)== g2)
          if(length(rindex2)!=0){
            cv <- cv + as.numeric(def_exp[rindex2,2])
          }
        }
        same_com[1,cj] <- cv
      }
      
    }
    exp_score[i,1] <- max(same_com[1,])
  }
  
  
  exp_score <- as.matrix(exp_score)
  exp_score[,3] <- exp_score[,1]  
  gene_score[,1] <- 1:nrow(gene_score)
  gene_score3 <- cbind(gene_score,exp_score)
  
  nor_genes <- matrix(0,length(gene_first),4)
  rownames(nor_genes) <- gene_first
  colnames(nor_genes) <- c("ncol","nnet","nexp","tscore")
  m_col <- max(gene_score3[,5])
  mi_col <- min(gene_score3[,5])
  m_net <- max(gene_score3[,6])
  mi_net <- min(gene_score3[,6])
  m_exp <- max(gene_score3[,7])
  mi_exp <- min(gene_score3[,7])
  
  for (i in 1:nrow(nor_genes)) {
    nor_genes[i,1] <- (gene_score3[i,5] - mi_col)/m_col
    va1 <- nor_genes[i,1]
    nor_genes[i,2] <- (gene_score3[i,6] - mi_net)/m_net
    va2 <- nor_genes[i,2]
    nor_genes[i,3] <- (gene_score3[i,7] - mi_exp)/m_exp
    va3 <- nor_genes[i,3]
    if(va1*va2*va3!=0){
      vl <- c(va1,va2,va3)
      nor_genes[i,4] <- harmonic.mean(vl)
      
    }
    
  }
  colnames(gene_score3) <- c("n","s1","s2","col1","col2","net","mrna","mirna","score")
  brca_score <- cbind(gene_score3,nor_genes)
  brca_score <- brca_score[order(brca_score[,13],decreasing = 'T'),]
  
  
  return(brca_score)
  
}


AdjMatrix <- Construct_PPINET(mul_edge_list,mul_gene_list,tf_gene)

brca_list <- IDGDD(AdjMatrix,G_mut,MIR,dif_module,0.4)
















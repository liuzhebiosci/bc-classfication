# install the required packages
# install.packages("ROCR")
# install.packages("e1071")

options(stringsAsFactors = FALSE)
options(help_type="html")
rm(list=ls())

library(e1071)
library(ROCR)
library(gplots)

setwd("/zheliu/DevGSEA/")

# function to read .gct .gmt and expression files
source("/zheliu/DevGSEA/code/core code/GSEA_Gct2Frame.R")
source("/zheliu/DevGSEA/code/core code/trim.R")
source("/zheliu/DevGSEA/code/core code/get_perf.R")

# read experiment data with Entrez gene IDs as rownames
raw.data<- GSEA.Gct2Frame("/zheliu/DevGSEA/input data/New data/new_exprs_entrez_id.gct")

# read clinical information
descrp.data<- read.delim("/zheliu/DevGSEA/input data/New data/descrp.txt", sep="\t", fill=T)

#read KEGG pathway data
kegg.gs=readLines("/zheliu/DevGSEA/input data/New data/kegg_pathway_mapped.gmt")
kegg.gs.anno=do.call(rbind, lapply(kegg.gs, function(x) trim(unlist(strsplit(x, "\t")))[c(1,2)]))
kegg.gs.gene=do.call(rbind, lapply(kegg.gs, function(x) trim(unlist(strsplit(x, "\t")))[-c(1,2)]))

#kegg.gs.gene[1:10, 1:10]
#kegg.gs.anno[1:10, ]

#read TF gene set data
TF.gs=readLines("/zheliu/DevGSEA/input data/New data/TF_gene_pathway.csv")
TF.gs.anno=do.call(rbind, lapply(TF.gs, function(x) trim(unlist(strsplit(x, "\t")))[c(1,2,3,4)]))
TF.gs.gene=do.call(rbind, lapply(TF.gs, function(x) trim(unlist(strsplit(x, "\t")))[-c(1,2,3,4)]))

#TF.gs.anno[1:10, ]
#TF.gs.gene[1:10, 1:10]

#gene level z-score normalization of the data
exprs.matrix<- round(t(apply(raw.data, 1, function(x) (x-mean(x))/sd(x))),4)
exprs.class<- descrp.data[, 3]

#construct pathway-TF gene set mapping relationship
pathway<- list()
pw.rank.entire<- list()

#TF.gs.anno[1:10,]
#kegg.gs.anno[1:10,]

#match pathway to transcription factor
TF.pathway.info=lapply(TF.gs.anno[,4], function(x) gsub(unlist(strsplit(x, "\\|")), pattern="\\.\\..*$", replacement="", ignore.case = FALSE))
pathway.list=unique(unlist(TF.pathway.info))
names(TF.pathway.info)=TF.gs.anno[,1]

#match tf to pathway
pathway.TF.info=lapply(lapply(pathway.list, grep, TF.pathway.info), function (x) names(TF.pathway.info[x]))
names(pathway.TF.info)=pathway.list

#match pathway to tf index
pathway.tf.map<- lapply(pathway.TF.info, match, TF.gs.anno[,1])
pathway.tf.map[1]

n.pathway=length(pathway.tf.map)
n.sample=ncol(exprs.matrix)

# load the data
# save.image("/zheliu/DevGSEA/result2019/imbalanced 50/input.Rdata")
# load("/zheliu/DevGSEA/result2019/imbalanced 50/input.Rdata")

n.sample<- length(exprs.class)
n.pathway<- length(pathway.tf.map)

sub.ratio<- 0.25

# create the results directory
dir.create(file.path(paste("/zheliu/DevGSEA/result2019/test_blind_imba_50_sub_", sub.ratio, sep="")), showWarnings = FALSE)
res.dir<- paste("/zheliu/DevGSEA/result2019/test_blind_imba_50_sub_", sub.ratio, sep="")

ptm<- proc.time()

test.perf.file=""
test.perf.cv.file=""
val.perf.cv.file=""

pos.ind<- which(exprs.class==1)
neg.ind<- which(exprs.class==-1)

n.pos.ind<- length(pos.ind)
n.neg.ind<- length(neg.ind)

# number of iteratioin
n.iter<- 50
n.cv<- 1

perf.file<- ""

for (k in 1:n.pathway){
  print(c("k", k))
  set.seed(1)
  
  # subsample ratio is 0.3 and the readjustment weight is 3.3
  adj.w<- round(1/sub.ratio, 3)
  
  # create a pathway folder
  dir.create(file.path(res.dir, kegg.gs.anno[k, 2]), showWarnings = FALSE)
  setwd(file.path(res.dir, kegg.gs.anno[k, 2]))
  
  pw.gs.gene.list<- list()
  pw.gs.anno.list<- list()
  
  # compute the length different between TF and KEGG
  len.dif<- length(TF.gs.gene[k, ])-length(kegg.gs.gene[k, ])
  
  # construct a list to store the pathway and associated TF gene sets
  pw.gs.gene.list[[1]]<- c(kegg.gs.gene[k, ], rep("", len.dif))
  pw.gs.anno.list[[1]]<- kegg.gs.anno[k, ]
  
  n.gs.pathway<- length(pathway.tf.map[[k]])+1
  
  for (j in 2:n.gs.pathway){
    pw.gs.gene.list[[j]] <- TF.gs.gene[pathway.tf.map[[k]][j-1], ]
    pw.gs.anno.list[[j]] <- TF.gs.anno[pathway.tf.map[[k]][j-1], ]
  }
  
  pw.gs.gene<- do.call(rbind, pw.gs.gene.list)
  pw.gs.anno<- do.call(rbind, pw.gs.anno.list)

  tryCatch({
    # testing dataset performance matrix
    val.perf.mat<- matrix(data=NA, nrow=1, ncol=6)
    test.perf.mat<- matrix(data=NA, nrow=1, ncol=6)
    
    # testing set classification result list
    test.class<- c()
    train.class<- c()
    val.class<- c()
    
    # testing set classification result list
    val.class.list<- list()
    test.class.list<- list()
    
    # training and testing dataset classification performance
    train.step.pred<- list()
    val.step.pred<- list()
    test.step.pred<- list()
    
    # test set step prediction results
    train.step.res<- list()
    val.setp.res<- list()
    test.step.res<- list()
    
    #training and testing set index 
    train.ind<- c()
    val.ind<- c()
    test.ind<- c()
    
    # parameter to control whether all of the data are from one class
    max.iter<- 0
    
    val.pos.ind<- seq(1, n.pos.ind)[floor(seq(1, n.pos.ind, 10))]
    val.neg.ind<- seq(1, n.neg.ind)[floor(seq(1, n.neg.ind, 10))]
    
    test.pos.ind<- seq(1, n.pos.ind)[floor(seq(2, n.pos.ind, 10))]
    test.neg.ind<- seq(1, n.neg.ind)[floor(seq(2, n.neg.ind, 10))]

    test.ind<- c(pos.ind[test.pos.ind], neg.ind[test.neg.ind])
    val.ind<- c(pos.ind[val.pos.ind], neg.ind[val.neg.ind])
    train.ind<- c(pos.ind[-c(test.pos.ind, val.pos.ind)], neg.ind[-c(test.neg.ind, val.neg.ind)])
    
    # subtract the respective experiments data and class labels
    train.exprs<- exprs.matrix[, train.ind]
    val.exprs<- exprs.matrix[, val.ind]
    test.exprs<- exprs.matrix[, test.ind]
    
    train.class<- exprs.class[train.ind]
    val.class<- exprs.class[val.ind]
    test.class<- exprs.class[test.ind]
    
    # put the class in the exprs.class in the respective cv
    test.class.list<- test.class
    val.class.list<- val.class
    
    cost<- 1/table(as.factor(train.class))
    ratio<- cost[2]/cost[1]
    
    alpha.list<- rep(NA, n.iter)
    err.list<- rep(NA, n.iter)
    class.sign<- rep(NA, n.iter)
    
    train.step.pred.mat<- matrix(data=NA, nrow=n.iter, ncol=6)
    val.step.pred.mat<- matrix(data=NA, nrow=n.iter, ncol=6)
    test.step.pred.mat<- matrix(data=NA, nrow=n.iter, ncol=6)
    
    model.list<- list()
    train.ratio.list<- list()
    
    n.train<- length(train.ind)
    n.test<- length(test.ind)
    
    w<- rep(1/n.train, n.train)
    
    misclass<- NA
    alpha<- NA
    
    for(t in 1:n.iter) {
      # subsampling the data with replacement
      n.sub.train<- round(sub.ratio*n.train)
      sample.ind<- sample(1:n.train, n.sub.train, replace=T, prob=w)
      w.sample<- w[sample.ind]
      sub.train.class<- train.class[sample.ind]
      
      # To avoid all of the samples are from one class, test whether all of the classes are the same.
      # If this is the case, then break the loop.
      if(abs(sum(sub.train.class))==n.sub.train){
        max.iter <- t-1
        print("all one class")
        break
      }
      
      # subtract the gene sets by repeatly indexing the TF gene sets.
      gs.ind<- na.omit(unlist(lapply(pw.gs.gene[rep_len(seq(1, n.gs.pathway), n.iter)[t], ], match, rownames(exprs.matrix))))
      gs.train.exprs<- train.exprs[gs.ind, ]
      
      # subtract the experiment data corresponding to the respective gene set 
      sub.gs.train.exprs<- train.exprs[gs.ind, sample.ind]
      
      #train model from dataset 1 and make prediction on the validation dataset
      # sub.train.weight<- 1/table(as.factor(sub.train.class))
      train.weight<- 1/table(as.factor(train.class))
      
      # a small value is added to circumvent duplicated data error
      dup.sample<- which(duplicated(colnames(sub.gs.train.exprs)))
      sub.gs.train.exprs[, dup.sample]<- sub.gs.train.exprs[, dup.sample]+seq(length(dup.sample))*1e-10
      colnames(sub.gs.train.exprs)<- seq(ncol(sub.gs.train.exprs))
      
      # svm parameters
      svm.cost=1
      svm.class.weigth=train.weight
      
      # train the model on the sub-sampled samples
      sub.train.model<- svm(t(sub.gs.train.exprs), as.factor(sub.train.class), kernel="linear", cost=svm.cost, class.weights=svm.class.weigth)
      
      # test the model on all samples 
      train.pred.res<- as.numeric(as.character(predict(sub.train.model, t(gs.train.exprs))))
      
      model.list[[t]]<- sub.train.model
      
      #boosting algorithm
      misclass<- abs(train.class-as.numeric(as.character(train.pred.res)))/2
      
      err<- t(w) %*% misclass
      err.list[[t]]<- err
      
      # assign class.sign as -1 when error is greater than 0.5, as +1 when error is less than 0.5
      if(err>0.5){
        class.sign[t]<- -1
        err<- 1-err
      }else{
        class.sign[t]<- 1
      }
      
      alpha<- 0.5*log((1-err)/(err+1e-15))
      alpha.list[t]<- alpha
      
      #weight the imbalaced class samples by cost adj.w
      w<- w*exp(-train.pred.res*train.class*alpha)*exp(train.class*adj.w*(1/n.iter)*log(sqrt(ratio), exp(1)))
      w<- w/sum(w)
      
      train.ratio.list[t]<- as.numeric(table(sub.train.class)[1]/length(sub.train.class))
      
      gs.val.exprs<- val.exprs[gs.ind, ]
      gs.test.exprs<- test.exprs[gs.ind, ]
      
      #calculate the training error when classifier is not null
      if(length(sub.train.model)!=0){
        train.step.pred[[t]]<- as.numeric(as.character(predict(sub.train.model, t(gs.train.exprs))))*class.sign[t]
        val.step.pred[[t]]<- as.numeric(as.character(predict(sub.train.model, t(gs.val.exprs))))*class.sign[t]
        test.step.pred[[t]]<- as.numeric(as.character(predict(sub.train.model, t(gs.test.exprs))))*class.sign[t]
      }
    }
    
    # use parameter max.iter to avoid all positive/negative samples
    if (max.iter!=0){
      max.iter <- min(n.iter, max.iter)
    }else{
      max.iter <- t
    }
    
    #calculate the training error and testing error
    train.step.pred.mat<- matrix(data = NA, nrow = max.iter, ncol = 6)
    colnames(train.step.pred.mat)<- c("auc", "accuracy", "fpr", "fnr", "fpn", "fnn")
    
    val.step.pred.mat<- matrix(data = NA, nrow = max.iter, ncol = 6)
    colnames(val.step.pred.mat)<- c("auc", "accuracy", "fpr", "fnr", "fpn", "fnn")
    
    test.step.pred.mat<- matrix(data = NA, nrow = max.iter, ncol = 6)
    colnames(test.step.pred.mat)<- c("auc", "accuracy", "fpr", "fnr", "fpn", "fnn")

    for (i in 1:max.iter) {
      train.step.res<- sign(alpha.list[1:i] %*% t(do.call(cbind, train.step.pred[1:i])))
      train.step.pred.mat[i, ]<- unlist(get.perf(as.numeric(as.character(train.step.res)), train.class))
      
      val.step.res<- sign(alpha.list[1:i] %*% t(do.call(cbind, val.step.pred[1:i])))
      val.step.pred.mat[i, ]<- unlist(get.perf(as.numeric(as.character(val.step.res)), val.class))
      
      test.step.res<- sign(alpha.list[1:i] %*% t(do.call(cbind, test.step.pred[1:i])))
      test.step.pred.mat[i, ]<- unlist(get.perf(as.numeric(as.character(test.step.res)), test.class))
    }
    
    test.perf.mat <- unlist(get.perf(as.numeric(as.character(test.step.res)), test.class))
    names(test.perf.mat)<- c("AUC", "accuracy", "fpr", "fnr", "fpn", "fnn")
    
    val.perf.mat <- unlist(get.perf(as.numeric(as.character(val.step.res)), val.class))
    names(val.perf.mat)<- c("AUC", "accuracy", "fpr", "fnr", "fpn", "fnn")
    ############################################################################################
    #plot the stepwise parameters at each iteration round, namely
    # error: err.list
    # alpha: alpha.list
    # training auc: train.step.pred.mat
    # testing auc: test.step.pred.mat
    # sampling ratio: train.ratio.list
    ############################################################################################
    pdf(file=paste("error.pdf", sep=""))
    plot(seq(1, max.iter), lty=1, err.list[1:max.iter], xlab="# of iteration", ylab="error", type="l")
    dev.off()
    
    pdf(file=paste("alpha.pdf",sep=""))
    plot(seq(1, max.iter), lty=1, alpha.list[1:max.iter], xlab="# of iteration", ylab="alpha", type="l")
    dev.off()
    
    pdf(file=paste("AUC.pdf",sep=""))
    plot(seq(1, n.iter), train.step.pred.mat[, 1], lty=1,  xlab="# of iteration", ylab="training/testing set AUC", ylim = c(0,1), type="l")
    lines(seq(1, n.iter), val.step.pred.mat[, 1], lty=2)
    lines(seq(1, n.iter), test.step.pred.mat[, 1], lty=3)
    legend("topleft", legend=c("training", "validation", "testing"), lty=1:3, cex=0.8)
    dev.off()
    
    pdf(file=paste("neg sample ratio.pdf",sep=""))
    plot(seq(1, max.iter), train.ratio.list[1:max.iter], xlab="# of iteration", ylab="neg ratio")
    dev.off()
    
    #write parameters to the files
    err.mat<- matrix(round(err.list, 3), ncol=1, nrow=length(err.list), byrow = TRUE, dimnames=list(c(seq(1, length(err.list))), c("error")))
    alpha.mat<- matrix(round(alpha.list, 3), ncol=1, nrow=length(alpha.list), byrow = TRUE, dimnames=list(c(seq(1, length(alpha.list))), c("alpha")))
    
    # write out the training and testing results
    write.table(t(c("iteration", "auc", "accuracy", "fpr", "fnr", "fpn", "fnn", "auc", "accuracy", "fpr", "fnr", "fpn", "fnn", "auc", "accuracy", "fpr", "fnr", "fpn", "fnn", "error", "alpha")), "results_cv.xls", append=T, quote=F, col.names=F, row.names=F, sep="\t")
    write.table(cbind(train.step.pred.mat, val.step.pred.mat, test.step.pred.mat, err.mat, alpha.mat), "results_cv.xls", append=T, quote=F, sep="\t", col.names=F)
    
    print(proc.time()-ptm)
  
    # compute the performance of validation and testing data split
    val.class.perf<- t(get.perf(as.numeric(unlist(val.step.res)), unlist(val.class)))
    test.class.perf<- t(get.perf(as.numeric(unlist(test.step.res)), unlist(test.class)))

    if(perf.file==""){
      perf.file<- paste(res.dir, "/perf.xls", sep="")
      write.table("\tAUC\taccuracy\tfpr\tfnr\tfpn\tfnn\tAUC\taccuracy\tfpr\tfnr\tfpn\tfnn", perf.file, quote=F, row.names=F, col.names=F)
    }
    write.table(cbind(kegg.gs.anno[k, 2], val.class.perf, test.class.perf), perf.file, append=T, quote=F, sep="\t", row.names=F, col.names=F)
    }, error = function(err){
    test.perf.mat="error here"
  })
}

par.mat<- cbind(c("subsampling ratio", "number of iterations", "svm prior", "svm cost"), c(ratio, n.iter, "train.weight", svm.cost))
write.table(par.mat, file=paste(res.dir, "/parameter.xls", sep=""), row.names=FALSE, col.names=F, sep="\t")

save.image(paste("/zheliu/DevGSEA/result2019/test_blind_imba_50_sub_", sub.ratio, "/test_set_blind.Rdata", sep=""))

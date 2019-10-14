options(stringsAsFactors = FALSE) 
options(help_type="html")
rm(list=ls())
sptm <- proc.time()

library(e1071)
library(ROCR)
library(gplots)

setwd("/C/Zhe Liu/DevGSEA/")

file.name<- c("Read Input File.R")
lapply(paste(getwd(), "code/core code", file.name, sep="/"), source)
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

#read in experiment data
#raw.data<- GSEA.Gct2Frame("./input data/New data/before_norm_exprs_entrez_id.gct")
raw.data<- GSEA.Gct2Frame("./input data/New data/new_exprs_entrez_id.gct")

#read in KEGG pathway data
kegg.gs=readLines("./input data/New data/kegg_pathway_mapped.gmt")
kegg.gs.gene=do.call(rbind, lapply(kegg.gs, function(x) trim(unlist(strsplit(x, "\t")))[-c(1,2)]))
kegg.gs.anno=do.call(rbind, lapply(kegg.gs, function(x) trim(unlist(strsplit(x, "\t")))[c(1,2)]))

#kegg.gs.gene[1:10, 1:10]
#kegg.gs.anno[1:10, ]

#read in TF gene set data
TF.gs=readLines("./input data/New data/TF_gene_pathway.csv")
TF.gs.anno=do.call(rbind, lapply(TF.gs, function(x) trim(unlist(strsplit(x, "\t")))[c(1,2,3,4)]))
TF.gs.gene=do.call(rbind, lapply(TF.gs, function(x) trim(unlist(strsplit(x, "\t")))[-c(1,2,3,4)]))

#TF.gs.anno[1:10, ]
#TF.gs.gene[1:10, 1:10]

#read in decription data
descrp.data<- read.table("./input data/New data/descrp.txt", sep="\t", head=T)

#Normalization of the data
entire.exprs <- t(apply(raw.data, 1, function(x) (x-mean(x))/sd(x)))
er.entire.class <- descrp.data[, 3]

#construct pathway-TF gene set mapping relationship
pathway<- list()
pw.rank.entire<- list()

#TF.gs.anno[1:10,]
#kegg.gs.anno[1:10,]

#match pathway to tf
TF.pathway.info=lapply(TF.gs.anno[,4], function(x) gsub(unlist(strsplit(x, "\\|")), pattern="\\.\\..*$", replacement="", ignore.case = FALSE))
pathway.list=unique(unlist(TF.pathway.info))
names(TF.pathway.info)=TF.gs.anno[,1]

#match tf to pathway
pathway.TF.info=lapply(lapply(pathway.list, grep, TF.pathway.info), function (x) names(TF.pathway.info[x]))
names(pathway.TF.info)=pathway.list

#match pathway to tf index
pathway.tf.map <- lapply(pathway.TF.info, match, TF.gs.anno[,1])
pathway.tf.map[1]

n.pathway=length(pathway.tf.map)
n.sample=ncol(entire.exprs)

#map KEGG gene set to exprs data
#index.exprs=function(gene.index){
#	entire.exprs[na.omit(unlist(lapply(gene.index[gene.index!=""], match, rownames(entire.exprs)))), ]
#}

###########################################################
#performance results
###########################################################
get.perf <- function(pred.res, class.labels){
	roc.pred=prediction(pred.res, class.labels)
	roc.pref=performance(roc.pred, "auc")
	auc=round(as.numeric(roc.pref@y.values), 3)
	
	#fp, tp, fn, tn of the prediction performance
	fp=roc.pred@fp[[1]][2]
	tp=roc.pred@tp[[1]][2]
	
	fn=roc.pred@fn[[1]][2]
	tn=roc.pred@tn[[1]][2]
	
	#classification accuracy 
	acc=round(((tp+tn)/sum(tp+tn+fp+fn)), 3)
	
	fpr=round((fp/(fp+tn)), 3)
	fnr=round((fn/(tp+fn)), 3)
	res <- list(auc, acc, fpr, fnr, fp, fn)
	names(res) <- c("AUC", "accuracy", "fpr", "fnr", "fpn", "fnn")
	return(res)
}

#save.image("/C/Zhe Liu/DevGSEA/result/NKI/test set blind/exprs.Rdata.Rdata")
#load(file="DevGSEA/result/NKI/test set blind/exprs.Rdata.Rdata")

setwd("/C/Zhe Liu")
description.data <- read.delim("./DevGSEA/dataset/expr data/raw data/descrp.txt", sep="\t", fill=T)

ptm <- proc.time()

n.cv <- 5

n.pathway=length(pathway.tf.map)
local.dir <- "/C/Zhe Liu/DevGSEA"
n.iter <- 50
para1 <- round(1/0.3, 3)

pos.ind <- which(er.entire.class==1)
neg.ind <- which(er.entire.class==-1)

pos.len <- length(pos.ind)
neg.len <- length(neg.ind)

#partition the data into training, validation and testing parts
train.pos.ind <- 1:round(1/3*length(pos.ind))
val.pos.ind <- round(1/3*length(pos.ind)):round(2/3*length(pos.ind))
test.pos.ind <- round(2/3*length(pos.ind)):length(pos.ind)

train.neg.ind <- 1:round(1/3*length(neg.ind))
val.neg.ind <- round(1/3*length(neg.ind)):round(2/3*length(neg.ind))
test.neg.ind <- round(2/3*length(neg.ind)):length(neg.ind)

train.ind <- c(train.pos.ind, train.neg.ind)
val.ind <- c(val.pos.ind, val.neg.ind)
test.ind <- c(test.pos.ind, test.neg.ind)

n.train <- length(train.ind)
train.exprs <- entire.exprs[, train.ind]
train.class <- er.entire.class[train.ind]
f.train <- as.factor(description.data[,9])[train.ind]

n.val <- length(val.ind)
val.exprs <- entire.exprs[, val.ind]
val.class <- er.entire.class[val.ind]
f.val <- as.factor(description.data[,9])[val.ind]

n.test <- length(test.ind)
test.exprs <- entire.exprs[, test.ind]
test.class <- er.entire.class[test.ind]
f.test <- as.factor(description.data[,9])[test.ind]

#put the classification performance in a list
combine.val.cl.pred <- list()
f1 <- as.factor(description.data[,9])

dir.create(file.path(local.dir, "result/NKI/test set blind"), showWarnings = FALSE)

for (k in 1:n.pathway){
#	k=1n.pathway
	combine.val.cl.pred[[k]] <- list()

	dir.create(file.path(local.dir, "result/NKI/test set blind", kegg.gs.anno[k, 2]), showWarnings = FALSE)
	setwd(file.path(local.dir, "/result/NKI/test set blind", kegg.gs.anno[k, 2]))
	
	set.seed(1)
	print(para1)
	val.auc.arr <- list()
	
	n.gs.pathway <- length(pathway.tf.map[[k]])+1
	
	pw.gs.gene.list <- list()
	pw.gs.anno.list <- list()
	
	len.dif <- length(TF.gs.gene[k, ])-length(kegg.gs.gene[k, ])
	
	pw.gs.gene.list[[1]] <- c(kegg.gs.gene[k, ], rep("", len.dif))
	pw.gs.anno.list[[1]] <- kegg.gs.anno[k, ]
	
	for (j in 2:n.gs.pathway){
		pw.gs.gene.list[[j]] <- TF.gs.gene[pathway.tf.map[[k]][j-1], ]
		pw.gs.anno.list[[j]] <- TF.gs.anno[pathway.tf.map[[k]][j-1], ]
	}
	
	pw.gs.gene <- do.call(rbind, pw.gs.gene.list)
	pw.gs.anno <- do.call(rbind, pw.gs.anno.list)
	
	print(proc.time()-ptm)
	
	tryCatch({
		print(c("k", k))
		
		val.perf.df <- matrix(data = NA, nrow = n.cv, ncol = 6)
		test.perf.df <- matrix(data = NA, nrow = n.cv, ncol = 6)

		max.iter <- 0

		cost <- 1/table(as.factor(train.class))
		ratio <- cost[2]/cost[1]
		
		alpha <- rep(NA,n.iter)
		errarr <- rep(NA,n.iter)
		
		cl.label <- rep(NA,n.iter)
		train.auc.arr <- matrix(data = NA, nrow = n.iter, ncol = 6)
		
		classifiers <- list()
		cl.pred <- list()
		train.ratio.arr <- list()
		
		n.train <- length(train.ind)
		n.val <- length(val.ind)
		
		w <- rep(1/n.train, n.train)
		
		#Initialize the subsampling of the training dataset
		f1.train <- f1[train.ind]
		sub.train.ind <- c()
		
		for (i in 1:nlevels(f1.train)){
			f1.length <- length(which(f1.train==levels(f1.train)[i]))
			sub.train.ind <- c(sub.train.ind, which(f1.train==levels(f1.train)[i])[1:round(0.3*f1.length)])
		}
		
		n.sub.train <- length(sub.train.ind)
		
		train.cl.pred <- list()
		val.cl.pred <- list()
		test.cl.pred <- list()
		
		misclass <- NA
		myalpha <- NA
		
		for(i in 1:n.iter) {
			#i=1
			sample.ind <- sample(1:n.train, n.sub.train, replace=F, prob=w)
			w.sample <- w[sample.ind]
			
			sub.train.class <- train.class[sample.ind]
			
			if(nlevels(as.factor(sub.train.class))==1){
				max.iter <- i-1
				break
			}
			
			gs.ind <- na.omit(unlist(lapply(pw.gs.gene[rep_len(seq(1,n.gs.pathway), n.iter)[i], ], match, rownames(entire.exprs))))
			
			gs.train.exprs <- train.exprs[gs.ind, ]
			sub.gs.train.exprs <- train.exprs[gs.ind, sample.ind]
			
			#train model from dataset 1 and make prediction on the validation dataset
			sub.train.weight <- 1/table(as.factor(sub.train.class))
			train.weight <- 1/table(as.factor(train.class))
			
			sub.train.model <- svm(t(sub.gs.train.exprs), as.factor(sub.train.class), kernel="linear", cost=1, class.weights=train.weight)#, class.weights=c("-1"=per.train.neg,"1"=per.train.pos))
			train.pred.res <- predict(sub.train.model, t(gs.train.exprs))
			classifiers[[i]] <- sub.train.model
			
			train.clas.res <- get.perf(as.numeric(as.character(train.pred.res)), train.class)
			train.clas.df <- do.call(rbind, train.clas.res)
			
			#boosting algorithm
			misclass <- abs(train.class-as.numeric(as.character(train.pred.res)))/2
			
			err <- (t(w) %*% misclass)
			errarr[[i]] <- err
			
			if(err>0.5){
				cl.label[i] <- -1
				err <- 1-err
			}else{
				cl.label[i] <- 1
			}
			
			myalpha <- 0.5*log((1-err)/(err+1e-15))
			alpha[i] <- myalpha
			
			#weight the imbalaced class samples by a cost
			w <- w*exp(myalpha * misclass)*exp(train.class*para1*(1/n.iter)*log(sqrt(ratio), exp(1)))
			w <- w/sum(w)
			
			train.ratio.arr[i] <- as.numeric(table(sub.train.class)[1]/length(sub.train.class))
			
			#calculate the training error
			if(length(classifiers[[i]])!=0){
				train.cl.pred[[i]] <- as.numeric(as.character(predict(classifiers[[i]], t(gs.train.exprs))))*cl.label[i]
			}
			
			gs.val.exprs <- val.exprs[gs.ind, ]
			gs.test.exprs <- test.exprs[gs.ind, ]
			
			if(length(classifiers[[i]])!=0){
				val.cl.pred[[i]] <- as.numeric(as.character(predict(classifiers[[i]], t(gs.val.exprs))))*cl.label[i]
				test.cl.pred[[i]] <- as.numeric(as.character(predict(classifiers[[i]], t(gs.test.exprs))))*cl.label[i]
			}
		}
		
		if (max.iter!=0){
			max.iter <- min(n.iter, max.iter)
		}else{
			max.iter <- n.iter
		}
		
		#calculate the training error
		train.auc.arr <- matrix(data = NA, nrow = max.iter, ncol = 6)
		
		for (i in 1:max.iter) {
			train.step.pred <- sign((alpha[1:i]) %*% t(do.call(cbind, train.cl.pred[1:i])))
			train.auc.arr[i, ] <- unlist(get.perf(as.numeric(as.character(train.step.pred)), train.class))
		}
		
		val.auc.arr <- matrix(data = NA, nrow = max.iter, ncol = 6)
		
		for (i in 1:max.iter) {
			val.step.pred <- sign((alpha[1:i]) %*% t(do.call(cbind, val.cl.pred[1:i])))
			val.auc.arr[i, ] <- unlist(get.perf(as.numeric(as.character(val.step.pred)), val.class))
		}
		
		test.auc.arr <- matrix(data = NA, nrow = max.iter, ncol = 6)
		
		for (i in 1:max.iter) {
			test.step.pred <- sign((alpha[1:i]) %*% t(do.call(cbind, test.cl.pred[1:i])))
			test.auc.arr[i, ] <- unlist(get.perf(as.numeric(as.character(test.step.pred)), test.class))
		}
		
		#plot error over iterations
		jpeg(filename = paste("error plot ", " para1 ", para1, ".jpg", sep=""))
		plot(seq(1, max.iter), errarr[1:max.iter])
		dev.off()
		
		#plot alpha over iterations
		jpeg(filename = paste("alpha plot ", " para1 ", para1, ".jpg",sep=""))
		plot(seq(1, max.iter), alpha[1:max.iter])
		dev.off()
		
		#plot the training auc
		jpeg(filename = paste("train auc plot ", " para1 ", para1 , ".jpg",sep=""))
		plot(seq(1, n.iter), train.auc.arr[, 1])
		dev.off()
		
		#plot the testing auc
		jpeg(filename = paste("val auc plot ", " para1 ", para1 , ".jpg",sep=""))
		plot(seq(1, n.iter), unlist(val.auc.arr[, 1]))
		dev.off()
		
		#plot sample ratio in training at each step
		jpeg(filename = paste("sample ratio", " para1 ", para1 , ".jpg",sep=""))
		plot(seq(1, max.iter), unlist(train.ratio.arr[1:max.iter]))
		dev.off()
		
		val.perf.df<- unlist(get.perf(as.numeric(as.character(val.step.pred)), val.class))
		test.perf.df<- unlist(get.perf(as.numeric(as.character(test.step.pred)), test.class))
		
		#write parameters to the files
		write(errarr, paste("err para1", " para1 ", para1 , ".txt", sep=""), append=T, ncolumns=1000)
		write(alpha, paste("alpha para1" , " para1 ", para1, ".txt", sep=""), append=T, ncolumns=1000)
		write(k, paste("val.auc para1", " para1 ", para1, ".txt", sep=""), append=T)
		write.table(val.auc.arr, paste("val.auc para1", " para1 ", para1 , ".txt", sep=""), append=T)
		write(k, paste("training.auc para1", " para1 ", para1, ".txt", sep=""), append=T)
		write.table(train.auc.arr, paste("training.auc para1", " para1 ", para1 , ".txt", sep=""), append=T)
		write(w, paste("weight para1", " para1 ", para1 , ".txt", sep=""), append=T, ncolumns=1000)	
	}, error = function(err){
		boost.perf.df="error here"
	})
	
	write.table(t(c(kegg.gs.anno[k,], val.perf.df)), paste("/C/Zhe Liu/DevGSEA/result/NKI/test set blind/val.res.df ", "para1 ",para1, ".txt", sep=""), append=T, col.names=F, row.names=F)
	write.table(t(c(kegg.gs.anno[k,], test.perf.df)), paste("/C/Zhe Liu/DevGSEA/result/NKI/test set blind/test.res.df ", "para1 ",para1, ".txt", sep=""), append=T, col.names=F, row.names=F)
}

write(c("sampling ratio", para1), paste("/C/Zhe Liu/DevGSEA/result/NKI/test set blind/parameter ", "para1 ", para1, ".txt", sep=""), append=T)
write(c("n.iteration", n.iter), paste("/C/Zhe Liu/DevGSEA/result/NKI/test set blind/parameter ", "para1 ", para1, ".txt", sep=""), append=T)
write(c("svm: balanced sub.train.class", "boost: viola balanced"), paste("/C/Zhe Liu/DevGSEA/result/NKI/test set blind/parameter ", "para1 ", para1, ".txt", sep=""), append=T)

save.image("/C/Zhe Liu/DevGSEA/result/NKI/test set blind/test set blind.Rdata")
#load("/C/Zhe Liu/DevGSEA/result/NKI/test set blind/combined weighted voting.Rdata")

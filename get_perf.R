# function to retrieve the performance results
get.perf<- function(pred, golden.st){
  roc.pred=prediction(pred, golden.st)
  roc.pref=performance(roc.pred, "auc")
  auc=round(as.numeric(roc.pref@y.values), 3)
  
  #fpn, tpn, fnn, tnn of the prediction performance
  fpn=roc.pred@fp[[1]][2]
  tpn=roc.pred@tp[[1]][2]
  
  fnn=roc.pred@fn[[1]][2]
  tnn=roc.pred@tn[[1]][2]
  
  #classification accuracy 
  acc=round(((tpn+tnn)/sum(tpn+tnn+fpn+fnn)), 3)
  
  #fpr and fnr
  fpr=round((fpn/(fpn+tnn)), 3)
  fnr=round((fnn/(tpn+fnn)), 3)
  
  res <- list(auc, acc, fpr, fnr, fpn, fnn)
  names(res) <- c("auc", "accuracy", "fpr", "fnr", "fpn", "fnn")
  return(res)
}

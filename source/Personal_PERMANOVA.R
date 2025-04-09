#Personal_PERMANOVA
#function to perform permanova analysis
Personal_PERMANOVA<-function(data.repeat,factors,response,algo,type,transform="4root",method='bray',binary=FALSE,only.intercept=FALSE,strata=NULL,nperm=999){
  #set.seed(123)
  fish.mat=Transformation_data(data.repeat,transform)
  #fish.mat=Transformation_data(Transformation_data(data.repeat,"root"),"Hellinger")
  perm <- how(nperm = nperm)
  if (is.null(strata)==FALSE){
    setBlocks(perm) <- with(data.repeat, fDate)
  }
  if (any(rowSums(fish.mat)==0)){
    data.repeat=data.repeat[-which(rowSums(fish.mat)==0),]
    fish.mat=fish.mat[-which(rowSums(fish.mat)==0),]
  }
  fish.dist=vegdist(fish.mat, method=method,binary=binary)
  if (type=="adonis2"){
    .form=as.formula(paste(response,"~", paste(factors, collapse=algo)))
    fish.div<-adonis2(.form, data=data.repeat, permutations = perm, method=method)
  }
  if (type=="adonis"){
    .form=as.formula(paste(response,"~", paste(factors, collapse=algo)))
    fish.div<-adonis(.form, data=data.repeat, permutations = perm, method=method)
  }
  if (only.intercept==TRUE){
    .form=as.formula(paste(response,"~1"))
    fish.div<-adonis(.form, data=data.repeat, permutations = perm, method=method)
  }
  print(fish.div)
  #print(densityplot(permustats(fish.div)))
  return(fish.div)
}
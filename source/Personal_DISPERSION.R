#Personal_DISPERSION
#function to perfrom permdisp and make plots
Personal_DISPERSION<-function(data.repeat,group){
  set.seed(123)
  fish.mat=Transformation_data(data.repeat,"4root")
  #fish.mat=Transformation_data(Transformation_data(N1,"root"),"Hellinger")
  if (any(rowSums(fish.mat)==0)){
    data.repeat=data.repeat[-which(rowSums(fish.mat)==0),]
    fish.mat=fish.mat[-which(rowSums(fish.mat)==0),]
  }
  fish.dist=vegdist(fish.mat, method='bray')
  dispersion<-betadisper(fish.dist, group=lapply(data.repeat[,group],as.character)[[1]])
  dispersion
  #print(anova(dispersion))
  boxplot(dispersion) #plot average distance to the centroid
  print(permutest(dispersion))
  plot(dispersion, hull=TRUE, ellipse=TRUE,main=group) #sd ellipse: Higher spread: Higher betadiversity  
  return(dispersion)
}
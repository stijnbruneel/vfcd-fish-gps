#Transformation_personal
#function to transform the biological data 
Transformation_data<-function(x,type){
  x=x[,c(which(colnames(x) %in% colnames(specieslist)))]
  x[is.na(x)==TRUE]=0
  if (type=="root"){
    x<-sqrt(x)
  }
  if (type=="4root"){
    x<-sqrt(sqrt(x))
  }
  if (type=="no transformation"){
    x<-x
  }
  if (type=="log"){
    x<-log(x+1)
  }
  if (type=="log10"){
    x<-log10(x+1)
  }
  if (type=="presence-absence"){
    x[x>0] <-1
  }
  if (type=="Hellinger"){
    x=decostand(x,"hellinger")
  }
  return(x)
}
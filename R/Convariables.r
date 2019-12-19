###################################################################
#######Develop conditional variable################################

GetConVariable=function(Y){
  convar=GetConVar(Y)
  return(convar)
}

ConVariable=function(y,X){
  X=as.matrix(X)
  Cyx=cov(y,X)
  Cx=cov(X)
  k=Cyx%*%solve(Cx)
  mux=apply(X,2,mean)
  x=t(X)-mux
  CY=y-t(k%*%x)
  return(CY)
}

GetConVar=function(Y){
  TraitNames=colnames(Y)
  ComNum=length(TraitNames)-1
  e=1
  y=Y[,-(1:ComNum)]
  for(i in 1:ComNum){
    if(i==ComNum)C=matrix(1:i)
    else C=combn(ComNum,i)
    nr=ncol(C)
    name=NULL
    for(k in 1:nr){
      for(j in 1:i){
        if(j==1){
          name[k]=TraitNames[C[j,k]]
        }
        else if(j>1){
          name[k]=paste(name[k],TraitNames[C[j,k]],sep="&")
        }
      }
      x=Y[,C[,k]]
      x=as.matrix(x)
      if(e==1)ConY=ConVariable(y,x)
      else if(e>1)ConY=cbind(ConY,ConVariable(y,x))
      e=e+1
    }
    if(i==1)ConNames=name
    else if(i>1)ConNames=c(ConNames,name)
  }
  ConNum=length(ConNames)

  colnames(ConY)=paste(TraitNames[ComNum+1],ConNames,sep="|")
  ConY=cbind(ConY,y)
  colnames(ConY)[ConNum+1]=TraitNames[ComNum+1]
  con=list(Y=Y,ConY=ConY)
  return(con)
}

ConVariable1=function(y,X){
  data=data.frame(y,X)
  mod=lm(y~X,data=data)
  return(mod$residuals)
}

GetConVar1=function(Y){
  #require(gtools)
  TraitNames=colnames(Y)
  ComNum=length(TraitNames)-1
  e=1
  y=Y[,-(1:ComNum)]
  for(i in 1:ComNum){
    C=combn(ComNum,i)
    nr=nrow(C)
    name=NULL
    for(k in 1:nr){
      for(j in 1:i){
        if(j==1){
          name[k]=TraitNames[C[k,j]]
        }
        else if(j>1){
          name[k]=paste(name[k],TraitNames[C[k,j]],sep="&")
        }
      }
      x=Y[,C[k,]]
      x=as.matrix(x)
      if(e==1)ConY=ConVariable1(y,x)
      else if(e>1)ConY=cbind(ConY,ConVariable1(y,x))
      e=e+1
    }
    if(i==1)ConNames=name
    else if(i>1)ConNames=c(ConNames,name)
  }
  ConNum=length(ConNames)

  colnames(ConY)=paste(TraitNames[ComNum+1],ConNames,sep="|")
  ConY=cbind(ConY,y)
  colnames(ConY)[ConNum+1]=TraitNames[ComNum+1]
  con=list(Y=Y,ConY=ConY)
  return(con)
}




##Functions used for linear mixed model analysis
## needed Dec 18, 2019

lmm.check=function (formula, data = list())
{
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    mk=length(gdata$U)
    xk=length(gdata$X)
    nx=numeric(xk)
    for(i in 1:xk){
       if(i==1)FixedNames=colnames(gdata$X[[i]])
       else FixedNames=c(FixedNames,colnames(gdata$X[[i]]))
       nx[i]=ncol(gdata$X[[i]])
    }
    nu=numeric(mk)
    for(i in 1:mk){
      if(i==1)RandomNames=colnames(gdata$U[[i]])
      else RandomNames=c(RandomNames,colnames(gdata$U[[i]]))
      nu[i]=ncol(gdata$U[[i]])
    }

    names(nx)=names(gdata$X)
    names(nu)=names(gdata$U)
    res=list(VarCompNum=(mk+1),VarCompNames=names(gdata$U),FixedCompNum=nx,RandomCompNum=nu,FixedEffectNames=FixedNames)
    return(res)
}

## ## needed Dec 18, 2019
lmm.mq.simu=function(formula, data = list(),v0,b0,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
    if (is.null(data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data)
    if(is.null(SimuNum))SimuNum=200
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(ALPHA))ALPHA=0.05
    result=genmod.simuold(gdata,v=v0,b=b0,SimuNum=SimuNum,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
    return(result)
}

## Needed Dec 18, 2019
ginv=function (X, tol = sqrt(.Machine$double.eps)){

    if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X)))
        stop("'X' must be a numeric or complex matrix")
    if (!is.matrix(X))
        X <- as.matrix(X)
    Xsvd <- svd(X)
    if (is.complex(X))
        Xsvd$u <- Conj(Xsvd$u)
    Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
    if (all(Positive))
        Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
    else if (!any(Positive))
        array(0, dim(X)[2L:1L])
    else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) *
        t(Xsvd$u[, Positive, drop = FALSE]))
}









## needed Dec 18, 2019

MINQUEVarPre=function(Qa,X,U,ML,MR,y,C=NULL,VI=NULL){
   svdinv=LargeInv
   #ML1=as.matrix(ML)
   #V0=ginv(ML)%*%MR ##Estimate variance components
   V0=svdinv(ML)%*%MR ##Estimate variance components
   v=as.vector(V0)
   mk=length(U)
   n=nrow(X)
   if(is.null(VI)){
      if(mk==0)VI=diag(1/v[1],n)
      else if(mk>0){
         id=which(v<=0)
         sumv=sum(v[id])
         smallvalue=sumv/10^3
         id=which(v<smallvalue)
         v[id]=smallvalue
         VI=Woodbury(v,U)
      }
   }
   #X0=X
   if(ncol(X)==1){
       HF=t(X[,1])%*%VI%*%X[,1]
       bhat=svdinv(HF)%*%(t(X[,1])%*%VI%*%y)  ##Estimate fixed effects
   }
   if(is.null(C)){
       HF=t(X)%*%VI%*%X
       bhat=svdinv(HF)%*%(t(X)%*%VI%*%y)
   }
   else{
      HF=t(X)%*%VI%*%X+C%*%t(C)
      bhat=ginv(HF)%*%(t(X)%*%VI%*%y)
   }
   #dim(X)
   #n=ncol(U[[1]])
   #mk=length(U)
   if(mk>0){
     k=numeric()
     for(i in 1:mk){
        nu=ncol(U[[i]])-1
        if(v[i]<=0){
          k[i]=0
          ##v[i]=0
        }
        else if(v[i]>0)k[i]=sqrt(nu*v[i]/MR[i])
     }
     PRE=numeric()
     vc=numeric()

     for(i in 1:mk){
        if(k[i]==0){
          es=rep(0,ncol(U[[i]]))
          PRE=c(PRE,es)
        }
        else if(k[i]>0){
          es<-k[i]*t(U[[i]])%*%Qa%*%y
          #es<-t(U[[i]])%*%Qa%*%y
          PRE=c(PRE,es)
        }
        #vc[i]=var(es)
     }
   }
   if(mk==0)result=list(v=v,b=bhat,Pre=NULL,HF=HF)
   else result=list(v=v,b=bhat,Pre=PRE,HF=HF)
   return(result)
}


## needed Dec 18, 2019

GetQa=function(au,U,X,C=NULL){
  svdinv=LargeInv
  mk=length(U)
  n=nrow(X)
  if(mk>0)VI=Woodbury(au,U)

  else if(mk==0)VI=diag(1/au,n)
  if(ncol(X)==1)VX=VI%*%X[,1]
  else VX=VI%*%X
  if(is.null(C))Qa=VI-VX%*%solve(t(X)%*%VX)%*%t(VX)
  else Qa=VI-VX%*%svdinv(t(X)%*%VX+C%*%t(C))%*%t(VX)
  res=list(Qa=Qa,VI=VI)
  return(res)

}

## Calculate the left side for MINQUE normal equations

## needed Dec 18, 2019

MINQUE_L=function(Qa,U){
  mk=length(U)
  ML=matrix(0, nrow=mk+1,ncol=mk+1)
  if(mk>0){
  for(i in 1:mk){
     for(j in i:mk){
        a1=t(U[[i]])%*%Qa%*%U[[j]]
        a=a1%*%t(a1)
        b=sum(diag(a))
        ML[i,j]=ML[j,i]=b[1]
     }
     a1=t(U[[i]])%*%Qa
     a=a1%*%t(a1)
     b=sum(diag(a))
     ML[i,(mk+1)]=ML[(mk+1),i]=b[1]
  }
  }

  a1=Qa
  a=a1%*%t(a1)
  b=sum(diag(a))
  ML[(mk+1),(mk+1)]=b[1]

  return(ML)
}

## needed Dec 18, 2019
MINQUE_R=function(Qa,U,y){
  mk=length(U)
  MR<-numeric()
  y=as.numeric(y)
  if(mk>0){
     for(i in 1:mk){
       #i=1
       t=t(y)%*%Qa%*%U[[i]]
       a=t%*%t(t)
       MR[i]=a[1,1]
     }
  }
  t=t(y)%*%Qa
  a=t%*%t(t)
  MR[mk+1]=a[1,1]
  return(MR)
}

##y=YD[,1]
## needed Dec 18, 2019

GetVarPre=function(Qa,U,X,ML,y){
   svdinv=LargeInv
   MR=MINQUE_R(Qa,U,y)
   V0=svdinv(ML)%*%MR


   v=as.vector(V0)
   mk=length(U)
   n=nrow(X)
   if(mk==0)VI=diag(1,n)
   else VI=Woodbury(v,U)
   X0=as.matrix(X)
   ##Estimate fixed effects by GLSE
   if(ncol(X)==1) bhat=solve(t(X[,1])%*%VI%*%X[,1])%*%(t(X[,1])%*%VI%*%y)  ## I don't know why, but it works this way
   else if(ncol(X)>1)bhat=solve(t(X0)%*%VI%*%X0)%*%(t(X0)%*%VI%*%y)

   #n=ncol(U[[1]])
   #mk=length(U)
   #v=as.vector(V0)
   if(mk>0){
     k=numeric()
     for(i in 1:mk){
       nu=ncol(U[[i]])-1
       if(v[i]<=0){
          k[i]=0
          ##v[i]=0
       }
       else if(v[i]>0)k[i]=sqrt(nu*v[i]/MR[i])
     }
     PRE=numeric()
     vc=numeric()

     for(i in 1:mk){
      if(k[i]==0){
          es=rep(0,ncol(U[[i]]))
          PRE=c(PRE,es)
      }
      else if(k[i]>0){
         es<-k[i]*t(U[[i]])%*%Qa%*%y
         #es<-t(U[[i]])%*%Qa%*%y
         PRE=c(PRE,es)
      }
      #vc[i]=var(es)
    }
   }
   if(mk==0)PRE=NULL
   result=list(v=v,b=bhat,Pre=PRE)
   return(result)
}



#######################################################################################
####### Function to Drop the columns of any matrix in which all elements are zero.#####
#######################################################################################

## Needed Dec 18, 2019

REML=function(U,X,y,ITMAX=NULL,...){
    #U=U1
    #y=yd
    if(is.null(ITMAX))ITMAX=5
    SMALLVALUE=0.001
    DIS=100
    it=1
    mk=length(U)
    v0=rep(1,mk+1)
    S=which(v0<=0)
    v0[S]=0

    while(DIS>SMALLVALUE&&it<=ITMAX){

        au=v0
        Qa=GetQa(au,U,X)$Qa
        ML=MINQUE_L(Qa,U)
        result=GetVarPre(Qa,U,ML,y)
        v1=result[[1]]
        S=which(v1<=0)
        v1[S]=0
        if(it==1)ITV=v1
        else if(it>1)ITV=cbind(ITV,v1)
        it=it+1
        d=v0-v1
        DIS=sqrt(mean(d^2)/sum(v0))
        v0=v1
    }
    Pre=result[[2]]
    reml=list(v1,Pre)
    return(reml)
}

## Needed Dec 18, 2019

Woodbury=function(au=NULL,U,R=NULL){
   svdinv=LargeInv
   mk=length(U)
   if(is.null(au))au=rep(1,mk+1)
   n=nrow(U[[1]])
   if(is.null(R)){
      V=diag(1/(au[mk+1]),n)
      for(i in 1:mk){
         a=au[i]
         if(a>0){
            c=ncol(U[[i]])
            m1=V%*%U[[i]]
            a=au[i]
            d=diag(1,c)
            vi=d+a*t(U[[i]])%*%m1
            vi=svdinv(vi)
            VI=V-a*m1%*%vi%*%t(m1)
            V=VI
         }
      }
      VI=V
   }
   else{
      V=R
      for(i in 1:mk){
         a=au[i]
         if(a>0)V=V+a*U[[i]]%*%t(U[[i]])
      }
      VI=svdinv(V)
   }
   return(VI)
}


## Needed Dec 18, 2019
SimuData=function(gdata,v,b,SimuNum=NULL){
   if(is.null(SimuNum))SimuNum=200
   gdata=SimuData_Old(gdata,v,b,SimuNum)
   return(gdata)
}

## Needed Dec 18, 2019
SimuData_Old=function(gdata,v,b,SimuNum=NULL){
  if(is.null(SimuNum))SimuNum=200
  U=gdata$U
  nk=length(gdata$X)
  if(nk==1)X=gdata$X[[1]]
  else{
    for(i in 1:nk){
       if(i==1)X=gdata$X[[i]]
       else X=cbind(X,gdata$X[[i]])
    }
  }
  mk=length(U)
  n=nrow(U[[1]])
  ct=numeric()
  for(i in 1:mk){
     ct[i]=ncol(U[[i]])
     if(i==1)BigU=U[[i]]
     else if(i>1)BigU=cbind(BigU,U[[i]])
  }
  #dim(X)
  mt=sum(ct)
  #b=as.matrix(b)
  RE0=matrix(0,nrow=mt,ncol=SimuNum)
  YS=matrix(0,nrow=n,ncol=SimuNum)
  for(i in 1:SimuNum){
     #i=1
     for(j in 1:mk){
        #j=1
        if(v[j]>0)a=rnorm(ct[j],0,sqrt(v[j]))
        if(v[j]<=0)a=rep(0,ct[j])
        if(j==1)RE=a
        else if(j>1)RE=c(RE,a)
     }
     RE0[,i]=RE
     #RE=as.matrix(RE)
     YS[,i]=X%*%b+BigU%*%RE+rnorm(n,0,sqrt(v[mk+1]))

  }
  #res=list(SimuY=YS,SimuEffect=RE0)
  gdata$Y=YS
  gdata$SimuEffect=RE0
  gdata$VC=rep(1,mk+1)
  return(gdata)
}




## Needed Dec 18, 2019
svdinv=function(A){
   return(LargeInv(A))
}



## A function used for Monte Carlo simulations
## Needed Dec 18, 2019
genmod.simuold=function(gdata,v,b,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL,...){
   if(is.null(SimuNum))SimuNum=200
   if(is.null(ALPHA))ALPHA=0.05
   gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
   if(is.null(JacNum))JacNum=10
   if(is.null(JacRep))JacRep=1

   jac=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep,LUP=1)
   ml=length(v)
   Var=jac$Var
   P=numeric(ml)
   vm=numeric(ml)
   for(i in 1:SimuNum){
      vm=vm+Var[[i]][,1]
      id=which(Var[[i]][,3]<=ALPHA)
      P[id]=P[id]+1
   }
   vm=vm/SimuNum
   P=P/SimuNum
   bias=vm-v
   a=data.frame(v,vm,bias,P)
   colnames(a)=c("True","Estimate","Bias","Power")
   rownames(a)=rownames(Var[[1]])
   res=list(P=a,ALPHA=ALPHA)
   return(res)
}


## Needed Dec 18, 2019
mixed.data=function(formula,data=list(),...){
   if(is.null(data))mq=lm.data(formula=formula)
   else mq=lm.data(formula=formula,data=data)
   Y=mq$Y
   X=mq$X
   U=mq$U
   xk=length(X)
   nx=0
   for(i in 1:xk)nx=nx+ncol(X[[i]])
   if(xk==1)C=NULL
   if(xk>1){
      C=matrix(0,nx,xk-1)
      e1=ncol(X[[1]])+1
      for(i in 2:xk){
         e2=e1+ncol(X[[i]])-1
         C[(e1:e2),(i-1)]=1
         e1=e2+1
      }
   }

   #for(i in 1:xk)X=cbind(X,mq$X[[i]])
   mk=length(U)
   VC=rep(1,(mk+1))
   if(class(Y)=="numeric")Y=data.frame(Y)
   result=list(Y=Y,X=X,U=U,VC=VC,C=C,Model=formula)
   result$call=match.call()
   class(result)="genmod.data"
   return(result)
}

## Needed Dec 18, 2019
genmod=function(gdata,au=NULL,LUP=NULL,...){
   Y=gdata$Y
   U=gdata$U
   C=gdata$C
   ml=length(gdata$U)+1
   if(is.null(au))au=rep(1,ml)
   xk=length(gdata$X)
   if(xk==1)X=gdata$X[[1]]
   if(xk>1)X=BigX(gdata$X)
   nx=ncol(X)
   #for(i in 1:nx)X[,i]=as.numeric(X[,i])
   VC=gdata$VC
   #Model=gdata$model
   #BigX=NULL
   #class(X[,3])
   #X=as.matrix(X)
   if(is.null(LUP))LUP=NULL
   #result=minque.default(X,U,Y,C,au)
   result=minque.default(X,U,Y,C,au,LUP)
   result$Var=result$Var*VC
   #result$call=match.call()
   class(result)="minque"
   return(result)
}

## Needed Dec 13, 2019
minque.default=function(X,U,Y,C=NULL,au=NULL,LUP=NULL){
  mk=length(U)
  n=nrow(U[[1]])
  if(is.null(au))au=rep(1,(mk+1))
  if(is.null(X))X=matrix(1,nrow=n,ncol=1)
  if(is.null(C))C=NULL
  xnames=colnames(X)
  if(mk==0)vnames="V(e)"
  else if(mk>0){
    rnames=NULL
    vnames=c(paste("V(",names(U),")",sep=""),"V(e)")
    for(i in 1:mk)rnames=c(rnames,colnames(U[[i]]))
  }

  aa=GetQa(au,U,X,C)
  Qa=aa$Qa
  VI=aa$VI
  ML=MINQUE_L(Qa,U)
  #if(class(Y)=="numeric")
  Y=data.frame(Y)
  TraitNum=ncol(Y)
  TraitNames=colnames(Y)
  Var=matrix(0,(mk+1),TraitNum)
  B=matrix(0,ncol(X),TraitNum)
  hf=list()
  if(mk>0)Pre=matrix(0,length(rnames),TraitNum)
  if(is.null(LUP))VI=NULL
  for(i in 1:TraitNum){
    MR=MINQUE_R(Qa,U,Y[,i])
    a=MINQUEVarPre(Qa,X,U,ML,MR,Y[,i],C,VI)
    Var[,i]=a[[1]]
    B[,i]=a[[2]]
    if(mk>0)Pre[,i]=a[[3]]
    hf[[i]]=a$HF
  }

  colnames(Var)=TraitNames
  rownames(Var)=vnames
  colnames(B)=TraitNames
  rownames(B)=xnames
  names(hf)=TraitNames
  if(mk>0){
    colnames(Pre)=TraitNames
    rownames(Pre)=rnames
  }
  if(mk>0)result=list(Var=data.frame(Var),
                      FixEffect=data.frame(B),
                      RandomEffect=data.frame(Pre),
                      HV=ML,
                      Qa=Qa,
                      HF=hf)
  if(mk==0)result=list(Var=data.frame(Var),
                       FixEffect=data.frame(B),
                       RandomEffect=NULL,
                       HV=ML,
                       Qa=Qa,
                       HF=hf)
  result$call=match.call()
  return(result)
}




#gdata=mod
print.genmod.data=function(gdata,...){
   cat("\nCall:\n")
   print(gdata$call)

   cat("\nThe information data are very complicated. We only provide summarized information as follows:\n")

   cat("\nSample size=", nrow(gdata$Y),"\n")
   cat("\nNo of traits=", ncol(gdata$Y),"\n")

   cat("\nNo of information matrices X for fixed effects is: ", length(gdata$X), "\n")
   mk=length(gdata$U)
   if(mk==0) cat("\nNo matrix for ramdom effects available\n")
   else if (mk>0) cat("\nNo of information matrices for random effects is:", length(gdata$U), "\n")

}



#minque=function(Y,formula,...) UseMethod("minque")


print.minque=function(mq,...){
   cat("Call: \n")
   print(mq$call)
   cat("\n")

   cat("\nMinque approach for variance components and for fixed and random effects\n")

   cat("\nVariance components:\n")
   print(mq$Var)
   cat("\n")

   cat("Fixed effects:\n")
   print(mq$FixEffect)
   cat("\n")

   cat("Random effects:\n")
   print(mq$RandomEffect)
   cat("\n")
}


print.minque.jack=function(mq,...){
   cat("Call: \n")
   print(mq$call)
   cat("\n")

   cat("\nJackknife technique with randomized grouping\n")

   cat("\nVariance components:\n")
   print(mq$Var)
   cat("\n")

   cat("Proportional variance components:\n")
   print(mq$PVar)
   cat("\n")

   cat("Fixed effects:\n")
   print(mq$FixedEffect)
   cat("\n")

   cat("Random effects:\n")
   print(mq$RandomEffect)
   cat("\n")
}

## Needed Dec 18, 2019
genmod.perm=function(gdata,au=NULL,PermNum=NULL,LUP=NULL,...){
    if(is.null(au))au=NULL
    if(is.null(PermNum))PermNum=200
    if(is.null(LUP))LUP=NULL
    #if(is.null(SIMUYES))result=genmod.oldperm(gdata,PermNum)
    #else
    result=genmod.oldperm(gdata,au,PermNum,LUP)
    return(result)

}


## Needed Dec 18, 2019
genmod.oldperm=function(gdata,au=NULL,PermNum,LUP=NULL,...){
   if(is.null(LUP))LUP=NULL
   if(is.null(au))au=NULL
   result0=mq(gdata)
   tn=ncol(gdata$Y)
   n=nrow(gdata$Y)
   Y0=gdata$Y
   gc=length(gdata$VC)
   mc=nrow(result0$RandomEffect)
   fc=nrow(result0$FixedEffect)
   VP=matrix(0,gc,tn)
   VT=list()
   ET=list()
   FT=list()
   for(i in 1:tn){
      vp=numeric(gc)
      ep=numeric(mc)
      fp=numeric(fc)
      y=Y0[,i]
      Y=matrix(0,n,PermNum)
      for(j in 1:PermNum)Y[,j]=sample(y)
      colnames(Y)=paste("R",1:PermNum,sep="")
      gdata$Y=Y
      result=mq(gdata)
      V=result$Var
      v0=result0$Var[,i]
      e0=result0$RandomEffect[,i]
      f0=result0$FixedEffect[,i]
      PE=result$RandomEffect
      FE=result$FixedEffect
      for(j in 1:gc){
         id=which(V[j,]>v0[j])
         vp[j]=length(id)/PermNum
      }
      for(j in 1:mc){
         if(e0[j]<0){
             id=which(PE[j,]<e0[j])
             ep[j]=length(id)/PermNum
         }
         if(e0[j]>=0){
             id=which(PE[j,]>e0[j])
             ep[j]=length(id)/PermNum
         }
      }
      for(j in 1:fc){
         if(f0[j]<0){
             id=which(FE[j,]<f0[j])
             fp[j]=length(id)/PermNum
         }
         if(f0[j]>=0){
             id=which(FE[j,]>f0[j])
             fp[j]=length(id)/PermNum
         }
      }

      vp[gc]=1-vp[gc]
      VT[[i]]=cbind(v0,vp)
      ET[[i]]=cbind(e0,ep)
      FT[[i]]=cbind(f0,fp)
      colnames(VT[[i]])=c("Est","Pvalue")
      rownames(VT[[i]])=rownames(result0$Var)

      colnames(ET[[i]])=c("Pre","Pvalue")
      rownames(ET[[i]])=rownames(result0$RandomEffect)

      colnames(FT[[i]])=c("Est","Pvalue")
      rownames(FT[[i]])=rownames(result0$FixedEffect)

  }
  names(VT)=colnames(result0$Var)
  names(ET)=colnames(result0$Var)

  res=list(Var=VT,RandomEffect=ET,FixedEffect=FT)
  return(res)
}

## Needed Dec 19, 2019
reml.perm=function(gdata,au=NULL,PermNum,LUP=NULL,...){
   if(is.null(LUP))LUP=NULL
   if(is.null(au))au=NULL
   result0=reml(gdata)
   tn=ncol(gdata$Y)
   n=nrow(gdata$Y)
   Y0=gdata$Y
   gc=length(gdata$VC)
   mc=nrow(result0$RandomEffect)
   fc=nrow(result0$FixedEffect)
   VP=matrix(0,gc,tn)
   VT=list()
   ET=list()
   FT=list()
   for(i in 1:tn){
      vp=numeric(gc)
      ep=numeric(mc)
      fp=numeric(fc)
      y=Y0[,i]
      Y=matrix(0,n,PermNum)
      for(j in 1:PermNum)Y[,j]=sample(y)
      colnames(Y)=paste("R",1:PermNum,sep="")
      gdata$Y=Y
      result=reml(gdata)
      V=result$Var
      v0=result0$Var[,i]
      e0=result0$RandomEffect[,i]
      f0=result0$FixedEffect[,i]
      PE=result$RandomEffect
      FE=result$FixedEffect
      for(j in 1:gc){
         id=which(V[j,]>v0[j])
         vp[j]=length(id)/PermNum
      }
      for(j in 1:mc){
         if(e0[j]<0){
             id=which(PE[j,]<e0[j])
             ep[j]=length(id)/PermNum
         }
         if(e0[j]>=0){
             id=which(PE[j,]>e0[j])
             ep[j]=length(id)/PermNum
         }
      }
      for(j in 1:fc){
         if(f0[j]<0){
             id=which(FE[j,]<f0[j])
             fp[j]=length(id)/PermNum
         }
         if(f0[j]>=0){
             id=which(FE[j,]>f0[j])
             fp[j]=length(id)/PermNum
         }
      }

      vp[gc]=1-vp[gc]
      VT[[i]]=cbind(v0,vp)
      ET[[i]]=cbind(e0,ep)
      FT[[i]]=cbind(f0,fp)
      colnames(VT[[i]])=c("Est","Pvalue")
      rownames(VT[[i]])=rownames(result0$Var)

      colnames(ET[[i]])=c("Pre","Pvalue")
      rownames(ET[[i]])=rownames(result0$RandomEffect)

      colnames(FT[[i]])=c("Est","Pvalue")
      rownames(FT[[i]])=rownames(result0$FixedEffect)

  }
  names(VT)=colnames(result0$Var)
  names(ET)=colnames(result0$Var)

  res=list(Var=VT,RandomEffect=ET,FixedEffect=FT)
  return(res)
}


#SIMUYES=1
## Needed Dec 19, 2019
genmod.jack=function(gdata,JacNum=NULL,JacRep=NULL,ALPHA=NULL,au=NULL,LUP=NULL,...){
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  if(is.null(ALPHA))ALPHA=0.05
  if(is.null(au))au=NULL
  if(is.null(LUP))LUP=NULL
  #if(is.null(SIMUYES))result=genmod.jack1(gdata=gdata,JacNum=JacNum,JacRep=JacRep,au=au,ALPHA=ALPHA)
  result=genmod.jack1(gdata,JacNum,JacRep,ALPHA,au,LUP)
}
## needed Dec 19, 2019
genmod0=function(gdata,au=NULL,LUP=NULL){
   gdata0=gdata
   if(is.null(au))au=NULL
   if(is.null(LUP))LUP=NULL
   res0=genmod(gdata0,au,LUP)
}

#JacNum=10
#JacRep=1
## Needed Dec 19, 2019
genmod.jack1=function(gdata,JacNum,JacRep,ALPHA,au=NULL,LUP=NULL,...){
   if(is.null(au))au=NULL
   if(is.null(LUP))LUP=NULL
   a=1-(1-ALPHA)^(1/(JacNum-2))
   t=qt(1-a/2,JacNum-1)
   Y=as.matrix(gdata$Y)
   TraitNum=ncol(Y)
   TraitNames=colnames(Y)
   n=nrow(Y)
   #colnames(Y)
   mk=length(gdata$U)
   xk=length(gdata$X)
   JB=rep(1:JacNum,length=n)
   gdata0=gdata
   ml=mk+1
   nx=0
   for(i in 1:xk)nx=nx+ncol(gdata0$X[[i]])

   res0=genmod0(gdata0,au=au,LUP=LUP)
   V0=res0$Var
   B0=res0$FixEffect
   P0=res0$RandomEffect
   if(mk>0){
      mt=nrow(P0)
      Str=rownames(P0)
   }
   VNames=rownames(V0)
   XNames=rownames(B0)

   ##############################################################################
   RES0=NULL
   VPower=numeric(ml)
   V0=matrix(0,TraitNum,ml)
   if(mk>0)P0=matrix(0,TraitNum,mt)
   #colnames(P0)
   Vp=matrix(0,TraitNum,ml)
   Cv=matrix(0,TraitNum,ml)
   Cv1=matrix(0,TraitNum,ml)
   Cv2=matrix(0,TraitNum,ml)
   Vm=matrix(0,TraitNum,ml)

   VCp=matrix(0,TraitNum,ml)
   VCv=matrix(0,TraitNum,ml)
   VCv1=matrix(0,TraitNum,ml)
   VCv2=matrix(0,TraitNum,ml)
   VCm=matrix(0,TraitNum,ml)

   if(mk>0){
      Pp=matrix(0,TraitNum,mt)
      Cp=matrix(0,TraitNum,mt)
      Cp1=matrix(0,TraitNum,mt)
      Cp2=matrix(0,TraitNum,mt)
      Pm=matrix(0,TraitNum,mt)
   }

   Cb1=matrix(0,TraitNum,nx)
   Cb2=matrix(0,TraitNum,nx)
   Cb=matrix(0,TraitNum,nx)
   Bm=matrix(0,TraitNum,nx)
   Bp=matrix(0,TraitNum,nx)


   VJR=matrix(0,nrow=ml,ncol=JacRep)
   JacVar=matrix(0,nrow=JacNum,ncol=ml)
   if(mk>0)JacPre=matrix(0,nrow=JacNum,ncol=mt)
   JacB=matrix(0,nrow=JacNum,ncol=nx)


   JV=matrix(0,nrow=ml,ncol=JacNum)
   if(mk>0)JP=matrix(0,nrow=mt,ncol=JacNum)
   #JB=matrix(0,nrow=nx,ncol=JacNum)

   df=JacNum-1

   JACVAR=matrix(0,TraitNum,ml)
   if(mk>0)JACPRE=matrix(0,TraitNum,mt)
   JACB=matrix(0,TraitNum,nx)

   JACKPRE=list()
   for(k in 1:JacRep){
      #k=1
      JB=sample(JB)
      U1=NULL
      X1=list()
      JacRes=NULL
      JAC1=NULL
      JAC2=NULL
      JAC3=NULL
      jackpre=list()

      for(i in 1:JacNum){
         #i=1
         index=which(JB==i)
         if(mk>0){
            for(j in 1:mk){
              m=gdata0$U[[j]]
              m=m[-index,]
              U1[[j]]=m
            }
         }
         gdata$U=U1
         names(gdata$U)=names(gdata0$U)
         for(j in 1:xk){
           m=gdata0$X[[j]]
           m=as.matrix(m[-index,])
           X1[[j]]=m
           colnames(X1[[j]])=colnames(gdata0$X[[j]])

         }
         gdata$X=X1
         #colnames(gdata$X[[1]])
         YD=as.matrix(Y[-index,])
         colnames(YD)=colnames(Y)
         gdata$Y=YD
         res=genmod(gdata,au=au,LUP=LUP)
         JACVAR=t(res$Var)
         JACB=t(res$FixEffect)
         if(mk>0){
            JACPRE=t(res$RandomEffect)
            jackpre[[i]]=JACPRE
         }
         JAC1[[i]]=as.matrix(JACVAR)
         if(mk>0)JAC2[[i]]=as.matrix(JACPRE)
         #if(nx==1)JACB
         JAC3[[i]]=as.matrix(JACB)


      }
      if(mk>0)JACKPRE[[k]]=jackpre
      for(s in 1:TraitNum){
          #s=1
          for(j in 1:JacNum){
             if(mk==0&&TraitNum==1)JacVar[j,]=JAC1[[j]]
             else JacVar[j,]=JAC1[[j]][s,]
             if(mk>0)JacPre[j,]=JAC2[[j]][s,]
             if(nx==1&&TraitNum==1)JacB[j,]=JAC3[[j]]
             else JacB[j,]=JAC3[[j]][s,]

          }
          for(j in 1:JacNum){
             index=which(JacVar[j,]<0)
             JacVar[j,index]=0
          }
          ####Calculate the variance components and their proportions from the jackknife estimates
          vt=as.vector(apply(JacVar,1,sum))
          for(i in 1:ml){
             m=mean(JacVar[,i])
             se=sqrt(var(JacVar[,i]))
             if(m<=0){
               m=0
               t1=0
               p=1.0
             }
             else if(m>0){
               t1=m/se
               p=pt(-abs(t1),df)
               p=1-(1-p)^(JacNum-2)
             }
             Vp[s,i]=Vp[s,i]+p
             Vm[s,i]=Vm[s,i]+m
             Cv[s,i]=Cv[s,i]+se
             m=mean(JacVar[,i]/vt)
             se=sqrt(var(JacVar[,i]/vt))
             if(m<=0){
               m=0
               t1=0
               p=1.0
             }
             else if(m>0){
               t1=m/se
               p=pt(-abs(t1),df)
               p=1-(1-p)^(JacNum-2)
             }
             VCp[s,i]=VCp[s,i]+p
             VCm[s,i]=VCm[s,i]+m
             VCv[s,i]=VCv[s,i]+se

          }
          #if(is.null(SIMUYES)){
          if(mk>0){
             for(i in 1:mt){
               m=mean(JacPre[,i])
               se=sqrt(var(JacPre[,i]))
               t1=m/se

               p=pt(-abs(t1),df)
               if(se<=0.0000001)p=1.0
               #p=a$p.value
               p=1-(1-p)^(JacNum-2)
               Pp[s,i]=Pp[s,i]+p
               Pm[s,i]=Pm[s,i]+m
               Cp[s,i]=Cp[s,i]+se
             }
          }
          #}
          for(i in 1:nx){
             m=mean(JacB[,i])
             se=sqrt(var(JacB[,i]))
             t1=m/se
             p=pt(-abs(t1),df)
             if(se<=0.0000001)p=1.0
             #p=a$p.value
             p=1-(1-p)^(JacNum-2)
             Bp[s,i]=Bp[s,i]+p
             Bm[s,i]=Bm[s,i]+m
             Cb[s,i]=Cb[s,i]+se
          }

      }
   }
   ### Variance components
   Vm=Vm/JacRep
   Vp=Vp/JacRep
   Cv=Cv/JacRep
   Cv1=Vm-t*Cv
   Cv2=Vm+t*Cv

   ### Proportional variance components
   VCm=VCm/JacRep
   VCp=VCp/JacRep
   VCv=VCv/JacRep
   VCv1=VCm-t*VCv
   VCv2=VCm+t*VCv
   #if(is.null(SIMUYES)){
   if(mk>0){
     ##Random effects
     Pm=Pm/JacRep
     Pp=Pp/JacRep
     Cp=Cp/JacRep
     Cp1=Pm-t*Cp
     Cp2=Pm+t*Cp
   }
   #}
   ##Fixed effwaects
   Bm=Bm/JacRep
   Bp=Bp/JacRep
   Cb=Cb/JacRep
   Cb1=Bm-t*Cb
   Cb2=Bm+t*Cb

   ##

   VAR=NULL
   PRE=NULL
   VC=NULL
   B=NULL
   CL=c(ALPHA/2,1-ALPHA/2)*100
   CL2=c("%LL","%UL")
   CL=paste(CL,CL2,sep="")
   CNames=c("Estimate","SE","PValue",CL)

   #CNames=c("Estimate","SE","PValue","2.5%CL","97.5%CU")
   for(i in 1: TraitNum){
      VAR[[i]]=data.frame(Vm[i,],Cv[i,],Vp[i,],Cv1[i,],Cv2[i,])
      VC[[i]]=data.frame(VCm[i,],VCv[i,],VCp[i,],VCv1[i,],VCv2[i,])
      if(mk>0)PRE[[i]]=data.frame(Pm[i,],Cp[i,],Pp[i,],Cp1[i,],Cp2[i,])
      B[[i]]=data.frame(Bm[i,],Cb[i,],Bp[i,],Cb1[i,],Cb2[i,])

      rownames(VAR[[i]])=VNames
      colnames(VAR[[i]])=CNames
      rownames(VC[[i]])=paste(VNames,"/VP",sep="")
      colnames(VC[[i]])=CNames
      #if(is.null(SIMUYES)){
        if(mk>0){
           rownames(PRE[[i]])=Str
           Rnames=CNames
           Rnames[1]="Pre"
           colnames(PRE[[i]])=Rnames

        }
      #}
      colnames(B[[i]])=CNames
      rownames(B[[i]])=XNames

   }

   names(VAR)=TraitNames
   names(VC)=TraitNames
   if(mk>0)names(PRE)=TraitNames
   names(B)=TraitNames
   ##a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
   #if(is.null(SIMUYES))a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,JackPre=JACKPRE)##,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
   #else a=list(Var=VAR,PVar=VC,FixedEffect=B)
   if(is.null(LUP))Prediction="AUP"
   else Prediction="LUP"
   #a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,JackPre=JACKPRE,Prediction=Prediction)
   a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Prediction=Prediction)

   #class(a)="jackknife.result"
   #a$call=match.call()
   #a$Var
   return(a)
}





lmm.mq <-
  function(formula,data=list()){
    if(is.null(data))gdata=mixed.data(formula)
    else gdata=mixed.data(formula,data)
    res=mq0(gdata)
    #res$FixedEffect=res$FixEffect
    return(res)
  }

lmm.reml <-
  function(formula,data=list(),criterion=NULL){
    if(is.null(criterion))criterion=1e-3
    if(is.null(data))gdata=mixed.data(formula)
    else gdata=mixed.data(formula,data)
    res=reml0(gdata,criterion)
    #res$FixedEffect=res$FixEffect
    return(res)
  }


lmm.mq.jack <-
  function(formula,data=list(),JacNum=NULL,JacRep=NULL){
    if(is.null(JacNum))JacNum=10
    if(is.null(JacRep))JacRep=1
    if(is.null(data))gdata=mixed.data(formula)
    else gdata=mixed.data(formula,data)
    res=genmod.jack(gdata,JacNum,JacRep)
    return(res)
  }
lmm.mq.perm <-
  function(formula,data=list(),PermNum=NULL){
    if(is.null(PermNum))PermNum=500
    ##if(is.null(JacRep))JacRep=1
    if(is.null(data))gdata=mixed.data(formula)
    else gdata=mixed.data(formula,data)
    res=genmod.perm(gdata=gdata,au=NULL,PermNum=PermNum)
    return(res)
  }

lmm.reml.jack <-
  function(formula,data=list(),criterion=NULL){
    if(is.null(criterion))criterion=1e-3
    if(is.null(data))gdata=mixed.data(formula)
    else gdata=mixed.data(formula,data)
    res=genmod.reml.jack(gdata,criterion=criterion)
    #res$FixedEffect=res$FixEffect
    return(res)
  }


##A function to generate a simulated data set
## Needed Dec 18, 2019
lmm.simudata=function(formula, data = list(),v,b,SimuNum=NULL){
  if (is.null(data))gdata = mixed.data(formula)
  else gdata = mixed.data(formula, data)
  if(is.null(SimuNum))SimuNum=50
  gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
  YS=gdata$Y
  return(YS)
}


#' Linear mixed model analysis with LMM approaches
#' @description n R function for linear mixed model analysis with REML and/or MINQUE approaches without resampling approaches.
#' @usage lmm(formula,data = list(), method = NULL, ALPHA = NULL)
#' @param formula	A user-defined linear mixed model formula used for data analysis
#' @param data Data frame. It is optional.
#' @param method The default linear mixed model approach is MINQUE. Users can choose both or one of two linear mixed model approaches, REML and MINQUE.
#' @param ALPHA a preset nominal probability level. The defualt value is 0.05.
#' @details A complete data set without missing points is needed. The model system includes
#' fixed effect and random effect components with separation of "|".
#' @return A list of lists of variance components,fixed effects, and random effects and their statistical tests.
#' @author Jixiang Wu <jixiang.wu@sdstate.edu>
#' @references Miller, R. G. 1974. The jackknife - a review. Biometrika, 61:1- 15.
#' @references Rao, C.R. 1971. Estimation of variance and covariance components-MINQUE theory. J Multiva Ana 1:19
#' @export
#'
lmm=function(formula,data = list(),method=NULL,ALPHA=NULL){
  if(is.null(method))method="minque"
  if(is.null(ALPHA))ALPHA=0.05
  if (is.null(data))data=NULL
  return(lmm1(formula,data=data,method=method,ALPHA=ALPHA))
}

lmm1=function(formula,data = list(),method=NULL,ALPHA=NULL){
  if(is.null(method))method="minque"
  if(is.null(ALPHA))ALPHA=0.05
  if (is.null(data))gdata = mixed.data(formula)
  else gdata = mixed.data(formula, data)
  if(method=="minque")RES=mq0(gdata)
  if(method=="reml")RES=reml0(gdata)
  mk=length(gdata$U)
  Uk=numeric()
  BigU=NULL
  if(mk>0){
    for(u in 1:mk){
      Uk[u]=ncol(gdata$U[[u]])
      BigU=cbind(BigU,gdata$U[[u]])
    }
  }
  mt=sum(Uk)
  xk=length(gdata$X)
  Xk=numeric()
  X=NULL
  for(i in 1:xk){
    Xk[i]=ncol(gdata$X[[i]])
    X=cbind(X,gdata$X[[i]])
  }

  TraitNum=ncol(gdata$Y)
  Residual=NULL
  Fitted=NULL
  for(j in 1:TraitNum){
    a=0
    re=RES$RandomEffect[[j]][,1]
    fe=RES$FixedEffect[[j]][,1]
    if(mk>0)a=a+BigU%*%re
    fit=a+X%*%fe
    res=gdata$Y[,j]-a
    res.mean=mean(res)
    fit=fit+res.mean
    res=res-res.mean
    Fitted=cbind(Fitted,fit)
    Residual=cbind(Residual,res)
  }
  colnames(Residual)=colnames(gdata$Y)
  colnames(Fitted)=colnames(gdata$Y)

  RES$Residual=Residual
  RES$Fitted=Fitted
  RES$method=method
  RES$ALPHA=ALPHA
  return(RES)
}


#' Linear mixed model analysis with LMM approaches
#' @description n R function for linear mixed model analysis with REML and/or MINQUE approaches with resampling approaches.
#' @usage lmm(formula,data = list(), method = NULL, ALPHA = NULL)
#' @param formula	A user-defined linear mixed model formula used for data analysis
#' @param data Data frame. It is optional.
#' @param JacNum A number of groups for jackknife. The default number is 10
#' @param JacReo A number of repeating jackkniffing. The default number is 1.Recommended number is 5.
#' @param method The default linear mixed model approach is MINQUE. Users can choose both or one of two linear mixed model approaches, REML and MINQUE.
#' @param ALPHA a preset nominal probability level. The defualt value is 0.05.
#' @details A complete data set without missing points is needed. The model system includes
#' fixed effect and random effect components with separation of "|".
#' @return A list of lists of variance components,fixed effects, and random effects and their statistical tests.
#' @author Jixiang Wu <jixiang.wu@sdstate.edu>
#' @references Miller, R. G. 1974. The jackknife - a review. Biometrika, 61:1- 15.
#' @references Rao, C.R. 1971. Estimation of variance and covariance components-MINQUE theory. J Multiva Ana 1:19
#' @export
#'

lmm.jack=function(formula,data = list(),method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
  if(is.null(method))method="minque"
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  if(is.null(ALPHA))ALPHA=0.05
  if(is.null(data))data=NULL
  res=lmm.jack1(formula,data=data,method=method,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
  return(res)
}

lmm.jack1=function(formula,data = list(),method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
  if(is.null(method))method="minque"
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  if(is.null(ALPHA))ALPHA=0.05
  if (is.null(data))gdata = mixed.data(formula)
  else gdata = mixed.data(formula, data)
  if(method=="minque")RES=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
  if(method=="reml")RES=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)

  mk=length(gdata$U)
  Uk=numeric()
  BigU=NULL
  if(mk>0){
    for(u in 1:mk){
      Uk[u]=ncol(gdata$U[[u]])
      BigU=cbind(BigU,gdata$U[[u]])
    }
  }
  mt=sum(Uk)

  xk=length(gdata$X)
  Xk=numeric()
  X=NULL
  for(i in 1:xk){
    Xk[i]=ncol(gdata$X[[i]])
    X=cbind(X,gdata$X[[i]])
  }

  TraitNum=ncol(gdata$Y)
  Residual=NULL
  Fitted=NULL
  for(j in 1:TraitNum){
    a=0
    re=RES$RandomEffect[[j]][,1]
    fe=RES$FixedEffect[[j]][,1]
    if(mk>0)a=a+BigU%*%re
    a=a+X%*%fe
    Fitted=cbind(Fitted,a)
    a=gdata$Y[,j]-a
    Residual=cbind(Residual,a)
  }
  colnames(Residual)=colnames(gdata$Y)
  colnames(Fitted)=colnames(gdata$Y)

  RES$residual=Residual
  RES$fitted=Fitted
  RES$method=method
  RES$ALPHA=ALPHA

  return(RES)
}

lmm.jack0=function(formula,lmlist=NULL){
  # data = list(),method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
  if(is.null(lmlist)){
    method=c("minque")
    JacNum=10
    JacRep=1
    ALPHA=0.05
    gdata = mixed.data(formula)
  }
  else{

    if(is.null(lmlist$method))method=c("minque")
    if(is.null(lmlist$JacNum))JacNum=10
    if(is.null(lmlist$JacRep))JacRep=1
    if(is.null(lmlist$ALPHA))ALPHA=0.05
    if (is.null(lmlist$data))gdata = mixed.data(formula)
    else gdata = mixed.data(formula, data=lmlist$data)
  }
  mn=length(method)
  RES=list()
  for(i in 1:mn){
    if(method[i]=="minque")RES[[i]]=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
    if(method[i]=="reml")RES[[i]]=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
  }
  mk=length(gdata$U)
  Uk=numeric()
  BigU=NULL
  for(u in 1:mk){
    Uk[u]=ncol(gdata$U[[u]])
    BigU=cbind(BigU,gdata$U[[u]])
  }
  mt=sum(Uk)
  xk=length(gdata$X)
  Xk=numeric()
  X=NULL
  for(i in 1:xk){
    Xk[i]=ncol(gdata$X[[i]])
    X=cbind(X,gdata$X[[i]])
  }
  TraitNum=ncol(gdata$Y)
  for(i in 1:mn){
    Residual=NULL
    Fitted=NULL
    for(j in 1:TraitNum){
      a=0
      re=RES[[i]]$RandomEffect[[j]][,1]
      fe=RES[[i]]$FixedEffect[[j]][,1]
      if(mk>0)a=a+BigU%*%re
      a=a+X%*%fe
      Fitted=cbind(Fitted,a)
      a=gdata$Y[,j]-a
      Residual=cbind(Residual,a)
    }
    colnames(Residual)=colnames(gdata$Y)
    colnames(Fitted)=colnames(gdata$Y)
    RES[[i]]$Residual=Residual
    RES[[i]]$Fitted=Fitted
  }

  names(RES)=method
  RES$ALPHA=ALPHA
  return(RES)
}

## Needed Dec 19
lmm.simu=function(formula,method=NULL,ALPHA=NULL){
  if(is.null(method))method="minque"
  if(is.null(ALPHA))ALPHA=0.05

  gdata = mixed.data(formula)

  if(method=="minque")RES=mq0(gdata)
  if(method=="reml")RES=reml0(gdata)

  #res=RES
  V=NULL
  Var=RES$Var
  SimuNum=length(Var)
  ml=length(Var[[1]][,1])
  P=numeric(ml)
  vm=numeric(ml)
  for(j in 1:SimuNum){
    vm=vm+Var[[j]][,1]
    V=cbind(V,Var[[j]][,1])
    id=which(Var[[j]][,4]<=ALPHA)
    P[id]=P[id]+1
  }
  SE=numeric(ml)
  for(j in 1:ml)SE[j]=sqrt(var(V[j,])/SimuNum)
  vm=vm/SimuNum
  P=P/SimuNum
  a=data.frame(vm,SE,P)
  colnames(a)=c("Estimate","SE","Power")

  rownames(a)=rownames(Var[[1]])

  res=list(Simu=a,method=method,alpha=ALPHA)
  #RES$Simu=a

  #RES$method=method
  #RES$ALPHA=ALPHA

  return(res)
}

## Needed Dec 2019
lmm.simu.jack=function(formula,method=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
  if(is.null(method))method="minque"
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  if(is.null(ALPHA))ALPHA=0.05

  gdata = mixed.data(formula)
  if(method=="minque")RES=genmod.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)
  if(method=="reml")RES=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)

  #res=RES
  V=NULL
  Var=RES$Var
  SimuNum=length(Var)
  ml=length(Var[[1]][,1])
  P=numeric(ml)
  vm=numeric(ml)
  for(j in 1:SimuNum){
    vm=vm+Var[[j]][,1]
    V=cbind(V,Var[[j]][,1])
    id=which(Var[[j]][,3]<=ALPHA)
    P[id]=P[id]+1
  }
  SE=numeric(ml)
  for(j in 1:ml)SE[j]=sqrt(var(V[j,])/SimuNum)
  vm=vm/SimuNum
  P=P/SimuNum
  a=data.frame(vm,SE,P)
  colnames(a)=c("Estimate","SE","Power")
  rownames(a)=rownames(Var[[1]])

  res=list(Simu=a,method=method,alpha=ALPHA)
  return(res)
}
lmm.reml.simu=function(formula, data = list(),v0,b0,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL){
  if (is.null(data))gdata = mixed.data(formula)
  else gdata = mixed.data(formula, data)
  if(is.null(SimuNum))SimuNum=200
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  if(is.null(ALPHA))ALPHA=0.05

  result=genmod.reml.simu(gdata,v=v0,b=b0,SimuNum=SimuNum,JacNum=JacNum,JacRep=JacRep,ALPHA=ALPHA)

  return(result)
}

genmod.reml.simu=function(gdata,v,b,SimuNum=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL,...){
  if(is.null(SimuNum))SimuNum=200
  if(is.null(ALPHA))ALPHA=0.05
  gdata=SimuData(gdata,v=v,b=b,SimuNum=SimuNum)
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1

  jac=genmod.reml.jack(gdata,JacNum=JacNum,JacRep=JacRep,LUP=1)
  ml=length(v)
  Var=jac$Var
  P=numeric(ml)
  vm=numeric(ml)
  for(i in 1:SimuNum){
    vm=vm+Var[[i]][,1]
    id=which(Var[[i]][,3]<=ALPHA)
    P[id]=P[id]+1
  }
  vm=vm/SimuNum
  P=P/SimuNum
  bias=vm-v
  a=data.frame(v,vm,bias,P)

  colnames(a)=c("True","Estimate","Bias","Power")
  rownames(a)=rownames(Var[[1]])
  res=list(P=a,ALPHA=ALPHA)

  return(res)
}






get_HR=function(gdata,Qa,v){
  mk=length(gdata$U)
  h=NULL
  dim(t(gdata$U[[1]]))
  n=ncol(Qa)
  #V=Diagonal(n)*v[mk+1]
  ##i=1
  #for(i in 1:mk){
  #    vv=gdata$U[[i]]%*%t(gdata$U[[i]])*v[i]
  #    V=V+vv
  #}
  #Q=Qa%*%V%*%Qa
  Q=Qa

  for(i in 1:mk){
    m=t(gdata$U[[i]])%*%Q%*%gdata$U[[i]]
    m=m*(v[i])^2
    #m=svdinv(m)
    if(i==1)h=diag(m)
    else h=c(h,diag(m))
  }
  summary(h)
  h=sqrt(h)
  return(h)
}
genmod.mq=function(gdata){
  res0=genmod(gdata)
  res=genmod(gdata,au=res0$Var[,1])
  ##Statistical tests for variance components
  h=diag(svdinv(res$HV))*2
  Chi_sq=(res0$Var[,1])^2/h
  P_value=(1-pchisq(Chi_sq,1))/2
  Var=data.frame(res0$Var,sqrt(h),Chi_sq,P_value)
  colnames(Var)=c("Est","SE","Chi_sq","P_value")
  res$Var=Var

  ##Statistical tests for fixed effects
  #res0$FixEffect[,1]
  h=diag(svdinv(res$HF[[1]]))
  h=sqrt(h)
  z=res0$FixEffect[,1]/h
  P_value=1-pchisq(z^2,1)
  attributes(res)
  FixEffect=data.frame(res0$FixEffect,h,z,P_value)
  colnames(FixEffect)=c("Est","SE","z_value","P_value")
  res$FixEffect=FixEffect

  ##Statistical tests for random effect
  mk=length(gdata$U)
  if(mk>0){

    h=get_HR(gdata,res$Qa,res0$Var[,1])
    z=res0$RandomEffect[,1]/h
    res0$RandomEffect[,1][1]/h[1]

    P_value=1-pchisq(z^2,1)
    RandomEffect=data.frame(res0$RandomEffect,h,z,P_value)
    colnames(RandomEffect)=c("Pre","SE","z_value","P_value")
    res$RandomEffect=RandomEffect
  }

  return(res)
}

genmod.reml=function(gdata,criterion=NULL){
  #svdinv=ginv
  if(is.null(criterion))DIFF=1e-3
  else DIFF=criterion
  d=100
  res=genmod(gdata)
  mk=length(gdata$U)
  au0=res$Var[,1]
  id=which(au0<=0)
  au0[id]=0
  its=1
  while(d>DIFF&&its<10){

    res=genmod(gdata,au=au0)
    (v1=res$Var[,1])
    id=which(v1<=0)
    v1[id]=0
    #print(v1)
    #cat("\n")

    d=sum(au0-v1)^2
    d=sqrt(d)
    d=d/sum(au0)
    au0=v1
    #print(d)
    #cat("\n")
    its=its+1
  }
  ##Statistical tests for variance components
  h=diag(svdinv(res$HV))*2
  Chi_sq=(res$Var[,1])^2/h
  P_value=(1-pchisq(Chi_sq,1))/2
  Var=data.frame(res$Var,sqrt(h),Chi_sq,P_value)
  colnames(Var)=c("Est","SE","Chi_sq","P_value")
  res$Var=Var

  ##Statistical tests for fixed effects
  res$FixEffect[,1]
  #h=diag(svdinv(res$HF[[1]]))
  h=diag(ginv(res$HF[[1]]))
  h=sqrt(h)
  z=res$FixEffect[,1]/h
  P_value=1-pchisq(z^2,1)
  attributes(res)
  FixEffect=data.frame(res$FixEffect,h,z,P_value)
  colnames(FixEffect)=c("Est","SE","z_value","P_value")
  res$FixEffect=FixEffect

  ##Statistical tests for random effect
  mk=length(gdata$U)
  if(mk>0){
    h=get_HR(gdata,res$Qa,res$Var[,1])
    z=res$RandomEffect[,1]/h
    res$RandomEffect[,1][1]/h[1]
    P_value=1-pchisq(z^2,1)
    RandomEffect=data.frame(res$RandomEffect,h,z,P_value)
    colnames(RandomEffect)=c("Pre","SE","z_value","P_value")
    res$RandomEffect=RandomEffect
  }

  return(res)
}
##This function is used without jackknife process
reml0=function(gdata,criterion=NULL){
  if(is.null(criterion))criterion=1e-3
  Y0=gdata$Y
  TraitNum=ncol(Y0)
  TraitNames=colnames(Y0)
  Var=list()
  RE=NULL
  FE=list()
  mk=length(gdata$U)
  for(i in 1:TraitNum){
    gdata$Y=as.matrix(Y0[,i])
    res=genmod.reml(gdata,criterion)
    Var[[i]]=res$Var
    if(mk>0)RE[[i]]=res$RandomEffect
    #else RE[[i]]=NULL
    FE[[i]]=res$FixEffect
  }

  names(Var)=TraitNames
  if(!is.null(RE))names(RE)=TraitNames
  names(FE)=TraitNames
  res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
  gdata$Y=Y0
  return(res)
}


##This function is used for jackknife process
reml=function(gdata,criterion=NULL){
  if(is.null(criterion))criterion=1e-3
  Y0=as.matrix(gdata$Y)
  TraitNum=ncol(Y0)
  TraitNames=colnames(Y0)
  Var=NULL
  RE=NULL
  FE=NULL
  mk=length(gdata$U)
  for(i in 1:TraitNum){
    gdata$Y=as.matrix(Y0[,i])
    res=genmod.reml(gdata,criterion)
    Var=cbind(Var,res$Var[,1])
    if(mk>0)RE=cbind(RE,res$RandomEffect[,1])
    FE=cbind(FE,res$FixEffect[,1])
  }
  res$FixEffect
  Var=as.matrix(Var)
  colnames(Var)=TraitNames
  rownames(Var)=rownames(res$Var)
  if(mk>0){
    RE=as.matrix(RE)
    colnames(RE)=TraitNames
    rownames(RE)=rownames(res$RandomEffect)
  }
  FE=as.matrix(FE)
  colnames(FE)=TraitNames
  rownames(FE)=rownames(res$FixEffect)
  res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
  gdata$Y=Y0
  return(res)
}

#A=NULL
#A[[1]]=c(1:2)
##This MINQUE function is used without jackknife process
mq0=function(gdata){
  Y0=as.matrix(gdata$Y)
  #Y0=as.matrix(Y0)
  TraitNum=ncol(Y0)
  TraitNames=colnames(Y0)
  Var=list()
  RE=NULL
  FE=list()
  #i=1
  mk=length(gdata$U)
  for(i in 1:TraitNum){
    gdata$Y=as.matrix(Y0[,i])
    res=genmod.mq(gdata)
    Var[[i]]=res$Var
    if(mk>0)RE[[i]]=res$RandomEffect
    #else RE[[i]]=NULL
    FE[[i]]=res$FixEffect
  }
  #res$FixEffect
  names(Var)=TraitNames
  if(!is.null(RE))names(RE)=TraitNames
  names(FE)=TraitNames
  res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
  gdata$Y=Y0
  return(res)
}


##This MINQUE function is used with jackknife
mq=function(gdata){
  Y0=gdata$Y
  TraitNum=ncol(Y0)
  TraitNames=colnames(Y0)
  Var=NULL
  RE=NULL
  FE=NULL
  for(i in 1:TraitNum){
    gdata$Y=as.matrix(Y0[,i])
    res=genmod(gdata)
    Var=cbind(Var,res$Var[,1])
    RE=cbind(RE,res$RandomEffect[,1])
    FE=cbind(FE,res$FixEffect[,1])
  }
  res$FixEffect
  Var=as.matrix(Var)
  colnames(Var)=TraitNames
  rownames(Var)=rownames(res$Var)
  RE=as.matrix(RE)
  colnames(RE)=TraitNames
  rownames(RE)=rownames(res$RandomEffect)
  FE=as.matrix(FE)
  colnames(FE)=TraitNames
  rownames(FE)=rownames(res$FixEffect)
  res=list(Var=Var,FixedEffect=FE,RandomEffect=RE)
  gdata$Y=Y0
  return(res)
}


genmod.reml.jack=function(gdata,criterion=NULL,JacNum=NULL,JacRep=NULL,ALPHA=NULL,au=NULL,LUP=NULL){
  if(is.null(criterion))criterion=1e-3
  if(is.null(JacNum))JacNum=10
  if(is.null(JacRep))JacRep=1
  if(is.null(ALPHA))ALPHA=0.05
  if(is.null(au))au=NULL
  if(is.null(LUP))LUP=NULL
  a=1-(1-ALPHA)^(1/(JacNum-2))
  t=qt(1-a/2,JacNum-1)
  Y=as.matrix(gdata$Y)
  TraitNum=ncol(Y)
  TraitNames=colnames(Y)
  n=nrow(Y)
  #colnames(Y)
  mk=length(gdata$U)
  xk=length(gdata$X)
  JB=rep(1:JacNum,length=n)
  gdata0=gdata
  ml=mk+1
  nx=0
  for(i in 1:xk)nx=nx+ncol(gdata0$X[[i]])
  res0=reml(gdata0,criterion)
  V0=res0$Var
  B0=res0$FixedEffect
  P0=res0$RandomEffect
  if(mk>0){
    mt=nrow(P0)
    Str=rownames(P0)
  }
  VNames=rownames(V0)
  XNames=rownames(B0)

  ##############################################################################
  RES0=NULL
  VPower=numeric(ml)
  V0=matrix(0,TraitNum,ml)
  if(mk>0)P0=matrix(0,TraitNum,mt)
  #colnames(P0)
  Vp=matrix(0,TraitNum,ml)
  Cv=matrix(0,TraitNum,ml)
  Cv1=matrix(0,TraitNum,ml)
  Cv2=matrix(0,TraitNum,ml)
  Vm=matrix(0,TraitNum,ml)

  VCp=matrix(0,TraitNum,ml)
  VCv=matrix(0,TraitNum,ml)
  VCv1=matrix(0,TraitNum,ml)
  VCv2=matrix(0,TraitNum,ml)
  VCm=matrix(0,TraitNum,ml)

  if(mk>0){
    Pp=matrix(0,TraitNum,mt)
    Cp=matrix(0,TraitNum,mt)
    Cp1=matrix(0,TraitNum,mt)
    Cp2=matrix(0,TraitNum,mt)
    Pm=matrix(0,TraitNum,mt)
  }

  Cb1=matrix(0,TraitNum,nx)
  Cb2=matrix(0,TraitNum,nx)
  Cb=matrix(0,TraitNum,nx)
  Bm=matrix(0,TraitNum,nx)
  Bp=matrix(0,TraitNum,nx)


  VJR=matrix(0,nrow=ml,ncol=JacRep)
  JacVar=matrix(0,nrow=JacNum,ncol=ml)
  if(mk>0)JacPre=matrix(0,nrow=JacNum,ncol=mt)
  JacB=matrix(0,nrow=JacNum,ncol=nx)


  JV=matrix(0,nrow=ml,ncol=JacNum)
  if(mk>0)JP=matrix(0,nrow=mt,ncol=JacNum)
  #JB=matrix(0,nrow=nx,ncol=JacNum)

  df=JacNum-1

  JACVAR=matrix(0,TraitNum,ml)
  if(mk>0)JACPRE=matrix(0,TraitNum,mt)
  JACB=matrix(0,TraitNum,nx)

  JACKPRE=list()
  for(k in 1:JacRep){
    #k=1
    JB=sample(JB)
    U1=NULL
    X1=list()
    JacRes=NULL
    JAC1=NULL
    JAC2=NULL
    JAC3=NULL
    jackpre=list()

    for(i in 1:JacNum){
      #i=1
      index=which(JB==i)
      if(mk>0){
        for(j in 1:mk){
          m=gdata0$U[[j]]
          m=m[-index,]
          U1[[j]]=m
        }
      }
      gdata$U=U1
      names(gdata$U)=names(gdata0$U)
      for(j in 1:xk){
        m=gdata0$X[[j]]
        m=as.matrix(m[-index,])
        X1[[j]]=m
        colnames(X1[[j]])=colnames(gdata0$X[[j]])

      }
      gdata$X=X1
      #colnames(gdata$X[[1]])
      YD=as.matrix(Y[-index,])
      colnames(YD)=colnames(Y)
      gdata$Y=YD
      res=reml(gdata,criterion)
      JACVAR=t(res$Var)
      JACB=t(res$FixedEffect)
      if(mk>0){
        JACPRE=t(res$RandomEffect)
        jackpre[[i]]=JACPRE
      }
      JAC1[[i]]=as.matrix(JACVAR)
      if(mk>0)JAC2[[i]]=as.matrix(JACPRE)
      #if(nx==1)JACB
      JAC3[[i]]=as.matrix(JACB)


    }
    if(mk>0)JACKPRE[[k]]=jackpre
    for(s in 1:TraitNum){
      #s=1
      for(j in 1:JacNum){
        if(mk==0&&TraitNum==1)JacVar[j,]=JAC1[[j]]
        else JacVar[j,]=JAC1[[j]][s,]
        if(mk>0)JacPre[j,]=JAC2[[j]][s,]
        if(nx==1&&TraitNum==1)JacB[j,]=JAC3[[j]]
        else JacB[j,]=JAC3[[j]][s,]

      }
      for(j in 1:JacNum){
        index=which(JacVar[j,]<0)
        JacVar[j,index]=0
      }
      ####Calculate the variance components and their proportions from the jackknife estimates
      vt=as.vector(apply(JacVar,1,sum))
      for(i in 1:ml){
        m=mean(JacVar[,i])
        se=sqrt(var(JacVar[,i]))
        if(m<=0){
          m=0
          t1=0
          p=1.0
        }
        else if(m>0){
          t1=m/se
          p=pt(-abs(t1),df)
          p=1-(1-p)^(JacNum-2)
        }
        Vp[s,i]=Vp[s,i]+p
        Vm[s,i]=Vm[s,i]+m
        Cv[s,i]=Cv[s,i]+se
        m=mean(JacVar[,i]/vt)
        se=sqrt(var(JacVar[,i]/vt))
        if(m<=0){
          m=0
          t1=0
          p=1.0
        }
        else if(m>0){
          t1=m/se
          p=pt(-abs(t1),df)
          p=1-(1-p)^(JacNum-2)
        }
        VCp[s,i]=VCp[s,i]+p
        VCm[s,i]=VCm[s,i]+m
        VCv[s,i]=VCv[s,i]+se

      }
      #if(is.null(SIMUYES)){
      if(mk>0){
        for(i in 1:mt){
          m=mean(JacPre[,i])
          se=sqrt(var(JacPre[,i]))
          t1=m/se

          p=pt(-abs(t1),df)
          if(se<=0.0000001)p=1.0
          #p=a$p.value
          p=1-(1-p)^(JacNum-2)
          Pp[s,i]=Pp[s,i]+p
          Pm[s,i]=Pm[s,i]+m
          Cp[s,i]=Cp[s,i]+se
        }
      }
      #}
      for(i in 1:nx){
        m=mean(JacB[,i])
        se=sqrt(var(JacB[,i]))
        t1=m/se
        p=pt(-abs(t1),df)
        if(se<=0.0000001)p=1.0
        #p=a$p.value
        p=1-(1-p)^(JacNum-2)
        Bp[s,i]=Bp[s,i]+p
        Bm[s,i]=Bm[s,i]+m
        Cb[s,i]=Cb[s,i]+se
      }

    }
  }
  ### Variance components
  Vm=Vm/JacRep
  Vp=Vp/JacRep
  Cv=Cv/JacRep
  Cv1=Vm-t*Cv
  Cv2=Vm+t*Cv

  ### Proportional variance components
  VCm=VCm/JacRep
  VCp=VCp/JacRep
  VCv=VCv/JacRep
  VCv1=VCm-t*VCv
  VCv2=VCm+t*VCv
  #if(is.null(SIMUYES)){
  if(mk>0){
    ##Random effects
    Pm=Pm/JacRep
    Pp=Pp/JacRep
    Cp=Cp/JacRep
    Cp1=Pm-t*Cp
    Cp2=Pm+t*Cp
  }
  #}
  ##Fixed effwaects
  Bm=Bm/JacRep
  Bp=Bp/JacRep
  Cb=Cb/JacRep
  Cb1=Bm-t*Cb
  Cb2=Bm+t*Cb

  ##

  VAR=NULL
  PRE=NULL
  VC=NULL
  B=NULL
  CL=c(ALPHA/2,1-ALPHA/2)*100
  CL2=c("%LL","%UL")
  CL=paste(CL,CL2,sep="")
  CNames=c("Estimate","SE","PValue",CL)

  #CNames=c("Estimate","SE","PValue","2.5%CL","97.5%CU")
  for(i in 1: TraitNum){
    VAR[[i]]=data.frame(Vm[i,],Cv[i,],Vp[i,],Cv1[i,],Cv2[i,])
    VC[[i]]=data.frame(VCm[i,],VCv[i,],VCp[i,],VCv1[i,],VCv2[i,])
    if(mk>0)PRE[[i]]=data.frame(Pm[i,],Cp[i,],Pp[i,],Cp1[i,],Cp2[i,])
    B[[i]]=data.frame(Bm[i,],Cb[i,],Bp[i,],Cb1[i,],Cb2[i,])

    rownames(VAR[[i]])=VNames
    colnames(VAR[[i]])=CNames
    rownames(VC[[i]])=paste(VNames,"/VP",sep="")
    colnames(VC[[i]])=CNames
    #if(is.null(SIMUYES)){
    if(mk>0){
      rownames(PRE[[i]])=Str
      Rnames=CNames
      Rnames[1]="Pre"
      colnames(PRE[[i]])=Rnames
    }
    #}
    colnames(B[[i]])=CNames
    rownames(B[[i]])=XNames

  }

  names(VAR)=TraitNames
  names(VC)=TraitNames
  if(mk>0)names(PRE)=TraitNames
  names(B)=TraitNames
  ##a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
  #if(is.null(SIMUYES))a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,JackPre=JACKPRE)##,Method=Method,JacNum=JacNum,JacRep=JacRep,VNames=VNames,Y=Y,X=X,U=U)
  #else a=list(Var=VAR,PVar=VC,FixedEffect=B)
  if(is.null(LUP))Prediction="AUP"
  else Prediction="LUP"
  #a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,JackPre=JACKPRE,Prediction=Prediction)
  #a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B,Prediction=Prediction)
  a=list(Var=VAR,PVar=VC,RandomEffect=PRE,FixedEffect=B)
  #class(a)="jackknife.result"
  #a$call=match.call()
  #a$Var
  return(a)
}

lmm.mq.perm <-
  function(formula,data=list(),PermNum=NULL){
    if(is.null(PermNum))PermNum=500
    if(is.null(data))gdata=mixed.data(formula)
    else gdata=mixed.data(formula,data)
    res=genmod.perm(gdata=gdata,au=NULL,PermNum=PermNum)
    return(res)
  }



## needed Dec 18, 2019
lmm.perm <-function(formula,data=list(),method=NULL,PermNum=NULL){
  if(is.null(PermNum))PermNum=100
  if(is.null(method))method="minque"
  if(is.null(data))gdata=mixed.data(formula)
  else gdata=mixed.data(formula,data)
  #mn=length(method)
  #RES=list()
  #for(i in 1:mn){
  if(method=="minque")RES=genmod.perm(gdata,au=NULL,PermNum=PermNum)
  if(method=="reml")RES=reml.perm(gdata,au=NULL,PermNum=PermNum)
  RES$method=method
  RES$PermNum=PermNum

  return(RES)
}



## what for?
lmm.ga=function(formula,data=list()){
  if(is.null(data))gdata=mixed.data(formula)
  else gdata=mixed.data(formula=formula,data=data)
  res=genmod.mq(gdata)
}



###A function to generate a list of information matrices #######
## needed dec 18, 2019
GenerateU=function(mat){
  names=colnames(mat)
  r=length(names)
  U=list(r)
  for(i in 1:r){
    U[[i]]=IndMatrix(mat[,i])
    colnames(U[[i]])=paste(names[i],"(",colnames(U[[i]]),")",sep="")
  }
  names(U)=names
  class(U)="Umatrix"
  return(U)
}

## needed dec 18, 2019
########################################################
## A very nice function to generate a indicator matrix##
IndMatrix <- function(v)
{
  nr=length(v)
  #stop("Sample size is too large, please contact the author <qgtools@gmail.com> for additional assistance")
  lv=sort(unique(v))
  #lv=unique(v)
  nc=length(lv)
  m=matrix(0,nr,nc)
  for(i in 1:nr){
    j=which(v[i]==lv)
    m[i,j]=1
  }
  colnames(m)=lv
  return(m)


}

##----------------------------------------------------##
## Commentted on March 6, 2017 -----------------------##
## A key function to generate Umatrix with a linear --##
## mixed model frame ---------------------------------##
## mixed-formula:  separate which effects are ranodm -##
## and which are fixed -------------------------------##

## Needed Dec 18, 2019

lm.data=function(formula,data=list(),...){
  a=mixed_formula(formula)
  model=a$model
  yv=a$yv
  mod1=a$mod1
  if(model==2)mod2=a$mod2
  if(is.null(data))a1=model_data(mod1,yv)
  else a1=model_data(mod1,yv,data)
  df1=a1$df
  y=a1$y
  data1=a1$data
  name1=a1$name

  if(model==2){
    a2=model_data(mod2,yv,data)
    df2=a2$df
    data2=a2$data
    name2=a2$name
  }
  if(model==1)df=df1
  if(model==2){
    e=1
    for(i in 1:length(name1)){
      #i=1
      id=which(name2==name1[i])
      if(length(id)>0){
        if(e==1)du_id=id
        else du_id=c(du_id,id)
        e=e+1
      }
    }
    if(e>1)name2=name2[-du_id]
    dfname1=colnames(df1)
    dfname2=colnames(df2)
    e=1
    for(i in 1:length(dfname1)){
      id=which(dfname2==dfname1[i])
      if(length(id)>0){
        if(e==1)id0=id
        else id0=c(id0,id)
        e=e+1
      }
    }
    df=data.frame(df1,df2[,-id0])
    colnames(df)=c(dfname1,dfname2[-id0])
    #head(df)
  }

  for(i in 1:ncol(df)){
    if(i==1)cls=class(df[,i])
    else cls=c(cls,class(df[,i]))
  }
  name0=colnames(df)
  if(length(name1)==1){
    dat=matrix(1,nrow(df),1)
    #colnames(dat)="mu"
    cls1="numeric"
  }
  if(length(name1)>1){
    d1=fixed_data(name0,name1,cls,df)
    dat1=d1$dat
    cls1=c("numeric",d1$cls)
    dat=data.frame(1,dat1)
    uname1=d1$uname
  }
  if(model==2){
    d2=random_data(name0,name2,cls,df)
    dat2=d2$dat
    cls2=d2$cls
    uname2=d2$uname
    dat=data.frame(dat,dat2)
  }

  #head(dat)
  #if(length(name1)>0)
  #dat=cbind(1,dat)
  dat=data.frame(dat)
  if(model==1){
    if(length(name1)==1)colnames(dat)="mu"
    else colnames(dat)=c("mu",uname1)
  }
  else if(model==2){
    if(length(name1)==1)colnames(dat)=c("mu",uname2)
    else if(length(name1)>1)colnames(dat)=c("mu",uname1,uname2)
  }

  for(i in 2:ncol(df)){
    if(cls[i-1]=="factor")dat[,i]=factor(dat[,i])
  }
  if(model==1)cls0=cls1
  else cls0=c(cls1,cls2)

  fid=length(cls1)
  id=which(cls0=="numeric")
  if(fid==1){
    X1=matrix(1,nrow(dat),1)
    colnames(X1)=colnames(dat)[1]
  }
  else{
    X0=dat[,1:fid]
    if(length(id)==1){
      X1=matrix(1,nrow(dat),1)
      colnames(X1)=colnames(dat)[1]
    }
    else X1=X0[,id]
  }
  X2=NULL
  if(fid>ncol(X1)){
    xmat=X0[,-id]
    if(ncol(X0)==(length(id)+1)){
      xmat=data.frame(xmat)
      colnames(xmat)=colnames(X0)[-id]
    }
    X2=GenerateU(xmat)
  }

  X=list()
  X4=matrix(0,nrow(X1),ncol(X1))
  for(i in 1:ncol(X1))X4[,i]=as.numeric(X1[,i])
  colnames(X4)=colnames(X1)
  X[[1]]=X4
  xk=length(X2)
  if(xk>0){
    for(i in 1:xk)X[[i+1]]=X2[[i]]
  }
  class(X)="Umatrix"
  names(X)=c("B0",names(X2))
  if(model==2){
    rmat=dat[,-(1:fid)]
    if(ncol(dat)==(fid+1)){
      rmat=data.frame(rmat)
      colnames(rmat)=colnames(dat)[fid+1]
    }
    U=GenerateU(rmat)
    class(U)="Umatrix"
  }
  #colnames(U[[1]])
  else U=NULL
  TraitNames=colnames(y)
  result=list(Y=y,dat=dat,X=X,U=U,class=cls0,TraitNames=TraitNames)
  return(result)
}

## Needed Dec 18, 2019
mixed_formula=function(formula,...){
  s1=toString(formula)
  s1=strsplit(s1,",")[[1]]
  mc=length(s1)

  s2=strsplit(s1[mc],"\\+")[[1]]
  id=grep("\\|",s2)
  if(length(id)==0)model=1
  else if(length(id)>0) model=2

  ##to get two model equations
  if(model==1)mod1=s1[mc]
  if(model==2){
    fullmod=strsplit(s1[mc],"\\|")[[1]]
    mod1=fullmod[1]
    mod2=fullmod[2]
  }

  ##to get names of Y used for model statement
  s=toString(formula)
  s=gsub("\\~,","",s)
  s=gsub(" ","",s)
  id1=grep("),",s)
  if(length(id1)==0){
    s2=strsplit(s,",")[[1]]
    yv=s2[1]
  }
  if(length(id1)>0){
    s2=strsplit(s,"),")[[1]]
    #length(s2)
    s2[1]=paste(s2[1],")",sep="")
    yv=s2[1]
  }
  if(model==1)res=list(mod1=mod1,yv=yv,model=model)
  else res=list(mod1=mod1,mod2=mod2,yv=yv,model=model)
  return(res)
}
#mod=mod1

## Needed Dec 18, 2019
model_data=function(mod,yv,data=NULL,...){
  #mod1=gsub(" ","",mod)
  #yv="Y"
  names=c(yv,mod)
  fm1=paste(names, collapse="~")
  fm1=as.formula(fm1)
  if(is.null(data))data=NULL
  if(is.null(data))df1=model.frame(fm1)
  else df1=model.frame(fm1,data)

  dat1=model.matrix(fm1,data)
  y=df1[,1]
  y=data.frame(y)
  if(ncol(y)==1)colnames(y)=yv
  name1=colnames(dat1)

  res=list(df=df1,data=dat1,y=y,name=name1)
  return(res)
}

## Needed DEc 18, 2019
random_data=function(name0,name2,cls,df,...){
  M=matrix(0,length(name0),length(name2))
  for(i in 1:length(name0)){
    id=grep(name0[i],name2)
    M[i,id]=1
  }
  for(i in 1:length(name2)){
    id=which(M[,i]==1)
    str=paste(name0[id],collapse=":")
    if(i==1)vname=str
    else vname=c(vname,str)
  }
  uname2=unique(vname)
  M2=matrix(0,length(name0),length(uname2))
  for(i in 1:length(uname2)){
    id=which(vname==uname2[i])
    if(length(id)>0)M2[,i]=M[,id[1]]

  }
  colnames(M2)=uname2
  cls2=NULL
  for(i in 1:length(uname2)){
    #i=5
    id=which(M2[,i]==1)
    a=CFactor(df,id)
    cls2=c(cls2,"factor")
    if(i==1)dat=a
    else dat=data.frame(dat,a)
  }
  if(length(uname2)==1){dat=data.frame(dat);colnames(dat)=uname2}
  #head(dat)
  a=list(dat=dat,cls=cls2,uname=uname2)
  return(a)
}

## Needed Dec 18, 2019
fixed_data=function(name0,name1,cls,df,...){
  M=matrix(0,length(name0),length(name1))
  for(i in 1:length(name0)){
    id=grep(name0[i],name1)
    M[i,id]=1
  }
  id=which(name1=="(Intercept)")
  if(length(id)>0){
    name1=name1[-id]
    M=M[,-id]
    if(length(name1)==1){
      M=data.frame(M)
      colnames(M)=name1
    }
    #class(M)

  }
  if(length(name1)==0){
    X1=rep(1,nrow(df))
    X1=data.frame(X1)
    cls1="numeric"
    dat=X1
  }
  else if(length(name1)>0){
    for(i in 1:length(name1)){
      id=which(M[,i]==1)
      str=paste(name0[id],collapse=":")
      if(i==1)vname=str
      else vname=c(vname,str)
    }
    uname1=unique(vname)
    M1=matrix(0,length(name0),length(uname1))
    for(i in 1:length(uname1)){
      id=which(vname==uname1[i])
      if(length(id)>0)M1[,i]=M[,id[1]]

    }
    colnames(M1)=uname1
    cls1=NULL
    for(i in 1:length(uname1)){
      #i=5
      id=which(M1[,i]==1)
      #if(length(id)==1)a=df[,id]
      #else if(length(id)>1){
      clst=cls[id]
      idn=which(clst=="numeric")
      idi=which(clst=="integer")
      idt=unique(c(idn,idi))
      if(length(idt)==length(id)){
        a=CNumeric(df,id)
        cls1=c(cls1,"numeric")
      }
      else {
        a=CFactor(df,id)
        cls1=c(cls1,"factor")
      }
      if(i==1)dat=a
      else dat=data.frame(dat,a)
    }
    if(length(uname1)==1){dat=data.frame(dat);colnames(dat)=uname1}
    #head(dat)
  }
  a=list(dat=dat,cls=cls1,uname=uname1)
  return(a)
}






mixed <- function(x, ...) UseMethod("mixed")
#a=mixed(y~1|G*B)
#formula=y~1|G*B







reg.formula=function(x,y){
  fmla <- as.formula(paste("y ~ ", paste(x, collapse= "+")))
  return(fmla)
}


print.Umatrix=function(U,...){
  mk=length(U)
  cat("These are very complicated matrices:\n")
  cat("We only provide summarized information:\n")
  for(i in 1:mk){
    cat("Dimensions for matrix", i,":", "row=", nrow(U[[i]]), "col=", ncol(U[[i]]), "\n")
  }
}
print.Matrix=function(X,...){
  row=nrow(X)
  if(row>=100){
    cat("This is very large data set:\n")
    cat("We only provide summarized information.\n")
    cat("The sample size is ", nrow(X), " and the number of variables is ", ncol(X), "\n")
  }
}


## Needed. Dec 18, 2019
CFactor=function(M,v){

  n=nrow(M)
  #m1=M[,v]
  #Factor=NULL
  #NAME=NULL
  r=length(v)
  ok=complete.cases(M[,v])
  for(i in 1:r){
    if(i==1){
      NAME=colnames(M)[v[i]]
      Factor=M[,v[i]]
    }
    else if(i>1){
      NAME=paste(NAME,colnames(M)[v[i]],sep=":")
      Factor=paste(Factor,M[,v[i]],sep=":")
    }
  }
  #F=matrix(0,nrow=n,ncol=1)
  #F[,1]=Factor
  #colnames(F)=NAME
  n=length(ok)
  id=which(ok==F)
  Factor[id]=NA
  return(Factor)
}

DropColumns = function(A){

  for(i in colnames(A)){
    x=sum(ifelse(A[,which(colnames(A)==i)]==0, 0, 1))

    if(x==0){
      A= A[, -which(colnames(A)==i)]
    }

  }
  return(A)

}

## Needed Dec 19, 2019
CNumeric=function(M,v){

  n=nrow(M)
  r=length(v)
  y=numeric(n)
  for(i in 1:r){
    if(i==1){
      NAME=colnames(M)[v[i]]
      a=M[,v[i]]
    }
    else if(i>1){
      NAME=paste(NAME,colnames(M)[v[i]],sep=":")
      a=a*M[,v[i]]
    }
  }
  y=a
  return(y)
}
#X=gdata$X
## Needed Dec 19, 2019
BigX=function(X){  ## X is a list of matrices

  X0=NULL
  xk=length(X)
  xnames=NULL
  #i=2
  for(i in 1:xk)X0=cbind(X0,as.matrix(X[[i]]))
  for(i in 1:xk)xnames=c(xnames,colnames(X[[i]]))
  colnames(X0)=xnames
  return(X0)
}

## Needed 19, 2019

LargeInv=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  else if(m<=2000)return(subinv1000(A))
  else if(m<=4000)return(subinv2000(A))
  else if(m<=8000)return(subinv4000(A))
  else if(m<=16000)return(subinv8000(A))
  else if(m<=32000)return(subinv16000(A))
  else if(m<=64000)return(subinv32000(A))
  else if(m<=128000)return(subinv64000(A))

  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }
  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=100){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }

  else if(p>500&&p<=1000){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000&&p<=2000){
    S10=subinv1000(S)
    P1=subinv1000(P-Q%*%S10%*%R)
  }
  else if(p>2000&&p<=4000){
    S10=subinv2000(S)
    P1=subinv2000(P-Q%*%S10%*%R)
  }
  else if(p>4000&&p<=8000){
    S10=subinv4000(S)
    P1=subinv4000(P-Q%*%S10%*%R)
  }
  else if (p>8000&&p<=16000){
    S10=subinv8000(S)
    P1=subinv8000(P-Q%*%S10%*%R)
  }
  else if (p>16000&&p<=32000){
    S10=subinv16000(S)
    P1=subinv16000(P-Q%*%S10%*%R)
  }
  else if (p>32000&&p<=64000){
    S10=subinv32000(S)
    P1=subinv32000(P-Q%*%S10%*%R)
  }
  else if (p>64000&&p<=128000){
    S10=subinv64000(S)
    P1=subinv64000(P-Q%*%S10%*%R)
  }
  else{
    S10=subinv128000(S)
    P1=subinv128000(P-Q%*%S10%*%R)
  }
  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)

}

subinv128000=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  else if(m<=2000)return(subinv1000(A))
  else if(m<=4000)return(subinv2000(A))
  else if(m<=8000)return(subinv4000(A))
  else if(m<=16000)return(subinv8000(A))
  else if(m<=32000)return(subinv16000(A))
  else if(m<=64000)return(subinv32000(A))
  else if(m<=128000)return(subinv64000(A))

  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }
  else if(p>500&&p<=1000){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000&&p<=2000){
    S10=subinv1000(S)
    P1=subinv1000(P-Q%*%S10%*%R)
  }
  else if(p>2000&&p<=4000){
    S10=subinv2000(S)
    P1=subinv2000(P-Q%*%S10%*%R)
  }
  else if(p>4000&&p<=8000){
    S10=subinv4000(S)
    P1=subinv4000(P-Q%*%S10%*%R)
  }
  else if(p>8000&&p<=16000){
    S10=subinv8000(S)
    P1=subinv8000(P-Q%*%S10%*%R)
  }
  else if(p>16000&&p<=32000){
    S10=subinv16000(S)
    P1=subinv16000(P-Q%*%S10%*%R)
  }
  else if(p>32000&&p<=64000){
    S10=subinv32000(S)
    P1=subinv32000(P-Q%*%S10%*%R)
  }
  else if(p>64000){
    S10=subinv64000(S)
    P1=subinv64000(P-Q%*%S10%*%R)
  }

  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}



subinv64000=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  else if(m<=2000)return(subinv1000(A))
  else if(m<=4000)return(subinv2000(A))
  else if(m<=8000)return(subinv4000(A))
  else if(m<=16000)return(subinv8000(A))
  else if(m<=32000)return(subinv16000(A))
  else if(m<=64000)return(subinv32000(A))



  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }
  else if(p>500&&p<=1000){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000&&p<=2000){
    S10=subinv1000(S)
    P1=subinv1000(P-Q%*%S10%*%R)
  }
  else if(p>2000&&p<=4000){
    S10=subinv2000(S)
    P1=subinv2000(P-Q%*%S10%*%R)
  }
  else if(p>4000&&p<=8000){
    S10=subinv4000(S)
    P1=subinv4000(P-Q%*%S10%*%R)
  }
  else if(p>8000&&p<=16000){
    S10=subinv8000(S)
    P1=subinv8000(P-Q%*%S10%*%R)
  }
  else if(p>16000&&p<=32000){
    S10=subinv16000(S)
    P1=subinv16000(P-Q%*%S10%*%R)
  }
  else if(p>32000){
    S10=subinv32000(S)
    P1=subinv32000(P-Q%*%S10%*%R)
  }


  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}




subinv32000=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  else if(m<=2000)return(subinv1000(A))
  else if(m<=4000)return(subinv2000(A))
  else if(m<=8000)return(subinv4000(A))
  else if(m<=16000)return(subinv8000(A))
  else if(m<=32000)return(subinv16000(A))



  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }
  else if(p>500&&p<=1000){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000&&p<=2000){
    S10=subinv1000(S)
    P1=subinv1000(P-Q%*%S10%*%R)
  }
  else if(p>2000&&p<=4000){
    S10=subinv2000(S)
    P1=subinv2000(P-Q%*%S10%*%R)
  }
  else if(p>4000&&p<=8000){
    S10=subinv4000(S)
    P1=subinv4000(P-Q%*%S10%*%R)
  }
  else if(p>8000&&p<=16000){
    S10=subinv8000(S)
    P1=subinv8000(P-Q%*%S10%*%R)
  }
  else if(p>16000){
    S10=subinv16000(S)
    P1=subinv16000(P-Q%*%S10%*%R)
  }

  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}




subinv16000=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  else if(m<=2000)return(subinv1000(A))
  else if(m<=4000)return(subinv2000(A))
  else if(m<=8000)return(subinv4000(A))
  else if(m<=16000)return(subinv8000(A))

  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }
  else if(p>500&&p<=1000){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000&&p<=2000){
    S10=subinv1000(S)
    P1=subinv1000(P-Q%*%S10%*%R)
  }
  else if(p>2000&&p<=4000){
    S10=subinv2000(S)
    P1=subinv2000(P-Q%*%S10%*%R)
  }
  else if(p>4000&&p<=8000){
    S10=subinv4000(S)
    P1=subinv4000(P-Q%*%S10%*%R)
  }
  else if(p>8000){
    S10=subinv8000(S)
    P1=subinv8000(P-Q%*%S10%*%R)
  }

  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}

subinv8000=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  else if(m<=2000)return(subinv1000(A))
  else if(m<=4000)return(subinv2000(A))
  else if(m<=8000)return(subinv8000(A))

  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }
  else if(p>500&&p<=1000){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000&&p<=2000){
    S10=subinv1000(S)
    P1=subinv1000(P-Q%*%S10%*%R)
  }
  else if(p>2000&&p<=4000){
    S10=subinv2000(S)
    P1=subinv2000(P-Q%*%S10%*%R)
  }
  else if(p>4000){
    S10=subinv4000(S)
    P1=subinv4000(P-Q%*%S10%*%R)
  }

  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}

subinv4000=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  else if(m<=2000)return(subinv1000(A))
  else if(m<=4000)return(subinv2000(A))


  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }
  else if(p>500&&p<=1000){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000&&p<=2000){
    S10=subinv1000(S)
    P1=subinv1000(P-Q%*%S10%*%R)
  }
  else if(p>2000){
    S10=subinv2000(S)
    P1=subinv2000(P-Q%*%S10%*%R)
  }


  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}




subinv2000=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  else if(m<=2000)return(subinv1000(A))

  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }
  else if(p>500&&p<=1000){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000){
    S10=subinv1000(S)
    P1=subinv1000(P-Q%*%S10%*%R)
  }

  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}

subinv1000=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  else if(m<=1000)return(subinv500(A))
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200&&p<=500){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }
  else if(p>500){
    S10=subinv500(S)
    P1=subinv500(P-Q%*%S10%*%R)
  }

  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}



subinv500=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  else if(m<=200)return(subinv100(A))
  else if(m<=500)return(subinv200(A))
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }
  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else if(p>50&&p<=100){
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  else if(p>100&&p<=200){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }
  else if(p>200){
    S10=subinv200(S)
    P1=subinv200(P-Q%*%S10%*%R)
  }


  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}

subinv200=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  else if(m<=200)return(subinv100(A))
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=100){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }

  else if(p>100){
    S10=subinv100(S)
    P1=subinv100(P-Q%*%S10%*%R)
  }

  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)

}

#A=vi
subinv100=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  #else if(m<=50)return(subinv20(A))
  #else if(m<=100)return(subinv50(A))
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=20){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else if(p>20&&p<=50){
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }

  else{
    S10=subinv50(S)
    P1=subinv50(P-Q%*%S10%*%R)
  }
  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)

}


subinv50=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  else if(m<=50)return(subinv20(A))
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  if(p<=4){
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
  }
  else{
    S10=subinv20(S)
    P1=subinv20(P-Q%*%S10%*%R)
  }


  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)

}

subinv20=function(A){
  sys=1
  m=ncol(A)
  if(m<=100)return(ginv(A))
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      if(A[i,j]!=A[j,i]){
        sys=0
        break
      }
    }
    if(sys==0)break
  }

  p=as.integer(m/2)
  s=m-p
  P=A[1:p,1:p]
  R=A[(p+1):m,1:p]
  Q=A[1:p,(p+1):m]
  S=A[(p+1):m,(p+1):m]

  S10=ginv(S)
  P1=ginv(P-Q%*%S10%*%R)
  Q1=-P1%*%Q%*%S10
  if(sys==0)R1=-(S10%*%R)%*%P1
  S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
  A1=cbind(P1,Q1)
  if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
  else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
  return(A1)
}





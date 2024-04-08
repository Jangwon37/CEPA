# Lin-Ying estimator for uni-censored data

library(survival);

linying = function(X,t1,t2) {
	tm = c(t1,t2);
	R = survLY2.fit(X[,1][,1],X[,2][,1],X[,1][,2],X[,2][,2],tm);
	return( R$surv.prob );
}

survLY2.fit=function(y1, y2, d1, d2, time.matrix=NULL){

ctime=pmax(y1, y2)
dc=1-d1*d2

fitc=survfit(Surv(ctime, dc)~1)


if(length(time.matrix)==0)
  {
   grd1=sort(unique(y1[d1==1]))
   grd2=sort(unique(y2[d2==1]))
   n1=length(grd1)
   n2=length(grd2)


   survp=matrix(0, n1, n2)
   for(i in 1:n1)
   for(j in 1:n2)
      {time0=max(grd1[i], grd2[j])
       weight=max(fitc$surv[fitc$time>time0])  
       survp[i,j]=mean((y1>=grd1[i])*(y2>=grd2[j]))/weight
       }

   result=list(t1.grid=grd1, t2.grid=grd2, surv.prob=survp)
   }

if(length(time.matrix)>0)
  {
   time.matrix=matrix(time.matrix, ncol=2)
   K=length(time.matrix[,1])
   
   survp=rep(0, K)
   for(i in 1:K)
      {t10=time.matrix[i,1]
       t20=time.matrix[i,2]

       time0=max(t10, t20)
       weight=max(fitc$surv[fitc$time>time0])  
       survp[i]=mean((y1>=t10)*(y2>=t20))/weight
      }


   result=data.frame(t1=time.matrix[,1], t2=time.matrix[,2], surv.prob=survp)
   }

   return(result)
}



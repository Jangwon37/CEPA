################################################################
#######(y1, y2):    observed bivariate survival times  #########
#######(d1, d2):    censoring indicators               #########
########time.matrix:two-column matrix, whose rows are 
########            pairs of time at which the survival
########            probability will be calculated.
########            If left as NULL, then the survival 
########            probability will be calculated for
########            all the possible grid points expanded
########            by observed bivariate survival times 
########Outputs are survival probability at given time #########                     
########If the survival probability is zero, then it   #########
########is not identifiable at that time point         #########
################################################################

dabrowska = function(y1, y2, d1, d2, time.matrix=NULL) {

	n=length(y1)

	index.order=order(y2)
	y1=y1[index.order]
	y2=y2[index.order]
	d1=d1[index.order]
	d2=d2[index.order]

	n1=length(unique(y1))
	n2=length(unique(y2))

	grd1=sort(unique(y1))
	grd2=sort(unique(y2))


	R=K11=K10=K01=matrix(0, n1, n2)

	for( j in 1:n2 ) {
		id=(1:n)[y2==grd2[j] & d2==1]
    		count=0
    		for( i in id ) count=count+(grd1<=y1[i]) 
    		K01[,j]=count
    	}


	for(i in 1:n1) {
		id=(1:n)[y1==grd1[i] & d1==1]
    		count=0
    		for( j in id ) count=count+(grd2<=y2[j]) 
    		K10[i,]=count
    	} 

	for(i in 1:n1) {
		id=(1:n)[y1==grd1[i] & d1==1]
    		count=0
    		for(j in id) count=count+(grd2==y2[j])*d2[j] 
    		K11[i,]=count
    	} 

	count.tot=rep(1:n2, table(y2))

	for( i in 1:n1 ) {
		id=(1:n)[y1>=grd1[i]]
    		number=as.vector(table(count.tot[id]))
    		index=unique(count.tot[id])
    		tot=length(id)
    		k=length(index)
    		final=rep((tot-c(0, cumsum(number[-k]))), c(index[1], index[-1]-index[-k]))
    		R[i,]=c(final, rep(0, n2-length(final)))
    	}   

	L=R*(R-K10-K01+K11)/(R-K10)/(R-K01)

	if( length(time.matrix)==0 ) {
		fit1=survfit(Surv(y1, d1)~1)
   		fit2=survfit(Surv(y2, d2)~1)

   		surv1=c(1, fit1$surv)
   		surv2=c(1, fit2$surv)
   		time1=c(0, fit1$time)
   		time2=c(0, fit2$time)

   		s1=rep(1, n1)
   		s2=rep(1, n2)
   		for( i in 1:n1 ) { s1[i]=min(surv1[time1<=grd1[i]]) }
   		for( i in 1:n2 ) { s2[i]=min(surv2[time2<=grd2[i]]) }
   
   		r.prod=apply(log(L), 1, cumsum)
   		r.prod=exp(apply(r.prod, 1, cumsum))
   		surv.prob=t(t(r.prod*s1)*s2)

   		surv.prob[is.na(surv.prob)]=0  

   		return(list(t1.grid=grd1, t2.grid=grd2, surv.prob=surv.prob))
   	}

	if( length(time.matrix)>0 ) {
		fit1=survfit(Surv(y1, d1)~1)
   		fit2=survfit(Surv(y2, d2)~1)

   		surv1=c(1, fit1$surv)
   		surv2=c(1, fit2$surv)
  		time1=c(0, fit1$time)
   		time2=c(0, fit2$time)

   		time.matrix=matrix(time.matrix, ncol=2)
   		K=length(time.matrix[,1])
   
   		surv.prob=rep(0, K)
   		for(i in 1:K) {
			t10=time.matrix[i,1]
       			t20=time.matrix[i,2]

			# added by junhee
			if( t10 < 0 ) t10 = 0;
			if( t20 < 0 ) t20 = 0;

       			index1=(1:n1)[grd1<=t10]
       			index2=(1:n2)[grd2<=t20]

       			rprod=exp(sum(log(L[index1, index2])))

       			s1=min(surv1[time1<=t10])
       			s2=min(surv2[time2<=t20])

       			surv.prob[i]=s1*s2*rprod

			# added by junhee
			if( (is.infinite(t10)&t10>0) | (is.infinite(t20)&t20>0) ) surv.prob[i] = 0;
      		}

		# added by junhee
		surv.prob[is.na(surv.prob)] = 0;

   		result=data.frame(t1=time.matrix[,1], t2=time.matrix[,2], surv.prob=surv.prob)
   		return(result)
	}
}
#tmps=rep(1,28)
#for(z in 1:27){
#  tmps[z+1]=tmps[z]*((n-sum(y1<=z+1))+sum((y1==z+1)*(d1==0)))/(n-sum((y1<=z)))
#}
#tmps=tmps[-1]

##############################################################################
# 2d OPT for censored data
#
# by Junhee Seok
##############################################################################

##############################################################################
# Required Libraries
##############################################################################

library(survival);

##############################################################################
# main OPT functions
##############################################################################

# opt2d
opt2d = function(X,A,max.iter=5,min.err=1E-3,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {

	fn.seed = sprintf('%d',round(runif(1)*10000));

	# options
	options = "";
	options = sprintf("%s -s %f-%f,%f-%f",options,A[1],A[2],A[3],A[4]);
	options = sprintf("%s -i %d -e %f",options,max.iter,min.err);
	options = sprintf("%s -d %d -n %f -r %f -a %f",options,max.depth,max.n,rho,alpha);
	if( int.part ) { options = sprintf("%s -I",options); }
	if( mixed.mode ) { options = sprintf("%s -M",options); }
	out.fn = paste("out",fn.seed,sep='.');
	options = sprintf("%s -o %s",options,out.fn);

	# data file
	d.out = cbind(X[,1][,1],X[,1][,2],X[,2][,1],X[,2][,2]);
	data.fn = paste("data",fn.seed,sep='.');
	write.table(d.out,file=data.fn,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');

	# command
	command = sprintf("cOPT/./copt %s %s",options,data.fn);
	system(command);

	# read output
	S = as.matrix(read.table(out.fn));
	system(sprintf('rm -rf *.%s',fn.seed));

	return( S );	
}

# detail
opt2d.detail = function(X,A,max.iter=5,min.err=1E-3,n.grid=11,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {
	S.iter = vector('list',1);
	Dp = Sp = NULL; 
	err = 1; n.iter = 1;
	while( n.iter <= max.iter & err > min.err ) {
		S = opt2d.one(X,A,Sp=Sp,max.depth=max.depth,max.n=max.n,rho=rho,alpha=alpha,int.part=int.part,mixed.mode=mixed.mode,verb=verb);
		D = opt2d.density_grid(S,seq(A[1],A[2],length=n.grid),seq(A[3],A[4],length=n.grid));
		if( !is.null(Dp) ) { err = sqrt(mean((D-Dp)^2)); }
		Sp = S; Dp = D;
		S.iter[[n.iter]] = S;
		cat('iter',n.iter,sprintf('%.5f',err),'\n');
		n.iter = n.iter+1;
	}

	return( list(S=S,S.iter=S.iter) );	
}

# opt2d once
opt2d.one = function(X,A,Sp=NULL,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {

	fn.seed = sprintf('%d',round(runif(1)*10000));

	# options
	options = "-O";
	options = sprintf("%s -s %f-%f,%f-%f",options,A[1],A[2],A[3],A[4]);
	options = sprintf("%s -d %d -n %f -r %f -a %f",options,max.depth,max.n,rho,alpha);
	if( !is.null(Sp) ) {
		fn = paste("prior",fn.seed,sep='.');
		write.table(Sp[,1:6],file=fn,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');
		options = sprintf("%s -p %s",options,fn);
	}
	if( int.part ) { options = sprintf("%s -I",options); }
	if( mixed.mode ) { options = sprintf("%s -M",options); }
	if( verb ) { options = sprintf("%s -v",options); }
	out.fn = paste("out",fn.seed,sep='.');
	options = sprintf("%s -o %s",options,out.fn);

	# data file
	d.out = cbind(X[,1][,1],X[,1][,2],X[,2][,1],X[,2][,2]);
	data.fn = paste("data",fn.seed,sep='.');
	write.table(d.out,file=data.fn,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');

	# command
	command = sprintf("./copt %s %s",options,data.fn);
	system(command);

	# read output
	S = as.matrix(read.table(out.fn));

	# remove tmp files
	system(sprintf('rm -rf *.%s',fn.seed));

	return( S );
}

##############################################################################
# partition functions
##############################################################################

# get partitions from a PHI table, inside handler
opt2d.partition = function(PHI,X,A,Aroi=NULL,int.part=FALSE,verb=FALSE) {
	if( is.null(Aroi) ) { Aroi = A; A[2] = A[4] = Inf; }
	N.ALL = nrow(X);
	S = opt2d.partition.sub(PHI,A,Aroi,N.ALL,int.part=int.part,verb=verb,depth=1);
	S[,6] = S[,6]/sum(S[,6]);
	S[,5] = S[,6]/(S[,2]-S[,1])/(S[,4]-S[,3]);
	o = order(S[,1],S[,2],S[,3],S[,4]);
	S = S[o,];
	return( S );
}

opt2d.partition.sub = function(PHI,A,Aroi,N.ALL,PARTITION=NULL,int.part=FALSE,verb=FALSE,depth=1) {

	if( verb ) { cat(rep("  ",depth-1),"enter: ",A,'\n'); }
    
	A1 = A[1:2]; A2 = A[3:4];
	Aroi1 = Aroi[1:2]; Aroi2 = Aroi[3:4];

        idx1 = PHI[,1]==A1[1] & PHI[,2]==A1[2];
	idx2 = PHI[,3]==A2[1] & PHI[,4]==A2[2];
	idx = idx1 & idx2;
	if( sum(idx) == 0 ) { 
		if( verb ) { cat(rep("  ",depth-1),"exit:  ",A,0,'\n'); }
                PARTITION = rbind(PARTITION,c(A,0,0,-1));
		return( PARTITION ); 
	}
        phi = PHI[idx,];

	d1 = d2 = -Inf;
	stop.partitioning = FALSE;
	if( phi[6] > phi[7] & !(A1[1]!=Aroi1[2]&is.infinite(A1[2])) & !(A2[1]!=Aroi2[2]&is.infinite(A2[2])) ) { stop.partitioning = TRUE; }
	if( verb ) { cat(rep("  ",depth-1),"mid:   ",A1,A2,d1,d2,phi[6:9],stop.partitioning,'\n'); }
	if( A1[1]==Aroi1[2] & is.infinite(A1[2]) ) { stop.partitioning = TRUE; }
	if( verb ) { cat(rep("  ",depth-1),"mid:   ",A1,A2,d1,d2,phi[6:9],stop.partitioning,'\n'); }
	if( A2[1]==Aroi2[2] & is.infinite(A2[2]) ) { stop.partitioning = TRUE; }
	if( verb ) { cat(rep("  ",depth-1),"mid:   ",A1,A2,d1,d2,phi[6:9],stop.partitioning,'\n'); }
	iidx1 = idx2 & (PHI[,1]==A1[1] | PHI[,2]==A1[2]) & PHI[,12]==phi[12]+1;
	if( sum(iidx1) > 0 ) { 
		d1 = setdiff(unique(as.vector(PHI[iidx1,1:2,drop=FALSE])),A1); 
		d1 = d1[ A1[1] < d1 & d1 < A1[2] ];
		if( length(d1) > 0 ) { d1 = min(d1); }
		else { d1 = -Inf; }
	}
	iidx2 = idx1 & (PHI[,3]==A2[1] | PHI[,4]==A2[2]) & PHI[,12]==phi[12]+1;
	if( sum(iidx2) > 0 ) { 
		d2 = setdiff(unique(as.vector(PHI[iidx2,3:4,drop=FALSE])),A2); 
		d2 = d2[ A2[1] < d2 & d2 < A2[2] ];
		if( length(d2) > 0 ) { d2 = min(d2); }
		else { d2 = -Inf; }
	}
	if( sum(iidx1) == 0 & sum(iidx2) == 0 ) { stop.partitioning = TRUE; }	

	if( verb ) { cat(rep("  ",depth-1),"mid:   ",A1,A2,d1,d2,phi[6:9],stop.partitioning,'\n'); }

	if( stop.partitioning ) {
		if( d1 == Aroi1[2] & d2 == Aroi2[2] ) {
                        PT1 = opt2d.partition.sub(PHI=PHI,A=c(A1[1],d1,A2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
                        PT1 = opt2d.partition.sub(PHI=PHI,A=c(d1,A1[2],A2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PT1,int.part=int.part,verb=verb,depth=depth+1);
                        PT2 = opt2d.partition.sub(PHI=PHI,A=c(A1,A2[1],d2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
                        PT2 = opt2d.partition.sub(PHI=PHI,A=c(A1,d2,A2[2]),Aro=Aroi,N.ALL=N.ALL,PARTITION=PT2,int.part=int.part,verb=verb,depth=depth+1);
			if( nrow(PT1) >= nrow(PT2) ) { PARTITION = PT1; } 
			else { PARTITION = PT2; }
		} else if( d1 == Aroi1[2] ) {
                        PARTITION = opt2d.partition.sub(PHI=PHI,A=c(A1[1],d1,A2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
                        PARTITION = opt2d.partition.sub(PHI=PHI,A=c(d1,A1[2],A2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
		} else if( d2 == Aroi2[2] ) {
                        PARTITION = opt2d.partition.sub(PHI=PHI,A=c(A1,A2[1],d2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
                        PARTITION = opt2d.partition.sub(PHI=PHI,A=c(A1,d2,A2[2]),Aro=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
		} else {
                	PARTITION = rbind(PARTITION,c(A,phi[11]/N.ALL/(A1[2]-A1[1])/(A2[2]-A2[1]),phi[11]/N.ALL,phi[12]));
		}
	} else {
		if( phi[8] == phi[9] & d1>0 & d2>0 ) {
                        PT1 = opt2d.partition.sub(PHI=PHI,A=c(A1[1],d1,A2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
                        PT1 = opt2d.partition.sub(PHI=PHI,A=c(d1,A1[2],A2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PT1,int.part=int.part,verb=verb,depth=depth+1);
                        PT2 = opt2d.partition.sub(PHI=PHI,A=c(A1,A2[1],d2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
                        PT2 = opt2d.partition.sub(PHI=PHI,A=c(A1,d2,A2[2]),Aro=Aroi,N.ALL=N.ALL,PARTITION=PT2,int.part=int.part,verb=verb,depth=depth+1);
			if( nrow(PT1) >= nrow(PT2) ) { PARTITION = PT1; } 
			else { PARTITION = PT2; }
                } else if( phi[8] >= phi[9] & d1>0 ) {
                        PARTITION = opt2d.partition.sub(PHI=PHI,A=c(A1[1],d1,A2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
                        PARTITION = opt2d.partition.sub(PHI=PHI,A=c(d1,A1[2],A2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
		} else if( phi[8] <= phi[9] & d2>0 ) {
                        PARTITION = opt2d.partition.sub(PHI=PHI,A=c(A1,A2[1],d2),Aroi=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
                        PARTITION = opt2d.partition.sub(PHI=PHI,A=c(A1,d2,A2[2]),Aro=Aroi,N.ALL=N.ALL,PARTITION=PARTITION,int.part=int.part,verb=verb,depth=depth+1);
		} else {
                	PARTITION = rbind(PARTITION,c(A,phi[11]/N.ALL/(A1[2]-A1[1])/(A2[2]-A2[1]),phi[11]/N.ALL,phi[12]));
		}
	}

	if( verb ) { cat(rep("  ",depth-1),"exit:  ",A,d1,d2,'\n'); }
        return( PARTITION );
}

##############################################################################
# probability handlers
##############################################################################

# probability in A, P(X in A) when A is finite
opt2d.p.area.sub = function(S,A) {
	P = 0; n = nrow(S);
	a1 = S[,1]; a1[a1<A[1]] = A[1];
	b1 = S[,2]; b1[b1>A[2]] = A[2];
	a2 = S[,3]; a2[a2<A[3]] = A[3];
	b2 = S[,4]; b2[b2>A[4]] = A[4];
	xx = b1-a1; xx[xx<0] = 0;
	yy = b2-a2; yy[yy<0] = 0;
	P = sum(xx*yy*S[,5]);

	return( P );
}

# probability in A, P(X in A)
opt2d.p.area = function(S,A) {
	if( all(is.finite(A)) ) { return( opt2d.p.area.sub(S,A) ); }
	if( is.infinite(A[1]) | is.infinite(A[3]) ) { return( 0 ); }

	# if A is an open space, we recalculate the density by replacing Inf to a large value
	t = S[,1:2]; r1 = range(t[is.finite(t)]);
	t = S[,3:4]; r2 = range(t[is.finite(t)]);
	L = 10*max(c(r1,r2));
	if( is.infinite(A[2]) ) { S[is.infinite(S[,2]),2] = L; A[2] = L; }
	if( is.infinite(A[4]) ) { S[is.infinite(S[,4]),4] = L; A[4] = L; }
	S[,5] = S[,6]/(S[,2]-S[,1])/(S[,4]-S[,3]);
	
	return( opt2d.p.area.sub(S,A) );
}

# survival probablity P(T1>=t1, T2>=t2)
opt2d.p = function(S,t1,t2) {
	return( opt2d.p.area(S,c(t1,Inf,t2,Inf)) );
}

# obtain a marginal distribution
opt2d.marginal = function(S,IND) {

	SM = NULL;
	if( IND == 1 ) {
		t = unique(sort(as.vector(S[,1:2])));
		for( i in 2:length(t) ) {
			p = opt2d.p.area(S,c(t[i-1],t[i],0,Inf));
			SM = rbind(SM,c(t[i-1],t[i],p/(t[i]-t[i-1]),p));
		}
	} else {
		t = unique(sort(as.vector(S[,3:4])));
		for( i in 2:length(t) ) {
			p = opt2d.p.area(S,c(0,Inf,t[i-1],t[i]));
			SM = rbind(SM,c(t[i-1],t[i],p/(t[i]-t[i-1]),p));
		}
	}

	return( SM );
}

##############################################################################
# plotting functions
##############################################################################

# plot partition
opt2d.plot.partition = function(S,X,xlab="",ylab="",main="",xlim=NULL,ylim=NULL) {

	if( is.null(xlim) ) { t = S[,1:2]; xlim = range(t[is.finite(t)]); } 
	if( is.null(ylim) ) { t = S[,3:4]; ylim = range(t[is.finite(t)]); } 
	large.a = 100*max(c(xlim,ylim));
	plot(0,0,xlim=xlim,ylim=ylim,col='white',xlab=xlab,ylab=ylab,main=main);
	for( i in 1:nrow(S) ) {
		A = S[i,];
		A[is.infinite(A)] = large.a;
		lines(c(A[2],A[2]),c(A[3],A[4]),col=2);
		lines(c(A[1],A[1]),c(A[3],A[4]),col=2);
		lines(c(A[1],A[2]),c(A[3],A[3]),col=2);
		lines(c(A[1],A[2]),c(A[4],A[4]),col=2);
	}
	points(X[,1][,1],X[,2][,1]);
	idx = X[,1][,2]==0;
	points(X[idx,1][,1],X[idx,2][,1],pch='|');
	idx = X[,2][,2]==0;
	points(X[idx,1][,1],X[idx,2][,1],pch=95,cex=1.3);
}

# plot one dimensional partition
opt2d.plot.marginal = function(S,X,xlab='Time',ylab='Probability',main='') {

	x = S[,1:2]; 
	r = range(x[is.finite(x)]);
	plot(survfit(X~1),conf.int=FALSE,xlim=c(r[1],r[2]*1.1),ylim=c(0,1.3),xlab=xlab,ylab=ylab,yaxt='n',main=main);
#	axis(2,at=seq(0,1,by=0.2));
	axis(1,at=r[2]*1.1,label=paste(round(max(r)*10)/10,'+',sep=""));

        x = 0; y = 1;
        S = rbind(c(0,S[1,1],0,0),S);
        for( i in 1:nrow(S) ) {
                x = c(x,S[i,2]);
                y = c(y,y[i]-S[i,3]*(S[i,2]-x[i]));
        }
	lines(x,y,type='l',col=4,lwd=2);

	x = S[,1:2]; 
	r = range(x[is.finite(x)]);
	x[is.infinite(x)] = r[2]*1.1;
	y = cbind(S[,4],S[,4]);
	lines(as.vector(t(x)),as.vector(t(y)),type='l',col=3,lty=1,lwd=2);

	x = unique(as.vector(S[,1:2])); x = x[is.finite(x)];
	for( i in 1:length(x) ) {
		lines(c(x[i],x[i]),c(-100,100),col=2,lty=2);
	}

	xx = X[,1]; yy = rep(1.2,nrow(X))+runif(nrow(X),-0.1,0.1);
	points(xx,yy);
	points(xx[X[,2]==0],yy[X[,2]==0],pch='|');
}

##############################################################################
# Misc.
##############################################################################

opt2d.density_grid = function(S,xgrid,ygrid) {

	xgrid = c(xgrid,Inf);
	ygrid = c(ygrid,Inf);

	D = matrix(rep(0,(length(xgrid)-1)*(length(ygrid)-1)),nrow=length(xgrid)-1);
	XN = rep("",length(xgrid)-1);
	YN = rep("",length(ygrid)-1);
	for( i in 2:length(xgrid) ) {
		for( j in 2:length(ygrid) ) {
			D[i-1,j-1] = opt2d.p.area(S,c(xgrid[i-1],xgrid[i],ygrid[j-1],ygrid[j]));
		}
	}
	for( i in 2:length(xgrid) ) { XN[i-1] = sprintf("[%.2f,%.2f)",xgrid[i-1],xgrid[i]); }
	for( i in 2:length(ygrid) ) { YN[i-1] = sprintf("[%.2f,%.2f)",ygrid[i-1],ygrid[i]); }
	row.names(D) = XN;
	colnames(D) = YN;
	
	return( D );
}


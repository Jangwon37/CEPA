##############################################################################
# 1d OPT for censored data
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

# opt1d
opt1d = function(X,A,max.iter=5,min.err=1E-3,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {

	fn.seed = sprintf('%d',round(runif(1)*10000));

	# options
	options = "";
	options = sprintf("%s -s %f-%f",options,A[1],A[2]);
	options = sprintf("%s -i %d -e %f",options,max.iter,min.err);
	options = sprintf("%s -d %d -n %f -r %f -a %f",options,max.depth,max.n,rho,alpha);
	if( int.part ) { options = sprintf("%s -I",options); }
	if( mixed.mode ) { options = sprintf("%s -M",options); }
	out.fn = paste("out",fn.seed,sep='.');
	options = sprintf("%s -o %s",options,out.fn);

	# data file
	d.out = cbind(X[,1],X[,2]);
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
opt1d.detail = function(X,A,max.iter=5,min.err=1E-3,n.grid=11,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {
	S.iter = vector('list',1);
	Dp = Sp = NULL; 
	err = 1; n.iter = 1;
	while( n.iter <= max.iter & err > min.err ) {
		S = opt1d.one(X,A,Sp=Sp,max.depth=max.depth,max.n=max.n,rho=rho,alpha=alpha,int.part=int.part,mixed.mode=mixed.mode,verb=verb);
		D = opt1d.density_grid(S,seq(A[1],A[2],length=n.grid));
		if( !is.null(Dp) ) { err = sqrt(mean((D-Dp)^2)); }
		Sp = S; Dp = D;
		S.iter[[n.iter]] = S;
		cat('iter',n.iter,sprintf('%.5f',err),'\n');
		n.iter = n.iter+1;
	}

	return( list(S=S,S.iter=S.iter) );	
}

# opt1d once
opt1d.one = function(X,A,Sp=NULL,max.depth=5,max.n=1,rho=0.5,alpha=0.5,int.part=FALSE,mixed.mode=FALSE,verb=FALSE) {

	fn.seed = sprintf('%d',round(runif(1)*10000));

	# options
	options = "-O";
	options = sprintf("%s -s %f-%f",options,A[1],A[2]);
	options = sprintf("%s -d %d -n %f -r %f -a %f",options,max.depth,max.n,rho,alpha);
	if( !is.null(Sp) ) {
		fn = paste("prior",fn.seed,sep='.');
		write.table(Sp[,1:4],file=fn,quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t');
		options = sprintf("%s -p %s",options,fn);
	}
	if( int.part ) { options = sprintf("%s -I",options); }
	if( mixed.mode ) { options = sprintf("%s -M",options); }
	if( verb ) { options = sprintf("%s -v",options); }
	out.fn = paste("out",fn.seed,sep='.');
	options = sprintf("%s -o %s",options,out.fn);

	# data file
	d.out = cbind(X[,1],X[,2]);
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
# probability handlers
##############################################################################

# probability in A, P(X in A) when A is finite
opt1d.p.area.sub = function(S,A) {
	P = 0; n = nrow(S);
	a1 = S[,1]; a1[a1<A[1]] = A[1];
	b1 = S[,2]; b1[b1>A[2]] = A[2];
	xx = b1-a1; xx[xx<0] = 0;
	P = sum(xx*S[,3]);

	return( P );
}

# probability in A, P(X in A)
opt1d.p.area = function(S,A) {
	if( all(is.finite(A)) ) { return( opt1d.p.area.sub(S,A) ); }
	if( is.infinite(A[1]) ) { return( 0 ); }

	# if A is an open space, we recalculate the density by replacing Inf to a large value
	t = S[,1:2]; r1 = range(t[is.finite(t)]);
	L = 10*max(r1);
	if( is.infinite(A[2]) ) { S[is.infinite(S[,2]),2] = L; A[2] = L; }
	S[,3] = S[,4]/(S[,2]-S[,1]);
	
	return( opt1d.p.area.sub(S,A) );
}

# survival probablity P(T>=t)
opt1d.p = function(S,t) {
	p = opt1d.p.area(S,c(t1,Inf));
	if( p > 1 ) p = 1;
	if( p < 0 ) p = 0;
	return(p);
}

##############################################################################
# plotting functions
##############################################################################

opt1d.plot = function(S,X,xlab='Time',ylab='Probability',main='') {

	x = S[,1:2]; 
	r = range(x[is.finite(x)]);
	plot(survfit(X~1),conf.int=FALSE,xlim=c(r[1],r[2]*1.1),ylim=c(0,1.3),xlab=xlab,ylab=ylab,yaxt='n',main=main);
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

opt1d.density_grid = function(S,xgrid) {

	xgrid = c(xgrid,Inf);

	D = rep(0,(length(xgrid)-1));
	XN = rep("",length(xgrid)-1);
	for( i in 2:length(xgrid) ) {
		D[i-1] = opt1d.p.area(S,c(xgrid[i-1],xgrid[i]));
	}
	for( i in 2:length(xgrid) ) { XN[i-1] = sprintf("[%.2f,%.2f)",xgrid[i-1],xgrid[i]); }
	names(D) = XN;
	
	return( D );
}


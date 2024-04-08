# library

####################################################################################################################
# Misc. Functiosn
####################################################################################################################

rowVars = function(E,na.rm=FALSE) {
	n = ncol(E);
	m = rowMeans(E,na.rm=na.rm);
	mm = rowMeans(E^2,na.rm=na.rm);
	v = (n*mm - n*m^2)/(n-1);
	return( v );
}

colVars = function(E,na.rm=FALSE) {
	E = t(E);
	n = ncol(E);
	m = rowMeans(E,na.rm=na.rm);
	mm = rowMeans(E^2,na.rm=na.rm);
	v = (n*mm - n*m^2)/(n-1);
	return( v );
}



############################################################################
# Statistics
############################################################################

t.test.mat = function(X,Y=NULL,equal.var=FALSE,report.t=FALSE,report.df=FALSE) {

	if( !is.null(Y) & equal.var==FALSE ) {	
		m1 = rowMeans(X,na.rm=TRUE); s1 = sqrt(rowVars(X,na.rm=TRUE)); n1 = ncol(X);
		m2 = rowMeans(Y,na.rm=TRUE); s2 = sqrt(rowVars(Y,na.rm=TRUE)); n2 = ncol(Y);
		s = sqrt( s1^2/n1 + s2^2/n2 );
		t = (m1-m2)/s;
		df = ( (s1^2/n1+s2^2/n2)^2 ) / ( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) );
		if( report.t ) { return( t ); }
		if( report.df ) { return( df ); }
		p = 2*(1-pt(abs(t),df));
	} else if( !is.null(Y) & equal.var==TRUE ) {
		m1 = rowMeans(X,na.rm=TRUE); s1 = sqrt(rowVars(X,na.rm=TRUE)); n1 = ncol(X);
		m2 = rowMeans(Y,na.rm=TRUE); s2 = sqrt(rowVars(Y,na.rm=TRUE)); n2 = ncol(Y);
		s = sqrt( ( (n1-1)*s1^2 + (n2-1)*s2^2 )/(n1+n2-2) );
		t = (m1-m2)/s/sqrt(1/n1+1/n2);	
		if( report.t ) { return( t ); }
		if( report.df ) { return( df ); }
		p = 2*(1-pt(abs(t),n1+n2-2));
	} else {
		m = rowMeans(X,na.rm=TRUE);
		s = sqrt(rowVars(X,na.rm=TRUE));
		t = m*sqrt(ncol(X))/s;
		if( report.t ) { return( t ); }
		if( report.df ) { return( df ); }
		p = 2*(1-pt(abs(t),ncol(X)-1));
	}
	
	return( p );
}

icdf = function(d,s,p,given.occurrence=FALSE) {
        if( given.occurrence ) p = 1-p*(1-s[length(s)]);
	ss = 1-cumsum(s); 
	kk = which(ss>p);
	if( length(kk) == 0 ) { 
		ii = -Inf;
	} else {
		ii = max(which(ss>p))+1;
	}
	if( is.na(ii) || is.infinite(ii) ) {
		d1 = min(d)-1; d2 = min(d);
		p1 = 1; p2 = ss[1];
	} else if( ii > length(ss) ) {
		return( NA );
	} else {
		d1 = d[ii-1]; d2 = d[ii]; p1 = ss[ii-1]; p2 = ss[ii];
	}
	m = (p-p1)*(d1-d2)/(p1-p2)+d1;
	return( m );
}

# median(x|y)
med.opt2d = function(S,given=1) {
        if( given == 2 ) S = S[,c(3,4,1,2,5,6)];
	d = p = NULL;
	tlist = 1:52;
	for( i in 1:length(tlist) ) {
		T1 = tlist[i];
		s = NULL;
		for( T2 in 2:52 ) s = c(s,opt2d.p.area(S,c(T1,T1+1,T2,T2+1)));
		s = c(s,opt2d.p.area(S,c(T1,T1+1,53,Inf)));
		p = c(p,sum(s)); s = s/sum(s);
		d = c(d,icdf(c(2:52,Inf)-T1+1,s,0.5));
	}
	return( cbind(tlist,d,p) );
}

############################################################################
# regression
############################################################################

reg.opt2d = function(S,response.var=2) {
        if( response.var == 1 ) S = S[,c(3,4,1,2,5,6)];

	x = y = 2:25;
	X =    matrix(rep(x,length(y)),ncol=length(y)); 
	Y = t( matrix(rep(y,length(x)),ncol=length(x)) ); 
	p.mat = matrix(rep(0,length(x)*length(y)),nrow=length(x));
	for( i in 1:length(x) ) {
		for( j in 1:length(y) ) p.mat[i,j] = opt2d.p.area(S,c(x[i],x[i]+1,y[j],y[j]+1));
	}
	p.mat = p.mat/sum(p.mat);
	px = rowSums(p.mat); py = colSums(p.mat);

	mx = sum(x*px); my = sum(y*py);
	mx2 = sum(x*x*px); mxy = sum(X*Y*p.mat);

	alpha = (mxy-mx*my)/(mx2-mx^2);
	beta = Y-alpha*X; mb = sum(p.mat*beta);

	d = Y-X; md = sum(p.mat*d); vd = sum(d*d*p.mat)-md^2;
	vx = mx2-mx^2; vb = sum(p.mat*beta*beta)-mb^2;
	#r = c((alpha-1)*mx/md,mb/md);
	#r = c( ((alpha-1)^2)*vx/vd,vb/vd);

	s = sqrt( (sum(beta*beta*p.mat)-sum(beta*p.mat)^2)/(mx2-mx^2) );
	z = (alpha-1)/s;
	t = (alpha)/s*sqrt(1875);

	return( c(alpha,sum(p.mat*beta),z,t) );

}

reg.opt3d = function(S,response.var=3) {
        if( response.var == 2 ) S = S[,c(1,2,5,6,3,4,7,8)];
        if( response.var == 1 ) S = S[,c(3,4,5,6,1,2,7,8)];

	x1 = x2 = y = 2:25;
	X1 =    matrix(rep(x1,length(x2)*length(y)),ncol=length(y)); 
	X2 =    matrix(rep(x2,each=length(x1)*length(y)),ncol=length(y)); 
	Y =  t( matrix(rep(y,each=length(x1)*length(x2)),ncol=length(x1)*length(x2)) ); 
	p.mat = matrix(rep(0,length(x1)*length(x2)*length(y)),ncol=length(y));
	for( i in 1:length(x1) ) {
		for( j in 1:length(x2) ) {
			for( k in 1:length(y) ) p.mat[(j-1)*length(x2)+i,k] = opt3d.p.area(S,c(x1[i],x1[i]+1,x2[j],x2[j]+1,y[k],y[k]+1));
		}
	}
	p.mat = p.mat/sum(p.mat);

	mx1 = sum(X1*p.mat); mx2 = sum(X2*p.mat); my = sum(Y*p.mat);
	XtX = matrix(c( sum(p.mat*(X1-mx1)^2), sum(p.mat*(X1-mx1)*(X2-mx2)), sum(p.mat*(X1-mx1)*(X2-mx2)), sum(p.mat*(X2-mx2)^2) ), nrow=2);
	XtY = c( sum(p.mat*(Y-my)*(X1-mx1)), sum(p.mat*(Y-my)*(X2-mx2)) );

	s1 = sqrt( (sum(beta*beta*p.mat)-sum(beta*p.mat)^2)/( sum(X1*X1*p.mat)-mx1^2 ) );
	s2 = sqrt( (sum(beta*beta*p.mat)-sum(beta*p.mat)^2)/( sum(X2*X2*p.mat)-mx1^2 ) );
	z1 = (alpha[1]-1)/s1;
	z2 = (alpha[2]-1)/s2;
	z = c(z1,z2);

	return( c(alpha,sum(p.mat*beta),z) );
}


############################################################################
# Mutual Information
############################################################################

# MI(x,y)
mi.opt2d = function(S) {
	tl = c(2:29,Inf);
	mi = 0;
	for( i in 1:(length(tl)-1) ) {
		px = opt2d.p.area(S,c(tl[i],tl[i+1],0,Inf));
		for( j in 1:(length(tl)-1) ) {
			py = opt2d.p.area(S,c(0,Inf,tl[j],tl[j+1]));
			pxy = opt2d.p.area(S,c(tl[i],tl[i+1],tl[j],tl[j+1]));
			if( px > 0 & py > 0 & pxy > 0 ) { mi = mi + pxy*log(pxy/px/py); }
		}
	}
	return(mi);
}

# MI(y-x,x)
mi.opt2d.diff = function(S,given=1,alpha=1,nbin=NULL,method='mi') {
        if( given == 2 ) S = S[,c(3,4,1,2,5,6)];

	x = y = 2:25;
	d.mat = matrix(rep(NA,length(x)*length(y)),ncol=length(y));
	for( i in 1:length(x) ) {
		for( j in 1:length(y) ) d.mat[i,j] = round( (y[j]-alpha*x[i])*10000 )/10000;
	}
	d = sort(unique(as.vector(d.mat)));
	if( !is.null(nbin) ) {
		dd = seq(min(d),max(d),length=nbin+1);
		d = dd[1:nbin];
	} else {
		d = seq(min(d),max(d),by=1);
	}

	p.mat = matrix(rep(0,length(d)*length(x)),nrow=length(d));
	pd = rep(0,length(d));
	for( i in 1:length(x) ) {
		for( j in 1:length(y) ) {
			pxy = opt2d.p.area(S,c(x[i],x[i]+1,y[j],y[j]+1));
			tmp.d = round( (y[j]-alpha*x[i])*10000 )/10000;
			iidx = which(d >= tmp.d);
			idx.d = length(d);
			if( length(iidx) > 0 ) idx.d = min(iidx);
			p.mat[idx.d,i] = p.mat[idx.d,i] + pxy;
		}
	}
	p.mat = p.mat/sum(p.mat);
	px = colSums(p.mat); pd = rowSums(p.mat);

	if( method == 'mi' ) {
		mi = 0;
		for( i in 1:nrow(p.mat) ) {
			for( j in 1:ncol(p.mat) ) {
				if( p.mat[i,j] > 0 ) mi = mi + p.mat[i,j]*log(p.mat[i,j]/pd[i]/px[j]);
			}
		}
		R = mi;
	} else {
	}

	return( R );
}

# MI(x,y|z)
mi.opt3d = function(S,given=3) {
	tl = c(2:29,Inf);
	mi = rep(0,length(tl)-1);
	names(mi) = c(as.character(2:28),'29+');

	if( given == 1 ) {
		S = S[,c(3,4,5,6,1,2,7,8)];
	} else if( given == 2 ) {
		S = S[,c(1,2,5,6,3,4,7,8)];
	}

	smi = 0;
	for( k in 1:(length(tl)-1) ) {
		pz = opt3d.p.area(S,c(0,Inf,0,Inf,tl[k],tl[k+1]));
		for( i in 1:(length(tl)-1) ) {
			px = opt3d.p.area(S,c(tl[i],tl[i+1],0,Inf,tl[k],tl[k+1]));
			for( j in 1:(length(tl)-1) ) {
				py = opt3d.p.area(S,c(0,Inf,tl[j],tl[j+1],tl[k],tl[k+1]));
				pxy = opt3d.p.area(S,c(tl[i],tl[i+1],tl[j],tl[j+1],tl[k],tl[k+1]));
				if( px > 0 & py > 0 & pxy > 0 ) { mi[k] = mi[k] + pxy*log(pxy/px/py); }
			}
		}
		smi = smi + mi[k]*pz;
	}

	return( list(mi=mi,smi=smi) );
}

# MI(x,y)
mi.2d = function(D,cutoff=0.001) {
	D[is.na(D)] = 0;
	D = D/sum(D); 
	D = D-cutoff;
	D[D<0] = 0;
	D = D/sum(D);
	mi = 0;
	for( i in 1:nrow(D) ) {
		for( j in 1:ncol(D) ) {
			px = sum(D[i,]); py = sum(D[,j]); pxy = D[i,j];
			if( px > 0 & py > 0 & pxy > 0 ) { mi = mi + pxy*log(pxy/px/py); }
		}
	}
	return( mi );
}

# COR(x,y)
cor.2d = function(D,N=10000) {
	D[is.na(D)] = 0;
	D = round(D/sum(D)*N);
	N = sum(D);
	V = matrix(rep(0,N*2),ncol=2); ii = 1;
	for( i in 1:nrow(D) ) {
		for( j in 1:ncol(D) ) {
			if( D[i,j] > 0 ) {
				n = D[i,j];
				idx = ii:(ii+n-1);
				V[idx,] = c(i,j);
				ii = ii+n
			}
		}
	}
	return( c(cor(V[,1],V[,2]),cor(V[,1],V[,2],method='spearman')) );
}

############################################################################
# Ordering
############################################################################

# 3-var sequence
get.prob.seq3d = function(SS,cnx,even.equal=FALSE) {
        tl = c(2:29);

	ss = SS; nn = cnx;
	Ppre0 = opt3d.p.area(ss,c(29,Inf,29,Inf,29,Inf));
	names(Ppre0) = sprintf('x%sx%sx%s',nn[1],nn[2],nn[3]);

	Ppre1 = rep(0,3);
	for( ii in 1:3 ) {
                if( ii == 1 ) { ss = SS[,c(1,2,3,4,5,6,7,8)]; nn = cnx[c(1,2,3)]; }
                if( ii == 2 ) { ss = SS[,c(3,4,1,2,5,6,7,8)]; nn = cnx[c(2,1,3)]; }
                if( ii == 3 ) { ss = SS[,c(5,6,1,2,3,4,7,8)]; nn = cnx[c(3,1,2)]; }
		Ppre1[ii] = opt3d.p.area(ss,c(2,29,29,Inf,29,Inf));
                names(Ppre1)[ii] = sprintf('%sx%sx%s',nn[1],nn[2],nn[3]);
	}

	Ppre2 = rep(0,6);
	for( ii in 1:6 ) {
                if( ii == 1 ) { ss = SS[,c(1,2,3,4,5,6,7,8)]; nn = cnx[c(1,2,3)]; }
                if( ii == 2 ) { ss = SS[,c(1,2,5,6,3,4,7,8)]; nn = cnx[c(1,3,2)]; }
                if( ii == 3 ) { ss = SS[,c(3,4,1,2,5,6,7,8)]; nn = cnx[c(2,1,3)]; }
                if( ii == 4 ) { ss = SS[,c(3,4,5,6,1,2,7,8)]; nn = cnx[c(2,3,1)]; }
                if( ii == 5 ) { ss = SS[,c(5,6,1,2,3,4,7,8)]; nn = cnx[c(3,1,2)]; }
                if( ii == 6 ) { ss = SS[,c(5,6,3,4,1,2,7,8)]; nn = cnx[c(3,2,1)]; }

                for( i in 1:(length(tl)-1) ) {
                        for( j in i:(length(tl)-1) ) {
				p = opt3d.p.area(ss,c(tl[i],tl[i+1],tl[j],tl[j+1],29,Inf));
				if( even.equal ) {
                                        if( i==j ) Ppre2[ii] = Ppre2[ii] + p/2;
                                        if( i<j ) Ppre2[ii] = Ppre2[ii] + p;
				} else {
					Ppre2[ii] = Ppre2[ii] + p;
				}
                        }
                }
                names(Ppre2)[ii] = sprintf('%s-%sx%s',nn[1],nn[2],nn[3]);
	}

        Ppre3 = rep(0,6);
        for( ii in 1:6 ) {
                if( ii == 1 ) { ss = SS[,c(1,2,3,4,5,6,7,8)]; nn = cnx[c(1,2,3)]; }
                if( ii == 2 ) { ss = SS[,c(1,2,5,6,3,4,7,8)]; nn = cnx[c(1,3,2)]; }
                if( ii == 3 ) { ss = SS[,c(3,4,1,2,5,6,7,8)]; nn = cnx[c(2,1,3)]; }
                if( ii == 4 ) { ss = SS[,c(3,4,5,6,1,2,7,8)]; nn = cnx[c(2,3,1)]; }
                if( ii == 5 ) { ss = SS[,c(5,6,1,2,3,4,7,8)]; nn = cnx[c(3,1,2)]; }
                if( ii == 6 ) { ss = SS[,c(5,6,3,4,1,2,7,8)]; nn = cnx[c(3,2,1)]; }

                for( i in 1:(length(tl)-1) ) {
                        for( j in i:(length(tl)-1) ) {
                                for( k in j:(length(tl)-1) ) {
                                        p = opt3d.p.area(ss,c(tl[i],tl[i+1],tl[j],tl[j+1],tl[k],tl[k+1]));
					if( even.equal ) {
                                        	if( i==j & j==k ) Ppre3[ii] = Ppre3[ii] + p/6;
                                        	if( i==j & j<k ) Ppre3[ii] = Ppre3[ii] + p/2;
                                        	if( i<j & j==k ) Ppre3[ii] = Ppre3[ii] + p/2;
                                        	if( i<j & j<k ) Ppre3[ii] = Ppre3[ii] + p;
					} else {
						Ppre3[ii] = Ppre3[ii] + p;
					}
                                }
                        }
                }
                names(Ppre3)[ii] = sprintf('%s-%s-%s',nn[1],nn[2],nn[3]);
        }

	Ppre = c(Ppre0,Ppre1,Ppre2,Ppre3);

        return( Ppre );
}

# E(z-y) according to E(y-x)
dist.opt3d = function(S,choice=123,cnx=c("","",""),method='seq') {

	if( choice == 123 ) { S = S[,c(1,2,3,4,5,6,7,8)]; nn = cnx[c(1,2,3)]; }
	if( choice == 132 ) { S = S[,c(1,2,5,6,3,4,7,8)]; nn = cnx[c(1,3,2)]; }
	if( choice == 213 ) { S = S[,c(3,4,1,2,5,6,7,8)]; nn = cnx[c(2,1,3)]; }
	if( choice == 231 ) { S = S[,c(3,4,5,6,1,2,7,8)]; nn = cnx[c(2,3,1)]; }
	if( choice == 312 ) { S = S[,c(5,6,1,2,3,4,7,8)]; nn = cnx[c(3,1,2)]; }
	if( choice == 321 ) { S = S[,c(5,6,3,4,1,2,7,8)]; nn = cnx[c(3,2,1)]; }

	tl = 2:29;
	dd = matrix(rep(0,53*53),ncol=53);
	for( i in 1:(length(tl)-1) ) {
		for( j in 1:(length(tl)-1) ) {
			for( k in 1:(length(tl)-1) ) {
				if( method == 'seq' ) {
					d1 = tl[j]-tl[i]; idx1 = d1+27;
					d2 = tl[k]-tl[j]; idx2 = d2+27;
				} else {
					d1 = tl[j]-tl[i]; idx1 = d1+27;
					d2 = tl[k]-tl[i]; idx2 = d2+27;
				}
				p = opt3d.p.area(S,c(tl[i],tl[i+1],tl[j],tl[j+1],tl[k],tl[k+1]));
				dd[idx1,idx2] = dd[idx1,idx2] + p;
			}
		}
	}

	dist.all = -26:26; dyz = pyz = rep(nrow(dd),0);
	for( i in 1:nrow(dd) ) {
		pyz[i] = sum(dd[i,]);
		dyz[i] = sum(dist.all*dd[i,]/sum(dd[i,]));
	}

	dist.all = 0:26; dyz2 = pyz2 = rep(nrow(dd),0);
	for( i in 27:nrow(dd) ) {
		pyz2[i] = sum(dd[i,27:53]);
		dyz2[i] = sum(dist.all*dd[i,27:53]/sum(dd[i,27:53]));
	}

	return( list(d=dyz,p=pyz,dd=dd,d2=dyz2,p2=pyz2) );
}

# get all potential sequences of N objects
get.allseq = function(N) {

	# prepare
	nN = length(N);
	R = T = NULL; 
	SEQ.ALL = NULL;
	for( i in 1:nN ) {
		seq.sub = rep(rep(1:nN,each=nN^(nN-i)),nN^(i-1));
		SEQ.ALL = cbind(SEQ.ALL,seq.sub);
	}
	nSEQ.ALL = apply(SEQ.ALL,1,function(x){return(length(unique(x)));});
	SEQ.ALL = SEQ.ALL[nSEQ.ALL==nN,];

	# no event
	r = rep("",nN); t = rep("",nN);
	for( i in 1:nN ) { r[i] = N[i]; t[i] = 'x'; } 
	R = rbind(R,r); T = rbind(T,t);

	for( ii in 1:nN ) {
		# there are nN-ii missing events
		n.miss = nN-ii;	
		ss = SEQ.ALL[,1:n.miss,drop=FALSE];
		nss = apply(ss,1,function(x){return(length(unique(x)));});
		ss = ss[nss==n.miss,,drop=FALSE];
		sss = rowSums(2^ss); usss = unique(sss);
		missed.all = ss[match(usss,sss),,drop=FALSE];

		# for not missed
		seq.all = NULL;
		for( i in 1:ii ) {
			seq.sub = rep(rep(1:ii,each=ii^(ii-i)),ii^(i-1));
			seq.all = cbind(seq.all,seq.sub);
		}
		nseq.all = apply(seq.all,1,function(x){return(length(unique(x)));});
		seq.all = seq.all[nseq.all==ii,,drop=FALSE];

		if( n.miss > 0 ) {
			for( jj in 1:nrow(missed.all) ) {
				missed = missed.all[jj,];
				not.missed = setdiff(1:nN,missed);
				for( kk in 1:nrow(seq.all) ) {
					r = N[not.missed[seq.all[kk,1]]]; t = "-";
					if( ncol(seq.all)>1 ) {
						for( jj in 2:ncol(seq.all) ) { r = c(r,N[not.missed[seq.all[kk,jj]]]); t = c(t,"-"); }
					}
					for( jj in 1:length(missed) ) { r = c(r,N[missed[jj]]); t = c(t,"x"); }
					R = rbind(R,r); T = rbind(T,t);
				}
			}
		} else {
			not.missed = 1:nN;
			for( kk in 1:nrow(seq.all) ) {
				r = N[not.missed[seq.all[kk,1]]]; t = "-";
				if( ncol(seq.all)>1 ) {
					for( jj in 2:ncol(seq.all) ) { r = c(r,N[not.missed[seq.all[kk,jj]]]); t = c(t,"-"); }
				}
				R = rbind(R,r); T = rbind(T,t);
			}
		}
	}

	RN = rep("",nrow(R));
	for( i in 1:length(RN) ) {
		r = R[i,1];
		if( T[i,1] == 'x' ) r = paste('x',R[i,1],sep='');
		for( j in 2:nN ) r = paste(r,T[i,j],R[i,j],sep='');
		RN[i] = r;
	}
	row.names(R) = row.names(T) = RN;

	return( list(R=R,T=T) );
}

# find all tripicate
find.all.triplicate = function(R,T) {
	O = NULL;
	for( i in 1:5 ) {
		for( j in (i+1):6 ) {
			for( k in (j+1):7 ) {
				if( T[i] == '-' ) {
					o = paste(R[i],T[j],R[j],T[k],R[k],sep='');
				} else {
					o = paste(T[i],R[i],T[j],R[j],T[k],R[k],sep='');
				}
				O = c(O,o);
			}
		}
	}
	return( O );
}

############################################################################
# Plotting Functions
############################################################################

# plot partition
my.opt2d.plot.partition = function(S,X,xlab="",ylab="",main="",xlim=NULL,ylim=NULL) {

	xlim = ylim = c(2,29);
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
my.opt2d.plot.marginal = function(S,X,xlab='Time',ylab='Probability',main='') {

	x = S[,1:2]; 
	r = range(x[is.finite(x)]);
	plot(survfit(X~1),conf.int=FALSE,xlim=c(1,30),ylim=c(0,1.3),xlab=xlab,ylab=ylab,yaxt='n',xaxt='n',main=main);
	axis(2,at=seq(0,1,by=0.2));
	axis(1,at=c(2,10,20,28,30),label=c('2','10','20','28','28+'));

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

# representive
plot_representive = function(S,X) {
	par(new=FALSE,bty='n');
	plot(0,0,xlim=c(10,11),xaxt='n',yaxt='n',xlab="",ylab="");
	par(fig=c(0.05,0.65,0.05,0.95),mar=c(4,4,2,2),bty='o',new=TRUE);
	my.opt2d.plot.partition(S,X,xlab=colnames(X)[1],ylab=colnames(X)[2]);
	par(fig=c(0.65,0.95,0.5,0.95),mar=c(2,2,2,1),bty='o',new=TRUE);
	my.opt2d.plot.marginal(opt2d.marginal(S,1),X[,1],main=colnames(X)[1]);
	par(fig=c(0.65,0.95,0.05,0.5),mar=c(2,2,2,1),bty='o',new=TRUE);
	my.opt2d.plot.marginal(opt2d.marginal(S,2),X[,2],main=colnames(X)[2]);
}

# 2d scatter plot
plot_2d_density = function(x,y,type='bubble',xlab="",ylab="",main="",scale=0.2,add=FALSE) {
	t = table(x,y);
	ux = sort(unique(x));	
	uy = sort(unique(y));	
	xx = rep(ux,length(uy));	
	yy = rep(uy,each=length(ux));
	if( type == 'bubble' ) {
		r = log2(sqrt(as.vector(t)/pi)+1)*scale;
		idx = r>1E-10;
		plot(0,0,xlim=c(1,53),ylim=c(1,53),col='white',xaxt='n',yaxt='n',xlab="",ylab="");
		points(xx[idx],yy[idx],cex=0.25*r[idx],pch=21,bg='black',col='black');
		lines(c(-100,100),c(-100,100),col='red',lwd=2,lty=2);
	} else if( type == 'ss' ) {
		smoothScatter(x,y,xlab=xlab,ylab=ylab,main=main);
	}
}

# plot event diagram
plot.event_diagram = function(score,med.t) {

	R = .4;
	nodes = as.character(names(med.t));
	par(mar=c(5,2,2,2));
	offset = c(-1.2,-0.5,0,0,0,0);

	# plot nodes
	m = (min(med.t)+max(med.t))/2; r = (max(med.t)-min(med.t))/2;
	x = med.t; y = r^2-(x-m)^2; y[y<0]=0; y = sqrt(y) + offset;
	plot(x,y,pch='.',col='white',ylim=c(min(y)-.3,max(y)+.3),xlim=c(2,11),bty='n',xaxt='n',yaxt='n',ylab='',xlab='Days Since Injury');
	axis(1,at=2:11);

	# circles
	symbols(x,y,circles=rep(R,length(x)),add=TRUE,inches=FALSE);
	text(x,y,nodes,cex=0.8);

	# arrows
	R = R*1.2;
	for( i in 1:nrow(score) ) {
		idx.from = which(nodes==as.character(score[i,1]));
		idx.to = which(nodes==as.character(score[i,2]));

		xf = x[idx.from]; xt = x[idx.to];
		yf = y[idx.from]; yt = y[idx.to];

		offset = 0.05; if( xt<xf ) offset = -offset;
		d = sqrt( (xt-xf)^2+(yt-yf)^2 );
		my.cos = (xt-xf)/d; my.sin = (yt-yf)/d;
		arrows( xf+R*my.cos+offset,yf+R*my.sin+offset,xt-R*my.cos+offset,yt-R*my.sin+offset,length=0.1,angle=20 );
	}

}



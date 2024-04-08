# library
library(survival);
source('cOPT/R/opt2d.R');
source('script99.misc.R');

# read data
RECOV.TIME = read.csv('../data/censored/floor_disease.time_sep.csv',header=TRUE);
RECOV.EVENT = read.csv('../data/censored/disease.event_sep.csv',header=TRUE);
X = NULL;
for( i in 1:7 ) {
	x = Surv(RECOV.TIME[,i+1],RECOV.EVENT[,i+1]);
	if( i == 1 ) { X = x; }
	else { X = cbind(as.data.frame(X),x); }
}
colnames(X) = colnames(RECOV.TIME)[2:8];
row.names(X) = RECOV.TIME[,1];
CNX = colnames(X);

# calculate pairwise distribution
if( TRUE ) {
	S2D.MAP = NULL;
	S2D.ALL = vector('list',1); ii = 1;
	for( i in 1:6 ) {
		for( j in (i+1):7 ) {
			A = c(1,52,1,52);
			S = opt2d(X[,c(i,j)],A,max.depth=12,max.n=1,int.part=TRUE);
			S2D.ALL[[ii]] = S; ii = ii+1;
			S2D.MAP = rbind(S2D.MAP,c(i,j));
		}
	}
	sS2D.MAP = matrix(CNX[S2D.MAP],ncol=2);
	save(S2D.ALL,S2D.MAP,file='../data/censored/MOF.RECOV_all.2D.RData');
}

# calculate triplicate-wise distribution
MAP = NULL;
for( i in 1:5 ) {
	for( j in (i+1):6 ) {
		for( k in (j+1):7 ) {
			MAP = rbind(MAP,c(i,j,k));
		}
	}
}


if( TRUE ) {
	load('../data/censored/MOF.RECOV_all.2D.RData');

	tl = c(1:52,Inf);
	for( ii in 1:nrow(S2D.MAP) ) {
		if( S2D.MAP[ii,2] != 7 ) { next; }
	
		S = NULL;
		for( i in 1:(length(tl)-1) ) {
			for( j in 1:(length(tl)-1) ) {
				if( j!=(length(tl)-1) & j<i ) { 
					p = 0;
				} else {
					p = opt2d.p.area(S2D.ALL[[ii]],c(tl[i],tl[i+1],tl[j],tl[j+1]));
				}
				S = rbind(S,c(tl[i],tl[i+1],tl[j],tl[j+1],0,p));
			}
		}
		S[,6] = S[,6]/sum(S[,6]);
		S[,5] = S[,6]/(S[,2]-S[,1])/(S[,4]-S[,3]);
	
		S2 = S2D.ALL[[ii]];
		for( i in 1:nrow(S2) ) S2[i,5] = opt2d.p.area(S,S2[i,1:4]);
		S2[,5] = S2[,6]/(S2[,2]-S2[,1])/(S2[,4]-S2[,3]);
		
		S2D.ALL[[ii]] = S2;
	}
	save(S2D.ALL,S2D.MAP,sS2D.MAP,file='../data/censored/MOF.RECOV_all.2D.Chop.RData');


}




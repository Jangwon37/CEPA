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
load('../data/censored/MOF.RECOV_all.2D.Chop.RData');

# pairwise precedence
PALL = NULL;
for( i in 1:nrow(S2D.MAP) ) {
        S = S2D.ALL[[i]]; p = rep(0,6);
        for( tl in 1:52 ) {
                p[1] = p[1] + opt2d.p.area(S,c(tl,tl+1,tl+1,Inf));
                p[2] = p[2] + opt2d.p.area(S,c(tl+1,Inf,tl,tl+1));
                p[3] = p[3] + opt2d.p.area(S,c(tl,tl+1,tl,tl+1));
                p[5] = p[5] + opt2d.p.area(S,c(tl,tl+1,tl+1,53));
                p[6] = p[6] + opt2d.p.area(S,c(tl+1,53,tl,tl+1));
        }
        p[4] = opt2d.p.area(S,c(53,Inf,53,Inf));
        PALL = rbind(PALL,p);
}

PRE = cbind(as.data.frame(sS2D.MAP[,1:2]),PALL);
colnames(PRE) = c('v1','v2','p1to2','p2to1','p1eq2','pnd','p1to2f','p2to1f');
PRE.RATIO1 = cbind(PRE[,1:2],PRE$p1to2f/(PRE$p1to2f+PRE$p2to1f)); colnames(PRE.RATIO1) = c('From','To','Pr');
PRE.RATIO2 = cbind(PRE[,c(2,1)],PRE$p2to1f/(PRE$p1to2f+PRE$p2to1f)); colnames(PRE.RATIO2) = c('From','To','Pr');
PRE.RATIO = rbind(PRE.RATIO1,PRE.RATIO2);

# plot figure4
setEPS();
# postscript(file='fig.pairwise_precedent.eps',width=8,height=8);
tiff(filename = '../result/censored/fig.pairwise_precedent_all.tiff', width = 2100, height = 2100, units = "px", res=300)
of = 0.05; dd = (1-2*of)/7;
par(fig=c(of,1-of,of,1-of),mar=c(0,0,0,0),bty='n',new=FALSE);
plot(0,0,xlim=c(11,11),col='white',xaxt='n',yaxt='n'); 
for( i in 1:7 ) {
	for( j in 1:7 ) {
		par(fig=c(of+dd*(i-1),of+dd*i,of+dd*(j-1),of+dd*j),mar=c(.1,.1,.1,.1),bty='o',new=TRUE);
		if( i == j ) { 
			plot(0,0,col='white',xaxt='n',yaxt='n'); 
			text(0,0,CNX[i],cex=1.75); 
		} else if( i > j ) {
			ii = which(PRE[,1]==CNX[i]&PRE[,2]==CNX[j]); ii1=7; ii2 = 8;
			if( length(ii) == 0 ) { ii = which(PRE[,2]==CNX[i]&PRE[,1]==CNX[j]); ii1=8; ii2=7; }
			p1 = PRE[ii,ii1]; p2 = PRE[ii,ii2]; 
			# if( i == 6 ) p1 = 0;
			ss = p1+p2; p1 = p1/ss; p2 = p2/ss;
			plot(0,0,col='white',xaxt='n',yaxt='n'); 
			text(0,0,sprintf('     %d%%\n%d%%   ',round(p1*100),round(p2*100)),cex=1.5); 
		} else {
			plot_2d_density(X[,i][,1],X[,j][,1],scale=1);
			if( i == 1 & j > 1 ) { axis(2,at=c(5,26,52),cex.axis=1); }
			if( i < 7 & j == 7 ) { axis(3,at=c(5,26,52),cex.axis=1); }
		}
	}
}
dev.off();


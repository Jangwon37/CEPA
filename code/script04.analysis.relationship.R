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
CNX.FULL = c('Infection','Eye','Respiratory','Digestive','Skin','Injury','Ear');


# median recovery
i = which(sS2D.MAP[,1]=='INFE' & sS2D.MAP[,2]=='RESP'); m1 = med.opt2d(S2D.ALL[[i]],given=2);
i = which(sS2D.MAP[,1]=='EYE' & sS2D.MAP[,2]=='RESP'); m2 = med.opt2d(S2D.ALL[[i]],given=2);
i = which(sS2D.MAP[,1]=='RESP' & sS2D.MAP[,2]=='DIGE'); m3 = med.opt2d(S2D.ALL[[i]],given=1);
i = which(sS2D.MAP[,1]=='RESP' & sS2D.MAP[,2]=='SKIN'); m4 = med.opt2d(S2D.ALL[[i]],given=1);
i = which(sS2D.MAP[,1]=='RESP' & sS2D.MAP[,2]=='INJU'); m5 = med.opt2d(S2D.ALL[[i]],given=1);
i = which(sS2D.MAP[,1]=='RESP' & sS2D.MAP[,2]=='EAR'); m6 = med.opt2d(S2D.ALL[[i]],given=1);
t = 1:52; m = cbind(m1[,2],m2[,2],m3[,2],m4[,2],m5[,2],m6[,2]); yy = NULL;
for( i in 1:6 ) {
	idx = t<52 & t>0 & is.finite(m[,i]);
	l = smooth.spline(t[idx],m[idx,i],df=6);
	yy = cbind(yy,predict(l,t)$y);
}


setEPS();
# postscript(file='fig.med_recov_since_resp_recov.eps',width=6,height=5);
tiff(filename = '../result/censored/fig.med_recov_since_resp_recov.tiff', width = 1400, height = 1000, units = "px", res = 300)
par(mar=c(4,4,2,2), cex.lab=0.8, cex.axis=0.8);
idx = t<21;
matplot(t[idx],m[idx,],cex=0,pch=1,ylab='Diagnosis Week Since Respiratory Diagnosis',xlab='Respiratory Diagnosis Week',col=c(1,2,3,4,5,6))
matlines(t[idx&t>0],yy[idx&t>0,],lty=1,lwd=2,col=c(1,2,3,4,5,6))
legend('topright',legend=c('Infection','Eye','Digestive','Skin','Injury','Ear'),lty=1,lwd=2,pch=1,col=c(1,2,3,4,5,6),bty='n', cex=0.7);
dev.off();


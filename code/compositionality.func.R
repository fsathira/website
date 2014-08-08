# compositionality.func.R
# =======================
# ReBoot is implemented in get.renorm.null.pval() when the options bootstrap=TRUE 
# and bootrenorm=TRUE (default).  We can go over the details of it when we meet.  
# There are lots of other functions that are there for historical reasons, 
# please ignore them.

library(gplots)
`%+%` = function(x,y) { paste(x,y,sep="") }

pos = function(x) {
	x[x<0] = 0
	return(x)
}

simulate.bug = function(type=c(1,2,3,4), noise=10, N=250, M=750, normalize=FALSE, baseline=2) {
	type = type[1]
	mm = data.frame(b1 = pos(1:N+10 + rnorm(N, 0, noise)))
	if (type==1) {
		mm$b2 = pos(rep(N/2, N) + rnorm(N, 0, noise))
		mm$b3 = pos((N:1+110)/30 + rnorm(N, 0, noise))
	} else if (type==2) {
		mm$b2 = pos((N:1+10)/3 + rnorm(N, 0, noise))
		mm$b3 = pos((110+1:N)/30 + rnorm(N, 0, noise))
	} else if (type==3) {
		mm$b2 = pos(rep(N/6, N) + rnorm(N, 0, noise))
		mm$b3 = pos((N:1+110)/30 + rnorm(N, 0, noise))
	} else if (type==4) {
		mm$b2 = pos(rep(20, N) + rnorm(N, 0, noise))
		mm$b3 = pos((110+1:N)/30 + rnorm(N, 0, noise))	
	} else if (type==5) {
		mm$b1 = pos(rep(N/2 + 10, N) + rnorm(N, 0, noise))
		mm$b2 = pos(rep(N/2, N) + rnorm(N, 0, noise))
		mm$b3 = pos(rep(10, N) + rnorm(N, 0, noise))
	} else if (type==6) {
		mm$b1 = pos((1:N)/7+10 + rnorm(N, 0, noise))
		mm$b2 = pos(rep(N/10 + 10, N) + rnorm(N, 0, noise))
		#mm$b3 = pos(rep(N/10 - 10, N) + rnorm(N, 0, noise))
		mm$b3 = pos((N:1+110)/10 + rnorm(N, 0, noise))
		mm$b4 = pos(rep(N/10, N) + rnorm(N, 0, noise))
	}
	if (type!=6) {
		mm$b4 = pos(rep(10, N) + rnorm(N, 0, noise))
	}
	if (M > 4) {
		for (i in (ncol(mm)+1):M) {
			mm[,i] = pos(rep(baseline, N) + rnorm(N, 0, noise))
		}
	}
	if (normalize) {
		return(normalize(mm))
	}
	return (mm)
}

# rankby is the column name to use as anchor
plot.bug = function(mm, raw=TRUE, plot.total=FALSE, write.file=TRUE, filename=NULL, rankby=NULL, decreasing=FALSE, main=NULL, legend=TRUE, ...) {
	if (raw) {
		if (is.null(main)) main = "Raw counts" 
		ylab = "Raw counts"
		if (is.null(filename)) filename = "composition.data.png"
	} else { 
		if (is.null(main)) main = "Relative abundance"
		ylab = "Relative abundance"
		if (is.null(filename)) filename = "composition.norm.data.png"
	}
	if (write.file) png(file=filename, width=1000, height=1000, pointsize=20)
	if (!is.null(rankby)) mm = mm[order(mm[,rankby],decreasing=decreasing),]
	if (plot.total) {
		total = apply(mm,1,sum)
		plot(total, col="gray", type='l', main=main, ylab=ylab, xlab="Subject Index", ylim=c(min(mm),max(total)), ...)
		lines(mm[,1], col=2)
	} else {
		if (raw) {
			plot(mm[,1], col=2, type='l', main=main, ylab=ylab, xlab="Subject Index", ylim=c(min(mm,na.rm=TRUE),max(mm,na.rm=TRUE)), ...)
		} else {
			plot(mm[,1], col=2, type='l', main=main, ylab=ylab, xlab="Subject Index", ylim=c(0,1), ...)
		}
	}
	lines(mm[,2], col=3, ...)
	lines(mm[,3], col=4, ...)
	lines(mm[,4], col=5, ...)
	if (legend) smartlegend("left", "top", colnames(mm)[1:4], col=2:5, lty=c(1,1,1,1))
	if (write.file) dev.off()
}

get.pval = function(x, y, lower.tail=TRUE, N.rand=10000, plot=FALSE, bootstrap=FALSE, ...) {
	this.cor = if (bootstrap) 0 else cor(x, y, use="complete.obs")
	if (bootstrap) lower.tail = !lower.tail
	rand.cor = rep(NA, N.rand)
	for (i in 1:N.rand) {
		if (bootstrap) {
			rand.idx = sample(1:length(x),replace=TRUE)
			rand.cor[i] = cor(x[rand.idx], y[rand.idx], use="complete.obs")
		} else {
			rand = sample(x, length(x))
			rand.cor[i] = cor(rand, y, use="complete.obs")
		}
	}
	if (lower.tail) {
		pval = (sum(this.cor > rand.cor) / N.rand)
	} else {
		pval = (sum(this.cor < rand.cor) / N.rand)
	}
	z.pval = pnorm(this.cor, mean=mean(rand.cor), sd=sd(rand.cor), lower.tail=lower.tail)
	if (plot) {
		hist(rand.cor, breaks="FD", xlim=c(-1,1), ...)
		abline(v=cor(x,y,use="complete.obs"), col="red")
		#smartlegend("left","top",c("sd="%+%sd(rand.cor),"pval="%+%pval,"z.pval="%+%z.pval),bty="n")
		smartlegend("left","top",c("sd="%+%sd(rand.cor),"pval="%+%pval),bty="n")
	}
	return(pval)
}

get.correlation = function(pp, method="variation.acomp.boot", filename=NULL) {
	cc = cor(pp[,1:4],use="complete.obs")
	bugs = colnames(pp)[1:4] # c("M1","M2","m1","m2")
	if (grepl("acomp",method)) {
		app = acomp(pp)
		acc = cor(app[,1:4],use="complete.obs")
	}
	if (is.null(filename)) filename = "composition.randhist." %+% method %+% ".png"
	png(file=filename, width=2000, height=1000, pointsize=20)
	par(mfrow=c(3,3))
	ss = matrix(NA, 4, 4, dimnames=list(bugs,bugs))
	ii = jj = 1
	for (i in bugs) {
		jj = 1
		for (j in bugs) {
			if (ii < jj) {
				print (i %+% "-" %+% j)
				if (method == "rand.boot") {
					ss[i,j] = ss[j,i] = get.pval(pp[,i], pp[,j], lower.tail=(cc[i,j]<0), plot=TRUE, bootstrap=TRUE, main=i %+% "-" %+% j)
				} else if (method == "rand.permute") {
					ss[i,j] = ss[j,i] = get.pval(pp[,i], pp[,j], lower.tail=(cc[i,j]<0), plot=TRUE, bootstrap=FALSE, main=i %+% "-" %+% j)
				} else if (method == "renorm.raw.permute") {
					ss[i,j] = ss[j,i] = get.renorm.pval(i, j, pp, N.rand=1000, lower.tail=(cc[i,j]<0), plot=TRUE, main=i %+% "-" %+% j)
				} else if (method == "renorm.boot") {
					ss[i,j] = ss[j,i] = get.renorm.pval(i, j, pp, N.rand=1000, lower.tail=(cc[i,j]<0), plot=TRUE, bootstrap=TRUE, main=i %+% "-" %+% j)
				} else if (method == "acomp.boot") {
					ss[i,j] = ss[j,i] = get.pval(app[,i], app[,j], lower.tail=(acc[i,j]<0), plot=TRUE, bootstrap=TRUE, main=i %+% "-" %+% j)
				} else if (method == "variation.acomp.boot") {
					tmp = 0 # this is a place holder for doing nothing
				} else if (method == "renorm.perm.boot") { # permute + renormalize the null, bootstrap the observed
					ss[i,j] = ss[j,i] = get.renorm.null.pval(i, j, pp, N.rand=1000, lower.tail=(cc[i,j]<0), plot=TRUE, bootstrap=TRUE, bootrenorm=FALSE, main=i %+% "-" %+% j)
				} else if (method == "renorm.perm.renorm.boot") { # permute + renormalize the null, bootstrap + renormalize the observed
					ss[i,j] = ss[j,i] = get.renorm.null.pval(i, j, pp, N.rand=1000, lower.tail=(cc[i,j]<0), plot=TRUE, bootstrap=TRUE, bootrenorm=TRUE, main=i %+% "-" %+% j)
				} else {
					dev.off()
					print ("Wrong method name")
					return()
				}
			} else {
				if (ii != 4 && jj != 1) plot(0,type='n', axes=F, xlab="", ylab="")
			}
			jj = jj + 1
		}
		ii = ii + 1
	}
	if (method == "variation.acomp.boot") {
		ss = get.var.pval(app, plot=TRUE)
	}
	dev.off()
}

normalize = function(data) {
	n = ncol(data)
	data.frame(diag(1/c(as.matrix(data)%*%rep(1,n))) %*% as.matrix(data))
}

perm.norm = function(bug1, bug2, data, bootstrap=FALSE) {
	n = ncol(data)
	m = nrow(data)
	if (bootstrap) {
		boot.idx = sample(1:m,replace=TRUE)
		data[,bug1] = data[boot.idx,bug1]
		data[,bug2] = data[boot.idx,bug2]
	} else {
		data[,bug1] = sample(data[,bug1], m)
		data[,bug2] = sample(data[,bug2], m)
	}
	renorm.data = normalize(data)
	names(renorm.data) = names(data)
	return(renorm.data)
}

get.renorm.pval = function(bug1, bug2, data, lower.tail=TRUE, N.rand=1000, bootstrap=FALSE, plot=FALSE, ...) {
	norm.data = normalize(data)
	this.cor = if (bootstrap) 0 else cor(norm.data[,bug1], norm.data[,bug2], use="complete.obs")
	if (bootstrap) lower.tail = !lower.tail
	rand.cor = rep(NA, N.rand)
	for (i in 1:N.rand) {
		rand.data = perm.norm(bug1, bug2, data, bootstrap=bootstrap)
		rand.cor[i] = cor(rand.data[,bug1], rand.data[,bug2], use="complete.obs")
	}
	if (lower.tail) {
		pval = (sum(this.cor > rand.cor) / N.rand)
	} else {
		pval = (sum(this.cor < rand.cor) / N.rand)
	}
	z.pval = pnorm(this.cor, mean=mean(rand.cor), sd=sd(rand.cor), lower.tail=lower.tail)
	if (plot) {
		hist(rand.cor, breaks="FD", xlim=c(-1,1), ...)
		abline(v=this.cor, col="red")
		smartlegend("left","top",c("sd="%+%sd(rand.cor),"pval="%+%pval,"z.pval="%+%z.pval),bty="n")
	}
	return(z.pval)
}

# bug1, bug2	names of the data column corresponding to microbes of interest
# data		data.frame of (relative) abundances
# lower.tail	boolean to identify which tail of the normal distribution to use to calculate p-value.
#			TRUE if null < boot, FALSE otherwise.
# N.rand		number of iteration
# bootstrap	degenerate, not used
# bootrenorm	if TRUE, the bootstrap distribution uses renormalized data
# plot		if TRUE, the null and bootstrap distributions are plotted.
get.renorm.null.pval = function(bug1, bug2, data, lower.tail=TRUE, N.rand=1000, bootstrap=FALSE, bootrenorm=TRUE, plot=FALSE, ...) {
	lower.tail = !lower.tail
	norm.data = normalize(data)
	null = rep(NA, N.rand)
	boot = rep(NA, N.rand)
	for (i in 1:N.rand) {
		null.data = perm.norm(bug1, bug2, data, bootstrap=FALSE)
		null[i] = cor(null.data[,bug1], null.data[,bug2], use="complete.obs")
		if (bootrenorm) {
			boot.data = perm.norm(bug1, bug2, data, bootstrap=TRUE)
			boot[i] = cor(boot.data[,bug1], boot.data[,bug2], use="complete.obs")
		} else {
			boot.idx = sample(1:nrow(data),replace=TRUE)
			boot[i] = cor(data[boot.idx,bug1], data[boot.idx,bug2], use="complete.obs")
		}
	}
	pval = t.test(null,boot,lower.tail=lower.tail)$p.value
	z.pval = pnorm(mean(null), mean=mean(boot), sd=sqrt((var(boot) + var(null))*0.5), lower.tail=lower.tail)
	if (plot) {
		hist(boot, breaks="FD", xlim=c(-1,1), ...)
		hist(null, breaks="FD", add=TRUE)
		abline(v=mean(boot), col="blue")
		abline(v=mean(null), col="red")
		smartlegend("right","top",c("Bootstrap SD: "%+%signif(sd(boot),5),"P-value: "%+%signif(z.pval,5)),bty="n")
	}
	return(z.pval)
}

## modify rawdata matrix to accommodate the higher level taxa
replace.taxa = function(taxon, fullname, rawdata) {
	bugs.idx = idBug(names(rawdata), taxon, findall=TRUE)
	newdata = rawdata[,!bugs.idx]
	newdata[,fullname] = apply(rawdata[,bugs.idx],1,sum)
	return(newdata)
}


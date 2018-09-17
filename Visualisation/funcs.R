
plot2DM <- function(wedge, name, depth, Ns, kill)
{
	wedge$Struct[wedge$Freq < kill] <- 0
	tst <- wedge$Struct/wedge$Freq
	tst[is.na(tst)] <- 0
	tst[tst > 1] <- 0

	plot(0,0, col=rgb(0,0,0,0), xlim=c(1,depth), ylim=c(0, 1), cex=1, xlab="l-mer size", 
		ylab="Mean K-mer Distinctness",
		main=paste(name, " l-mer Distinctness with ", Ns, " Ns", sep=""))

	colrs = c("black", "darkred", "darkblue", "purple")

	for(i in 0:max(wedge[,2]))
	{
		ts2 <- tst[wedge[,2] == i]
		xp <- wedge[wedge[,2] == i, 1]
		points(xp, ts2, col=colrs[i+1], type="l", lwd=2, lty=i+1)
	}

	legs <- list()
	legs[[1]] <- "N=0"
	legs[[2]] <- c("N=0", "N=1")
	legs[[3]] <- c("N=0", "N=1", "N=2")
	legs[[4]] <- c("N=0", "N=1", "N=2", "N=3")
	
	legend("topleft", bty="n", col=colrs[1:(1+max(wedge[,2]))], lty=c(1:(1+max(wedge[,2]))), lwd=c(rep(2, 1+max(wedge[,2]))), 
		legend=legs[[max(wedge[,2]) + 1]])
}

plot_line_fams <- function(wedge, name, depth, Ns, inc)
{
	wedge <- trimtable(wedge, depth)
	
	tst <- wedge[,4]
	tst[is.na(tst)] <- 0
	tst[tst > 1] <- 0
	wedge <- wedge[tst > 0,]
	cols <- colTransSD(wedge)
	ypos = findYpos(wedge)

	a2 <- makeRedSD(wedge)

	tst <- tst[tst > 0]
	
	plot(0,0, col=rgb(0,0,0,0), xlim=c(0,1), ylim=c(1, depth), ylab="K-size w/flattened seed lengths and N-mask", 
		xlab="Mean l-mer Distinctness",
		main=paste(name, " l-mer Distinctness with ", Ns, " Ns", sep=""))

	drawSLines_fam(wedge)
	
	symbols(tst, ypos,cex=1, circles=c(a2), inches=inc, fg=cols, bg=cols, add=TRUE)

	legend("bottomright", title="Color Coding", col=c("grey", "grey", "red", "green", "blue"), pch=c(21, 19, 19, 19, 19), cex=0.9,
		pt.cex=c(2,1,2,2,2),
		legend=c("(alpha) K-Norm Freq", "(scale) Inverse Weighted SD", "Absolute Frequency", "Left Seed", "Right Seed"), bty="n")	
 
}

plot_lineM <- function(wedge, name, depth, Ns, inc)
{
	wedge <- trimtable(wedge, depth)
	
	tst <- wedge[,4]/wedge[,6]
	tst[is.na(tst)] <- 0
	tst[tst > 1] <- 0
	wedge <- wedge[tst > 0,]
	cols <- colTrans(wedge)
	ypos = findYpos(wedge)

	tst <- tst[tst > 0]
	
	plot(0,0, col=rgb(0,0,0,0), xlim=c(0,1), ylim=c(1, depth), ylab="K-size w/flattened seed lengths and N-mask", 
		xlab="Mean l-mer Distinctness",
		main=paste(name, " l-mer Distinctness with ", Ns, " Ns", sep=""))

	drawSLines(wedge)
	
	symbols(tst, ypos,cex=1, circles=c(rep(1,length(ypos))), inches=inc, fg=cols, bg=cols, add=TRUE)

	legend("bottomright", title="Color Coding", col=c("grey", "red", "green", "blue"), pch=c("*", "*", "*", "*"), cex=0.9,
		pt.cex=c(2,2,2),
		legend=c("(alpha) K-Norm Freq", "Absolute Frequency", "Left Seed", "Right Seed"), bty="n")	
 
}
drawSLines_fam <- function(x)
{
	for(i in 0:max(x[,2]))
	{
		for(j in 0:max(x[,3]))
		{
			set2draw <- x[(x[,1] - x[,2]) == i & x[,3] == j,]
			if(nrow(set2draw) > 0)
			{
				tst <- set2draw[,4]
				tst[is.na(tst)] <- 0
				set2draw <- set2draw[tst > 0,]
				tst <- tst[tst > 0]
				if(length(tst) > 0)
				{
					if(i > 1000)
					{
						prec <- x[x[,1] == i & (x[,1] - x[,2]) == j & x[,3] == max(c(j-1, 0)),]
						if(nrow(prec) > 0)
						{
							cont <- prec[,4]/prec[,6]
							if(is.na(cont))
							{
								cont <- 0
							}
							if(cont > 0)
							{
								tst <- c(cont, tst)
								set2draw <- rbind(prec, set2draw)
							}
						}
					}
					ypos <- findYpos(set2draw)
					ypos <- ypos[tst > 0.01]
					tst <- tst[tst > 0.01]
					points(tst, ypos, type="l", col=rgb(0.1,0.1,0.1,0.02))
				}
			}
		}
	}
}

drawSLines <- function(x)
{
	for(i in 0:max(x[,2]))
	{
		for(j in 0:max(x[,3]))
		{
			set2draw <- x[(x[,1] - x[,2]) == i & x[,3] == j,]
			if(nrow(set2draw) > 0)
			{
				tst <- set2draw[,4]/set2draw[,6]
				tst[is.na(tst)] <- 0
				set2draw <- set2draw[tst > 0,]
				tst <- tst[tst > 0]
				if(length(tst) > 0)
				{
					if(i > 1000)
					{
						prec <- x[x[,1] == i & (x[,1] - x[,2]) == j & x[,3] == max(c(j-1, 0)),]
						if(nrow(prec) > 0)
						{
							cont <- prec[,4]/prec[,6]
							if(is.na(cont))
							{
								cont <- 0
							}
							if(cont > 0)
							{
								tst <- c(cont, tst)
								set2draw <- rbind(prec, set2draw)
							}
						}
					}
					ypos <- findYpos(set2draw)
					ypos <- ypos[tst > 0.01]
					tst <- tst[tst > 0.01]
					points(tst, ypos, type="l", col=rgb(0.1,0.1,0.1,0.05))
				}
			}
		}
	}
}

plot_mats <- function(wedge, name, depth, Ns)
{
	wedge <- trimtable(wedge, depth)
	
	tst <- wedge[,4]/wedge[,6]
	tst[is.na(tst)] <- 0
	wedge <- wedge[tst > 0,]
	cols <- colTrans(wedge)
	ypos = findYpos(wedge)

	tst <- tst[tst > 0]
	
	plot(tst, ypos,pch="*", cex=1, col=cols, xlim=c(0,1), ylim=c(1,depth),
		ylab="K-size w/flattened seed lengths and N-mask", xlab="Mean K-mer Distinctness",
		main=paste(name, " DNA l-mer Distinctness with ", Ns, " Ns", sep=""))

	legend("topleft", title="Color Scheme", col=c("grey", "red", "green", "blue"), pch=c("*", "*", "*", "*"), cex=0.9,
		pt.cex=c(2,2,2),
		legend=c("(alpha) K-Norm Freq", "Absolute Frequency", "Left Seed", "Right Seed"), bty="n")	
 
}

plotExtra <- function(wedge, inc)
{
	wedge <- trimtable(wedge, max(wedge[,1]))
	wedge$Struct <- pmin(wedge$Struct, wedge$Freq)
	wedge$ShapeN[wedge$Shape == 0] <- 0
	wedge$Shape <- pmax(wedge$Shape, 1)

	pareto <- wedge$ShapeN/wedge$Shape

	pair <- cbind(wedge[,1], pareto)

	rs <- makeReds2(pair)
	
	symbols(wedge$Ns, wedge$Depth, bg=rgb(rs, rep(0.0, length(rs)), rs, rs), squares=(wedge$Struct/pmax(wedge$Freq,1)), inches=inc, xlab="N Count", ylab="K-size", main="N-Mask Scale", xlim=c(-0.5,3.5), ylim=c(0,max(wedge[,1])+1), bty="n")
	#legend("topleft")

	#pareto
}

makeReds2 <- function(x)
{
	xr <- vector()
	n <- max(x[,1])

	for(i in 1:n)
	{
		x1 <- x[x[,1] == i,]
		xb <- x1[,2]

		if(length(xb) > 0)
		{
			xb <- xb - min(xb)

			if(max(xb) > 0)
			{
				xb <- xb/max(xb)
			}
			else
			{
				xb <- rep(0, length(xb))
			}

			xr <- c(xr, xb)
		}
	}

	xr
}

makeReds <- function(x, n)
{
	xr <- vector()

	for(i in 1:n)
	{
		x1 <- x[x[,1] == i,]
		xb <- x1[,5]

		if(length(xb) > 0)
		{
			if(max(xb) > 0)
			{
				xb <- pmax(rep(0.1, length(xb)), xb/max(xb))
			}
			else
			{
				xb <- rep(0, length(xb))
			}

			xr <- c(xr, xb)
		}
	}

	xr
}

makeRedSD <- function(x, n)
{
	xr <- vector()

	xb <- x[,6]

	#xb <- -log( pmax(rep(0.1, length(xb)), xb/max(xb)) )

	xb <- 1-pmax(rep(0.02, length(xb)), xb/max(xb))

	xb
}


colTrans <- function(x)
{
	r <- pmax(0, (log2(x[,5]/(median(x[,5]))*2))+3)

	r <- pmin(r, 10)
	r <- r/10
	g <- (x[,1] - x[,2])/x[,1]
	b <- (x[,2] - x[,3])/x[,2]

	r[is.na(r)] <- 0
	g[is.na(g)] <- 0
	b[is.na(b)] <- 0

	a <- makeReds(x, 29)

	rgb(r,g,b,a)
}

colTransSD <- function(x)
{
	r <- pmax(0, (log2(x[,5]/(median(x[,5]))*2))+3)

	r <- pmin(r, 10)
	r <- r/10
	g <- (x[,1] - x[,2])/x[,1]
	b <- (x[,2] - x[,3])/x[,2]

	r[is.na(r)] <- 0
	g[is.na(g)] <- 0
	b[is.na(b)] <- 0

	a1 <- makeReds(x, 29)
	a2 <- makeRedSD(x,29)

	#a2 <- 1-a2

	rgb(r,g,b, a1)
}

findYpos <- function(x)
{
	x[,1] + (x[,2] * (0.5/x[,1])) + (x[,3] * (0.5/x[,1]))
}

trimtable <- function(x, d)
{
	x[x[,1] < d,]
	#x[x[,2] != 0,]
}

extractBase <- function(x)
{
	x[x[,2] == 1 & x[,3] == 0,]
}

extractSingle <- function(x)
{
	x[x[,3] == (x[,2] - 1),]
}

extractRest <- function(x)
{
	x[x[,3] < (x[,2] - 1) & x[,2] > 1,]
}

library(plot3D)

findZpos <- function(x)
{
	(x[,2] * (0.5/x[,1])) + (x[,3] * (0.5/x[,1])) * 10
}

plot_line_fam3D <- function(wedge, name, depth, Ns, inc)
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

	drawSLines_fam3D(wedge)
	
	symbols(tst, ypos,cex=1, circles=c(a2), inches=inc, fg=cols, bg=cols, add=TRUE)

	legend("bottomright", title="Color Coding", col=c("grey", "grey", "red", "green", "blue"), pch=c(21, 19, 19, 19, 19), cex=0.9,
		pt.cex=c(2,1,2,2,2),
		legend=c("(alpha) K-Norm Freq", "(scale) Inverse Weighted SD", "Absolute Frequency", "Left Seed", "Right Seed"), bty="n")	
 
}

plot_line_3D <- function(wedge, name, depth, Ns, inc)
{
	wedge <- trimtable(wedge, depth)
	
	tst <- wedge[,4]/wedge[,6]
	tst[is.na(tst)] <- 0
	tst[tst > 1] <- 0
	wedge <- wedge[tst > 0,]
	cols <- colTrans(wedge)
	ypos <- findYpos(wedge)
	zpos <- findZpos(wedge)

	tst <- tst[tst > 0]
	
	scatter3D(0, 0, 0, pch=19, cex=0.6, colvar=NULL, col=rgb(1,1,1,0), add=FALSE, zlim=c(0, max(zpos)), xlim=c(0,1), ylim=c(1, depth), ylab="K-size w/flattened seed lengths and N-mask", 
		xlab="Mean l-mer Distinctness", zlab="extention by depth-uniqueness",
		main=paste(name, " l-mer Distinctness with ", Ns, " Ns", sep=""))

	drawSLines3D(wedge)

	scatter3D(tst, ypos, zpos, pch=19, cex=0.6, colvar=NULL, col=cols, add=TRUE)

	#symbols(tst, ypos,cex=1, circles=c(rep(1,length(ypos))), inches=inc, fg=cols, bg=cols, add=TRUE)

	legend("bottomright", title="Color Coding", col=c("grey", "red", "green", "blue"), pch=c(19, 19, 19, 19), cex=0.9,
		pt.cex=c(2,2,2,2),
		legend=c("(alpha) K-Norm Freq", "Absolute Frequency", "Left Seed", "Right Seed"), bty="n")	
 
}

drawSLines_fam3D <- function(x)
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
					zpos <- findZpos(set2draw)
					ypos <- ypos[tst > 0.01]
					zpos <- zpos[tst > 0.01]
					tst <- tst[tst > 0.01]
					if(length(tst) > 0 & length(tst) == length(zpos) & length(zpos) == length(ypos))
					{
						lines3D(tst, ypos, zpos, col=rgb(0.3,0.3,0.3), add=TRUE)
					}
				}
			}
		}
	}
}

drawSLines3D <- function(x)
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
					zpos <- findZpos(set2draw)
					ypos <- ypos[tst > 0.01]
					zpos <- zpos[tst > 0.01]
					tst <- tst[tst > 0.01]
					if(length(tst) > 0 & length(tst) == length(zpos) & length(zpos) == length(ypos))
					{
						lines3D(tst, ypos, zpos, col=rgb(0.3,0.3,0.3), add=TRUE)
					}
				}
			}
		}
	}
}
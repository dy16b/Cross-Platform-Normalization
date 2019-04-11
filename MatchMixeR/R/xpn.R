######################################################################
#Copyright Jason Rudy & Faramarz Valafar 2009-2010

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
######################################################################

#Implementation of the xpn cross platform normalization algorithm.

#'@param X gene expression matrix
#'@param Y gene expression matrix
#'@export
xpn = function(platform1.data,platform2.data,K=10,L=4,p1.names=0,p2.names=0,gene.cluster="kmeans",assay.cluster="kmeans",corr="pearson", iterations=30, skip.match=FALSE){
	#If K or L is not a single value, it is taken as a list of possible values.

	#Match names
	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)
	x <- input$x
	y <- input$y
	input <- NA

	#Remove the medians
	x_demed = x - rowMedians(as.matrix(x))
	y_demed = y - rowMedians(as.matrix(y))

	#Get the dimensions
	nx=dim(x)[2]
	ny=dim(y)[2]
	mx=dim(x)[1]
	my=dim(y)[1]


	#Create the combined dataframe for clustering purposes
	combined <- cbind(x_demed,y_demed)

	#Detect K and L if necessary
	K = detect_K(combined,K,corr)
	L = detect_L(x_demed,y_demed,L,corr)


	xout = 0
	yout = 0


	for(iter in 1:iterations){
	cat("XPN iteration",iter,"out of",iterations,"\n")

	#Do the assay clustering
	assayclusters = xpnassaycluster(x_demed, y_demed, L, assay.cluster, corr)

	#Do the gene clustering
	geneclusters = xpngenecluster(combined, K, gene.cluster, corr)

	#Estimate the XPN parameters
	xparams = xpnparamsp(x,geneclusters,assayclusters[1:nx])
	yparams = xpnparamsp(y,geneclusters,assayclusters[(nx+1):(nx+ny)])

	#Calcuate the weighted averages
	nor = 1/(nx + ny)
	Aav = try((1/(xparams$nblock+yparams$nblock))*(xparams$nblock*xparams$A + yparams$nblock*yparams$A))
	bav = nor*(nx*xparams$b + ny*yparams$b)
	cav = nor*(nx*xparams$c + ny*yparams$c)
	sigmaav = sqrt(nor*(nx*xparams$s2 + ny*yparams$s2))
	sigmax = sqrt(xparams$s2)
	sigmay = sqrt(yparams$s2)


	#Calculate the expanded A
	expAx = xpnclusterexpand(xparams$A,geneclusters,assayclusters[1:nx])
	expAy = xpnclusterexpand(yparams$A,geneclusters,assayclusters[(nx+1):(nx+ny)])

	#Calculate the residuals
	epsilonx = as.matrix(x) - (as.vector(xparams$b) * expAx + kronecker(ones(1,nx),xparams$c))
	epsilony = as.matrix(y) - (as.vector(yparams$b) * expAy + kronecker(ones(1,ny),yparams$c))

	#Calculate the expanded average A
	expAavx = xpnclusterexpand(Aav,geneclusters,assayclusters[1:nx])
	expAavy = xpnclusterexpand(Aav,geneclusters,assayclusters[(nx+1):(nx+ny)])

	#Calculate the output values
	xout = xout + (1/iterations)*((as.vector(bav) * expAavx) + kronecker(ones(1,nx),cav) + as.vector(sigmaav/sigmax) * epsilonx)
	yout = yout + (1/iterations)*((as.vector(bav) * expAavy) + kronecker(ones(1,ny),cav) + as.vector(sigmaav/sigmay) * epsilony)

	}#end of the enclosing for loop

	#Put the rownames back in and convert to data frames
	xout = as.data.frame(xout,row.names=rownames(x))
	yout = as.data.frame(yout,row.names=rownames(y))

	#All done!
	return(list(x=xout,y=yout))
}

xpnparamsp = function(x, geneclusters, assayclusters, method){
	#x should be a dataframe or hashdataframe

	#Get K and L
	K=max(geneclusters)
	L=max(assayclusters)

	#Get number of assays and genes
	numassays=length(x)
	numgenes=length(rownames(x))

	#Number of assays in each cluster
	nj = matrix(table(as.factor(assayclusters)),1,L)

	#Set up the output variables
	A = matrix(nrow=K,ncol=L)
	nblock = matrix(nrow=K,ncol=L)
	b = matrix(nrow=numgenes, ncol=1)
	c = matrix(nrow=numgenes, ncol=1)
	sigma = matrix(nrow=numgenes, ncol=1)
	mu = matrix(nrow=numgenes,ncol=L)


	#For each gene cluster, estimate the XPN parameters
	for (a in 1:K){

		#Get the logical indexing vector for gene cluster a
		geneinds = geneclusters==a
		numgenesa = sum(as.numeric(geneinds))

		#Get the data for this gene cluster
		xa = as.matrix(x[geneinds,])

		#Get the number of genes in this cluster
		Ga = dim(xa)[1]
		S = dim(xa)[2]

		#For each assay cluster and each gene, get the (relatively) unconstrained
		#MLE for mu
		mua = matrix(nrow=Ga,ncol=L)
		for (bb in 1:L){
			assayinds = assayclusters==bb
			xab = matrix(xa[,assayinds],numgenesa,sum(as.numeric(assayinds)))
			mua[,bb]=t(t(rowMeans(xab)))
		}


		#Calculate sigma2
		expmua = xpnclusterexpand(mua, colclusters = assayclusters)
		sigma2 = xa - expmua
		sigma2 = sigma2*sigma2
		sigma2 = xpnclustercollapse(sigma2,colclusters = assayclusters)

		solna = XPN_MLE_C(mua,sigma2,nj)

		ba = solna$b
		Aa = solna$A
		ca = solna$c
		s2a = solna$s2

		#Write the values for this gene group into the larger arrays
		sigma[geneinds,] = s2a
		mu[geneinds,] = mua
		A[a,] = Aa
		nblock[a,] = nj
		b[geneinds,] = ba
		c[geneinds,] = ca
	}

	return(list(A=A,nblock=nblock,b=b,c=c,s2=sigma,mu=mu))

}


xpnclustercollapse = function(aMatrix, rowclusters = 1:(dim(aMatrix)[1]), colclusters = 1:(dim(aMatrix)[2])){
	#This undoes the work of xpnclusterexpand

	l = dim(aMatrix)[1]
	w = dim(aMatrix)[2]

	#collapse the rows
	rowidx = split(1:l,rowclusters)
	lc = length(rowidx)
	collapse1 = matrix(NA,lc,w)
	for ( i in 1:lc){
		collapse1[i,] =  colMeans(matrix(aMatrix[rowidx[[i]],],length(rowidx[[i]]),w))
	}

	#collapse the columns
	colidx = split(1:w,colclusters)
	wc = length(colidx)
	collapse2 = matrix(NA,lc,wc)
	for ( i in 1:wc){
		collapse2[,i] =  rowMeans(matrix(collapse1[,colidx[[i]]],lc,length(colidx[[i]])))
	}

	return(collapse2)

}


xpnclusterexpand = function(aMatrix, rowclusters = 1:(dim(aMatrix)[1]), colclusters = 1:(dim(aMatrix)[2]), noise = 0){
	#This is basically a fancy repmat.  rowclusters and colclusters
	#as produced by pam.  aMatrix must have the same number of rows
	#and columns as there are row and column clusters.

	l = dim(aMatrix)[1]
	w = dim(aMatrix)[2]

	#Get the dimensions of the output
	m = length(rowclusters)
	n = length(colclusters)

	#Get the number of row and column clusters
	nrowclust = max(rowclusters)
	ncolclust = max(colclusters)

	#initialize the vertically expanded matrix
	vert = zeros(m,dim(aMatrix)[2])

	#Do the vertical expansion
	for (i in 1:m){
		vert[i,]=aMatrix[rowclusters[i],] + (noise*randn(1,w))
	}

	#initialize the fully expanded matrix
	output = matrix(NA,m,n)

	#Do the horizontal expansion
	for (j in 1:n){
		output[,j] = vert[,colclusters[j]] + (noise*randn(m,1))
	}

	#That is it!
	return(output)
}

xpnassaycluster = function(x,y,L,cluster="kmeans",corr="pearson"){

	if(is.numeric(cluster) | is.factor(cluster)){
		return(as.numeric(cluster))
	}

	#The numbers of assays
	nx = dim(x)[2]
	ny = dim(y)[2]

	#The number of genes
	m = dim(x)[1]

	#Compute the two correlation matrices if necessary
	if (cluster == 'pam' | (length(L)>1)){
		xdiss = 1 - cor(as.matrix(x), method = corr)
		ydiss = 1 - cor(as.matrix(y), method = corr)
	}

	#Determine L if a range has been given
	if(length(L)>1){
		Lx=pamk(xdiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
		Ly = pamk(ydiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
		L = max(Lx,Ly)
		cat(as.character(L), "assay clusters detected...\n")
	}

	#Do the clustering
	if (cluster == "classic"){

		done = FALSE
		while (!done){
			xyclust <- ckmeans(t(cbind(x,y)),centers=L, iter.max=20, distance=corr)
			#print(xyclust$cluster)
			xclust <- xyclust$cluster[1:nx]
			yclust <- xyclust$cluster[(nx+1):(nx+ny)]
			#print(xclust)
#			print(yclust)
#			print(nlevels(as.factor(xclust)))
#			print(nlevels(as.factor(yclust)))
			done <- nlevels(as.factor(xclust)) == L & nlevels(as.factor(yclust)) == L
			if(!done){
				cat("Not all platforms contained all assay clusters.  Trying again (only \"classic\" mode has this problem)...\n")
			}
		}
		return(xyclust$cluster)
	}
	else if (cluster == "pam"){
		xclust = pam(xdiss, k=L, diss=TRUE, keep.diss=FALSE, keep.data=FALSE,cluster.only=TRUE)
		yclust = pam(ydiss, k=L, diss=TRUE, keep.diss=FALSE, keep.data=FALSE,cluster.only=TRUE)
	}else if (cluster == "kmeans"){

		#This loop ensures there are no empty clusters.
		done = FALSE
		while (!done){
			xclust = ckmeans(t(x), centers=L, iter.max = 20, distance = corr)
			yclust = ckmeans(t(y), centers=L, iter.max = 20, distance = corr)

			done = !xclust$empty & !yclust$empty
		}
		xclust = xclust$cluster
		yclust = yclust$cluster

	} else if (cluster == "flexclust"){

		distCor1 = function (x, centers)
		{
    		z <- matrix(0, nrow(x), ncol = nrow(centers))
    		for (k in 1:nrow(centers)) {
        		z[, k] <- 1 - cor(t(x), centers[k, ], method=corr)
    		}
    		z
		}
		#This loop ensures there are no empty clusters.
		done = FALSE
		while (!done){
		xclust = kcca(t(x), k=L, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE)
		yclust = kcca(t(y), k=L, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE)

			done = (dim(xclust@clusinfo)[1] == L) &(dim(yclust@clusinfo)[1] == L)
		}
		xclust = xclust@cluster
		yclust = yclust@cluster
	}else if(cluster == "random"){
			xclust = sample(c(1:L,sample(1:L,nx-L,replace = TRUE)),nx,replace=FALSE)
			yclust = sample(c(1:L,sample(1:L,ny-L,replace = TRUE)),ny,replace=FALSE)
	}

	#Compute cluster averages
	xave = matrix(NA,m,L)
	yave = matrix(NA,m,L)
	for (i in 1:L){
		xinds = xclust==i
		yinds = yclust==i
		xave[,i] = rowMeans(as.matrix(x[,xinds],nrow=m,ncol=sum(as.numeric(xinds))))
		yave[,i] = rowMeans(as.matrix(y[,yinds],nrow=m,ncol=sum(as.numeric(yinds))))
	}

	#Compute the cluster correlation matrix
	clustercor = cor(xave,yave, method = corr)

	#Map clusters
	xtoymap = matrix(NA,L,1)
	for (i in 1:L){
		highest = which.max(clustercor)
		xind = highest%%L

		yind = highest%/%L + 1

		if (xind==0){
			xind = L
			yind = yind-1
		}
		xtoymap[xind]=yind
		clustercor[xind,] = -2
		clustercor[,yind] = -2
	}

	#Change the x clusters using the map
	newxclust = xclust
	for (i in 1:L){
		newxclust[xclust==i] = xtoymap[i]
	}
	xclust = newxclust

	#Return the clustering in a combined vector
	return(c(xclust,yclust))
}


xpngenecluster = function(x,K,cluster="pam",corr="pearson"){

	#Dimensions
	m = dim(x)[1]
	n = dim(x)[2]

	#Compute the dissimilarity matrix
	if (cluster=="pam"|(length(K)>1)){
		xdiss = 1 - cor(t(x),method=corr)
	}

	#Determine K if a range was given
	if(length(K)>1 ){
		pamklist=pamk(xdiss,krange=K,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)
		K = pamklist$nc
		genepamobj = pamklist$pamobject
	}else{
		genepamobj = NULL
	}


	#Do the gene clustering
	if(cluster=="pam"){
		if (!is.null(genepamobj)){
			geneclusters = genepamobj$clustering
		}else{
			geneclusters = pam(xdiss, k=K, diss=TRUE, keep.diss=FALSE, keep.data=FALSE, cluster.only=TRUE)
		}
	}else if(cluster=="kmeans"){
		#This loop ensures there are no empty clusters.
		done = FALSE
		while (!done){
			geneclusters = ckmeans(x,centers=K,iter.max=1000,distance=corr)
			done = !geneclusters$empty
		}
		geneclusters = geneclusters$cluster
	}else if(cluster=="flexclust"){


		distCor1 = function (x, centers) #change
		{
    		z <- matrix(0, nrow(x), ncol = nrow(centers))
    		for (k in 1:nrow(centers)) {
        		z[, k] <- 1 - cor(t(x), centers[k, ], method=corr)
    		}
    		z
		}
		#This loop ensures there are no empty clusters.
		done = FALSE
		while (!done){
			geneclusters = kcca(x, k=K, family=kccaFamily(dist=distCor1,cent="centMean"), simple=TRUE) #change
			#geneclusters = ckmeans(x,centers=K,iter.max=1000,distance=corr)
			done = (dim(geneclusters@clusinfo)[1]==K)#change
		}
		geneclusters = geneclusters@cluster
	}else if (cluster == "random"){
		geneclusters = sample(c(1:K,sample(1:K,m-K,replace = TRUE)),m,replace=FALSE)
	}else{
		cat("Unknown gene clustering method:",cluster,"\n")
	}

	return(geneclusters)

}


XPN_MLE_C = function(xbar,sigma2,nj){


	n = sum(nj)
	I = dim(xbar)[1]
	J = dim(xbar)[2]

	A = zeros(1,J)
	b = ones(I,1)
	c = zeros(I,1)
	s2 = ones(I,1)

	cout = .C("XPN_MLE_C_C",as.double(xbar),A=as.double(A),b=as.double(b),c=as.double(c),s2=as.double(s2),as.double(sigma2),as.integer(I),as.integer(J),as.integer(n),as.integer(nj))

	A = matrix(cout$A,1,J)
	b = matrix(cout$b,I,1)
	c = matrix(cout$c,I,1)
	s2 = matrix(cout$s2,I,1)

	return(list(A=A,b=b,c=c,s2=s2))
}

XPN_MLE = function(xbar,sigma2,nj){
	#This function is never used.  It serves to illustrate the method of maximum likelihood estimation used by XPN.  It has been replaced by XPN_MLE_C, which calls a more efficient C function.  The inputs and outputs of these functions are identical.

	n = sum(nj)
	I = dim(xbar)[1]
	J = dim(xbar)[2]

	A = zeros(1,J)
	b = ones(I,1)
	c = zeros(I,1)
	s2 = ones(I,1)

	old = c(A, b, c, s2)*0
	iter = 0
	while(sum((c(A, b, c, s2)-old)^2)>(1e-16)*max(old)){

		print(iter<-iter+1)
		print(sum((c(A, b, c, s2)-old)^2))
		old = c(A, b, c, s2)

		c = matrix(rowSums((xbar-repmat(b,1,J)*repmat(A,I,1))*repmat(nj,I,1))/n,I,1)

		if(sum(b)<0){
			b = -1*b;
		}

		A = matrix(colSums(repmat(b,1,J)*(xbar-repmat(c,1,J))/repmat(s2,1,J)/sum(b^2/s2)),1,J)

		A = A - mean(A)
		A = A * sqrt(J/sum(A^2))

		b = matrix(rowSums(repmat(A,I,1)*(xbar-repmat(c,1,J))*repmat(nj,I,1))/sum(A^2 * nj),I,1)


		s2 = matrix(rowSums(((xbar-repmat(c,1,J)-repmat(A,I,1)*repmat(b,1,J))^2 + sigma2)*repmat(nj,I,1))/n,I,1)



		s2[s2==0] = 2.2251e-308



	}

	return(list(A=A,b=b,c=c,s2=s2))
}



detect_L = function(x,y,L,corr="pearson"){

	#Determine L if a range has been given
	if(length(L)>1 & length(L) != (dim(x)[2]+dim(y)[2])){
		xdiss = 1 - cor(as.matrix(x), method = corr)
		ydiss = 1 - cor(as.matrix(y), method = corr)
		Lx=pamk(xdiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
		Ly = pamk(ydiss,krange=L,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)$nc
		L = max(Lx,Ly)
		cat(as.character(L), "assay clusters detected...\n")
	}
	return(L)
}

detect_K = function(x,K,corr="pearson"){

	#Determine K if a range was given
	if(length(K)>1 & length(K) != dim(x)[1]){
		xdiss = 1 - cor(t(x),method=corr)
		pamklist=pamk(xdiss,krange=K,diss=TRUE, keep.diss=FALSE, keep.data=FALSE)
		K = pamklist$nc
		genepamobj = pamklist$pamobject
	}else{
		genepamobj = NULL
	}
	return(K)
}

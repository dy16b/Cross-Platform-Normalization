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

sepelimdwd = function(Xp,Xn,penalty, useSparse=TRUE){
	#Xp and Xn are matrices.  penalty is a scalar.  This is an adaptation of the Matlab function of the same name, written by J. S. Marron, available at https://genome.unc.edu/pubsup/dwd/.

	flag <- 0
	#Dimensions of the data
	dp <- dim(Xp)[1]
	np <- dim(Xp)[2]
	dn <- dim(Xn)[1]
	nn <- dim(Xn)[2]
	if (dn != dp) {stop('The dimensions are incomapatible.')}
	d <- dp

	#Dimension reduction in HDLSS setting.
	XpnY <- as.matrix(cbind(Xp,-1*Xn))
	XpnY11 <- XpnY[1,1]
	n <- np + nn
	if(d>n){
		qrfact = qr(XpnY)
		Q = qr.Q(qrfact)
		R = qr.R(qrfact)
		RpnY = R
		dnew = n
	}else{
		RpnY = XpnY
		dnew = d
	}

	y = ones(np + nn,1)
	y[(np+1):(np+nn),1] = -1
	ym = y[2:n,1]

	# nv is the number of variables (eliminating beta)
	# nc is the number of constraints
	nv = 1 + dnew + 4*n
	nc = 2*n
	#Set up the block structure, constraint matrix, rhs, and cost vector

	blk = list()
	blk$type = character()
	blk$size = list()
	blk$type[1] = 's'
	blk$size[[1]] = cbind(dnew+1,3*ones(1,n))
	blk$type[2] = 'l'
	blk$size[[2]] = n


	Avec = list()
	A = zeros(nc,nv-n)
	col1 = RpnY[,1]
	A[1:(n-1),2:(dnew+1)] = t(RpnY[,2:n] - col1%*%t(ym))
	A[1:(n-1),seq(dnew+5,dnew+1+3*n,3)]	= -1*speye(n-1)
	A[1:(n-1),seq(dnew+6,dnew+2+3*n,3)]	= speye(n-1)
	A[1:(n-1),dnew+2] = ym
	A[1:(n-1),dnew+3] = -1*ym
	A[n,1] = 1
	A[(n+1):(n+n),seq(dnew+4,dnew+3+3*n,3)] = speye(n)


	Avec[[1]] = t(A)
	Avec[[2]] = (rbind(cbind(-1*ym,speye(n-1)),zeros(1+n,n)))
	b = rbind(zeros(n-1,1),ones(1+n,1))


	C = list()
	c = zeros(nv-n,1)
	c[seq(dnew+2,dnew+1+3*n,3),1] = ones(n,1)
	c[seq(dnew+3,dnew+2+3*n,3),1] = ones(n,1)


	C[[1]] = c
	C[[2]] = penalty*ones(n,1)


	#Solve the SOCP problem

##################CLSOCP############################
	CL_K <- unlist(blk$size)
	CL_qlen <-sum(blk$size[[1]])
	CL_llen <-sum(blk$size[[2]])
	CL_type <- c(rep('q',length(blk$size[[1]])),rep('l',length(blk$size[[2]])))
	CL_A <- cbind(t(Avec[[1]]),Avec[[2]])
	CL_c <- rbind(C[[1]],C[[2]])


	soln <- CLSOCP::socp(CL_A,b,CL_c,CL_K,CL_type,gamma_fac=.3,sigma0 = .1,use_sparse=useSparse)


	X1 <- soln$x[1:CL_qlen]
	X2 <- soln$x[(CL_qlen+1):(CL_qlen+CL_llen)]
	lambda <- soln$y
#####################################################


	# Compute the normal vector w and constant term beta.

	barw = X1[2:(dnew+1)]
	if (d>n){
		w = Q %*% barw
	}else{
		w = barw
	}
	beta = X1[dnew + 2] - X1[dnew + 3] - X2[1] - t(col1)%*%barw
	normw = norm(w)
	if (normw < 1 - 1e-3){
		print(normw)
	}
	normwm1 = 0
	if (normw > 1 - 1e-3){
		w = w/normw
		normwm1 = norm(w) - 1
		beta = beta/normw
	}


	# Compute the minimum of the supposedly positive
	# and the maximum of the supposedly negative residuals.
	# Refine the primal solution and print its objective value.

	residp = t(Xp) %*% w + beta[1] #optimization
	residn = t(Xn) %*% w + beta[1] #optimization
	minresidp = min(residp)
	maxresidn = max(residn)
	res = t(XpnY) %*% w + beta[1] * y
	rsc = 1/sqrt(penalty)
	xi = -1* res + rsc[1]
	xi[xi<0] <- 0
	totalviolation = sum(xi)
	minresidpmod = min(residp + xi[1:np])
	maxresidnmod = max(residn - xi[(np+1):n])
	minxi = min(xi)
	maxxi = max(xi)
	resn = res + xi
	rresn = 1 / resn
	primalobj = penalty * sum(xi) + sum(rresn)
#	print(primalobj)


	#Compute the dual solution alp and print its objective value.
	alp = zeros(n,1)
	lambda1 = lambda[1:(n-1)]
	alp[1] = -1*t(ym)%*%lambda1
	alp[2:n] = lambda1
	alp = alp * (as.numeric(alp>0))
	sump = sum(alp[1:np])
	sumn = sum(alp[(np+1):n])
	sum2 = (sump + sumn)/2
	alp[1:np] = (sum2/sump)*alp[1:np]
	alp[(np+1):n] = (sum2/sumn)*alp[(np+1):n]
	maxalp = max(alp)
	if (maxalp > penalty | maxxi > 1e-3){
		alp = (penalty[1]/maxalp)*alp
	}
	minalp = min(alp)
	p = RpnY%*%alp
	eta = -1*norm(p)
	gamma = 2*sqrt(alp)
	dualobj = eta + sum(gamma)
#	print(dualobj)



	#dualgap is a measure of the accuracy of the solution
	dualgap = primalobj - dualobj
#	print(dualgap)

	if (dualgap > 1e-4){
		flag = -1
	}


	return(list(w=w,beta=beta,residp=residp,residn=residn,alp=alp,totalviolation=totalviolation,dualgap=dualgap,flag=flag))

}


DWD1SM = function(trainp,trainn,threshfact = 100, useSparse = TRUE){

	np = dim(trainp)[2]
	nn = dim(trainn)[2]
#	vpwdist2 = numeric(np*nn)
	vpwdist2x <- fields::rdist(t(trainp),t(trainn))
#	for (ip in 1:np){
#		vpwdist2[((ip-1)*nn+1):(ip*nn)] <- colSums((trainp[,ip] - trainn)^2) #optimization
#	}
#	medianpwdist2 = median(vpwdist2)
	medianpwdist2 = median(vpwdist2x)^2

	penalty = threshfact / medianpwdist2
	sepelimout = sepelimdwd(trainp,trainn,penalty,useSparse=useSparse)
	w = sepelimout$w
	flag = sepelimout$flag
	if (flag == -1){
		cat("Inaccurate solution!\n")
	}
	if (flag == -2){
		cat("Infeasible or unbounded optimization problem!\n")
	}
	dirvec = w/norm(w)
	return(dirvec)
}

#'@param X gene expression matrix
#'@param Y gene expression matrix
#'@export
dwd = function(platform1.data, platform2.data, platform1.train = NULL, platform2.train=NULL, p1.names=0, p2.names=0,p1.train.names=0, p2.train.names=0, skip.match=FALSE, use.sparse = TRUE) {

	#Match names
	if(is.null(platform1.train) & is.null(platform2.train)){
		input = processplatforms(list(x=platform1.data,y=platform2.data), namesvec = c(p1.names, p2.names), skip.match=skip.match)



		#Get normal vector for DWD adjustment
		dirvec = DWD1SM(input[[1]],input[[2]],useSparse=use.sparse)
		#Project the data
		vprojp = t(input[[1]]) %*% dirvec
		vprojn = t(input[[2]]) %*% dirvec
		meanprojp = mean(vprojp)
		meanprojn = mean(vprojn)
		output = list()
		p1.adjust <- -1 * meanprojp * dirvec
		p2.adjust <- -1 * meanprojn * dirvec
		names(p1.adjust) <- rownames(input[[1]])
		names(p2.adjust) <- rownames(input[[2]])
		output$x = input[[1]] + p1.adjust
		output$y = input[[2]] + p2.adjust
		output$p1.adjust <- p1.adjust
		output$p2.adjust <- p2.adjust

	}else{
		input = processplatforms(list(x=platform1.data,y=platform2.data,platform1.train,platform2.train), namesvec = c(p1.names, p2.names, p1.train.names, p2.train.names), skip.match=skip.match)


		#Get normal vector for DWD adjustment
		dirvec = DWD1SM(input[[3]],input[[4]],useSparse=use.sparse)

		#Project the data
		vprojp = t(input[[3]]) %*% dirvec
		vprojn = t(input[[4]]) %*% dirvec

		meanprojp = mean(vprojp)
		meanprojn = mean(vprojn)

		output = list()
		p1.adjust <- -1 * meanprojp * dirvec
		p2.adjust <- -1 * meanprojn * dirvec
		names(p1.adjust) <- rownames(input[[1]])
		names(p2.adjust) <- rownames(input[[2]])
		output$x = input[[1]] + p1.adjust
		output$y = input[[2]] + p2.adjust
		output$p1.adjust <- p1.adjust
		output$p2.adjust <- p2.adjust
	}

	return(output)

}

norm = function(aMatrix){
#Returns the largest singular value of aMatrix
	o = svd(aMatrix,nu=0,nv=0)
	return(o$d[1])

}




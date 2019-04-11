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

cross.platform.normalize=function(method,...){
	if(is.character(method)){
		method = get(method)
	}
	method(...)
}


processplatforms = function(datalist, namesvec=NULL, skip.match=FALSE){
	#Convert data from various formats to the proper format for use 
	#with all the crossnorm normalization functions
	
	for(i in 1:length(datalist)){
		if(is.matrix(datalist[[i]])){
			datalist[[i]] <- as.data.frame(datalist[[i]])
		}
	}
	
	if (is.null(namesvec)){
		namesvec <- numeric(length(datalist))
		for (i in 1:length(datalist)){
			namesvec[i] <- 0
		}
	}
	
	#Put the row names in their places
	for (i in 1:length(namesvec)){
		if(namesvec[i] != 0){
			rownames(datalist[[i]]) = datalist[[i]][,namesvec[i]]
			datalist[[i]] = datalist[[i]][,-1*namesvec[i],drop=FALSE]
		}	
	}

	if(!skip.match){
		#Create the common genes list
		commongenes <- rownames(datalist[[1]])
		for (i in 2:length(datalist)){
			commongenes <- intersect(commongenes,rownames(datalist[[i]]))
		}
	
	
		#Put it all together
		for (i in 1:length(datalist)){
			datalist[[i]] <- datalist[[i]][commongenes,,drop=FALSE]
		}
	}
	return(datalist)
}


assert = function(aBool){
	if (!aBool){
		stop("Assert failed!")
	}	
}
	
repmat = function(mat,m,n){
	out = mat
	if(m>1){
	for (i in 2:m){
		out = rbind(out,mat)
	}
	}
	if(n>1){
	mat = out
	for (j in 2:n){
		out = cbind(out,mat)
	}
	}
	return(out)
}
ones = function(m,n){
	return(matrix(1,nrow=m,ncol=n))
}

zeros = function(m,n){
	return(matrix(0,nrow=m,ncol=n))
}

speye = function(m){
#This function is not actually sparse yet.  Needs to be fixed.

	return(diag(m))
	
}



rowMedians = function(aMatrix, na.rm = FALSE){
#Also works on data.frames
	m = dim(aMatrix)[1]
	ret = numeric(m)
	for (i in 1:m){
		ret[i] = mean(median(aMatrix[i,],na.rm=na.rm),na.rm=TRUE)
	}
	return(ret)
}

randn = function(m,n){
	return(matrix(rnorm(m*n),m,n))	
}

colRanks = function(aMatrix){
    return(apply(aMatrix,2,rank))
}

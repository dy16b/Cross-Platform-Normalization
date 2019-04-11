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


######################################################################
#This file contains code adapted from the WebArrayDB program,
#available from http://www.webarraydb.org/webarray/index.html  Please
#cite appropriately when using these functions.  Type
#citation("CONOR") at the R prompt for details.

## Implementation of the "Gene Quantile" cross platform
## normalization algorithm, as available from webarraydb.
## The function gq() is a wrapper for normalizeGQ, which
## was provided by the authors of webarraydb.
#'
#'fit mixed effect regression model
#'@param X gene expression matrix
#'@param Y gene expression matrix
#'@export
gq = function(platform1.data, platform2.data, p1.names=0, p2.names=0, skip.match=FALSE){
	#This function is basically a wrapper for normalizeGQ

	#Match names
	input = processplatforms(list(x=platform1.data,y=platform2.data),namesvec = c(p1.names, p2.names), skip.match=skip.match)

	#Prepare for normalizeGQ
	combined = cbind(input$x,input$y)
	pf = c(seq(1,1,length.out=dim(input$x)[2]),seq(2,2,length.out=dim(input$y)[2]))

	#Call normalizeGQ
	ngq = normalizeGQ(combined,pf)

	#Split the results and return
	out=split(seq(pf),pf)
	out[[1]] = ngq[,out[[1]]]
	out[[2]] = ngq[,out[[2]]]
	names(out) <- c("x","y")
	return(out)
}


normalizeGQ <- function(M, pf, ...) {
	#This function was provided by Xiao-Qin Xia, one of the authors of
	#webarraydb.
    # modified MRS
    # M is the data matrix
    # pf is the vector to specify the platform for each column of M.
    idx <- split(seq(pf), pf)
    if (length(pf)<=1) return(M)
    imax <- which.max(sapply(idx, length)) # for reference
    ref_med <- apply(M[, idx[[imax]]], 1, function(x) median(x, na.rm=TRUE))
    ref_med_srt <- sort(ref_med)
    idx[imax] <- NULL
    lapply(idx, function(i) {
         MTMP <- sapply(i, function(x) ref_med_srt[rank(M[,x])]);
         M[,i] <<- MTMP - apply(MTMP, 1, median) + ref_med
         } )
    invisible(M)
}


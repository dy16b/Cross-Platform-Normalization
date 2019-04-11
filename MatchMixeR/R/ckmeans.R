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

ckmeans = function (x, centers, iter.max = 10, nstart = 1, distance = "pearson")
{


	if (distance == "pearson"){
		distance = 1
	}else if(distance == "euclidean"){
		distance = 2
	}else if(distance == "manhattan"){
		distance = 3
	}else if(distance == "absolutepearson"){
		distance = 4
	}
	nmeth = "lloyd"

	#Normalization consistent with flexclust
	x = x/sqrt(rowSums(x^2))


    do_one <- function(nmeth) {
      #' @useDynLib MatchMixeR ckmeans_c
            Z <- .C(ckmeans_c, as.double(x), as.integer(m),
                as.integer(ncol(x)), centers = as.double(centers),
                as.integer(k), c1 = integer(m), iter = as.integer(iter.max),
                nc = integer(k), wss = double(k), as.integer(distance))
            if (Z$iter > iter.max) {
            	warning("did not converge in ",
                iter.max, " iterations", call. = FALSE)
                Z$converged = FALSE
            }else{
                Z$converged = TRUE
            }
            if (any(Z$nc == 0)) {
            	warning("empty cluster: try a better set of initial centers",
                call. = FALSE)
                Z$empty <- TRUE
            }else{
                Z$empty <- FALSE
            }
			Z
    }
    x <- as.matrix(x)
    m <- nrow(x)
    if (missing(centers))
        stop("'centers' must be a number or a matrix")

    if (length(centers) == 1L) {
        k <- centers
        if (nstart == 1)
            centers <- x[sample.int(m, k), , drop = FALSE]
        if (nstart >= 2 || any(duplicated(centers))) {
            cn <- unique(x)
            mm <- nrow(cn)
            if (mm < k)
                stop("more cluster centers than distinct data points.")
            centers <- cn[sample.int(mm, k), , drop = FALSE]
        }
    }
    else {
        centers <- as.matrix(centers)
        if (any(duplicated(centers)))
            stop("initial centers are not distinct")
        cn <- NULL
        k <- nrow(centers)
        if (m < k)
            stop("more cluster centers than data points")
    }
    if (iter.max < 1)
        stop("'iter.max' must be positive")
    if (ncol(x) != ncol(centers))
        stop("must have same number of columns in 'x' and 'centers'")
    Z <- do_one(nmeth)
    if (nstart >= 2 && !is.null(cn)) {
        best <- sum(Z$wss)
        for (i in 2:nstart) {
            centers <- cn[sample.int(mm, k), , drop = FALSE]
            ZZ <- do_one(nmeth)
            if ((z <- sum(ZZ$wss)) < best) {
                Z <- ZZ
                best <- z
            }
        }
    }
    centers <- matrix(Z$centers, k)
    dimnames(centers) <- list(1L:k, dimnames(x)[[2L]])
    cluster <- Z$c1
    if (!is.null(rn <- rownames(x)))
        names(cluster) <- rn
    out <- list(cluster = cluster, centers = centers, withinss = Z$wss,
        size = Z$nc, empty = Z$empty, converged = Z$converged)
    class(out) <- "kmeans"
    out
}

/*
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
*/
#include <R.h>
#include <math.h>
//#include "modreg.h" /* for declarations for registration */

double pearson(double *x, double *cen,  int p, int i, int j, int n, int k);
double euclidean(double *x, double *cen, int p, int i, int j, int n, int k);
double manhattan(double *x, double *cen, int p, int i, int j, int n, int k);
double absolutepearson(double *x, double *cen, int p, int i, int j, int n, int k);

void ckmeans_c(double *x, int *pn, int *pp, double *cen, int *pk, int *cl, int *pmaxiter, int *nc, double *wss, int *pdistance )
{
    int n = *pn, k = *pk, p = *pp, maxiter = *pmaxiter, distance = *pdistance;
    int iter, i, j, c, it, inew = 0;
    double best, dd, tmp;
    Rboolean updated;
	
    for(i = 0; i < n; i++) cl[i] = -1;
    for(iter = 0; iter < maxiter; iter++) {
	updated = FALSE;
	for(i = 0; i < n; i++) {
	    /* find nearest centre for each point */
	    best = R_PosInf;
	    for(j = 0; j < k; j++) {
		dd = 0.0;
		
		switch (distance){
			case 1: dd = pearson(x,cen,p, i, j, n, k);
				break;
			case 2: dd = euclidean(x,cen,p, i, j, n, k);
				break;
			case 3: dd = manhattan(x,cen,p, i, j, n, k);
				break;
			case 4: dd = absolutepearson(x,cen,p, i, j, n, k);
				break;
			default: error("Unkown distance measure provided to ckmeans.");
		}
		if(dd < best) {
		    best = dd;
		    inew = j+1;
		}
	    }
	    if(cl[i] != inew) {
		updated = TRUE;
		cl[i] = inew;
	    }
	}
	if(!updated) break;
	/* update each centre */
	for(j = 0; j < k*p; j++) cen[j] = 0.0;
	for(j = 0; j < k; j++) nc[j] = 0;
	for(i = 0; i < n; i++) {
	    it = cl[i] - 1; nc[it]++;
	    for(c = 0; c < p; c++) cen[it+c*k] += x[i+c*n];
	}
	for(j = 0; j < k*p; j++) cen[j] /= nc[j % k];
    }

    *pmaxiter = iter + 1;
    for(j = 0; j < k; j++) wss[j] = 0.0;
    for(i = 0; i < n; i++) {
	it = cl[i] - 1;
	for(c = 0; c < p; c++) {
	    tmp = x[i+n*c] - cen[it+k*c];
	    wss[it] += tmp * tmp;
	}
    }
}

double pearson(double *x, double *cen,  int p, int i, int j, int n, int k){
	double top = 0;
	double bottom1 = 0;
	double bottom2 = 0;
	int c;
	for (c = 0; c < p; c++){
		top += x[i+n*c] * cen[j+k*c];
		bottom1 += x[i+n*c]*x[i+n*c];
		bottom2 += cen[j+k*c]*cen[j+k*c];
	}
	return(1 - top/(sqrt(bottom1)*sqrt(bottom2)));
}






double euclidean(double *x, double *cen, int p, int i, int j, int n, int k){
	double tmp, dd = 0;
	int c;
	for(c = 0; c < p; c++) {
	    tmp = x[i+n*c] - cen[j+k*c];
	    dd += tmp * tmp;
	}
	return(dd);
}


double manhattan(double *x, double *cen, int p, int i, int j, int n, int k){
	double dd = 0;
	int c;
	for(c = 0; c < p; c++) {

	    dd += abs(x[i+n*c] - cen[j+k*c]);
	}
	return(dd);
}


double absolutepearson(double *x, double *cen, int p, int i, int j, int n, int k){
	
	return(abs(pearson(x,cen,p, i, j, n, k)));

}

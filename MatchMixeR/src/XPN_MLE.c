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
#include <float.h>
#include <stdio.h>


void XPN_MLE_C_C(double *xbar, double *A, double *b, double *c, double *s2, double *sigma2, int *pI, int *pJ, int *pn, int *nj){
	int I = *pI, J = *pJ, n = *pn;
	double *old;
	double maxold, tmp, sumb, sumA, sumA2, meanA, Anorm, sumnjA2, sumb2overs2, change = 1, radius = 0.0000000000000001;
	int changed = 1, oldlength = J+3*I;
	old = malloc((oldlength)*sizeof(double));
	
	//Copy starting values
	int i,j;
	maxold = 0;
	for(i=0;i<J;i++){
		old[i] = A[i];
		if(old[i] > maxold){
			maxold = old[i];
		}
	}
	for(i=0;i<I;i++){
		old[i + J] = b[i];
		if(old[i+J] > maxold){
			maxold = old[i+J];
		}
	}
	for(i=0;i<I;i++){
		old[i + J+I] = c[i];
		if(old[i+J+I] > maxold){
			maxold = old[i+J+I];
		}
	}
	for(i=0;i<I;i++){
		old[i + J+2*I] = s2[i];
		if(old[i+J+2*I] > maxold){
			maxold = old[i+J+2*I];
		}
	}

	

	//Until we stabilize
	int iter = 0;
	while(changed && (iter < 100)){
		iter++;
	//	Rprintf("%d ",iter);
		fflush(stdout);


	
		//update c
		for (i=0;i<I;i++){
			c[i] = 0;
			for (j=0;j<J;j++){
				c[i] += (xbar[i + j*I] - b[i] * A[j]) * nj[j];
			}
			c[i] = c[i]/n;
		//	Rprintf("c[%d]: %f\n",i,sumb);
		}
		
		//Correct b if needed
		sumb =0;
		for (i=0;i<I;i++){
			sumb += b[i];
			
		}
		
		if (sumb < 0) {
			for(i=1;i<I;i++){
				b[i] = -1*b[i];
			}
			sumb = -1*sumb;
		}

		
		//Update A
		sumb2overs2 = 0;
		for(i=0;i<I;i++){
			sumb2overs2 += b[i]*b[i]/s2[i];
		}
		sumA = 0;
		for(j=0;j<J;j++){
			A[j] = 0;
			for(i=0;i<I;i++){
				A[j] += b[i]*(xbar[i+j*I]-c[i])/s2[i];
			}
			A[j] = A[j]/sumb2overs2;
			sumA += A[j];
		//	Rprintf("A[%d]: %f\n",j,A[j]);
		}
		meanA = sumA/J;
		for(j=0;j<J;j++){
			A[j] = A[j] - meanA;
		}
		sumA2 = 0;
		for(j=0;j<J;j++){
			sumA2 += A[j]*A[j];
		}
		Anorm = sqrt(J/sumA2);
		for(j=0;j<J;j++){
			A[j] = A[j]*Anorm;
		}

		
		//Update b
		sumnjA2 = 0;
		for(j=0;j<J;j++){
			sumnjA2 += nj[j]*A[j]*A[j];
		}
		for (i=0;i<I;i++){
			b[i] = 0;
			for (j=0;j<J;j++){
				b[i] += A[j] * (xbar[i + I*j] - c[i]) * nj[j];
			}
			b[i] = b[i]/sumnjA2;
		//	Rprintf("b[%d]: %f\n",i,b[i]);
		}
		
		//update s2
		for (i=0;i<I;i++){
			s2[i] = 0;
			for (j=0;j<J;j++){
		//		Rprintf("Going into s2 calculation, xbar[%d,%d] = %f, c[%d] = %f, A[%d] = %f, b[%d] = %f n = %d\n", i,j,xbar[i + I*j] ,i, c[i],j,A[j],i,b[i],n);
				tmp = (xbar[i + I*j] - c[i] - A[j]*b[i]);
				s2[i] += (tmp*tmp + sigma2[i+j*I])*nj[j];
			}
			s2[i] = s2[i]/n;
			if (s2[i] == 0){
		//		Rprintf("OMG!!!\n");
				s2[i] = DBL_MIN;
			}
		//	Rprintf("s2[%d]: %f\n",i,s2[i]);
		}

		//Calculate change
		change = 0;
		for(i=0;i<J;i++){
			tmp = old[i] - A[i];
		//	Rprintf("%f ",tmp);
			change += tmp*tmp;
		}
		for(i=0;i<I;i++){
			tmp = old[i + J] - b[i];
		//	Rprintf("%f ",tmp);
			change += tmp * tmp;
		}
		for(i=0;i<I;i++){
			tmp = old[i + J + I] - c[i];
		//	Rprintf("%f ",tmp);
			change += tmp*tmp;
		}
		for(i=0;i<I;i++){
			tmp = old[i + J + 2*I] - s2[i];
		//	Rprintf("%f ",tmp);
			change += tmp*tmp;
		}
		
		changed = (change > radius*maxold);
		//printf("%f ",change);
		
		//Copy previous values
		maxold = 0;
		for(i=0;i<J;i++){
			old[i] = A[i];
			if(old[i] > maxold){
				maxold = old[i];
			}
		}
		for(i=0;i<I;i++){
			old[i + J] = b[i];
			if(old[i+J] > maxold){
				maxold = old[i+J];
			}
		}
		for(i=0;i<I;i++){
			old[i + J+I] = c[i];
			if(old[i+J+I] > maxold){
				maxold = old[i+J+I] ;
			}
		}
		for(i=0;i<I;i++){
			old[i + J+2*I] = s2[i];
			if(old[i+J+2*I] > maxold){
				maxold = old[i+J+2*I];
			}
		}

	}	

	free(old);

}

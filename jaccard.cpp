//
//  Jaccard.cpp
//
//
//
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Jac_dist(IntegerVector x, IntegerVector y, int kx, int ky)
{
    int i,j,k,l=x.size();
    float v1,v2,v3,r1;
    NumericMatrix jac(kx,ky);
    for (i=0;i<kx;i++){
        for (j=0;j<ky;j++){
            v1=v2=v3=0.0;
            for (k=0;k<l;k++){
                if (x[k]==i && y[k]!=j)
                    v1+=1.0; //in cluster i but not in cluster j
                if (x[k]!=i && y[k]==j)
                    v2+=1.0; //in cluster j but not in cluster i
                if (x[k]==i && y[k]==j)
                    v3+=1.0; //in both clusters
            }
            if (v3+v1+v2==0.0) r1=1.0; else r1=(v1+v2)/(v3+v1+v2);
            jac(i,j)=r1;
        }
    }
    return jac;
}









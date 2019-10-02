/*
 * lakemaggiore.cpp
 *
 *  Created on: 7/feb/2014
 *      Author: MatteoG
 */

#include "lakemaggiore.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <math.h>

using namespace std;
#define K pow(10,9)

lakemaggiore::lakemaggiore() {
    // TODO Auto-generated constructor stub
    theta = utils ::loadVector("../data/maggiore_theta.txt", 7); //vector containing natural storage discharge relation parameters
    theta1 = utils ::loadVector("../data/maggiore_theta1.txt", 2); // vector containing maximum storage discharge relation parameters where s <= -103311800
    theta2 = utils ::loadVector("../data/maggiore_theta2.txt", 2); // vector containing minimum storage discharge relation parameters where s <= -103311800
    smax = utils ::loadVector("../data/s_max.txt", 365);
}

lakemaggiore::~lakemaggiore() {
    // TODO Auto-generated destructor stub
}

double lakemaggiore::storageToLevel(double s){
    double h;
    if(lsv_rel.size()>0){
        h = utils::interp_lin(lsv_rel[1],lsv_rel[0],s);
    }else{
        h = s/A;
    }
    // double h = 4.23836E-09 * s + 0.13051;
    return h;
}

double lakemaggiore::levelToStorage(double h){
    double s;
    if(lsv_rel.size()>0){
        s = utils::interp_lin(lsv_rel[0],lsv_rel[1],h);
    }else{
        s = h*A;
    }
    //double s = 2.35639E+08 * h - 2.94332E+07;
    return s;
}

double lakemaggiore::levelToSurface(double h){

	double S ;
	if (h > 0){
		S = (0.13000 * pow(h,2) + 4.65000 * h + 207.51000)*1000000 ;
	} else {
		S = A;
	}
    return S;
}


double lakemaggiore::s_max(int cday){

  double s = 0.0;

  s = smax[cday-1];
	return s;
}


double lakemaggiore::min_release(double s, int cday){
    double q = 0.0;
    double DMV = minEnvFlow[cday-1];
    if(s <= -105420000){
        q = 0.0;
    }else if (s <= -103311800) {
        q =  theta2[0]*s*K +  theta2[1];
    }else if (s <= s_max(cday)) {
        q = DMV;
    }else{
        //q = 0.0567 * pow(h,4) + 0.6130 * pow(h,3) + 58.0684*pow(h,2) + 232.4234*h + 302.5700;
        q = theta[0]*(s/K)*(s/K)*(s/K)*(s/K)*(s/K)*(s/K) +theta[1]*(s/K)*(s/K)*(s/K)*(s/K)*(s/K) + theta[2]*(s/K)*(s/K)*(s/K)*(s/K) + theta[3]*(s/K)*(s/K)*(s/K)+ theta[4]*(s/K)*(s/K)+ theta[5]*(s/K) + theta[6];
    }
    return q;
}

double lakemaggiore::max_release(double s, int cday){

    double q = 0.0;
    if(s <= -105420000){
        q = 0.0;
    }else if (s <= -103311800) {
        q = theta1[0]*s/K + theta1[1];
    }else{
        //q = 0.0567 * pow(h,4) + 0.6130 * pow(h,3) + 58.0684*pow(h,2) + 232.4234*h + 302.5700;
        q = theta[0]*(s/K)*(s/K)*(s/K)*(s/K)*(s/K)*(s/K) +theta[1]*(s/K)*(s/K)*(s/K)*(s/K)*(s/K) + theta[2]*(s/K)*(s/K)*(s/K)*(s/K) + theta[3]*(s/K)*(s/K)*(s/K)+ theta[4]*(s/K)*(s/K)+ theta[5]*(s/K) + theta[6];
    }
    return q;
}

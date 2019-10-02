/*
 * pwLinear.cpp
 *
 *  Created on: 31/oct/2013
 *      Author: MatteoG
 */

#include "pwLinear.h"

using namespace std;


pwLinear::pwLinear()
{
}

pwLinear::~pwLinear(){
}


pwLinear::pwLinear(unsigned int pM, unsigned int pK, unsigned int pN){
    M = pM;
    K = pK;
    N = pN;
}


void pwLinear::setParameters(double* pTheta){

    // piecewise linear parameters
    unsigned int count = 0;
    for(unsigned int i=0; i<N-1; i++){
        param.a.push_back( pTheta[count] );
        count = count+1;
        param.b.push_back( pTheta[count] );
        count = count+1;
        param.c.push_back( pTheta[count] );
        count = count+1;
    }
    param.a.push_back( pTheta[count] );
    count = count+1;
    param.b.push_back( pTheta[count] );
    param.c.push_back( -999 );

}


void pwLinear::clearParameters(){
    param.a.clear();
    param.b.clear();
    param.c.clear();
}

vector<double> pwLinear::get_output(vector<double> input){

    double y = 0.0;
    for(unsigned int i=0; i<N-1; i++){
        if( input[0] <= param.c[i] ){
            y = param.a[i] + param.b[i]*input[0] ;
        }
    }
    if( input[0] > param.c[N-2] ){
        y = param.a[N] + param.b[N]*input[0] ;
    }
    vector<double> yy;
    yy.push_back( y );

    return yy;
}

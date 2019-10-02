/*
 * maggiore_panperduto.cpp
 */


#include "maggiore_panperduto.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <algorithm>

using namespace std;

maggiore_panperduto::maggiore_panperduto(){

	a_ic = utils ::loadVector("../data/Maggiore_Panperduto/a_ic.txt", 365);
	a_rec = utils ::loadVector("../data/Maggiore_Panperduto/a_rec.txt", 365);
	a_vc = utils ::loadVector("../data/Maggiore_Panperduto/a_vc.txt", 365);
 	q_turbigo = utils ::loadVector("../data/Maggiore_Panperduto/q_turbigo.txt", 365);
 	max_rid_ic = utils ::loadVector("../data/Maggiore_Panperduto/max_rid_ic.txt", 365);
 	max_rid_rec = utils ::loadVector("../data/Maggiore_Panperduto/max_rid_rec.txt", 365);
 	max_rid_vc = utils ::loadVector("../data/Maggiore_Panperduto/max_rid_vc.txt", 365);
 	q_min_rec = utils ::loadVector("../data/Maggiore_Panperduto/q_min_rec.txt", 365);
 	q_min_vc = utils ::loadVector("../data/Maggiore_Panperduto/q_min_vc.txt", 365);
 	q_ope_rec = utils ::loadVector("../data/Maggiore_Panperduto/q_ope_rec.txt", 365);
 	q_ope_vc = utils ::loadVector("../data/Maggiore_Panperduto/q_ope_vc.txt", 365);
 	q_thre_ic = utils ::loadVector("../data/Maggiore_Panperduto/q_thre_ic.txt", 365);
 	q_thre_rec = utils ::loadVector("../data/Maggiore_Panperduto/q_thre_rec.txt", 365);
 	DMV = utils ::loadVector("../data/MEF_maggioreTV.txt", 365);

    // TODO Auto-generated constructor stub
}

maggiore_panperduto::~maggiore_panperduto(){
    // TODO Auto-generated destructor stub
}

vector<vector<double> > maggiore_panperduto::Distribution_rule( vector<double> q, vector<double> d, vector<double> dREC, vector<double> dVC, vector<double> dIC, vector<double> doy){

double hlp, deficit;
//vector <double> hlp, deficit;
vector<double> q_c;
vector<vector<double > > q_flow, q_alloc;

for( int i=0; i<q.size(); i++){

//cout << "i=" << i << endl;

// Compute flow avaiable
hlp = q[i]- DMV[doy[i]-1];
//hlp.push_back( q[i]- DMV[doy[i]-1] );
//Compute current deficit
deficit = d[doy[i]-1] -hlp;
//deficit.push_back( d[doy[i]-1] - hlp[i] );


//1Condition for normal case
if (deficit <= 0)
{ q_c = Normal_Condition( dREC[doy[i]-1], dVC[doy[i]-1], dIC[doy[i]-1]);
	//cout << "NC" << endl;
	//id[i] = (double) 1;
}
//2Condition for Hydropower reduction
else if (( deficit > 0 ) && ( deficit <=  (max_rid_ic[doy[i]-1] + max_rid_rec[doy[i]-1])))
{ q_c = Hydro_reduction ( hlp, d[doy[i]-1], dREC[doy[i]-1], dVC[doy[i]-1], dIC[doy[i]-1], max_rid_rec[doy[i]-1], max_rid_ic[doy[i]-1], q_thre_rec[doy[i]-1], q_thre_ic[doy[i]-1] ,q_ope_rec[doy[i]-1], deficit );
	//cout << "HR" << endl;
}
//3Condition for irr reduction
else if ( (  deficit >  ( max_rid_ic[doy[i]-1]+max_rid_rec[doy[i]-1] ) ) && ((deficit<= max_rid_rec[doy[i]-1]+max_rid_ic[doy[i]-1]+max_rid_vc[doy[i]-1])) )
{ q_c = Irrig_reduction( hlp, d[doy[i]-1], dREC[doy[i]-1], dVC[doy[i]-1], dIC[doy[i]-1], max_rid_rec[doy[i]-1], max_rid_ic[doy[i]-1]);
	//cout << "IR" << endl;
}
//4Conflict or extreme case
else if (  deficit > max_rid_ic[doy[i]-1]+max_rid_rec[doy[i]-1]+max_rid_vc[doy[i]-1])
{ q_c = Conflict ( hlp, d[doy[i]-1], dREC[doy[i]-1], dVC[doy[i]-1], dIC[doy[i]-1], q_turbigo[doy[i]-1], a_rec[doy[i]-1], a_ic[doy[i]-1], a_vc[doy[i]-1], q_min_vc[doy[i]-1], q_ope_vc[doy[i]-1], q_min_rec[doy[i]-1], q_ope_rec[doy[i]-1], max_rid_rec[doy[i]-1] );
	//cout << "CO" << endl;
}

else { q_c.push_back( (double) -999);
			q_c.push_back( (double) -999);
			q_c.push_back( (double) -999);
			q_c.push_back( (double) -999);
			//cout << "W" << endl;
}

q_flow.push_back( q_c );
}

//cout << "i=" << q_flow[364][2]<< endl;
//cout << "deficit[4751]=" << deficit[4751] << endl;
//cout << "hlp[4751]=" << hlp[4751] << endl;
return q_flow;
}

vector<double> maggiore_panperduto::Normal_Condition( double re_dem, double vc_dem, double ic_dem){

vector<double> q;
q.push_back( re_dem );
q.push_back( vc_dem );
q.push_back( ic_dem );
q.push_back( (double) 1 );

return q;
}

vector<double> maggiore_panperduto::Hydro_reduction ( double q, double tot_dem, double re_dem, double vc_dem, double ic_dem, double maxRe, double maxIc, double qThreR, double qThreI, double opeR, double deficit){
vector<double> qout;
if (deficit <= qThreR ) {
  // solo il CRE riduce i prelievi, il Villoresi e il CI prendono la loro domanda
  qout.push_back( re_dem-deficit );
  qout.push_back( vc_dem );
  qout.push_back( ic_dem );
	qout.push_back( (double) 21 );
}

else if ( ( deficit > qThreR ) && ( deficit <= ( qThreR+ qThreI ) ) ){
  // caso 1b: anche il CI riduce i prelievi: il CRE riduce fino ad un massimo e il ci riduce della differenza
  qout.push_back( re_dem-qThreR );
  qout.push_back( vc_dem );
  qout.push_back( ic_dem -  (deficit-qThreR) );
	qout.push_back( (double) 22 );
}

else if ( ( deficit > ( qThreR + qThreI ) && (deficit <=  ( maxRe + qThreI  )))){
  // caso 1c: il CRE riduce ulteriormente i prelievi: il CI riduce i prelievi e il CRE riduce della differenza ma se qR< q operativa GLI ASSEGNAMO MIN ASS PER PESCI
		double qV, qI, qR, i;
		qI = ic_dem - qThreI;
		qR = re_dem - deficit + qThreI;
		i = 232;
		if ( qR < opeR ){
			qR = re_dem - maxRe;
			qI = ic_dem - deficit + maxRe;
			i = 231;
		}
		qV = vc_dem;
		qout.push_back( qR );
		qout.push_back( qV );
		qout.push_back( qI );
		qout.push_back( i );

}

else if( deficit > ( maxRe + qThreI ) ) {
//il CI riduce ulteriormente il prelievo: il cre riduce la portata fino al minimo e il ci della differenza
qout.push_back( re_dem - maxRe );
qout.push_back( vc_dem );
qout.push_back( ic_dem - ( deficit - maxRe ) );
qout.push_back( (double) 24 );

}

return qout;
}


vector<double> maggiore_panperduto::Irrig_reduction( double q, double tot_dem, double re_dem, double vc_dem, double ic_dem, double maxRe, double maxIc ){
vector<double> qout;
qout.push_back( re_dem - maxRe );
qout.push_back( q - ((re_dem - maxRe)+(ic_dem - maxIc)) );
qout.push_back( ic_dem - maxIc );
qout.push_back( (double) 3 );
return qout;
}

vector<double> maggiore_panperduto::Conflict (double q, double tot_dem, double re_dem, double vc_dem, double ic_dem, double qTur, double aR, double aI, double aV, double qminV, double opeV, double qminR, double opeR, double maxRe){
vector<double> qout;
double qV, qCRE, qIC, villo2elena, elena2indu ;

if (q < (double) 25) {
		qout.push_back( (double) 0);
		qout.push_back( (double) 0);
		qout.push_back( (double) 0);
		qout.push_back( (double) 5);
		//cout << "EX" << endl;
}

else {
// Il CV prende la sua % di portata meno la portata di Turbigo, ma se e' meno della portata minima trasferisce al regina elena
villo2elena = (double) 0 ;
if ( (q-qTur)*aV < qminV ) {
	villo2elena = (q-qTur)*aV;
	qV = (double) 0 ; }
else { qV = (q-qTur)*aV; }
if( ((q-qTur)*aV < opeV) && ( (q-qTur)*aV >= qminV ) ) {
	villo2elena = (q-qTur)*aV-qminV;
	qV = qminV;
}
else if ((q-qTur)*aV > vc_dem) {
	villo2elena = (q-qTur)*aV - vc_dem;
	qV = vc_dem;
}
// il CRE prende la sua % di portata + eccesso di CV
elena2indu = (double) 0 ;
if ( ((q-qTur)*aR + villo2elena) < qminR) {
	elena2indu = (q-qTur)*aR + villo2elena;
	qCRE = (double) 0 ; }
else { qCRE = ((q-qTur)*aR + villo2elena) ;}
if ( (((q-qTur)*aR + villo2elena) < opeR) && ( ((q-qTur)*aR + villo2elena) >= qminR ) ) {
	qCRE = re_dem - maxRe; }
else if ( ((q-qTur)*aR + villo2elena) > re_dem ) {
	elena2indu = ((q-qTur)*aR + villo2elena) - re_dem;
	qCRE = re_dem; }
// il CI prende la sua % di portata piu' eccesso di CRE
if ( ((q-qTur)*aI + elena2indu + qTur) > ic_dem ) {
	qIC = ic_dem; }
else { qIC = ((q-qTur)*aI + elena2indu + qTur); }

qout.push_back( qCRE );
qout.push_back( qV );
qout.push_back( qIC );
qout.push_back( (double) 4 );
}

return qout;
}

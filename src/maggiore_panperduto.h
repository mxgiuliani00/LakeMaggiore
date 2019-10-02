/*
 * maggiore_panperduto.h
 *
 */

#ifndef maggiore_panperduto_H_
#define maggiore_panperduto_H_

#include <math.h>
#include <vector>
#include <string>
#include "utils.h"


namespace std{



class maggiore_panperduto {
public:
	maggiore_panperduto();
	virtual ~maggiore_panperduto();


	vector<vector<double> > Distribution_rule( vector<double> q, vector<double> d ,vector<double> dREC, vector<double> dVC, vector<double> dIC, vector<double> doy);
	vector<double> Normal_Condition ( double re_dem, double vc_dem, double ic_dem);
	vector<double> Hydro_reduction ( double q, double tot_dem, double re_dem, double vc_dem, double ic_dem, double maxRe, double maxIc, double qThreR, double qThreI, double opeR, double deficit);
	vector<double> Irrig_reduction( double q, double tot_dem, double re_dem, double vc_dem, double ic_dem, double maxRe, double maxIc );
	vector<double> Conflict (double q, double tot_dem, double re_dem, double vc_dem, double ic_dem, double qTur, double aR, double aI, double aV, double qminV, double opeV, double qminR, double opeR, double maxRe);




protected:
	vector<double> a_ic, a_rec, a_vc, q_turbigo;
	vector<double> max_rid_ic, max_rid_rec, max_rid_vc;
	vector<double> q_min_rec, q_min_vc;
	vector<double> q_ope_rec, q_ope_vc;
	vector<double> q_thre_ic, q_thre_rec;
	vector<double> DMV ;
};


}
#endif

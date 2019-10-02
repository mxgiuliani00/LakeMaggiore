/*
 * model_lakemaggiore.h
 *
 *  Created on: 20/may/2016
 *      Author: MatteoG
 */

#ifndef model_lakemaggiore_H_
#define model_lakemaggiore_H_

#include <math.h>
#include <vector>
#include <string>
#include "lake.h"
#include "lakemaggiore.h"
#include "catchment.h"
#include "param_function.h"
#include "rbf.h"
#include "ann.h"
#include "pwLinear.h"
#include "ncRBF.h"
#include "maggiore_panperduto.h"

namespace std{

class model_lakemaggiore {
public:
        model_lakemaggiore();
        virtual ~model_lakemaggiore();

        model_lakemaggiore(string filename);
        void clear_model_lakemaggiore();

        /**
        * number of objectives and variables
        */
        int getNobj();
        int getNvar();

        /**
         * function to perform the optimization:
         *  - var = decision variables
         *  - obj = objectives
        */
        void evaluate(double* var, double* obj);

protected:
        void readFileSettings(string filename);

        // problem setting
        int Nsim;               // number of simulation (1=deterministic, >1 MC)
        int NN;                 // dimension of the stochastic ensemble
        int T;                  // period
        int integStep;          // integration timestep = number of subdaily steps
        int H;                  // simulation horizon
        int Nobj;               // number of objectives
        int Nvar;               // number of variables
        int initDay;            // first day of simulation
        vector<double> doy_file;     // day of the year (it includes leap years, otherwise doy is computed runtime in the simulation)

        // catchment
        catchment_param Maggiore_catch_param;
        catchment* MaggioreCatchment;


        // reservoir: Lake Maggiore
        reservoir_param Maggiore_param;
        lake* LakeMaggiore;

        // diversion:panperduto
        maggiore_panperduto* PanperdutoDiv;

        // operating policy
        pFunction_param p_param;
        param_function* mPolicy;

        // objective function data

        // flood
        int warmup;                                 // number of days of warmup before obj calculation starts
        vector<vector<double> > level_areaFlood;    // level (cm) - flooded area (m2)
        double h_flo;                               // flooding threshold

        // agricultural
        vector<double> demand;                          // total downstream demand (m3/s)
        vector<double> rec_demand;                      // Regina Elena demand (m3/s)
        vector<double> vc_demand;                       // Villoresi demand (m3/s)
        vector<double> ic_demand;                       // Canale Industriale demand (m3/s)

        // hydropower production
        vector<double> w_rec_hydro;                   // Regina Elena hydropower demand
        vector<double> AN_demand;                     // Alto Novarese demand
        vector<double> hydro;                         // hydropower parameters
        vector<vector< double> > pp_Enel;             // Enel power plants parameters
        vector< vector < double> > pp_Sesia;          // Est Sesia power plants parameters

        // initial condition

        double inflow00;    // previous day inflow
        double Qmio;		// previous day release
        double hSC;			// h Sesto Calende


        // results



        /**
         * function to perform the simulation over the scenario ps
         */
         vector<double> simulate(int ps);

        /**
         * Functions to compute the objective functions:
         **/
        double floodDays(vector<double> h, vector<double> r, double Q0, double h_flo);
        vector<double> g_hFlood(vector<double> h, vector<double> r, double Q0, double h_flo);
        vector<vector<double> > g_deficit(vector<double> dR, vector<double> dV, vector<double> dI, vector<double> qR, vector<double> qV, vector<double> qI);
        vector<vector<double> > g_deficit_hydro(vector<double> qR,vector<double> w_rec_hydro, vector<double> qI, vector<double> r, vector<double> Nov_dem, vector<double> hydro_param, vector<vector<double> >pp_Enel_param, vector<vector<double> > pp_Sesia_param, vector<double> doy );
        vector<double> g_deficitTOT(vector<double> q, vector<double> w, vector<double> doy);

        /**
         * function to convert lake level into level at Sesto Calende (and vice-versa)
         **/
        double SCtoLakeLevel(double hSC, double qMio);
        double lakeLevelToSC(double hL, double r);


};

}

#endif /* model_lakemaggiore_H_ */

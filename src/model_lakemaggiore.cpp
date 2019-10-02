/*
 * model_lakemaggiore.cpp
 *
 *  Created on: 20/may/2016
 *      Author: MatteoG
 */


#include "model_lakemaggiore.h"
#include "utils.h"
#include <fstream>
#include <iostream>
#include <string>
#include <stdlib.h>

using namespace std;

#define PI 3.14159265

model_lakemaggiore::model_lakemaggiore() {
    // TODO Auto-generated constructor stub
}

model_lakemaggiore::~model_lakemaggiore() {
    // TODO Auto-generated destructor stub
}

model_lakemaggiore::model_lakemaggiore(string filename){

    readFileSettings(filename);
    // create catchment
    catchment* mc = new catchment(Maggiore_catch_param);
    MaggioreCatchment = mc;
    // create lake
    lake* ml = new lakemaggiore();
    LakeMaggiore = ml;
    LakeMaggiore->setEvap(0);
    Maggiore_param.minEnvFlow.filename = "../data/MEF_maggiore.txt";
    Maggiore_param.minEnvFlow.row = T;
    LakeMaggiore->setMEF(Maggiore_param.minEnvFlow);
    LakeMaggiore->setSurface(207510000);
    LakeMaggiore->setInitCond(Maggiore_param.initCond);
    Maggiore_param.lsv_rel.filename = "../data/maggiore_levelVolume.txt";
    Maggiore_param.lsv_rel.row = 2;
    Maggiore_param.lsv_rel.col = 111;
    LakeMaggiore->setLSV_Rel(Maggiore_param.lsv_rel);    

	  //create panperduto
	  maggiore_panperduto* md= new maggiore_panperduto();
	  PanperdutoDiv = md;


    // policy
    switch (p_param.tPolicy) {
    case 1: // RBF policy
    {
        param_function* mp1 = new rbf(p_param.policyInput,p_param.policyOutput,p_param.policyStr);
        mPolicy = mp1;
        break;
    }
    case 2: // ANN
    {
        param_function* mp2 = new ann(p_param.policyInput,p_param.policyOutput,p_param.policyStr);
        mPolicy = mp2;
        break;
    }
    case 3: // piecewise linear policy
    {
        param_function* mp3 = new pwLinear(p_param.policyInput,p_param.policyOutput,p_param.policyStr);
        mPolicy = mp3;
        break;
    }
    case 4:
    {
        param_function* mp4 = new ncRBF(p_param.policyInput,p_param.policyOutput,p_param.policyStr);
        mPolicy = mp4;
        break;
    }
    default:
        break;
    }
    // min-max policy input
    mPolicy->setMaxInput(p_param.MIn); mPolicy->setMaxOutput(p_param.MOut);
    mPolicy->setMinInput(p_param.mIn); mPolicy->setMinOutput(p_param.mOut);


    // OBJECTIVES
    h_flo = 1.73 ; // hsc threshold
    demand = utils::loadVector("../data/maggioreDemand.txt", 365);
    rec_demand = utils::loadVector("../data/maggioreREC_Demand.txt", 365);
    vc_demand = utils::loadVector("../data/maggioreVC_Demand.txt", 365);
    ic_demand = utils::loadVector("../data/maggioreIC_Demand.txt", 365);
    AN_demand = utils::loadVector("../data/AltoNovarese_demand.txt", 365);
    hydro = utils::loadVector("../data/hydro_parameters.txt", 4);
    pp_Enel = utils::loadMatrix("../data/hydro_powerplant_Enelparam.txt", 4, 3); //matrix that includes the ENEL power plants parameters (minimum turb flow, hydraulic haed..)
    pp_Sesia = utils::loadMatrix("../data/hydro_powerplant_EstSesiaparam.txt", 6, 3); //matrix that includes the Consorzio Est Sesia power plants parameters (minimum turb flow, hydraulic haed..)
    w_rec_hydro = utils::loadVector("../data/w_rec_hyd.txt", 365);

}

void model_lakemaggiore::clear_model_lakemaggiore(){
    delete LakeMaggiore;
    delete mPolicy;
    delete MaggioreCatchment;
}


int model_lakemaggiore::getNobj() {
	return Nobj;
}

int model_lakemaggiore::getNvar() {
        return Nvar;
}


void model_lakemaggiore::evaluate(double* var, double* obj){

    // set CONTROL POLICY
 mPolicy->setParameters(var);

    vector<double> J ;

    if(Nsim < 2){ // single simulation
        J = simulate(0);
        for(unsigned int i=0; i<Nobj; i++){
            obj[i] = J[i];
        }
    }else{ 
        // MC simulation to be implemented
    }

   mPolicy->clearParameters();

}



vector<double> model_lakemaggiore::simulate(int ps){


    // INITIALIZATION: storage, level, decision, release
    vector<double> s (H+1,-999) ;
    vector<double> h (H+1,-999) ;
    vector<double> u (H,-999) ;
    vector<double> r (H+1,-999) ;
    vector<double> doy (H,-999) ;
    vector<double> rPan (H,-999);

    // simulation variables
    double qIn, qIn_1, Q;   // daily inflow (today and yesterday)
    vector<double> sh_rh;   // storage and release resulting from hourly integration
    vector<double> uu;      // decision vector
    vector<double> input;   // policy input vector
    vector<double> JJ;      // objective vector

    // initial condition

    Q = Qmio;
    double hSC = LakeMaggiore->getInitCond();
    h[0] = SCtoLakeLevel( hSC, Q );
    s[0] = LakeMaggiore->levelToStorage(h[0]);
    qIn_1 = inflow00 ;

    //cout << "### start simulation ###" << endl;
    // Run simulation:
    for(unsigned int t = 0; t < H; t++){
        // day of the year
        doy[t] =  (initDay+t-1)%T+1;  // day of the year
        
        // inflows
        qIn = MaggioreCatchment->getInflow(t, ps);
        // compute decision - standard
        input.push_back(sin(2*PI*doy[t])/365) ;
        input.push_back(cos(2*PI*doy[t])/365) ;
        input.push_back( h[t] );
        
        uu = mPolicy->get_NormOutput(input);
        u[t] = uu[0]; // single release decision

        // hourly integration of mass balance
        sh_rh = LakeMaggiore->integration(integStep,t,s[t],u[t],qIn,doy[t],ps);

        // assignment of daily values
        s[t+1] = sh_rh[0];
        r[t+1] = sh_rh[1];
        h[t+1] = LakeMaggiore->storageToLevel(s[t+1]);
        qIn_1 = qIn;

        // clear subdaily values
        sh_rh.clear();
        input.clear();
        uu.clear();
    }
    rPan = r;
    rPan.erase(rPan.begin());
    
	// compute panperduto allocation
  vector<vector<double> > qAlloc = PanperdutoDiv -> Distribution_rule( rPan, demand, rec_demand, vc_demand, ic_demand, doy );
  vector<double> qReginaElena, qVilloresi, qIndustriale, idx;

  for ( int i=0; i<H; i++){
    qReginaElena.push_back(qAlloc[i][0]);
    qVilloresi.push_back(qAlloc[i][1]);
    qIndustriale.push_back(qAlloc[i][2]);
    idx.push_back(qAlloc[i][3]);
  }

    // compute objectives
    h.erase(h.begin()); r.erase(r.begin());
    h.erase(h.begin(),h.begin()+warmup);
    r.erase(r.begin(),r.begin()+warmup);
    doy.erase(doy.begin(),doy.begin()+warmup);

    // flood days h vs deficit2
    vector<double> g_flo;
    vector<double> g_irr;
    vector<vector<double> > g_hydro;
    //g_flo = g_hFlood(h, r, Q, h_flo);
    int NYears = H/T; // number of years

    g_irr = g_deficitTOT(r, demand, doy);
    
    g_hydro = g_deficit_hydro (qReginaElena, w_rec_hydro, qIndustriale, r, AN_demand, hydro, pp_Enel, pp_Sesia, doy);
    vector<double> def_ENEL, def_ESTS;
    def_ENEL = g_hydro[0];
    def_ESTS = g_hydro[1];
    
    //JJ.push_back( utils::computeMean(g_flo) );
    JJ.push_back( floodDays( h, r, Q, h_flo )/NYears );
    JJ.push_back( utils::computeMean(g_irr) );

    return JJ;
}


vector<double> model_lakemaggiore::g_deficitTOT(vector<double> q, vector<double> w, vector<double> doy){

  vector<double> gt;
  double d, qdiv;
  for(unsigned int i=0; i<q.size(); i++){
      qdiv = q[i] - LakeMaggiore->getMEF(doy[i]-1);
      if( qdiv<0.0 ){
          qdiv = 0.0;
      }
      d = w[doy[i]-1] - qdiv;
      if( d < 0.0 ){
          d = 0.0;
      }
      gt.push_back( d*d );
  }
  return gt;

}


double model_lakemaggiore::floodDays(vector<double> h, vector<double> r, double Q0, double h_flo){

    r.insert(r.begin(),Q0);
    double c=0.0, hSC = 0.0;
    for(unsigned int i=0; i<h.size(); i++){
        hSC = lakeLevelToSC( h[i], r[i]) ;
        if(hSC>h_flo){
            c=c+1.0;
        }
    }
    return c;
}

vector<double> model_lakemaggiore::g_hFlood( vector<double> h, vector<double> r, double Q0, double h_flo) {

    r.insert(r.begin(),Q0);
    vector<double> h_flood;
    double hSC ;
    for (unsigned int i = 0; i<h.size(); i++){
        hSC = lakeLevelToSC( h[i], r[i]) ;
        if (hSC > h_flo) {
            h_flood.push_back(hSC-h_flo);
        }else{
            h_flood.push_back(0.0);
        }
    }

    return h_flood;
}

// given Regina Elena, Villoresi and Industriale inflows it computes the deficit for each canal
vector<vector<double> > model_lakemaggiore::g_deficit(vector<double> dR, vector<double> dV, vector<double> dI, vector<double> qR, vector<double> qV, vector<double> qI){
    vector <double> g_irrR, g_irrV, g_irrI;
    double R = 0.0, V = 0.0, I = 0.0;
    for(unsigned int i = 0; i<dR.size(); i ++){
    R = dR[i]-qR[i];
    g_irrR.push_back(R);
    V = dV[i]-qV[i];
    g_irrV.push_back(V);
    I = dI[i]-qI[i];
    g_irrI.push_back(I);
    }

    vector<vector<double> > deficit;
    deficit.push_back(g_irrR);
    deficit.push_back(g_irrV);
    deficit.push_back(g_irrI);

    return deficit;
}



vector<vector<double> > model_lakemaggiore::g_deficit_hydro(vector<double> qR, vector<double> w_rec_hydro, vector<double> qI, vector<double> r, vector<double> Nov_dem, vector<double> hydro_param, vector<vector<double> >pp_Enel_param, vector<vector<double> > pp_Sesia_param, vector<double> doy){
    // define all the parameters needed_ENEL
    double gam, eta, fi, g;
    double q_min_PT, q_max_PT, q_min_V, q_max_V, q_min_T, q_max_T, q_min_TS, q_max_TS, H_PT, H_V, H_T, H_TS;
    // define all the parameters needed_SESIA
    double q_min_SES, q_max_CAV, H_CAV, H_URI, H_CID, H_MONT, H_TERD, H_VERS;
    // set parameters values
    gam = hydro_param[0];
    eta = hydro_param[1];
    fi = hydro_param[2];
    g = hydro_param[3];

    q_min_PT = pp_Enel_param[0][1]; q_max_PT = pp_Enel_param[0][2]; q_min_V= pp_Enel_param[1][1]; q_max_V= pp_Enel_param[1][2]; q_min_T= pp_Enel_param[2][1]; q_max_T= pp_Enel_param[2][2]; q_min_TS= pp_Enel_param[3][1]; q_max_TS= pp_Enel_param[3][2]; H_PT= pp_Enel_param[0][0]; H_V= pp_Enel_param[1][0]; H_T= pp_Enel_param[2][0]; H_TS= pp_Enel_param[3][0];
    q_min_SES = pp_Sesia_param[0][1]; q_max_CAV = pp_Sesia_param[0][2]; H_CAV = pp_Sesia_param[0][0]; H_URI = pp_Sesia_param[1][0]; H_CID = pp_Sesia_param[2][0]; H_MONT = pp_Sesia_param[3][0]; H_TERD = pp_Sesia_param[4][0]; H_VERS = pp_Sesia_param[5][0];

    //compute effective flow for each power plant
    vector<double> qPT, qV, qT, qTS, qSESIA;
    double q_Porto_Torre, q_Sesia_pp;

    // ENEL Porto Torre
    for(unsigned int i = 0; i<r.size(); i ++){
     q_Porto_Torre = r[i] - qR[i];
     if ( q_Porto_Torre< q_min_PT) { qPT.push_back( (double) 0 );}
     else if ( q_Porto_Torre< q_max_PT ) { qPT.push_back(q_Porto_Torre) ;}
     else { qPT.push_back(q_max_PT) ;}
    }

    //ENEL Vizzola
    for(unsigned int i = 0; i<r.size(); i ++){
     if ( qI[i]< q_min_V) { qV.push_back( (double) 0 );}
     else if ( qI[i]< q_max_V ) { qV.push_back(qI[i]) ;}
     else { qV.push_back(q_max_V) ;}
    }

    //ENEL Tornavento
    for(unsigned int i = 0; i<r.size(); i ++){
     if ( qI[i]< q_min_T) { qT.push_back( (double) 0 );}
     else if ( qI[i]< q_max_T ) { qT.push_back(qI[i]) ;}
     else { qT.push_back(q_max_T) ;}
    }

    //ENEL Turbigo Superiore
    for(unsigned int i = 0; i<r.size(); i ++){
     if ( qI[i]< q_min_TS) { qTS.push_back( double(0)) ;}
     else if ( qI[i]< q_max_TS ) { qTS.push_back(qI[i]) ;}
     else { qTS.push_back(q_max_TS) ;}
    }

    //EST SESIA POWER PLANTS (we need just one of them given that every est sesia power plant has the same qmin,qmax)

    for(unsigned int i = 0; i<r.size(); i ++){
     q_Sesia_pp = qR[i] - Nov_dem[i]; // Alto Novarese demand is subtracted to the flow avaiable for the est sesia power plants
     if ( q_Sesia_pp< q_min_SES) { qSESIA.push_back( (double) 0 );}
     else if ( q_Sesia_pp< q_max_CAV ) { qSESIA.push_back(q_Sesia_pp );}
     else { qSESIA.push_back(q_max_CAV) ;}
    }

    //compute deficit_hydro
    vector <double> g_hydro_ENEL, g_hydro_SESIA;
    double enel=0.0, est_sesia=0.0;
    double deficitPT = 0.0, deficitV = 0.0, deficitT = 0.0, deficitTS = 0.0;
    double deficit_est_s = 0.0;
    // verify when a deficit occurs (ENEL power plants)
    for(unsigned int i = 0; i<r.size(); i ++){
    deficitPT = q_max_PT - qPT[i];
    deficitV = q_max_V - qV[i];
    deficitT = q_max_T - qT[i];
    deficitTS = q_max_TS - qTS[i];
    // compute the difference between energy that could be produced and the one produced (ENEL power plants)
    enel = gam * eta * g * fi *( H_PT*deficitPT + H_V*deficitV + H_T*deficitT + H_TS *deficitTS );
    g_hydro_ENEL.push_back(enel);
    // verify when a deficit occurs (Est sesia power plants, nb: the estSesia hydo demand is time-varying)
    deficit_est_s = w_rec_hydro[doy[i]-1] - qSESIA[i];
    if (deficit_est_s < 0) {deficit_est_s=0;}
    // compute the difference between energy that could be produced and the one produced (Est sesia power plants)
    est_sesia =  gam * eta * g * fi * deficit_est_s*( H_CAV+ H_URI + H_CID + H_MONT + H_TERD + H_VERS );
    g_hydro_SESIA.push_back(est_sesia);
    }

    //save deficit
    vector<vector<double> > deficit_hydro;
    //cout << "g_hydro_Enel="<< g_hydro_ENEL.size() << endl;
    //cout << "g_hydro_Est="<< g_hydro_SESIA.size() << endl;
    deficit_hydro.push_back(g_hydro_ENEL);
    deficit_hydro.push_back(g_hydro_SESIA);
    return deficit_hydro;

}


double model_lakemaggiore::SCtoLakeLevel(double hSC, double qMio){

    double h = hSC + 78.74E-6 * qMio -0.01045;
    return h;
}

double model_lakemaggiore::lakeLevelToSC(double hL, double r){

    double hSC = hL - 78.74E-6 * r +0.01045;
    return hSC;
}





void model_lakemaggiore::readFileSettings(string filename){

    ifstream in;
    string sJunk = "";

    in.open(filename.c_str(), ios_base::in);
    if(!in)
    {
        cout << "The input file specified: " << filename << " could not be found!" << endl;
        exit(1);
    }

    // PROBLEM SETTING
    //Look for the <NUM_SIM> key
    while (sJunk != "<NUM_SIM>")
    {
        in >> sJunk;
    }
    in >> Nsim;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <DIM_ENSEMBLE> key
    while (sJunk != "<DIM_ENSEMBLE>")
    {
        in >> sJunk;
    }
    in >> NN;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <PERIOD> key
    while (sJunk != "<PERIOD>")
    {
        in >> sJunk;
    }
    in >> T;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <INTEGRATION> key
    while (sJunk != "<INTEGRATION>")
    {
        in >> sJunk;
    }
    in >> integStep;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <SIM_HORIZON> key
    while (sJunk != "<SIM_HORIZON>")
    {
        in >> sJunk;
    }
    in >> H;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <NUM_OBJ> key
    while (sJunk != "<NUM_OBJ>")
    {
        in >> sJunk;
    }
    in >> Nobj;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <NUM_VAR> key
    while (sJunk != "<NUM_VAR>")
    {
        in >> sJunk;
    }
    in >> Nvar;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <WARMUP> key
    while (sJunk != "<WARMUP>")
    {
        in >> sJunk;
    }
    in >> warmup;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <DOY> key
    string file0;
    int temp=0;
    while (sJunk != "<DOY>")
    {
        in >> sJunk;
    }
    in >> temp;
    if(temp>0){
        initDay = temp;
    }else{
        in.ignore(1000,'\n');
        in >> file0;
        doy_file = utils::loadVector(file0,H);
    }
    //Return to the beginning of the file
    in.seekg(0, ios::beg);


    // CATCHMENT MODEL
    //Look for the <CATCHMENT> key
    while (sJunk != "<CATCHMENT>")
    {
        in >> sJunk;
    }
    in >> Maggiore_catch_param.CM;
    in.ignore(1000,'\n');
    in >> Maggiore_catch_param.inflow_file.filename;
    Maggiore_catch_param.inflow_file.row = NN;
    Maggiore_catch_param.inflow_file.col = H;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    // MODEL INITIALIZATION
    //Look for the <INIT_CONDITION> key
    while (sJunk != "<INIT_CONDITION>")
    {
        in >> sJunk;
    }
    in >> Maggiore_param.initCond;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <INIT_INFLOW> key
    while (sJunk != "<INIT_INFLOW>")
    {
        in >> sJunk;
    }
    in >> inflow00;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

	//Look for the <INIT_QMIO> key
    while (sJunk != "<INIT_QMIO>")
    {
        in >> sJunk;
    }
    in >> Qmio;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <POLICY_CLASS> key
    while (sJunk != "<POLICY_CLASS>")
    {
        in >> sJunk;
    }
    in >> p_param.tPolicy;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <NUM_INPUT> key
    double i1, i2;
    while (sJunk != "<NUM_INPUT>")
    {
        in >> sJunk;
    }
    in >> p_param.policyInput;
    in.ignore(1000,'\n');
    //Loop through all of the input data and read in this order:
    for (int i=0; i<p_param.policyInput; i++)
    {
        in >> i1 >> i2;
        p_param.mIn.push_back(i1);
        p_param.MIn.push_back(i2);
        in.ignore(1000,'\n');
    }
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Look for the <NUM_OUTPUT> key
    double o1, o2;
    while (sJunk != "<NUM_OUTPUT>")
    {
        in >> sJunk;
    }
    in >> p_param.policyOutput;
    in.ignore(1000,'\n');
    //Loop through all of the input data and read in this order:
    for (int i=0; i<p_param.policyOutput; i++)
    {
        in >> o1 >> o2;
        p_param.mOut.push_back(o1);
        p_param.MOut.push_back(o2);
        in.ignore(1000,'\n');
    }
    //Return to the beginning of the file
    in.seekg(0, ios::beg);


    //Look for the <POLICY_STRUCTURE> key
    while (sJunk != "<POLICY_STRUCTURE>")
    {
        in >> sJunk;
    }
    in >> p_param.policyStr;
    //Return to the beginning of the file
    in.seekg(0, ios::beg);

    //Close the input file
    in.close();

    return;

}

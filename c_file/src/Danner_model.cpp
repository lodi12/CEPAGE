/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   HH_model.cpp
 * Author: teo
 *
 * Created on July 12, 2017, 10:33 AM
 */

#include "../inc/Danner_model.hpp"

Danner_model::Danner_model()
{
    this->nx = 2;
    this->C = 0;
    this->gNaP = 0;
    this->ENa = 0;
    this->gl = 0;
    this-> El = 0;
    this->VhalfM = 0;
    this->km = 0;
    this->Vhalfh = 0;
    this->kh = 0;
    this->tau0 = 0;
    this->tauMax = 0;
    this->VhalfTau = 0;
    this->kTau = 0;
    this->Di = 0;
    this->gSynE = 0;
    this->EsynE = 0;
}

Danner_model::Danner_model(double C, double gNaP, double ENa, double gl, double  El, double VhalfM, double km, double Vhalfh, double kh, double tau0, double tauMax, double VhalfTau, double kTau, double Di, double gSynE, double EsynE)
{
    this->nx = 2;
    this->C = C;
    this->gNaP = gNaP;
    this->ENa = ENa;
    this->gl = gl;
    this->El = El;
    this->VhalfM = VhalfM;
    this->km = km;
    this->Vhalfh = Vhalfh;
    this->kh = kh;
    this->tau0 = tau0;
    this->tauMax = tauMax;
    this->VhalfTau = VhalfTau;
    this->kTau = kTau;
    this->Di = Di;
    this->gSynE = gSynE;
    this->EsynE = EsynE;
}


void Danner_model::getXdot(double t, double *x, double *xdot,double Iext)
{
    double nx = this->nx;
    double C = this->C;
    double gNaP = this->gNaP;
    double ENa = this->ENa;
    double gl = this->gl;
    double El = this->El;
    double VhalfM = this->VhalfM;
    double km = this->km;
    double Vhalfh = this->Vhalfh;
    double kh = this->kh;
    double tau0 = this->tau0;
    double tauMax = this->tauMax;
    double VhalfTau = this->VhalfTau;
    double kTau = this->kTau;
    double Di = this->Di;
    double gSynE = this->gSynE;
    double EsynE = this->EsynE;
    
    double V = x[0];
    double h = x[1];
    double m = 1.0/(1.0+exp((V-VhalfM)/km));
    double tauH = tau0 +((tauMax-tau0)/cosh((V-VhalfTau)/kTau));
    double hInf = 1.0/(1+exp((V-Vhalfh)/kh));
    double Idrive = gSynE*(V-EsynE)*Di;
    xdot[1] = 1.0/tauH*(hInf-h);
    xdot[0] = 1.0/C*(-gl*(V-El)-gNaP*m*h*(V-ENa)-Idrive+Iext);
    
}

Danner_model::~Danner_model() {
}


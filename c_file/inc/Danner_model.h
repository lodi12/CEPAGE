/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   Danner_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:33 AM
 */

#ifndef DANNER_MODEL_H
#define DANNER_MODEL_H

#include "neuron_model.h"
#include <math.h>

class Danner_model : public neuron_model {
    
private:
    double C;
    double gNaP;
    double ENa;
    double gl;
    double  El;
    double VhalfM;
    double km;
    double Vhalfh;
    double kh;
    double tau0;
    double tauMax;
    double VhalfTau;
    double kTau;
    double Di;
    double gSynE;
    double EsynE;
    
public:
    Danner_model();
    Danner_model(double C, double gNaP, double ENa, double gl, double  El, double VhalfM, double km, double Vhalfh, double kh, double tau0, double tauMax, double VhalfTau, double kTau, double Di, double gSynE, double EsynE);
    void getXdot(double t, double *x, double *xdot,double Iext);
    
    virtual ~Danner_model();
    
};

#endif /* DANNER_MODEL_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ML_model.cpp
 * Author: teo
 * 
 * Created on July 12, 2017, 10:21 AM
 */

#include "../inc/ML_model.h"

ML_model::ML_model() 
{
    this->nx = 2;
    this->CM = 0;
    this->gCa = 0;
    this->V3 = 0;
    this->V4 = 0;
    this->phi = 0;
    this->gl = 0;
    this->Vl = 0;
    this->I = 0;
}

ML_model::ML_model(double CM,double gCa,double V3,double V4,double phi,double gl,double Vl,double I)
{
    this->nx = 2;
    this->CM = CM;
    this->gCa = gCa;
    this->V3 = V3;
    this->V4 = V4;
    this->phi = phi;
    this->gl = gl;
    this->Vl = Vl;
    this->I = I;
}

void ML_model::getXdot(double t, double *x, double *xdot,double *Iext)
{
double gCa = this->gCa;
double I = this->I;
double V3 = this->V3;
 double V4 = this->V4;
double phi = this->phi;
double Vl = this->Vl;
double gl = this->gl;

double CM = this->CM;

const double gK = 8;;
const double VCa = 120;
const double Vk = -80;
const double V1 = -1.2;
const double V2 = 18;

double Minf = 0.5*(1 + tanh((x[0] - V1 )/V2 ));
double Ninf = 0.5*(1 + tanh((x[0] - V3 )/V4 ));
double lamda_N = phi*cosh((x[0] - V3 )/(2*V4) );

xdot[1] = lamda_N*(Ninf-x[1]);
xdot[0] = (-gl*(x[0] - Vl) - gCa*Minf*(x[0] - VCa) - gK*x[1]*(x[0] - Vk ) + I + Iext[0])/CM;

}

ML_model::~ML_model() 
{
}


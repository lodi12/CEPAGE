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

#include "../inc/HH_model.hpp"

HH_model::HH_model() 
{
    this->nx = 3;
    this->gna = 0; 
    this->ENa = 0; 
    this->gk2 = 0;
    this->Ek = 0; 
    this->gl = 0;
    this->El = 0; 
    this->tNa = 0;
    this->tk2 = 0;
    this->C = 0; 
    this->Iapp = 0;
    this->VshiftK2 = 0; 
}

HH_model::HH_model(const HH_model &h)
{
    this->nx = h.nx;
    this->gna = h.gna; 
    this->ENa = h.ENa; 
    this->gk2 = h.gk2;
    this->Ek = h.Ek; 
    this->gl = h.gl;
    this->El = h.El; 
    this->tNa = h.tNa;
    this->tk2 = h.tk2;
    this->C = h.C; 
    this->Iapp = h.Iapp;
    this->VshiftK2 = h.VshiftK2; 
}

HH_model::HH_model(double gna,double ENa,double gk2,double Ek,double gl,double El,double tNa,double tk2,double C,double Iapp,double VshiftK2)
{
    this->nx = 3;
    this->gna = gna; 
    this->ENa = ENa; 
    this->gk2 = gk2;
    this->Ek = Ek; 
    this->gl = gl;
    this->El = El; 
    this->tNa = tNa;
    this->tk2 = tk2;
    this->C = C; 
    this->Iapp = Iapp;
    this->VshiftK2 = VshiftK2; 
}

HH_model *HH_model::clone() const
{
    return new HH_model(*this);
}


void HH_model::getXdot(double t, double *x, double *xdot,double Iext)
{
double gna = this->gna;
double ENa = this->ENa;
double gk2 = this->gk2;
double Ek = this->Ek;
double gl = this->gl;
double El = this->El;
double tNa = this->tNa;
double tk2 = this->tk2;
double C = this->C;
double Iapp = this->Iapp;
double VshiftK2 = this->VshiftK2;


double hInf = 1/(1+exp(500*(x[0]+0.0333)));
double nInf = 1/(1+exp(-150*(x[0]+0.0305)));
double mInf = 1/(1+exp(-83*(x[0]+0.018+VshiftK2)));

double INa = gna*nInf*nInf*nInf*x[1]*(x[0]-ENa);
double Ik2 = gk2*x[2]*x[2]*(x[0]-Ek);
double Il = gl*(x[0]-El);

xdot[0] = (-INa-Ik2-Il-Iapp+Iext)/C;
xdot[1] = (hInf- x[1])/tNa;
xdot[2] = (mInf- x[2])/tk2;

}

HH_model::~HH_model() {
}


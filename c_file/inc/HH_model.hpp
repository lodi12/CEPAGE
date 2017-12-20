/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HH_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:33 AM
 */

#ifndef HH_MODEL_H
#define HH_MODEL_H

#include "neuron_model.hpp"
#include <math.h>

class HH_model : public neuron_model {
    
private:
    double gna; 
    double ENa; 
    double gk2;
    double Ek; 
    double gl;
    double El; 
    double tNa;
    double tk2;
    double C; 
    double Iapp;
    double VshiftK2; 
    
public:
    HH_model();
    HH_model(double gna,double ENa,double gk2,double Ek,double gl,double El,double tNa,double tk2,double C,double Iapp,double VshiftK2);
    
    void getXdot(double t, double *x, double *xdot,double Iext);
    
    virtual ~HH_model();

};

#endif /* HH_MODEL_H */


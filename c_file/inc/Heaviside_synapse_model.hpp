/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   Heaviside_synapse_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:21 AM
 */

#ifndef HEAVISIDE_SYNAPSE_MODEL_H
#define HEAVISIDE_SYNAPSE_MODEL_H

#include "synapse_model.hpp"
#include <math.h>

class Heaviside_synapse_model : public synapse_model {
    
private:
    double theta;
    
public:
    Heaviside_synapse_model();
    Heaviside_synapse_model(double theta);
    
    void getXdot(double t, double *x, double *xdot,double Vpre);
    double getActivation(double *x, double Vpre);
    
    virtual ~Heaviside_synapse_model();

};


#endif /* HEAVISIDE_SYNAPSE_MODEL_H */


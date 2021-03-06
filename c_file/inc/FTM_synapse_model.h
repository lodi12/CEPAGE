/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FTM_synapse_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:21 AM
 */

#ifndef FTM_SYNAPSE_MODEL_H
#define FTM_SYNAPSE_MODEL_H

#include "synapse_model.h"
#include <math.h>
class FTM_synapse_model : public synapse_model {
    
private:
    double nu;
    double theta;
    
public:
    FTM_synapse_model();
    FTM_synapse_model(double nu,double theta);
    
    void getXdot(double t, double *x, double *xdot);
    double getActivation(double Vpre);
    
    virtual ~FTM_synapse_model();

};


#endif /* FTM_SYNAPSE_MODEL_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FTM_synapse_model.cpp
 * Author: picio
 * 
 * Created on 14 maggio 2017, 19.31
 */

#include "../inc/FTM_synapse_model.h"

FTM_synapse_model::FTM_synapse_model()
{
    this->nx = 0;
    this->nu = 0;
    this->theta = 0;
}

FTM_synapse_model::FTM_synapse_model(double nu, double theta)
{
    this->nx = 0;
    this->nu = nu;
    this->theta = theta;
}


double FTM_synapse_model::getActivation(double Vpre)
{
    double nu = this->nu;
    double theta = this->theta;    
    
    return 1.0/(1.0+exp(-nu*(Vpre-theta)));
}


void FTM_synapse_model::getXdot(double t, double *x, double *xdot)
{
    
}

FTM_synapse_model::~FTM_synapse_model() {
}


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

#include "../inc/Heaviside_synapse_model.hpp"

Heaviside_synapse_model::Heaviside_synapse_model()
{
    this->nx = 0;
    this->theta = 0;
}

Heaviside_synapse_model::Heaviside_synapse_model(double theta)
{
    this->nx = 0;
    this->theta = theta;
}


double Heaviside_synapse_model::getActivation(double *x,double Vpre)
{
    double theta = this->theta;    
    
    if(Vpre > theta)
        return 1.0;
    else
        return 0.0;
}


void Heaviside_synapse_model::getXdot(double t, double *x, double *xdot,double Vpre)
{
    
}

Heaviside_synapse_model::~Heaviside_synapse_model() {
}


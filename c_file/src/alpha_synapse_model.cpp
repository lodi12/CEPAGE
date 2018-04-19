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

#include "../inc/alpha_synapse_model.hpp"

alpha_synapse_model::alpha_synapse_model()
{
    this->nx = 1;
    this->alpha = 0;
    this->beta = 0;
    this->nu = 0;
    this->theta = 0;
}

alpha_synapse_model::alpha_synapse_model(const alpha_synapse_model &other)
{
    this->nx = other.nx;
    this->alpha = other.alpha;
    this->beta = other.beta;
    this->nu = other.nu;
    this->theta = other.theta;
}

alpha_synapse_model::alpha_synapse_model(double alpha, double beta, double nu, double theta)
{
    this->nx = 1;
    this->alpha = alpha;
    this->beta = beta;
    this->nu = nu;
    this->theta = theta;
}


alpha_synapse_model *alpha_synapse_model::clone() const
{
    return new alpha_synapse_model(*this);
}


double alpha_synapse_model::getActivation(double *x,double Vpre)
{
    return x[0];
}


void alpha_synapse_model::getXdot(double t, double *x, double *xdot,double Vpre)
{
    double nu = this->nu;
    double theta = this->theta;    
    double alpha = this->alpha;
    double beta = this->beta;
    
    xdot[0] = alpha*(1.0-x[0])/(1.0+exp(-nu*(Vpre-theta)))-beta*x[0];
}

alpha_synapse_model::~alpha_synapse_model() {
}


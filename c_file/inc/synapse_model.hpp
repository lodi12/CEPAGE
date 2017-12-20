/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   neuron_model.h
 * Author: picio
 *
 * Created on 14 maggio 2017, 19.28
 */

#ifndef SYNAPSE_MODEL_H
#define SYNAPSE_MODEL_H


class synapse_model
{        
    
protected:    
    int nx;
    
public:
    virtual double getActivation(double Vpre) = 0;
    virtual void getXdot(double t, double *x, double *xdot) = 0;
    virtual bool getResetConditions(double *x);
    virtual void resetStates(double *x);
    int getnx();
};

#endif /* NEURON_MODEL_H */


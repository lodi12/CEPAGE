/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   CPG.h
 * Author: picio
 *
 * Created on 16 maggio 2017, 18.41
 */

#ifndef CPG_H
#define CPG_H

#include "dynSys.hpp"
#include "neuron_model.hpp"
#include "synapse_model.hpp"

#include "stdlib.h"
#include "math.h"

class CPG : public dynSys {
private:
    int N;
    neuron_model **neuroni;
    double *g_in;
    double *g_ex;
    double *g_el;
    double EsynIn;
    double EsynEx; 

    synapse_model **inhActivation;
    synapse_model **excActivation;
    
    int *firstState;
    
    
public:
    CPG();
    CPG(int N,neuron_model **neuroni, double *g_in, double *g_ex, double *g_el,double EsynIn, double EsynEx, synapse_model **inhActivation, synapse_model **excActivation ); 
    
    void getXdot(double t, double *x, double *xdot,double Iext);
    bool getResetConditions(double *x);
    void resetStates(double *x);
    
    void getFirstIndex(int *firstIndex);
    
    virtual ~CPG();


};

#endif /* CPG_H */


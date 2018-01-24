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


#include <vector>

#include "dynSys.hpp"
#include "neuron_model.hpp"
#include "synapse_model.hpp"

#include "stdlib.h"
#include "math.h"

using namespace std;

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
    
    vector<int> **neuronDelaysIndex;
    vector<int> **inhSynDelaysIndex;
    vector<int> **excSynDelaysIndex;
 
    int *NneuronDelaysIndex;
    int *NinhSynDelaysIndex;
    int *NexcSynDelaysIndex;
    
    int *firstState;
    
    
public:
    CPG();
    CPG(int N,neuron_model **neuroni, double *g_in, double *g_ex, double *g_el,double EsynIn, double EsynEx, synapse_model **inhActivation, synapse_model **excActivation , double networkDelays[]); 
    
    /* Xprec is a Ndelays x Nstates matrix */
    void getXdot(double t, double *x, double *xdot,double Iext,double **Xprec);
    void getXdot(double t, double *x, double *xdot,double Iext);
    
    bool getResetConditions(double *x);
    
    void resetStates(double *x);
    
    void getFirstIndex(int *firstIndex);
    
    virtual ~CPG();


};

#endif /* CPG_H */


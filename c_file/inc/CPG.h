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

#include "dynSys.h"
#include "neuron_model.h"
#include "synapse_model.h"

#include "stdlib.h"
#include "math.h"



typedef struct SynStruct
{
    int i;
    int j;
    synapse_model *activationFunction;
    double g;
    
    SynStruct(int i, int j, synapse_model *activationFunction, double g) : i(i), j(j), activationFunction(activationFunction), g(g) {}
    
} t_SynStruct;


class CPG : public dynSys {
    
    
    
private:
    int N;
    int Ninh;
    int Nexc;
    neuron_model **neuroni;
    
    double *g_el;
    
    double EsynIn;
    double EsynEx;
    
    t_SynStruct **inhSyn;
    t_SynStruct **excSyn;
    
    int *firstState;
    
    
public:
    
    CPG();
    CPG(int N,neuron_model **neuroni, double *g_el,double EsynIn, double EsynEx, t_SynStruct **inhSyn, t_SynStruct **excSyn );
    
    void getXdot(double t, double *x, double *xdot,double *Iext);
    bool getResetConditions(double *x);
    void resetStates(double *x);
    
    void getFirstIndex(int *firstIndex);
    
    virtual ~CPG();
    
    
};

#endif /* CPG_H */


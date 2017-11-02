/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   CPG.cpp
 * Author: picio
 *
 * Created on 16 maggio 2017, 18.41
 */

#include "../inc/CPG.h"

CPG::CPG()
{
    this->N = 0;
    this->neuroni = (neuron_model **)malloc(0);
    this->g_in = (double *)malloc(0);
    this->g_ex = (double *)malloc(0);
    this->g_el = (double *)malloc(0);
    this->EsynIn = 0;
    this->EsynEx = 0;
    this->inhActivation = (synapse_model **)malloc(0);
    this->excActivation = (synapse_model **)malloc(0);
    
    this->firstState = (int *)malloc(0);
}

CPG::CPG(int N,neuron_model **neuroni, double *g_in, double *g_ex, double *g_el,double EsynIn, double EsynEx, synapse_model **inhActivation, synapse_model **excActivation )
{
    int i,index;
    this->N = N;
    this->neuroni = neuroni;
    this->EsynIn = EsynIn;
    this->EsynEx = EsynEx;

    
    this->g_in = (double *)malloc(N*N*sizeof(double));
    this->g_ex = (double *)malloc(N*N*sizeof(double));
    this->g_el = (double *)malloc(N*N*sizeof(double));
    
    for(i=0;i<N*N;i++)
    {
        this->g_in[i] = g_in[i];
        this->g_ex[i] = g_ex[i];
        this->g_el[i] = g_el[i];
    }
    

    this->inhActivation = inhActivation;
    this->excActivation = excActivation;

    
    this->firstState = (int *)malloc(N*sizeof(int));
    
    index = 0;
    for(i=0;i<N;i++)
    {
        this->firstState[i] = index;
        index += (this->neuroni[i])->getnx();
    }
    
}


void CPG::getXdot(double t, double *x, double *xdot,double Iext)
{
    int N = this->N;
    int i,j,matrInd;
    neuron_model **neuroni = this->neuroni;
    
    double *g_in = this->g_in;
    double *g_ex = this->g_ex;
    double *g_el = this->g_el;
    double EsynIn = this->EsynIn;
    double EsynEx = this->EsynEx;
    
    synapse_model **inhActivation = this->inhActivation;
    synapse_model **excActivation = this->excActivation;
    
    int *Vindex = this->firstState;
    
    double Isyn,Vi,Vj;
    
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        Vi = x[Vindex[i]];
        
        for(j=0;j<N;j++)
        {
            matrInd = i*N+j;
            
            Vj = x[Vindex[j]];
            if(g_in[matrInd] != 0)
                Isyn += g_in[matrInd]*(EsynIn-Vi)*(inhActivation[matrInd]->getActivation(Vj));
            
            if(g_ex[matrInd] != 0)   
                Isyn += g_ex[matrInd]*(EsynEx-Vi)*(excActivation[matrInd]->getActivation(Vj));
            
            if(g_el[matrInd] != 0)
                Isyn += g_el[matrInd]*(Vj-Vi);
            
        }
        
        neuroni[i]->getXdot(t,x+Vindex[i],xdot+Vindex[i],Isyn);
    }
    
}


bool CPG::getResetConditions(double *x)
{
    int N = this->N;
    int i;
    neuron_model **neuroni = this->neuroni;
    int *Vindex = this->firstState;
    
    for(i=0;i<N;i++)
    {
        if((neuroni[i])->getResetConditions(x+Vindex[i]))
        {
            return true;
            break;
        }
    }
    
    return false;
}

void CPG::resetStates(double *x)
{
    int N = this->N;
    int i;
    neuron_model **neuroni = this->neuroni;
    int *Vindex = this->firstState;
    
    for(i=0;i<N;i++)
    {
        if((neuroni[i])->getResetConditions(x+Vindex[i]))
        {
            (neuroni[i])->resetStates(x+Vindex[i]);
        }
    }
    return;
}


CPG::~CPG() 
{
    free(this->g_in);
    free(this->g_ex);
    free(this->g_el);
 
}


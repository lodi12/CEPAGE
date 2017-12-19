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
    this->Ninh = 0;
    this->Nexc = 0;
    
    this->neuroni = (neuron_model **)malloc(0);

    this->g_el = (double *)malloc(0);
    
    this->inhSyn = new SynStruct*[0];
    this->excSyn = new SynStruct*[0];
    
    this->EsynIn = 0;
    this->EsynEx = 0;

    
    this->firstState = (int *)malloc(0);
}

CPG::CPG(int N,neuron_model **neuroni, double *g_el,double EsynIn, double EsynEx, t_SynStruct **inhSyn, t_SynStruct **excSyn )
{
    int i,index;
    
    this->N = N;
    this->Ninh = sizeof(inhSyn)/sizeof(t_SynStruct *);
    this->Nexc = sizeof(excSyn)/sizeof(t_SynStruct *);
    
    
    this->neuroni = neuroni;
    
    this->EsynIn = EsynIn;
    this->EsynEx = EsynEx;

    
    this->g_el = (double *)malloc(N*N*sizeof(double));
    
    for(i=0;i<N*N;i++)
    {
        this->g_el[i] = g_el[i];
    }
    

    this->inhSyn = inhSyn;
    this->excSyn = excSyn;

    
    this->firstState = (int *)malloc((N+Ninh+Nexc)*sizeof(int));
    
    index = 0;
    for(i=0;i<N;i++)
    {
        this->firstState[i] = index;
        index += (this->neuroni[i])->getnx();
    }
    
    for(i=0;i<Ninh;i++)
    {
        this->firstState[N+i] = index;
        index += ((this->inhSyn[i])->activationFunction)->getnx();
    }
    
    for(i=0;i<Nexc;i++)
    {
        this->firstState[N+N*N+i] = index;
        index += ((this->excSyn[i])->activationFunction)->getnx();
    }
    
    
}


void CPG::getXdot(double t, double *x, double *xdot,double *Iext)
{
    int i,j,matrInd;
    int inhIndex, excIndex;
    double Isyn,Vi,Vj;
    
    int N = this->N;
    int Ninh = this->Ninh;
    int Nexc = this->Nexc;
    
    t_SynStruct *tmp;

    
    neuron_model **neuroni = this->neuroni;
    
    double *g_el = this->g_el;
    double EsynIn = this->EsynIn;
    double EsynEx = this->EsynEx;
    
    
    t_SynStruct **inhSyn = this->inhSyn;
    t_SynStruct **excSyn = this->excSyn;
    
    int *Vindex = this->firstState;
    

    
    /* Compute neurons differentials */
    
    inhIndex = 0;
    excIndex = 0;
    
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        Vi = x[Vindex[i]];
        
        if(inhIndex < Ninh)
        {
            while(inhSyn[inhIndex]->i == i)
            {
                Vj = x[Vindex[inhSyn[inhIndex]->j]];
                Isyn += (inhSyn[inhIndex]->g)*(EsynIn-Vi)*((inhSyn[inhIndex]->activationFunction)->getActivation(Vj));
                inhIndex++;
            }
        }
        
        if(excIndex < Nexc)
        {
            while(excSyn[excIndex]->i == i)
            {
                Vj = x[Vindex[excSyn[excIndex]->j]];
                Isyn += (excSyn[excIndex]->g)*(EsynIn-Vi)*((excSyn[excIndex]->activationFunction)->getActivation(Vj));
                excIndex++;
            }
        }
        
        for(j=0;j<N;j++)
        {
            matrInd = i*N+j;
            
            
           /* if(g_in[matrInd] != 0)
                Isyn += g_in[matrInd]*(EsynIn-Vi)*(inhActivation[matrInd]->getActivation(Vj));
            
            if(g_ex[matrInd] != 0)   
                Isyn += g_ex[matrInd]*(EsynEx-Vi)*(excActivation[matrInd]->getActivation(Vj));
            */
            if(g_el[matrInd] != 0)
            {
                Vj = x[Vindex[j]];
                Isyn += g_el[matrInd]*(Vj-Vi);
            }
            
        }
        

        neuroni[i]->getXdot(t,x+Vindex[i],xdot+Vindex[i],&Isyn);
    }
    
    
    /* Compute synapses */
    for(i=0;i<Ninh;i++)
    {
        (inhSyn[i]->activationFunction)->getXdot(t,x+Vindex[N+i],xdot+Vindex[N+i]);
    }
   
    for(i=0;i<Nexc;i++)
    {
        (excSyn[i]->activationFunction)->getXdot(t,x+Vindex[N+Ninh+i],xdot+Vindex[N+Ninh+i]);
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
    
    for(i=0;i<Ninh;i++)
    {
        if((inhSyn[i]->activationFunction)->getResetConditions(x+Vindex[N+i]))
        {
            return true;
            break;
        }
    }
    for(i=0;i<Nexc;i++)
    {
        if ((excSyn[i]->activationFunction)->getResetConditions(x+Vindex[N+Ninh+i]))
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
    
    
    for(i=0;i<Ninh;i++)
    {
        if((inhSyn[i]->activationFunction)->getResetConditions(x+Vindex[N+i]))
        {
            (inhSyn[i]->activationFunction)->resetStates(x+Vindex[N+i]);
        }
    }
    for(i=0;i<Nexc;i++)
    {
        if ((excSyn[i]->activationFunction)->getResetConditions(x+Vindex[N+Ninh+i]))
        {
            (excSyn[i]->activationFunction)->resetStates(x+Vindex[N+Ninh+i]);
        }
    }
    
    return;
}


void CPG::getFirstIndex(int *firstIndex)
{
    int i;
    for(i=0;i<N+2*N*N;i++)
        firstIndex[i] = this->firstState[i];
}

CPG::~CPG() 
{
    
    free(this->g_el);
 
}


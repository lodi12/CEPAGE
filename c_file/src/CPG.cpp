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
#include <stdio.h>
#include "../inc/CPG.hpp"

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

CPG::CPG(int N,neuron_model **neuroni, double *g_in, double *g_ex, double *g_el,double EsynIn, double EsynEx, synapse_model **inhActivation, synapse_model **excActivation , double networkDelays[])
{
    int i,index,j,k;
    
    double *tmpVect;
    int Ndelaytmp;
    
    double diffMin;
    double iMin;
    
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
    
    
    
    
    this->Ndelay = sizeof(networkDelays)/sizeof(double);
    
        
    this->delays = (double *)malloc(Ndelay*sizeof(double));
    
    for(i=0;i<Ndelay;i++)
        this->delays[i] = networkDelays[i];
    
    this->firstState = (int *)malloc(N*sizeof(int)+2*N*N*sizeof(int));
    
    
    
    tmpVect = (double *)malloc(Ndelay*sizeof(double));
    
    this->neuronDelaysIndex = new vector<int>*[N];
    this->NneuronDelaysIndex = (int *)malloc(N*sizeof(int));
    
    index = 0;
    for(i=0;i<N;i++)
    {
        this->firstState[i] = index;
        index += (this->neuroni[i])->getnx();
        
        
        (this->neuroni[i])->getDelays(tmpVect);
        Ndelaytmp = (this->neuroni[i])->getNdelay();
        
        neuronDelaysIndex[i] = new vector<int>();
        for(j=0;j<Ndelaytmp;j++)
        {
            iMin = 0;
            diffMin = abs(networkDelays[0] - tmpVect[j]);
            for(k=1;k<Ndelay;k++)
            {
                if(diffMin > abs(networkDelays[k] - tmpVect[j]))
                {
                    diffMin = abs(networkDelays[k] - tmpVect[j]);
                    iMin = k;
                }
            }
                
            neuronDelaysIndex[i]->push_back(iMin);
        }
        
        
        NneuronDelaysIndex[i] = Ndelaytmp;
        
        
    }
    
    
    this->inhSynDelaysIndex = new vector<int>*[N*N];
    this->NinhSynDelaysIndex = (int *)malloc(N*N*sizeof(int));
    
    for(i=0;i<N*N;i++)
    {
        this->firstState[N+i] = index;
        index += (this->inhActivation[i])->getnx();
        
        
        (this->inhActivation[i])->getDelays(tmpVect);
        Ndelaytmp = (this->inhActivation[i])->getNdelay();
        
        inhSynDelaysIndex[i] = new vector<int>();
        for(j=0;j<Ndelaytmp;j++)
        {
            iMin = 0;
            diffMin = abs(networkDelays[0] - tmpVect[j]);
            for(k=1;k<Ndelay;k++)
            {
                if(diffMin > abs(networkDelays[k] - tmpVect[j]))
                {
                    diffMin = abs(networkDelays[k] - tmpVect[j]);
                    iMin = k;
                }
            }
                
            inhSynDelaysIndex[i]->push_back(iMin);

        }
        
        
        NinhSynDelaysIndex[i] = Ndelaytmp;
        
    }
    
    this->excSynDelaysIndex = new vector<int>*[N*N];
    this->NexcSynDelaysIndex = (int *)malloc(N*N*sizeof(int));
    
    for(i=0;i<N*N;i++)
    {
        this->firstState[N+N*N+i] = index;
        index += (this->excActivation[i])->getnx();
        
        
        (this->excActivation[i])->getDelays(tmpVect);
        Ndelaytmp = (this->excActivation[i])->getNdelay();
        
        excSynDelaysIndex[i] = new vector<int>();
        for(j=0;j<Ndelaytmp;j++)
        {
            iMin = 0;
            diffMin = abs(networkDelays[0] - tmpVect[j]);
            for(k=1;k<Ndelay;k++)
            {
                if(diffMin > abs(networkDelays[k] - tmpVect[j]))
                {
                    diffMin = abs(networkDelays[k] - tmpVect[j]);
                    iMin = k;
                }
            }
                
            excSynDelaysIndex[i]->push_back(iMin);
        }
        
        NexcSynDelaysIndex[i] = Ndelaytmp;
        
    }
    
    
    free(tmpVect);
    
}


void CPG::getXdot(double t, double *x, double *xdot,double Iext,double **Xprec)
{   
    int N = this->N;
    int i,j,k,matrInd;
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
    
    
    int maxNdelay = this->Ndelay;
    int *delays_;
    
    
    double *Vjprec;
    double **Xprec_;

    Xprec_ = (double**)malloc(maxNdelay*sizeof(double *));
    
    Vjprec = (double *)malloc(maxNdelay*sizeof(double));
    
    
    /* Compute neurons differentials */
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        Vi = x[Vindex[i]];
        
        for(j=0;j<N;j++)
        {
            matrInd = i*N+j;
            
            Vj = x[Vindex[j]];
            if(g_in[matrInd] != 0)
            {
                
                for(k=0;k<NinhSynDelaysIndex[matrInd];k++)
                {
                    Xprec_[k] = Xprec[inhSynDelaysIndex[matrInd]->at(k)]+Vindex[N+matrInd];
                    Vjprec[k] = Xprec[inhSynDelaysIndex[matrInd]->at(k)][Vindex[j]];
                }
                
                Isyn += g_in[matrInd]*(EsynIn-Vi)*(inhActivation[matrInd]->getActivation(x+Vindex[N+matrInd],Vj,Vjprec));
                inhActivation[matrInd]->getXdot(t,x+Vindex[N+matrInd],xdot+Vindex[N+matrInd],Vj,Xprec_);
            }

            if(g_ex[matrInd] != 0)
            {
                
                for(k=0;k<NexcSynDelaysIndex[matrInd];k++)
                {
                    Xprec_[k] = Xprec[excSynDelaysIndex[matrInd]->at(k)]+Vindex[N+N*N+matrInd];
                    Vjprec[k] = Xprec[excSynDelaysIndex[matrInd]->at(k)][Vindex[j]];
                }
                
                excActivation[matrInd]->getXdot(t,x+Vindex[N+N*N+matrInd],xdot+Vindex[N+N*N+matrInd],Vj,Xprec_);
                Isyn += g_ex[matrInd]*(EsynEx-Vi)*(excActivation[matrInd]->getActivation(x+Vindex[N+N*N+matrInd],Vj,Vjprec));
            }
            
            if(g_el[matrInd] != 0)
                Isyn += g_el[matrInd]*(Vj-Vi);
            
        }
        
     
        for(k=0;k<NneuronDelaysIndex[i];k++)
        {
            Xprec_[k] = Xprec[neuronDelaysIndex[i]->at(k)]+Vindex[i];
        }

        neuroni[i]->getXdot(t,x+Vindex[i],xdot+Vindex[i],Isyn,Xprec_);
    }
    
    free(Xprec_);
    free(Vjprec);
    
    
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
    
    /* Compute neurons differentials */
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        Vi = x[Vindex[i]];
        
        for(j=0;j<N;j++)
        {
            matrInd = i*N+j;
            
            Vj = x[Vindex[j]];
            if(g_in[matrInd] != 0)
            {
                Isyn += g_in[matrInd]*(EsynIn-Vi)*(inhActivation[matrInd]->getActivation(x+Vindex[N+matrInd],Vj));
                inhActivation[matrInd]->getXdot(t,x+Vindex[N+matrInd],xdot+Vindex[N+matrInd],Vj);
            }
            
            if(g_ex[matrInd] != 0)
            {
                excActivation[matrInd]->getXdot(t,x+Vindex[N+N*N+matrInd],xdot+Vindex[N+N*N+matrInd],Vj);
                Isyn += g_ex[matrInd]*(EsynEx-Vi)*(excActivation[matrInd]->getActivation(x+Vindex[N+N*N+matrInd],Vj));
            }
            
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
    double Vj;
    for(i=0;i<N;i++)
    {
        Vj = x[Vindex[i%N]];
        if((neuroni[i])->getResetConditions(x+Vindex[i]))
        {
            return true;
            break;
        }
    }
    
    for(i=0;i<N*N;i++)
    {
        Vj = x[Vindex[i%N]];
        if(((inhActivation[i])->getResetConditions(x+Vindex[N+i],Vj)) || ((excActivation[i])->getResetConditions(x+Vindex[N+N*N+i],Vj)))
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
    double Vj;
    
    for(i=0;i<N;i++)
    {
        if((neuroni[i])->getResetConditions(x+Vindex[i]))
        {
            (neuroni[i])->resetStates(x+Vindex[i]);
        }
    }
    
    
    for(i=0;i<N*N;i++)
    {
        
        Vj = x[Vindex[i%N]];
        if((inhActivation[i])->getResetConditions(x+Vindex[N+i],Vj))
        {
            (inhActivation[i])->resetStates(x+Vindex[N+i],Vj);
        }
        
        if((excActivation[i])->getResetConditions(x+Vindex[N+N*N+i],Vj))
        {
            (excActivation[i])->resetStates(x+Vindex[N+N*N+i],Vj);
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
    int i;
    
    free(this->g_in);
    free(this->g_ex);
    free(this->g_el);
    
    free(this->delays);
    
    for(i=0;i<N;i++)
    {
        delete(neuronDelaysIndex[i]);
    }
    delete[] (neuronDelaysIndex);
    free(this->NneuronDelaysIndex);
    
    
    
    for(i=0;i<N*N;i++)
    {
        delete(inhSynDelaysIndex[i]);
    }
    delete[] (inhSynDelaysIndex);
    free(this->NinhSynDelaysIndex);
    
    
    for(i=0;i<N*N;i++)
    {
        delete(excSynDelaysIndex[i]);
    }
    delete[] (excSynDelaysIndex);
    free(this->NexcSynDelaysIndex);
    
}


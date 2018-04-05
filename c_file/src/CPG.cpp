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

#include "../inc/CPG.hpp"

CPG::CPG()
{
    this->N = 0;
    this->neuroni = (neuron_model **)malloc(0);
    
    
    this->firstState = (int *)malloc(0);
    
    
    NneuronDelaysIndex = (int *)malloc(0);
    NinhSynDelaysIndex = (int *)malloc(0);
    NexcSynDelaysIndex = (int *)malloc(0);
    
    this->Ninh = 0;
    this->Nexc = 0;
    this->Nel = 0;
    
    firstState = (int *)malloc(0);
    
}

CPG::CPG(int N,neuron_model **neuroni, int Ninh, int Nexc, int Nel, synStruct_t **inhSyn, synStruct_t **excSyn , synStruct_t **elSyn, int Ndelay, double networkDelays[])
{
    int i,index,j,k;
    
    double *tmpVect;
    int Ndelaytmp;
    
    double diffMin;
    double iMin;
    
    synapse_model *activation;
    
    this->N = N;
    
    this->neuroni = new neuron_model*[N];
    
    for(i=0;i<N;i++)
        this->neuroni[i] = neuroni[i];
    
    this->Ninh = Ninh;
    this->Nexc = Nexc;
    this->Nel = Nel;
    
    
    this->inhSyn = new synStruct*[Ninh];
    this->excSyn = new synStruct*[Nexc];
    this->elSyn = new synStruct*[Nel];
    
    for(i=0;i<Ninh;i++)
        this->inhSyn[i] = inhSyn[i];
    for(i=0;i<Nexc;i++)
        this->excSyn[i] = excSyn[i];
    for(i=0;i<Nel;i++)
        this->elSyn[i] = elSyn[i];
    
    
    this->Ndelay = Ndelay;//((int)sizeof(networkDelays))/sizeof(double);
    
    this->delays = networkDelays;
    
    this->firstState = (int *)malloc(N*sizeof(int)+Ninh+Nexc*sizeof(int));
    
    
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
 
    
    this->inhSynDelaysIndex = new vector<int>*[Ninh];
    this->NinhSynDelaysIndex = (int *)malloc(Ninh*sizeof(int));
    
    for(i=0;i<Ninh;i++)
    {
        this->firstState[N+i] = index;
        
        activation = this->inhSyn[i]->activation;
        
        index += activation->getnx();
        
        activation->getDelays(tmpVect);
        Ndelaytmp = activation->getNdelay();
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
    

    
    this->excSynDelaysIndex = new vector<int>*[Nexc];
    this->NexcSynDelaysIndex = (int *)malloc(Nexc*sizeof(int));
    
    for(i=0;i<Nexc;i++)
    {
        this->firstState[N+Ninh+i] = index;
        
        activation = this->excSyn[i]->activation;
        
        index += activation->getnx();
        
        
        activation->getDelays(tmpVect);
        Ndelaytmp = activation->getNdelay();
        
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
    
    int Ninh = this->Ninh;
    int Nexc = this->Nexc;
    int Nel = this->Nel;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->excSyn;
    
    int *Vindex = this->firstState;
    
    double Isyn,Vi,Vj,g,Esyn;
    
    synapse_model *a;
    
    int maxNdelay = this->Ndelay;
    int *delays_;
    
    
    double *Vjprec;
    double **Xprec_;
    
    
    int iInh,iExc,iEl;
    
    
    Xprec_ = (double**)malloc(maxNdelay*sizeof(double *));
    
    Vjprec = (double *)malloc(maxNdelay*sizeof(double));
    
    iInh = 0;
    iExc = 0;
    iEl = 0;
    
    /* Compute neurons differentials */
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        Vi = x[Vindex[i]];
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {
                j = inhSyn[iInh]->j;
                g = inhSyn[iInh]->g;
                Esyn = inhSyn[iInh]->Esyn;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                
                
                for(k=0;k<NinhSynDelaysIndex[matrInd];k++)
                {
                    Xprec_[k] = Xprec[inhSynDelaysIndex[iInh]->at(k)]+Vindex[N+iInh];
                    Vjprec[k] = Xprec[inhSynDelaysIndex[iInh]->at(k)][Vindex[j]];
                }
                
                Isyn += g*(Esyn-Vi)*(a->getActivation(x+Vindex[N+iInh],Vj,Vjprec));
                a->getXdot(t,x+Vindex[N+matrInd],xdot+Vindex[N+iInh],Vj,Xprec_);
                
                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                g = excSyn[iExc]->g;
                Esyn = excSyn[iExc]->Esyn;
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                
                for(k=0;k<NexcSynDelaysIndex[matrInd];k++)
                {
                    Xprec_[k] = Xprec[excSynDelaysIndex[iExc]->at(k)]+Vindex[N+Ninh+iExc];
                    Vjprec[k] = Xprec[excSynDelaysIndex[iExc]->at(k)][Vindex[j]];
                }
                
                a->getXdot(t,x+Vindex[N+N*N+matrInd],xdot+Vindex[N+Ninh+iExc],Vj,Xprec_);
                Isyn += g*(Esyn-Vi)*(a->getActivation(x+Vindex[N+Ninh+iExc],Vj,Vjprec));
                
                iExc++;
            }
        }
        
        while(iEl < Nel)
        {
            if (elSyn[iEl]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                g = excSyn[iExc]->g;
                
                Vj = x[Vindex[j]];
                Isyn += g*(Vj-Vi);
                
                iEl++;
            }
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
    
    int Ninh = this->Ninh;
    int Nexc = this->Nexc;
    int Nel = this->Nel;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->excSyn;
    
    int *Vindex = this->firstState;
    
    double Isyn,Vi,Vj,g,Esyn;
    
    synapse_model *a;
    
    int iInh,iExc,iEl;
    
    /* Compute neurons differentials*/
    
    iInh = 0;
    iExc = 0;
    iEl = 0;
    
    for(i=0;i<N;i++)
    {
        Isyn = 0;
        Vi = x[Vindex[i]];
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {
                j = inhSyn[iInh]->j;
                g = inhSyn[iInh]->g;
                Esyn = inhSyn[iInh]->Esyn;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                
                Isyn += g*(Esyn-Vi)*(a->getActivation(x+Vindex[N+iInh],Vj));
                a->getXdot(t,x+Vindex[N+matrInd],xdot+Vindex[N+iInh],Vj);
                
                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                g = excSyn[iExc]->g;
                Esyn = excSyn[iExc]->Esyn;
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                
                a->getXdot(t,x+Vindex[N+N*N+matrInd],xdot+Vindex[N+Ninh+iExc],Vj);
                Isyn += g*(Esyn-Vi)*(a->getActivation(x+Vindex[N+Ninh+iExc],Vj));
                
                iExc++;
            }
        }
        
        while(iEl < Nel)
        {
            if (elSyn[iEl]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                g = excSyn[iExc]->g;
                
                Vj = x[Vindex[j]];
                Isyn += g*(Vj-Vi);
                
                iEl++;
            }
        }
        neuroni[i]->getXdot(t,x+Vindex[i],xdot+Vindex[i],Isyn);
    }
}

bool CPG::getResetConditions(double *x)
{
    int N = this->N;
    int i,j;
    neuron_model **neuroni = this->neuroni;
    int *Vindex = this->firstState;
    double Vj;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->excSyn;
    
    synapse_model *a;
    
    int iInh,iExc;
    
    iInh = 0;
    iExc = 0;
    
    
    
    for(i=0;i<N;i++)
    {
        
        if((neuroni[i])->getResetConditions(x+Vindex[i]))
        {
            return true;
            break;
        }
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {
                j = inhSyn[iInh]->j;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                
                if (a->getResetConditions(x+Vindex[N+iInh],Vj))
                    return true;
                
                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                if (a->getResetConditions(x+Vindex[N+Ninh+iExc],Vj))
                    return true;
                
                iExc++;
            }
        }
    }
    return false;
}


void CPG::resetStates(double *x)
{
    int N = this->N;
    int i,j;
    neuron_model **neuroni = this->neuroni;
    int *Vindex = this->firstState;
    double Vj;
    
    synStruct **inhSyn = this->inhSyn;
    synStruct **excSyn = this->excSyn;
    synStruct **elSyn = this->excSyn;
    
    synapse_model *a;
    
    int iInh,iExc;
    
    iInh = 0;
    iExc = 0;
    
    for(i=0;i<N;i++)
    {
        if((neuroni[i])->getResetConditions(x+Vindex[i]))
        {
            (neuroni[i])->resetStates(x+Vindex[i]);
        }
        
        while(iInh < Ninh)
        {
            if (inhSyn[iInh]->i > i)
            {
                break;
            }
            else
            {
                j = inhSyn[iInh]->j;
                a = inhSyn[iInh]->activation;
                Vj = x[Vindex[j]];
                if (a->getResetConditions(x+Vindex[N+iInh],Vj))
                    a->resetStates(x+Vindex[N+iInh],Vj);
                
                iInh++;
            }
        }
        
        while(iExc < Nexc)
        {
            if (excSyn[iExc]->i > i)
            {
                break;
            }
            else
            {
                j = excSyn[iExc]->j;
                a = excSyn[iExc]->activation;
                Vj = x[Vindex[j]];
                
                if (a->getResetConditions(x+Vindex[N+Ninh+iExc],Vj))
                    a->resetStates(x+Vindex[N+Ninh+iExc],Vj);
                
                iExc++;
            }
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
    
    free(this->firstState);
        
    for(i=0;i<N;i++)
        delete(this->neuroni[i]);
    delete[](neuroni);
    
    for(i=0;i<Ninh;i++)
        delete(this->inhSyn[i]);
    delete[](inhSyn);
    
    for(i=0;i<Nexc;i++)
        delete(this->excSyn[i]);
    delete[](excSyn);
    
    for(i=0;i<Nel;i++)
        delete(this->elSyn[i]);
    delete[](elSyn);
}


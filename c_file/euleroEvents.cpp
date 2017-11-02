#include "mex.h"
#include "math.h"
#include "vectorField.h"
#include <stdio.h>
/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    mwSize nStep;   /* Number of step */
    double dt;      /* Integration step */
    mwSize Nstati;       /* Number of state variables*/
    mwSize i,j;
    
    double *dx;

    double *oldState;
    double *currentState;
    double currentT;
    
    double **eventMatrix;
    
    double *phiOut;
    
    double Vth;
    
    double period;
    
    double stopThreshold;
    
    double diff,tmp;
    
    mwSize *nEvent;

    mwSize minEventNumber;
    
    mwSize nx_;
    
    mwSize N;
    
        dynSys *vectorField;

        
    /* check for proper number of arguments */
    if(nrhs!=7) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 7 Input required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 1 Output required.");
    }
    
     initVectorField(&vectorField);
  
    
    Nstati =  mxGetScalar(prhs[0]);   
    nStep = mxGetScalar(prhs[1]);    
    dt = mxGetScalar(prhs[2]);
    oldState = mxGetPr(prhs[3]);
    Vth = mxGetScalar(prhs[4]);
    stopThreshold = mxGetScalar(prhs[5]);
    nx_ =  mxGetScalar(prhs[6]);
    
    N = Nstati/nx_;
    
    currentState = (double *)mxMalloc(Nstati*sizeof(double));
    dx = (double *)mxMalloc(Nstati*sizeof(double));
    nEvent = (mwSize *)mxMalloc(N*sizeof(mwSize));
    eventMatrix = (double **)mxMalloc(N*sizeof(double *));   
   
    for(j=0;j<N;j++)
    {
        eventMatrix[j] = (double *)mxMalloc(0);
        nEvent[j] = 0;
    }
    
    for(j=0;j<Nstati;j++)
    {
        currentState[j] = 0;
    }

    currentT = 0;
    
    for(i=1;i<nStep;i++)
    {
        
        
        
        vectorField->getXdot(0,oldState,dx,0);
        currentT += dt;

        for(j=0;j<Nstati;j++)
            currentState[j] = oldState[j]+dt*dx[j];
        
        
        if (vectorField->getResetConditions(currentState))
        {
            vectorField->resetStates(currentState);
        }

        
        for(j=0;j<Nstati;j+=nx_)
        {
            if((oldState[j] < Vth) && (currentState[j] > Vth))
            {
                nEvent[j/nx_]++;
                eventMatrix[j/nx_] = (double *)mxRealloc(eventMatrix[j/nx_],nEvent[j/nx_]*sizeof(double));
                eventMatrix[j/nx_][nEvent[j/nx_]-1] = currentT;
            }
        }

        for(j=0;j<Nstati;j++)
        {
            oldState[j] = currentState[j];
        }
        /*diff = 0;

        for(j=0;j<Nstati;j++)
        {
            tmp = oldState[j] - currentState[j];
            if (tmp < 0)
                tmp *=-1;
            diff+= tmp;
            oldState[j] = currentState[j];
        }
        
        

        if(diff < stopThreshold)
            break;*/
    }
    
    
    
     minEventNumber = nEvent[0];

     for(i=1;i<N;i++)
         if(minEventNumber > nEvent[i])
             minEventNumber = nEvent[i];
    
    plhs[0] = mxCreateDoubleMatrix(minEventNumber,N-1,mxREAL);
    phiOut = mxGetPr(plhs[0]);
    if(minEventNumber > 1)
    {
        period = eventMatrix[0][minEventNumber-1] - eventMatrix[0][minEventNumber-2];

        for(i=0;i<N-1;i++)
        {
            for(j=0;j<minEventNumber;j++)
            {
               phiOut[i*minEventNumber+j] =  ((eventMatrix[i+1][j]-eventMatrix[0][j])/period);
            }
        }
    }
    
    mxFree(dx);
    mxFree(currentState);
    mxFree(nEvent);
    for(j=0;j<N;j++)
    {
        mxFree(eventMatrix[j]);
    }
    mxFree(eventMatrix);
    
}

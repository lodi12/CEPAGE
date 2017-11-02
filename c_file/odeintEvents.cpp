#include "mex.h"
#include "math.h"
#include "vectorField.h"
#include <boost/serialization/array_wrapper.hpp>
#include <boost/numeric/odeint.hpp>
#include <vector>

typedef std::vector< double > state_type;

mwSize Nstati;       /* Number of state variables*/
double *xC;
double *dxC;

double Vth;


mwSize *nEvent;

mwSize minEventNumber;

mwSize nx_;

state_type oldState;

double oldT;

dynSys *vectorField;

/* The rhs of x' = f(x) */
void cppVectorialField( const state_type &x , state_type &dxdt , const double  t  )
{
    mwSize i;
    for(i=0;i<Nstati;i++)
        xC[i] = x[i];
    
    
    vectorField->getXdot(t,xC,dxC,0);
    
    for(i=0;i<Nstati;i++)
        dxdt[i] = dxC[i];
}

struct push_back_events
{
   double **eventMatrix;
    
    push_back_events(double **eventMatrix_) : eventMatrix(eventMatrix_) { }
    
    void operator()( state_type &x , double t )
    {
        mwSize i,j;
        for(i=0;i<Nstati;i++)
            xC[i] = x[i];


        if (vectorField->getResetConditions(xC))
        {
            vectorField->resetStates(xC);
        }

        for(i=0;i<Nstati;i++)
            x[i] = xC[i];
        
        for(j=0;j<Nstati;j+=nx_)
        {
            if((oldState[j] < Vth) && (x[j] > Vth))
            {
                nEvent[j/nx_]++;
                eventMatrix[j/nx_] = (double *)mxRealloc(eventMatrix[j/nx_],nEvent[j/nx_]*sizeof(double));
                eventMatrix[j/nx_][nEvent[j/nx_]-1] = (Vth-oldState[j])/(x[j]-oldState[j])*(t-oldT)+oldT;
            }
        }
        copy(x.begin(), x.end(), oldState.begin());
        oldT = t;
    }
};

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    using namespace std;
    using namespace boost::numeric::odeint;
    
    mwSize nStep;   /* Number of step */
    double dt;      /* Integration step */
    mwSize i,j;
    double *ptr;
    double Tfinale;
    
    
    double *x0;     /* Inital condition */
    double *xOutMatrix; /* Out vector */
    double *dx;
    
     double *phiOut;   
    double period;

    /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 4 Input required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:sim:nrhs","Error! 1 Output required.");
    }
    
    Nstati =  mxGetScalar(prhs[0]);
    Tfinale = mxGetScalar(prhs[1]);
    dt = mxGetScalar(prhs[2]);
    x0 = mxGetPr(prhs[3]);
    Vth = mxGetScalar(prhs[4]);
    nx_ =  mxGetScalar(prhs[5]);
    state_type xinit(Nstati);
    state_type x(Nstati);
    
    double **eventMatrix;
    
    xC = (double *)(mxMalloc(Nstati*sizeof(double)));
    dxC = (double *)(mxMalloc(Nstati*sizeof(double)));
    
    oldT = 0;
    nEvent = (mwSize *)mxMalloc(Nstati/nx_*sizeof(mwSize));
    eventMatrix = (double **)mxMalloc(Nstati/nx_*sizeof(double *));
    
    initVectorField(&vectorField);
    
    for(i=0;i<Nstati;i++)
    {
        xinit[i] = x0[i];
        oldState.push_back(x0[i]);
    }
    
    for(j=0;j<Nstati/nx_;j++)
    {
        eventMatrix[j] = (double *)mxMalloc(0);
        nEvent[j] = 0;
    }
    

    typedef runge_kutta_cash_karp54< state_type > error_stepper_type;
    
    
    typedef controlled_runge_kutta< error_stepper_type > controlled_stepper_type;
    controlled_stepper_type controlled_stepper;
    
    size_t steps = integrate_adaptive( controlled_stepper , cppVectorialField , xinit , 0.0 , Tfinale, dt ,
            push_back_events(eventMatrix));
    minEventNumber = nEvent[0];

     for(i=1;i<(Nstati/nx_);i++)
         if(minEventNumber > nEvent[i])
             minEventNumber = nEvent[i];

    mxFree(xC);
    mxFree(dxC);
    
    plhs[0] = mxCreateDoubleMatrix(minEventNumber,(Nstati/nx_)-1,mxREAL);
    phiOut = mxGetPr(plhs[0]);

    if(minEventNumber > 1)
    {
        period = eventMatrix[0][minEventNumber-1] - eventMatrix[0][minEventNumber-2];

        for(i=0;i<(Nstati/nx_)-1;i++)
        {
            for(j=0;j<minEventNumber;j++)
            {
               phiOut[i*minEventNumber+j] =  ((eventMatrix[i+1][j]-eventMatrix[0][j])/period);
            }
        }
    }
    
    for(j=0;j<(Nstati/nx_);j++)
    {
        mxFree(eventMatrix[j]);
    }
    mxFree(eventMatrix);
    
}

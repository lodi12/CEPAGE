/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FN_relaxation_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 9:48 AM
 */

#ifndef FN_RELAXATION_MODEL_H
#define FN_RELAXATION_MODEL_H



#include "neuron_model.hpp"
#include <math.h>
class FN_relaxation_model : public neuron_model {
    
private:
    double eps;
    double I;
    
public:
    FN_relaxation_model();
    FN_relaxation_model(double I,double eps);
    
    void getXdot(double t, double *x, double *xdot,double Iext);
    
    virtual ~FN_relaxation_model();

};

#endif /* FN_RELAXATION_MODEL_H */


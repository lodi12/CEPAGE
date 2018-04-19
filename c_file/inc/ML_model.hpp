/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   ML_model.h
 * Author: teo
 *
 * Created on July 12, 2017, 10:21 AM
 */

#ifndef ML_MODEL_H
#define ML_MODEL_H

#include "neuron_model.hpp"
#include <math.h>
class ML_model : public neuron_model {
    
private:
    double CM;
    double gCa;
    double V3;
    double V4;
    double phi;
    double gl;
    double Vl;
    double I;
    
public:
    ML_model();
    ML_model(const ML_model &m);
    ML_model(double CM,double gCa,double V3,double V4,double phi,double gl,double Vl,double I);
    
    ML_model *clone() const;
    
    void getXdot(double t, double *x, double *xdot,double Iext);
    
    virtual ~ML_model();

};


#endif /* ML_MODEL_H */


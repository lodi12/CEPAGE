/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   HR_model.h
 * Author: picio
 *
 * Created on 14 maggio 2017, 19.31
 */

#ifndef HR_MODEL_H
#define HR_MODEL_H

#include "neuron_model.hpp"

class HR_model : public neuron_model {
    
private:
    double b;
    double I;
    double mu;
    double s;
    double x_rest;
    
public:
    HR_model();
    HR_model(double b, double I, double mu, double s, double x_rest);
    
    void getXdot(double t, double *x, double *xdot,double Iext);
    
    virtual ~HR_model();

};

#endif /* HR_MODEL_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   FN_relaxation_model.cpp
 * Author: teo
 * 
 * Created on July 12, 2017, 9:48 AM
 */

#include "../inc/FN_relaxation_model.hpp"


FN_relaxation_model::FN_relaxation_model() 
{
    this->nx = 2;
    this->I = 0;
    this->eps = 0;
}

FN_relaxation_model::FN_relaxation_model(double I,double eps) 
{
    this->nx = 2;
    this->I = I;
    this->eps = eps;
}

FN_relaxation_model::FN_relaxation_model(const FN_relaxation_model &f)
{
    this->nx = f.nx;
    this->I = f.I;
    this->eps = f.eps;
}

FN_relaxation_model *FN_relaxation_model::clone() const
{
    return new FN_relaxation_model(*this);
}

void FN_relaxation_model::getXdot(double t, double *x, double *xdot,double Iext)
{
    double eps = this->eps;
    double I = this->I;
    
    xdot[0] = x[0]-x[0]*x[0]*x[0]+I-x[1]+Iext;
    xdot[1] = eps*(1/(1+exp(-10*x[0]))-x[1]);
    
    
}

FN_relaxation_model::~FN_relaxation_model() {
}


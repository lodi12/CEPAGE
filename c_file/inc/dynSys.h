/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   dynSys.h
 * Author: teo
 *
 * Created on July 12, 2017, 11:27 AM
 */

#ifndef DYNSYS_H
#define DYNSYS_H

#include "stdlib.h"

class dynSys 
{
public:
    virtual void getXdot(double t, double *x, double *xdot,double Iext) = 0;
    virtual bool getResetConditions(double *x) = 0;
    virtual void resetStates(double *x) = 0;

    virtual void getFirstIndex(int *firstIndex);
    

};

#endif /* DYNSYS_H */


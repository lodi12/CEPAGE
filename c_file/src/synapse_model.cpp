/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   synapse_model.cpp
 * Author: picio
 * 
 * Created on 14 maggio 2017, 19.28
 */

#include "../inc/synapse_model.h"


int synapse_model::getnx() 
{
    return this->nx;
}

bool synapse_model::getResetConditions(double *x)
{
    return false;
}

void synapse_model::resetStates(double *x)
{
    return;
}
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   neuron_model.cpp
 * Author: picio
 * 
 * Created on 14 maggio 2017, 19.28
 */

#include "../inc/neuron_model.hpp"

int neuron_model::getnx() 
{
    return this->nx;
}

bool neuron_model::getResetConditions(double *x)
{
    return false;
}

void neuron_model::resetStates(double *x)
{
    return;
}



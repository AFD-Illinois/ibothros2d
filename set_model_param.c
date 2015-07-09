/*
 
 set_model_param.c
 
 
 Created by John McCann on 6/23/13.
 Copyright (c) 2013 Gammie Labs. All rights reserved.
 
 File Description:
 Reads in all the model parameters which have been so conveniently placed in
 one file.
 
 */

#include "decs.h"
#include "defs.h"

#define stringMAX 121

void set_model_param(void){
    
    FILE *fp;
    
    fp = fopen("model_param.dat","r");
    if(fp==NULL){
        printf("Failed to open model_param.dat\n");
        exit(0);
    }
    
    char buffer[stringMAX];
    fgets(buffer, stringMAX, fp); //Reads Instructions, does nothing
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &MBH);
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &M_unit);
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &Dsource);
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &TP_OVER_TE);
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &THETAE_MAX);
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &pin);
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &FRAC);
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &Gmin);
    fgets(buffer, stringMAX, fp);
    sscanf(buffer, "%*s %*s %lf", &Gmax);
    
    //printf("MBH = %.2e, M_unit = %.2e, Dsource = %.2e, \npin = %.2e, FRAC = %.2e, Gmin = %.2e, Gmax = %.2e \n", MBH, M_unit, Dsource, pin, FRAC, Gmin, Gmax);
    
}
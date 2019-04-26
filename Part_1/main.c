//
//  main.c
//  Assignment_6
//
//  Created by Jake Ursetta on 3/13/17.
//  Copyright © 2017 Jake Ursetta. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <libgen.h>
#include <sys/stat.h>
#include <sys/types.h>

#define pi acos(-1.0)

void deriv(double (*SCstate)[4], double (*M_State)[4], double (*E_state)[4],double G, double m_E, double m_M, double m_S){
    double dSE = pow(pow((*E_state)[0] - (*SCstate)[0],2) + pow((*E_state)[1] - (*SCstate)[1],2),.5);
    double dME = pow(pow((*E_state)[0] - (*M_State)[0],2) + pow((*E_state)[1] - (*M_State)[1],2),.5);
    double dSM = pow(pow((*SCstate)[0] - (*M_State)[0],2) + pow((*SCstate)[1] - (*M_State)[1],2),.5);
    
    double SC0Hold = (*SCstate)[2];
    double SC1Hold = (*SCstate)[3];
    double M0Hold  = (*M_State)[2];
    double M1Hold  = (*M_State)[3];
    double E0Hold  = (*E_state)[2];
    double E1Hold  = (*E_state)[3];
    
    (*SCstate)[2] = G*m_E*((*E_state)[0] - (*SCstate)[0])/pow(dSE,3) + G*m_M*((*M_State)[0] - (*SCstate)[0])/pow(dSM,3);
    (*SCstate)[3] = G*m_E*((*E_state)[1] - (*SCstate)[1])/pow(dSE,3) + G*m_M*((*M_State)[1] - (*SCstate)[1])/pow(dSM,3);
    (*M_State)[2] = G*m_E*((*E_state)[0] - (*M_State)[0])/pow(dME,3) + G*m_S*((*SCstate)[0] - (*M_State)[0])/pow(dSM,3);
    (*M_State)[3] = G*m_E*((*E_state)[1] - (*M_State)[1])/pow(dME,3) + G*m_S*((*SCstate)[1] - (*M_State)[1])/pow(dSM,3);
    (*E_state)[2] = G*m_M*((*M_State)[0] - (*E_state)[0])/pow(dME,3) + G*m_S*((*SCstate)[0] - (*E_state)[0])/pow(dSE,3);
    (*E_state)[3] = G*m_M*((*M_State)[1] - (*E_state)[1])/pow(dME,3) + G*m_S*((*SCstate)[1] - (*E_state)[1])/pow(dSE,3);
    
    (*SCstate)[0] = SC0Hold;
    (*SCstate)[1] = SC1Hold;
    (*M_State)[0] = M0Hold;
    (*M_State)[1] = M1Hold;
    (*E_state)[0] = E0Hold;
    (*E_state)[1] = E1Hold;
    
    return;
}



int rungeKutta(double h, double (*earthPosVel)[4], double (*spacePosVel)[4], double (*moonPosVel)[4], double (*varArray)[6], double *t, double *tmin, int plot, const char *saveName){

    //printf("%s\n",saveName);
    
    FILE *outputFile;
    outputFile = fopen(saveName,"w");
    
    double m_E  = (*varArray)[0];
    double m_M  = (*varArray)[1];
    double m_S  = (*varArray)[2];
    double r_E  = (*varArray)[3];
    double r_M  = (*varArray)[4];
    double G    = (*varArray)[5];
    
    double dSE;
    double dME;
    double dSM;
    
    double c1 = 1;
    double c2 = 1;
    double c3 = 1;
    
    double spaceStateTemp[4];
    double lunarStateTemp[4];
    double earthStateTemp[4];
    double K[12][5];
    
    for (int i = 0; i<8; i++) {
        K[i][0] = 0;
    }
    
    double hval[4];
    hval[0] = 0;
    hval[1] = h/2;
    hval[2] = h/2;
    hval[3] = h;
    
    
    
    //int counter = 0;
    
    while (c1 > 0 && c2 > 0 && c3 > 0 && *t <= *tmin ) {
        
        for (int i = 1; i <= 4; i++) {
            spaceStateTemp[0] = (*spacePosVel)[0] + hval[i-1]*K[0][i-1];
            spaceStateTemp[1] = (*spacePosVel)[1] + hval[i-1]*K[1][i-1];
            spaceStateTemp[2] = (*spacePosVel)[2] + hval[i-1]*K[2][i-1];
            spaceStateTemp[3] = (*spacePosVel)[3] + hval[i-1]*K[3][i-1];
            lunarStateTemp[0] = (*moonPosVel)[0]  + hval[i-1]*K[4][i-1];
            lunarStateTemp[1] = (*moonPosVel)[1]  + hval[i-1]*K[5][i-1];
            lunarStateTemp[2] = (*moonPosVel)[2]  + hval[i-1]*K[6][i-1];
            lunarStateTemp[3] = (*moonPosVel)[3]  + hval[i-1]*K[7][i-1];
            earthStateTemp[0] = (*earthPosVel)[0] + hval[i-1]*K[8][i-1];
            earthStateTemp[1] = (*earthPosVel)[1] + hval[i-1]*K[9][i-1];
            earthStateTemp[2] = (*earthPosVel)[2] + hval[i-1]*K[10][i-1];
            earthStateTemp[3] = (*earthPosVel)[3] + hval[i-1]*K[11][i-1];
            
            deriv(&spaceStateTemp, &lunarStateTemp, &earthStateTemp, G, m_E, m_M, m_S);
            
            K[0][i]  = spaceStateTemp[0];
            K[1][i]  = spaceStateTemp[1];
            K[2][i]  = spaceStateTemp[2];
            K[3][i]  = spaceStateTemp[3];
            K[4][i]  = lunarStateTemp[0];
            K[5][i]  = lunarStateTemp[1];
            K[6][i]  = lunarStateTemp[2];
            K[7][i]  = lunarStateTemp[3];
            K[8][i]  = earthStateTemp[0];
            K[9][i]  = earthStateTemp[1];
            K[10][i] = earthStateTemp[2];
            K[11][i] = earthStateTemp[3];
        }

        (*spacePosVel)[0] = (*spacePosVel)[0] + (h/6)*(K[0][1]  + 2*K[0][2]  + 2*K[0][3]  + K[0][4]);
        (*spacePosVel)[1] = (*spacePosVel)[1] + (h/6)*(K[1][1]  + 2*K[1][2]  + 2*K[1][3]  + K[1][4]);
        (*spacePosVel)[2] = (*spacePosVel)[2] + (h/6)*(K[2][1]  + 2*K[2][2]  + 2*K[2][3]  + K[2][4]);
        (*spacePosVel)[3] = (*spacePosVel)[3] + (h/6)*(K[3][1]  + 2*K[3][2]  + 2*K[3][3]  + K[3][4]);
        
        (*moonPosVel)[0] =  (*moonPosVel)[0]  + (h/6)*(K[4][1]  + 2*K[4][2]  + 2*K[4][3]  + K[4][4]);
        (*moonPosVel)[1] =  (*moonPosVel)[1]  + (h/6)*(K[5][1]  + 2*K[5][2]  + 2*K[5][3]  + K[5][4]);
        (*moonPosVel)[2] =  (*moonPosVel)[2]  + (h/6)*(K[6][1]  + 2*K[6][2]  + 2*K[6][3]  + K[6][4]);
        (*moonPosVel)[3] =  (*moonPosVel)[3]  + (h/6)*(K[7][1]  + 2*K[7][2]  + 2*K[7][3]  + K[7][4]);
        
        (*earthPosVel)[0] = (*earthPosVel)[0] + (h/6)*(K[8][1]  + 2*K[8][2]  + 2*K[8][3]  + K[8][4]);
        (*earthPosVel)[1] = (*earthPosVel)[1] + (h/6)*(K[9][1]  + 2*K[9][2]  + 2*K[9][3]  + K[9][4]);
        (*earthPosVel)[2] = (*earthPosVel)[2] + (h/6)*(K[10][1] + 2*K[10][2] + 2*K[10][3] + K[10][4]);
        (*earthPosVel)[3] = (*earthPosVel)[3] + (h/6)*(K[11][1] + 2*K[11][2] + 2*K[11][3] + K[11][4]);
        
        dSE  = pow(pow((*earthPosVel)[0] - (*spacePosVel)[0],2) + pow((*earthPosVel)[1] - (*spacePosVel)[1],2),.5);
        dME  = pow(pow((*earthPosVel)[0] - (*moonPosVel)[0],2)  + pow((*earthPosVel)[1] - (*moonPosVel)[1],2),.5);
        dSM  = pow(pow((*spacePosVel)[0] - (*moonPosVel)[0],2) + pow((*spacePosVel)[1] - (*moonPosVel)[1],2),.5);
        
        *t += h;
        
        if (plot == 1){
            fprintf(outputFile, "%g %g %g %g %g " , *t, (*spacePosVel)[0], (*spacePosVel)[1], (*spacePosVel)[2], (*spacePosVel)[3]);
            fprintf(outputFile, "%g %g %g %g " ,        (*moonPosVel)[0] , (*moonPosVel)[1] , (*moonPosVel)[2] , (*moonPosVel)[3]);
            fprintf(outputFile, "%g %g %g %g\n",        (*earthPosVel)[0], (*earthPosVel)[1], (*earthPosVel)[2], (*earthPosVel)[3]);
        }
        
        c1 = dSE - r_E;
        c2 = dSM - r_M;
        c3 = 2*dME - dSE;
        
    }
    
    
    fclose(outputFile);
    
    if (c1 <= 0 && *t <= *tmin){
        *tmin = *t;
        return 1;
    }
    else{
        return 0;
    }
    
}



int main(int argc, const char * argv[]) {
    
    if(argc != 4){
        printf("Incorrect Number of Command Line Arguments\n");
        return 1;
    }
    
    dirname((char*)argv[0]);
    
    double objective = atof(argv[1]);
    double clearance = atof(argv[2]);
    double tol = atof(argv[3]);
    
    int path_Len  = (int)strlen(argv[0]);
    int objecLen  = (int)strlen(argv[1]);
    int clearLen  = (int)strlen(argv[2]);
    int tolerLen  = (int)strlen(argv[3]);
    int outputLen = path_Len + objecLen + clearLen + tolerLen + 17;
    
    char path_String[path_Len];
    char objecString[objecLen];
    char clearString[clearLen];
    char tolerString[tolerLen];
    char outputString[outputLen];
    
    strcpy(path_String, dirname((char*)argv[0]));
    strcpy(objecString, argv[1]);
    strcpy(clearString, argv[2]);
    strcpy(tolerString, argv[3]);
    

    
    for (int i = 0; i < clearLen; i++) {
        if(!strncmp(&clearString[i],".",1)){
            clearString[i] = 'p';
        }
    }
    
    for (int i = 0; i < tolerLen; i++) {
        if(!strncmp(&tolerString[i],".",1)){
            tolerString[i] = 'p';
        }
    }
    
    strcpy(outputString, path_String);
    strcat(outputString, "/Output");
    mkdir(outputString,07777);
    strcat(outputString, "/Optimum_");
    strcat(outputString, objecString);
    strcat(outputString, "_");
    strcat(outputString, clearString);
    strcat(outputString, "_");
    strcat(outputString, tolerString);
    
    
    
    //Initial Conditions//
    // Masses, Radius, and constants
    double m_E = 5.97219E24;    // (kg)
    double m_M = 7.34767309E22; // (kg)
    double m_S = 28833;              // (kg)
    double r_E = 6371000;            // (m)
    double r_M = 1737100 + clearance;            // (m)
    double G   = 6.674E-11;     // (N*m^2/kg^2)
    
    // Space Craft
    double dES    = 340000000;       // (m)
    double thetaS = 50*pi/180;       // (rad)
    double Vs     = 1000;            // (m/s)
    double Xs0    = dES*cos(thetaS); // (m)
    double Ys0    = dES*sin(thetaS); // (m)
    double Vsx    = Vs*cos(thetaS);  // (m/s)
    double Vsy    = Vs*sin(thetaS);  // (m/s)
    
    // Moon
    double dEM    = 384403000;                           // (m)
    double thetaM = 42.5*pi/180;                         // (rad)
    double Vm     = sqrt(G*pow(m_E,2)/((m_E+m_M)*dEM));  // (m/s)
    double Xm0    = dEM*cos(thetaM);                     // (m)
    double Ym0    = dEM*sin(thetaM);                     // (m)
    double Vmx    = -Vm*sin(thetaM);                     // (m/s)
    double Vmy    = Vm*cos(thetaM);                      // (m/s)
    
    double Vsxnew;
    double Vsynew;
    
    // Pos, Vel and Var arrays
    double earthPosVel[4] = {0,0,0,0};
    double spacePosVel[4] = {Xs0, Ys0, Vsx, Vsy};
    double lunarPosVel[4] = {Xm0, Ym0, Vmx, Vmy};
    double varArray[6]    = {m_E, m_M, m_S, r_E, r_M, G};
    double *earthPointer = earthPosVel;
    double *spacePointer = spacePosVel;
    double *moonPointer  = lunarPosVel;
    double t = 0;
    double tmin = INFINITY;
    
    double timeTemp[2];
    
    
    if (objective == 2){
        for(double vel = 0; vel <= 100; vel += tol/2){
            for (double ang = 0; ang <= 2*pi; ang += pi/40) {
                t = 0;
                Vsxnew = Vsx + vel*cos(ang);
                Vsynew = Vsy + vel*sin(ang);
                
                *(spacePointer + 0) = Xs0;
                *(spacePointer + 1) = Ys0;
                *(spacePointer + 2) = Vsxnew;
                *(spacePointer + 3) = Vsynew;
                
                *(moonPointer + 0) = Xm0;
                *(moonPointer + 1) = Ym0;
                *(moonPointer + 2) = Vmx;
                *(moonPointer + 3) = Vmy;
                
                *(earthPointer + 0) = 0;
                *(earthPointer + 1) = 0;
                *(earthPointer + 2) = 0;
                *(earthPointer + 3) = 0;
                
                if(rungeKutta(30, &earthPosVel, &spacePosVel, &lunarPosVel, &varArray, &t, &tmin, (int)0, (const char*)outputString)){
                    timeTemp[0] = vel;
                    timeTemp[1] = ang;
                }
            }
        }
        t = 0;
        tmin = INFINITY;
        Vsxnew = Vsx + timeTemp[0]*cos(timeTemp[1]);
        Vsynew = Vsy + timeTemp[0]*sin(timeTemp[1]);
        
        *(spacePointer + 0) = Xs0;
        *(spacePointer + 1) = Ys0;
        *(spacePointer + 2) = Vsxnew;
        *(spacePointer + 3) = Vsynew;
        
        *(moonPointer + 0) = Xm0;
        *(moonPointer + 1) = Ym0;
        *(moonPointer + 2) = Vmx;
        *(moonPointer + 3) = Vmy;
        
        *(earthPointer + 0) = 0;
        *(earthPointer + 1) = 0;
        *(earthPointer + 2) = 0;
        *(earthPointer + 3) = 0;
        
        rungeKutta(30, &earthPosVel, &spacePosVel, &lunarPosVel, &varArray, &t, &tmin, (int)1, (const char*)outputString);
        printf("Optimized Velocity and Angle for Objective 2:\nVelocity: %g (m/s)\nAngle: %g (deg)\nTime: %g (s)\n", timeTemp[0], timeTemp[1]*180/pi, t);
        return 0;
    }
    
    else{
        for(double vel = 0; vel <= 100; vel += tol/2){
            for (double ang = 0; ang <= 2*pi; ang += pi/40) {
                t = 0;
                Vsxnew = Vsx + vel*cos(ang);
                Vsynew = Vsy + vel*sin(ang);
                
                *(spacePointer + 0) = Xs0;
                *(spacePointer + 1) = Ys0;
                *(spacePointer + 2) = Vsxnew;
                *(spacePointer + 3) = Vsynew;
                
                *(moonPointer + 0) = Xm0;
                *(moonPointer + 1) = Ym0;
                *(moonPointer + 2) = Vmx;
                *(moonPointer + 3) = Vmy;
                
                *(earthPointer + 0) = 0;
                *(earthPointer + 1) = 0;
                *(earthPointer + 2) = 0;
                *(earthPointer + 3) = 0;
                
                if(rungeKutta(30, &earthPosVel, &spacePosVel, &lunarPosVel, &varArray, &t, &tmin, (int)0, (const char*)outputString)){
                    t = 0;
                    Vsxnew = Vsx + vel*cos(ang);
                    Vsynew = Vsy + vel*sin(ang);
                    
                    *(spacePointer + 0) = Xs0;
                    *(spacePointer + 1) = Ys0;
                    *(spacePointer + 2) = Vsxnew;
                    *(spacePointer + 3) = Vsynew;
                    
                    *(moonPointer + 0) = Xm0;
                    *(moonPointer + 1) = Ym0;
                    *(moonPointer + 2) = Vmx;
                    *(moonPointer + 3) = Vmy;
                    
                    *(earthPointer + 0) = 0;
                    *(earthPointer + 1) = 0;
                    *(earthPointer + 2) = 0;
                    *(earthPointer + 3) = 0;
                    
                    rungeKutta(30, &earthPosVel, &spacePosVel, &lunarPosVel, &varArray, &t, &tmin, (int)1, (const char*)outputString);
                    printf("Optimized Velocity and Angle for Objective 1:\nVelocity: %g (m/s)\nAngle: %g (deg)\n", vel, ang*180/pi);
                    return 0;
                }
            }
        }
    }
    
    
    return 0;
}

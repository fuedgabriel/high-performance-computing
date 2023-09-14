#include <iostream>

using namespace std;

void HLLC(double rL,double ruL,double rEL,double rR,double ruR,double rER,double* Fx, double gamma) {
    double uL = ruL / rL;
    double uR = ruR / rR;
    double pL = (gamma-1.0) * (rEL - 0.5 * rL * uL * uL);
    double pR = (gamma-1.0) * (rER - 0.5 * rR * uR * uR);
    double aL = sqrt(gamma*pL/rL);
    double aR = sqrt(gamma*pR/rR);

    double SL = min(uL-aL,uR-aR);
    double SR = max(uL+aL,uR+aR);

    double SS = (pR-pL+ruL*(SL-uL)-ruR*(SR-uR))/(rL*(SL-uL)-rR*(SR-uR));

    if (SL>0.0) {
        Fx[0] = ruL;
        Fx[1] = ruL*uL+pL;
        Fx[2] = uL*(rEL+pL);
        return;
    } else if (SS>0.0) {
        double rsL = rL * (SL-uL)/(SL-SS);
        double rusL = rsL * SS;
        double rEsL = rsL * (rEL / rL + (SS-uL)*(SS+pL/rL/(SL-uL)));
        Fx[0] = ruL + SL * (rsL - rL);
        Fx[1] = (ruL*uL+pL) + SL * (rusL-ruL);  
        Fx[2] = uL*(rEL+pL) + SL * (rEsL-rEL);
        return;
    } else if (SR>0.0) {
        double rsR = rR * (SR-uR)/(SR-SS);
        double rusR = rsR * SS;
        double rEsR = rsR * (rER / rR + (SS-uR)*(SS+pR/rR/(SR-uR)));
        Fx[0] = ruR + SR * (rsR - rR);
        Fx[1] = (ruR*uR+pR) + SR * (rusR-ruR);
        Fx[2] = uR*(rER+pR) + SR * (rEsR-rER);
        return;
    } else {
        Fx[0] = ruR;
        Fx[1] = ruR*uR+pR;
        Fx[2] = uR*(rER+pR);
        return;
    }
}
//######################################################################

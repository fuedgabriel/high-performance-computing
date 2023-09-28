#include <iostream>
#include <fstream>
#include <cmath>
#include "problems.h"
#include <chrono>
using namespace std;

int dims;
int grid[3];
int N,Cycle;
double Time,EndTime;
double* UVector;
double* nUVector;
double* Flux;
double* X;
double dX;
double* P;
int* neighbour;
double Cgama;

class Timer
{
public:
    void start()
    {
        m_StartTime = std::chrono::system_clock::now();
        m_bRunning = true;
    }
    
    void stop()
    {
        m_EndTime = std::chrono::system_clock::now();
        m_bRunning = false;
    }
    
    double elapsedMilliseconds()
    {
        std::chrono::time_point<std::chrono::system_clock> endTime;
        
        if(m_bRunning)
        {
            endTime = std::chrono::system_clock::now();
        }
        else
        {
            endTime = m_EndTime;
        }
        
        return std::chrono::duration_cast<std::chrono::milliseconds>(endTime - m_StartTime).count();
    }
    
    double elapsedSeconds()
    {
        return elapsedMilliseconds() / 1000.0;
    }

private:
    std::chrono::time_point<std::chrono::system_clock> m_StartTime;
    std::chrono::time_point<std::chrono::system_clock> m_EndTime;
    bool                                               m_bRunning = false;
}; 



//######################################################################
int initialize(int Dims,int *Grid) {
    int N = 0;
    switch(Dims) {
        case 1:
            N = Grid[0];
            break;
        case 2:
            N = Grid[0]*Grid[1];
            break;
        case 3:
            N = Grid[0]*Grid[1]*Grid[2];
            break;
        default:
            cout << " Error in Dims " << endl;
            exit(0);
            return -1;
            break;
    }
    //-- 
    return N;
}
//######################################################################
void Problem1() {
    dims = 1;
    grid[0] = 100;
    //-- initialization
    N = initialize(dims,grid);

    UVector = new double [N*(dims+2)] ;
    nUVector = new double [N*(dims+2)] ;
    Flux = new double [(N+1)*(dims+2)] ;
    P = new double [N] ;
    X = new double [N*dims] ;
    neighbour  = new int [N*2*dims] ;

    //-------------------------------
    double Length = 1.0;
    Cgama = 1.4;
    EndTime = 0.25;

    dX = Length / grid[0];
    for (int i=0;i<grid[0];i++) {
        X[i] = Length / (grid[0]-1.0) * i + dX / 2;
        if (X[i]<Length/2.0) {
            //1.0,0.0,1.0
            UVector[i*3+0] = 1.0;
            UVector[i*3+1] = 0.0;
            UVector[i*3+2] = 1.0 / 0.4;
            P[i] = 1.0;
        } else {
            //0.125,0,0.1
            UVector[i*3+0] = 0.125;
            UVector[i*3+1] = 0.0;
            UVector[i*3+2] = 0.1 / 0.4;
            P[i] = 0.1;
        }
        neighbour[2*dims*i+0] = 0;
        neighbour[2*dims*i+1] = 0;
    }
    neighbour[0+0] = 1;
    neighbour[2*dims*(grid[0]-1)+1] = 1;
    //-------------------------------

    return ;
}
//######################################################################
void Output(int C,double T,int Nt, double* U, double* P,double* Pos) {
    char fileName[150];
    sprintf(fileName, "Output/Cycle_%06d.dat", C);
    ofstream myfile(fileName);

    myfile << "# Time = " << T << "  Cycle = " << C << endl; 
    myfile << "# X,UVector,P" << endl;

    for (int n = 0; n<Nt;n++) {
        #pragma opm parallel
        {
        #pragma omp for
        for (int j=0;j<dims;j++)
            myfile << Pos[n*dims+j] << ",";
        #pragma omp for
        for (int j=0;j<2+dims;j++)
            myfile << U[n*3+j] << ",";        
        
        myfile << P[n] << endl;       
        } 
    }
    myfile.close();
}
//######################################################################
double CFL(int Nt,double* U,double CT,double CET) {
    double maxVelocity = -1000.0; 
    for (int n = 0;n<Nt;n++){
        double c = sqrt(Cgama*P[n]/U[n*(2+dims)+0]);
        for (int d=0;d<dims;d++){
            double Vel = U[n*(2+dims)+1+d] / U[n*(2+dims)+0];
            maxVelocity = max(maxVelocity,abs(Vel+c));
            maxVelocity = max(maxVelocity,abs(-Vel+c));
        }
    }
    
    double ndT =  dX / maxVelocity * 0.1;
    if (CT+ndT>CET) ndT = CET - CT + 1.0e-8;
    return ndT;
}
//######################################################################
void HLLC(double rL,double ruL,double rEL,double rR,double ruR,double rER,double* Fx) {
    double uL = ruL / rL;
    double uR = ruR / rR;
    double pL = (Cgama-1.0) * (rEL - 0.5 * rL * uL * uL);
    double pR = (Cgama-1.0) * (rER - 0.5 * rR * uR * uR);
    double aL = sqrt(Cgama*pL/rL);
    double aR = sqrt(Cgama*pR/rR);

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
void boundaryCondition(int TYP,double rC,double ruC,double rEC,double* r,double* ru,double* rE) {
    if (TYP==1) {
        *r =rC;
        *ru =-ruC;
        *rE = rEC;
    } else {
        cout << "error in BND" << endl;
        exit(0);
    }
}
//######################################################################
//######################################################################
int main() {
    Timer timer;

    timer.start();
    //-- initialization
    Problem1();
    Time = 0.0;
    Cycle = 0;
    double dT = CFL(N,UVector,Time,EndTime);
    //-- initializtion
    Output(Cycle,Time,N,UVector,P,X);

    while (Time<EndTime && Cycle < 1000) {
        Time = Time + dT;
        //-- calc Flux
        for(int n=0;n<N;n++) {
            int LEFT = neighbour[n*dims*2+0];
            int RIGHT = neighbour[n*dims*2+1];

            double rL,ruL,rEL;
            double rC = UVector[n*(2+dims)+0],ruC = UVector[n*(2+dims)+1],rEC = UVector[n*(2+dims)+2];
            double rR,ruR,rER;

            if (LEFT==0) {
                //#n-1 is the left
                rL  = UVector[(n-1)*(2+dims)+0];
                ruL = UVector[(n-1)*(2+dims)+1];
                rEL = UVector[(n-1)*(2+dims)+2];
            } else {
                //# the left is boundary
                boundaryCondition(LEFT,rC,ruC,rEC,&rL,&ruL,&rEL);
            }
            //-- Flux n is LC
            HLLC(rL,ruL,rEL,rC,ruC,rEC,Flux+n*(2+dims));


            if (RIGHT==0) {
                //-- n+1 is the right [DO NOTHING]
            } else {
                //-- the right is boundary
                boundaryCondition(RIGHT,rC,ruC,rEC,&rR,&ruR,&rER);
                //-- Flux n+1 is CR
                HLLC(rC,ruC,rEC,rR,ruR,rER,Flux+(n+1)*(2+dims));
            }
        }
        //-- update new UVector
        #pragma omp parallel for collapse(2)
        for (int n=0;n<N;n++) 
            for (int i=0;i<2+dims;i++) {
                UVector[n*(2+dims)+i] = UVector[n*(2+dims)+i] - dT/dX * (Flux[(n+1)*(2+dims)+i]-Flux[n*(2+dims)+i]);            
            
            //cout << n << "  " << UVector[n*(2+dims)+0] << "  " << UVector[n*(2+dims)+1] << " @ " << UVector[n*(2+dims)+2] << endl;
            double re = UVector[n*(2+dims)+(1+dims)] - 0.5 * UVector[n*(2+dims)+1] * UVector[n*(2+dims)+1] / UVector[n*(2+dims)+0];
            P[n] = (Cgama-1.0)*re;

        }
        //-- CFL
        dT = CFL(N,UVector,Time,EndTime);
        //cout << Cycle << "  " << Time << "  " << dT << endl;
        Cycle++;
        if (Cycle%1==0) Output(Cycle,Time,N,UVector,P,X);
    }

    Output(Cycle,Time,N,UVector,P,X);

    timer.stop();
    cout << "\nSeconds: " << timer.elapsedSeconds() << "\n";
    cout << "Milliseconds: " << timer.elapsedMilliseconds() << "\n";


    delete [] UVector;
    delete [] nUVector;
    delete [] Flux;
    delete [] P;
    delete [] X;
    delete [] neighbour;
    return 0;
}
//######################################################################

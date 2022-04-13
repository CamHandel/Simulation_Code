#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <math.h>
#include <vector>

using namespace std;

int main()
{
    string line,fName;
    int    id,n,t,tTime,nAtoms,yStep,currentTimeStep,typ;
    double delta,sine,sine1,forceY;

    double kB = 1.38064852e-23;
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength/refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;
   
    
    //****  INPUT *******************************************************************
    
    int nTimeSteps=100000;
    
    int tSkip= 1; // dump output
    
    double deltaT=2.0; // time step
    double amp=9.0;
    double Time_period=62100;

    //*******************************************************************************

    ifstream data("../forceY.txt",ios::in);

    ofstream dataFile("time_vs_energyInput.txt",ios::out);

    double energyInput = 0.0;
    
    for(t=0;t<nTimeSteps;t++)
    {
      delta = t*deltaT;
      sine = amp*sin((2*3.14159265359 / Time_period) * delta);
      sine1 = amp*sin((2*3.14159265359 / Time_period) * (delta- deltaT));
        
      data >> currentTimeStep >> forceY;

      energyInput = energyInput + forceY*(sine1 -sine);

      dataFile<< deltaT*t*tSkip <<'\t'
                << energyInput
                << endl;
                                      
              
            getline(data,line);
        }
        
        

    return 0;
}

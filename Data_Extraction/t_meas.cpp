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
    double x,y,z,vx,vy,vz;
    double Lx,Ly,Lz,xLo,xHi,yLo,yHi,zLo,zHi;
    double vxAvgSol,vyAvgSol,vzAvgSol,vxAvgWall,vyAvgWall,vzAvgWall; 
    int    nMolSol, nMolWall;
    double keSol,keWall,keMacroSol,volume;

    double kB = 1.38064852e-23;
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength/refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;
   
    
    //****  INPUT *******************************************************************

    double miWall = 195.084;// mass of one molecule
    double miSol = 195.084;
    
    int nTimeSteps=10000; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    
    int tSkip= 10; // dump output
    
    double deltaT=2.0; // time step

    //*******************************************************************************
    
    vector<double> velX;    
    vector<double> velY;
    vector<double> velZ;
    vector<double> type;






    ifstream data("../dump_meas.lammpstrj",ios::in);

    ofstream dataFile("1_ke_vs_time.txt",ios::out);
    
    for(t=0;t<nTimeSteps;t++)
    {
        for(n=1;n<10;n++)
        {
            if(n == 2)
            {
                data >> currentTimeStep;
                
                cout << "currentTimeStep = " << currentTimeStep
                    << "; t = " << t << " [ " 
                    << 100*float(t+1)/float(nTimeSteps) 
                    << "% ]" << endl;    
            }

            if(n == 4)
            {
                data >> nAtoms;
//                 cout << "nAtoms = " << nAtoms << endl; 
            }
            
            if(n == 6)
            {
                data >> xLo >> xHi;
            }
            if(n == 7)
            {
                data >> yLo >> yHi;
            }
            if(n == 8)
            {
                data >> zLo >> zHi;
            }                      
              
            getline(data,line);
        }
        
	/*       if (t == 0)
        {
            Lx = xHi-xLo;
            Ly = yHi-yLo;
            Lz = zHi-zLo;
	    } */
        
        vxAvgSol = 0;
        vyAvgSol = 0;
	vzAvgSol = 0;
	vxAvgWall = 0;
        vyAvgWall = 0;
	vzAvgWall = 0;
	keMacroSol=0.0;
        keSol=0.0;
        keWall=0.0;
        nMolSol = 0;
	nMolWall = 0;
	velX.resize(nAtoms, 0.0);
	velY.resize(nAtoms, 0.0);
	velZ.resize(nAtoms, 0.0);
	type.resize(nAtoms, 0.0);

        for(n=0;n<nAtoms;n++)
        {
            data>>id>>typ>>x>>y>>z>>vx>>vy>>vz;
            
            if(typ == 1) // wall
            {
                nMolWall++;
                vxAvgWall += vx;
                vyAvgWall += vy;
                vzAvgWall += vz;
                velX[n] = vx;
		velY[n] = vy;
		velZ[n] = vz;
		type[n] = typ;
            }


	    if(typ == 2) // sol
	    {

	        nMolSol++;
	        vxAvgSol += vx;
                vyAvgSol += vy;
                vzAvgSol += vz;
		velX[n] = vx;
		velY[n] = vy;
		velZ[n] = vz;
                type[n] = typ;
            
            }
        }  
        
    //       volume = Lx*Ly*Lz;
        

        vxAvgSol = vxAvgSol/nMolSol;
        vyAvgSol = vyAvgSol/nMolSol;
	vzAvgSol = vzAvgSol/nMolSol;
	vxAvgWall = vxAvgWall/nMolWall;
        vyAvgWall = vyAvgWall/nMolWall;
	vzAvgWall = vzAvgWall/nMolWall;

	keMacroSol = 0.5*miSol*(vxAvgSol*vxAvgSol + vyAvgSol*vyAvgSol + vzAvgSol*vzAvgSol); 



	 for(n=0;n<nAtoms;n++)
	 {   
            if(type[n] == 1) // wall
            {
  
	      keWall +=0.5*miWall*((velX[n]-vxAvgWall)*(velX[n]-vxAvgWall) + (velY[n]-vyAvgWall)*(velY[n]-vyAvgWall) + (velZ[n]-vzAvgWall)*(velZ[n]-vzAvgWall));
		
            }


	    if(type[n] == 2) // sol
	    {

	      keSol+=0.5*miSol*((velX[n]-vxAvgSol)*(velX[n]-vxAvgSol) + (velY[n]-vyAvgSol)*(velY[n]-vyAvgSol) + (velZ[n]-vzAvgSol)*(velZ[n]-vzAvgSol));
            
            }
        }

	 
        
        
        getline(data,line); 
        
        dataFile<< deltaT*t*tSkip <<'\t'
                << keMacroSol <<'\t'
                << keSol <<'\t'
                << keWall
                << endl;                   
    }
     
    


    return 0;
}

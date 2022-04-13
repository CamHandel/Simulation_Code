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
    double keSol, keWall,keMacroSol,volume;
    double VelocitySol, VelocityWall, VelocityRel;

    double kB = 1.38064852e-23;
    double refMass = 1.66054e-27;
    double refLength = 1e-10;
    double refTime = 1e-15;
    double refVelocity = refLength/refTime; // 1 Angstrom / 1 fs
    double refPressure = 101325;
   
    
    //****  INPUT *******************************************************************

    double miWall = 195.084;// mass of one molecule
    double miSol = 195.084;
    
    int nTimeSteps= 1000; // CHANGE, use command line: grep -o 'TIMESTEP' dump_meas.lammpstrj | wc -l
    // 1000 100 or 400 and 5
    int tSkip= 100; // dump output
    
    double deltaT=2.0; // time step
    
    cout <<  "The number of timesteps:" << nTimeSteps*tSkip << endl;

    //*******************************************************************************
    
    vector<double> velX;    
    vector<double> velY;
    vector<double> velZ;
    vector<double> type;
    
    vector<double> posX;    
    vector<double> posY;
    vector<double> posZ;
    
    vector<double> cornerY;






    ifstream data("../dump_meas.lammpstrj",ios::in);

    ofstream rotDataFile("rotation.txt",ios::out);
    ofstream velDataFile("velocity.txt",ios::out);
    ofstream kinDataFile("kinetic.txt",ios::out);
    ofstream cornerDataFile("corners.txt",ios::out);
    
    
    for(t=0;t<nTimeSteps;t++)
    {
    	// Reads through the start of the file gathering simulation data
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
	
	posX.resize(2, 0.0);
	posY.resize(2, 0.0);
	posZ.resize(2, 0.0);
	
	cornerY.resize(9, 0.0);
	
	VelocitySol = 0.0;
	VelocityWall = 0.0;
	VelocityRel = 0.0;
	

	
	// loops through atom data from a single timestep
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
                
                vx *= pow(10,5);
                vy *= pow(10,5);
                vz *= pow(10,5);
                // (1.0/4184.0 * 6.0221 * pow(10.0,23.0)) *
                //  * (vx*vx + vy*vy + vz*vz)) 0.5*miSol*
                keSol +=  (1.0/4184.0 * 6.0221 * pow(10.0,23.0)) * (0.5 * miSol * 1.6605 * pow(10,-27) * (vx*vx + vy*vy + vz*vz));

            }
            
            if(id == 10423){
            	posX[0] = x;
               posY[0] = y;
               posZ[0] = z;
            }
            
            if(id == 10884){
            	posX[1] = x;
               posY[1] = y;
               posZ[1] = z;
            }
            
            // corners added
            if(id == 10531){
            	cornerY[0] = y;
            }
            
            if(id == 10884){
            	cornerY[1] = y;
            }
            
            if(id == 10971){
            	cornerY[2] = y;
            }
            
            if(id == 10952){
            	cornerY[3] = y;
            }
            
            if(id == 10512){
            	cornerY[4] = y;
            }
            
            if(id == 10369){
            	cornerY[5] = y;
            }
            
            if(id == 10378){
            	cornerY[6] = y;
            }
            
            if(id == 10865){
            	cornerY[7] = y;
            }
            
            if(id == 3888){
            	cornerY[8] = y;
            }
            
        }  
        
    //       volume = Lx*Ly*Lz;

        vxAvgSol = vxAvgSol/nMolSol;
        vyAvgSol = vyAvgSol/nMolSol;
	vzAvgSol = vzAvgSol/nMolSol;
	vxAvgWall = vxAvgWall/nMolWall;
        vyAvgWall = vyAvgWall/nMolWall;
	vzAvgWall = vzAvgWall/nMolWall;

	VelocityRel = vyAvgSol-vyAvgWall;
	
	keMacroSol = 1.0/4184.0 * 6.0221 * pow(10,23) * (0.5*miSol*(1.6605*pow(10.0,-27.0)) * (vxAvgSol*vxAvgSol + vyAvgSol*vyAvgSol + vzAvgSol*vzAvgSol)); 
	
	/*

	 if(type[n] == 2) // sol
	 {

	      keSol+=0.5*miSol*((velX[n]-vxAvgSol)*(velX[n]-vxAvgSol) + (velY[n]-vyAvgSol)*(velY[n]-vyAvgSol) + (velZ[n]-vzAvgSol)*(velZ[n]-vzAvgSol));
            
         }
         
         */

	 
        
        
        getline(data,line); 
        
        velDataFile<< deltaT*t*tSkip <<'\t'
       	<< VelocityRel <<'\t'
       	<< vyAvgSol <<'\t'
       	<< vyAvgWall
                << endl;
                
        rotDataFile<< deltaT*t*tSkip <<'\t'
        	 << posX[0] <<'\t'
        	 << posX[1] <<'\t'
                << posY[0] <<'\t'
                << posY[1] <<'\t'
                << posZ[0] <<'\t'
                << posZ[1]
                << endl;                   
    
    	kinDataFile<< deltaT*t*tSkip <<'\t'
        	 << keSol
                << endl;
                
        cornerDataFile << deltaT*t*tSkip <<'\t'
        	<< cornerY[0] <<'\t'
        	<< cornerY[1] <<'\t'
        	<< cornerY[2] <<'\t'
        	<< cornerY[3] <<'\t'
        	<< cornerY[4] <<'\t'
        	<< cornerY[5] <<'\t'
        	<< cornerY[6] <<'\t'
        	<< cornerY[7] <<'\t'
        	<< cornerY[8]
        	<< endl;
    
    }
     
    


    return 0;
}

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TGeoManager.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoVolume.h>
#include <TGeoBBox.h>
#include <TGeoTube.h>
#include <TGeoPcon.h>
#include <TGeoHalfSpace.h>
#include <TGeoMatrix.h>
#include <TGeoCompositeShape.h>
#include <TGraph.h>
#include <TFile.h>
#include <TMath.h>
#include <TTree.h>

#include "ComponentComsol.hh"
#include "ViewField.hh"
#include "ViewFEMesh.hh"
#include "MediumMagboltz.hh"
#include "Sensor.hh"
#include "AvalancheMicroscopic.hh"
#include "AvalancheMC.hh"
#include "Random.hh"
#include "Plotting.hh"
#include "ViewSignal.hh"

#include <time.h>
using namespace Garfield;
using namespace std;


void ReadGas(const char *filename);

int main(){

    // ReadGas("Ar-CO2(90-10).gas");
    
    for(int i=11;i>0;i--){
        char filename[50];
        sprintf( filename, "%s%d%s%d%s", "Ar-CO2(", 100-i, "-", i, ").gas");
        ReadGas(filename);
    }
    
    return 0;
}

void ReadGas(const char *filename) {
    MediumMagboltz *gas;
    gas = new MediumMagboltz();

    gas->LoadGasFile(filename);

    const double DriftLength = 20; // unit cm
    double ex = 10;
    double vx, vy, vz;
    double t[5], dt[5], dvx[5], dex[5], zRes[5];
	double dl, dt;
	double alpha; 
	double eta;

	while(ex < 2000) {
        gas->ElectronVelocity  (ex, 0, 0, 0, 0, 0, vx, vy, vz);
        gas->ElectronDiffusion (ex, 0, 0, 0, 0, 0, dl, dt);
        gas->ElectronTownsend  (ex, 0, 0, 0, 0, 0, alpha);
        gas->ElectronAttachment(ex, 0, 0, 0, 0, 0, eta);

        vx = vx * 1000; // unit cm/us
        t[0] = DriftLength / vx; // DriftLengthcm drift zone; unit us
        dt[0] = 0;
        dex[0] = 0;
        dvx[0] = 0;
        zRes[0] = 0;


        dex[1] = ex + 0.1; // Electric field fluctuation 0.1V/cm
        gas->ElectronVelocity  (dex[1], 0, 0, 0, 0, 0, dvx[1], vy, vz);

        dex[2] = ex + 1; // Electric field fluctuation 1V/cm
        gas->ElectronVelocity  (dex[2], 0, 0, 0, 0, 0, dvx[2], vy, vz);

        dex[3] = ex * 1.01; // Electric field fluctuation 1%
        gas->ElectronVelocity  (dex[3], 0, 0, 0, 0, 0, dvx[3], vy, vz);

        dex[4] = ex * 1.1; // Electric field fluctuation 10%
        gas->ElectronVelocity  (dex[4], 0, 0, 0, 0, 0, dvx[4], vy, vz);
        
        for(int i=1;i<5;i++){
            dvx[i] = dvx[i] * 1000; // unit cm/us
            t[i] = DriftLength / dvx[i]; 
            dt[i] = (t[i]-t[0])/1000; // unit ns
            zRes[i] = (dvx[i]*t[0] - DriftLength) * 10; // unit mm
            dvx[i] = dvx[i] - vx;
        }
        

        cout << "E = " << ex << "V/cm, vx = " << vx << "cm/us. Under 1V/cm fluctuation, v_uncertain = " << dvx[2] << "cm/us, tRes = " << dt[2] << "ns, zRes = " << zRes[2] << "mm" << endl;

        if(ex < 400) ex+=10;
        else if(ex < 1000) ex+=50;
        else ex+=500;
    }

}

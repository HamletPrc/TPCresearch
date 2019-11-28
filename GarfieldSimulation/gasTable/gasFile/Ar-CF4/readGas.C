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

    double ex = 10;
	double vx, vy, vz;
	double dl, dt;
	double alpha; 
	double eta;

	while(ex < 10000) {
        gas->ElectronVelocity  (ex, 0, 0, 0, 0, 0, vx, vy, vz);
        gas->ElectronDiffusion (ex, 0, 0, 0, 0, 0, dl, dt);
        gas->ElectronTownsend  (ex, 0, 0, 0, 0, 0, alpha);
        gas->ElectronAttachment(ex, 0, 0, 0, 0, 0, eta);

	    // cout << "E = " << ex << "V/cm,   Vx =" << vx * 1000 << "cm/us; alpha = " << alpha << "/cm" <<endl;
        cout << "E = " << ex << "V/cm,   Vx = " << vx * 1000 << " cm/us,    Dt = " << dt << " some unit,    Dl = " << dl << " some unit" << endl;
	// cout << "E = " << ex << "V/cm,   Dl =" << dl << "um/sqrt(cm)" << endl;

        if(ex < 400) ex+=10;
        else if(ex < 1000) ex+=50;
        else ex+=500;
    }

}

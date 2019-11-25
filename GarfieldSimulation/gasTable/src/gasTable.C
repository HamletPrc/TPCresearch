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

void GenGasTable(const char *gas1, double rat1, const char *gas2, double rat2, const char *filename);
void GenGasTable(const char *gas1, double rat1, const char *gas2, double rat2, const char *gas3, double rat3, const char *filename);
void GasSet(MediumMagboltz *gas);
void GasPrintOut(MediumMagboltz *gas);


int main(){

    // for(proportion=10;proportion>0;proportion--){
    //     char filename[20];
    //     sprintf( filename, "%s%d%s", "Ar_iC4H10_", proportion,".gas");
    //     GenGasTable("argon", 100-proportion, "iC4H10", proportion, filename);
    // }

    GenGasTable("argon", 90, "CH4", 10, "H2O", 0.000070, "P10+70ppm_H2O.gas");

    // GenGasTable("argon", 90, "CH4", 10, "P10.gas");
        
    return 0;
}




void GenGasTable(const char *gas1, double rat1, const char *gas2, double rat2, const char *filename) {
    MediumMagboltz *gas;
    gas = new MediumMagboltz();

    gas->SetComposition(gas1, rat1, gas2, rat2);
    GasSet(gas);
    gas->WriteGasFile(filename);

    GasPrintOut(gas);
}


void GenGasTable(const char *gas1, double rat1, const char *gas2, double rat2, const char *gas3, double rat3, const char *filename)
{
    MediumMagboltz *gas;
    gas = new MediumMagboltz();

    gas->SetComposition(gas1, rat1, gas2, rat2, gas3, rat3);
    GasSet(gas);
    gas->WriteGasFile(filename);

    GasPrintOut(gas);
}


void GasSet(MediumMagboltz *gas) {
    gas->SetTemperature(293.15);
    gas->SetPressure(760);  //In Huang XueFeng's work, he uses 0.8atm=0.8*760mmHg gas-pressure
    gas->EnableDebugging();
    gas->Initialise();
    gas->DisableDebugging();
    
    gas->SetMaxElectronEnergy(200.); //eV
    gas->SetMaxPhotonEnergy(200.);
    gas->EnableEnergyRangeAdjustment(true);
    gas->EnableAnisotropicScattering();
    
    gas->Initialise(true);

	int nE, nB, nA;              // number of E/B/Angles between E and B
	double eMin, eMax;           // min and max of E [V/cm]
	double bMin, bMax;           // min and max of B [T]
	double aMin, aMax;           // min and max of angle [rad]
	bool logE;                   // use evenly spaced (false) or logarithmically spaced (true)

	nE = 35;
	eMin = 10; 
	eMax = 10000;
	logE = true;

	nB = 1;
	bMin = 0;
	bMax = 0;

	nA = 1;
	aMin = 0;
	aMax =0;
	
	gas->SetFieldGrid(eMin, eMax, nE, logE, bMin, bMax, nB, aMin, aMax, nA);

    // Set the Penning transfer efficiency.
    const double rPenning = 0.57;
    const double lambdaPenning = 0.;
    gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
    // Generate the gas table file for electron drift
    gas->GenerateGasTable(1,false);
}

void GasPrintOut(MediumMagboltz *gas) {
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
        cout << "E = " << ex << "V/cm,   Dl =" << dl << "um/sqrt(cm)" << endl;

        if(ex < 100) ex+=10;
        else if(ex < 1000) ex+=50;
        else ex+=500;
    }
}

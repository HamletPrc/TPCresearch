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

void DefineGas(const char *gas1, double rat1, const char *gas2, double rat2, const char *filename);

int main(){
    int proportion;
    for(proportion=10;proportion>0;proportion--){
        char filename[20];
		sprintf( filename, "%s%d%s", "Ar_iC4H10_", proportion,".gas");
        DefineGas("argon", 100-proportion, "iC4H10", proportion, filename);
    }
        
    return 0;
}




void DefineGas(const char *gas1, double rat1, const char *gas2, double rat2, const char *filename)
{
    MediumMagboltz *gas;
    gas = new MediumMagboltz();

    cout<<"?"<<endl;

    if(!gas->LoadGasFile(filename)) 
    {
        gas->SetComposition(gas1, rat1, gas2, rat2);
        //gas->SetComposition(gas1, rat1);
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
        // Set the Penning transfer efficiency.
        const double rPenning = 0.57;
        const double lambdaPenning = 0.;
        gas->EnablePenningTransfer(rPenning, lambdaPenning, "ar");
        // Generate the gas table file for electron drift
        gas->GenerateGasTable(1,false);
	    gas->WriteGasFile(filename);
    }

    double ex;
	double vx, vy, vz;
	double dl, dt;
	double alpha; 
	double eta;

    for(int i=1; i<100; i++)
	{
	    ex = i*10;
	    gas->ElectronVelocity  (ex, 0, 0, 0, 0, 0, vx, vy, vz);
            gas->ElectronDiffusion (ex, 0, 0, 0, 0, 0, dl, dt);
            gas->ElectronTownsend  (ex, 0, 0, 0, 0, 0, alpha);
            gas->ElectronAttachment(ex, 0, 0, 0, 0, 0, eta);

	    cout<<"E = "<<ex<<"V/cm,   Vx ="<<vx*100<<"cm/us; alpha = "<<alpha<<"/cm"<<endl;
    }
}
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
    
    for(int i=10;i>0;i--){
        char filename[50];
        sprintf( filename, "%s%d%s%d%s", "Ar-iC4H10(", 100-i, "-", i, ").gas");
        ReadGas(filename);
    }
    
    return 0;
}

void ReadGas(const char *filename) {
    MediumMagboltz *gas;
    gas = new MediumMagboltz();

    gas->LoadGasFile(filename);

    double ex = 10.0, dex1, dex2, dex3, dex4;
    double vx, vy, vz, dvx1, dvx2, dvx3, dvx4;
    double t, t1, t2, t3, t4;
    double dl, dt;
    double alpha; 
    double eta;

    while(ex < 10000.0) {
    gas->ElectronVelocity  (ex, 0, 0, 0, 0, 0, vx, vy, vz);
    gas->ElectronDiffusion (ex, 0, 0, 0, 0, 0, dl, dt);
    gas->ElectronTownsend  (ex, 0, 0, 0, 0, 0, alpha);
    gas->ElectronAttachment(ex, 0, 0, 0, 0, 0, eta);
	
    // cout << "E = " << ex << "V/cm,   Vx =" << vx * 1000 << "cm/us; alpha = " << alpha << "/cm" <<endl;
    // cout << "E = " << ex << "V/cm,   Vx = " << vx * 1000 << " cm/us,    Dt = " << dt << " some unit,    Dl = " << dl << " some unit" << endl;
    // cout << "E = " << ex << "V/cm,   Dl =" << dl << "um/sqrt(cm)" << endl;


    dex1 = ex + 0.1;
    gas->ElectronVelocity  (dex1, 0, 0, 0, 0, 0, dvx1, vy, vz);
    //dvx1 = dvx1 - vx;
	
    dex2 = ex + 1;
    gas->ElectronVelocity  (dex2, 0, 0, 0, 0, 0, dvx2, vy, vz);
    //dvx2 = dvx2 - vx;

    dex3 = ex * 1.01;
    gas->ElectronVelocity  (dex3, 0, 0, 0, 0, 0, dvx3, vy, vz);
    //dvx3 = dvx3 - vx;
	
    dex4 = ex * 1.1;
    gas->ElectronVelocity  (dex4, 0, 0, 0, 0, 0, dvx4, vy, vz);
    //dvx4 = dvx4 - vx;
	
    //cout << "E = " << ex << "V/cm,   Vx = " << vx * 1000 << " cm/us,  dvx1 = " << dvx1 * 1000 << " cm/us,  dvx2 =  " << dvx2 * 1000 << " cm/us, dvx3 = " << dvx3 * 1000 << " cm/us, dvx4 = " << dvx4 * 1000 << " cm/us " << endl;


    vx = vx * 1000;
    dvx1 = dvx1 * 1000;
    dvx2 = dvx2 * 1000;
    dvx3 = dvx3 * 1000;
    dvx4 = dvx4 * 1000;
    
    t = 20 / vx; // unit us
    t1 = 20 / dvx1; t1 = (t1-t)*1000; // unit ns
    t2 = 20 / dvx2; t2 = (t2-t)*1000;
    t3 = 20 / dvx3; t3 = (t3-t)*1000;
    t4 = 20 / dvx4; t4 = (t4-t)*1000;

    cout << "E = " << ex << "V/cm,  Vx = " << vx << "cm/us, to drift 20cm, t = " << t << "us, 0,1V/cm fluctuation, dt = " << t1 << "ns, 1V/cm fluctuation, dt = " << t2 << "ns, \%1 fluctuation, dt = " << t3 << "ns, \%10 fluctuation, dt = " << t4 << "ns" << endl;

    if(ex < 400) ex+=10;
    else if(ex < 1000) ex+=50;
    else ex+=500;
    }

    
/*    ex = 50;
    gas->ElectronVelocity  (ex, 0, 0, 0, 0, 0, vx, vy, vz);
    cout << "E = " << ex << "V/cm, vx = " << vx * 1000 << endl;

    ex = 50.1;
    gas->ElectronVelocity  (ex, 0, 0, 0, 0, 0, vx, vy, vz);
    cout << "E = " << ex << "V/cm, vx = " << vx * 1000 << endl;
    
    ex = 51;
    gas->ElectronVelocity  (ex, 0, 0, 0, 0, 0, vx, vy, vz);
    cout << "E = " << ex << "V/cm, vx = " << vx * 1000 << endl;

    ex = 55;
    gas->ElectronVelocity  (ex, 0, 0, 0, 0, 0, vx, vy, vz);
    cout << "E = " << ex << "V/cm, vx = " << vx * 1000 << endl;
*/
}

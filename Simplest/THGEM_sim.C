//------------------------------------------------------------------//
// Author: Qian LIU
//------------------------------------------------------------------//


#include <iostream>
#include <fstream>
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

void ReadFiledMap(TString);
void DefineGas(const char *gas1, double rat1, const char *gas2, double rat2);
void PlotField(ComponentComsol* fm);
void PlotDriftLine(AvalancheMC* drift, ViewDrift* driftView);

ComponentComsol *fm;
MediumMagboltz    *gas;

// Dimensions of the GEM, in !!!! CM !!!!, note that in Ansys it's in MM !
const double pitch  = 0.015;
const double kapton = 0.04;
const double metal  = 0.001;
const double diam   = 0.01;
const double induct = 0.2;
const double drift  = 0.2;
//const double truepitch = 0.1*0.01;

// Dimensions for area selection for DriftView & Sensor etc.
const double xmin = -pitch*3/2;
const double xmax =  pitch*3/2;
const double ymin = -pitch*3/2;
const double ymax =  pitch*3/2;
const double zmin = -induct;
const double zmax =  drift+kapton+2*metal;

char Path[100] = "./";
TString RootFile = TString(Path)+TString("/result.root");
//------------------------------------------------------------------//
// Main
//------------------------------------------------------------------//
int main(int argc, char * argv[]) {
    
    TApplication app("app", &argc, argv);
    plottingEngine.SetDefaultStyle();
    
    // Control flag
    const bool debug = true;
    const bool plotField = false;
    const bool plotDrift = true; 
    const bool plotDriftLine = false;
    const bool plotEFDrift = true;

    //const bool plotGeo = false;
    const bool plotHistogram = true;
    
    //the prim ionization Z pozition. Now I set it to induction field to study IONs!!
    const double zprim  = drift/2+kapton+2*metal;
    //const double zprim  = drift/2+kapton/2+metal;
    //const double zprim = 0.2;
    
	const double tStart =0;
	const double tStop =100;
	const int nSteps =100;
	const double tStep =(tStop - tStart) /nSteps;
	const std::string label = "readout";
    //-----------------
    // Start simulation
    ReadFiledMap(TString(Path));
    
    DefineGas("argon", 90, "C4H10", 10);
	//DefineGas("argon", 90, "isopropanol", 10);
	//DefineGas("DME", 100, "neon", 0);
    
    //---------------------------------------------------------------
    //!!!!! AT LEAST RUN THIS ONCE TO CHECK THE ELECTRIC FIELD !!!!!!
    //!!!!!           IS CORRECTLY SIMULATED BY ANSYS          !!!!!!
    //---------------------------------------------------------------
    // for example, check the anode Z position and potential.
    if(plotField) {PlotField(fm); app.Run(kTRUE);}
    cout << "111" << endl; 
	//fm->SetWeightingField("WPOT.lis", label);
    //------------------------------------------------------
    // Create the sensor. Initialize
    Sensor* sensor = new Sensor();
    sensor->AddComponent(fm);
	//sensor->AddElectrode(fm, label);
	sensor->SetTimeWindow(-1, tStep, 200);
    //sensor->SetArea(-truepitch, -truepitch, zmin, truepitch, truepitch, zmax);	// create a sensor area
    
    sensor->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
	sensor->ClearSignal();
	//sensor->EnableDebugging();
    //------------------------------------------------------
    // Now avalanche and track electron
    AvalancheMicroscopic* aval = new AvalancheMicroscopic();
    aval->SetSensor(sensor);
	aval->EnableSignalCalculation();
    // Now track ions
    AvalancheMC* drift = new AvalancheMC();
    drift->SetSensor(sensor);
    drift->SetDistanceSteps(2.e-4);
    
    //------------------------------------------------------
    // plot drift line.
    ViewDrift* driftView = new ViewDrift();
    if (plotDrift) {
        driftView->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
        // Plot every 10 collisions (in microscopic tracking).
        aval->SetCollisionSteps(10);
        aval->EnablePlotting(driftView);
        drift->EnablePlotting(driftView);
	 
        if(plotDriftLine) { PlotDriftLine(drift, driftView); app.Run(kTRUE);}
    }
  


    //------------------------------------------------------
    // define Histograms according to your needs
    int nBinsGain = 500;
    double gmin =   0.;
    double gmax = 10000.;
    TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons", nBinsGain, gmin, gmax);
    TH1F* hIons      = new TH1F("hIons",      "Number of ions",      nBinsGain, gmin, gmax);
    
    int nBinsChrg = 100;
    TH1F* hChrgE = new TH1F("hChrgE", "Electrons on plastic", nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);
    TH1F* hChrgI = new TH1F("hChrgI", "Ions on plastic", nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);
	TH1F* hVeloc = new TH1F("hVeloc", "velocity of electron", 100, gmin, 100);
	TH2F* hposi = new TH2F("hposi", "position of electron", 100, xmin, xmax, 100, ymin, ymax);
	//TH1F* hYposi = new TH1F("hYposi", "y position of electron", 100, ymin, ymax);
    TGraph* hElecZpos = new TGraph();
    TGraph* hIonZpos  = new TGraph();
    TGraph* hIonCath= new TGraph();
    TGraph* hIonTop= new TGraph();
    TGraph* hIonKap= new TGraph();
    TGraph* hIonBot= new TGraph();
    hElecZpos->SetName("hElecZpos");
    hIonZpos->SetName("hIonZpos");
    hIonCath->SetName("hIonCath");
    hIonTop->SetName("hIonTop");
    hIonKap->SetName("hIonKap");
    hIonBot->SetName("hIonBot");
    
    double sumIonsTotal = 0.;
    double sumIonsDrift = 0.;
    double sumIonsPlastic = 0.;
    
    double sumElectronsTotal = 0.;
    double sumElectronsPlastic = 0.;
    double sumElectronsUpperMetal = 0.;
    double sumElectronsLowerMetal = 0.;
    double sumElectronsTransfer = 0.;
    double sumElectronsOther = 0.;
    
    int nions = 0;
    int ntop  = 0;
    int nbot  = 0;
    int nkap  = 0;
    int ncath = 0;


	ViewSignal* signalView = new ViewSignal();
    //------------------------------------------------------
    // Start to do the simulation.
    //------------------------------------------------------
    // The goal of this simulation is to study the ion transparency through THGEM.
    // So the ions are put at induction field, so they can drift upwards.
    // In garfield I haven't figured out how to put ions sperately, so I put electrons.
    // When an electron is put, an ion is also generated in pires, that's how I got ions done.
    //------------------------------------------------------
    const int nEvents   = 100;
	clock_t start, stop;
    start = clock();
    for (int i = nEvents; i--;) {
        if (debug && i%100==0) std::cout << i << "/" << nEvents << "\n";
        // Randomize the initial position.zprim
        // no smear is needed here.
        //double smear = truepitch; 
		double smear = diam /2 *3;
        double x0 = - smear + 2*RndmUniform() * smear;
		//double x0 = 0;
        double y0 = 0.;
        double z0 = zprim;
        double t0 = 0.;
        double e0 = 0.1;
        cout << i << " 217" << endl;
	//aval->EnableDebugging();
        //x0 = x0 - smear + RndmUniform() * smear;
	aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);	// Calculate an electron avalanche
		int ne = 0, ni = 0;  //ne, ni represents the number of electrons and ions
        cout << "222" << endl;
	aval->GetAvalancheSize(ne, ni);
        //std::cout<<"Avalance "<<i<<" - electrons = "<<ne<<", ions = "<<ni<<"\n";
        hElectrons->Fill(ne);
        hIons->Fill(ni);
        const int np = aval->GetNumberOfElectronEndpoints();
        double xe1, ye1, ze1, te1, e1;
        double xe2, ye2, ze2, te2, e2;
        double xi1, yi1, zi1, ti1;
        double xi2, yi2, zi2, ti2;
        int status;


        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1,
                                      xe2, ye2, ze2, te2, e2, status);
            hElecZpos->SetPoint(int(sumElectronsTotal), xe2 * 1.e4, ze2 * 1.e4);
            sumElectronsTotal += 1.;
            if (ze2 > -kapton / 2. && ze2 < kapton / 2.) {
                hChrgE->Fill(ze2 * 1.e4);
                sumElectronsPlastic += 1.;
            } else if (ze2 >= kapton / 2. && ze2 <= kapton  / 2. + metal) {
                sumElectronsUpperMetal += 1.;
            } else if (ze2 <= -kapton / 2. && ze2 >= -kapton / 2. - metal) {
                sumElectronsLowerMetal += 1.;
            } else if (ze2 < -kapton / 2. - metal) {
                sumElectronsTransfer += 1.;
            } else {
                sumElectronsOther += 1.;
            }
            drift->DriftIon(xe1, ye1, ze1, te1);
            drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1,
                                  xi2, yi2, zi2, ti2, status);
            
            if (zi1 >= -zprim) {
                hIonZpos->SetPoint(nions++, xi2 * 1.e4, zi2 * 1.e4);
                sumIonsTotal += 1.;
                if (zi2 > -kapton / 2. && zi2 < kapton / 2.) {
                    hChrgI->Fill(zi2 * 1.e4);
                    sumIonsPlastic += 1.;
                    hIonKap->SetPoint(nkap++, xi2 * 1.e4, zi2 * 1.e4);
                } else if (zi2 > kapton/2.+metal) {
                    sumIonsDrift += 1.;
                    hIonCath->SetPoint(ncath++, xi2 * 1.e4, zi2 * 1.e4);
                } else if (zi2 <= -kapton/2. && zi2 >= -kapton/2.-metal) {
                    hIonTop->SetPoint(ntop++, xi2 * 1.e4, zi2 * 1.e4);
                } else if (zi2 >= kapton/2. && zi2 <= kapton/2.+metal) {
                    hIonBot->SetPoint(nbot++, xi2 * 1.e4, zi2 * 1.e4);
                }
            }
			//std::cout<<"end z is "<<ze2<<" timestart "<<te1<<" tend "<<te2<<" energy "<<e2<<std::endl;
			double velocity = sqrt(pow(ze2-ze1,2)+pow(ye2-ye1,2)+pow(xe2-xe1,2))/(te2-te1)*10000;
			hVeloc->Fill(velocity);
			if (ze2 <= -0.4){
			//std::cout<<xe2<<ye2<<std::endl;
			hposi->Fill(xe2, ye2, e2);
			}
			//std::cout<<velocity<<std::endl;
        }
		/*signalView->SetSensor(sensor);
		TCanvas* c1step = new TCanvas();
		signalView->SetCanvas(c1step);
		signalView->PlotSignal(label); 
		signalView->PlotSignal(label,false,false,true);
       // double signalinte = 0;
		for (int i =0; i<200; i++){
            //std::cout<<"time is "<<i<<"signal is "<<sensor->GetIonSignal(label,i)<<std::endl;
            signalinte += sensor->GetElectronSignal(label,i);
            }
        std::cout<<"-------------------- "<<signalinte<<std::endl;
		sensor->ClearSignal();
		std::cout<<"======="<<nEvents<<std::endl;
*/
    }
	
    cout << "295" << endl; 
    double fFeedback = 0.;
    if (sumIonsTotal > 0.) fFeedback = sumIonsDrift / sumIonsTotal;
    std::cout << "Fraction of ions drifting back: " << fFeedback << "\n";
    const double neMean = hElectrons->GetMean();
    const double niMean = hIons->GetMean();
    std::cout << "Mean number of electrons: " << neMean << "\n";
    std::cout << "Mean number of ions: " << niMean << "\n";
    std::cout << "Mean number of electrons on plastic: " << sumElectronsPlastic / nEvents << "\n";
    std::cout << "Mean number of ions on plastic: "      << sumIonsPlastic / nEvents << "\n";
    
    const double fUpperMetal = sumElectronsUpperMetal / sumElectronsTotal;
    const double fPlastic = sumElectronsPlastic / sumElectronsTotal;
    const double fLowerMetal = sumElectronsLowerMetal / sumElectronsTotal;
    const double fTransfer = sumElectronsTransfer / sumElectronsTotal;
    const double fOther = sumElectronsOther / sumElectronsTotal;
    std::cout << "Electron endpoints:\n";
    std::cout << "    upper metal: " << fUpperMetal * 100. << "%\n";
    std::cout << "    plastic:     " << fPlastic * 100. << "%\n";
    std::cout << "    lower metal: " << fLowerMetal * 100. << "%\n";
    std::cout << "    transfer:    " << fTransfer * 100. << "%\n";
    std::cout << "    other:       " << fOther * 100. << "%\n";

    TFile *f = new TFile(RootFile.Data(), "recreate");
    f->cd();
	
    hElectrons->Write();
    hIons->Write();
    hChrgE->Write();
    hChrgI->Write();
	hVeloc->Write();
	hposi->Write();
	//hYposi->Write();
    hElecZpos->Write();
    hIonZpos->Write();
    hIonCath->Write();
    hIonTop->Write();
    hIonKap->Write();
    hIonBot->Write();
    f->Close();
    
	stop = clock();
    cout<<" Time consuming: "<<(double)(stop - start)/CLOCKS_PER_SEC<<endl;
    TCanvas* cD = new TCanvas();
    if (plotDrift) {
        driftView->SetCanvas(cD);
        driftView->Plot();
        cD->SaveAs("Drift.eps");
        if (plotEFDrift) {
            ViewFEMesh* meshView = new ViewFEMesh();
            meshView->SetComponent(fm);
            meshView->SetPlane(0., -1, 0., 0., 0., 0.);
			meshView->SetArea(xmin, zmin, zmin, xmax, zmax, zmax);
			//std::cout<<"it come in "<<std::endl;
            meshView->SetFillColor(2,3); //matid=2 is FR4, painted as green.
			//std::cout<<"it come in 2 "<<std::endl;
            meshView->SetFillMesh(true);
			//std::cout<<"it come in 3 "<<std::endl;
            meshView->SetViewDrift(driftView);
			//std::cout<<"it come in 4"<<std::endl;
            TCanvas* cM = new TCanvas();
            meshView->SetCanvas(cM);
            meshView->Plot();

	    cM->SaveAs("EFMesh.eps");
        }
    }
    
    if (plotHistogram) {
        TCanvas* cH = new TCanvas("cH", "Histograms", 800, 700);
        cH->Divide(2, 2);
        cH->cd(1);
        hElecZpos->SetMarkerSize(0.5);
        hElecZpos->SetMarkerStyle(8);
        hElecZpos->Draw("ap");
        hElectrons->Draw();
        cH->cd(2);
        hIons->Draw();
        cH->cd(3);
        hChrgE->Draw();
        cH->cd(4);
        hChrgI->Draw();
	cH->SaveAs("Histogram.eps");
    }
    
    app.Run(kTRUE);
}




//------------------------------------------------------------------//
// Read ANSYS 12 output
//------------------------------------------------------------------//
void ReadFiledMap(TString path)
{
    std::cout<<" Entering "<<path<<std::endl;
    // Load the field map.
	fm = new ComponentComsol();
    TString mfile = path+TString("/dat/mesh.mphtxt");
    TString dfile = path+TString("/dat/dielectrics_article.dat");
    TString ffile = path+TString("/dat/Field.txt");

    fm->Initialise(mfile.Data(), dfile.Data(), ffile.Data());
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
	fm->SetMagneticField(0.,0.,0.);
    fm->PrintRange();				// print some information about the cell dimention
}



//------------------------------------------------------------------//
// define gas compounents
//------------------------------------------------------------------//
void DefineGas(const char *gas1, double rat1, const char *gas2, double rat2)
{
    gas = new MediumMagboltz();
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
    // Load the ion mobilities.
    gas->LoadIonMobility("/work/src/programs/garfield/garfield/Data/IonMobility_Ar+_Ar.txt");
    


    // Associate the gas with the corresponding field map material.
    const int nMaterials = fm->GetNumberOfMaterials();
    for (int i = 0; i < nMaterials; ++i) {
        const double eps = fm->GetPermittivity(i);	 	// Get permittivity (episilon)
        if (fabs(eps - 1.) < 1.e-3) fm->SetMedium(i, gas);	// set gas as medium if eps is arround 1
    }
    fm->PrintMaterials();
}


//------------------------------------------------------------------//
// 画电场，电力线等图
//------------------------------------------------------------------//
void PlotField(ComponentComsol* fm)
{
    
    ViewField* fieldView = new ViewField();
    fieldView->SetComponent(fm);
    
    fieldView->SetNumberOfSamples1d(2000); 	      //default is 1000, the density of the plotting grid
    fieldView->SetNumberOfSamples2d(100,100);	      //default is 200, 200. 二维图的nX, nY
    fieldView->SetNumberOfContours(50); 	      //default is 100, 等高线的次数
    //fieldView->SetVoltageRange(-3000., 3000.);    //the min/max of potential in the Component is used as default, but it can also be defined by user.
    //fieldView->SetElectricFieldRange(0,20000);	  //default is (0, 100). it must be defined by user for electric field.
    //fieldView->SetWeightingFieldRange(0,10);	  //it must be defined by user for weighting field.

     //fieldView->SetPlane(0.866, -0.5, 0., 0., 0., 0.);	      //(fx,fy,fz,x0,y0,z0), (fx,fy,fz) 指定的是平面的法线方向。
/*	fieldView->SetPlane(0, -1, 0, 0., 0., 0.);	
     //fieldView->Rotate(TMath::Pi()*3./2.);                     //指定切出来的平面，在做图的时候选择的角度。比如有时候出来的是Y:Z(Y是横轴，Z是纵轴），可以选择Pi/2来交换成Z:Y
   fieldView->SetArea(-pitch / 2., zmin, pitch / 2., zmax);  //指定这个平面所画的区域大小
    TCanvas* cF1 = new TCanvas();
    fieldView->SetCanvas(cF1);
    fieldView->PlotContour("v");                            // v/p/phi/volt/voltage/pot/potential 都是画电势, e/field 画电场, or ex/ey/ez         
*/
    fieldView->SetPlane(0, 0, 0, 0., 0., 0.2);
    cout << "OK" << endl;
    fieldView->SetArea(-pitch, -pitch, pitch, pitch);  //指定这个平面所画的区域大小
    TCanvas* cF2 = new TCanvas();
    fieldView->SetCanvas(cF2);
    fieldView->PlotContour("e"); 
    cF2->SaveAs("Field2.eps");

    TCanvas* cF3 = new TCanvas();
    fieldView->SetCanvas(cF3);
    fieldView->PlotProfile(0, 0, -1*(kapton/2+metal+induct), 0, 0, kapton/2+metal+drift, "v");
    cF3->SaveAs("Field3.eps");
    cout << "480" << endl; 
/*
    TCanvas* cF4 = new TCanvas();
    fieldView->SetCanvas(cF4);
    fieldView->PlotProfile(0, 0, -1*(kapton/2+metal+induct), 0, 0, kapton/2+metal+drift, "ex");

    TCanvas* cF5 = new TCanvas();
    fieldView->SetCanvas(cF5);
    fieldView->PlotProfile(0, 0, -1*(kapton/2+metal+induct), 0, 0, kapton/2+metal+drift, "ey");
   
    TCanvas* cF6 = new TCanvas();
    fieldView->SetCanvas(cF6);
    fieldView->PlotProfile(0, 0, -1*(kapton/2+metal+induct), 0, 0, kapton/2+metal+drift, "ez");

    */
    return;
    
}


//------------------------------------------------------------------//
// Plot drift lines
//------------------------------------------------------------------//
void PlotDriftLine(AvalancheMC* drift, ViewDrift* driftView)
{ 
    drift->DisableDiffusion();
    const int Nbin=50;
    for(int i=0;i<=Nbin;i++) {  
        double xtmp = xmin+i*(xmax-xmin)/Nbin;
        drift->DriftIon(xtmp,0,zmin,0);
        //   drift->DriftIon(xtmp,0,-0.0120,0);
    }
    
    TCanvas* cD = new TCanvas();
    driftView->SetCanvas(cD);
    driftView->Plot();
     
    ViewFEMesh* meshView = new ViewFEMesh();
    meshView->SetComponent(fm);
    meshView->SetPlane(0., -1, 0., 0., 0., 0.);
    meshView->SetArea(xmin, ymin, zmin, xmax, ymax, zmax);
    meshView->SetFillColor(2,3); //matid=2 is FR4, painted as green.
    meshView->SetFillMesh(true);  
    meshView->SetViewDrift(driftView);
    TCanvas* cM = new TCanvas();
    meshView->SetCanvas(cM); 
    meshView->Plot();

    cD->SaveAs("DriftLine.eps");
    cM->SaveAs("Mesh.eps");
}


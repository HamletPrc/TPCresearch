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

void ReadFiledMap(TString);
void DefineGas(const char *gas1, double rat1, const char *gas2, double rat2, const char *filename);
void PlotField(ComponentComsol* fm);
void PlotDriftLine(AvalancheMC* drift, ViewDrift* driftView);

ComponentComsol *fm;
MediumMagboltz    *gas;

// Dimensions of the GEM, in !!!! CM !!!!, note that in Ansys it's in MM !
const double pitch  = 0.1;
const double kapton = 0.04;
const double metal  = 0.001;
const double diam   = 0.04;
const double induct = 0.2;
const double drift  = 0.2;

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
    const bool plotField = true;
    const bool plotDrift = true; 
    const bool plotDriftLine = false;
    const bool plotEFDrift = true;
    const bool plotHistogram = false;
    const bool plotGain = true;
    const bool plotIon = true;

    //the prim ionization Z pozition. Now I set it to induction field to study IONs!!
    const double zprim  = zmax-0.01;
    
    const double tStart =0;
    const double tStop =100;
    const int nSteps =100;
    const double tStep =(tStop - tStart) /nSteps;
    const std::string label = "readout";
    //-----------------
    // Start simulation
    ReadFiledMap(TString(Path));
    
    DefineGas("argon", 90, "iC4H10", 10, "Ar_90_iC4H10_10.gas");
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
    double gmax = 5000.;
    TH1F* hElectrons = new TH1F("hElectrons", "Number of electrons", nBinsGain, gmin, gmax);
    TH1F* hIons      = new TH1F("hIons",      "Number of ions",      nBinsGain, gmin, gmax);
    TH1F* hGainAbs = new TH1F("hGainAbs", "Absolute gain", nBinsGain, gmin, gmax);
    TH1F* hGainEff = new TH1F("hGainEff", "Effective gain", nBinsGain, gmin, gmax);
    
    int nBinsChrg = 100;
    TH1F* hChrgE = new TH1F("hChrgE", "Electrons on plastic", nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);
    TH1F* hChrgI = new TH1F("hChrgI", "Ions on plastic", nBinsChrg, -0.5e4 * kapton, 0.5e4 * kapton);
    TH1F* hVeloc = new TH1F("hVeloc", "velocity of electron", 100, gmin, 100);
    TH1F* hIonInitX = new TH1F("hIonInitX", "Production x-coordinate of ions", 40, xmin, xmax);
    TH1F* hIonInitY = new TH1F("hIonInitY", "Production y-coordinate of ions", 40, ymin, ymax);
    TH1F* hIonInitR = new TH1F("hIonInitR", "Production r-coordinate of ions", 40, 0, diam);
    TH1F* hIonInitPhi = new TH1F("hIonInitPhi", "Production phi-coordinate of ions", 40, -3.14,3.14);
    TH1F* hIonInitZ = new TH1F("hIonInitZ", "Production z-coordinate of ions", 40, -0.05,kapton+0.05);
    TH2F* hposi = new TH2F("hposi", "position of electron", 100, xmin, xmax, 100, ymin, ymax);
    TH2F* hIonProduce = new TH2F("hIonProduce", "Production positions of ions", 100, xmin, xmax, 100, zmin, zmax);
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
    // wmz wants more information on detection effeciency and effective gain, so define these variables
    bool findElectron=false;
    double singleElectronGain=0.;
    double singleElectronAnode=0.;
    double nDetects=0.;
    //
    
    int nions = 0;
    int ntop  = 0;
    int nbot  = 0;
    int nkap  = 0;
    int ncath = 0;
    int nOther = 0.;


    double xe1, ye1, ze1, te1, e1;
    double xe2, ye2, ze2, te2, e2;
    double xi1, yi1, zi1, ti1,ri1,phi1;
    double xi2, yi2, zi2, ti2;
    int status;


    TFile *data = new TFile("result0.root","recreate");
    TTree *tdata = new TTree("tdata","TTree for result");
    tdata->Branch("xe1",&xe1);
    tdata->Branch("ye1",&ye1);
    tdata->Branch("ze1",&ze1);
    tdata->Branch("xe2",&xe2);
    tdata->Branch("ye2",&ye2);
    tdata->Branch("ze2",&ze2);
    tdata->Branch("xi1",&xi1);
    tdata->Branch("yi1",&yi1);
    tdata->Branch("zi1",&zi1);
    tdata->Branch("xi2",&xi2);
    tdata->Branch("yi2",&yi2);
    tdata->Branch("zi2",&zi2);

    //------------------------------------------------------
    // Start to do the simulation.
    //------------------------------------------------------
    // The goal of this simulation is to study the ion transparency through THGEM.
    // So the ions are put at induction field, so they can drift upwards.
    // In garfield I haven't figured out how to put ions sperately, so I put electrons.
    // When an electron is put, an ion is also generated in pires, that's how I got ions done.
    //------------------------------------------------------
    const int nEvents   = 1;
    clock_t start, stop;
    start = clock();
    for (int i = nEvents; i--;) {
        if (debug && i%10==0) std::cout << i << "/" << nEvents << "\n";
        // Randomize the initial position.zprim
        // no smear is needed here
		//double smear = diam /2 *3;
        //double x0 = - smear + 2*RndmUniform() * smear;
		double x0 = 0;
        double y0 = 0.;
        double z0 = zprim;
        double t0 = 0.;
        double e0 = 0.1;
	    //aval->EnableDebugging();
		int ne = 0, ni = 0;  //ne, ni represents the number of electrons and ions

        aval->AvalancheElectron(x0, y0, z0, t0, e0, 0., 0., 0.);	// Calculate an electron avalanche
	    aval->GetAvalancheSize(ne, ni);

        //std::cout<<"Avalance "<<i<<" - electrons = "<<ne<<", ions = "<<ni<<"\n";
        hElectrons->Fill(ne);
        hIons->Fill(ni);
        const int np = aval->GetNumberOfElectronEndpoints();

        findElectron=false;
        singleElectronGain=np;
        singleElectronAnode=0;
        hGainAbs->Fill(singleElectronGain);

        for (int j = np; j--;) {
            aval->GetElectronEndpoint(j, xe1, ye1, ze1, te1, e1,
                                    xe2, ye2, ze2, te2, e2, status);
            hElecZpos->SetPoint(int(sumElectronsTotal), xe2 * 1.e4, ze2 * 1.e4);
            sumElectronsTotal += 1.;
            if (ze2 > metal && ze2 < metal + kapton) {
                hChrgE->Fill(ze2 * 1.e4);
                sumElectronsPlastic += 1.;
            } else if (ze2 >= kapton + metal && ze2 <= kapton + 2*metal) {
                sumElectronsUpperMetal += 1.;
            } else if (ze2 >= 0 && ze2 <= metal) {
                sumElectronsLowerMetal += 1.;
            } else if (ze2 < 0) {
                sumElectronsTransfer += 1.;
                findElectron=true;
                singleElectronAnode += 1.;
            } else {
                sumElectronsOther += 1.;
            } 

	        drift->DriftIon(xe1, ye1, ze1, te1);
            drift->GetIonEndpoint(0, xi1, yi1, zi1, ti1,
                                  xi2, yi2, zi2, ti2, status);
            hIonZpos->SetPoint(nions++, xi2 * 1.e4, zi2 * 1.e4);
            sumIonsTotal += 1.;
            tdata->Fill();
                
            if (zi2 > metal && zi2 < kapton + metal) {
                hChrgI->Fill(zi2 * 1.e4);
                sumIonsPlastic += 1.;
                hIonKap->SetPoint(nkap++, xi2 * 1.e4, zi2 * 1.e4);
            } else if (zi2 > kapton + 2*metal) {
                sumIonsDrift += 1.;
                hIonCath->SetPoint(ncath++, xi2 * 1.e4, zi2 * 1.e4);
            } else if (zi2 >= kapton +metal && zi2 >= kapton + 2*metal) {
                    hIonTop->SetPoint(ntop++, xi2 * 1.e4, zi2 * 1.e4);
            } else if (zi2 >= 0  && zi2 <= metal) {
                    hIonBot->SetPoint(nbot++, xi2 * 1.e4, zi2 * 1.e4);
            } else nOther += 1.;
            
            // The production information of ions
            hIonProduce->Fill(xi1,zi1);
            hIonInitX->Fill(xi1);
            hIonInitY->Fill(yi1);
            hIonInitZ->Fill(zi1);
            ri1 = sqrt(xi1*xi1+yi1*yi1);  hIonInitR->Fill(ri1);
            phi1 = asin(xi1/ri1); hIonInitPhi->Fill(phi1);

			//std::cout<<"end z is "<<ze2<<" timestart "<<te1<<" tend "<<te2<<" energy "<<e2<<std::endl;
			double velocity = sqrt(pow(ze2-ze1,2)+pow(ye2-ye1,2)+pow(xe2-xe1,2))/(te2-te1)*10000;
			hVeloc->Fill(velocity);
			if (ze2 <= -0.4){
			//std::cout<<xe2<<ye2<<std::endl;
			hposi->Fill(xe2, ye2, e2);
			}
			//std::cout<<velocity<<std::endl;
        }   // end of the "for" loop of one event!

        if(findElectron) {
            nDetects += 1.;
            hGainEff->Fill(singleElectronAnode);
        }
    
    } // end of the "for" loop of all events!
    	
    cout << "295" << endl;
    double fDetect = 0.;
    double fFeedback = 0.;
    if(sumElectronsTotal>0.) fDetect = nDetects / nEvents *100;
    if (sumIonsTotal > 0.) fFeedback = sumIonsDrift / sumIonsTotal;
    std::cout << "Number of detection: " << nDetects << endl;
    std::cout << "Fraction of detection of events: " << fDetect << "%" << endl;
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

    const double IonUpperMetal  = ntop / sumIonsTotal;
    const double IonPlastic =  sumIonsPlastic / sumIonsTotal;
    const double IonLowerMetal =  nbot / sumIonsTotal;
    const double IonOther = nOther /sumIonsTotal;
    std::cout << "Ion endpoints:\n";
    std::cout << "    upper metal: " << IonUpperMetal * 100. << "%\n";
    std::cout << "    plastic:     " << IonPlastic * 100. << "%\n";
    std::cout << "    lower metal: " << IonLowerMetal * 100. << "%\n";
    std::cout << "    other:       " << IonOther * 100. << "%\n";


    tdata->Write();
    data->Close();
    

    TFile *f = new TFile(RootFile.Data(), "recreate");
    f->cd();
	
    hElectrons->Write();
    hIons->Write();
    hChrgE->Write();
    hChrgI->Write();
	hVeloc->Write();
	hposi->Write();
    hIonProduce->Write();
    hIonInitX->Write();
    hIonInitY->Write();
    hIonInitZ->Write();
    hIonInitR->Write();
    hIonInitPhi->Write();
	//hYposi->Write();
    hElecZpos->Write();
    hIonZpos->Write();
    hIonCath->Write();
    hIonTop->Write();
    hIonKap->Write();
    hIonBot->Write();
    hGainAbs->Write();
    hGainEff->Write();
    f->Close();
    
	stop = clock();
    cout<<" Time consuming: "<<(double)(stop - start)/CLOCKS_PER_SEC<<endl;
    TCanvas* cD = new TCanvas();
    if (plotDrift) {
        driftView->SetCanvas(cD);
        driftView->Plot();
        cD->SaveAs("output/xyDiffusion.eps");
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

	    cM->SaveAs("output/DriftLine.eps");
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
	cH->SaveAs("output/Histogram.eps");
    }
    
    if(plotGain) {
        TCanvas* cG = new TCanvas("cG", "Gains", 800, 700);
        cG->Divide(2,1);
        cG->cd(1);
        hGainAbs->Draw();
        cG->cd(2);
        hGainEff->Draw();
		cG->SaveAs("output/Gain.eps");
    }

    if(plotIon){
        TCanvas* cI1 = new TCanvas("cI1", "Production positions of ions", 800, 700);
        cI1->Divide(1,1);
        hIonProduce->Draw();
        TCanvas* cI2 = new TCanvas("cI2", "Production positions of ions", 800, 700);
        cI2->Divide(3,1);
        cI2->cd(1);
        hIonInitX->Draw();
        hIonInitX->SetTitle("Ion x-coordinate");
        hIonInitX->GetTitle();
        cI2->cd(2);
        hIonInitY->Draw();
        cI2->cd(3);
        hIonInitZ->Draw();
        TCanvas* cI3 =  new TCanvas("cI3", "Production positions of ions", 800, 700);
        cI3->Divide(3,1);
        cI3->cd(1);
        hIonInitR->Draw();
        cI3->cd(2);
        hIonInitPhi->Draw();
        cI3->cd(3);
        hIonInitZ->Draw();
		cI1->SaveAs("output/Ion1.eps");
		cI2->SaveAs("output/Ion2.eps");
		cI3->SaveAs("output/Ion3.eps");
    }
    
    app.Run(kTRUE);
}




//------------------------------------------------------------------//
// Read COMSOL5.4 output
//------------------------------------------------------------------//
void ReadFiledMap(TString path)
{
    std::cout<<" Entering "<<path<<std::endl;
    // Load the field map.
	fm = new ComponentComsol();
    TString mfile = path+TString("/dat/meshNew.mphtxt");
    TString dfile = path+TString("/dat/dielectrics_article.dat");
    TString ffile = path+TString("/dat/FieldNew.txt");

    fm->Initialise(mfile.Data(), dfile.Data(), ffile.Data());
    fm->EnableMirrorPeriodicityX();
    fm->EnableMirrorPeriodicityY();
	fm->SetMagneticField(0.,0.,0.);
    fm->PrintRange();				// print some information about the cell dimention
}



//------------------------------------------------------------------//
// define gas compounents
//------------------------------------------------------------------//
void DefineGas(const char *gas1, double rat1, const char *gas2, double rat2, const char *filename)
{

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
    // Load the ion mobilities
    gas->LoadIonMobility("/work/src/programs/garfield/garfield/Data/IonMobility_Ar+_Ar.txt");	
	
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

    
    cout << "OK" << endl;
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
    cF2->SaveAs("output/xyField.eps");

    TCanvas* cF3 = new TCanvas();
    fieldView->SetCanvas(cF3);
    fieldView->PlotProfile(0, 0, -1*(kapton/2+metal+induct), 0, 0, kapton/2+metal+drift, "v");
    cF3->SaveAs("output/FieldConverge.eps");
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

    TCanvas* cF7 = new TCanvas();
    double z_print = -0.001;
    fieldView->SetCanvas(cF7);
    fieldView->PlotProfile(-diam*2, 0, z_print, diam*2, 0, z_print, "ez");
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

    cD->SaveAs("output/Drift.eps");
    cM->SaveAs("output/FieldLine.eps");
}


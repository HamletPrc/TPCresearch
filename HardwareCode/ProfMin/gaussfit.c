void gaussfit(){
    const int bins = 50;
    const double xmin = -1.85;
    const double xmax = -1.35;

    TFile* file=new TFile("result_root/0V.root");
    TTree* t1=(TTree*)file->Get("calibconst");
    TH1F* currenthist=new TH1F("currenthist","The current at 0V", bins, xmin, xmax);
    
    double c;
    t1->SetBranchAddress("current", &c);
    const int entries=t1->GetEntries();

    for (int i=0; i<=entries; ++i) {
        t1->GetEntry(i);
        currenthist->Fill(c);
    }
    currenthist->Scale(bins/( (xmax-xmin)* currenthist->Integral()) );

    TF1 *fit=new TF1("fit","1/(sqrt(2*pi)*[1])*exp(-(x-[0])*(x-[0])/(2*[1]*[1]))", xmin, xmax);
    fit->SetParameters(-1.6,0.04);
    currenthist->Fit(fit);
	currenthist->Draw();
	
}

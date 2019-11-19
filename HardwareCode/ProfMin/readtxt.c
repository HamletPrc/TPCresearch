void readtxt(){
	TString dataFile = Form("data_used/3000V_stable.txt");
	TTree *datatree = new TTree("calibconst","calibconst");
	datatree->ReadFile(dataFile,"current/D:Temperature/D:Wet/D");
	datatree->SaveAs("root_result/3000V.root");
}

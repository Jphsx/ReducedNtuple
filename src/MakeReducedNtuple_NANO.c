#ifndef MAKEREDUCE
#define MAKEREDUCE
// C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <ostream>
#include <istream>
#include <stdio.h>
#include <dirent.h>
#include <vector>

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TList.h>
#include <TLeafElement.h>
#include <TLorentzVector.h>

//#include "ReducedNtuple.hh"
//#include "SUSYNANOBase.hh"
#include "prod2018MC_reducedNANO_Muon.h"

//MET classes
#include "prod2018MC_reducedNANO_MET.h"
#include "prod2017MC_reducedNANO_MET.h"
#include "prod2016MC_reducedNANO_MET.h"

//trigger classes
#include "prod2018MC_reducedNANO_Triggers.h"
#include "prod2017MC_reducedNANO_Triggers.h"
#include "prod2016MC_reducedNANO_Triggers.h"

using namespace std;
using std::vector;

//template function to deal with various tselector
template<class selectortype>
void produceReducedTree(selectortype& selector, std::string ofilename){
	//copy branches to output file	
	auto ofile =  new TFile(ofilename.c_str(), "RECREATE");
//	auto reducedTree = selector.fReader.GetTree()->CloneTree();
	auto reducedTree = selector.fChain->CloneTree();
	//auto reducedTree = selector.fChain->CloneTree();
	//reducedTree->Write();
	//reducedTree->CopyEntries(selector.fReader.GetTree());
	reducedTree->Write();
	ofile->Write();
	ofile->Close();
}

/// Main function that runs the analysis algorithm on the
/// specified input files
int main(int argc, char* argv[]) {

  /// Gets the list of input files and chains
  /// them into a single TChain
  char inputFileName[400];
  char inputListName[400];
  char inputFolderName[400];
  char outputFileName[400];
  char TreeName[400];
  char DataSet[400];
  char FileTag[400];
  char EventCount[400];

  char SelectorClassName[400];

  bool DO_FILE = false;
  bool DO_LIST = false;
  bool DO_FOLDER = false;
  bool DO_TREE = false;
  bool DO_SMS = false;

  if ( argc < 4 ){
    cout << "Error at Input: please specify an input file name, a list of input ROOT files and/or a folder path"; 
    cout << " , an output filename, and a selector class name:" << endl; 
    cout << "  Example:      ./MakeReducedNtuple_NANO.x -ifile=input.root -ofile=output.root -dataset=dataset_name -filetag=sample_tag"  << endl;
    cout << "  Example:      ./MakeReducedNtuple_NANO.x -ilist=input.list -ofile=output.root -dataset=dataset_name -filetag=sample_tag"  << endl;
    cout << "  Example:      ./MakeReducedNtuple_NANO.x -ifold=folder_path -ofile=output.root -dataset=dataset_name -filetag=sample_tag -tree=treename -eventcount=event_count --sms" << endl;
    cout << " additional tags for object based reduced tree: -selector=TSelector_ClassName "<<endl; 
    return 1;
  }
  for (int i=0;i<argc;i++){
    if (strncmp(argv[i],"-ifile",6)==0){
      sscanf(argv[i],"-ifile=%s",  inputFileName);
      DO_FILE = true;
    }
    if (strncmp(argv[i],"-ilist",6)==0){
      sscanf(argv[i],"-ilist=%s",  inputListName);
      DO_LIST = true;
    }
    if (strncmp(argv[i],"-ifold",6)==0){
      sscanf(argv[i],"-ifold=%s",  inputFolderName);
      DO_FOLDER = true;
    }
    if (strncmp(argv[i],"-tree",5)==0){
      sscanf(argv[i],"-tree=%s",  TreeName);
      DO_TREE = true;
    } 
    if (strncmp(argv[i],"-selector",9)==0){
      sscanf(argv[i],"-selector=%s", SelectorClassName); 
    }
    if (strncmp(argv[i],"-ofile",6)==0) sscanf(argv[i],"-ofile=%s", outputFileName);
    if (strncmp(argv[i],"-dataset",8)==0)   sscanf(argv[i],"-dataset=%s", DataSet);
    if (strncmp(argv[i],"-filetag",8)==0)   sscanf(argv[i],"-filetag=%s", FileTag);
    if (strncmp(argv[i],"-eventcount",11)==0)   sscanf(argv[i],"-eventcount=%s", EventCount);
    if (strncmp(argv[i],"--sms",5)==0)  DO_SMS = true;
  }


  gROOT->ProcessLine("#include <vector>");

  vector<string> filenames;

  char Buffer[500];
  char MyRootFile[2000];  

  if(DO_FOLDER){
    DIR *dpdf;
    struct dirent *epdf;
    dpdf = opendir(inputFolderName);
    if(dpdf == NULL){
      cout << "ERROR: " << inputFolderName << " is not a directory" << endl;
      return 1;
    }
    string dname(inputFolderName);
    while ((epdf = readdir(dpdf))){
      if(string(epdf->d_name).find(".root") == string::npos)
	continue;
      string name = dname+"/"+string(epdf->d_name);
      filenames.push_back(name);
    }
  }

  if(DO_LIST){
    ifstream *inputFile = new ifstream(inputListName);
    while( !(inputFile->eof()) ){
      inputFile->getline(Buffer,500);
      if (!strstr(Buffer,"#") && !(strspn(Buffer," ") == strlen(Buffer))){
	sscanf(Buffer,"%s",MyRootFile);
	filenames.push_back(MyRootFile);
      }
    }
    inputFile->close();
    delete inputFile;
  }

  if(DO_FILE){
    filenames.push_back(inputFileName);
  }

  TChain* chain;
  if(DO_TREE)
    chain = (TChain*) new TChain(TreeName);
  else
    chain = (TChain*) new TChain("Events");
  
  int Nfile = filenames.size();
  for(int i = 0; i < Nfile; i++){
    chain->Add(filenames[i].c_str());
    cout << "   Adding file " << filenames[i] << endl;
  }


  //use std::string because who would want to use char array
  std::string _selectorClassName(SelectorClassName);
  std::string _ofilename(outputFileName);
  //create appropriate selector class
  if(_selectorClassName.compare("prod2018MC_reducedNANO_Muon") == 0){
	std::cout<<"Using selector: "<< _selectorClassName <<std::endl;
	prod2018MC_reducedNANO_Muon s(chain);
	produceReducedTree(s,_ofilename);	
  }

//MET CLASSES
  if(_selectorClassName.compare("prod2018MC_reducedNANO_MET") == 0){
	std::cout<<"Using selector: "<< _selectorClassName <<std::endl;
	prod2018MC_reducedNANO_MET s(chain);
	produceReducedTree(s,_ofilename);
  }
  if(_selectorClassName.compare("prod2017MC_reducedNANO_MET") == 0){
	std::cout<<"Using selector: "<< _selectorClassName <<std::endl;
	prod2017MC_reducedNANO_MET s(chain);
	produceReducedTree(s,_ofilename);
  }
  if(_selectorClassName.compare("prod2016MC_reducedNANO_MET") == 0){
	std::cout<<"Using selector: "<< _selectorClassName <<std::endl;
	prod2016MC_reducedNANO_MET s(chain);
	produceReducedTree(s,_ofilename);
  }
// Trigger classes
  if(_selectorClassName.compare("prod2018MC_reducedNANO_Triggers") == 0){
	std::cout<<"Using selector: "<< _selectorClassName <<std::endl;
	prod2018MC_reducedNANO_Triggers s(chain);
	produceReducedTree(s,_ofilename);
  }
  if(_selectorClassName.compare("prod2017MC_reducedNANO_Triggers") == 0){
	std::cout<<"Using selector: "<< _selectorClassName <<std::endl;
	prod2017MC_reducedNANO_Triggers s(chain);
	produceReducedTree(s,_ofilename);
  }
  if(_selectorClassName.compare("prod2016MC_reducedNANO_Triggers") == 0){
	std::cout<<"Using selector: "<< _selectorClassName <<std::endl;
	prod2016MC_reducedNANO_Triggers s(chain);
	produceReducedTree(s,_ofilename);
  }


 // ReducedNtuple<SUSYNANOBase>* ntuple = new ReducedNtuple<SUSYNANOBase>(chain);

 // ntuple->AddLabels(string(DataSet),string(FileTag));
 // ntuple->AddEventCountFile(string(EventCount));

  //if(DO_SMS)
   // ntuple->DoSMS();

 // ntuple->WriteNtuple(string(outputFileName));

  //delete ntuple;
 
  return 0;

}
#endif

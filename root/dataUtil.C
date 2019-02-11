

TChain *getFriendedFromDir(int zip=4,string dir="YBe",string selection="")
{
  ostringstream mergezip,calibzip,event;
  mergezip << "rqDir/zip" << zip;
  calibzip << "rrqDir/calibzip" << zip;
  event << "rqDir/eventTree";
  cout << mergezip.str() << "\t" << calibzip.str() << "\t" << event.str() << endl;
  TChain *merged = getChainFromSeries(dir+"/merge"+selection,mergezip.str());
  TChain *calib = getChainFromSeries(dir+"/calib"+selection,calibzip.str());
  TChain *eventtree = getChainFromSeries(dir+"/merge"+selection,event.str());

  eventtree->AddFriend(merged);
  eventtree->AddFriend(calib);

  return eventtree;

}
TChain *getChainFromSeries(string series="07140524_0040",string tree="rqDir/zip1",TChain *old=NULL)
{
  //parse out everything earlier than the last /
  string dir="./";
  size_t found = series.rfind("/");
  if(found != string::npos){
    dir = series.substr(0,found);
    series = series.substr(found+1,series.size()-found);
  }

  //get the files
  vector<string> files;
  getRQFiles(files,series,dir);

  //loop through the files and add them to a TChain
  TChain *ch;
  if(!old)
    ch = new TChain(Form("%s",tree.c_str()));
  else
    ch = old;
  for(int i=0;i<files.size();i++){
    ch->Add(files[i].c_str());
  }

  return ch;
}
TChain *getChainFromList(string list="cf_biasp",string tree="rqDir/zip1")
{
  //make sure we have the env var
  char *dd = getenv("CDMSBATSFILELISTS");
  string fileListDir = string(dd);

  //open the file list
  string filename = fileListDir+"/"+list;
  ifstream seriesfile(filename.c_str(),ios::in);
  if(!seriesfile){
    cout << "ERROR!(getChainFromList): list file cannot be opened.." << endl;
    return NULL;
  }
  cout << filename << endl;

  //loop through the list and gather the series names
  string thisseries;
  seriesfile >> thisseries;
  TChain *ch=NULL;
  while(!seriesfile.eof()){
    cout << thisseries << endl;
    ch = getChainFromSeries(thisseries,tree,ch);
    seriesfile >> thisseries;
  }

  return ch;
}
TChain *getChain(string seriesorlist="cf_biasp",string tree="rqDir/zip1")
{
  //get the list var
  char *dd = getenv("CDMSBATSFILELISTS");
  string fileListDir = string(dd);

  //initialize the returned pointer
  TChain *ch = NULL;

  //open the file list
  string filename = fileListDir+"/"+seriesorlist;
  ifstream seriesfile(filename.c_str(),ios::in);
  if(!seriesfile){
    //assume it's a series
    seriesfile.close();
    ch = getChainFromSeries(seriesorlist,tree,ch);
  }
  else{
    //assume it's a list
    seriesfile.close();
    ch = getChainFromList(seriesorlist,tree);
  }
  
  return ch;
}

#include <google/sparse_hash_set>
#include <iostream>
#include <fstream>
#include <tr1/functional>
#include <map>
#include <sstream>
#include <stdlib.h> 
#include <math.h>
#include <gsl/gsl_randist.h>
#include <boost/lexical_cast.hpp>
#include <iomanip>
using namespace std;

class commonUtil{
    
    char delimiter ;
public:
    commonUtil(char seperator){
	delimiter=seperator;
    }
    commonUtil(){
	delimiter=',';
    }
    vector<string> digestLine(string line) {
	vector<string> elems;
	string token="";
	size_t i=0;
	while(i<line.length()){
	    if (line[i]==delimiter){
		elems.push_back(token);
		token="";
	    }
	    else
		token=token+line[i];
	    i++;
	}
	elems.push_back(token);
	return elems;
    }
    string ConvertToString (float number){
	std::ostringstream buff;
	buff<< fixed <<setprecision(2)<<number;
	return buff.str();   
    }
    double strToDouble(string str){
	double value;
	try
	{
	    value = boost::lexical_cast<double>(str);
	}
	catch (boost::bad_lexical_cast const&)
	{
	    value = 0;
	}
	return value;
    }
    float strToFloat(string str){
	return (atof(str.c_str()));
    }
    bool is_digits(const char &c)
    {
	if (c=='0'||c=='1'||c=='2'||c=='3'||c=='4'||c=='5'||c=='6'||c=='7'||c=='8'||c=='9'||c=='.' )
	    return true;
	return false;
    }
    bool cheakKmerHelth(string &kmer, size_t kmerSize){
	if (kmer.length()!=kmerSize)
	    return false;
	for (size_t i=0;i<kmerSize;i++){
	    if( (kmer[i]!='A'&& kmer[i]!='T'&&kmer[i]!='C'&&kmer[i]!='G'&&kmer[i]!='a'&&kmer[i]!='t'&&kmer[i]!='c'&&kmer[i]!='g'))
		return false;
	}
	for (size_t i=0;i<kmerSize;i++){
	    if(kmer[i]=='a')
		kmer[i]='A';
	    if (kmer[i]=='t')
		kmer[i]='T';
	    if (kmer[i]=='c')
		kmer[i]='C';
	    if (kmer[i]=='g')
		kmer[i]='G';
	}
	
	return true;
    }
    string reverseComplement(string kmer){
	string result="";
	for (int i=kmer.length()-1;i>=0;i--){
	    switch(kmer[i]){
		case 'A':
		    result=result+'T';
		    break;
		case 'T':
		    result=result+'A';
		    break;
		case 'C':
		    result=result+'G';
		    break;
		case 'G':
		    result=result+'C';
		    break;
	    }
	}
	return result;
    }
};
class jellyfishUtil{
    typedef google::sparse_hash_set<std::string> HashSet;
    typedef HashSet::const_iterator HashSetIt;
    typedef map< int , int > freqVec;
    typedef freqVec::const_iterator freqVecIt; 
    commonUtil util; 
    size_t kmerSize;
    size_t numThreads;
    string Mem;
    string FileName;
    char delimiter ;
    bool fillSet(HashSet &Set, string fileName)
    {  
	ifstream  readFastaFileStream;
	readFastaFileStream.open(fileName.c_str());
	if(!readFastaFileStream.good()){
	    std::cerr << "Error opening file"<< fileName << std::endl;
	    return false;
	}
	string line;
	while (!readFastaFileStream.eof())
	{
	    readFastaFileStream>>line;
	    if (line[0] != '>' ){
		string kmer =line;                                     //line.substr (0,kmerSize-1);
		if (!util.cheakKmerHelth(kmer,kmerSize ))
		    continue;
		Set.insert(kmer);
	    }
	    
	}
	//Set.insert( reverseComplement("gtctggatgactggacccacaacaaaagcat"));
	cout<<Set.size()<<endl;
	readFastaFileStream.close();
	return true;
	
    }
    void generateKmerFile(string inputFileName, string outputFileName){
	string tempFile="tempKmerFreq";
	string command="jellyfish count -t "+ConvertToString(numThreads)+" -C -m "+ConvertToString( kmerSize)+" -s "+Mem+"M -o "+tempFile +" "+ inputFileName;
	
	cout<<command<<endl;
	system(command.c_str()); 
	command="jellyfish dump "+tempFile +"* > "+ outputFileName;
	cout<<command<<endl;
	system(command.c_str());
	command="rm "+tempFile+"*";
	system(command.c_str());
	
    }

public:
    HashSet genomeKmerSet;
    
    jellyfishUtil(){
        FileName="genome.fasta";
	kmerSize=31;
	numThreads=4;
	Mem="100M";
	delimiter = ',';
    }

     bool lookup(const string word)
    {
	HashSetIt it=genomeKmerSet.find(word);
	if (it != genomeKmerSet.end())
	    return true;
	else
	{
	    it=genomeKmerSet.find(util.reverseComplement( word));
	    if (it != genomeKmerSet.end())
		return true;
	}
	return false;
    }
    void fillSetwithFileContent(){
	string outputFileName="temp.fasta";
	generateKmerFile(FileName,outputFileName);
	//put all previously extracted kmers from genome to the hash set
	fillSet(genomeKmerSet,outputFileName);
	string command="rm "+outputFileName;
	system(command.c_str());
    }
};
class nodeAnalysis{
    size_t kmerSize;
    jellyfishUtil jellyfish;
    string nodeFileName;
    commonUtil util;
    float erronousNodeThreshold;
public:
    nodeAnalysis(){
	kmerSize=31;
	jellyfish.fillSetwithFileContent();
	erronousNodeThreshold=.1;
	nodeFileName="nodes.stage3";
    }
    
    bool traverseNodeFile(){
	ifstream  readFastaFileStream;
	readFastaFileStream.open(nodeFileName.c_str());
	size_t numofcorrectKmerInNode;
	if(!readFastaFileStream.good()){
	    std::cerr << "Error opening file"<< nodeFileName << std::endl;
	    return false;
	}
	double numOfcorreckmerInNode=0;
	double numOferronouskmerInNode=0;
	string line;
	string infoLine;
	size_t j=1;
	size_t maxEsize=0;
	double avgOfEsoze=0;
	double numberOfincorrect=0;
	while (getline(readFastaFileStream,line ))
	{
	    j++;
	    if (j%2==1 ){
		numOfcorreckmerInNode=0;
		numOferronouskmerInNode=0;
		int i=kmerSize;
		while(i<line.length()+1){
		    string kmer =line.substr (i-kmerSize,kmerSize);
		    i++;
		    if (!util.cheakKmerHelth(kmer,kmerSize ))
			continue;
		    if (jellyfish.lookup(kmer)){
			numOfcorreckmerInNode++;
		    }
		    else
		    {
			numOferronouskmerInNode++;
		    }
		}
		if (numOferronouskmerInNode/(numOferronouskmerInNode+numOfcorreckmerInNode)>.1){
		    numberOfincorrect++;
		if (maxEsize<i-kmerSize)
		    maxEsize=i-kmerSize;
		avgOfEsoze=(avgOfEsoze*(numberOfincorrect-1)+(i-kmerSize))/numberOfincorrect;
		//cout<<line<<endl;
		//cout <<line.length()<<endl;
		//cout<<"marginal length of  node		:"<<i-kmerSize<<endl;
		//cout<<"number of correctKmer in node	:"<<numOfcorreckmerInNode<<endl;
		//cout<<"number of inccorrect kmer in node:"<<numOferronouskmerInNode<<endl;
		//cout <<"info line			:"<<infoLine<<endl;
		//cout<<"****************************************************************************************************"<<endl;
		}

	    }
	    else{
		infoLine=line;
	    }

	}
	cout <<"max error size is 		:"<<maxEsize<<endl;
	cout <<"num of Incorrect nodes		:"<<numberOfincorrect<<endl;
	cout <<"avg of Incorrect nodes size	:"<<avgOfEsoze<<endl;
	return true;
	
    }
    
};
class CutOffAnalysis
{
public:
    CutOffAnalysis(){
	
    }
    commonUtil util; 
    float calKmerFreGainValue(float cutOffVal, string fileToProcess, bool print)
    {
	double FP = 0,  FN = 0, TP = 0, TN = 0;
	float index = 0,   numOfcorrect = 0,  numOfIncorrect = 0;
	float gainValue;
	
	ifstream  frequencyStream;
	frequencyStream.open(fileToProcess.c_str());
	if(!frequencyStream.good()){
	    std::cerr << "Error opening file "<< fileToProcess << std::endl;
	    return false;
	}
	
	string line;
	frequencyStream>>line;
	while (!frequencyStream.eof())
	{
	    frequencyStream>>line;
	    vector<string> tokens= util.digestLine(line);
	    index =util.strToDouble(tokens[0]);
	    numOfcorrect= util.strToDouble(tokens[2]);
	    numOfIncorrect= util.strToDouble( tokens[3]);
	    if (index>=cutOffVal) {
		TN =TN+ numOfcorrect;
		FN =FN+numOfIncorrect;
		
	    }else{
		TP=TP+numOfIncorrect;
		FP=FP+numOfcorrect;
	    }
	}
	frequencyStream.close();
	if ((FN+TP)!=0)
	    gainValue = (TP-FP)/(FN+TP);
	else
	    return 0;
	if (print){
	    
	    cout<<"**************************************************************************************"<<endl;
	    cout << " gain value for cut off is "+ util.ConvertToString(cutOffVal)+" is ("<<100*((double)(TP-FP)/(double)(TP+FN))<< "%)"<<endl;
	    cout<< "TP:	"<<TP<<"	TN:	"<<TN<<"	FP:	"<<FP<<"	FN:	"<<FN<<endl;
	    cout << "Sensitivity: ("<<100*((double)TP/(double)(TP+FN))<<"%)"<<endl;
	    cout<<"Specificity: ("<<100*((double)TN/(double)(TN+FP))<<"%)"<<endl;
	    
	}
	return gainValue;
    }
    float getBestGain(string fileToProcess){
	float max=0;
	float cur=0;
	float bestIndex=0;
	float c=1;
	while (c<50){
	    cur=calKmerFreGainValue(c, fileToProcess, false);
	    if (cur>max){
		max=cur;
		bestIndex=c;
	    }
	    c=c+.1;
	}
	calKmerFreGainValue(bestIndex ,fileToProcess, true);
	//calKmerFreGainValue(bestIndex,fileToProcess);
	cout<<"********************"<<endl;
	cout <<"The best gain value is: "<<max*100<< "% This is found for cut of value equal to :"<<bestIndex<< endl;
	return cur;
    }
};

class ExpMaxClustering{
    
private:
    float perErronousClusterMean;
    float perCorrectClusterMean;
    float curErronousClusterMean;
    float curCorrectClusterMean;
    string erronousNodesFileName;
    string correctNodesFileName;
    std::string frequencyFileName;
    double intersectionPoint;
    float divergenceThreshold;
    commonUtil util;
    CutOffAnalysis cutOff;    
    void initialization(string FileName, string firstClusterFileName, string secondClusterFileName,float divergenceValue,float initialMeanOfFirstCluster,float initialMeanOfSecondCluster  ){
	divergenceThreshold=divergenceValue;
	perCorrectClusterMean=initialMeanOfSecondCluster;
	perErronousClusterMean=initialMeanOfFirstCluster;
	curErronousClusterMean=perErronousClusterMean+divergenceThreshold*2;
	curCorrectClusterMean=perCorrectClusterMean+divergenceThreshold*2;
	frequencyFileName=FileName;
	erronousNodesFileName=firstClusterFileName;//"erronousCluster.dat";
	correctNodesFileName=secondClusterFileName;//"correctCluster.dat";
    }
    
    float calculateMean(string clusterFileName){
	ifstream  frequencyStream;
	frequencyStream.open(clusterFileName.c_str());
	string line;
	double sumOfAll=0;
	double totalNum=0;
	float avg=0;
	double num=0;
	double index=0;
	//cout<<clusterFileName<<endl;
	while (!frequencyStream.eof())
	{
	    frequencyStream>>line;
	    vector<std::string> tokens= util.digestLine(line);
	    index =util.strToDouble(tokens[0]);
	    num= util.strToDouble(tokens[1]);
	    //cout<<num<<endl;
	    sumOfAll=sumOfAll+index*num;
	    avg=(avg*totalNum+index*num)/(totalNum+num);
	    totalNum=totalNum+num;
	}
	//avg=sumOfAll/totalNum;
	frequencyStream.close();
	return avg;
    }
    bool isErroneous(double index){
	float erronousProb=gsl_ran_poisson_pdf(index,curErronousClusterMean);
	float correctProb=gsl_ran_poisson_pdf(index,curCorrectClusterMean);
	if (erronousProb/correctProb>1)
	    return true;
	return false;
	
    }
    bool classifier(){
	double  numOfcorrect = 0,  numOfIncorrect = 0, sumOfBoth=0;
	double index;
	ifstream  frequencyStream;
	ofstream correctNodeStream;
	ofstream erronousNodeStream;
	correctNodeStream.open(correctNodesFileName.c_str(), ios::out);
	erronousNodeStream.open(erronousNodesFileName.c_str(), ios::out);
	frequencyStream.open(frequencyFileName.c_str());
	correctNodeStream.setf(ios_base::fixed);
	erronousNodeStream.setf(ios_base::fixed);
	if(!frequencyStream.good()){
	    std::cerr << "Error opening file "<< frequencyFileName << std::endl;
	    return false;
	}
	string line;
	frequencyStream>>line;
	while (!frequencyStream.eof())
	{
	    frequencyStream>>line;
	    vector<string> tokens= util.digestLine(line);
	    index =util.strToDouble(tokens[0]);
	    numOfcorrect= util.strToDouble(tokens[2]);
	    numOfIncorrect= util.strToDouble( tokens[3]);
	    //cout<<fixed<<setprecision(2)<<numOfIncorrect<<endl;
	    sumOfBoth=numOfcorrect+numOfIncorrect;
	    if (isErroneous(index))
		erronousNodeStream<<fixed<<setprecision(2)<<index<<"," <<util.ConvertToString(sumOfBoth)<<","<<setprecision(4)<< gsl_ran_poisson_pdf(index,curErronousClusterMean)<<endl;
	    else
		correctNodeStream<<fixed<<setprecision(2)<<index<<","<<fixed <<util.ConvertToString(sumOfBoth)<<","<<setprecision(4)<<gsl_ran_poisson_pdf(index,curCorrectClusterMean) <<endl;
	}
	correctNodeStream.close();
	erronousNodeStream.close();
	frequencyStream.close();
	return true;
    }
    void updateParameter(){
	perCorrectClusterMean=curCorrectClusterMean;
	perErronousClusterMean=curErronousClusterMean;
	curErronousClusterMean=calculateMean(erronousNodesFileName);
	curCorrectClusterMean=calculateMean(correctNodesFileName);
    }
    
    
    void findIntersectionPoint(){
	double e=2.718281;
	double c=curErronousClusterMean/curCorrectClusterMean;
	intersectionPoint=(curErronousClusterMean-curCorrectClusterMean)* (log(e)/log(c));
	cout<<"The intersection of these curves is :"<<intersectionPoint<<endl;
    }
    
    
public:
    ExpMaxClustering(int argc, char ** args){
	parseInputArguments(argc,args);
    }
    ExpMaxClustering(string fileName){
	initialization(fileName, "erronousNodes.dat","correctNodes.dat",.01,1, 50 );
    }
    ExpMaxClustering(){
	initialization("scov_001.dat", "erronousNodes.dat","correctNodes.dat",.01,1, 50 );
    }
    void doClassification(){
	size_t i=0;
	while(fabs( perErronousClusterMean-curErronousClusterMean)>divergenceThreshold ||fabs( perCorrectClusterMean-curCorrectClusterMean)>divergenceThreshold){
	    //cout<<"correctMean :	"<<curCorrectClusterMean<<endl;
	    //cout<<"erronousMean:	"<<curErronousClusterMean<<endl;
	    //cout<<"next round------------------------"<<endl;
	    classifier();
	    updateParameter();
	    i++;
	    if (i>20)
		break;
	}
	cout<<"Mean of nodes in correct nodes cluster	:	"<<curCorrectClusterMean<< endl;
	cout<<"Mean of nodes in erronous nodes cluster	:	"<<curErronousClusterMean<<endl;
	findIntersectionPoint();
	cutOff.calKmerFreGainValue(intersectionPoint, frequencyFileName,true);
	cutOff.getBestGain(frequencyFileName);
    }
    void parseInputArguments(int argc, char ** args){
	string firstClusterFileName="erronousNodes.dat";
	string secondClusterFileName="correctNodes.dat"; 
	double divergenceValue=.01;
	double initialMeanOfFirstCluster=1; 
	double initialMeanOfSecondCluster=50;
	string frequency_fileName="";
	if (argc<3){
	    cout<<"wrong way of using this software"<<endl;
	    printUsage();
	    exit(EXIT_SUCCESS);
	}
	
	for (int i = 2; i < argc; i++) {
	    string arg(args[i]);
	    if (arg.empty())
		continue;              // this shouldn't happen
		
		if ((arg == "-h") || (arg == "--help")) {
		    printUsage();
		    exit(EXIT_SUCCESS);
		}else if ((arg == "-of") || (arg == "--firstOutput")) {
		    i++;
		    if (i < argc){
			firstClusterFileName = args[i]; 
		    }
		    
		}else if ((arg == "-os") || (arg == "--secondOutput")) {
		    i++;
		    if (i < argc){
			secondClusterFileName =args[i];
		    }
		    
		}else if ((arg == "-d") || (arg == "--divergenceValue")) {
		    i++;
		    if (i < argc){
			divergenceValue =util.strToDouble(args[i]);
		    }
		    
		}
		else if ((arg == "-ifm") || (arg == "--firstInitialMean")) {
		    i++;
		    if (i < argc){
			initialMeanOfFirstCluster =util.strToDouble( args[i]);
		    }
		    
		}
		else if ((arg == "-ism") || (arg == "--secondInitialMean")) {
		    i++;
		    if (i < argc){
			initialMeanOfSecondCluster =util.strToDouble( args[i]);
		    }   
		}
		else {        // it must be an input file
		    frequency_fileName = args[i];
		    initialization( frequency_fileName,  firstClusterFileName, secondClusterFileName, divergenceValue, initialMeanOfFirstCluster, initialMeanOfSecondCluster );
		    doClassification();
		}
	}
	
	
    }
    void printUsage(){
	cout<<"This command do the clustering operation with Expectation Maximization method."<<endl;
	cout<<"There are two known clusters"<<endl;
	cout <<"Distribution of both clusters is Poisson."<<endl;
	cout<<"Usage: 	BrowniePlus ExpectationMaximization [file_options] frequency_fileName"<<endl;
	cout<<"options"<<endl;
	cout<<"-of	name of first ouptut cluster"<<endl;
	cout<<"-os	name of second ouptut cluster"<<endl;
	cout<<"-ifm	the initial Mean vlaue for the first cluster"<<endl;
	cout<<"-ism	the initial Mean value for the second cluster"<<endl;
	cout<<"-d	divergence Value"<<endl;
    }
    
};


typedef google::sparse_hash_set<std::string> HashSet;
typedef HashSet::const_iterator HashSetIt;
typedef map< int , int > freqVec;
typedef freqVec::const_iterator freqVecIt; 

freqVec correctKmers;
freqVec erronouskmers;

string ReadsFileName;
string genomeFileName;
string plotOutFileName;
string frequencyDatFileName;
string nodeCoveDatFile;
string outputfileNameforSampling;
string lowFrequentFileName;
int lowFrequentthreshold=10;
float samplinPer;
size_t kmerSize=31;
size_t numThreads=4;
string Mem="200M";
char delimiter = ',';
string ConvertToString (float number){
    std::ostringstream buff;
    buff<<number;
    return buff.str();   
}
float strToFloat(string str){
    return (atof(str.c_str()));
}
bool cheakKmerHelth(string &kmer){
    if (kmer.length()!=kmerSize)
	return false;
    for (size_t i=0;i<kmerSize;i++){
	if( (kmer[i]!='A'&& kmer[i]!='T'&&kmer[i]!='C'&&kmer[i]!='G'&&kmer[i]!='a'&&kmer[i]!='t'&&kmer[i]!='c'&&kmer[i]!='g'))
	    return false;
    }
    for (size_t i=0;i<kmerSize;i++){
	if(kmer[i]=='a')
	    kmer[i]='A';
	if (kmer[i]=='t')
	    kmer[i]='T';
	if (kmer[i]=='c')
	    kmer[i]='C';
	if (kmer[i]=='g')
	    kmer[i]='G';
    }
    
    return true;
}
string reverseComplement(string kmer){
    string result="";
    for (int i=kmer.length()-1;i>=0;i--){
	switch(kmer[i]){
	    case 'A':
		result=result+'T';
		break;
	    case 'T':
		result=result+'A';
		break;
	    case 'C':
		result=result+'G';
		break;
	    case 'G':
		result=result+'C';
		break;
	}
    }
    return result;
}
bool lookup(const HashSet &Set,
	    const string word)
{
    HashSetIt it=Set.find(word);
    if (it != Set.end())
	return true;
    else
    {
	it=Set.find(reverseComplement( word));
	if (it != Set.end())
	    return true;
    }
    return false;
}
bool fillSet(HashSet &Set, string fileName)
{  
    ifstream  readFastaFileStream;
    readFastaFileStream.open(fileName.c_str());
    if(!readFastaFileStream.good()){
	std::cerr << "Error opening file"<< fileName << std::endl;
	return false;
    }
    string line;
    while (!readFastaFileStream.eof())
    {
	readFastaFileStream>>line;
	if (line[0] != '>' ){
	    string kmer =line;                                     //line.substr (0,kmerSize-1);
	    if (!cheakKmerHelth(kmer))
		continue;
	    Set.insert(kmer);
	}
	
    }
    //Set.insert( reverseComplement("gtctggatgactggacccacaacaaaagcat"));
    cout<<Set.size()<<endl;
    readFastaFileStream.close();
    return true;
    
}
void generateKmerFile(string inputFileName, string outputFileName){
    string tempFile="tempKmerFreq";
    string command="jellyfish count -t "+ConvertToString(numThreads)+" -C -m "+ConvertToString( kmerSize)+" -s "+Mem+"M -o "+tempFile +" "+ inputFileName;
    cout<<command<<endl;
    system(command.c_str()); 
    command="jellyfish dump "+tempFile +"* > "+ outputFileName;
    cout<<command<<endl;
    system(command.c_str());
    command="rm "+tempFile+"*";
    system(command.c_str());
    
}

bool separateTrueFalseKmers(HashSet& genomeKmerSet, string inputFileName){
    cout<<"seperating true and false kmers  starting soon........."<<endl;
    ifstream  readFastaFileStream;
    readFastaFileStream.open(inputFileName.c_str());
    
    if(!readFastaFileStream.good()){
	std::cerr << "Error opening file"<< inputFileName << std::endl;
	return false;
    }
    string line;
    string numStr;
    int num;
    string kmer;
    int max=0;
    int i=1;
    while (getline(readFastaFileStream,line ))
    {
	if (i%100000==0)
	{
	    cout << "Processing kmer number " << i/2<<"\r";
	    cout.flush();
	}
	i++;
	if (line[0] == '>' ){
	    numStr =line.substr (1,line.length());
	    num= atoi(numStr.c_str());
	    if (num>max)
		max=num;
	}
	else
	{
	    kmer =line;                                            //(line.substr (0,kmerSize-1));
	    if (!cheakKmerHelth(kmer))
		continue;
	    if (lookup(genomeKmerSet,kmer))
		correctKmers[num]=correctKmers[num]+1;
	    
	    else
		erronouskmers[num]=erronouskmers[num]+1;
	}
	
    }
    readFastaFileStream.close();   
    cout<<"Writing true and false kmer frequency in file ..."<<endl;
    ofstream notCorrectedReadF;
    notCorrectedReadF.open(frequencyDatFileName.c_str(), ios::out);
    notCorrectedReadF<<"index1"<<delimiter<<"index2"<<delimiter<<"numOfcorrect" <<delimiter<<"numOfIncorrect"<<endl;
    i=1;
    while(i<max)
    {
	int numOfcorrect=correctKmers[i];
	int numOfIncorrect=erronouskmers[i];
	notCorrectedReadF<<i<<delimiter<<i+.5<<delimiter<<numOfcorrect<<delimiter<<numOfIncorrect<<endl;
	i++;
    }
    notCorrectedReadF.close();
    return true;
}
bool separateTrueFalseKmers(HashSet& genomeKmerSet, string inputFileName, string outputCorrectKmerFile){
    ofstream outputCorrectStream;
    outputCorrectStream.open(outputCorrectKmerFile.c_str(), ios::out);
    
    cout<<"seperating true and false kmers and Writing low freqent kmer in seperate file starting soon........."<<endl;
    ifstream  readFastaFileStream;
    readFastaFileStream.open(inputFileName.c_str());
    
    if(!readFastaFileStream.good()){
	std::cerr << "Error opening file"<< inputFileName << std::endl;
	return false;
    }
    string line;
    string numStr;
    int num;
    string kmer;
    int max=0;
    int i=1;
    while (getline(readFastaFileStream,line ))
    {
	if (i%100000==0)
	{
	    cout << "Processing kmer number " << i/2<<"\r";
	    cout.flush();
	}
	i++;
	if (line[0] == '>' ){
	    numStr =line.substr (1,line.length());
	    num= atoi(numStr.c_str());
	    if (num>max)
		max=num;
	}
	else
	{
	    kmer =line;                                            //(line.substr (0,kmerSize-1));
	    if (!cheakKmerHelth(kmer))
		continue;
	    if (lookup(genomeKmerSet,kmer)){
		correctKmers[num]=correctKmers[num]+1;
		if (num<10){
		    outputCorrectStream<<">"+numStr<<endl;
		    outputCorrectStream<<kmer<<endl;
		}
	    }
	}
	
    }
    readFastaFileStream.close();   
    cout<<"Writing true and false kmer frequency in file ..."<<endl;
    ofstream notCorrectedReadF;
    notCorrectedReadF.open(frequencyDatFileName.c_str(), ios::out);
    notCorrectedReadF<<"index1"<<delimiter<<"index2"<<delimiter<<"numOfcorrect" <<delimiter<<"numOfIncorrect"<<endl;
    i=1;
    while(i<max)
    {
	int numOfcorrect=correctKmers[i];
	int numOfIncorrect=erronouskmers[i];
	notCorrectedReadF<<i<<delimiter<<i+.5<<delimiter<<numOfcorrect<<delimiter<<numOfIncorrect<<endl;
	i++;
    }
    outputCorrectStream.close();
    notCorrectedReadF.close();
    return true;
}
void makeGnuPlot(){
    string command="gnuplot -e  \"filename='"+frequencyDatFileName+"'\" -e \"outputfile='"+plotOutFileName+"'\" plot.dem";
    cout<<command<<endl;
    system(command.c_str());
}


void makePlot(bool gnuplot, bool printlowFrequentkmer ){
    HashSet genomeKmerSet;
    // generate all kmers from genom file
    string outputFileName="genomekmerFreq.fasta";
    generateKmerFile(genomeFileName,outputFileName);
    //put all previously extracted kmers from genome to the hash set
    fillSet(genomeKmerSet,outputFileName);
    string command="rm "+outputFileName;
    system(command.c_str());
    //extract all kmers from read file
    outputFileName="readkmerFreq.fasta";
    generateKmerFile(ReadsFileName,outputFileName);
    if (!printlowFrequentkmer)
	separateTrueFalseKmers(genomeKmerSet,outputFileName );
    else
	separateTrueFalseKmers(genomeKmerSet,outputFileName, lowFrequentFileName );
    if (gnuplot)
	makeGnuPlot();
    command="rm "+outputFileName;
    system(command.c_str());
    //gnuplot -e  "kmerfile='1/frequency.dat'" -e  "filename='1/scov_001.dat'" -e "outputfile='1/mixCovKmer1.pdf'" -e "xlimit=50" -e "ylimit=10000" -e "xtic=2"  mixplot.dem
}
vector<string> digestLine(string line) {
    vector<string> elems;
    string token="";
    size_t i=0;
    while(i<line.length()){
	if (line[i]==delimiter){
	    elems.push_back(token);
	    token="";
	}
	else
	    token=token+line[i];
	i++;
    }
    elems.push_back(token);
    return elems;
}

bool is_digits(const char &c)
{
    if (c=='0'||c=='1'||c=='2'||c=='3'||c=='4'||c=='5'||c=='6'||c=='7'||c=='8'||c=='9'||c=='.' )
	return true;
    return false;
}
float calKmerFreGainValue(float cutOffVal, string fileToProcess){
    double FP = 0,  FN = 0, TP = 0, TN = 0;
    float index = 0,   numOfcorrect = 0,  numOfIncorrect = 0;
    float gainValue;
    
    ifstream  frequencyStream;
    frequencyStream.open(fileToProcess.c_str());
    if(!frequencyStream.good()){
	std::cerr << "Error opening file "<< fileToProcess << std::endl;
	return false;
    }
    
    string line;
    frequencyStream>>line;
    while (!frequencyStream.eof())
    {
	frequencyStream>>line;
	vector<string> tokens= digestLine(line);
	index =strToFloat(tokens[0]);
	numOfcorrect= strToFloat(tokens[2]);
	numOfIncorrect= strToFloat( tokens[3]);
	if (index>=cutOffVal) {
	    TN =TN+ numOfcorrect;
	    FN =FN+numOfIncorrect;
	    
	}else{
	    TP=TP+numOfIncorrect;
	    FP=FP+numOfcorrect;
	}
    }
    frequencyStream.close();
    if ((FN+TP)!=0)
	gainValue = (TP-FP)/(FN+TP);
    else
	return 0;
    /*cout<<"**************************************************************************************"<<endl;
     *    cout << " gain value for cutoff-"+ ConvertToString(cutOffVal)+" is ("<<100*((double)(TP-FP)/(double)(TP+FN))<< "%)"<<endl;
     *    cout<< "TP:	"<<TP<<"	TN:	"<<TN<<"	FP:	"<<FP<<"	FN:	"<<FN<<endl;
     *    cout << "Sensitivity: ("<<100*((double)TP/(double)(TP+FN))<<"%)"<<endl;
     *    cout<<"Specificity: ("<<100*((double)TN/(double)(TN+FP))<<"%)"<<endl;*/
    return gainValue;
}

float getBestGain(string fileToProcess){
    float max=0;
    float cur=0;
    float bestIndex=0;
    float c=1;
    while (c<200){
	cur=calKmerFreGainValue(c, fileToProcess);
	if (cur>max){
	    max=cur;
	    bestIndex=c;
	}
	c=c+.5;
    }
    //calKmerFreGainValue(bestIndex,fileToProcess);
    cout<<"********************"<<endl;
    cout <<"The best gain value is: "<<max*100<< "% This is found for cut of value equal to :"<<bestIndex<< endl;
    return cur;
}
void makeSampleReadFile(string readFileNameToSamle,float percentage)
{
    
    ifstream readsFile;
    ofstream sampleReadFile;
    string perStr=ConvertToString(percentage);
    cout<<perStr<<endl;
    sampleReadFile.open(outputfileNameforSampling.c_str(), ios::out);
    readsFile.open(readFileNameToSamle.c_str(),ios::in);
    string firstLine="", secodnLine="",thirdLine="",forthLine="";
    while(getline(readsFile,firstLine )&&getline(readsFile,secodnLine )&& getline(readsFile,thirdLine )&&getline(readsFile,forthLine )) {
	int r = ((rand()%100));
	if (r>percentage)
	    continue;
	sampleReadFile<<firstLine<<endl;
	sampleReadFile<<secodnLine<<endl;
	sampleReadFile<<thirdLine<<endl;
	sampleReadFile<<forthLine<<endl;
    }
    readsFile.close();
    sampleReadFile.close();
}
void testroutin(string fileName){
    /*cout<<"seperating true and false kmers starting soon........."<<endl;
     *  ifstream  nodeFileStream;
     *  nodeFileStream.open(fileName.c_str());
     *  if(!nodeFileStream.good()){
     *    std::cerr << "Error opening file"<< fileName << std::endl;
     *    return false;
}
while(getline(nodeFileStream,firstLine )&&getline(nodeFileStream,secodnLine )) {
    for (int i=0;i<secodnLine.length()-kmerSize;i++){
	string kmer=string.substr(i,kmerSize);
	
}
}
nodeFileStream.close()*/
    
}

void printUsage()
{ 
    cout<<"BrowniePlus"<<endl;
    cout<<"contact mahdi.heydari@ibcn.ugent.be"<<endl;
    cout<<"Usage: 	BrowniePlus<command> [options]"<<endl;
    cout<<"Command:	comparison	Evaluation of error correction method"<<endl;
    cout<<"Command:	sampling	Make smaple of original read file"<<endl;
    cout<<"Command:	kmerAnalysis	Analysing kmer frequency and node coverage frequency"<<endl;  
    cout<<"Command:	clustering 	Clustering based on Expectation Maximization method"<<endl;
}
void printUsageForSampling(){
    cout<<"Usage: 	BrowniePlus sampling [file_options] file"<<endl;
    cout<<"options"<<endl;
    cout<<"-o	outpuFileName"<<endl;
    cout<<"-p	percentage of read file which you want to make a smaple of it (defult is 5)"<<endl;
}

void parseCommandForSampling(int argc, char ** args){
    if (argc<3){
	cout<<"wrong way of using this software"<<endl;
	printUsageForSampling();
	exit(EXIT_SUCCESS);
    }
    samplinPer=5;  
    bool outputFileNotGiven=true;;
    for (int i = 2; i < argc; i++) {
	string arg(args[i]);
	if (arg.empty())
	    continue;              // this shouldn't happen
	    
	    if ((arg == "-h") || (arg == "--help")) {
		printUsageForSampling();
		exit(EXIT_SUCCESS);
	    }else if ((arg == "-p") || (arg == "--percentage")) {
		i++;
		if (i < argc){
		    samplinPer = strToFloat(args[i]);  
		}
		
	    }else if ((arg == "-o") || (arg == "--output")) {
		i++;
		if (i < argc){
		    outputfileNameforSampling =args[i];
		    outputFileNotGiven=false;
		}
		
	    } else {        // it must be an input file
		ReadsFileName = args[i];
	    }
    }
    if (outputFileNotGiven)
	outputfileNameforSampling="S_"+ConvertToString( samplinPer)+"_"+ReadsFileName;
    makeSampleReadFile(ReadsFileName,samplinPer);
    
}
void printUsageCommandForKmerAnalysis(){
    cout<<"Usage: 	BrowniePlus kmerAnalysis [file_options] readFile -g genomeFile -ff frequencyFileNmae -cf nodeKmerCovFileName "<<endl;
    cout<<"options"<<endl;
    cout<<"-k	Kmer size, defult value is 31"<<endl;
    cout<<"-m	Amount of memory dedicated to jellyfish, defult value is 400M"<<endl;  
    cout<<"-d	Delimiter for Writing frequency result , defult value is ,"<<endl;
    cout<<"-pf	plot file name for making gnuplot frequency histogram"<<endl;  
    cout<<"-nj	Frequency file is ready, jellyfish is not needed "<<endl;
    cout<<"-np	Gnuplot diagram not needed"<<endl;  
}
void parseCommandForKmerAnalysis(int argc, char ** args){
    
    bool jellyfish=true;
    bool gnuplot=true;
    bool doclustering=true;
    bool lowFrequentgiven=false;
    bool genomeFileIsGiven=false;
    bool nodeKmerCovFileIsGiven=false;
    bool frequenyFileIsGive=false;
    for (int i = 2; i < argc; i++) {
	string arg(args[i]);
	if (arg.empty())
	    continue;              // this shouldn't happen
	    if ((arg == "-h") || (arg == "--help")) {
		printUsageCommandForKmerAnalysis();
		exit(EXIT_SUCCESS);
	    } else if ((arg == "-k") || (arg == "--kmersize")) {
		i++;
		if (i < argc)
		    kmerSize = atoi(args[i]);
	    } else if ((arg == "-t") || (arg == "--threads")) {
		i++;
		if (i < argc)
		    numThreads = atoi(args[i]);
	    } else if ((arg == "-g") || (arg == "--genomeFile")) {
		i++;
		if (i < argc){
		    genomeFileName = args[i];
		    genomeFileIsGiven=true;
		}
	    }else if ((arg == "-m") || (arg == "--memoryForJellyfish")) {
		i++;
		if (i < argc)
		    Mem = args[i];
	    }
	    else if ((arg == "-d") || (arg == "--delimiter")) {
		i++;
		if (i < argc)
		    delimiter = args[i][0];
	    } 
	    else if ((arg == "-pf") || (arg == "--plotFileName")) {
		i++;
		if (i < argc)
		    plotOutFileName = args[i];
	    } 
	    else if ((arg == "-ff") || (arg == "--frequenyFileName")) {
		i++;
		if (i < argc){
		    frequencyDatFileName = args[i];
		    frequenyFileIsGive=true;
		}
	    }
	    
	    else if ((arg == "-lf") || (arg == "--lowFrequencykmer")) {
		i++;
		if (i < argc){
		    lowFrequentFileName = args[i];
		    lowFrequentgiven=true;
		}
	    }
	    
	    else if ((arg == "-cf") || (arg == "--nodeKmerCovFileName")) {
		i++;
		if (i < argc){
		    nodeCoveDatFile = args[i];
		    nodeKmerCovFileIsGiven=true;
		}
	    }
	    else if ((arg == "-nj") || (arg == "--jellyfishNotNedded")) {
		jellyfish=false;
	    }
	    else if ((arg == "-np") || (arg == "--gnuPlotNotNedded")) {
		gnuplot=false;
		
	    }
	    else if ((arg == "-nc") || (arg == "--clusteringNotNedded")) {
		doclustering=false;
		
	    }
	    else {
		
		ReadsFileName = args[i];
	    }
	    
    }
    if (genomeFileIsGiven&&frequenyFileIsGive&&nodeKmerCovFileIsGiven)
    {
	if (jellyfish)
	    makePlot(gnuplot, lowFrequentgiven);	
	cout <<"Finding best gain for kmer coverage ................"<<endl;
	getBestGain(frequencyDatFileName);
	cout <<"Finding best gain for node kmer coverage ................"<<endl;
	getBestGain(nodeCoveDatFile);
	if (doclustering){
	    cout<<".................................."<<endl<<endl;    
	    cout <<"Finding cut off value based on clustering ................"<<endl;
	    ExpMaxClustering exp(nodeCoveDatFile);
	    exp.doClassification();
	}
    }
    else{
	printUsageCommandForKmerAnalysis();
    }
    
    
}
void parseCommandForComparison(int argc, char ** args){
    /* for (int i = 2; i < argc; i++) {
     *    string arg(args[i]);
     *    
     *    if (arg.empty())
     *      continue;              // this shouldn't happen
     *      
     *      if ((arg == "-h") || (arg == "--help")) {
     *	printUsage();
     *	exit(EXIT_SUCCESS);
} else if ((arg == "-g") || (arg == "--genomeFile")) {
    i++;
    if (i < argc)
	genomeFileName = atoi(args[i]);
} else if ((arg == "-o") || (arg == "--output")) {
    i++;
    if (i < argc)
	outputFilename = args[i];
} else {        // it must be an input file
    inputFilename = args[i];
}
}*/
}
void parseCommandLineArguments(int argc, char ** args)
{
    if (argc<2){
	cout<<"wrong way of using this software"<<endl;
	printUsage();
	exit(EXIT_SUCCESS);
    }
    string command=args[1];
    if (command=="comparison"){
	parseCommandForComparison(argc,args);
    }else if (command=="clustering"){
	ExpMaxClustering exp;
	exp.parseInputArguments(argc,args);
    }
    else if (command=="sampling"){
	parseCommandForSampling(argc,args);
    }else if (command=="kmerAnalysis"){
	parseCommandForKmerAnalysis(argc,args);
    }else {
	cout<<"wrong way of using this software"<<endl;
	printUsage();
	exit(EXIT_SUCCESS);
    }
    
}
int main(int argc, char ** args)
{
    cout<<"Welcom to BrowniePlus"<<endl;
    cout << "----------------------------------------------------\n" << endl;
    /*if (argc != 6) {
     *    cout << "ERROR: Wrong number of arguments!\n";
     *    cout << "Usage: ./BrowniePlus <reads.fasq> <genome.fasta> <plotOutFileName> <frequencyDatFileName> <nodeCoveDatFile>" << endl << flush;
     *    cout << "this software produces kmer distribution plot.\n Jellyfish and gnuplot is used in this software.\nPlease look at the readme file for more information." << endl << flush;
     *    return 0;
    }
    */
    //"how to compile:"
    //g++ -lgsl  hashSet.cpp -o browniePlus.exe -lm
    parseCommandLineArguments(argc,args);
    //nodeAnalysis node;
    //node.traverseNodeFile("nodes.stage3");
    cout<<"bye bye, ...... see you soon"<<endl;
}

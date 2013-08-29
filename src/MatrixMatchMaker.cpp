// MatrixMatchMaker
// Version 2010.10.31
// (c) 2008-2010, Robert L. Charlebois and Elisabeth R. M. Tillier

// Do not distribute

#include "classCmdLineArgParser.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>

#include "mmm_algorithm.h"
#include "timer.h"
#include "sighand.h"
#include "MMML.h"

#include "Alignment.hpp"
#include "DistanceMatrix.hpp"

#define ALIGNMENT_CACHE_SIZE 10
#define DISTANE_MATRIX_CACHE_SIZE 1000

using namespace std;


struct CmpTreesByWeight
{
	bool operator()(const WeightedProteinPairVec& a, const WeightedProteinPairVec& b)
	{
		return a.second > b.second;
	}
};


static void printUsage1(const string& pname);
static void printUsage2(const string& pname);
static void printUsage(const string& pname);
static void printVerboseUsage(const string& pname);
static int readDistFile(const string& which, const string& distfileName, vector<double>& distances, vector<string>& names,
						vector<string>& taxLabels, bool useTaxInfo);
static size_t findBiggestPossibleMatch(const vector<string>& taxLabels1, const vector<string>& taxLabels2,
									   const bool useTaxInfo);
//static void processTaxonLabels(vector<string>& taxLabels1, vector<string>& taxLabels2, string& reqTaxName,
//							   vector<int>& taxIndices1, vector<int>& taxIndices2, vector<string>& taxLabels, int& reqTaxIndex);
static void processTaxonLabels(vector<string> taxLabels, const string reqTaxName, vector<int>& taxIndices,
										 vector<string>& taxLabelsGlobal, int& reqTaxIndexGlobal);
static string itoa(int number);
static void showStats(vector<int>& taxIndices1, vector<int>& taxIndices2, vector<string>& taxLabels);
static int calcTaxonCombinations(vector<int>& taxIndices1, vector<int>& taxIndices2);
static string makeStateSaveFileName(const string& first, const string& second);
static void processDistances(vector<double>& distances, double cutoffLower, double cutoffUpper, bool logDistances);
static void handleErrors(bool ignoreErrors, ofstream& outfile);
static void rearrangeMatrix(vector<string>& names, vector<string>& taxLabels, vector<double>& distances, const string reqTaxName);



static unsigned int getDistMatrix(const string& which, const string& distfileNamePrefix, const string& distfileName,
									  list <DistanceMatrix>& distMatricesCache, list <Alignment>& alignmnentsCache,
									  vector <double>& distances, vector <string>& names, vector <string>& taxLabels,
									  bool useTaxInfo, bool aln2pmb);
static Alignment const& getAlignment(const string& inFilepathPrefix, const string& inFileName, list <Alignment>& alignmnentsCache);
static bool find_in_alignments_cache(list <Alignment>& cache, string const& name);
static bool find_in_dist_matrices_cache(list <DistanceMatrix>& cache, string const& name);
static Alignment getSlice(stringstream& distFilenameSS, Alignment const& alignment);


int main(int argc, const char* argv[])
{
	// Read arguments from the command line:
	string distfileName1, distfileName2, batchFileName, outfileName, reqTaxName, distfileNamesPrefix, stateFileName;
	int strictSize = -1;
	int minSize = 0;
	int maxTrees = numeric_limits<int>::max();
	int rangeStart = -1;
	int rangeEnd = numeric_limits<int>::max();
	bool tabulateTop = false;
	bool useTaxInfo = false;
	bool reqTax = false;
	bool tabDelimitedTable = false;
	bool brief = false;
	bool verbose = false;
	bool doShowStats = false;
	bool outputZeros = false;
	bool stateRestore = false;
	bool exactTime = false;
	bool useHardware = false;
	double allow;
	// MMML extensions: New 2010.02.28: RLC:
	string treeFileName, ignoreTaxaFileName;
	int minBootstrapSupport = 1;
	bool MMMLmode = false, listTaxa = false, topCoverage = false, LverboseMode = false; 
	// New 2010.02.28: RLC.

	
	double cutoffLower = 0.000011; //protdist's threshold of 0.000010 fits here
	double cutoffUpper = 1.0e+10;
	bool logDistances = false;
	
	bool noHeader = false;
	bool onlyHeader = false;
	bool onlyMP = false;
	bool quickMP = false;
	bool ignoreErrors = false;
	
	bool aln2pmb = false;
	
	{
		CmdLineArgParser options(argc, argv);
		doShowStats = options.parse("-stats");
		if (options.parse("-help")) printVerboseUsage(options.programName());
		
		bool haveDist1 = options.parse("-d1", &distfileName1);
		bool haveDist2 = options.parse("-d2", &distfileName2);
		bool haveBatch = options.parse("-b", &batchFileName);

		stateRestore = options.parse("-restore", &stateFileName);

		if (!doShowStats || !haveDist1 || !haveDist2)
		{
			if (!haveBatch && (!haveDist1 || !haveDist2) ||
				!options.parse("-o", &outfileName) || !options.parse("-a", &allow))
			{
				printUsage(options.programName());
			}
		}

		if (haveBatch && stateRestore)
		{
			printUsage(options.programName());
		}

		brief = options.parse("-brief"); // optional; edit New 2010.03.01: RLC, fixed by abezginov
		tabDelimitedTable = (batchFileName.empty() && options.parse("-table")); // optional, precluded if -b is selected
		if (!batchFileName.empty() && !distfileName1.empty()) printUsage(options.programName());
		tabulateTop = options.parse("-top"); // optional;
		if (tabulateTop && !tabDelimitedTable) printUsage(options.programName());
		useTaxInfo = options.parse("-u"); // optional
		options.parse("-sz", &strictSize); // optional
		if (strictSize != -1 && strictSize < 3) {
			cout << "strictSize must be 3 or more." << endl;
			exit(1);
		}
		if (options.parse("-minsize", &minSize) && minSize < 3) // optional
		{
			cout << "minSize must be 3 or more." << endl;
			exit(1);
		}
		options.parse("-maxtrees", &maxTrees); // optional
		verbose = options.parse("-v"); // optional
		options.parse("-req", &reqTaxName); // optional
		doShowStats = options.parse("-stats");
		options.parse("-p", &distfileNamesPrefix);
		outputZeros = options.parse("-z");
		exactTime = options.parse("-exacttime");
		useHardware = options.parse("-hard");
		options.parse("-lower", &cutoffLower);
		options.parse("-upper", &cutoffUpper);
		logDistances = options.parse("-log"); 
		ignoreErrors = options.parse("-ignoreerrors");
		
		
		reqTax = !reqTaxName.empty();
		if (reqTax && !useTaxInfo)
		{	cout << "\nError: Must use taxon info when requiring a taxon." << endl;
			exit(1);
		}
		
		if (reqTax)
			cout << "Warning: -req option is currently broken (may return lower scores than it should). Use with caution." << endl;
			
			
		// MMML extensions: New 2010.02.28: RLC:
		options.parse("-Lt", &treeFileName); // optional
		if (!treeFileName.empty()) {
			if (brief || outputZeros) printUsage(options.programName()); // incompatible options
			MMMLmode = true;
			// Note: outfileName now refers to MMML output
			options.parse("-Li", &ignoreTaxaFileName); // optional
			options.parse("-Lb", &minBootstrapSupport); // optional
			listTaxa = options.parse("-List"); // optional
			topCoverage = options.parse("-LtopCov"); // optional
			LverboseMode = options.parse("-Lv"); // optional
		}
		// New 2010.02.28: RLC.
			
		// abezgino 2011.01.21
		if (brief)
		{	quickMP = options.parse("-quickmp");
			onlyMP = options.parse("-onlymp");
			if (onlyMP) outputZeros = true;
		}
		
		
		if (brief || MMMLmode)
		{	noHeader = options.parse("-noheader"); 
			onlyHeader = options.parse("-onlyheader");
			if (onlyHeader && noHeader)
			{	cout << "Incompatible options: -noheader -onlyheader" << endl;
				printUsage(options.programName());
			}
		}
		
		
		const int foundRangeStart = options.parse("-rangeStart", &rangeStart); // optional
		const int foundRangeEnd = options.parse("-rangeEnd", &rangeEnd); // optional
		if (foundRangeStart ^ foundRangeEnd) printUsage(options.programName());
		
		aln2pmb = options.parse("-aln2pmb"); 
		
		if (!options.empty())
		{	cout << "Unrecognized option:" << endl;
			options.print();
			cout << endl;
			printUsage(options.programName());
		}
	}
	//

	//populate the list of pairs of distance matrices to run
	vector<pair<string, string> > theDistFileNames;
	if (batchFileName.empty())
	{
		theDistFileNames.push_back(pair<string, string>(distfileName1, distfileName2));
	} 
	else
	{
		ifstream batchFile(batchFileName.c_str());
		if (!batchFile.is_open())
		{
			cout << "\nError: cannot open file " << batchFileName << endl;
			exit(1);
		}
		string thePair;
		string::size_type p;
		getline(batchFile, thePair, '\n');
		while (batchFile)
		{	p = thePair.find('\t');
			if (p == string::npos) break;
			theDistFileNames.push_back(pair<string, string>(thePair.substr(0, p), thePair.substr(p + 1)));
			getline(batchFile, thePair, '\n');
		}
		batchFile.close();
	}

	// MMML extensions: New 2010.02.28: RLC:
	const set<string> ignoreTaxa(identifyTaxaToIgnore(ignoreTaxaFileName));
	
	const pair<TreeType, map<string, set<string> > > theClades(
	retrieveRefCladesFromFile(treeFileName, minBootstrapSupport));
	const TreeType& refSets = theClades.first; // for convenience
	const map<string, set<string> >& refClades = theClades.second;
	// These are empty if we're not in MMML mode...
	// MMML extensions: New 2010.02.28: RLC.

	ofstream outfile, tablefile; // New 2010.01.23: RLC edit
	if (MMMLmode || (brief && outputZeros) || !noHeader) // New 2010.02.28: RLC, added MMMLmode file opening
	{	outfile.open(outfileName.c_str());
	}
	
	// MMML extensions: New 2010.02.28: RLC:
	if (MMMLmode)
	{	if (LverboseMode) writeTree(outfile, refSets, treeFileName);
	}
	// MMML extensions: New 2010.02.28: RLC.
	
	//abezgino 2011.01.21 : print header when in batch mode
	if (!noHeader)
	{	if (MMMLmode) writeHeader(outfile, topCoverage, listTaxa);
		if (brief)
		{	outfile << "#Matrix1\tMatrix2\ttaxonCombinations";
			if (quickMP) outfile <<"\tquickMP";
			outfile << "\tMP";	
			if (!onlyMP) outfile << "\ttime\tscore\tRMSD\ttrees";
			outfile << endl;
		}
			
		if (onlyHeader)
		{	outfile.close();
			exit (0);
		}
	}

	if (useHardware)
	{
		if (!CHardAlgorithm::initHardware())
		{
			cout << "Could not initialize max clique hardware" << endl;
			return -1;
		}
	}

	list <DistanceMatrix> distMatricesCache;
	list <Alignment> alignmentsCache;
	
// Batch loop:
	vector<pair<string, string> >::const_iterator batchIt, batchItEnd = theDistFileNames.end();
	for (batchIt = theDistFileNames.begin(); batchIt != batchItEnd; ++batchIt)
	{
		// Open output file
		if (!outfile.is_open())
		{
			if (rangeStart >= 0)
			{
				string::size_type lastDot = outfileName.find_last_of('.');
				if (lastDot == string::npos)
				{
					outfile.open((outfileName + '_' + itoa(rangeStart) + '-' + itoa(rangeEnd)).c_str());
				}
				else
				{
					outfile.open((outfileName.substr(0, lastDot) + '_' + itoa(rangeStart) + '-' + itoa(rangeEnd) +
						outfileName.substr(lastDot)).c_str());
				}
			}
			else
			{
				outfile.open(outfileName.c_str());
			}
		}

		// Print matrix file names to stdout
		cout << batchIt->first << '\t' << batchIt->second; 

		// Brief format output: start each row with the two matrix file names
		if(brief && outputZeros) 
		{	
			outfile << batchIt->first << '\t' << batchIt->second;
		}

		// Read and process the distance matrices:
		vector<double> distances1, distances2;
		vector<string> names1, names2, taxLabels1, taxLabels2, taxLabels;

		//const int num1 = readDistFile("1:", distfileNamesPrefix + batchIt->first, distances1, names1, taxLabels1, true);
		//const int num2 = readDistFile("2:", distfileNamesPrefix + batchIt->second, distances2, names2, taxLabels2, true);
		
		int num1, num2;
		try
		{	num1 = getDistMatrix("1:", distfileNamesPrefix, batchIt->first, distMatricesCache, alignmentsCache, distances1, names1, taxLabels1, true, aln2pmb);
			num2 = getDistMatrix("2:", distfileNamesPrefix, batchIt->second, distMatricesCache, alignmentsCache, distances2, names2, taxLabels2, true, aln2pmb);
		} catch (exception &e)
		{	cout << "\tError: Could not obtain at least one distance matrix: " << e.what() << endl;
			cerr << "\tError: Could not obtain at least one distance matrix: " << e.what() << endl;
			if(brief && outputZeros) 
			{	outfile << endl;
			}
			continue;
		}
		
		if (reqTax)
		{	rearrangeMatrix(names1, taxLabels1, distances1, reqTaxName);
			//rearrangeMatrix(names2, taxLabels2, distances2, reqTaxName);
			
			//reqTax=false; //debug - check whether simply rearranging the matrices affects the results
		}

		// Assign unique integers to taxon labels
		vector<int> taxIndices1, taxIndices2;
		int reqTaxIndex;
		processTaxonLabels(taxLabels1, reqTaxName, taxIndices1, taxLabels, reqTaxIndex);
		processTaxonLabels(taxLabels2, reqTaxName, taxIndices2, taxLabels, reqTaxIndex);
		
		// showstats: just show information about matrices and do nothing
		if (doShowStats)
		{
			showStats(taxIndices1, taxIndices2, taxLabels);
			continue;
		}
	
		processDistances(distances1, cutoffLower, cutoffUpper, logDistances);
		processDistances(distances2, cutoffLower, cutoffUpper, logDistances);

		// Initialize the algorithm
		CMMMAlgorithm* algo = new CHardAlgorithm(allow, distances1, distances2,
				taxIndices1, taxIndices2, useTaxInfo, reqTaxIndex, 
				reqTax, strictSize, maxTrees, minSize, useHardware);

		// Calculate match potential and taxon combinations
		int taxonCombinations = calcTaxonCombinations(taxIndices1, taxIndices2);
		
		
		//abezgino 2011.01.21: calculate both match potentials (should not add too much overhead)
		int quickMatchPotential = findBiggestPossibleMatch(taxLabels1, taxLabels2, useTaxInfo);
		int matchPotential = algo->calcMatchPotential();
		
		//abezgino 2011.01.21
		//debug code to ensure that new matchPotential is never larger than quickMatchPotential
		if (matchPotential > quickMatchPotential)
		{	matchPotential = quickMatchPotential;
			#ifdef DEBUG
			cout << batchIt->first << '\t' << batchIt->second << "\tError: matchPotential=" << matchPotential << " > quickMatchPotential=" << quickMatchPotential;
			cerr << batchIt->first << '\t' << batchIt->second << "\tError: matchPotential=" << matchPotential << " > quickMatchPotential=" << quickMatchPotential;
			handleErrors(ignoreErrors, outfile);
			#endif
	
		}
		// Brief format output: print taxon combinations and match potential before analysis begins
		if(outputZeros && brief)
		{
			//abezgino: print known information about current combination before the actual analysis. Make taxonCombinations optional? 
			//abezgino 2011.01.21
			outfile << '\t' << taxonCombinations;
			if (quickMP) outfile << '\t' << quickMatchPotential;
			outfile << '\t' << matchPotential;
			if (onlyMP)
			{	cout << endl;
				outfile << endl;
				continue;
			} else
			{	cout << flush;
				outfile << flush;
			}
		}
		
		
		
		// Check that the required taxon name, if specified, is in both matrices:
		if (reqTax && (find(taxLabels1.begin(), taxLabels1.end(), reqTaxName) == taxLabels1.end() ||
							find(taxLabels2.begin(), taxLabels2.end(), reqTaxName) == taxLabels2.end()
							)
			)
		{	if (brief && outputZeros) outfile << "\t0\t0\t0\t0" << endl; //format: time	score	rmsd	trees
			delete algo;
			continue; // New 2010.01.23: RLC bug fix
		}

		// Conduct the analysis:
		Timeval startTime, endTime;
		WeightedProteinPairVecVec resultVec;

		int wkRangeStart = max(0, min(rangeStart, num1));
		int wkRangeEnd = min(rangeEnd, num1);
		int maxScore = 0;

		timer_gettime(&startTime);
	
		if (stateRestore)
		{
			if (!algo->stateRestore(stateFileName))
			{
				cout << "\tError loading state file " << stateFileName << endl;
				exit(1);
			}
		}
		
		string saveFile = makeStateSaveFileName(batchIt->first, batchIt->second);
		sighand_init(algo, saveFile);
		
		try
		{
			maxScore = algo->launchAnalysis(wkRangeStart, wkRangeEnd);	
			
			if (maxScore < 0)
			{
				sighand_restore();
				cout << "\tLast reached index: " << algo->lastReachedIndex() << endl;
				delete algo;
				continue;
			}
			algo->calcScores(resultVec);
			algo->printBonusStats(outfile, brief);
		}
		catch (bad_alloc&)
		{
			cout << "\tOut of memory (bad_alloc)" << endl;
			sighand_restore();
			delete algo;
			continue;
		}
		
		timer_gettime(&endTime);

		double runtime = timer_diff(&startTime, &endTime);
		if (!exactTime)
		{
			runtime = floor(runtime + 0.5);
		}

		cout << "\tTime: " << runtime;

		if (maxTrees >= 0 && algo->maxTreesReached())
		{
			cout << "\tWarning: result was truncated by -maxTrees to " << maxTrees;
		}
		cout << endl;
		
		sighand_restore();
		delete algo;

		if (brief && outputZeros && maxScore == 0) 
		{	
			outfile << '\t' << runtime << "\t0\t0\t0" << endl; //format:time	score	rmsd	trees
			continue;
		}
		
		// Generate the output files:
		if ((strictSize < 3 && maxScore > 0) || (maxScore == strictSize))
		{
			if (maxScore < 3)
			{
				//cout << "Warning: maxScore=" << maxScore << "< 3" << flush;
				cerr << batchIt->first << '\t' << batchIt->second << "\tError: maxScore=" << maxScore << " < 3";
				handleErrors(ignoreErrors, outfile);
			}
			
			const size_t numTrees = resultVec.size();
			
			sort(resultVec.begin(), resultVec.end(), CmpTreesByWeight());
			
			// Debug code: returned trees should not be larger than biggest possible match
			if (maxScore > matchPotential)
			{
				cerr << batchIt->first << '\t' << batchIt->second << "\tError: maxScore=" << maxScore << " > matchPotential=" << matchPotential;
				handleErrors(ignoreErrors, outfile);
			}
			//assert(maxScore <= matchPotential);
			
			// Output summary
			if (brief)
			{
				if (!outputZeros)
				{	
					outfile << batchIt->first << '\t' << batchIt->second << '\t' << taxonCombinations; 
					if (quickMP) outfile << '\t' << quickMatchPotential;
					outfile << '\t' << matchPotential;
				}
				
				double rmsd = 0.0;
				if (!resultVec.empty())
				{
					rmsd = resultVec.front().second;
				}
				outfile << '\t' << runtime << '\t' << maxScore << '\t' << rmsd << '\t' << numTrees << endl;
				if (!tabDelimitedTable) continue; //nothing more to do in brief mode except when in table mode. Do we actually need table mode at all now?
			}
			else if (!MMMLmode) // New 2010.02.28: RLC
			{
				outfile << "infile 1: " << batchIt->first << '\n'; // New 2010.01.23: RLC, to support MMML
				outfile << "infile 2: " << batchIt->second << '\n'; // New 2010.01.23: RLC, to support MMML
				outfile << "maxScore = " << maxScore << " (max possible score = " << matchPotential << ")\n";
				outfile << "number of matching trees = " << numTrees << '\n';
			}

			// MMML extensions: New 2010.02.28: RLC:
			MMMdata* theMMMdata = 0;
			if (MMMLmode) {
				theMMMdata = new MMMdata(batchIt->first, batchIt->second, distfileNamesPrefix, numTrees, ignoreTaxa, refClades, listTaxa);
				if (theMMMdata->bad()) {
					cout << "Error constructing MMMdata object for " << batchIt->first << " - " << batchIt->second << endl;
					delete theMMMdata;
					continue;
				}
			}
			// MMML extensions: New 2010.02.28: RLC.

			// Output the trees
			WeightedProteinPairVecVec::const_iterator vmIt, vmItEnd = resultVec.end();
			ProteinPairVec::const_iterator mIt, mItEnd;
			int t = 1;
			vector<vector<float> > counts(num1, vector<float>(num2));
			for (vmIt = resultVec.begin(); vmIt != vmItEnd; ++vmIt, ++t)
			{
				mItEnd = vmIt->first.end();
				// New 2010.02.28: RLC, modified to feed MMML where desired
				if (MMMLmode)
				{	set<string> theTaxa1, theTaxa2;
					for (mIt = vmIt->first.begin(); mIt != mItEnd; ++mIt)
					{	theTaxa1.insert(names1[mIt->first].substr(2));
						theTaxa2.insert(names2[mIt->second].substr(2));
					}
					theMMMdata->processSubmatricesFromMMM(theTaxa1, theTaxa2, t, (topCoverage ? 0 : &outfile));
				} else
				{	// New 2010.02.28: RLC.
					double theWtScore = vmIt->second;
					if (!brief) outfile << "Submatrix " << t << ": weighted score = " << theWtScore << '\n'; // New 2010.01.23: RLC edit
					for (mIt = vmIt->first.begin(); mIt != mItEnd; ++mIt)
					{	if (!brief) outfile << names2[mIt->second] << "->" << names1[mIt->first] << '\n';
						if (t == 1 || !tabulateTop) counts[mIt->first][mIt->second] += 1.0F;
					}
				} // New 2010.02.28: RLC
			}
			// MMML extensions: New 2010.02.28: RLC:
			if (MMMLmode)
			{	theMMMdata->produceSummaryOutput(outfile);
				delete theMMMdata;
			}
			// MMML extensions: New 2010.02.28: RLC.

			if (tabDelimitedTable)
			{
				if (!tablefile.is_open()) // New 2010.01.23: RLC bug fix: as it was, the file kept getting overwritten
				{
					tablefile.open((outfileName + ".tab").c_str());
				}
				tablefile << numTrees << '\t' << maxScore << "\n----";
				for (t = 0; t < num1; ++t) 
				{
					tablefile << '\t' << names1[t];
				}
				
				const size_t numTreeNorm = tabulateTop ? 1 : numTrees;
				for (int u = 0; u < num2; ++u)
				{
					tablefile << '\n' << names2[u];
					for (t = 0; t < num1; ++t) tablefile << '\t' << counts[t][u] / numTreeNorm;
				}
				tablefile << endl;
				// New 2010.01.23: RLC bug fix: moved closing of file out of loop lest it be repeatedly overwritten
			}
		}
	}

	if (useHardware)
	{
		CHardAlgorithm::closeHardware();
	}

	tablefile.close(); // New 2010.01.23: RLC bug fix
	outfile.close();
	
	return 0;
}
//



static void printUsage1(const string& pname)
{
	cout << "Usage 1: " << pname <<
		" -d1 distfile1 -d2 distfile2 -o outfile -a allowance [-brief] " <<
		"[-table [-top]] [-u] [-req reqTaxName] [-sz strictSize] [-rangeStart start -rangeEnd end] " <<
		"[-lower cutoffLower] [-upper cutoffUpper] [-aln2pmb] " <<
		"[-v] [-restore stateFileName]" << endl;
}
//

static void printUsage2(const string& pname)
{
	cout << "Usage 2: " << pname <<
		" -b batchfile -o outfile -a allowance [-brief] " <<
		"[-u] [-req reqTaxName] [-sz strictSize] [-rangeStart start -rangeEnd end] " <<
		"[-lower cutoffLower] [-upper cutoffUpper] [-aln2pmb] " <<
		"[-v]" << endl;
}
//

static void printUsage(const string& pname)
{
	printUsage1(pname);
	cout << "...or..." << endl;
	printUsage2(pname);
	cout << "...or..." << endl;
	cout << "Usage: " << pname << " -help" << endl;
	cout << "MMML extensions, incompatible with -brief or -z: -Lt rootedTreeFileName [-LtopCov] [-Li ignoreTaxaFileName] " <<
		"[-Lb minBootstrapSupport] [-List] [-Lv]" << endl;
	exit(1);
}
//

static void printVerboseUsage(const string& pname)
{
	cout << endl;
	printUsage1(pname);
	cout << "...or..." << endl;
	printUsage2(pname);
	cout << "\nRequired parameters:\n" << endl;
	cout << "-d1 distfile1 (if -b batchfile is not specified)" << endl;
	cout << "\tSpecifies the file name of the first distance matrix." << endl;
	cout << "-d2 distfile2 (if -b batchfile is not specified)" << endl;
	cout << "\tSpecifies the file name of the second distance matrix." << endl;
	cout << "-b batchfile (if -d1 distfile1 -d2 distfile2 are not specified)" << endl;
	cout << "\tSpecifies the name of a file containing tab-delimited pairs of distance matrix file names, one pair per line." << endl;
	cout << "-o outfile" << endl;
	cout << "\tThe desired name of the output file." << endl;
	cout << "-a allowance" << endl;
	cout << "\tSpecifies a value between 0.0 and 1.0, with smaller values requiring more precise matching, and "
		  << "larger values tolerating looser matching." << endl;
	cout << "\nOptional parameters:\n" << endl;
	cout << "-brief" << endl;
	cout << "\tOutputs a single tab-delimited line of output for the pair of distance matrices." << endl;
	cout << "-table (valid only when -b batchfile is not specified)" << endl;
	cout << "\tProduces a tab-delimited matrix that lists the "
		  << "fit-weighted proportion of matches for each entry from the two distance matrices, "
		  << "over all common submatrices. For example, if seq1 matched seq2 in two of ten discovered submatrices, with "
		  << "a score of 0.8 and 0.9, respectively, then the entry for seq1_seq2 would be (0.8 + 0.9) / 10 = 0.17" << endl;
	cout << "-top (valid only when -table is specified)" << endl;
	cout << "\tRestricts the -table output to just the best-scoring submatrix." << endl;
	cout << "-u" << endl;
	cout << "\tOnly attempts to match entries from the pair of distance matrices that have the same tag in their name. A tag "
		  << "ends with a pipe (\"|\") character, e.g. Ecoli in Ecoli|DnaA, or Hsapiens in Hsapiens|myc" << endl;
	cout << "-req reqTaxName" << endl;
	cout << "\tRequires that a submatrix include at least one entry with the specified reqTaxName, e.g. Hsapiens." << endl;
	cout << "-sz strictSize" << endl;
	cout << "\tReturns common submatrices of size strictSize, if any. Does not attempt to find larger submatrices." << endl;

	cout << "-p distfileNamesPrefix" << endl;
	cout << "\tThe prefix is appended to the path of all distance matrices." << endl;
	cout << "-z" << endl;
	cout << "\tForces the output of zeros for combinations that have not reached the minSize (to be used in 'brief' mode)." << endl;
	cout << "-minsize minSize" << endl;
	cout << "\tReturns common submatrices of at least minSize, if any. Will not return smaller submatricies. Use this option to save memory with large matrices." << endl;
	cout << "-maxtrees maxTrees" << endl;
	cout << "\tLimits the number of returned trees. Saves memory." << endl;
	cout << "-rangeStart start -rangeEnd end" << endl;
	cout << "\tFor distributed computing, launch this program multiple times with non-overlapping ranges, [start, end), "
		  << "start >= 0 and end <= distfile1's matrix size (auto-adjusted when -b is specified). Incompatible with -rndOrder." << endl;
	cout << "-quickmp" << endl;
	cout << "\tIf specified, the quickly but naively estimated match potential is included in batch output." << endl;
	cout << "-v" << endl;
	cout << "\tOutputs the size of the largest common submatrix found so far, to stdout." << endl;
	cout << "-onlymp" << endl;
	cout << "\tDoes not perform MMM analysis, only prints MatchPotential (brief and MMML modes only)." << endl;
	cout << "-restore stateFileName" << endl;
	cout << "\tResumes work from a previously interrupted MMM session" << endl;
	cout << "-exacttime" << endl;
	cout << "\tDo not round running time to nearest second, in brief mode" << endl;
	cout << "-upper cutoffUpper" << endl;
	cout << "\tAll distances above cutoffUpper are set to zero." << endl;
	cout << "-lower cutoffLower" << endl;
	cout << "\tAll distances below cutoffLower are set to zero." << endl;
	cout << "-log" << endl;
	cout << "\tReplaces the distances with their logs: dst=ln(dst)+1 if dst > 1." << endl;
	cout << "-noheader" << endl;
	cout << "\tOmits the header from the output file (brief and MMML modes only)." << endl;
	cout << "-onlyheader" << endl;
	cout << "\tPrints the header to the output file and exits (brief and MMML modes only)." << endl;
	cout << "-ignoreerrors" << endl;
	cout << "\tDoes not abort the run even in the case of program errors (intended for debugging only)." << endl;
	cout << "-aln2pmb" << endl;
	cout << "\tEnable recognition of .aln files as FASTA alignments and processing them into PMB matrices." << endl;
	cout << "\tSimply specify the alignment filename in place of distance matrix file. This works in both command-line and batch mode." << endl;
	cout << "\tA named alignment slice can be defined as comma-separated ranges of columns added after the filename." << endl;
	cout << "\tWhen using multiple ranges to define a slice, the ranges must go in increasing order, and must not overlap." << endl;
	cout << "\tAn \"inverse\"slice (remove the slice itself and keep all other columns) can be specified by adding 'r' after the range. Do not use with multiple ranges." << endl;
	cout << "\tFormat: <filename>|[range name]:<start>-<end>[r][,<start>-<end> ...]" << endl;
	cout << "\tExamples:" << endl;
	cout << "\t\talignment1.aln" << endl;
	cout << "\t\talignment2.aln|pfam12820:373-546" << endl;
	cout << "\t\talignment3.aln|disorder:16-90,150-150,430-519" << endl;
	cout << "\t\talignment4.aln|:300-400r" << endl;
	cout << "MMML extensions:" << endl;
	cout << "\t-Lt rootedTreeFileName (required)" << endl;
	cout << "\t\tFile name for a Newick style tree file encompassing all of the taxa to be encountered in the MMM run." << endl;
	cout << "\t-LtopCov (optional)" << endl;
	cout << "\t\ttop-coverage style output." << endl;
	cout << "\t-Li ignoreTaxaFileName (optional)" << endl;
	cout << "\t\tfile listing taxon names to ignore in the MMML analysis." << endl;
	cout << "\t-Lb minBootstrapSupport (optional, default=1)" << endl;
	cout << "\t\tan integer between 0 and 100 representing the threshold below which to collapse branches in parsing the treefile." << endl;
	cout << "\t-List (optional)" << endl;
	cout << "\t\tif specified, lists the taxa in each clade in the output." << endl;
	cout << "\t-Lv (optional)" << endl;
	cout << "\t\tif specified, outputs the clade sets." << endl;
	exit(1);
}
//

static int readDistFile(const string& which, const string& distfileName, vector<double>& distances,
						vector<string>& names, vector<string>& taxLabels, bool useTaxInfo)
{
	ifstream distfile(distfileName.c_str());
	if (!distfile.is_open()) {
		cout << "\nError: cannot open file " << distfileName << endl;
		exit(1);
	}
	int num;
	distfile >> num;
	names.resize(num);
	taxLabels.resize(num);
	distances.clear();
	distances.reserve(num*num);
	string theName;
	double theDist;
	int i, j;
	for (i = 0; i < num; ++i) {
		distfile >> theName;
		if (useTaxInfo) {
			taxLabels[i] = theName.substr(0, theName.find('|')); // If not found, selects the whole name
		}
		names[i] = which + theName;
		for (j = 0; j < num; ++j) {
			distfile >> theDist;
			distances.push_back(max(1e-20,theDist)); // Don't allow negative or zero distances
		}
	}
	distfile.close();
	return num;
}
//

//abezgino
//rearrange the matrices to ensure that the sequences for the reqTax are placed at the top
static void rearrangeMatrix(vector<string>& names, vector<string>& taxLabels, vector<double>& distances, const string reqTaxName)
{	
	int size = taxLabels.size();
	vector<int> map; //mapping of original indices to the rearranged ones
	map.reserve(size);
			
	unsigned int firstNonReqTaxIdx = 0;
	
	//rearrange taxLabels and names
	for (unsigned int i = 0; i < size; i++)
	{	map.push_back(i);
		if (taxLabels[i] == reqTaxName)
		{	swap(map[i], map[firstNonReqTaxIdx]);
			swap(taxLabels[i], taxLabels[firstNonReqTaxIdx]);
			swap(names[i], names[firstNonReqTaxIdx]);
			firstNonReqTaxIdx++;
		}
	}
	
	//now rearrange distances
	//use the new in-place approach
	for (unsigned int i = 0; i < size; i++)
	{	for (unsigned int j = 0; j < size; j++)
		{	//make sure that each element gets swapped only once
			if(map[i] > i || (map[i] == i && map[j] > j))
			{ swap(distances[i*size+j], distances[map[i]*size+map[j]]);
			}
		}
	}
	
	return;
}
//

//abezgino
//convert taxon labels to integers for easier manipulation
static void processTaxonLabels(vector<string> taxLabels, const string reqTaxName, vector<int>& taxIndices,
										 vector<string>& taxLabelsGlobal, int& reqTaxIndexGlobal)
{	for (unsigned int i = 0; i < taxLabels.size(); i++)
	{
		bool found = false;
		
		unsigned int j;
		for (j = 0; j < taxLabelsGlobal.size(); j++)
		{	if (taxLabels[i] == taxLabelsGlobal[j])
			{	found = true;
				taxIndices.push_back(j);
				break;
			}
		}
		
		if (!found)
		{	//since the above loop was allowed to run to the end then j == taxLabelsGlobal.size();
			taxLabelsGlobal.push_back(taxLabels[i]);
			taxIndices.push_back(j);
			if (reqTaxName == taxLabels[i])
				reqTaxIndexGlobal = j;
		}
	}
}
//

size_t findBiggestPossibleMatch(const vector<string>& taxLabels1, const vector<string>& taxLabels2,
								const bool useTaxInfo)
{
	const size_t num1 = taxLabels1.size();
	const size_t num2 = taxLabels2.size();

	if (!useTaxInfo) return min(num1, num2);

	map<string, pair<size_t, size_t> > taxCounts;
	size_t i;

	for (i = 0; i < num1; ++i)
	{
		++taxCounts[taxLabels1[i]].first;
	}
	for (i = 0; i < num2; ++i)
	{
		++taxCounts[taxLabels2[i]].second;
	}

	size_t answer = 0;
	map<string, pair<size_t, size_t> >::const_iterator mIt, mItEnd = taxCounts.end();
	for (mIt = taxCounts.begin(); mIt != mItEnd; ++mIt)
	{
		answer += min(mIt->second.first, mIt->second.second);
	}
	return answer;
}
//

static string itoa(int number)
{
	ostringstream ost;
	ost << number;
	return ost.str();
}
//

static void showStats(vector<int>& taxIndices1, vector<int>& taxIndices2, vector<string>& taxLabels)
{
	cout << "Matrix 1 size: " << taxIndices1.size() << endl;
	cout << "Matrix 2 size: " << taxIndices2.size() << endl;
	cout << "Number of unique taxons: " << taxLabels.size() << endl;

	unsigned int compat = calcTaxonCombinations(taxIndices1, taxIndices2);
	unsigned int total = taxIndices1.size() * taxIndices2.size();
	float ratio = (float)compat/(float)total;

	cout << "Compatible/total pairings: " << compat << '/' << total << " (" << ratio << ")" << endl;
}
//

static int calcTaxonCombinations(vector<int>& taxIndices1, vector<int>& taxIndices2)
{
	//abezgino: calculate the number of combinations for all taxon names.
	map<int, int> taxCounts1, taxCounts2;
	int taxonCombinations = 0;

	for (vector<int>::iterator it = taxIndices1.begin();  it < taxIndices1.end(); it++ )
		taxCounts1[*it]++;

	for (vector<int>::iterator it = taxIndices2.begin();  it < taxIndices2.end(); it++ )
		taxCounts2[*it]++;

	for (map<int, int>::iterator it = taxCounts1.begin();  it != taxCounts1.end(); it++ )
	{	
		int currentTaxCount1 = it->second;
		int currentTaxCount2 = taxCounts2[it->first];

		taxonCombinations += currentTaxCount1 * currentTaxCount2;
	}

	return taxonCombinations;
}
//

static string makeStateSaveFileName(const string& first, const string& second)
{
	string::size_type slashpos1, slashpos2;
	string result;

	slashpos1 = first.find_last_of("/\\") + 1;
	slashpos2 = second.find_last_of("/\\") + 1;
	result = first.substr(slashpos1) + '_' + second.substr(slashpos2) + ".sav";

	return result;
}
//

//process a raw distance matrixs 
//for now only remove distances above and below specified thresholds
static void processDistances(vector<double>& distances, double cutoffLower, double cutoffUpper, bool logDistances)
{	bool errorReported = false;
	for (unsigned int i = 0; i < distances.size(); i++ )
	{
		//cout << distances[i] <<" \t";
		if (distances[i] < cutoffLower) distances[i] = 1.0e-20;
		if (distances[i] > cutoffUpper) distances[i] = 1.0e-20;
		if (abs(distances[i] - 0.000010) < 1.0e-20 && !errorReported)
		{	cout << "\tWarning: distance of 0.000010 found. Forgot to zero protdist distances?";
			errorReported = true;
		}
		//cout << distances[i] <<"\n";
		if (logDistances && distances[i] > 1)
		{	distances[i] = log(distances[i])+1;
		}
	}
	return;
}
//


static void handleErrors(bool ignoreErrors, ofstream& outfile)
{	if (ignoreErrors)
	{	cerr << "\tIgnored ..." << endl;
	}
	else
	{	cerr << endl;
		cout << endl;
		if (outfile.is_open())
		{	outfile << endl;
			outfile.close();
		}
		exit(1);
	}
	return;
}
//

static unsigned int getDistMatrix(const string& which, const string& distFilenamePrefix, const string& distFilename,
									  list <DistanceMatrix>& distMatricesCache, list <Alignment>& alignmnentsCache,
									  vector <double>& distances, vector <string>& names, vector <string>& taxLabels,
									  bool useTaxInfo, bool aln2pmb)
{
	//parameters
	bool quickPmb=true, printPmb=false, printSlice=false; 
	
	//try to find the matrix in the cache
	if(!find_in_dist_matrices_cache(distMatricesCache, distFilename))
	{
		#ifdef DEBUG
		cout << "\nWarning: matrix [" << distFilename << "] is not cached. Trying to read from file. Cache size=" << distMatricesCache.size() << '/' << DISTANE_MATRIX_CACHE_SIZE << endl;
		#endif
		distMatricesCache.push_front(DistanceMatrix());
		DistanceMatrix& distanceMatrixTmp = distMatricesCache.front();
		distanceMatrixTmp.set_name(distFilename);
		
		
		stringstream distFilenameSS(distFilename);
		string filename;	
	
		getline(distFilenameSS, filename, '|');
		bool subAln = !distFilenameSS.eof();
		
		string extension = filename.substr(filename.find_last_of('.')+1);
		
		if (aln2pmb && extension == "aln")
		{
			#ifdef DEBUG
			cout << "\tWarning: extension [" << extension <<"] found. Calculating PMB distnce matrix from alignment ..." << endl;
			#endif
			
			Alignment const& alignment = getAlignment(distFilenamePrefix, filename, alignmnentsCache);
			#ifdef DEBUG
			cout << "Alignment length: " << alignment.get_nOfColumns() << endl;
			//cout << alignment.get_sequence(1) << endl;
			#endif
			if(subAln)
			{	Alignment alignmentSlice = getSlice(distFilenameSS, alignment);
				#ifdef DEBUG
				cout << "Slice length: " << alignmentSlice.get_nOfColumns() << endl;
				//cout << alignmentSlice.get_sequence(1) << endl;
				#endif
				if(printSlice)
				{	string filepath = distFilenamePrefix+distFilename+".slice";
					ifstream file(filepath.c_str()); 
					if (!file)
						alignmentSlice.print_to_file(filepath);
					else
						file.close();
				}
				distanceMatrixTmp.create_from_alignment(alignmentSlice, quickPmb);
			} else
			{	distanceMatrixTmp.create_from_alignment(alignment, quickPmb);
			}
			
			if (printPmb)
			{	string filepath = distFilenamePrefix+distFilename+".pmbq";
				ifstream file(filepath.c_str()); 
				if (!file)
					distanceMatrixTmp.print_to_file(filepath);
				else
					file.close();
			}
		} else
		{
			#ifdef DEBUG
			if(extension == "aln")
				cout << "\tWarning: extension [" << extension <<"] found. Forgot to add -aln2pmb option?" << endl;
			#endif
			distanceMatrixTmp.read_from_file(distFilenamePrefix+distFilename);
		}
		
		if(distMatricesCache.size() > DISTANE_MATRIX_CACHE_SIZE)
		{	
			#ifdef DEBUG
			cout << "\tWarning: cache has too many elements. Removing the last one [" << distMatricesCache.back().get_name() << "]" << endl;
			#endif
			distMatricesCache.pop_back();
		}
	}
	
	//at this point the first element of cache list must be the matrix.
	
	DistanceMatrix& distanceMatrix = distMatricesCache.front();
	
	//fill the output variables
	unsigned int nOfSequences = distanceMatrix.get_nOfSequences();
	names.resize(nOfSequences);
	taxLabels.resize(nOfSequences);
	for (unsigned int i = 0; i < nOfSequences; i++)
	{
		string header = distanceMatrix.get_header(i);
		if (useTaxInfo)
		{	taxLabels[i] = header.substr(0, header.find('|')); // If not found, selects the entire header
		}
		names.at(i) = which + header;
	}
	distanceMatrix.get_distances_as_1d(distances);
	
	return nOfSequences; //can be zero?
}
//

static Alignment const& getAlignment(const string& inFilepathPrefix, const string& inFilename, list <Alignment>& alignmnentsCache)
{	
	//if(! Alignment::find_in_cache(alignmnentsCache, inFilename))
	if(!find_in_alignments_cache(alignmnentsCache, inFilename))
	{
		#ifdef DEBUG
		cout << "\tWarning: alignment [" << inFilename << "] not cached. Trying to read from file. Cache size=" << alignmnentsCache.size() << '/' << ALIGNMENT_CACHE_SIZE << endl;
		#endif
		alignmnentsCache.push_front(Alignment());
		alignmnentsCache.front().set_name(inFilename);
		alignmnentsCache.front().read_alignment(inFilepathPrefix+inFilename, -1); //even if this fails, the alignment name is already cached, so the file won't be attempted to get re-read several times over.

		if(alignmnentsCache.size() > ALIGNMENT_CACHE_SIZE)
		{
			#ifdef DEBUG
			cout << "\tWarning: cache has too many elements. Removing the last one [" << alignmnentsCache.back().get_name() << "]" << endl;
			#endif
			alignmnentsCache.pop_back();
		}
	}
	//at this point the front element of cache list must be the alignment of interest.
	
	return alignmnentsCache.front(); 
}
//


static bool find_in_alignments_cache(list <Alignment>& cache, string const& name)
{	bool found = false;
	unsigned int i = 0;
	for(list <Alignment>::iterator it = cache.begin(); it != cache.end() && !found; it++)
	{	i++;
		if (it->get_name() == name)
		{	found = true;
			#ifdef DEBUG
			cout << "\talignment [" << name << "] is cached at position " << i << " of " << cache.size() << endl;
			#endif
			//move this element to the front (as the most recently used)
			cache.splice(cache.begin(), cache, it);
		}
	}
	return found;
}
//

static bool find_in_dist_matrices_cache(list <DistanceMatrix>& cache, string const& name)
{	bool found = false;
	unsigned int i = 0;
	for(list <DistanceMatrix>::iterator it = cache.begin(); it != cache.end() && !found; it++)
	{	i++;
		if (it->get_name() == name)
		{	found = true;
			#ifdef DEBUG
			cout << "\tmatrix [" << name << "] is cached at position " << i << " of " << cache.size() << endl;
			#endif
			//move this element to the front (as the most recently used)
			cache.splice(cache.begin(), cache, it);
		}
	}
	return found;
}
//



static Alignment getSlice(stringstream& distFilenameSS, Alignment const& alignment)
{	string name;
	getline(distFilenameSS, name, ':');
	
	#ifdef DEBUG
	cout << "\tWarning: a sub-alignment found: name=[" << name << "]:";
	#endif
	
	vector<pair<int, int > > ranges;
	
	while (!distFilenameSS.eof())
	{	int start=0, end=0;
		distFilenameSS >> start;
		if(distFilenameSS.peek() != '-')
			throw runtime_error("could not parse range for sub-alignment: " + name);
		distFilenameSS.ignore(1);
		distFilenameSS >> end;

		#ifdef DEBUG
		cout << " [" << start <<"]-[" << end << "]";
		#endif
		
		if(!distFilenameSS.eof() && distFilenameSS.peek() == 'r') //"reverse-slice" mode, where only the slice itself is removed
		{	distFilenameSS.ignore(1);
			#ifdef DEBUG
			cout << 'r';
			#endif
			if (start > 1)
				ranges.push_back(pair<int, int>(1,start-1));
			if (end < alignment.get_nOfColumns())
				ranges.push_back(pair<int, int>(end+1,alignment.get_nOfColumns()));
		} else
		{	ranges.push_back(pair<int, int>(start,end));

		}
		
		if(!distFilenameSS.eof() && distFilenameSS.peek() == ',')
			distFilenameSS.ignore(1);
	}
	#ifdef DEBUG
	cout << endl;
	#endif
	if (distFilenameSS.fail() || ranges.size() < 1)
		throw runtime_error("could not parse sub-alignment: " + name);
	return alignment.get_slices(ranges);
}
//
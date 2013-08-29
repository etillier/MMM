#ifndef _MMM_ALGO_H_
#define _MMM_ALGO_H_

#include <vector>
#include <string>
#include <map>
#include <list>
#include <cmath>
#include <queue>
using namespace std;

typedef pair<int,int> ProteinPair;
typedef vector<ProteinPair> ProteinPairVec;
typedef vector<ProteinPairVec> ProteinPairVecVec;

typedef pair<ProteinPairVec,double> WeightedProteinPairVec;
typedef vector<WeightedProteinPairVec> WeightedProteinPairVecVec;


class CMMMAlgorithm
{
public:
	CMMMAlgorithm() { }

	CMMMAlgorithm
	(
		const double allow,
		const vector<double>& distances1, const vector<double>& distances2,
		const vector<int>& taxLabels1, const vector<int>& taxLabels2, const bool useTaxInfo,
		const int& reqTaxName, const bool reqTax, 
		const int strictSize,
		const int minSize,
		const int maxTrees
	)
		: _allow(allow), _distances1(distances1), _distances2(distances2),
		_taxLabels1(taxLabels1), _taxLabels2(taxLabels2), _useTaxInfo(useTaxInfo),
		_reqTaxName(reqTaxName), _reqTax(reqTax),
		_strictSize(strictSize), _minSize(minSize), _maxTrees(maxTrees)
	{
		_num1 = (int)sqrt((float)distances1.size());
		_num2 = (int)sqrt((float)distances2.size());
	}

	virtual ~CMMMAlgorithm() { };
	
	virtual int launchAnalysis
	(
		const int wkRangeStart,
		const int wkRangeEnd
	) = 0;

	virtual void calcScores(WeightedProteinPairVecVec& out) = 0;
	virtual bool stateRestore(string& fileName) = 0;
	virtual void stateSave(string& fileName) = 0;
	virtual bool maxTreesReached() = 0;
	virtual int lastReachedIndex() = 0;
	virtual int calcMatchPotential() = 0;
	virtual void printBonusStats(ostream& str, bool brief) = 0;

protected:
	inline double getDist1(unsigned int a, unsigned int b)
	{
		return _distances1[b*_num1+a];
	}

	inline double getDist2(unsigned int w, unsigned int x)
	{
		return _distances2[x*_num2+w];
	}

	double _allow;
	int _num1;
	int _num2;
	vector<double> _distances1;
	vector<double> _distances2;
	vector<int> _taxLabels1;
	vector<int> _taxLabels2;
	bool _useTaxInfo;
	int _reqTaxName;
	bool _reqTax;
	int _strictSize;
	int _minSize;
	int _maxTrees;
};



class CHardAlgorithm : public CMMMAlgorithm
{
public:
	CHardAlgorithm
	(
		const double allow,
		const vector<double>& distances1, const vector<double>& distances2,
		const vector<int>& taxLabels1, const vector<int>& taxLabels2, const bool useTaxInfo,
		const int& reqTaxName, const bool reqTax, 
		const int strictSize, const int maxTrees, const int minSize,
		const bool useHW
	);

	~CHardAlgorithm();

	static bool initHardware();
	static void closeHardware();

	virtual int launchAnalysis
	(
		const int wkRangeStart,
		const int wkRangeEnd
	);

	virtual void calcScores(WeightedProteinPairVecVec& out);
	virtual bool stateRestore(string& fileName);
	virtual void stateSave(string& fileName);
	virtual bool maxTreesReached();
	virtual int lastReachedIndex();
	virtual int calcMatchPotential();
	virtual void printBonusStats(ostream& str, bool brief);

protected:
	// Structures
	struct CNLEntry
	{
		int dest_a;
		int dest_b;
		double ratio;
	};

	struct CProblemDesc
	{
		int j;
		int start_idx;
	};

	struct CSortEntries
	{
		bool operator() (const CNLEntry& left, const CNLEntry& right)
		{
			return left.ratio < right.ratio;
		}
	};

	struct CSortByDegree
	{
		CSortByDegree(const vector<int>& degs) : _degs(degs) { }
		const vector<int>& _degs;
		bool operator () (const int& l, const int& r) { return (_degs[l] > _degs[r]); }
	};

	typedef vector<CNLEntry> CNumberLine;

	// Functions
	bool validProteinPair(int a, int b);
	bool satisfiesReqTax(int a);
	void buildNumberLine(int i_a, int i_b);
	void buildProblemSet(int i_a, int i_b);
	bool forwardCompatible(int i_a, int j_a, int j_b, int k_a, int k_b, 
		const double& r_min, const double& r_max);
	bool buildMatrix(int min_a, const double& r_min, const double& r_max);
	bool sortVertices();
	void appendDataBuffer(int* header_idx);
	void maxCliquesEmu();
	void maxCliquesHard();
	void copyCliques(int i_a, int i_b);
	void dumpSet(ofstream& file);

	// Member variables
	ProteinPairVecVec _resultVec;
	int _maxScore;
	double _D;
	bool _noTrees;
	bool _useHW;
	int _probNumColorClasses;

	// Pre-allocated data structures used in algorithm
	CNumberLine _numberLine;
	vector<CProblemDesc> _setProbs;
	vector<int> _setVerts;
	vector<int> _probVerts;
	vector<unsigned int> _setInputBuffer;
	vector<unsigned short> _setOutputBuffer;
	vector<char> _probMatrix;
	vector<int> _probOrder;
	vector<int> _probDegree;
	vector<int> _probColorClass;
	vector<int> _probColorClassSize;

	// Used in software max clique algorithm
	vector<int> _mcMSPV;
	vector<int> _mcVStack;
	vector<int> _mcCliqueStack;
	vector<int> _mcPtrStack;
};




#endif

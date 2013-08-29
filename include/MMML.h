// MMML.h
// Version 2010.02.28
// (c) 2010, Robert L. Charlebois and Elisabeth R.M. Tillier

#include "classDistanceMatrixFile.h"
#include "treeManip.h"
#include <iosfwd>

typedef float MMML_float_type;

class MMML_DistanceMatrixFile : public DistanceMatrixFile<MMML_float_type> {
private:
public:
	MMML_DistanceMatrixFile(const std::string& matrixFileName) :
		DistanceMatrixFile<MMML_float_type>(matrixFileName) { }
	~MMML_DistanceMatrixFile() { }
	
	std::set<std::string> taxonNames(const std::set<std::string>& toIgnore) const;
	std::set<std::string> taxonNamesFromSeqsEquivalentTo(const std::string& label) const;
};

inline std::string stripPath(const std::string& withPath)
{
	const std::string::size_type p = withPath.find_last_of('/');
	return (p == std::string::npos) ? withPath : withPath.substr(p + 1);
}

class MMMdata {
private:
	const std::map<std::string, std::set<std::string> >& refClades_;
	const std::set<std::string>& ignoreTaxa_;
	std::vector<MMML_DistanceMatrixFile> matrix_;
	std::pair<std::set<std::string>, std::set<std::string> > union_mTaxa_;
		// union of the list of taxa with and without zero-distance equivalents for each submatrix
	std::pair<std::set<std::string>, std::set<std::string> > top0Coverage_mTaxa_;
		// the list of taxa with and without zero-distance equivalents for the submatrix having the best zero-extended coverage
	std::set<std::string> matchPotentialTaxa_;
	std::string matrixInfoLines_;
	std::string matchPotentialTaxonOutput_;
	float top0Coverage_;
	std::size_t top0CoverageSubmatrixIdx_;
	const std::size_t numTrees_;
	bool listTaxa_;
	bool bad_;
	
	void generateOutputText(void);
public:
	MMMdata(const std::string& sourceFile1, const std::string& sourceFile2, const std::string& distMatrixFilePath, 
		const std::size_t numTrees, const std::set<std::string>& ignoreTaxa,
		const std::map<std::string, std::set<std::string> >& refClades, bool listTaxa);
	~MMMdata() { }

	bool bad() const { return bad_; }
	const std::pair<std::set<std::string>, std::set<std::string> >& itsUnionTaxa() const { return union_mTaxa_; }
	const std::pair<std::set<std::string>, std::set<std::string> >& itsTopCoverageTaxa() const { return top0Coverage_mTaxa_; }
	float itsTopCoverage() const { return top0Coverage_; }
	std::size_t itsTopCoverageSubmatrixIdx() const { return top0CoverageSubmatrixIdx_; }
	const std::string& itsMatrixInfoLines() const { return matrixInfoLines_; }
	const std::string& itsMatchPotentialTaxonOutput() const { return matchPotentialTaxonOutput_; }
	const std::set<std::string>& itsMatchPotentialTaxa() const { return matchPotentialTaxa_; }
	std::set<std::string> matrixTaxonNames(std::size_t idx, const std::set<std::string>& toIgnore) const {
		return matrix_.at(idx).taxonNames(toIgnore);
	}
	std::size_t matrixSize(std::size_t idx) const {
		return matrix_.at(idx).size();
	}
	const std::string& matrixName(std::size_t idx) const {
		return matrix_.at(idx).fileName();
	}
	std::string matrixNameNoPath(std::size_t idx) const {
		return stripPath(matrix_.at(idx).fileName());
	}

	//bool failsSanityCheck(const std::string& mmmFileName);
	//void processSubmatrices(std::ifstream& infile, std::string& theLine);
	//void processSubmatricesWithOutput(std::ifstream& infile, std::ofstream& outfile, std::string& theLine);
	void processSubmatricesFromMMM(const std::set<std::string>& theSet1, const std::set<std::string>& theSet2,
		int subtreeIdx, std::ofstream* outfile);
	void produceSummaryOutput(std::ofstream& outfile);	
	//std::pair<std::set<std::string>, std::set<std::string> > getNextSubtreeFromFile(std::ifstream& infile,
	//	bool* done, std::string& theLine);
};

void writeTree(std::ofstream& outfile, const TreeType& refSets, const std::string& treeFileName);
void writeHeader(std::ofstream& outfile, bool topCoverage, bool listTaxa);
std::string cladeRestriction(const std::map<std::string, std::set<std::string> >& refClades, const std::set<std::string>& taxa, bool listMembers);
std::set<std::string> identifyTaxaToIgnore(const std::string& ignoreTaxaFileName);
std::pair<TreeType, std::map<std::string, std::set<std::string> > >
	retrieveRefCladesFromFile(const std::string& treeFileName, int minBootstrapSupport);

// MMML.cpp
// Version 2010.03.01
// (c) 2010, Robert L. Charlebois and Elisabeth R.M. Tillier

#include "MMML.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>

static const std::string infile1Label("infile 1: ");
static const std::string infile2Label("infile 2: ");
static const std::string submatrixLabel("Submatrix ");
	
std::map<std::string, std::set<std::string> >::const_iterator
	commonCladeCode(const std::map<std::string, std::set<std::string> >& refClades, const std::set<std::string>& taxa);
float computeCladeCoverage(const std::map<std::string, std::set<std::string> >& refClades, const std::set<std::string>& taxa);

std::set<std::string> MMML_DistanceMatrixFile::taxonNames(const std::set<std::string>& toIgnore) const
{
	std::set<std::string> theNames;
	std::vector<std::string>::const_iterator nIt, nItEnd = names_.end();
	for (nIt = names_.begin(); nIt != nItEnd; ++nIt) {
		const std::string taxon(nIt->substr(0, nIt->find('|')));
		if (toIgnore.find(taxon) == toIgnore.end()) {
			theNames.insert(taxon);
		}
	}
	return theNames;
}

std::set<std::string> MMML_DistanceMatrixFile::taxonNamesFromSeqsEquivalentTo(const std::string& label) const
{
	static const MMML_float_type nearZero = 1.0e-10;
	std::set<std::string> answer;
	answer.insert(label.substr(0, label.find('|')));
		// Instead of relying on it being discovered as zero self-distance below
	std::vector<std::string>::const_iterator nIt = std::find(names_.begin(), names_.end(), label); // to find index of label
	if (nIt != names_.end()) {
		const std::size_t idx = std::distance(names_.begin(), nIt);
		// Now find all of the names with zero distance to this sequence and extract their taxon labels:
		const std::vector<MMML_float_type>& v = distances_[idx];
		std::size_t numElements = v.size();
		for (std::size_t i = 0; i < numElements; ++i) {
			if (i != idx && v[i] < nearZero && v[i] >= 0.0) answer.insert(names_[i].substr(0, names_[i].find('|')));
		}
	} else {
		std::cout << "Warning: " << label << " not found in the distance matrix!" << std::endl;
		for (nIt = names_.begin(); nIt != names_.end(); ++nIt) {
			std::cout << *nIt << std::endl;
		}
		std::exit(1);
	}
	return answer;
}

MMMdata::MMMdata(const std::string& sourceFile1, const std::string& sourceFile2, const std::string& distMatrixFilePath, 
 const size_t numTrees, const std::set<std::string>& ignoreTaxa, const std::map<std::string, std::set<std::string> >& refClades,
 bool listTaxa) :
	refClades_(refClades),
	ignoreTaxa_(ignoreTaxa),
	matchPotentialTaxonOutput_("\t\t\t"),
	top0Coverage_(0.0F),
	top0CoverageSubmatrixIdx_(0),
	numTrees_(numTrees),
	listTaxa_(listTaxa),
	bad_(false)
{
	matrix_.push_back(MMML_DistanceMatrixFile(distMatrixFilePath + sourceFile1));
	matrix_.push_back(MMML_DistanceMatrixFile(distMatrixFilePath + sourceFile2));
	bad_ |= (!matrix_[0].hasData() || !matrix_[1].hasData());
	if (bad_) {
		if (!matrix_[0].hasData()) std::cout << "Error: problem with " << sourceFile1 << std::endl;
		if (!matrix_[1].hasData()) std::cout << "Error: problem with " << sourceFile2 << std::endl;
	} else {
		generateOutputText();
	}
}

void MMMdata::generateOutputText(void)
{
	// Generate matrixInfoLines_ and matchPotentialOutput_ for output
	const std::set<std::string> sourceFile1Taxa(matrixTaxonNames(0, ignoreTaxa_));
	const std::set<std::string> sourceFile2Taxa(matrixTaxonNames(1, ignoreTaxa_));
	if (listTaxa_) matchPotentialTaxonOutput_ += "\t\t";
		// default tabs in case the data aren't available
	const std::string matrix1output = cladeRestriction(refClades_, sourceFile1Taxa, false);
	const std::string matrix2output = cladeRestriction(refClades_, sourceFile2Taxa, false);
	std::set_intersection(sourceFile1Taxa.begin(), sourceFile1Taxa.end(), sourceFile2Taxa.begin(),
		sourceFile2Taxa.end(), std::inserter(matchPotentialTaxa_, matchPotentialTaxa_.begin()));
	std::map<std::string, std::set<std::string> >::const_reverse_iterator cItBegin = refClades_.rbegin();
	if (matchPotentialTaxa_.empty()) {
		std::cout << "***Error*** No taxa in common between the two source distance matrices" << std::endl;
	} else if (!std::includes(cItBegin->second.begin(), cItBegin->second.end(), matchPotentialTaxa_.begin(), matchPotentialTaxa_.end())) {
		std::cout << "***Error*** Taxa in common between the two source distance matrices are not in the supplied tree" << std::endl;
	} else {
		matchPotentialTaxonOutput_ = '\t' + cladeRestriction(refClades_, matchPotentialTaxa_, listTaxa_);
	}
	std::ostringstream ost;
	ost << //tabFlankedInfileNameNoPath has been removed from the output. Can now discard altogether?
		matrixNameNoPath(0) << '\t' << matrixSize(0) << '\t' << matrix1output << '\t' <<
		matrixNameNoPath(1) << '\t' << matrixSize(1) << '\t' << matrix2output;
	matrixInfoLines_ = ost.str();
}

void MMMdata::processSubmatricesFromMMM(const std::set<std::string>& theSet1, const std::set<std::string>& theSet2,
	int subtreeIdx, std::ofstream* outfile)
{
	std::set<std::string> theSet, theSetZeroExt;
	std::set<std::string>::const_iterator sIt, sItEnd = theSet1.end();
	for (sIt = theSet1.begin(); sIt != sItEnd; ++sIt) {
		theSet.insert(sIt->substr(0, sIt->find('|')));
		const std::set<std::string> zeroExtended(matrix_.at(0).taxonNamesFromSeqsEquivalentTo(*sIt));
		theSetZeroExt.insert(zeroExtended.begin(), zeroExtended.end());
	}
	sItEnd = theSet2.end();
	for (sIt = theSet2.begin(); sIt != sItEnd; ++sIt) {
		theSet.insert(sIt->substr(0, sIt->find('|')));
		const std::set<std::string> zeroExtended(matrix_.at(1).taxonNamesFromSeqsEquivalentTo(*sIt));
		theSetZeroExt.insert(zeroExtended.begin(), zeroExtended.end());
	}
	std::set<std::string>::const_iterator igIt, igItEnd = ignoreTaxa_.end();
	for (igIt = ignoreTaxa_.begin(); igIt != igItEnd; ++igIt) {
		theSet.erase(*igIt);
		theSetZeroExt.erase(*igIt);
	}
	if (!theSet.empty()) {
		union_mTaxa_.first.insert(theSet.begin(), theSet.end());
		union_mTaxa_.second.insert(theSetZeroExt.begin(), theSetZeroExt.end());
		if (outfile) {
			*outfile << subtreeIdx << '\t' << numTrees_ << '\t' << matrixInfoLines_ << '\t';
			*outfile << cladeRestriction(refClades_, theSet, listTaxa_) << '\t';
			*outfile << cladeRestriction(refClades_, theSetZeroExt, listTaxa_) << matchPotentialTaxonOutput_ << std::endl;
		} else {
			float theCoverage = computeCladeCoverage(refClades_, theSetZeroExt);
			if (theCoverage > top0Coverage_) {
				top0Coverage_ = theCoverage;
				top0Coverage_mTaxa_.first = theSet;
				top0Coverage_mTaxa_.second = theSetZeroExt;
				top0CoverageSubmatrixIdx_ = static_cast<std::size_t>(subtreeIdx);
			}
		}
	}
}

void MMMdata::produceSummaryOutput(std::ofstream& outfile)
{
	if (top0CoverageSubmatrixIdx_ > 0) {
		outfile << top0CoverageSubmatrixIdx_ << '\t' << numTrees_ <<  '\t' << matrixInfoLines_ << '\t';
		outfile << cladeRestriction(refClades_, top0Coverage_mTaxa_.first, listTaxa_) << '\t';
		outfile << cladeRestriction(refClades_, top0Coverage_mTaxa_.second, listTaxa_) << matchPotentialTaxonOutput_ << '\t';
		outfile << cladeRestriction(refClades_, union_mTaxa_.first, listTaxa_) << '\t';
		outfile << cladeRestriction(refClades_, union_mTaxa_.second, listTaxa_) << '\t';
		outfile << top0Coverage_ << std::endl;
	} else  {
		outfile << "Union" << '\t' << numTrees_ <<  '\t' << matrixInfoLines_ << '\t';
		outfile << cladeRestriction(refClades_, union_mTaxa_.first, listTaxa_) << '\t';
		outfile << cladeRestriction(refClades_, union_mTaxa_.second, listTaxa_) << matchPotentialTaxonOutput_ << std::endl;
	}
}


void writeTree(std::ofstream& outfile, const TreeType& refSets, const std::string& treeFileName)
{
	std::cout << "\nThe reference tree " << treeFileName << std::endl;
	outfile << "The reference tree " << treeFileName << std::endl;
	printSets(refSets, outfile);
	outfile << "\n\n";
}

void writeHeader(std::ofstream& outfile, bool topCoverage, bool listTaxa)
{
	outfile << "#Submatrix\tTotal submatrices\t";
	outfile << "Matrix 1 name\tMatrix 1 size\tMatrix 1 clade|included|total\t";
	outfile << "Matrix 2 name\tMatrix 2 size\tMatrix 2 clade|included|total\t";
	outfile << "Intersection submatrix clade|included|total\t";
	if (listTaxa) outfile << "Clade members included\tClade members excluded\t";
	outfile << "Intersection submatrix clade (zero-extended)|included|total\t";
	if (listTaxa) outfile << "Clade members included\tClade members excluded\t";
	outfile << "Match potential clade|included|total";
	if (listTaxa) outfile << "\tClade members included\tClade members excluded";
	if (topCoverage) {
		outfile << "\tUnion submatrix clade|included|total";
		if (listTaxa) outfile << "\tClade members included\tClade members excluded";
		outfile << "\tUnion submatrix clade (zero-extended)|included|total";
		if (listTaxa) outfile << "\tClade members included\tClade members excluded";
		outfile << "\tZero-extended coverage";
	}
	outfile << std::endl;
}

std::map<std::string, std::set<std::string> >::const_iterator
	commonCladeCode(const std::map<std::string, std::set<std::string> >& refClades, const std::set<std::string>& taxa)
{
	std::map<std::string, std::set<std::string> >::const_iterator bestIt, cIt, cItBegin = refClades.begin(), cItEnd = refClades.end();
	// Find the smallest subtree that includes all of the taxa in matchPotentialTaxa:
	const std::size_t taxaSize = taxa.size();
	std::size_t cladeSize, smallest = 2000000000;
	for (bestIt = cItBegin, cIt = cItBegin; cIt != cItEnd; ++cIt) {
		cladeSize = cIt->second.size();
		if (cladeSize < smallest && cladeSize >= taxaSize && std::includes(cIt->second.begin(), cIt->second.end(), taxa.begin(), taxa.end())) {
			smallest = cladeSize;
			bestIt = cIt;
		}
	}
	if (smallest == 2000000000) {
		std::cout << "***Error: taxa detected that are not in the tree" << std::endl;
	}
	return bestIt;
}

float computeCladeCoverage(const std::map<std::string, std::set<std::string> >& refClades, const std::set<std::string>& taxa)
{
	std::map<std::string, std::set<std::string> >::const_iterator bestIt = commonCladeCode(refClades, taxa);
	return static_cast<float>(taxa.size()) / static_cast<float>(bestIt->second.size());
}

std::string cladeRestriction(const std::map<std::string, std::set<std::string> >& refClades, const std::set<std::string>& taxa, bool listMembers)
{
	std::map<std::string, std::set<std::string> >::const_iterator bestIt = commonCladeCode(refClades, taxa);
	std::string included, excluded;
	if (listMembers) {
		std::set<std::string>::const_iterator sIt, sItEnd = bestIt->second.end();
		for (sIt = bestIt->second.begin(); sIt != sItEnd; ++sIt) {
			if (taxa.find(*sIt) != taxa.end()) {
				if (!included.empty()) included += ", ";
				included += *sIt;
			} else {
				if (!excluded.empty()) excluded += ", ";
				excluded += *sIt;
			}
		}
	}
	std::ostringstream ost;
	ost << bestIt->first << '|' << taxa.size() << '|' << bestIt->second.size();
	if (listMembers) ost << '\t' << included << '\t' << excluded;
	return ost.str();
}

std::set<std::string> identifyTaxaToIgnore(const std::string& ignoreTaxaFileName)
{
	std::set<std::string> ignoreTaxa;
	if (!ignoreTaxaFileName.empty()) {
		std::ifstream ignoreFile(ignoreTaxaFileName.c_str());
		if (!ignoreFile.is_open()) {
			std::cout << "Cannot open " << ignoreTaxaFileName << std::endl;
			std::exit(1);
		}
		std::string ignoreTaxon;
		ignoreFile >> ignoreTaxon;
		while (ignoreFile) {
			ignoreTaxa.insert(ignoreTaxon);
			ignoreFile >> ignoreTaxon;
		}
	}
	return ignoreTaxa;
}

std::pair<TreeType, std::map<std::string, std::set<std::string> > >
	retrieveRefCladesFromFile(const std::string& treeFileName, int minBootstrapSupport)
{
	std::pair<TreeType, std::map<std::string, std::set<std::string> > > theClades;
	if (!treeFileName.empty()) {
		// Read the reference tree file into the string refTree:
		std::ifstream treeFile(treeFileName.c_str());
		if (!treeFile.is_open()) {
			std::cout << "Cannot open tree file " << treeFile << std::endl;
			std::exit(1);
		}
		const std::string refTree((std::istreambuf_iterator<char>(treeFile)), std::istreambuf_iterator<char>());
		treeFile.close();

		// Parse the reference tree into clade sets:
		TreeType& refSets = theClades.first;
		parseTree(refTree, refSets, 'r');
		// Collapse nodes with less support than minBootstrapSupport:
		collapseWeakNodes(refSets, minBootstrapSupport); // default collapses nodes with zero support.
		// Generate the clade lists:
		std::map<std::string, std::set<std::string> >& refClades = theClades.second;
		TreeType::const_iterator tIt, tItEnd = refSets.end();
		for (tIt = refSets.begin(); tIt != tItEnd; ++tIt) {
			refClades[tIt->first] = elaborateSet(tIt->first, refSets);
		}
	}
	return theClades;
}

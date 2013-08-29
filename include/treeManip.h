// treeManip.h
// Version 2010.01.09
// (c) 2004-2005, NeuroGadgets Inc.
// (c) 2004-2010, Author: Robert L. Charlebois

// Building and improving upon HorizStory code:

// MacLeod, D., R.L. Charlebois, W.F. Doolittle and E. Bapteste. 2005. Deduction of probable events of lateral gene
// transfer through comparison of phylogenetic trees by recursive consolidation and rearrangement. BMC Evolutionary
// Biology 5, 27.

#ifndef PHYLOGENETIC_TREE_MANIP_H
#define PHYLOGENETIC_TREE_MANIP_H

#include <ctime>
#include <iosfwd>
#include <map>
#include <set>
#include <string>
#include <vector>

struct TreeBranches {
	std::set<std::string> branches_;
	int bootstrapValue_;
	double branchLength_;
	
	TreeBranches() { }
	TreeBranches(const std::set<std::string>& theBranches, int theBootValue, double theBranchLength) :
		branches_(theBranches),
		bootstrapValue_(theBootValue),
		branchLength_(theBranchLength)
		{ }
};

struct StringTetrad {
	std::string first;
	std::string second;
	std::string third;
	std::string fourth;
	StringTetrad(const std::string& a, const std::string& b, const std::string& c, const std::string& d) :
		first(a), second(b), third(c), fourth(d) { }
};

struct PrintPathOutput {
	int a;
	int b;
	int c;
	PrintPathOutput() : a(0), b(0), c(0) { }
};

typedef std::map<std::string, TreeBranches> TreeType;
typedef std::pair<TreeType, TreeType> TreePair;

void recursiveRearrange(TreePair theTrees, std::size_t bestScore, std::vector<StringTetrad> thePath,
	std::vector<std::vector<StringTetrad> >& vecPaths, std::size_t& shortestPathSoFar,
	std::time_t startTime, std::time_t timeLimit);
PrintPathOutput printPath(const std::vector<StringTetrad>& thePath,
	std::map<std::string, int>& theSummary, std::ofstream& outFile, bool verboseMode);
std::string computePartnerLGT(const std::string& s, const std::string& from,
	const std::string& to, TreeType& theTree, TreeType& theRefTree);
std::string writeTree(const std::string& origTree, const std::string& s, TreeType& theTree);
std::string writeTreeRecursor(const std::string& origTree, const std::string& s, TreeType& theTree);
std::string elaborateSetName(const std::string& s, TreeType& theTree);
std::set<std::string> elaborateSet(const std::string& s, TreeType& theTree);
void printSets(const TreeType& theSets, std::ofstream& outFile);
void printVerticalSets(const TreeType& theSets, std::map<std::string, int>& theSummary,
	std::ofstream& outFile, bool verboseMode);
void findVerticalSetsInString(const std::string& theString, std::set<std::string>& theSets);
long parseTree(std::string t, TreeType& m, char treeType);
void collapseWeakNodes(TreeType& theSets, const int minSupport);
void pruneSequence(TreeType& theTree, const std::string& toRemove);
void consolidateRefTopologies(TreeType& theRefTree, TreeType& theTestTree);
void consolidateSimilarities(TreeType& theRefTree, TreeType& theTestTree);
void moveItem(TreeType& theTree, const std::string& toKey, const std::string& fromKey,
	const std::string& theItem);
std::size_t numLeavesInTree(const TreeType& theTree);
bool isLeafSequence(const std::string& theName);
void trimMultipleUnderscores(std::string& theString);
void findLeafSets(const TreeType& theTree, std::set<std::string>& theLeafSet);

void printOneSet(const TreeType& theSets, const std::string& setName, std::ofstream& outFile);
int computeNumberOfNodesFromRoot(TreeType& theTree, const std::string& name, int count);

#endif

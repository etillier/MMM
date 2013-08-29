// treeManip.cpp
// Version 2010.01.09
// (c) 2004-2005, NeuroGadgets Inc.
// (c) 2004-2010, Author: Robert L. Charlebois

// Building and improving upon HorizStory code:

// MacLeod, D., R.L. Charlebois, W.F. Doolittle and E. Bapteste. 2005. Deduction of probable events of lateral gene
// transfer through comparison of phylogenetic trees by recursive consolidation and rearrangement. BMC Evolutionary
// Biology 5, 27.

#include "treeManip.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cctype>

void recursiveRearrange(TreePair theTrees, std::size_t bestScore, std::vector<StringTetrad> thePath,
	std::vector<std::vector<StringTetrad> >& vecPaths, std::size_t& shortestPathSoFar,
	std::time_t startTime, std::time_t timeLimit)
{
	if (std::time(0) - startTime > timeLimit) return;
	if (bestScore <= 2) {
		vecPaths.push_back(thePath);
		if (thePath.size() < shortestPathSoFar) {
			std::cout << "\n  shortest path length found so far: " << thePath.size() << std::endl;
			std::cout << "   ." << std::flush;
		} else if (thePath.size() == shortestPathSoFar) std::cout << '.' << std::flush;
		shortestPathSoFar = std::min(shortestPathSoFar, thePath.size());
		return;
	}
	if (thePath.size() >= shortestPathSoFar) return; // Optimization.
	const std::size_t lastBestScore = bestScore;
	size_t theScore;
	TreeType tempTree, tempRefTree;
	TreeType::iterator a, b;
	std::set<std::string>::iterator s;
	std::vector<std::pair<TreeType, TreeType> > vecBest;
	std::vector<StringTetrad> vecMoveFromTo;
	std::vector<std::size_t> vecScore;
	TreeType& refSets = theTrees.first;
	TreeType& testSets = theTrees.second;
	for (a = testSets.begin(); a != testSets.end(); ++a) {
		for (b = testSets.begin(); b != testSets.end(); ++b) {
			for (s = b->second.branches_.begin(); s != b->second.branches_.end(); ++s) {
				if (isLeafSequence(*s)) {
					tempTree = testSets;
					tempRefTree = refSets;
					moveItem(tempTree, a->first, b->first, *s);
					if (tempTree.size() > 1UL) {
						consolidateSimilarities(tempRefTree, tempTree);
						theScore = numLeavesInTree(tempTree);
					} else {
						theScore = 1;
					}
					if (theScore < lastBestScore) {
						vecBest.push_back(TreePair(tempRefTree, tempTree));
						vecScore.push_back(theScore);
						vecMoveFromTo.push_back(StringTetrad(*s, b->first, a->first,
							computePartnerLGT(*s, b->first, a->first, testSets, refSets)));
					}
				}
			}
		}
	}

	const std::size_t n = vecBest.size();
	std::vector<StringTetrad> tempPath;
	for (std::size_t i = 0; i < n; ++i) {
		tempPath = thePath;
		tempPath.push_back(vecMoveFromTo[i]);
		recursiveRearrange(vecBest[i], vecScore[i], tempPath, vecPaths, shortestPathSoFar,
			startTime, timeLimit);
	}
}

PrintPathOutput printPath(const std::vector<StringTetrad>& thePath,
	std::map<std::string, int>& theSummary, std::ofstream& outFile, bool verboseMode)
{
	PrintPathOutput phnstCounts;
	for (std::vector<StringTetrad>::const_iterator i = thePath.begin(); i != thePath.end(); ++i) {
		if (verboseMode) std::cout << "Moving " << i->first << " from " << i->second <<
			" to within " << i->third << std::endl;
		outFile << "Moving " << i->first << " from " << i->second <<
			" to within " << i->third << std::endl;
		std::string theLGT = i->fourth + ", " + i->first;
		if (i->fourth.find("*\'") != std::string::npos) {
			++phnstCounts.c;
		} else if (i->fourth.find('*') != std::string::npos) {
			++phnstCounts.a;
		} else if (i->fourth.find('\'') != std::string::npos) {
			++phnstCounts.b;
		}
		if (verboseMode) std::cout << "  Putative LGT: { " << theLGT << " }" << std::endl;
		outFile << "  Putative LGT: { " << theLGT << " }" << std::endl;
		++theSummary[theLGT];
	}
	if (verboseMode) std::cout << std::endl;
	outFile << std::endl;
	return phnstCounts;
}

std::string computePartnerLGT(const std::string& s, const std::string& from,
 const std::string& to, TreeType& theTree, TreeType& theRefTree) {
	// result is <from> minus <s>, appended with a prime if <to> is included within <from>.
	// Also, result is appended with an asterisk if <from> is not a set in theRefTree.
	const std::set<std::string>& fromSet = theTree[from].branches_;
	if (fromSet.empty()) { // Just for debug:
		std::cout << "computePartnerLGT(), set is empty" << std::endl;
		std::exit(1);
	}
	std::set<std::string>::const_iterator iter = fromSet.begin();
	std::set<std::string>::const_iterator iterEnd = fromSet.end();
	std::string tempResult, tempString;
	std::set<std::string> tempSet;
	bool hasPlus = false;
	while (iter != iterEnd) {
		if (*iter != s) {
			if (!tempResult.empty()) {
				tempResult += " + ";
				hasPlus = true;
			}
			tempString = elaborateSetName(*iter, theTree);
			tempResult += tempString;
			tempSet.insert(tempString);
		}
		++iter;
	}
	std::string result;
	if (hasPlus) {
		result = "[" + tempResult + "]";
	} else {
		result = tempResult;
	}
	// Append an asterisk if tempResult is not in theRefTree:
	if (tempResult.find("+") != std::string::npos) {
		bool notInRefTree = true;
		for (TreeType::iterator ti = theRefTree.begin(); ti != theRefTree.end(); ++ti) {
			if (ti->second.branches_ == tempSet) {
				notInRefTree = false;
				break;
			}
		}
		if (notInRefTree) result += '*';
	} // If there is no '+', then it's a single leaf sequence, hence is in theRefTree.
	// Append a prime if <to> is included within <from>:
	if (fromSet.find(to) != fromSet.end()) result += '\'';
	return result;
}

std::string writeTree(const std::string& origTree, const std::string& s, TreeType& theTree)
{
	// launch using s = the name of the last set of theTree
	std::string answer(writeTreeRecursor(origTree, s, theTree));
	// trim the last bit off:
	std::string::size_type p = answer.find_last_of(')');
	answer[p + 1] = ';';
	for (std::string::size_type i = p + 2; i < answer.length(); ++i) {
		answer[i] = ' ';
	}
	return answer;
}

std::string writeTreeRecursor(const std::string& origTree, const std::string& s, TreeType& theTree)
{
	std::string result;
	if (isLeafSequence(s)) {
		result = s;
		std::string::size_type p = origTree.find(s);
		if (p == std::string::npos) {
			std::cout << "Error: " << s << " not found in original tree" << std::endl;
			std::exit(1);
		}
		p = origTree.find(':', p);
		if (p == std::string::npos) {
			std::cout << "Error: bad format in original tree" << std::endl;
			std::exit(1);
		}
		while (origTree[p] != ',' && origTree[p] != ')') {
			if (origTree[p] != '\n') result += origTree[p];
			++p;
		}
	} else {
		const std::set<std::string>& theSet = theTree[s].branches_;
		if (theSet.size() > 1) result += '(';
		std::set<std::string>::const_iterator iter = theSet.begin();
		std::set<std::string>::const_iterator iterEnd = theSet.end();
		while (iter != iterEnd) {
			result += writeTreeRecursor(origTree, *iter, theTree);
			if (++iter != iterEnd) result += ',';
		}
		if (theSet.size() > 1) {
			std::ostringstream ost;
			ost << ")'" << theTree[s].bootstrapValue_ << "':" << theTree[s].branchLength_;
			result += ost.str();
		}
	}
	return result;
}

std::string elaborateSetName(const std::string& s, TreeType& theTree)
{
	if (isLeafSequence(s)) return s;
	const std::set<std::string>& theSet = theTree[s].branches_;
	std::string result;
	if (theSet.size() > 1) result += '[';
	std::set<std::string>::const_iterator iter = theSet.begin();
	std::set<std::string>::const_iterator iterEnd = theSet.end();
	while (iter != iterEnd) {
		result += elaborateSetName(*iter, theTree);
		if (++iter != iterEnd) result += " + ";
	}
	if (theSet.size() > 1) result += ']';
	return result;
}

std::set<std::string> elaborateSet(const std::string& s, TreeType& theTree)
{
	std::set<std::string> result;
	if (isLeafSequence(s)) {
		result.insert(s);
	} else {
		const std::set<std::string>& theSet = theTree[s].branches_;
		std::set<std::string>::const_iterator iter = theSet.begin();
		std::set<std::string>::const_iterator iterEnd = theSet.end();
		while (iter != iterEnd) {
			std::set<std::string> temp(elaborateSet(*iter++, theTree));
			result.insert(temp.begin(), temp.end());
		}
	}
	return result;
}
	
inline std::string itoa(long number)
{
	std::ostringstream ost;
	ost << std::setw(3) << std::setfill('0') << number;
	return ost.str();
}

void printSets(const TreeType& theSets, std::ofstream& outFile)
{
	const bool toFile = outFile.is_open();
	TreeType::const_iterator mIter;
	std::set<std::string>::const_iterator sIter, sIterEnd;
	for (mIter = theSets.begin(); mIter != theSets.end(); ++mIter) {
		std::cout << mIter->first << " = |" << mIter->second.bootstrapValue_ << ':' <<
			mIter->second.branchLength_ << "| ";
		if (toFile) {
			outFile << mIter->first << " = |" << mIter->second.bootstrapValue_ << ':' <<
				mIter->second.branchLength_ << "| ";
		}
		sIter = mIter->second.branches_.begin();
		sIterEnd = mIter->second.branches_.end();
		while (sIter != sIterEnd) {
			std::cout << *sIter;
			if (toFile) outFile << *sIter;
			if (++sIter == sIterEnd) break;
			std::cout << " + ";
			if (toFile) outFile << " + ";
		}
		std::cout << std::endl;
		if (toFile) outFile << std::endl;
	}
	std::cout << std::endl;
	if (toFile) outFile << std::endl;
}

void printVerticalSets(const TreeType& theSets, std::map<std::string, int>& theSummary,
	std::ofstream& outFile, bool verboseMode)
{
	TreeType::const_iterator mIter;
	std::set<std::string>::const_iterator sIter, sIterEnd;
	std::vector<std::string> v;
	std::vector<std::string>::iterator viter;
	std::string::const_iterator s;
	for (mIter = theSets.begin(); mIter != theSets.end(); ++mIter) {
		sIter = mIter->second.branches_.begin();
		sIterEnd = mIter->second.branches_.end();
		for (; sIter != sIterEnd; ++sIter) {
			for (s = sIter->begin(); s != sIter->end(); ++s) {
				if (*s == '(') {
					for (viter = v.begin(); viter != v.end(); ++viter) {
						*viter += *s;
					}
					v.push_back(std::string());
				} else if (*s == ')') {
					if (v.empty()) { // for debug
						std::cout << "printVerticalSets(), vector is empty!" << std::endl;
						std::exit(1);
					}
					for (viter = v.begin(); viter + 1 != v.end(); ++viter) {
						*viter += *s;
					}
					if (verboseMode) std::cout << " (" << v.back() << ")" << std::endl;
					outFile << " (" << v.back() << ")" << std::endl;
					++theSummary[v.back()];
					v.pop_back();
				} else {
					for (viter = v.begin(); viter != v.end(); ++viter) {
						*viter += *s;
					}
				}
			}
		}
	}
}

void findVerticalSetsInString(const std::string& theString, std::set<std::string>& theSets)
{
	// This has much code in common with printVerticalSets(), so could refactor...
	std::vector<std::string> v;
	std::vector<std::string>::iterator viter;
	std::string::const_iterator s;
	for (s = theString.begin(); s != theString.end(); ++s) {
		if (*s == '(') {
			for (viter = v.begin(); viter != v.end(); ++viter) {
				*viter += *s;
			}
			v.push_back(std::string());
		} else if (*s == ')') {
			if (v.empty()) { // for debug
				std::cout << "printVerticalSetsFromString(), vector is empty!" << std::endl;
				std::exit(1);
			}
			for (viter = v.begin(); viter + 1 != v.end(); ++viter) {
				*viter += *s;
			}
			theSets.insert(v.back());
			v.pop_back();
		} else {
			for (viter = v.begin(); viter != v.end(); ++viter) {
				*viter += *s;
			}
		}
	}
}

long parseTree(std::string t, TreeType& m, char treeType)
{
	std::string::size_type p, pos1;
	std::vector<std::string> nameVec;
	std::string taxName, setname, stringy;
	std::string setbasename("set");
	std::string::size_type pos2 = t.find(')', 0);
	double branchLen = 0.0;
	int bootValue;
	long count = 0;
	do {
		pos1 = t.rfind('(', pos2);
		p = pos1;
		taxName = "";
		while (t[++p] != ':') if (t[p] != '\n') taxName += t[p];
		nameVec.assign(1, taxName);
		p = t.find(',', p);
		while (p != std::string::npos && p < pos2) {
			taxName = "";
			while (t[++p] != ':') if (t[p] != '\n') taxName += t[p];
			nameVec.push_back(taxName);
			p = t.find(',', p);
		}
		bootValue = 0;
		if (t[++pos2] == '\'') ++pos2;
		if (std::isdigit(t[pos2])) {
			stringy = t[pos2];
			while (t[++pos2] != ':' && t[pos2] != '\'') stringy += t[pos2];
			bootValue = int(std::atof(stringy.c_str()) + 0.5); // round up
		} else {
			bootValue = 100; // nothing was specified
		}
		if (t[pos2] == '\'') ++pos2;
		
		// Read the branch length, another indication of support:
		if (t[pos2] != ';') { // if not the last set
			while (t[pos2] != ':') ++pos2; // in case there's a label other than a bootstrap value
			p = pos2;
			stringy.clear();
			while (t[++p] != ',' && t[p] != ')') stringy += t[p];
			branchLen = std::atof(stringy.c_str());
			if (branchLen <= 0.0) bootValue = 0;
		}
		
		setname = setbasename + itoa(count++) + treeType;
		t.replace(pos1, pos2 - pos1, setname);
		std::set<std::string> theBranches(nameVec.begin(), nameVec.end());
		m[setname] = TreeBranches(theBranches, bootValue, branchLen);
		pos2 = t.find(')', 0);
	} while (pos2 != std::string::npos);
	return count;
}

void collapseWeakNodes(TreeType& theSets, const int minSupport)
{
	TreeType::iterator i, j;
	std::set<std::string>::iterator foundIter;
	bool doContinue, didFindIter;
	do {
		doContinue = false;
		for (i = theSets.begin(); i != theSets.end(); ++i) {
			if (i->second.bootstrapValue_ < minSupport) {
				// This node's branches need to be moved up to the parent,
				// then this node must be deleted. Find the parent:
				didFindIter = false;
				for (j = theSets.begin(); j != theSets.end(); ++j) {
					foundIter = j->second.branches_.find(i->first);
					if (foundIter != j->second.branches_.end()) {
						// j is the parent: remove foundIter, then add i's set to j's:
						j->second.branches_.erase(foundIter);
						j->second.branches_.insert(i->second.branches_.begin(), i->second.branches_.end());
						// And delete i:
						theSets.erase(i);
						didFindIter = true;
						break;
					}
				}
				doContinue |= didFindIter;
				break; // Start over from the beginning.
			}
		}
	} while (doContinue);
}

void pruneSequence(TreeType& theTree, const std::string& toRemove)
{
	// <toRemove> may be a compound sequence, so start by parsing it:
	std::vector<std::string> vecToRemove;
	std::string temp;
	for (std::string::const_iterator k = toRemove.begin(); k != toRemove.end(); ++k) {
		if (*k != '[' && *k != '(') {
			if (*k == ',' || *k == ']' || *k == ')') {
				if (!temp.empty()) {
					vecToRemove.push_back(temp);
					temp.clear();
				}
			} else {
				temp += *k;
			}
		}
	}
	if (!temp.empty()) vecToRemove.push_back(temp);
	std::set<std::string>::iterator foundIter;
	TreeType::iterator i;
	std::vector<std::string>::const_iterator vi;
	for (vi = vecToRemove.begin(); vi != vecToRemove.end(); ++vi) {
		for (i = theTree.begin(); i != theTree.end(); ++i) {
			i->second.branches_.erase(*vi);
			if (i->second.branches_.size() < 2) {
				// We can simplify the tree:
				bool foundSetName = false;
				for (TreeType::iterator j = theTree.begin(); j != theTree.end(); ++j) {
					foundIter = j->second.branches_.find(i->first);
					if (foundIter != j->second.branches_.end()) {
						j->second.branches_.erase(foundIter);
						j->second.branches_.insert(*(i->second.branches_.begin()));
						foundSetName = true;
						break;
					}
				}
				if (foundSetName) theTree.erase(i);
				break; // There is only one copy, so we're done with this name.
			}
		}
	}
}

void consolidateRefTopologies(TreeType& theRefTree, TreeType& theTestTree)
{
	if (theRefTree.size() == 1UL || theTestTree.size() == 1UL) return;
	TreeType::iterator i, j, k;
	std::set<std::string>::const_iterator sIter;
	std::set<std::string>::iterator foundIter;
	bool doContinue, foundSetName;
	do {
		doContinue = false;
		for (i = theRefTree.begin(); i != theRefTree.end(); ++i) {
			for (j = theTestTree.begin(); j != theTestTree.end(); ++j) {
				if (i->second.branches_ == j->second.branches_) {
					// Given that reference sets and test sets are named differently,
					// the only way for them to match is if they are identical terminal nodes.
					// Replace the entries containing i->first and j->first with the catenated
					// contents of i->second.branches_, then delete *i and *j:
					std::string theCat("(");
					sIter = i->second.branches_.begin();
					while (true) {
						theCat += *sIter;
						if (++sIter != i->second.branches_.end()) {
							theCat += ',';
						} else {
							theCat += ')';
							break;
						}
					}
					foundSetName = false;
					for (k = theRefTree.begin(); k != theRefTree.end(); ++k) {
						foundIter = k->second.branches_.find(i->first);
						if (foundIter != k->second.branches_.end()) {
							k->second.branches_.erase(foundIter);
							k->second.branches_.insert(theCat);
							foundSetName = true;
							break;
						}
					}
					if (foundSetName) theRefTree.erase(i);
					foundSetName = false;
					for (k = theTestTree.begin(); k != theTestTree.end(); ++k) {
						foundIter = k->second.branches_.find(j->first);
						if (foundIter != k->second.branches_.end()) {
							k->second.branches_.erase(foundIter);
							k->second.branches_.insert(theCat);
							foundSetName = true;
							break;
						}
					}
					if (foundSetName) theTestTree.erase(j);
					doContinue = true;
					break;
				}
			}
			if (doContinue) break; // Start over from the beginning.
		}
	} while (doContinue && theRefTree.size() > 1UL && theTestTree.size() > 1UL);
}

void consolidateSimilarities(TreeType& theRefTree, TreeType& theTestTree)
{
	if (theRefTree.size() == 1UL || theTestTree.size() == 1UL) return;
	std::vector<std::string> v;
	std::vector<std::string>::const_iterator vIter;
	std::set<std::string>::iterator foundIter;
	TreeType::iterator i, j, k;
	bool doContinue, foundSetName;
	do {
		doContinue = false;
		for (i = theRefTree.begin(); i != theRefTree.end(); ++i) {
			for (j = theTestTree.begin(); j != theTestTree.end(); ++j) {
				v.clear();
				std::set_intersection(i->second.branches_.begin(), i->second.branches_.end(),
					j->second.branches_.begin(), j->second.branches_.end(), std::back_inserter(v));
				if (v.size() > 1UL) {
					// The sets *i and *j can be simplified:
					std::string theCat("[");
					vIter = v.begin();
					while (true) {
						theCat += *vIter;
						// Remove the item from *i and from *j:
						i->second.branches_.erase(*vIter);
						j->second.branches_.erase(*vIter);
						if (++vIter != v.end()) {
							theCat += ',';
						} else {
							theCat += ']';
							break;
						}
					}
					// Add the new name to *i and to *j, or to their parents:
					if (i->second.branches_.empty()) {
						foundSetName = false;
						for (k = theRefTree.begin(); k != theRefTree.end(); ++k) {
							foundIter = k->second.branches_.find(i->first);
							if (foundIter != k->second.branches_.end()) {
								k->second.branches_.erase(foundIter);
								k->second.branches_.insert(theCat);
								foundSetName = true;
								break;
							}
						}
						if (foundSetName) theRefTree.erase(i);
					} else {
						i->second.branches_.insert(theCat);
					}
					if (j->second.branches_.empty()) {
						foundSetName = false;
						for (k = theTestTree.begin(); k != theTestTree.end(); ++k) {
							foundIter = k->second.branches_.find(j->first);
							if (foundIter != k->second.branches_.end()) {
								k->second.branches_.erase(foundIter);
								k->second.branches_.insert(theCat);
								foundSetName = true;
								break;
							}
						}
						if (foundSetName) theTestTree.erase(j);
					} else {
						j->second.branches_.insert(theCat);
					}
					doContinue = true;
					break;
				}
			}
			if (doContinue) break; // Start over from the beginning.
		}
	} while (doContinue && (theRefTree.size() > 1UL || theTestTree.size() > 1UL));
}

void moveItem(TreeType& theTree, const std::string& toKey, const std::string& fromKey,
	const std::string& theItem)
{
	theTree[fromKey].branches_.erase(theItem);
	theTree[toKey].branches_.insert(theItem);
	if (theTree[fromKey].branches_.size() < 2) {
		// That node needs to be collapsed:
		const std::string leftover = *(theTree[fromKey].branches_.begin());
		std::set<std::string>::iterator foundIter;
		bool foundSetName = false;
		for (TreeType::iterator k = theTree.begin(); k != theTree.end(); ++k) {
			foundIter = k->second.branches_.find(fromKey);
			if (foundIter != k->second.branches_.end()) {
				k->second.branches_.erase(foundIter);
				k->second.branches_.insert(leftover);
				foundSetName = true;
				break;
			}
		}
		if (foundSetName) theTree.erase(fromKey);
	}
}

std::size_t numLeavesInTree(const TreeType& theTree)
{
	std::size_t num = 0UL;
	std::set<std::string>::const_iterator j;
	for (TreeType::const_iterator i = theTree.begin(); i != theTree.end(); ++i) {
		for (j = i->second.branches_.begin(); j != i->second.branches_.end(); ++j) {
			if (isLeafSequence(*j)) ++num;
		}
	}
	return num;
}

bool isLeafSequence(const std::string& theName)
{
	return (theName.length() < 4 || theName[0] != 's' || theName[1] != 'e' || theName[2] != 't' ||
		!std::isdigit(theName[3]));
}

void trimMultipleUnderscores(std::string& theString)
{
	std::string::size_type p = theString.find("__");
	while (p != std::string::npos) {
		theString.replace(p, 2, "");
		p = theString.find("__");
	}
	p = theString.find("_:");
	while (p != std::string::npos) {
		theString.replace(p, 2, ":");
		p = theString.find("_:");
	}
}

void findLeafSets(const TreeType& theTree, std::set<std::string>& theLeafSet)
{
	std::set<std::string>::const_iterator j;
	for (TreeType::const_iterator i = theTree.begin(); i != theTree.end(); ++i) {
		for (j = i->second.branches_.begin(); j != i->second.branches_.end(); ++j) {
			if (isLeafSequence(*j)) theLeafSet.insert(*j);
		}
	}
}

void printOneSet(const TreeType& theSets, const std::string& setName, std::ofstream& outFile)
{
	TreeType::const_iterator mIter = theSets.find(setName);
	if (mIter == theSets.end()) {
		outFile << "Error: set not found.\n";
	} else {
		std::set<std::string>::const_iterator sIter, sIterEnd;
		sIter = mIter->second.branches_.begin();
		sIterEnd = mIter->second.branches_.end();
		while (sIter != sIterEnd) {
			outFile << *sIter;
			if (++sIter == sIterEnd) break;
			outFile << " + ";
		}
		outFile << std::endl;
	}
}

int computeNumberOfNodesFromRoot(TreeType& theTree, const std::string& name, int count)
{
	TreeType::iterator i;
	for (i = theTree.begin(); i != theTree.end(); ++i) {
		if (i->second.branches_.find(name) != i->second.branches_.end()) {
			if (i->first == theTree.rbegin()->first) return count;
			else return computeNumberOfNodesFromRoot(theTree, i->first, count + 1);
		}
	}
	return -1;
}

// classCmdLineArgParser.cpp
// Version 2005.10.18
// © 2001-2005, Robert L. Charlebois

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "classCmdLineArgParser.h"
#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <iostream>
using namespace std;

CmdLineArgParser::CmdLineArgParser(int argc, const char* argv[])
{
	theArgs_.reserve(argc);
	for (int i = 0; i < argc; ++i) {
		theArgs_.push_back(std::string(argv[i]));
	}
}

bool CmdLineArgParser::parse(const std::string& key, std::string* theString, bool needArg)
{
	std::vector<std::string>::iterator vi = std::find(theArgs_.begin(), theArgs_.end(), key);
	if (vi == theArgs_.end()) {
		if (needArg) throw 0;
		else return 0; // key not found.
	}
	*vi = ""; //This key has been successfully recognized. Remove
	if (++vi == theArgs_.end()) throw -1; // No argument specified: an error.
	*theString = std::string(*vi);
	*vi = ""; //This argument has been successfully recognized. Remove
	
	if ((*theString)[0] == '-' &&
		(theString->length() < 2 || (*theString)[1] == '-' || std::isalpha((*theString)[1]))) throw -1;
		// The argument is missing: an error.
		//cout << "found: " << key << "=" << *theString << endl; //debug code
	return true; // Success.
}

bool CmdLineArgParser::parse(const std::string& key, int* theInt, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theInt = std::atoi(stringRep.c_str());
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, unsigned int* theUInt, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theUInt = (unsigned int)std::atol(stringRep.c_str());
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, long long* theLongLong, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theLongLong = std::atol(stringRep.c_str()); // •• atol is limited ••
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, unsigned long long* theULongLong, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theULongLong = (unsigned long long)std::atol(stringRep.c_str()); // •• atol is limited ••
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, long* theLong, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theLong = std::atol(stringRep.c_str());
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, unsigned long* theULong, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theULong = (unsigned long)std::atol(stringRep.c_str());
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, short* theShort, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theShort = (short)std::atoi(stringRep.c_str());
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, unsigned short* theUShort, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theUShort = (unsigned short)std::atoi(stringRep.c_str());
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, double* theDouble, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theDouble = std::atof(stringRep.c_str());
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, float* theFloat, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theFloat = (float)std::atof(stringRep.c_str());
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, bool* theBool, bool needArg)
{
	std::string stringRep;
	bool resultCode = parse(key, &stringRep, needArg);
	if (resultCode == true) *theBool = (stringRep[0] == 'T');
	return resultCode;
}

bool CmdLineArgParser::parse(const std::string& key, bool needArg)
{
	std::vector<std::string>::iterator found = std::find(theArgs_.begin(), theArgs_.end(), key);
	bool result = (found != theArgs_.end());
	if (needArg && !result) throw 0; // Why anyone would specify needArg here is a mystery.
	if (result)
	{	
		*found = ""; //This key has been successfully recognized. Remove
		//cout << "found: " << key << endl;  //debug code
	}
	return result;
}

// abezgino 2011.01.21
//prints command-line options that have not been yet processed
void CmdLineArgParser::print()
{
	for (std::vector<std::string>::iterator it=theArgs_.begin()+1; it!=theArgs_.end(); it++)
		if (!(*it).empty()) std::cout << *it << " ";
	std::cout << std::endl;
	return;
}

//returns whether there no options that have not been recognized yet
bool CmdLineArgParser::empty()
{
	for (std::vector<std::string>::iterator it=theArgs_.begin()+1; it!=theArgs_.end(); it++)
		if (! (*it).empty()) return false;
	return true;
}

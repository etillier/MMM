// classCmdLineArgParser.h
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


#ifndef CLASS_CMDLINEARGPARSER_H
#define CLASS_CMDLINEARGPARSER_H

#include <string>
#include <vector>

class CmdLineArgParser {
private:
	std::vector<std::string> theArgs_;
public:
	CmdLineArgParser(int argc, const char* argv[]);
	~CmdLineArgParser() { }

	const std::string& programName() const { return theArgs_.front(); }

	bool parse(const char* key, std::string* theString, bool needArg = false)
		{ return parse(std::string(key), theString, needArg); }	
	bool parse(const char* key, int* theInt, bool needArg = false)
		{ return parse(std::string(key), theInt, needArg); }	
	bool parse(const char* key, unsigned int* theUInt, bool needArg = false)
		{ return parse(std::string(key), theUInt, needArg); }	
	bool parse(const char* key, long long* theLongLong, bool needArg = false)
		{ return parse(std::string(key), theLongLong, needArg); }	
	bool parse(const char* key, unsigned long long* theULongLong, bool needArg = false)
		{ return parse(std::string(key), theULongLong, needArg); }	
	bool parse(const char* key, long* theLong, bool needArg = false)
		{ return parse(std::string(key), theLong, needArg); }	
	bool parse(const char* key, unsigned long* theULong, bool needArg = false)
		{ return parse(std::string(key), theULong, needArg); }	
	bool parse(const char* key, short* theShort, bool needArg = false)
		{ return parse(std::string(key), theShort, needArg); }	
	bool parse(const char* key, unsigned short* theUShort, bool needArg = false)
		{ return parse(std::string(key), theUShort, needArg); }	
	bool parse(const char* key, double* theDouble, bool needArg = false)
		{ return parse(std::string(key), theDouble, needArg); }	
	bool parse(const char* key, float* theFloat, bool needArg = false)
		{ return parse(std::string(key), theFloat, needArg); }	
	bool parse(const char* key, bool* theBool, bool needArg = false)
		{ return parse(std::string(key), theBool, needArg); }	
	bool parse(const char* key, bool needArg = false)
		{ return parse(std::string(key), needArg); }	

	bool parse(const std::string& key, std::string* theString, bool needArg = false);	
	bool parse(const std::string& key, int* theInt, bool needArg = false);
	bool parse(const std::string& key, unsigned int* theUInt, bool needArg = false);
	bool parse(const std::string& key, long long* theLongLong, bool needArg = false);
	bool parse(const std::string& key, unsigned long long* theULongLong, bool needArg = false);
	bool parse(const std::string& key, long* theLong, bool needArg = false);
	bool parse(const std::string& key, unsigned long* theULong, bool needArg = false);
	bool parse(const std::string& key, short* theShort, bool needArg = false);
	bool parse(const std::string& key, unsigned short* theUShort, bool needArg = false);
	bool parse(const std::string& key, double* theDouble, bool needArg = false);
	bool parse(const std::string& key, float* theFloat, bool needArg = false);
	bool parse(const std::string& key, bool* theBool, bool needArg = false); // Assumes T or F.
	bool parse(const std::string& key, bool needArg = false);
	
	void print();
	bool empty();
};

#endif

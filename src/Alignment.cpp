#include "Alignment.hpp"

//Read an alignment in fasta format. If a non-negative gapThreshold is supplied then
//sequences with the proportion of gaps exceeding this value will not be included
void Alignment::read_alignment(string const& inFilename, double gapThreshold)
{	
	ifstream inFile(inFilename.c_str(), ios::in);
	if(!inFile.is_open())
	{	//cerr << "Error: cannot read file: " << inFilename << endl;
		//abort();
		throw runtime_error("cannot read file: " + inFilename);
	}
	
	char buffer[BUFFER_SIZE] = {0};
	bool sequenceInProgress = false;
	string headerCurrent = "";
	string sequenceCurrent = "";
	int nOfGapsCurrent = 0;
	
	while(!inFile.eof())
	{	inFile.getline(buffer,sizeof(buffer));
		if(inFile.fail() && !inFile.eof()) //can happen if the buffer is too small to fit the line
		{	//cerr << "Error: could not read line"<< endl;
			//abort();
			throw runtime_error("could not read line");
		}
		
		//remove trailing carriage return character (applicable to files in DOS format)
		if (buffer[strlen(buffer)-1] == '\r')
		{	buffer[strlen(buffer)-1] == '\0';
		}
		
		if (buffer[0] == '>') //started reading a new sequence
		{	//need to store the previous sequence (unless this is the first one)
			if (sequenceInProgress && (gapThreshold < 0.0 || (double) nOfGapsCurrent / sequenceCurrent.length() < gapThreshold))
			{	add_sequence(headerCurrent, sequenceCurrent);
			}
			
			//initialize the next sequence
			sequenceInProgress = true;
			headerCurrent = buffer+1;//skip the initial '>' character
			sequenceCurrent = "";
			nOfGapsCurrent = 0;
		} else
		{	if (sequenceInProgress)
			{	nOfGapsCurrent += prepare_sequence(buffer);
				sequenceCurrent+=buffer;
			}
		}
	}
	inFile.close();
	//don't forget to add the last sequence 
	if (sequenceInProgress && (gapThreshold < 0 || (double) nOfGapsCurrent / sequenceCurrent.length() < gapThreshold))
	{	add_sequence(headerCurrent, sequenceCurrent);
	}
	
	//a sequence alignment must have at least 2 sequences
	if(nOfSequences_ > 1)
		valid_ = true;
}
//

void Alignment::add_sequence(string header, string sequence)
{	if (nOfColumns_ < 0)
	{	nOfColumns_ = sequence.length();
	} else if (nOfColumns_ != sequence.length())
	{	stringstream err;
		err << "inconsistent number of columns: " << nOfColumns_ << " != " << sequence.length() << " for sequence: " << header;
		//abort();
		throw runtime_error(err.str());
	}
		
	headers_.push_back(header);
	sequences_.push_back(sequence);
	nOfSequences_++;
}
//

//ensures that only 20 AA letters and '-' are used in this sequence
//returns the number of gaps
int Alignment::prepare_sequence(char *buffer)
{	int nOfGaps = 0;
	
	for(int i=0; i < strlen(buffer); i++)
	{	if(isupper(buffer[i]))
		{	switch (buffer[i])
			{	case 'U':
					cerr << "Warning: 'U' (selenocysteine) will be converted to 'C' cysteine" << endl;
					buffer[i] = 'C';
					break;
				case 'O':
					cerr << "Warning: 'O' (pyrrolysine) will be converted to 'K' lysine" << endl;
					buffer[i] = 'K';
					break;
				case 'B':
				case 'J':
				case 'Z':
					cerr << "Warning: '" << buffer[i] << "' (non-standard residue) will be converted to 'X' unknown" << endl;
					buffer[i] = 'X';
				case 'X': //unknown AA
					nOfGaps++; //treated as a gap
					break;
				
				//{	cerr << "Error: unexpected FASTA character: " << buffer[i] << endl;
				//	exit(1);
				//}
			}
		} else
		{	switch (buffer[i])
			{	case '.': //abezgino: why is dot consireded a gap? OK. pfam format
				case '_':
					buffer[i] == '-'; //standartize gaps
				case '-':
					nOfGaps++;
					break;
				default:
				{	stringstream err;
					err << "unexpected FASTA character: [" << buffer[i] << "]";
					throw runtime_error(err.str());
				}
			}
		}
	}
	return nOfGaps;
}
//

//appends sequences from another MSA. Only the sequences with matching taxon IDs are processed
Alignment Alignment::append_alignment(Alignment const& other)
{	validate();
	typedef map <string, int> TaxonMap;
	
	//cerr << "Appending alignment\n";
	TaxonMap taxon2index;
		
	//store the list of taxons found in the first alignment 
	for(int i=0; i<get_nOfSequences(); i++)
	{	string taxon = get_taxon(i);
		//todo: abort if mapping already exists, i.e. we can't handle paralogs
		taxon2index[taxon]=i;
	}
	
	Alignment alignment(get_nOfColumns());//label the alignment as joined. This will help to recover the original aminoacid positions
	
	//go through the sequences of the second alignment, and if the taxon matches a taxon of a sequence in the first alignment then join the two sequences 
	for(int j=0; j<other.get_nOfSequences(); j++)
	{	string taxon = other.get_taxon(j);
		TaxonMap::iterator it = taxon2index.find(taxon);
		if(it != taxon2index.end())
		{	int i = it->second; //the index of the corresponding sequence in the first alignment
			
			string headerJoined = get_header(i) + "+" + other.get_header(j); //concatenate the headers
			
			string sequenceJoined = get_sequence(i) + other.get_sequence(j);
			alignment.add_sequence(headerJoined, sequenceJoined);
		}
	}
	
	
	//cerr << flush;
	
	return alignment;
}
//

Alignment Alignment::get_slices(vector<pair<int, int > > const& ranges) const
{	validate();
	Alignment alignmentSlice;
	for(int i=0; i<get_nOfSequences(); i++)
	{	string sequenceSlice = "";
		int min=1;
		for(int j=0; j<ranges.size(); j++)
		{	int const& start = ranges.at(j).first;
			int const& end   = ranges.at(j).second;
			
			if(start < min || start > end || end > get_nOfColumns())
			{	stringstream err;
				err << "invalid slice range: " << start << "-" << end << "(";
				for(int k=0; k<ranges.size(); k++)
				{	if (k>0)
						 err << ",";
					err << ranges.at(k).first << "-" << ranges.at(k).second;
				}
				err << ") for alignment 1-" << get_nOfColumns();
				throw runtime_error(err.str());
			}
			
			min=end+1;
			int raw_start = start-1;
			int raw_length = end-start+1;
			sequenceSlice += get_sequence(i).substr(raw_start, raw_length);
		}
		
		size_t nonGapFound = sequenceSlice.find_first_not_of("-");
		if (nonGapFound == string::npos)
		{	cout << "Warning: only gaps found within range (";
			for(int k=0; k<ranges.size(); k++)
			{	if (k>0)
					cout << ",";
				cout << ranges.at(k).first << "-" << ranges.at(k).second;
			}
			cout << ") Skipping sequence:\n" << get_header(i) << endl;
		} else
		{	alignmentSlice.add_sequence(get_header(i), sequenceSlice);
		}
	}
		
	if (alignmentSlice.nOfSequences_ > 0)
		alignmentSlice.valid_=true;
	if (alignmentSlice.nOfSequences_ == 1)
		cout << "Warning: this alignment slice only has one sequence" << endl;
	return alignmentSlice;
}
//

void Alignment::print_to_file(string const& outFilepath) const
{	validate();
	ofstream outFile(outFilepath.c_str());
	if (!outFile.is_open())
	{	throw runtime_error("cannot create file: " + outFilepath);
	}
		
	for (int i=0; i<nOfSequences_; i++)
	{	outFile << '>' << headers_.at(i) << '\n' << sequences_.at(i) << '\n';
	}
		
	outFile.close();
}
//
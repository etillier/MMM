#include "ProteinGroup.hpp"

void ProteinGroup::validate(void) const
{	if (! valid_)
	{	stringstream err;
		//err << "attempting to use an invalid " << typeid(this).name();
		err << "attempting to use an invalid ProteinGroup";
		throw runtime_error(err.str());
	}
}

//extracts taxon from a header. Everything preceeding the '|' character is considered to be the taxon
string ProteinGroup::get_taxon(int index) const
{	validate();
	string const& header = headers_.at(index);
	size_t pos = header.find('|');
	if (pos == string::npos)
	{	//cerr << "Error: could not find taxon in the header: " << header << endl;
		//exit (1);
		stringstream err;
		err << "could not extract a taxon from the header: " << header;
		throw runtime_error(err.str());
	}
	string taxon = header.substr(0,pos); //extract the taxon from the header
	return taxon;
}
//


//bool ProteinGroup::find_in_cache(list <ProteinGroup>& cache, string const& name)
//{	bool found = false;
//	//unsigned int ii = 0;
//	for(list <ProteinGroup>::iterator it = cache.begin(); it != cache.end() && !found; it++)
//	{	//ii++;
//		if (it->get_name() == name)
//		{	found = true;
//			//move this alignment to the front (as the most recently used)
//			cache.splice(cache.begin(), cache, it);
//		}
//	}
//	return found;
//}
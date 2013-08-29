#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <stack>
using namespace std;

#if defined (_MSC_VER)
#	define FORCEINLINE		__forceinline
#elif defined (__GNUG__)
#	define FORCEINLINE		__inline
#else
#	define FORCEINLINE		inline
#endif

#include "bonus.h"
#include "profiler.h"
#include "hw_iface.h"
#include "mmm_algorithm.h"

#ifdef EDGE_STATS
static unsigned long long int s_total_edges = 0;
static unsigned long long int s_total_verts = 0;
static unsigned long long int s_total_sq_verts = 0;
#endif

// Lower limit of what '0' is for matrix entries 
const double TOO_SMALL = -16.118095650958319788125940182791; //ln(1.0e-7)

//
// Funcs
//

CHardAlgorithm::CHardAlgorithm
(
	const double allow,
	const vector<double>& distances1, const vector<double>& distances2,
	const vector<int>& taxLabels1, const vector<int>& taxLabels2, const bool useTaxInfo,
	const int& reqTaxName, const bool reqTax, 
	const int strictSize, const int maxTrees, const int minSize,
	bool useHW
)
: CMMMAlgorithm(allow, distances1, distances2, taxLabels1, taxLabels2, useTaxInfo,
reqTaxName, reqTax, strictSize, minSize, maxTrees), _useHW(useHW)
{
	// Take the logs of all the distances
	int ptr = 0;
	for (int i = 0; i < _num1; i++)
	{
		for (int j = 0; j < _num1; j++)
		{
			// Diagonal entries get treated as zero/bad distances
			if (i == j || _distances1[ptr] <= 0.0)
			{
				_distances1[ptr] = TOO_SMALL;
			}
			else
			{
				double log_dist = log(_distances1[ptr]);
				_distances1[ptr] = (log_dist < TOO_SMALL)? TOO_SMALL : log_dist;
			}

			ptr++;
		}
	}

	// Do the same for distance matrix 2/B
	ptr = 0;
	for (int i = 0; i < _num2; i++)
	{
		for (int j = 0; j < _num2; j++)
		{
			if (i == j || _distances2[ptr] <= 0.0)
			{
				_distances2[ptr] = TOO_SMALL;
			}
			else
			{
				double log_dist = log(_distances2[ptr]);
				_distances2[ptr] = (log_dist < TOO_SMALL)? TOO_SMALL : log_dist;
			}

			ptr++;
		}
	}

	// Convert the -a tolerance factor into the logged distance measure 'D'
	_D = log( (1.0 + _allow) / (1.0 - _allow) );

	_noTrees = (_maxTrees == 0) || (_maxTrees == 1);

	// Allocate data structures
	_numberLine.reserve(1024);
	_setVerts.reserve(1024);
	_setProbs.reserve(1024);
	_setInputBuffer.reserve(131072);

	_probVerts.reserve(1024);
	_probMatrix.reserve(65536);
	_probOrder.reserve(1024);
	_probDegree.reserve(1024);
	_probColorClass.reserve(65536);
	_probColorClassSize.reserve(1024);

	_setOutputBuffer.reserve(16384);
	_mcVStack.resize(16384);
	_mcCliqueStack.resize(32);
	_mcPtrStack.resize(128);
	_mcMSPV.reserve(1024);
}

CHardAlgorithm::~CHardAlgorithm()
{
}


bool CHardAlgorithm::initHardware()
{
	return hw_init();
}


void CHardAlgorithm::closeHardware()
{
	hw_close();
}


bool CHardAlgorithm::stateRestore(string& fileName)
{
	cout << "State restore not supported." << endl;
	return false;
}


void CHardAlgorithm::stateSave(string& fileName)
{
	cout << "State save not supported." << endl;
}


bool CHardAlgorithm::maxTreesReached()
{
	return (int)_resultVec.size() == _maxTrees;
}


int CHardAlgorithm::lastReachedIndex()
{
	return _num1;
}


bool CHardAlgorithm::validProteinPair(int a, int b)
{
	return (!_useTaxInfo || (_taxLabels1[a] == _taxLabels2[b]));
}

bool CHardAlgorithm::satisfiesReqTax(int a)
{
	return (!_reqTax || (_taxLabels1[a] == _reqTaxName));
}

void CHardAlgorithm::printBonusStats(ostream& str, bool brief)
{
	#ifdef CACHE_SIM
	if (brief) str << hw_sim_cache_get_hit_ratio() << '\t';
	else str << "Cache hit/miss ratio: " << hw_sim_cache_get_hit_ratio() << endl;
	#endif

	#ifdef HW_COUNTERS
	if (_useHW)
	{
		unsigned long long int total, wu, usage;
		hw_dbg_get_counters(&total, &wu, &usage);
		
		double total_time = (double)(total) / HW_TICKS_PER_S;
		double working_time = (double)(wu) / HW_TICKS_PER_S;
		double usage_time = (double)(usage) / HW_TICKS_PER_S;

		double util_frac = usage_time / (working_time + 1e-7);
	
		if (brief) str << total_time << '\t' << working_time << '\t' << util_frac << '\t';
		else 
		{
			str << "HW total time: " << total_time << endl;
			str << "HW solving time: " << working_time << endl;
			str << "HW utilization: " << util_frac << endl;
		}
	}
	#endif

	#ifdef EDGE_STATS
	double density = (double)s_total_edges / (double)s_total_sq_verts;
	if (brief) str << _num1 << '\t' << _num2 << '\t' << s_total_verts << '\t' << '\t' << density << '\t';
	else 
	{
		str << "N: " << _num1 << " M: " << _num2 << endl;
		str << "Total vertices: " << s_total_verts << endl;
		str << "Average density: " << density << endl;
	}
	#endif

	#ifdef PROFILE
	if (brief)
	{
		float clique_time = _useHW ? PROFILER_OUTPUT_ZONE_TIME(maxCliquesHard) :
			PROFILER_OUTPUT_ZONE_TIME(maxCliquesEmu);
		float setup_time = PROFILER_OUTPUT_ZONE_TIME(buildProblemSet) +
			PROFILER_OUTPUT_ZONE_TIME(buildNumberLine) +
			PROFILER_OUTPUT_ZONE_TIME(copyCliques);

		str << setup_time << '\t' << clique_time << '\t';
	}
	else str << PROFILER_OUTPUT_FLAT_STRING() << endl;
	#endif
}


FORCEINLINE bool CHardAlgorithm::forwardCompatible(int min_a, int i_a, int i_b, int j_a, int j_b, const double& r_min, const double& r_max)
{
	double dist_a = getDist1(i_a,j_a);
	double dist_b = getDist2(i_b,j_b);
	double ratio = dist_a - dist_b;

	if (dist_a == TOO_SMALL || dist_b == TOO_SMALL ||
		ratio < r_min || ratio >= r_max)
	{
		return false;
	}
	else if (ratio == r_min && (i_a < min_a || j_a < min_a))
	{
		// Fix for duplicate ratios=minRatio
		return false;
	}

	return true;
}


void CHardAlgorithm::buildNumberLine(int i_a, int i_b)
{
	PROFILE_FUNC();

	_numberLine.clear();

	CNLEntry e;
	for (int j_a = 0; j_a < _num1; j_a++)
	{
		double dist_a = getDist1(i_a,j_a);
		if (dist_a == TOO_SMALL)
			continue;

		for (int j_b = 0; j_b < _num2; j_b++)
		{
			double dist_b = getDist2(i_b,j_b);
			if (dist_b == TOO_SMALL)
				continue;

			if (!validProteinPair(j_a, j_b))
				continue;

			e.dest_a = j_a;
			e.dest_b = j_b;
			e.ratio = dist_a - dist_b;

			_numberLine.push_back(e);
		}
	}

	// Sort the numberline by ascending ratio
	std::sort(_numberLine.begin(), _numberLine.end(), CSortEntries());
}


bool CHardAlgorithm::buildMatrix(int min_a, const double& r_min, const double& r_max)
{
	PROFILE_FUNC();

	bool emptyMatrix = true;

	int nverts = (int)_probVerts.size();
	_probMatrix.resize(nverts*nverts);
	_probDegree.resize(nverts);

	// Initialize degrees
	for (int i = 0; i < nverts; i++)
	{
		_probDegree[i] = 0;
	}

	#ifdef EDGE_STATS
	s_total_verts += nverts;
	s_total_sq_verts += nverts*(nverts+1) / 2;
	#endif

	// Build matrix and calculate degrees
	for (int i_idx = 0; i_idx < nverts; i_idx++)
	{
		const CNLEntry& i = _numberLine[_probVerts[i_idx]];
		int i_a = i.dest_a;
		int i_b = i.dest_b;

		// Set diagonal to 0
		_probMatrix[i_idx*(nverts+1)] = 0;

		for (int j_idx = i_idx+1; j_idx < nverts; j_idx++)
		{
			const CNLEntry& j = _numberLine[_probVerts[j_idx]];
			int j_a = j.dest_a;
			int j_b = j.dest_b;

			bool conn = forwardCompatible(min_a, i_a, i_b, j_a, j_b, r_min, r_max);
			if (conn)
			{
				emptyMatrix = false;
				_probMatrix[i_idx*nverts + j_idx] = 1;
				_probMatrix[j_idx*nverts + i_idx] = 1;
				_probDegree[i_idx]++;
				_probDegree[j_idx]++;

				#ifdef EDGE_STATS
				s_total_edges++;
				#endif
			}
			else
			{
				_probMatrix[i_idx*nverts + j_idx] = 0;
				_probMatrix[j_idx*nverts + i_idx] = 0;
			}
		}
	}

	// Early exit: matrix is all zeros, which means the
	// largest cliques/matches will be of size 1/3 respectively.
	if (emptyMatrix && _maxScore > 3)
		return false;

	return true;
}


bool CHardAlgorithm::sortVertices()
{
	PROFILE_FUNC();

	int nverts = (int)_probVerts.size();

	// Initialize order and color arrays, and determine max degree
	_probColorClassSize.resize(nverts);
	_probOrder.resize(nverts);
	int sufDegreeVertices = 0;
	int num_colors = 0;

	for (int i = 0; i < nverts; i++)
	{
		_probOrder[i] = i;
		_probColorClassSize[i] = 0;

		bool cond = _noTrees?
			_probDegree[i]+1 + 2 > _maxScore :
			_probDegree[i]+1 + 2 >= _maxScore;

		if (cond)
			sufDegreeVertices++;
	}

	// Early exit: not enough sufficient-degree vertices to match or beat max clique size
	bool cond = _noTrees?
		sufDegreeVertices + 2 <= _maxScore :
		sufDegreeVertices + 2 < _maxScore;

	if (cond)
		return false;

	// Sort vertices by degree
	std::sort(_probOrder.begin(), _probOrder.end(), CSortByDegree(_probDegree));

	// Do greedy vertex coloring
	_probColorClass.resize(nverts*nverts);

	for (int i = 0; i < nverts; i++)
	{
		int ii = _probOrder[i];

		bool found = false;
		for (int c = 0; c < num_colors; c++)
		{
			found = true;

			// See if ii connects to anything in this color class
			int cc_size = _probColorClassSize[c];
			for (int v = 0; v < cc_size; v++)
			{
				int jj = _probColorClass[c*nverts + v];
				if (_probMatrix[ii*nverts + jj])
				{
					found = false;
					break;
				}
			}

			// If not, assign it that color class and we're done
			if (found)
			{
				_probColorClass[c*nverts + cc_size] = ii;
				_probColorClassSize[c] = cc_size+1;
				break;
			}
		}

		// No color class will take us. Create a new one.
		if (!found)
		{
			_probColorClassSize[num_colors] = 1;
			_probColorClass[num_colors*nverts + 0] = ii;
			num_colors++;
		}
	}

	// Early exit: number of color classes is smaller than max clique size
	cond = _noTrees?
		num_colors + 2 <= _maxScore :
		num_colors + 2 < _maxScore;

	if (cond)
		return false;

	// Sort by decreasing vertex color class, and by decreasing degree within each color class

	int order_ptr = 0;
	for (int i = num_colors-1; i >= 0; i--)
	{
		int col_ptr = i*nverts;
		for (int j = 0; j < _probColorClassSize[i]; j++)
		{
			_probOrder[order_ptr++] = _probColorClass[col_ptr++];
		}
	}

	_probNumColorClasses = num_colors;

	return true;
}


void CHardAlgorithm::appendDataBuffer(int* header_idx)
{
	PROFILE_FUNC();

	int nverts = _probVerts.size();
	int nbits = 0;
	unsigned int val = 0;

	// First value (header) is the number of vertices in this problem.
	// It also contains a bit that says whether this is the last problem of the set.
	// That's unknown right now, so the header will be modified at the very end.
	// We return a pointer to this header to let that happen.
	*header_idx = (int)_setInputBuffer.size();
	_setInputBuffer.push_back(nverts);

	// Insert the matrix next, literally bit by bit
	for (int i = 0; i < nverts; i++)
	{
		int ii = _probOrder[i];
		for (int j = 0; j < nverts; j++)
		{
			int jj = _probOrder[j];

			val >>= 1;
			val |= _probMatrix[ii*nverts + jj] * 0x80000000;
			nbits++;

			if (nbits == 32)
			{
				nbits = 0;
				_setInputBuffer.push_back(val);
				val = 0;
			}
		}
	}

	if (nbits != 0)
	{
		val >>= (32 - nbits);
		_setInputBuffer.push_back(val);
	}
}


void CHardAlgorithm::buildProblemSet(int i_a, int i_b)
{
	PROFILE_FUNC();

	_setVerts.clear();
	_setProbs.clear();

	int start_idx = 0;
	CProblemDesc prob_desc;

	// Initialize data buffer with current max clique size
	_setInputBuffer.clear();

	unsigned int global_header = std::max(0, _maxScore - 2);
	if (_noTrees) global_header |= 0x80000000;
	_setInputBuffer.push_back(global_header);

	int last_header_idx = -1;

	int numberLine_size = (int)_numberLine.size();
	for (int j = 0; j < numberLine_size; j++)
	{
		int j_a = _numberLine[j].dest_a;
		int j_b = _numberLine[j].dest_b;

		// Only consider lexicographically-later vertices than i, so that
		// the i and j loops together visit all edges only once
		if (j_a < i_a)
			continue;

		// Set min/max ratio
		double r_min = _numberLine[j].ratio;
		double r_max = r_min + _D;

		// Gather forward-compatible vertices from the numberline
		_probVerts.clear();
		for (int k = j+1; k < numberLine_size; k++)
		{
			int k_a = _numberLine[k].dest_a;
			int k_b = _numberLine[k].dest_b;

			if (_numberLine[k].ratio >= r_max)
				break;

			if (!forwardCompatible(i_a, j_a, j_b, k_a, k_b, r_min, r_max))
				continue;

			_probVerts.push_back(k);
		}

		// Early exit: not enough vertices in problem to match/beat maximum
		int probVerts_size = (int)_probVerts.size();
		bool not_enough_verts = _noTrees ?
			probVerts_size + 2 <= _maxScore :
			probVerts_size + 2 < _maxScore;

		if (probVerts_size == 0 || not_enough_verts)
			continue;

		// Too many vertices for hardware?
		if (_useHW && probVerts_size > HW_MAX_VERTS)
		{
			cout << "Warning: too many vertices. Falling back to software mode." << endl;
			_useHW = false;
		}

		// Build adjacency matrix
		bool worthyProblem = buildMatrix(i_a, r_min, r_max);
		if (!worthyProblem)
			continue;

		// Sort vertices
		worthyProblem = sortVertices();
		if (!worthyProblem)
			continue;

		// Copy ordered problem vertices into set vertex buffer
		for (int v = 0; v < probVerts_size; v++)
		{
			_setVerts.push_back(_probVerts[_probOrder[v]]);
		}

		// Write to problem input data buffer
		appendDataBuffer(&last_header_idx);

		// Create an entry for this problem
		prob_desc.j = j;
		prob_desc.start_idx = start_idx;
		start_idx += probVerts_size;

		_setProbs.push_back(prob_desc);
	}

	if (_useHW && _setProbs.size() > HW_MAX_PROBS)
	{
		cout << "Warning: too many problems in set. Falling back to software mode." << endl;
		_useHW = false;
	}

	// If at least one problem exists, then there exists a last problem in this set.
	// Modify the last problem's header to include the bit that says "this is the last problem"
	if (last_header_idx >= 0)
	{
		_setInputBuffer[last_header_idx] |= 0x80000000;
	}
}


void CHardAlgorithm::copyCliques(int i_a, int i_b)
{
	PROFILE_FUNC();

	int data_ptr = 0;

	// Get number of cliques and clique size from the output buffer
	unsigned int n_cliques = _setOutputBuffer[data_ptr++];
	unsigned int maxsize = _setOutputBuffer[data_ptr++];

	// Maxscore is maxsize + 2 (the i and j vertices chosen)
	int new_maxscore = maxsize + 2;
	if (new_maxscore < _maxScore)
	{
		return;
	}
	else if (new_maxscore > _maxScore)
	{
		_maxScore = new_maxscore;
		_resultVec.clear();
	}

	ProteinPairVec matchVec(_maxScore);
	matchVec[0] = ProteinPair(i_a, i_b);

	unsigned int cliques_to_copy = n_cliques;
	if (_maxTrees >= 0)
	{
		//cliques_to_copy = std::min(cliques_to_copy, (unsigned int)_resultVec.size() - (unsigned int)_maxTrees);
		cliques_to_copy = std::min(cliques_to_copy, (unsigned int) (_maxTrees - _resultVec.size())); //abezgino 12-01-12: attempted fix for the bug where maxtrees parameter gets ignored. Alex can't remember the actual reason, so for now we consider that the bug is fixed.
	}

	for (unsigned int i = 0; i < cliques_to_copy; i++)
	{
		// Each clique contains: problem number, then vertex IDs.
		// Get the problem number, and then use it as an index into
		// setVerts, which contains the mapping from hardware-land 
		// vertex IDs to the numberLine array that contains all the
		// unique vertices for the problem set.

		int prob_no = _setOutputBuffer[data_ptr++];
		int prob_offs = _setProbs[prob_no].start_idx;
		int j_vert = _setProbs[prob_no].j;

		matchVec[1].first = _numberLine[j_vert].dest_a;
		matchVec[1].second = _numberLine[j_vert].dest_b;

		for (unsigned int j = 0; j < maxsize; j++)
		{
			int v_offs = _setOutputBuffer[data_ptr++];
			int v = _setVerts[prob_offs + v_offs];

			matchVec[j+2].first = _numberLine[v].dest_a;
			matchVec[j+2].second = _numberLine[v].dest_b;
		}

		_resultVec.push_back(matchVec);
	}
}


void CHardAlgorithm::maxCliquesHard()
{
	PROFILE_FUNC();

	if (!hw_write(_setInputBuffer))
	{
		cout << "Error writing to hardware" << endl;
		exit(1);
	}
	
	if (!hw_read(_setOutputBuffer))
	{
		cout << "Error reading from hardware" << endl;
		exit(1);
	}
}


int CHardAlgorithm::launchAnalysis
(
	const int wkRangeStart, 
	const int wkRangeEnd
)
{
	PROFILER_DESTROY();
	PROFILE_BEGIN(MMM);

	_maxScore = _minSize;
	_resultVec.clear();

	#ifdef DUMP_TESTBENCH
	ofstream tbfile("tb_data.txt");
	unsigned int n_sets = 0;
	#endif

	#ifdef CACHE_SIM
	hw_sim_cache_init();
	#endif

	#ifdef HW_COUNTERS
	hw_dbg_clear_counters();
	#endif

	#ifdef EDGE_STATS
	s_total_edges = 0;
	s_total_verts = 0;
	s_total_sq_verts = 0;
	#endif

	// Outer loop - iterate over all protein pairs (vertices)
	for (int i_a = wkRangeStart; i_a < wkRangeEnd ; ++i_a)
	{	if (!satisfiesReqTax(i_a)) //abezgino: only need to test one protein, since reqTax is meaningless when taxon matching is not enforced
			break; // *** if matrix 1 has all of its reqTaxName at the top, then there's nothing left to do
		for (int i_b = 0 ; i_b < _num2; ++i_b)
		{
			if (!validProteinPair(i_a, i_b))
				continue;
			
			// Build the numberline for protein pair i
			buildNumberLine(i_a, i_b);

			// Build a problem set
			buildProblemSet(i_a, i_b);

			if (_setProbs.empty())
				continue;

			_useHW? maxCliquesHard() : maxCliquesEmu();

			#ifdef DUMP_TESTBENCH
			dumpSet(tbfile);
			n_sets++;
			#endif

			copyCliques(i_a, i_b);
		} // i_b
	} // i_a

	#ifdef DUMP_TESTBENCH
	tbfile.close();
	cout << "Dumped " << n_sets << " sets\n";
	#endif

	PROFILE_END();
	PROFILER_UPDATE();

	return _maxScore;
}


void CHardAlgorithm::dumpSet(ofstream& file)
{
	for (unsigned int i = 0; i < _setInputBuffer.size(); i++)
	{
		file << _setInputBuffer[i] << endl;
	}

	for (unsigned int i = 0; i < _setOutputBuffer.size(); i++)
	{
		file << _setOutputBuffer[i] << endl;
	}
}


void CHardAlgorithm::calcScores(WeightedProteinPairVecVec& out)
{
	PROFILE_FUNC();

	unsigned int ncliques = _resultVec.size();
	out.reserve(ncliques);
	
	for (ProteinPairVecVec::iterator i = _resultVec.begin(); i != _resultVec.end(); ++i)
	{
		ProteinPairVec& clique = *i;

		unsigned int cliquesize = clique.size();
		unsigned int cliquesize_c2 = ( (cliquesize) * (cliquesize-1) / 2 );
		double score = 0.0;
		double mean = 0.0;

		// first find the mean ratio between the two cliques.
		for (unsigned int ai = 0; ai < cliquesize; ai++)
		{
			int a = clique[ai].first;
			int w = clique[ai].second;

			for (unsigned int bi = ai+1; bi < cliquesize; bi++)
			{
				int b = clique[bi].first;
				int x = clique[bi].second;
				double ab = getDist1(a,b);
				double wx = getDist2(w,x);

				mean += ab - wx;
			}
		}

		// there are (n choose 2) matching edge pairs
		mean /= (double)cliquesize_c2;

		// then get the std dev of the ratios
		for (unsigned int ai = 0; ai < cliquesize; ai++)
		{
			int a = clique[ai].first;
			int w = clique[ai].second;

			for (unsigned int bi = ai+1; bi < cliquesize; bi++)
			{
				int b = clique[bi].first;
				int x = clique[bi].second;
				double ab = getDist1(a,b);
				double wx = getDist2(w,x);

				double val = ab - wx;
				val = val - mean;
				score += val*val;
			}
		}

		// the worst-case stddev happens when all the ratios lie at the extremes of an interval
		// of size D, and we normalize according to this
		double scoreMax = (double)(cliquesize_c2 * (_D/2.0)*(_D/2.0));
		score = 1.0 - _allow*sqrt(score / scoreMax);
		out.push_back(WeightedProteinPairVec(clique, score));
	}
}


void CHardAlgorithm::maxCliquesEmu()
{
	PROFILE_FUNC();

	// Assuming the clique buffer can hold 65536 entries, this lookup table
	// says how many cliques of each size can be held in the buffer.
	// It includes the problem number value that's part of each clique.
	//
	// Note: this is only here in order to have identical behavior to the hardware.
	// It should totally go away.
	static const unsigned int max_buf_cliques[32] = {65535, 32767, 21117, 16383, 13106, 10921, 9361, 
		8191, 7280, 6552, 5956, 5460, 5041, 4680, 4368, 4095, 3854, 3639, 3448, 3275, 3119, 2977, 2848,
		2729, 2620, 2519, 2426, 2339, 2258, 2183, 2113, 2047};

	// Points to current problem in input buffer
	unsigned int data_ptr = 1;

	// Current max clique problem in this set
	unsigned int prob_no = 0;
	
	// Initialize global max size and number of cliques
	// Leave first two words in output buffer for storing
	// the final maxsize and number of cliques.
	unsigned int global_maxsize = _setInputBuffer[0] & 0xFFFF;
	unsigned int num_cliques = 0;
	bool nocliques = (_setInputBuffer[0] & 0x80000000) != 0;
	_setOutputBuffer.resize(2);

	while (true)
	{
		// Get number of vertices in this problem.
		// Bit 31 is on if this is the last problem in the set.
		unsigned int nverts = _setInputBuffer[data_ptr++];
		bool last_problem = (nverts & 0x80000000) != 0;
		nverts &= 0x7FFFFFFF;
		
		// This is the number of 32-bit words needed to hold the adjacency matrix
		// assuming 1 bit per entry.
		unsigned int nwords = (nverts*nverts + 31) / 32;

		// Initialize MaxSizePerVertex and per-problem maxsize
		_mcMSPV.resize(nverts);
		unsigned int local_maxsize = 0;

		// Clear the simulated cache
		#ifdef CACHE_SIM
		hw_sim_cache_clear();
		#endif

		// Go through the vertices in reverse order
		for (int startv = nverts - 1; startv >= 0; startv--)
		{
			// Initialize clique stack
			unsigned int cur_size = 0;

			// Initialize pointer stack
			unsigned int depth = 0;

			// Initialize first vspan. It contains vertices startv to nverts-1
			unsigned int end_ptr = 0;
			unsigned int v_ptr = 0;
			unsigned int last_v_ptr = 1;
			for (unsigned int v = startv; v < nverts; v++)
			{	if (end_ptr >= _mcVStack.size())
				{
					#ifdef DEBUG
					cout << "Warning: _mcVStack: end_ptr " << end_ptr << " >= " << _mcVStack.size() << " Resizing ..." << endl;
					#endif
					_mcVStack.resize(_mcVStack.size()*2);
				}
				_mcVStack[end_ptr++] = v;
			}

			bool found_larger = false;

			// Recursive processing
			while(true)
			{
				// While the vspan is not empty
				while (v_ptr < last_v_ptr)
				{
					// First vertex of the vspan is called 'v'
					unsigned int v = _mcVStack[v_ptr];
					// It gets included in the clique
					
					
					if (cur_size >= _mcCliqueStack.size())
					{
						#ifdef DEBUG
						cout << "Warning: _mcCliqueStack: cur_size " << cur_size << " >= " << _mcCliqueStack.size() << " Resizing ..." << endl;
						#endif
						_mcCliqueStack.resize(_mcCliqueStack.size()*2);
					}
					_mcCliqueStack[cur_size++] = v;
					//_mcCliqueStack.at(cur_size++) = v;
					
					if (cur_size > local_maxsize)
						local_maxsize = cur_size;

					// The rest of the vspan now needs to be intersected
					// with the neighbour set of v, creating a new vspan
					unsigned int write_ptr = end_ptr;
					unsigned int v_ptr_new = end_ptr;
					unsigned int last_v_ptr_new = end_ptr;
					unsigned int n_children = 0;
					for (unsigned int u_ptr = v_ptr + 1; u_ptr < end_ptr; u_ptr++)
					{
						unsigned int u = _mcVStack[u_ptr];

						// Lookup matrix[v][u]
						unsigned int bit_addr = v*nverts + u;
						unsigned int word_addr = bit_addr >> 5;
						unsigned int word = _setInputBuffer[data_ptr + word_addr];
						unsigned int bit_offs = bit_addr & 0x1F;
						bool isconn = (word & (1 << bit_offs)) != 0;
					
						#ifdef CACHE_SIM
						hw_sim_cache_access(bit_addr);
						#endif

						if (isconn)
						{
							// If '1', then u gets added to the child vspan.
							n_children++;
							
							if (write_ptr >= _mcVStack.size())
							{
								#ifdef DEBUG
								cout << "Warning: _mcVStack: write_ptr " << write_ptr << " >= " << _mcVStack.size() << " Resizing ..." << endl;
								#endif
								_mcVStack.resize(_mcVStack.size()*2);
							}
							_mcVStack[write_ptr++] = u;
							
							// However, it may fail the branch and bound test. This
							// means it's not fit to have child vspans because
							// (current clique size) + (u plus the most children it could
							// possibly have), and therefore can't make big enough cliques.

							bool mspv_passed = nocliques?
								cur_size + _mcMSPV[u] > local_maxsize :
								cur_size + _mcMSPV[u] >= local_maxsize;

							if (mspv_passed)
								last_v_ptr_new++;
						}
					}

					#ifdef DUMP_TESTBENCH
					if (write_ptr >= HW_MAX_STACKSIZE)
					{
						cout << "Warning: Stack overflow will happen in hardware" << endl;
					}
					#endif

					// If some of the other vspan vertices survived the culling (creating a
					// nonempty child vspan), and if there's enough of them to be
					// worth it, then recurse downwards.
					bool enough_children = nocliques?
						n_children > 0 && cur_size + n_children > local_maxsize :
						n_children > 0 && cur_size + n_children >= local_maxsize;

					if (enough_children)
					{
						// Save current vstack pointers on the pointer stack, so this
						// vspan can be resumed later.
						
						if (depth*4+2 >= _mcPtrStack.size())
						{
							#ifdef DEBUG
							cout << "Warning: _mcPtrStack: depth " << depth << "*4+2 >= " << _mcPtrStack.size() << " Resizing ..." << endl;
							#endif
							_mcPtrStack.resize(_mcPtrStack.size()*2);
						}
						_mcPtrStack[depth*4 + 0] = v_ptr + 1;
						_mcPtrStack[depth*4 + 1] = end_ptr;
						_mcPtrStack[depth*4 + 2] = last_v_ptr;
						
						depth++;

						// Set new pointers
						v_ptr = v_ptr_new;
						last_v_ptr = last_v_ptr_new;
						end_ptr = write_ptr;
					}
					else
					{
						// No downward recursion possible, v is the last vertex
						// to be added to this clique. Record this clique if it
						// big enough.	
						if (cur_size >= global_maxsize)
						{
							// New largest clique found. Empty current list.
							if (cur_size > global_maxsize)
							{
								global_maxsize = cur_size;
								num_cliques = 0;
								_setOutputBuffer.resize(2);

								found_larger = true;
							}

							bool store_clique = nocliques? (num_cliques == 0) :
								(num_cliques < max_buf_cliques[global_maxsize]);

							if (store_clique)
							{
								// Stuff problem# + vertex IDs into output buffer
								_setOutputBuffer.push_back(prob_no);
								for (unsigned int c = 0; c < cur_size; c++)
								{
									_setOutputBuffer.push_back(_mcCliqueStack[c]);
								}

								num_cliques++;
							}
						}

						// In the nocliques case, once we find a clique larger than the previous
						// best size (and it's always guaranteed to be exactly 1 more vertex than the previous),
						// then we can stop looking any further.
						if (nocliques && found_larger)
							break;

						// Remove 'v' from the clique
						cur_size--;

						// Give the next vertex in the vspan the chance to be 'v'
						v_ptr++;
					}
				}

				// Break out of main recursion loop - see above.
				if (nocliques && found_larger)
					break;

				// All possibilities for this vspan have been tried. Time to exit one
				// level from recursion.
				if (depth == 0)
				{
					// At the top? Done this problem.
					break;
				}
				else
				{
					depth--;
					cur_size--;

					// Restore pointers from the ptr stack to resume
					// processing the parent vspan.
					v_ptr = _mcPtrStack[depth*4 + 0];
					end_ptr = _mcPtrStack[depth*4 + 1];
					last_v_ptr = _mcPtrStack[depth*4 + 2];
				}
			}

			// Store the size of the max cliques attainable when just
			// considering vertices startv to nverts-1
			_mcMSPV[startv] = local_maxsize;
		}

		// Advance to next max clique problem
		prob_no++;
		data_ptr += nwords;

		if (last_problem)
			break;
	}

	// Store the number of maximum cliques and their size in the output stream.
	_setOutputBuffer[0] = num_cliques;
	_setOutputBuffer[1] = global_maxsize;
}


int CHardAlgorithm::calcMatchPotential()
{
	PROFILE_FUNC();
	
	// Build vertices
	ProteinPairVec vertices(_num1 * _num2);
	int nverts = 0;
	for (int a = 0; a < _num1; a++)
	{
		for (int b = 0; b < _num2; b++)
		{
			if (validProteinPair(a,b))
			{
				vertices[nverts].first = a;
				vertices[nverts].second = b;
				nverts++;
			}
		}
	}

	_probMatrix.resize(nverts*nverts);
	_probDegree.resize(nverts);

	// Initialize degrees
	for (int i = 0; i < nverts; i++)
	{
		_probDegree[i] = 0;
	}

	// Build matrix and calculate degrees
	for (int i_idx = 0; i_idx < nverts; i_idx++)
	{
		int i_a = vertices[i_idx].first;
		int i_b = vertices[i_idx].second;

		// Set diagonal to 0
		_probMatrix[i_idx*(nverts+1)] = 0;

		for (int j_idx = i_idx+1; j_idx < nverts; j_idx++)
		{
			int j_a = vertices[j_idx].first;
			int j_b = vertices[j_idx].second;

			bool conn = getDist1(i_a, j_a) != TOO_SMALL && getDist2(i_b, j_b) != TOO_SMALL;
			if (conn)
			{
				_probMatrix[i_idx*nverts + j_idx] = 1;
				_probMatrix[j_idx*nverts + i_idx] = 1;
				_probDegree[i_idx]++;
				_probDegree[j_idx]++;
			}
			else
			{
				_probMatrix[i_idx*nverts + j_idx] = 0;
				_probMatrix[j_idx*nverts + i_idx] = 0;
			}
		}
	}

	// Sort vertices
	_probVerts.resize(nverts); // hack: sortVertices uses _probVerts.size
	_maxScore = 0;
	sortVertices();

	return _probNumColorClasses;	
}

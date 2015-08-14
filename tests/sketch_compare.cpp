#include "SequenceCodec.hpp"
#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "iomanip.hpp"
#include "create_smatrix.hpp"

#include "bio/fastx_iterator.hpp"
#include "bio/sequence_compare.hpp"
#include "jaz/hash.hpp"
#include "jaz/string.hpp"



template <typename Hash>
void shingles(const AppConfig& opt, const std::string& s, Hash h, std::vector<uint64_t>& sh) {
    int l = s.size();
    int k = opt.kmer;

    SequenceCodec sc(opt.is_dna);
    sh.resize(l - k + 1);

    for (int i = 0; i < l - k + 1; ++i) {
	// sh[i] = h(sc.code(s.substr(i, k)));
	sh[i] = h(s.substr(i, k));
    }
} // shingles


//XORshift
uint64_t xorshift64(uint64_t x) {
    x ^= x >> 12; // a
    x ^= x << 25; // b
    x ^= x >> 27; // c
    return x * static_cast<uint64_t>(2685821657736338717ULL);
}

//XORShift on shingles
void xorsketches (std::vector<uint64_t>& shingles, int NF, std::vector<uint64_t>& sketch) {
    sketch.resize(NF);
    for(int i = 0; i < NF; i++) {
        for(int s = 0; s < shingles.size(); s++) {
            shingles[s] = xorshift64(shingles[s]);
        }
        sketch[i]= *std::min_element(shingles.begin(), shingles.end());
    }
}//xorsketches

//resize sketch to NF
//partial sort to find 10th element
//copy the 10th element from shingles to sketch
//sort sktech
void sketches(std::vector<uint64_t>& shingles, int NF, std::vector<uint64_t>& sketch) {
    sketch.resize(NF);
    std::nth_element(shingles.begin(), shingles.begin() + NF, shingles.end());
    std::copy(shingles.begin(), shingles.begin() + NF, sketch.begin());
    std::sort(sketch.begin(), sketch.end());
} // sketches


void sketches2(const AppConfig& opt, std::vector<uint64_t>& shingles, int NF, std::vector<uint64_t>& sketch) {
    sketch.clear();
    for (int i = 0; i < shingles.size(); ++i) {
	if (shingles[i] % opt.mod == NF) sketch.push_back(shingles[i]);
    }
    std::sort(sketch.begin(), sketch.end());
} // sketches


template <typename Ostream>
void process(AppConfig opt, const SequenceList& SL, Ostream& of) {

    //NF = length of static sequence sketch
    int NF = opt.mod;
    jaz::murmur264 h;

    //Sketch vectors
    //contains all hashed sketches
    //unsigned 64 bit ints
    std::vector<uint64_t> sketch0;
    std::vector<uint64_t> sketch1;

    //Kmer vectors
    //contains all shingles from sketches in vectors after repeated hashing
    //unsigned 64 bit ints
    std::vector<uint64_t> shingles0;
    std::vector<uint64_t> shingles1;

    //Set n to size of input sequences
    int n = SL.seqs.size();
    n >>= 1;

    //comp provides exact shared kmer count
    bio::kmer_fraction comp(opt.kmer);

    //Compute semi global alignment    
    bio::scoring_matrix sm;
    int f, g;
    create_smatrix(opt.gaps, opt.is_dna, sm, f, g);
    bio::semi_global_alignment sgalign(sm, f, g);

    //Used for correlation
    double S0 = 0.0;
    double S1 = 0.0;
    double S2 = 0.0;
    double S01 = 0.0;
    double S02 = 0.0;
    double S0_2 = 0.0;
    double S1_2 = 0.0;
    double S2_2 = 0.0;

    //tuples of 64 bit int sketch and vector of sequence ID's it's contained in
    std::map<uint64_t, std::vector<int>> sketch2seq;

    //iterate through sequece strings and ID all kmers
    for (int i = 0; i < n; ++i) {
	const std::string& s0 = SL.seqs[2 * i].s;
	const std::string& s1 = SL.seqs[2 * i + 1].s;

	if ((s0.size() < opt.kmer) || (s1.size() < opt.kmer)) continue;

    //generate shingles for input sequences
    shingles(opt, s0, h, shingles0);
    shingles(opt, s1, h, shingles1);
    
    //sketches(shingles0, NF, sketch0);
    //sketches(shingles1, NF, sketch1);
	
    //sketches2(opt, shingles0, 0, sketch0);
    //sketches2(opt, shingles1, 0, sketch1);

    xorsketches(shingles0, NF, sketch0);
    xorsketches(shingles1, NF, sketch1);

	/*
	  std::copy(sketch0.begin(), sketch0.end(), std::ostream_iterator<uint64_t>(std::cout, " "));
	  std::cout << std::endl;

	  std::copy(sketch1.begin(), sketch1.end(), std::ostream_iterator<uint64_t>(std::cout, " "));
	  std::cout << std::endl;
	*/

	for (int j = 0; j < sketch0.size(); ++j) sketch2seq[sketch0[j]].push_back(2 * i);
	for (int j = 0; j < sketch1.size(); ++j) sketch2seq[sketch1[j]].push_back(2 * i + 1);

    //calculate containment
    int is = jaz::intersection_size(sketch0.begin(), sketch0.end(), sketch1.begin(), sketch1.end());
    double contain = (1.0 * is) / std::min(sketch0.size(), sketch1.size());

    //calculate jaccard
    int jc = jaz::intersection_size(sketch0.begin(), sketch0.end(), sketch1.begin(), sketch1.end());
    double  jacc = (1.0 * jc / ((sketch0.size() + sketch1.size())));

    //calculate actual shared kmer fraction
    auto ki = comp(s0, s1);
    double ksim = 1.0 * boost::get<0>(ki) / std::min(boost::get<1>(ki), boost::get<2>(ki));

    //calculate semi_global_alignment
    auto sga = sgalign(s0, s1);
    
    S0 += ksim;
    S1 += contain;
    S2 += jacc;
    S01 += (ksim * contain);
    S02 += (ksim * jacc);
    S0_2 += (ksim * ksim);
    S1_2 += (contain * contain);
    S2_2 += (jacc * jacc);

    of << contain << " " << jacc << " " << boost::get<0>(sga) << " " << boost::get<1>(sga) << " " << boost::get<2>(sga) << std::endl;

    }//for (int i = 0; i < n; ++i)

    //for (auto iter = sketch2seq.begin(); iter != sketch2seq.end(); ++iter) {
	//   of << iter->second.size() << std::endl;
    //}


    //correlation for jaccard
    double cor0 = (S02 - ((S0 * S2) / n)) / std::sqrt((S0_2 - ((S0 * S0) / n)) * (S2_2 - ((S2 * S2) /n )));
    //correlation for containment
    double cor1 = (S01 - ((S0 * S1) / n)) / std::sqrt((S0_2 - ((S0 * S0) / n)) * (S1_2 - ((S1 * S1) / n)));
    

    //#pragma omp critical
    //of << opt.kmer << " " << opt.mod << " " << cor << std::endl;
}; // process


int main(int argc, char* argv[]) {
    std::map<std::string, std::string> conf;

    bool res = false;
    int pos = 0;

    boost::tie(res, pos) = jaz::parse_argv(argc, argv, conf);
    if (res == false) return -1;

    AppConfig opt;
    std::string err = "";

    boost::tie(res, err) = opt.set(conf);
    if (res == false) {
	std::cout << "error: " << err << std::endl;
	return -1;
    }

    AppLog log;
    Reporter report(std::cout, std::cout);

    SequenceList SL;

    std::ifstream f(opt.input.c_str());
    if (!f) return -1;

    pos = 0;
    bio::fasta_input_iterator<> fi(f), end;

    Sequence s;

    for (; fi != end; ++fi, ++pos) {
	s.id = pos;
	s.s = fi->second;
	SL.seqs.push_back(s);
    }

    SL.N = SL.seqs.size();

    f.close();

    std::ofstream of(opt.output.c_str());

    /*
#pragma omp parallel
    {

#pragma omp single nowait
	for (int k = 6; k < 20; ++k) {
	    for (int m = 8; m < 25; ++m) {

#pragma omp task private(opt)
		{
		    opt.kmer = k;
		    opt.mod = m;
		    process(opt, SL, of);
		}

	    } // for m
	} // for k
#pragma omp barrier

    }
    */

    //set variables
    opt.kmer = 10;
    opt.mod = 10;

    //call process with opt, sequence list, and outputstream
    //opt = OPTIONS, stores all variables for run
    //ie. kmer, mod, iterations, t...
    process(opt, SL, of);

    //oputstream
    of.close();

    return 0;
} // main

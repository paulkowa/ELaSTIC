#include "SequenceCodec.hpp"
#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "iomanip.hpp"

#include "bio/fastx_iterator.hpp"
#include "bio/sequence_compare.hpp"
#include "jaz/hash.hpp"


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
    const int NF = 10;
    jaz::murmur264 h;

    std::vector<uint64_t> sketch0;
    std::vector<uint64_t> sketch1;

    std::vector<uint64_t> shingles0;
    std::vector<uint64_t> shingles1;

    int n = SL.seqs.size();
    n >>= 1;

    bio::kmer_fraction comp(opt.kmer);

    double S0 = 0.0;
    double S1 = 0.0;
    double S01 = 0.0;
    double S0_2 = 0.0;
    double S1_2 = 0.0;

    std::map<uint64_t, std::vector<int>> sketch2seq;

    for (int i = 0; i < n; ++i) {
	const std::string& s0 = SL.seqs[2 * i].s;
	const std::string& s1 = SL.seqs[2 * i + 1].s;

	if ((s0.size() < opt.kmer) || (s1.size() < opt.kmer)) continue;

	shingles(opt, s0, h, shingles0);
	shingles(opt, s1, h, shingles1);

	// sketches(shingles0, NF, sketch0);
	// sketches(shingles1, NF, sketch1);

	sketches2(opt, shingles0, 0, sketch0);
	sketches2(opt, shingles1, 0, sketch1);

	/*
	  std::copy(sketch0.begin(), sketch0.end(), std::ostream_iterator<uint64_t>(std::cout, " "));
	  std::cout << std::endl;

	  std::copy(sketch1.begin(), sketch1.end(), std::ostream_iterator<uint64_t>(std::cout, " "));
	  std::cout << std::endl;
	*/

	for (int j = 0; j < sketch0.size(); ++j) sketch2seq[sketch0[j]].push_back(2 * i);
	for (int j = 0; j < sketch1.size(); ++j) sketch2seq[sketch1[j]].push_back(2 * i + 1);

	int is = jaz::intersection_size(sketch0.begin(), sketch0.end(), sketch1.begin(), sketch1.end());
	double msim = (1.0 * is) / std::min(sketch0.size(), sketch1.size());

	auto ki = comp(s0, s1);
	double ksim = 1.0 * boost::get<0>(ki) / std::min(boost::get<1>(ki), boost::get<2>(ki));

	S0 += ksim;
	S1 += msim;
	S01 += (ksim * msim);
	S0_2 += (ksim * ksim);
	S1_2 += (msim * msim);

	// of << ksim << " " << msim << std::endl;
    }

    for (auto iter = sketch2seq.begin(); iter != sketch2seq.end(); ++iter) {
	of << iter->second.size() << std::endl;
    }

    double cor = (S01 - ((S0 * S1) / n)) / std::sqrt((S0_2 - ((S0 * S0) / n)) * (S1_2 - ((S1 * S1) / n)));
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

    opt.kmer = 10;
    opt.mod = 10;
    process(opt, SL, of);

    of.close();

    return 0;
} // main

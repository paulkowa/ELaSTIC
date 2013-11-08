#include "SequenceDB.hpp"
#include "config_log.hpp"
#include "iomanip.hpp"
#include "make_shingles.hpp"

#include "bio/fastx_iterator.hpp"


inline unsigned int hash_sketch(const sketch_id& si) {
    uint64_t x = si.sketch;
    x >>= 3;
    return (x ^ (x >> 10) ^ (x >> 20)) & 0xFFFFFFFF;
} // hash_sketch


int main(int argc, char* argv[]) {
    std::map<std::string, std::string> conf;

    bool res = false;
    int pos = 0;

    boost::tie(res, pos) = jaz::parse_argv(argc, argv, conf);
    if (res == false) return -1;

    AppConfig opt;
    std::string err = "";

    boost::tie(res, err) = opt.set(conf);
    if (res == false) return -1;

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

    std::cout << opt << std::endl;
    std::cout << "sequences: " << SL.N << std::endl;

    jaz::murmur264 hash;
    std::vector<shingle_list_type> shingles;

    boost::tie(res, err) = make_shingles(opt, log, report, SL, hash, shingles);

    int M = 0;
    int p = 1024;

    std::vector<unsigned int> hist(p, 0);

    for (int i = 0; i < shingles.size(); ++i) {
	for (int j = 0; j < shingles[i].size(); ++j) {
	    if (shingles[i][j] % opt.mod == M) {
		sketch_id sid = make_sketch_id(shingles[i][j], i, 0);
		hist[hash_sketch(sid) % p]++;
	    }
	}
    }

    for (int i = 0; i < p; ++i) {
	std::cerr << i << " " << hist[i] << std::endl;
    }

    return 0;
} // main

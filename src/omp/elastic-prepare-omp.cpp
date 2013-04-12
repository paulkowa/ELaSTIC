/***
 *  $Id$
 **
 *  File: elastic-prepare-omp.cpp
 *  Created: May 20, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012-2013 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying LICENSE.
 *
 *  This file is part of ELaSTIC.
 */

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include "SequenceCodec.hpp"
#include "config.hpp"
#include "iomanip.hpp"
#include "tools.hpp"

#include <arpa/inet.h>

#include <boost/lexical_cast.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/bzip2.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/tuple/tuple.hpp>

#include <bio/fastx_iterator.hpp>

#include <jaz/hash.hpp>
#include <jaz/parameters.hpp>
#include <jaz/boost/files.hpp>
#include <jaz/science/disjoint_set.hpp>


std::istream* open_stream(const fs::path& name,
			  std::ifstream& fs, boost::iostreams::filtering_istream& cs) {
    std::istream* is = 0;

    if ((name.extension() != ".bz2") && (name.extension() != ".gz")) {
	fs.open(name.string().c_str());
	if (fs) is = &fs;
    } else {
	fs.open(name.string().c_str(), std::ios_base::binary);
	if (fs) {
	    if (name.extension() == ".gz") cs.push(boost::iostreams::gzip_decompressor());
	    else cs.push(boost::iostreams::bzip2_decompressor());
	    cs.push(fs);
	}
	is = &cs;
    }

    return is;
} // open_stream


struct AppConfig {
    AppConfig() {
	input = "";
	output = "";
	dna = true;
	length = 100;
	clean = true;
	group = true;
	kmer1 = 20;
	kmer2 = 15;
    } // AppConfig

    static void usage() {
	std::cout << "Usage: " << ELASTIC_PREPARE_SHORT << " --input name --output name [options...]\n";
	std::cout << "\n";
	std::cout << "Options:\n";
	std::cout << "  --input name          read input from this file/directory\n";
	std::cout << "  --output name         write output to files with this prefix\n";
	std::cout << "  --type {nt|aa}        set input sequence type (default nt)\n";
	std::cout << "  --length size         remove sequences shorter than this length (default 100)\n";
	std::cout << "  --clean {0|1}         remove sequences with missing bases (default 1)\n";
	std::cout << "  --group {0|1}         cluster very high similarity sequences (default 1)\n";
	std::cout << "\n";
    } // usage

    template <typename Container>
    std::pair<bool, std::string> set(const Container& conf) {
	std::string val;

	// check major options
	if (jaz::check_option(conf, "input", val) == false) {
	    return std::make_pair(false, "missing input parameter");
	}

	input = val;

	if (jaz::check_option(conf, "output", val) == false) {
	    return std::make_pair(false, "missing output parameter");
	}

	output = val;

	if (jaz::check_option(conf, "type", val) == true) {
	    if ((val != "nt") && (val != "aa")) {
		return std::make_pair(false, "incorrect sequence type");
	    }
	    if (val == "aa") {
		dna = false;
		kmer1 = 5;
		kmer2 = 3;
	    }
	}

	try {

	    // check hidden options
	    if (jaz::check_option(conf, "kmer1", val) == true) {
		kmer1 = boost::lexical_cast<unsigned short int>(val);
		if ((kmer1 < 5) || (kmer1 > 31)) {
		    return std::make_pair(false, "incorrect size of kmer1");
		}
	    }

	    if (jaz::check_option(conf, "kmer2", val) == true) {
		kmer2 = boost::lexical_cast<unsigned short int>(val);
		if ((kmer2 < 3) || (kmer2 > kmer1)) {
		    return std::make_pair(false, "incorrect size of kmer2");
		}
	    }

	    // check other options
	    if (jaz::check_option(conf, "clean", val) == true) {
		clean = boost::lexical_cast<unsigned int>(val);
	    }

	    if (jaz::check_option(conf, "group", val) == true) {
		group = boost::lexical_cast<unsigned int>(val);
	    }

	    if (jaz::check_option(conf, "length", val) == true) {
		length = boost::lexical_cast<unsigned short int>(val);
		if (group == true) {
		    if ((length < kmer1) || (length - kmer1 < kmer2)) {
			return std::make_pair(false, "incorrect minimal sequence length");
		    }
		}
	    }

	} catch (boost::bad_lexical_cast& ex) {
	    return std::make_pair(false, "incorrect argument(s)");
	}

	return std::make_pair(true, "");
    } // set

    std::string input;
    std::string output;
    unsigned short int length;
    bool dna;
    bool clean;
    bool group;
    unsigned short int kmer1;
    unsigned short int kmer2;

    friend std::ostream& operator<<(std::ostream& os, const AppConfig& opt) {
	os << "input = " << opt.input << "\n";
	os << "output = " << opt.output << "\n";
	os << "dna = " << opt.dna << "\n";
	os << "length = " << opt.length << "\n";
	os << "clean = " << opt.clean << "\n";
	os << "group = " << opt.group << "\n";
	os << "kmer1 = " << opt.kmer1 << "\n";
	os << "kmer2 = " << opt.kmer2;
	return os;
    } // operator<<
}; // struct AppConfig


struct AppLog {
    AppLog() : argv(), wtime(0), input(0), extracted(0), groups(0) {
	time_t t;
	time(&t);
	date = ctime(&t);
    } // AppLog

    std::string date;
    std::string argv;
    double wtime;
    unsigned int input;
    unsigned int extracted;
    unsigned int groups;

    friend std::ostream& operator<<(std::ostream& os, const AppLog& log) {
	os << "execution date: " << log.date;
	os << "program version: " << ELASTIC_PREPARE_SHORT << " " << ELASTIC_PREPARE_VERSION << "\n";
	os << "program options: " << log.argv << "\n";
	os << "walltime used: " << log.wtime << "\n";
	os << "input sequences: " << log.input << "\n";
	os << "extracted sequences: " << log.extracted << "\n";
	os << "output groups: " << log.groups << "\n";
	return os;
    } // operator<<
}; // struct AppLog


struct SequenceDesc {
    uint64_t sketch;
    std::string name;
    unsigned short int size;
}; // struct SequenceDesc

class sequence_compare {
public:
    explicit sequence_compare(std::vector<SequenceDesc>& sdesc) : sdesc_(sdesc) { }

    bool operator()(unsigned int lhs, unsigned int rhs) const {
	return sdesc_[lhs].size < sdesc_[rhs].size;
    } // operator()

private:
    std::vector<SequenceDesc>& sdesc_;

}; // class sequence_compare


typedef bio::fasta_input_iterator<>::value_type sequence_type;


class is_clean {
public:
    is_clean(bool check, unsigned int size, bool dna = true) : check_(check), size_(size) {
	std::memset(sig_, 0, 256);
	if (dna == true) {
	    sig_['A'] = sig_['a'] = 1;
	    sig_['C'] = sig_['c'] = 1;
	    sig_['G'] = sig_['g'] = 1;
	    sig_['T'] = sig_['t'] = 1;
	} else {
	    sig_['A'] = sig_['a'] = 1;
	    sig_['C'] = sig_['c'] = 1;
	    sig_['D'] = sig_['d'] = 1;
	    sig_['E'] = sig_['e'] = 1;
	    sig_['F'] = sig_['f'] = 1;
	    sig_['G'] = sig_['g'] = 1;
	    sig_['H'] = sig_['h'] = 1;
	    sig_['I'] = sig_['i'] = 1;
	    sig_['K'] = sig_['k'] = 1;
	    sig_['L'] = sig_['l'] = 1;
	    sig_['M'] = sig_['m'] = 1;
	    sig_['N'] = sig_['n'] = 1;
	    sig_['P'] = sig_['p'] = 1;
	    sig_['Q'] = sig_['q'] = 1;
	    sig_['R'] = sig_['r'] = 1;
	    sig_['S'] = sig_['s'] = 1;
	    sig_['T'] = sig_['t'] = 1;
	    sig_['V'] = sig_['v'] = 1;
	    sig_['W'] = sig_['w'] = 1;
	    sig_['Y'] = sig_['y'] = 1;
	}
    } // is_clean

    bool operator()(const sequence_type& s) const {
	unsigned int n = s.second.size();
	if (n < size_) return false;
	if (check_ == true) {
	    for (unsigned int i = 0; i < n; ++i) {
		if (sig_[s.second[i]] == 0) return false;
	    }
	}
	return true;
    } // operator

private:
    char sig_[256];

    bool check_;
    unsigned int size_;

}; // class is_clean


typedef std::vector<unsigned int> int_list_type;

struct Cluster {
    unsigned int id;
    int_list_type seqs;
}; // struct Cluster


inline bool id_compare(const Cluster& lhs, const Cluster& rhs) {
    return lhs.id < rhs.id;
} // id_compare

inline bool seq_compare(const Cluster& lhs, const Cluster& rhs) {
    return lhs.seqs[0] < rhs.seqs[0];
} // seq_compare


void clean_sequence(std::string& s, bool dna = true) {
    const char NT[] = "ACGT";
    const char AA[] = "ACDEFGHIKLMNPQRSTVWY";

    unsigned int l = s.size();
    for (unsigned int i = 0; i < l; ++i) {
	s[i] = std::toupper(s[i]);
	if (dna == true) {
	    if (std::find(NT, NT + sizeof(NT), s[i]) == false) s[i] = 'T';
	} else {
	    if (std::find(AA, AA + sizeof(AA), s[i]) == false) s[i] = 'A';
	}
    }
} // clean_sequence


uint64_t super_shingle(const std::string& s, unsigned int kmer1, unsigned int kmer2) {
    // stage 1
    unsigned int l = s.size();
    std::vector<uint64_t> sh1(l - kmer1 + 1, 0);

    jaz::murmur264 hash;

    for (unsigned int i = 0; i < l - kmer1 + 1; ++i) {
	sh1[i] = hash(std::string(s.begin() + i, s.begin() + i + kmer1));
    }

    std::sort(sh1.begin(), sh1.end());

    // stage 2
    l = sh1.size();
    std::vector<uint64_t> sh2(l - kmer2 + 1, 0);

    unsigned int vsz = kmer2 * sizeof(uint64_t);

    for (unsigned int i = 0; i < l - kmer2 + 1; ++i) {
	const char* v = reinterpret_cast<const char*>(&sh1[i]);
	sh2[i] = hash(v, vsz);
    }

    return *std::min_element(sh2.begin(), sh2.end());
} // super_shingle


void welcome() {
    std::cout << ELASTIC_PREPARE_FULL << " Version " << ELASTIC_PREPARE_VERSION << "\n";
    std::cout << ELASTIC_PREPARE_COPYRIGHT << "\n";
    std::cout << "\n";
} // welcome


std::pair<bool, std::string> run(const AppConfig& opt, AppLog& log, Reporter& report) {
    double t0 = get_time();

    // pre-open output files
    std::ofstream fmap((opt.output + ".emap").c_str());
    if (!fmap) return std::make_pair(false, "unable to create " + opt.output + ".emap");

    std::ofstream fdel((opt.output + ".edel").c_str());
    if (!fdel) return std::make_pair(false, "unable to create " + opt.output + ".edel");

    std::ofstream fseq((opt.output + ".eseq").c_str(), std::ios_base::binary);
    if (!fseq) return std::make_pair(false, "unable to create " + opt.output + ".eseq");

    std::ofstream fidx((opt.output + ".eidx").c_str(), std::ios_base::binary);
    if (!fidx) return std::make_pair(false, "unable to create " + opt.output + ".eidx");

    std::ofstream flog((opt.output + ".eplog").c_str());
    if (!flog) return std::make_pair(false, "unable to create " + opt.output + ".eplog");

    // get files to process
    std::cout << step << "scanning " << opt.input << " for input files..." << std::endl;

    std::vector<fs::path> files;
    std::vector<fs::path> cfiles; // clean files

    if (jaz::files(opt.input, std::back_inserter(files)) == false) {
	return std::make_pair(false, "unable to scan " + opt.input);
    }

    if (files.empty() == true) return std::make_pair(false, "no files to process");

    report << info << "found " << files.size() << " input file(s)" << std::endl;
    report << step << "extracting sequences..." << std::endl;

    // GET INPUT SEQUENCES
    const unsigned int SBUF_SIZE = 512 * 1024;

    std::vector<sequence_type> seqs; // sequences to process
    std::vector<std::string> delseq; // removed sequences

    std::vector<SequenceDesc> sdesc;

    bool is_last = false;

    for (unsigned int i = 0; i < files.size(); ++i) {
	if (fs::is_regular_file(files[i]) == false) continue;

	std::ifstream fs;
	boost::iostreams::filtering_istream cs;

	std::istream* is = open_stream(files[i], fs, cs);

	if (is == 0) {
	    report.critical << warning << "unable to open " << files[i].string() << ", ignoring" << std::endl;
	    continue;
	}

	cfiles.push_back(files[i]);

	bio::fasta_input_iterator<> fi(*is), end;

	for (;;) {
	    if (is_last == false) {
		seqs.push_back(*fi);
		++fi;
	    }

	    // if buffer is full or we are done processing
	    if ((seqs.size() == SBUF_SIZE) || (is_last == true)) {
		// remove dirty sequences
		unsigned int m = seqs.size();

		std::vector<sequence_type>::iterator iter;
		iter = std::stable_partition(seqs.begin(), seqs.end(), is_clean(opt.clean, opt.length, opt.dna));

		for (std::vector<sequence_type>::iterator j = iter; j != seqs.end(); ++j) {
		    delseq.push_back(j->first);
		}

		seqs.erase(iter, seqs.end());

		// update map
		unsigned int pos = sdesc.size();
		unsigned int n = seqs.size();
		sdesc.resize(pos + n);

		report << info << "valid " << n << " out of " << m << " sequences" << std::endl;
		log.input += m;
		log.extracted += n;

#pragma omp parallel for schedule(dynamic) shared(opt, seqs, sdesc, pos)
		for (unsigned int j = 0; j < n; ++j) {
		    sdesc[pos + j].name = seqs[j].first;
		    sdesc[pos + j].size = seqs[j].second.size();
		    if (opt.group == true) {
			sdesc[pos + j].sketch = super_shingle(seqs[j].second, opt.kmer1, opt.kmer2);
		    }
		}

		seqs.clear();
	    } // if

	    // trick to handle the last buffer
	    if (is_last == true) break;

	    if (fi == end) {
		if (i + 1 == files.size()) is_last = true;
		else break;
	    }

	} // for fi
    } // for i

    std::vector<sequence_type>().swap(seqs);

    // PERFORM CLUSTERING
    if (sdesc.empty() == true) return std::make_pair(false, "no sequences found");

    unsigned int n = sdesc.size();

    report << step << "indexing/clustering " << n << " sequences..." << std::endl;

    // union find clustering
    std::vector<unsigned int> uf(n);
    jaz::set_make(&uf[0], n);

    if (opt.group == true) {
	std::map<uint64_t, int_list_type> sketch2ids;
	for (unsigned int i = 0; i < n; ++i) sketch2ids[sdesc[i].sketch].push_back(i);

	std::map<uint64_t, int_list_type>::iterator iter(sketch2ids.begin());

	for (; iter != sketch2ids.end(); ++iter) {
	    unsigned int l = iter->second.size();
	    for (unsigned int i = 1; i < l; ++i) {
		jaz::set_union(&uf[0], iter->second[0], iter->second[i]);
	    }
	} // for iter
    } // if opt.group

    // extract clusters
    Cluster c;
    std::vector<Cluster> cluster;

    cluster.reserve(n);

    for (unsigned int i = 0; i < n; ++i) {
	c.id = jaz::set_find(&uf[0], i);
	std::vector<Cluster>::iterator cfirst, clast;
	boost::tie(cfirst, clast) = std::equal_range(cluster.begin(), cluster.end(), c, id_compare);
	if (cfirst == clast) cfirst = cluster.insert(cfirst, c);
	cfirst->seqs.push_back(i);
    } // for i

    log.groups = cluster.size();
    report << info << "created " << cluster.size() << " output groups" << std::endl;

    // find longest sequence for each cluster
    // and move it to head
    std::vector<Cluster>::iterator iter(cluster.begin());

    for (; iter != cluster.end(); ++iter) {
	int_list_type& lst = iter->seqs;
	std::swap(*lst.begin(), *std::max_element(lst.begin(), lst.end(), sequence_compare(sdesc)));
    }

    // now we can sort cluster based on the representative sequence
    // and assign ids to sequences representing clusters
    std::sort(cluster.begin(), cluster.end(), seq_compare);

    unsigned int sid = 0;

    std::map<std::string, unsigned int> name2id;

    for (iter = cluster.begin(); iter != cluster.end(); ++iter, ++sid) {
	bool res;
	std::map<std::string, unsigned int>::iterator niter;

	boost::tie(niter, res) = name2id.insert(std::make_pair(sdesc[iter->seqs[0]].name, sid));

	if (res == false) return std::make_pair(false, "sequence with this name " + niter->first + " exists");
    } // for iter

    // WRITE OUTPUT
    report << step << "writing output files..." << std::endl;

    // sequence map
    unsigned int cid = 0;

    for (iter = cluster.begin() ; iter != cluster.end(); ++iter, ++cid) {
	fmap << cid << " " << iter->seqs.size() << std::endl;
	for (unsigned int i = 0; i < iter->seqs.size(); ++i) {
	    fmap << sdesc[iter->seqs[i]].name << std::endl;
	}
    }

    fmap.close();

    // deleted sequences
    std::copy(delseq.begin(), delseq.end(), std::ostream_iterator<std::string>(fdel, "\n"));
    fdel.close();

    // create index and store sequences
    SequenceCodec sc(opt.dna);

    unsigned int pos = 0;
    std::vector<uint16_t> index(name2id.size(), 0);

    for (unsigned int i = 0; i < cfiles.size(); ++i) {
	std::ifstream fs;
	boost::iostreams::filtering_istream cs;

	std::istream* is = open_stream(cfiles[i], fs, cs);
	if (is == 0) return std::make_pair(false, "file " + cfiles[i].string() + " disappeared");

	bio::fasta_input_iterator<> fi(*is), end;

	for (; fi != end; ++fi) {
	    std::map<std::string, unsigned int>::iterator iter = name2id.find(fi->first);
	    if (iter != name2id.end()) {
		uint32_t id = iter->second;
		id = htonl(id);

		std::string s = fi->second;

		clean_sequence(s, opt.dna);
		s = sc.code(s);

		fseq.write(reinterpret_cast<char*>(&id), sizeof(id));
		fseq.write(s.c_str(), s.size());

		index[pos] = s.size() + sizeof(id);
		pos++;
	    }
	}
    } // for i

    fseq.close();

    // store index
    std::transform(index.begin(), index.end(), index.begin(), htons);
    fidx.write(reinterpret_cast<char*>(&index[0]), index.size() * sizeof(index[0]));
    fidx.close();

    log.wtime = get_time() - t0;

    // write final log
    flog << log;
    flog << "config:" << std::endl;
    flog << opt;
    flog.close();

    return std::make_pair(true, "");
} // run


int main(int argc, char* argv[]) {
    welcome();

    // get parameters
    std::map<std::string, std::string> conf;

    if (argc == 1) {
	AppConfig::usage();
	return 0;
    }

    bool res = false;
    int pos = 0;
    std::string err = "";

    boost::tie(res, pos) = jaz::parse_argv(argc, argv, conf);

    if (res == false) {
	if (pos == -1) {
	    AppConfig::usage();
	    std::cout << error << "incorrect command line arguments\n";
	    return -1;
	} else {
	    std::cout << error << "incorrect command line argument " << argv[pos] << "\n";
	    return -1;
	}
    }

    // create config and log
    AppConfig opt;
    boost::tie(res, err) = opt.set(conf);

    if (res == false) {
	std::cout << error << err << "\n";
	return -1;
    }

    AppLog log;
    log.argv = jaz::join(' ', argv + 1, argv + argc);

    Reporter report(std::cout, std::cout);

    boost::tie(res, err) = run(opt, log, report);

    if (res == false) {
	std::cout << error << err << "\n";
	return -1;
    }

    std::cout << "time: " << log.wtime << std::endl;
    std::cout << "done!" << std::endl;
    std::cout << std::endl;

    return 0;
} // main

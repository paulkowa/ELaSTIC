/***
 *  $Id$
 **
 *  File: SequenceMap.hpp
 *  Created: May 20, 2012
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 *  Copyright (c) 2012 Jaroslaw Zola
 *  Distributed under the MIT License.
 *  See accompanying file LICENSE_MIT.txt.
 *
 *  This file is part of ELaSTIC.
 */

#ifndef SEQUENCE_MAP_HPP
#define SEQUENCE_MAP_HPP

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/tuple/tuple.hpp>


struct SequenceGroup {
    unsigned int id;
    // name[0] is the name of the representative sequence
    std::vector<std::string> name;
}; // struct SequenceGroup

inline bool operator==(const SequenceGroup& lhs, const SequenceGroup& rhs) {
    return (lhs.id == rhs.id);
} // operator==

inline bool operator<(const SequenceGroup& lhs, const SequenceGroup& rhs) {
    return (lhs.id < rhs.id);
} // operator<


class SequenceMap {
public:
    typedef std::vector<SequenceGroup>::iterator iterator;
    typedef std::vector<SequenceGroup>::const_iterator const_iterator;


    iterator begin() { return seqmap_.begin(); }

    const_iterator begin() const { return seqmap_.begin(); }

    iterator end() { return seqmap_.end(); }

    const_iterator end() const { return seqmap_.end(); }

    // std::pair<iterator, bool> insert(const SequenceGroup& obj) {
    // 	iterator first, last;

    // 	boost::tie(first, last) = std::equal_range(seqmap_.begin(), seqmap_.end(), obj);
    // 	if (first != last) return std::make_pair(first, false);

    // 	return std::make_pair(seqmap_.insert(first, obj), true);
    // } // insert

    const_iterator find(unsigned int id) const {
	SequenceGroup key;
	key.id = id;

	const_iterator first, last;

	boost::tie(first, last) = std::equal_range(seqmap_.begin(), seqmap_.end(), key);
	if (first == last) return end();

	return first;
    } // find


    bool read(const std::string& name) {
	std::ifstream f(name.c_str());
	if (!f) return false;

	seqmap_.clear();

	std::string buf;
	SequenceGroup obj;

	while (!f.eof()) {
	    f >> obj.id;

	    unsigned int size = 0;
	    f >> size;

	    if (!f.eof() && !f) return false;
	    seqmap_.push_back(obj);

	    std::getline(f, buf);

	    seqmap_.back().name.resize(size);
	    for (unsigned int i = 0; i < size; ++i) {
		std::getline(f, seqmap_.back().name[i]);
	    }
	} // while

	f.close();

	std::sort(seqmap_.begin(), seqmap_.end());
	seqmap_.erase(std::unique(seqmap_.begin(), seqmap_.end()), seqmap_.end());

	return true;
    } // read


    bool write(const std::string& name) const {
	std::ofstream f(name.c_str());
	if (!f) return false;

	unsigned int n = seqmap_.size();

	for (unsigned int i = 0; i < n; ++i) {
	    f << seqmap_[i].id << " " << seqmap_[i].name.size() << std::endl;
	    std::copy(seqmap_[i].name.begin(), seqmap_[i].name.end(),
		      std::ostream_iterator<std::string>(f, "\n"));
	}

	f.close();

	return true;
    } // write


private:
    std::vector<SequenceGroup> seqmap_;

}; // class SequenceMap

#endif // SEQUENCE_MAP_HPP

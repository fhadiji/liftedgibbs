/*
    Copyright (C) 2010
    Babak Ahmadi [babak dot ahmadi at iais dot fraunhofer dot de]
    Fabian Hadiji [fabian dot hadiji at iais dot fraunhofer dot de]
    Kristian Kersting (coordination) [kristian dot kersting at iais dot fraunhofer dot de]

    STREAM Project at
        Fraunhofer IAIS, Sankt Augustin, Germany, and
        KDML, Unversity of Bonn, Germany

    This file is part of libSTREAM.

    libSTREAM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    libSTREAM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, see <http://www.gnu.org/licenses/>.
*/
/*
 * FactorGraphProb.cpp
 *
 *  Created on: May 26, 2010
 *      Author: Fabian Hadiji
 */
#include <STREAM/FactorGraph.h>

using namespace dai;
using namespace std;

namespace stream {

template<>
void FactorGraph<CFactor,Var>::initCounts() {
	_counts.reserve( nrVars() );
	_edgeCounts.reserve( nrVars() );

	vector<vector<int> > varCounts;
	for (size_t i=0; i<nrVars(); i++) {
		varCounts.clear();
		varCounts.reserve(nbV(i).size());
		foreach(const Neighbor tmpFac, nbV(i)) {
			varCounts.push_back(vector<int>(nbF(tmpFac).size(),0));
			size_t pos = find(factor(tmpFac).sigma().begin(), factor(tmpFac).sigma().end(), tmpFac.dual) - factor(tmpFac).sigma().begin();
			varCounts.back()[pos] = 1;
		}
		_edgeCounts.push_back(varCounts);
		_counts.push_back(varCounts);
	}
}

template<>
FactorGraph<CFactor, Var>::FactorGraph(const CFactorGraph & cfg) : G(), _vars(), _factors() {
	create(cfg.factors());

	_counts.reserve( nrVars() );
	_edgeCounts.reserve( nrVars() );
	vector<vector<int> > varCounts;
	for (size_t i=0; i<nrVars(); i++) {
		varCounts.clear();
		varCounts.reserve(nbV(i).size());
		foreach(const Neighbor tmpFac, nbV(i)) {
			varCounts.push_back(vector<int>(cfg.factor(tmpFac).sigma().size(),0));
		}
		_edgeCounts.push_back(varCounts);
		_counts.push_back(varCounts);
	}

	map<size_t,int>::const_iterator countsIter;
	for (size_t i=0; i<cfg.nrFactors(); i++) {
		foreach (const Neighbor tmpVar, nbF(i)) {
			_edgeCounts[tmpVar][tmpVar.dual][tmpVar.iter] = cfg.factor(i).counts()[tmpVar.iter].size();
			for (countsIter = cfg.factor(i).counts()[tmpVar.iter].begin(); countsIter != cfg.factor(i).counts()[tmpVar.iter].end(); ++countsIter) {
				_counts[tmpVar][tmpVar.dual][countsIter->first] = countsIter->second;
			}
		}
	}

}

// taken from CFactorGraph implementation
template <>
void FactorGraph<CFactor, Var>::readFromFile(const char *filename) {
	CFactorGraph cfg;
	cfg.readFromFile(filename);
	*this = FactorGraph<CFactor, Var>(cfg.factors());
}

}

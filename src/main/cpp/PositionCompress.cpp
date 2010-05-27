/*
    Copyright (C) 2009
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
 * PositionCompress.cpp
 *
 *  Created on: Oct 2, 2009
 *      Author: Fabian Hadiji
 */
#include "STREAM/CBP/Compress.h"

namespace stream {

using namespace std;
using namespace dai;


void PositionCompress::initFacColors() {
	_facSigs = vector<Signature>(_cfg.nrFactors());
	vector<size_t> dims;
	vector<Real> correctPotential;
	boost::hash<vector<Real> > hashVector;
	for (size_t i=0; i<_cfg.nrFactors(); i++) {
		dims.clear();
		dims.reserve(_cfg.factor(i).vars().size());
		for (vector<Var>::const_iterator iter=_cfg.factor(i).vars().begin(); iter<_cfg.factor(i).vars().end(); iter++) {
			dims.push_back(iter->states());
		}
		Permute permute(dims, _cfg.factor(i).sigma());

		vector<Real> correctPotential(_cfg.factor(i).states());
		for (size_t k=0; k<_cfg.factor(i).states(); k++) {
			correctPotential[permute.convert_linear_index(k)] = _cfg.factor(i)[k];
		}

		_facSigs[i] = hashVector(correctPotential);
		_facInbox.push_back(std::vector<size_t>(_cfg.nbF(i).size() + 1) );
	}
}

size_t PositionCompress::facColor(size_t i) {
	for (size_t j=0; j< _cfg.nbF(i).size(); j++) {
		_facInbox[i][j] = _varSigs[_cfg.nbF(i,_cfg.getSigma(i)[j])];
	}
	_facInbox[i][_facInbox[i].size() - 1] = _facSigs[i];
	return hashVector(_facInbox[i]);
}


size_t PositionCompress::varColor(size_t i) {
	foreach (const dai::BipartiteGraph::Neighbor tmpFac, _cfg.nbV(i)) {
		vector<size_t> hashAndPos(2);
		size_t pos = find(_cfg.factor(tmpFac).sigma().begin(), _cfg.factor(tmpFac).sigma().end(), tmpFac.dual) - _cfg.factor(tmpFac).sigma().begin();

		hashAndPos[0] = _facSigs[tmpFac];
		hashAndPos[1] = pos;
		_varInbox[i][tmpFac.iter] = hashVector(hashAndPos);
	}
	_varInbox[i][_varInbox[i].size() - 1] = _varSigs[i];
	sort(_varInbox[i].begin(), _varInbox[i].end());

	return hashVector(_varInbox[i]);
}

}

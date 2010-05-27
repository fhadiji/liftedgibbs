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
 * Compress.cpp
 *
 *  Created on: Sep 29, 2009
 *      Author: Fabian Hadiji
 */

#include <STREAM/CBP/Compress.h>

namespace stream {

using namespace std;
using namespace dai;

void CompressInterface::initVarColors() {
	_varSigs = vector<Signature>(_cfg.nrVars(),0);
	for (size_t i=0; i<_cfg.vars().size(); i++) {
		map<long, size_t>::const_iterator evidence = _cfg.evidence().find(_cfg.var(i).label());
		if ( evidence != _cfg.evidence().end()) {
			_varSigs[i] = evidence->second;
		}
		_varInbox.push_back(vector<size_t>(_cfg.nbV(i).size() + 1));
	}
}

void CompressInterface::initFacColors() {
	_facSigs = vector<Signature>(_cfg.nrFactors());
	vector<Real> potential;
	boost::hash<vector<Real> > hashVector;
	for (size_t i=0; i<_cfg.nrFactors(); i++) {
		bool allTrue = false;
		for (size_t j=0; j<i; j++)  {
			// muss set _facColorVec[i] properly
			if (_cfg.factor(i).states() != _cfg.factor(j).states()) {
				continue;
			}

			allTrue = true;
			for (size_t k=0; k<_cfg.factor(i).states(); k++) {
				if(_cfg.factor(i)[k] != _cfg.factor(j)[k]) {
					allTrue = false;
					break;
				}
			}

			if (allTrue) {
				_facSigs[i] = _facSigs[j];
				break;
			}
		}

		if (!allTrue) {
			potential.clear();
			potential.reserve(_cfg.factor(i).states());
			for (size_t k=0; k<_cfg.factor(i).states(); k++) {
				potential.push_back(_cfg.factor(i)[k]);
			}
			_facSigs[i] = hashVector(potential);
		}
		_facInbox.push_back(std::vector<size_t>(_cfg.nbF(i).size() + 1) );
	}
}

void CompressInterface::execute() {

	if( _verbose >= 1 ) {
		cout << "Called " << name() << " execute method...";
	}
	if (_verbose >= 3) {
		cout << endl;
	}

	double tic = toc();

	init();

	if (_initialVarColors.size() > 0) {
		assert(_initialVarColors.size() == _varSigs.size());
		vector<size_t> testDict;
		for (size_t i=0; i<_initialVarColors.size(); i++) {
			_varSigs[i] = _initialVarColors[i];
		}
	}

	_iters = 0;
	do {
		iterate();
		_iters++;
		if( _verbose >= 3 ) {
			cout << name() << "::execute:  #varClusters " << _varClustering.size() << ", #facClusters " << _facClustering.size() << " after " << _iters << " passes (" << toc() - tic << " secs.)" <<  endl;
		}
	} while (!hasConverged());

	if( _verbose >= 1 ) {
		cout << "Finished in " << _iters << " iterations (" << toc() - tic << " secs.)" << endl;
	}
}

void CompressInterface::init() {
	if( _verbose >= 3 ) {
		cout << name() << "::init::initFacColors... ";
	}
	initFacColors();
	if( _verbose >= 3 ) {
		cout << endl;
	}
	if( _verbose >= 4 ) {
		dai::operator <<(cout, _facSigs) << endl;
	}

	if( _verbose >= 3 ) {
		cout << name() << "::init::initVarColors... ";
	}
	initVarColors();
	if( _verbose >= 3 ) {
		cout << endl;
	}
	if( _verbose >= 4 ) {
			dai::operator <<(cout, _varSigs) << endl;
	}

	_iters = 0;
}

void CompressInterface::iterate() {
	map<size_t, vector<size_t> >::iterator mapIter;
	_varClusterCountOld = _varClustering.size();
	_facClusterCountOld = _facClustering.size();

	// determine new colors of variables
	_varClustering.clear();
	for (size_t i=0; i<_cfg.nrVars(); i++) {
		size_t hashVal = varColor(i);
		mapIter = _varClustering.find(hashVal);
		if (mapIter == _varClustering.end()) {
			_varClustering[hashVal] = vector<size_t>(1, i);
		} else {
			mapIter->second.push_back(i);
		}
		_varSigs[i] = hashVal;
	}
	if( _verbose >= 4 ) {
		dai::operator <<(cout, _varInbox);
		cout << endl;
	}

	// determine new colors of factors
	_facClustering.clear();
	for (size_t i=0; i<_cfg.nrFactors(); i++) {
		size_t hashVal = facColor(i);
		mapIter = _facClustering.find(hashVal);
		if (mapIter == _facClustering.end()) {
			_facClustering[hashVal] = vector<size_t>(1, i);
		} else {
			mapIter->second.push_back(i);
		}
		_facSigs[i] = hashVal;
	}
	if( _verbose >= 4 ) {
		dai::operator <<(cout, _facInbox);
		cout << endl;
	}
}

bool CompressInterface::hasConverged() {
	if (_varClustering.size() == _varClusterCountOld /*&& _facClustering.size() == _facClusterCountOld*/) {
		return true;
	}
	return false;
}

size_t CompressInterface::facColor(size_t i) {
	for (size_t j=0; j< _cfg.nbF(i).size(); j++) {
		_facInbox[i][j] = _varSigs[_cfg.nbF(i,j)];
	}
	_facInbox[i][_facInbox[i].size() - 1] = _facSigs[i];
	return hashVector(_facInbox[i]);
}

size_t CompressInterface::varColor(size_t i) {
	foreach (const dai::BipartiteGraph::Neighbor tmpFac, _cfg.nbV(i)) {
		vector<size_t> hashAndPos(2);
		hashAndPos[0] = _facSigs[tmpFac];
		hashAndPos[1] = tmpFac.dual;
		_varInbox[i][tmpFac.iter] = hashVector(hashAndPos);

	}
	_varInbox[i][_varInbox[i].size() - 1] = _varSigs[i];
	sort(_varInbox[i].begin(), _varInbox[i].end());

	return hashVector(_varInbox[i]);
}

void CompressInterface::createFacClustering() {
	_facColorVec = vector<size_t>(_cfg.nrFactors());
	size_t cluster = 0;
	for (map<size_t, vector<size_t> >::const_iterator mapIter=_facClustering.begin(); mapIter!=_facClustering.end(); mapIter++ ) {
		for (size_t i=0; i<mapIter->second.size(); i++) {
			if (i==0) {
				_facRepr[cluster] = mapIter->second[i];
			}
			_facColorVec[mapIter->second[i]] = cluster;
		}
		cluster++;
	}
}

void CompressInterface::createVarClustering() {
	_varColorVec = vector<size_t>(_cfg.nrVars());
	size_t cluster = 0;
	for (map<size_t, vector<size_t> >::const_iterator mapIter=_varClustering.begin(); mapIter!=_varClustering.end(); mapIter++ ) {
		for (size_t i=0; i<mapIter->second.size(); i++) {
			if (i==0) {
				_varRepr[cluster] = _cfg.var(mapIter->second[i]);
			}
			_varColorVec[mapIter->second[i]] = cluster;
		}
		cluster++;
	}
}

CFactorGraph CompressInterface::compress(const CFactorGraph& cfg) {
	if( _verbose >= 3 ) {
		cout << "setCfg:" << endl;
	}
	setCfg(cfg);
	if( _verbose >= 3 ) {
		cout << "execute:" << endl;
	}
	execute();
	if( _verbose >= 3 ) {
		cout << "createCFactorGraph:" << endl;
	}
	return createCFactorGraph();
}

CFactorGraph CompressInterface::compress(const CFactorGraph& cfg, const vector<size_t> &varColors) {
	_initialVarColors = varColors;
	return compress(cfg);
}

CFactorGraph CompressInterface::createCFactorGraph() {

	double tic = toc();

	createVarClustering();
	if( _verbose >= 3 ) {
		cout << name() << "::createCFactorGraph()->createFacClustering took " << toc() - tic << " secs." <<  endl;
	}
	createFacClustering();
	if( _verbose >= 3 ) {
		cout << name() << "::createCFactorGraph()->createFacClustering took " << toc() - tic << " secs." <<  endl;
	}

	// create lifted fg here
	vector<CFactor> superFacs;
	superFacs.reserve(_facRepr.size());

	// clusterIdx => facIdx
	for (map<size_t,size_t>::iterator facIter = _facRepr.begin(); facIter != _facRepr.end(); facIter++) {
		VarSet superVarSet;

		foreach (const dai::BipartiteGraph::Neighbor &tmpVar, _cfg.nbF(facIter->second)) {
			if (!superVarSet.contains(_varRepr[_varColorVec[tmpVar]])) {
				superVarSet |= _varRepr[_varColorVec[tmpVar]];
			}
		}

		CFactor superFac = CFactor(superVarSet, _cfg.factor(facIter->second).p());
		superFac.sigma() = _cfg.factor(facIter->second).sigma();
		superFac.position() = _cfg.factor(facIter->second).position();

		superFac.counts() = createCounts(facIter->second, superVarSet);
		superFacs.push_back(superFac);
	}
	if( _verbose >= 3 ) {
		cout << name() << "::createCFactorGraph()->created superFacs: " << toc() - tic << " secs." <<  endl;
	}
	return CFactorGraph(superFacs);
}

std::vector<std::map<size_t, int> > CompressInterface::createCounts(size_t &gndFactor, VarSet &superVarSet) {
	// create zero entries for each position
	map<long, map<size_t, int> > countMap;
	foreach (const dai::BipartiteGraph::Neighbor &tmpVar, _cfg.nbF(gndFactor)) {
		Var liftedVar = _varRepr[_varColorVec[tmpVar]];
		size_t pos = find(_cfg.factor(gndFactor).sigma().begin(), _cfg.factor(gndFactor).sigma().end(), tmpVar.iter) - _cfg.factor(gndFactor).sigma().begin();
		countMap[liftedVar.label()][pos] = 0;
	}

	vector<map<size_t, int> > counts;
	for (vector<Var>::const_iterator iter = superVarSet.begin(); iter < superVarSet.end(); iter++) {
		foreach(const dai::BipartiteGraph::Neighbor gndFac, _cfg.nbV(_cfg.findVar(*iter))) {
			if (_facRepr[_facColorVec[gndFac]] == gndFactor) {
				size_t pos = find(_cfg.factor(gndFac).sigma().begin(), _cfg.factor(gndFac).sigma().end(), gndFac.dual) - _cfg.factor(gndFac).sigma().begin();
				countMap[iter->label()][pos]++;
			}
		}
		counts.push_back(countMap[iter->label()]);
	}
	return counts;
}

void CompressInterface::rearrangeCnfFactor(CFactor &cnfClause) {
	State state(cnfClause.vars());
	vector<size_t> (cnfClause.vars().size());
	vector<size_t> pos,neg,sigma;
	for (size_t i=0; i<cnfClause.states(); i++) {
		if (cnfClause[i] <= 1) {
			for (vector<Var>::const_iterator varIter=cnfClause.vars().begin(); varIter!=cnfClause.vars().end();varIter++) {
				if (state(*varIter) == 0) {
					pos.push_back(varIter - cnfClause.vars().begin());
				} else {
					neg.push_back(varIter - cnfClause.vars().begin());
				}
			}
			break;
		} else if (cnfClause[i] == 0) {
			//TODO handle factors which contain zero states
			// in such cases the possible zero state could have been clamped;
			// possible solutions:
			// i) restore factor from backup (requires to clamp factors accordingly)
			// ii) cache zero states in advance
			return;
		}
		state++;
	}
	sigma.insert(sigma.begin(),pos.begin(), pos.end());
	sigma.insert(sigma.begin() + pos.size(),neg.begin(), neg.end());
	cnfClause.sigma() = sigma;
}

}

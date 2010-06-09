/*
 * LiftedGibbs.cpp
 *
 *  Created on: Mar 30, 2010
 *      Author: stream
 */
#include <STREAM/Gibbs/LiftedGibbs.h>
#include <dai/util.h>
#include <stdio.h>

namespace stream {

using namespace std;
using namespace dai;

string LiftedGibbs::_name = "LiftedGibbs";

void LiftedGibbs::setProperties( const PropertySet &opts ) {
	DAI_ASSERT( opts.hasKey("iters") );
	props.iters = opts.getStringAs<size_t>("iters");

	if( opts.hasKey("burnin") )
		props.burnin = opts.getStringAs<size_t>("burnin");
	else
		props.burnin = 0;

	if( opts.hasKey("verbose") )
		props.verbose = opts.getStringAs<size_t>("verbose");
	else
		props.verbose = 0;
}

string LiftedGibbs::printProperties() const {
	stringstream s( stringstream::out );
	s << "[";
	s << "iters=" << props.iters << ",";
	s << "burnin=" << props.burnin << ",";
	s << "verbose=" << props.verbose << "]";
	return s.str();
}

void LiftedGibbs::construct() {
	_var_counts.clear();
	_var_counts.reserve( _fg.gndVarToSuperVar().size() );
	for( size_t i = 0; i < _fg.gndVarToSuperVar().size(); i++ ) {
		_var_counts.push_back( _count_t( _fg.var(_fg.reprV(i)).states(), 0 ) );
	}

	_sample_count = 0;

	_states.clear();
	_states.resize( _fg.gndVarToSuperVar().size(), 0 );
	randomizeStates();

	_schedule.reserve(_fg.gndVarToSuperVar().size());
	for (size_t i=0; i<_fg.gndVarToSuperVar().size(); i++) {
		_schedule.push_back(i);
	}

	_configMapping.reserve(_fg.nrFactors());
	vector<vector<vector<size_t> > > configs;
	for (size_t i=0; i<_fg.nrFactors(); i++) {
		configs.clear();
		configs.reserve(_fg.factor(i).states());
		for (size_t j=0; j<_fg.factor(i).states(); j++) {
			configs.push_back(vector<vector<size_t> >());
		}
		_configMapping.push_back(configs);
	}


}

Real LiftedGibbs::run() {
	if( props.verbose >= 1 )
		cout << "Starting " << identify() << "...";
	if( props.verbose >= 3 )
		cout << endl;

	double tic = toc();

	cout << "schedule: ";
	dai::operator <<(cout, _schedule) << endl;

	for( size_t iter = 0; iter < props.iters; iter++ ) {
		for (vector<size_t>::const_iterator j=_schedule.begin(); j<_schedule.end(); j++) {
			// j is an idx for a ground variable
			resampleVar( *j );
		}
		updateCounts();
	}

	cout << "iters: " << _sample_count++ << endl;
	for( size_t i = 0; i < _var_counts.size(); i++ ) {
		cout << i << " [" << _var_counts[i][0] << "," << _var_counts[i][1] << "] = " << _var_counts[i][0] + _var_counts[i][1] << endl;
	}

	if( props.verbose >= 3 ) {
		for( size_t i = 0; i < _fg.nrVars(); i++ ) {
			cout << "belief for variable " << _fg.var(i) << ": " << beliefV(i) << endl;
			cout << "counts for variable " << _fg.var(i) << ": " << Prob( _var_counts[i] ) << endl;
		}
	}

	if( props.verbose >= 1 )
		cout << _name << "::run:  ran " << props.iters << " passes (" << toc() - tic << " clocks)." << endl;

	return 0.0;
}

void LiftedGibbs::randomizeStates() {
	for( size_t i = 0; i < _states.size(); i++ ) {
		_states[i] = rnd( _fg.var(_fg.reprV(i)).states() );
	}
}

inline void LiftedGibbs::resampleVar( size_t gndVarIdx ) {
	if (_sample_count == 0) {
		cout << "sample var " << gndVarIdx << endl;
		cout << "states: ";
		dai::operator <<(cout, _states) << endl;
	}
	size_t oldState = _states[gndVarIdx];
	_states[gndVarIdx] = getVarDist(gndVarIdx).draw();
	if (oldState != _states[gndVarIdx]) {
		updateConfigs( gndVarIdx );
	}
}

void LiftedGibbs::updateConfigs(size_t gndVarIdx) {
	vector<vector<vector<size_t> > > newConfigs;
	size_t j,skip;
	size_t superVarIdx = _fg.reprV(gndVarIdx);
	vector<size_t>::const_iterator pos;
	vector<vector<size_t> >::reverse_iterator rPos;

	foreach(const BipartiteGraph::Neighbor& tmpFac, _fg.nbV(superVarIdx)) {
		if (_sample_count == 0) {
			dai::operator <<(cout, _configMapping[tmpFac]) << endl;
		}

		newConfigs.assign(_fg.factor(tmpFac).states(), vector<vector<size_t> >());
		for (j=0; j<_fg.factor(tmpFac).states(); j++) {
			for ( rPos=_configMapping[tmpFac][j].rbegin() ; rPos < _configMapping[tmpFac][j].rend(); ++rPos ) {
				pos = find(rPos->begin(), rPos->end(), gndVarIdx);
				if (pos != rPos->end()) {
					skip = getFactorEntryDiff(pos - rPos->begin());

					if (_sample_count == 0) {
						cout << gndVarIdx << " found in ";
						dai::operator <<(cout, *rPos) << endl;
					}

					if (_states[gndVarIdx] == 0) {
						newConfigs[j - skip].push_back(*rPos);
					} else {
						newConfigs[j + skip].push_back(*rPos);
					}

					_configMapping[tmpFac][j].erase(rPos.base()-1);
				}
			}
		}

		if (_sample_count == 0) {
			dai::operator << (cout, _configMapping[tmpFac]);
			cout << " + ";
			dai::operator << (cout, newConfigs) << endl;
		}

		for (j=0; j<_fg.factor(tmpFac).states(); j++) {
			if (newConfigs[j].size() > 0) {
				_configMapping[tmpFac][j].insert(_configMapping[tmpFac][j].end(),newConfigs[j].begin(),newConfigs[j].end());
			}
		}

		if (_sample_count == 0) {
			dai::operator << (cout, _configMapping[tmpFac]);
		}

	}
}

Prob LiftedGibbs::getVarDist( size_t gndVarIdx ) {
	Prob iGivenMB( 2, 1.0 );
	size_t skip, entry;
	vector<size_t>::iterator pos;
	size_t superVarIdx = _fg.reprV(gndVarIdx);
	size_t j, k, l,posNum;
	vector<size_t> posAndCount;
	bool possiblePos;
	unsigned int ix;
	if (_sample_count == 0) {
		cout << "state is " << _states[gndVarIdx] << endl;
		cout << "Neighbors of " << gndVarIdx << endl;
	}
	foreach(const BipartiteGraph::Neighbor& tmpFac, _fg.nbV(superVarIdx)) {
		if (_sample_count == 0) {
			cout << "\t" << tmpFac << ": " << _fg.factor(tmpFac) << endl << "\t";
			dai::operator <<(cout, _fg.count(superVarIdx,tmpFac.iter)) << endl;
		}
		for (j=0; j<_configMapping[tmpFac].size(); j++) {
			if (false) {
				for (k=0; k<_configMapping[tmpFac][j].size(); k++) {
					pos = find(_configMapping[tmpFac][j][k].begin(), _configMapping[tmpFac][j][k].end(), gndVarIdx);
					if (pos != _configMapping[tmpFac][j][k].end()) {
						skip = getFactorEntryDiff(pos - _configMapping[tmpFac][j][k].begin());
						entry = j - (_states[gndVarIdx] * skip);
						if (_sample_count == 0) {
							cout << "\t\t" << pos - _configMapping[tmpFac][j][k].begin() << " ";
							dai::operator <<(cout, _configMapping[tmpFac][j][k]);
							printf (" [%f|%f]", _fg.factor(tmpFac)[entry], _fg.factor(tmpFac)[entry + skip]);
							cout << endl;
						}
						for (l=0; l<2; l++) {
							iGivenMB[l] *= _fg.factor(tmpFac)[entry];
							entry += skip;
						}
					}
				}
			} else {
				possiblePos = false;
				ix = 1;
				for (k=0;k<_fg.count(superVarIdx,tmpFac.iter).size(); k++) {
					if (_fg.count(superVarIdx,tmpFac.iter)[k] > 0) {
						if (j & ix && _states[gndVarIdx] == 1 or !(j & ix) && _states[gndVarIdx] == 0 ) {
							possiblePos = true;
							break;
						}
					}
					ix <<= 1;
				}
				if (!possiblePos) {
					continue;
				}

				posAndCount.assign(_fg.count(superVarIdx,tmpFac.iter).size(), 0);
				for (k=0; k<_configMapping[tmpFac][j].size(); k++) {
					pos = find(_configMapping[tmpFac][j][k].begin(), _configMapping[tmpFac][j][k].end(), gndVarIdx);
					if (pos != _configMapping[tmpFac][j][k].end()) {
						posNum = pos - _configMapping[tmpFac][j][k].begin();
						posAndCount[posNum]++;
					}
				}

				for (size_t k=0; k<posAndCount.size(); k++) {
					if (posAndCount[k] > 0) {
						skip = getFactorEntryDiff(k);
						entry = j - (_states[gndVarIdx] * skip);
						if (posAndCount[k] == 1) {
							for (l=0; l<2; l++) {
								iGivenMB[l] *= _fg.factor(tmpFac)[entry];
								entry += skip;
							}
						} else {
							for (l=0; l<2; l++) {
								iGivenMB[l] *= pow(_fg.factor(tmpFac)[entry], posAndCount[k]);
								entry += skip;
							}
						}
					}
				}
			}
		}
	}

	if (_sample_count == 0) {
		cout << "iGivenMB: " << iGivenMB << endl;
	}
	iGivenMB.normalize();
	return iGivenMB;
}

inline size_t LiftedGibbs::getFactorEntry(const vector<size_t>& groundMB) {
}

inline size_t LiftedGibbs::getFactorEntryDiff( size_t pos ) {
//	skip = 1;
//	for (l=0; l<rPos->size(); l++) {
//		if (gndVarIdx == (*rPos)[l]) {
//			break;
//		} else {
//			skip *= 2;
//		}
//	}
	return pow(2, pos);
}

Factor LiftedGibbs::beliefV( size_t superVarIdx ) const {
	_count_t accumulatedCounts(_fg.var(superVarIdx).states(), 0);
	for (size_t i=0; i<_var_counts.size(); i++) {
		if (_fg.reprV(i) == superVarIdx) {
			for (size_t j=0; j<_var_counts[i].size(); j++) {
				accumulatedCounts[j] += _var_counts[i][j];
			}
		}
	}
	return Factor( _fg.var(superVarIdx), accumulatedCounts ).normalized();
}

void LiftedGibbs::updateCounts() {
	_sample_count++;
	if( _sample_count > props.burnin ) {
		for( size_t i = 0; i < _states.size(); i++ ) {
			_var_counts[i][_states[i]]++;
		}
	}
}

void LiftedGibbs::createConfigMapping() {
	map<Var, size_t> mapping;
	size_t linearIdx;
	vector<size_t> varIdxs;
	for (size_t i=0; i<_fg.nrFactors(); i++) {
		varIdxs.clear();
		varIdxs.reserve(_fg.nbF(i).size());
		mapping.clear();
		foreach (const BipartiteGraph::Neighbor& tmpVar, _fg.nbF(i)) {
			varIdxs.push_back(tmpVar);
			mapping[_fg.var(tmpVar)] = _states[tmpVar];
		}
		linearIdx = _fg.factor(i).vars().calcState(mapping);
		_configMapping[i][linearIdx].push_back(varIdxs);
	}
}

void LiftedGibbs::createConfigMapping(CFactorGraph &cfg) {
	map<Var, size_t> mapping;
	size_t linearIdx;
	size_t prod;
	vector<size_t> varIdxs;
	dai::operator <<(cout, _states) << endl;

	for (size_t i=0; i<cfg.nrFactors(); i++) {
		cout << cfg.factor(i) << endl;
		cout << _fg.factor(_fg.reprF(i)) << endl;
		varIdxs.assign(cfg.nbF(i).size(), 0);

		prod = 1;
		linearIdx = 0;

		for (size_t j=0; j<_fg.factor(_fg.reprF(i)).sigma().size(); j++) {
			// IMPORTANT:
			// we have to use the idxs which represent the ground variables in the super factor
			size_t posInSuper = find(_fg.factor(_fg.reprF(i)).sigma().begin(),_fg.factor(_fg.reprF(i)).sigma().end(), j) - _fg.factor(_fg.reprF(i)).sigma().begin();
			size_t posInGround = cfg.factor(i).sigma()[posInSuper];
			linearIdx += prod * _states[cfg.nbF(i, posInGround)];
			prod *= cfg.var(cfg.nbF(i, posInGround)).states();

			varIdxs[posInSuper] = cfg.nbF(i, cfg.factor(i).sigma()[j]);
		}

//		for (size_t j=0; j<cfg.getSigma(i).size(); j++) {
//
//			cout << "position for " << cfg.nbF(i,j) << "=>"<< cfg.var(cfg.nbF(i,j)) << " in ground is " << cfg.getSigma(i)[j] << endl;
//			size_t posInSuperFac = find(cfg.getSigma(_fg.reprF(i)).begin(), cfg.getSigma(_fg.reprF(i)).end(), cfg.getSigma(i)[j]) - cfg.getSigma(_fg.reprF(i)).begin();
//			cout << "and is represented by " << cfg.nbF(_fg.reprF(i), posInSuperFac) << "=>" << cfg.var(cfg.nbF(_fg.reprF(i), posInSuperFac)) << endl;
//			mapping[cfg.var(_fg.nbF(_fg.reprF(i), posInSuperFac))] = _states[cfg.nbF(i,j)];
////			mapping[cfg.var(_fg.reprV(varIdxs[j]))] = _states[varIdxs[j]];
//		}
		dai::operator <<(cout, varIdxs) << endl;
		cout << "linearIdx: " << linearIdx << endl;
		_configMapping[_fg.reprF(i)][linearIdx].push_back(varIdxs);
	}
}

}

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
 * CFactorGraph.cpp
 *
 *  Created on: Jul 13, 2009
 *      Author: Fabian Hadiji
 *
 * Code originally copied from <src/factorgraph.cpp> from libDAI
 */

#include <STREAM/CBP/CFactorGraph.h>
#include <STREAM/FactorGraph.h>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <fstream>
#include <string>
#include <algorithm>
#include <functional>
#include <dai/util.h>
#include <dai/exceptions.h>

namespace stream {

using namespace dai;
using namespace std;

CFactorGraph::CFactorGraph(const std::vector<CFactor>& factors) :
	G(), _backup(), _verbose(1) {
	// add factors, obtain variables
	set<Var> varset;
	_factors.reserve(factors.size());
	size_t nrEdges = 0;
	for (vector<CFactor>::const_iterator p2 = factors.begin(); p2 != factors.end(); p2++) {
		_factors.push_back(*p2);
		copy(p2->vars().begin(), p2->vars().end(), inserter(varset,
				varset.begin()));
		nrEdges += p2->vars().size();
	}

	// add vars
	_vars.reserve(varset.size());
	for (set<Var>::const_iterator p1 = varset.begin(); p1 != varset.end(); p1++)
		_vars.push_back(*p1);

	// create graph structure
	constructGraph(nrEdges);
}

CFactorGraph::CFactorGraph(const std::vector<dai::Factor>& factors) :
	G(), _backup(), _verbose(1) {
	// add factors, obtain variables
	set<Var> varset;
	_factors.reserve(factors.size());
	size_t nrEdges = 0;
	for (vector<Factor>::const_iterator p2 = factors.begin(); p2 != factors.end(); p2++) {
		_factors.push_back(CFactor(*p2));
		copy(p2->vars().begin(), p2->vars().end(), inserter(varset,
				varset.begin()));
		nrEdges += p2->vars().size();
	}

	// add vars
	_vars.reserve(varset.size());
	for (set<Var>::const_iterator p1 = varset.begin(); p1 != varset.end(); p1++)
		_vars.push_back(*p1);

	// create graph structure
	constructGraph(nrEdges);
}

void CFactorGraph::constructGraph(size_t nrEdges) {
	// create a mapping for indices
	hash_map<size_t, size_t> hashmap;

	for (size_t i = 0; i < vars().size(); i++)
		hashmap[var(i).label()] = i;

	// create edge list
	vector<Edge> edges;
	edges.reserve(nrEdges);
	for (size_t i2 = 0; i2 < nrFactors(); i2++) {
		const VarSet& ns = factor(i2).vars();
		for (VarSet::const_iterator q = ns.begin(); q != ns.end(); q++)
			edges.push_back(Edge(hashmap[q->label()], i2));
	}

	// create bipartite graph
	G.construct(nrVars(), nrFactors(), edges.begin(), edges.end());
}

void CFactorGraph::makeCavity(unsigned i, bool backup) {
	// fills all Factors that include var(i) with ones
	map<size_t, CFactor> newFacs;
	foreach( const Neighbor &I, nbV(i) )
	// for all neighboring factors I of i
newFacs	[I] = CFactor(factor(I).vars(), 1.0);
	setFactors( newFacs, backup );
}

void CFactorGraph::writeToFile( const char *filename ) const {
	ofstream outfile;
	outfile.open( filename );
	if( outfile.is_open() ) {
		outfile << *this;
		outfile.close();
	} else
	DAI_THROW(CANNOT_WRITE_FILE);
}

void CFactorGraph::printDot( std::ostream &os ) const {
	os << "graph G {" << endl;
	os << "node[shape=circle,width=0.4,fixedsize=true];" << endl;
	for( size_t i = 0; i < nrVars(); i++ )
	os << "\tv" << var(i).label() << ";" << endl;
	os << "node[shape=box,width=0.3,height=0.3,fixedsize=true];" << endl;
	for( size_t I = 0; I < nrFactors(); I++ )
	os << "\tf" << I << ";" << endl;
	for( size_t i = 0; i < nrVars(); i++ )
	foreach( const Neighbor &I, nbV(i) ) // for all neighboring factors I of i
	os << "\tv" << var(i).label() << " -- f" << I << ";" << endl;
	os << "}" << endl;
}

void CFactorGraph::clamp( const Var & n, size_t i, bool backup ) {
	assert( i <= n.states() );

	// store evidence of variable
	_evidence[n.label()] = i;

	// Multiply each factor that contains the variable with a delta function
	CFactor delta_n_i(n,0.0);
	delta_n_i[i] = 1.0;

	map<size_t, CFactor> newFacs;
	// For all factors that contain n
	for( size_t I = 0; I < nrFactors(); I++ ) {
		if( factor(I).vars().contains( n ) ) {
			// Multiply it with a delta function
			newFacs[I] = factor(I) * delta_n_i;
			newFacs[I].sigma() = factor(I).sigma();
			newFacs[I].counts() = factor(I).counts();
		}
	}

	setFactors( newFacs, backup );

	return;
}

// Method to clamp mln factors, makes following assumtions:
// - variables are binary
// - state with weight 1.0 is the False state
// - any other weight is seen as true
void CFactorGraph::clampMlnFactor(size_t I, bool clampState, float weight, bool backup ) {
	CFactor newFac(factor(I));

	for( size_t state = 0; state < factor(I).states(); state++){
		if ((fabs(factor(I)[state] - 1.0) < 1e-8) != clampState)
			newFac[state] = weight;
		else
			newFac[state] = factor(I)[state];
	}

	newFac.sigma() = factor(I).sigma();
	newFac.counts() = factor(I).counts();

	setFactor( I, newFac, backup );
	return;
}



void CFactorGraph::backupFactor( size_t I ) {
	map<size_t,CFactor>::iterator it = _backup.find( I );
	if( it != _backup.end() ){
		DAI_THROW( MULTIPLE_UNDO );
	}
	_backup[I] = factor(I);
}

void CFactorGraph::backupFactors( const std::set<size_t> & facs ) {
	for( std::set<size_t>::const_iterator fac = facs.begin(); fac != facs.end(); fac++ )
	backupFactor( *fac );
}

void CFactorGraph::backupFactors( const VarSet &ns ) {
	for( size_t I = 0; I < nrFactors(); I++ )
	if( factor(I).vars().intersects( ns ) )
	backupFactor( I );
}


void CFactorGraph::restoreFactor( size_t I ) {
    map<size_t,CFactor>::iterator it = _backup.find( I );
    if( it != _backup.end() ) {
        setFactor(I, it->second);
        _backup.erase(it);
    }
}


void CFactorGraph::restoreFactors( const VarSet &ns ) {
    map<size_t,CFactor> facs;
    for( map<size_t,CFactor>::iterator uI = _backup.begin(); uI != _backup.end(); ) {
        if( factor(uI->first).vars().intersects( ns ) ) {
            facs.insert( *uI );
            _backup.erase(uI++);
        } else
            uI++;
    }
    setFactors( facs );
}


void CFactorGraph::restoreFactors() {
    setFactors( _backup );
    _backup.clear();
}


void CFactorGraph::readFromFile( const char *filename ) {
	ifstream infile;
	infile.open( filename );
	if( infile.is_open() ) {
		infile >> *this;
		infile.close();
	} else {
		DAI_THROW(CANNOT_READ_FILE);
	}
}

/*
vector<CFactor> readFactorsFromFile(const char *filename) {
	ifstream infile;
	infile.open( filename );
	if( infile.is_open() ) {
		infile >> *this;
		infile.close();
	} else {
		DAI_THROW(CANNOT_READ_FILE);
	}
}
*/

/// Reads a FactorGraph from an input stream
istream& operator >> (istream& is, CFactorGraph& cfg) {
	vector<CFactor> facs;
	CFactorGraph::readFactors(&facs, is, cfg._verbose);
	cfg = CFactorGraph(facs);
	return is;
}

void CFactorGraph::readFactors (vector<CFactor>* facs, istream& is, long verbose) {
	try {

		size_t nr_Factors;
		string line;

		while( (is.peek()) == '#' ) {
			getline(is,line);
		}
		is >> nr_Factors;
		if( is.fail() ) {
			DAI_THROW(INVALID_FACTORGRAPH_FILE);
		}
		if( verbose >= 2 ) {
			cout << "Reading " << nr_Factors << " factors..." << endl;
		}

		// empty line after number of factors
		getline (is,line);
		if( is.fail() ) {
			DAI_THROW(INVALID_FACTORGRAPH_FILE);
		}
		map<long,size_t> vardims;
		for( size_t I = 0; I < nr_Factors; I++ ) {
			if( verbose >= 3 ) {
				cout << "Reading factor " << I << "..." << endl;
			}

			size_t nrSuperVars = 0;
			// empty lines between factors
			while (is.peek() == '\n') {
				getline(is,line);
			}

			while( (is.peek()) == '#' )  {
				// get #
				is.get();
				// number of super vars is given
				if (is.peek() == '!') {
					// get !
					is.get();
					is >> nrSuperVars;
					if( verbose >= 3 ) {
						cout << "  nrSuperVars: " << nrSuperVars << endl;
					}
				}
				getline(is,line);
			}

			size_t nr_members;
			while( (is.peek()) == '#' )  {
				getline(is,line);
			}
			is >> nr_members;
			if( verbose >= 3 ) {
				cout << "  nr_members: " << nr_members << endl;
			}
			vector<long> labels;
			for( size_t mi = 0; mi < nr_members; mi++ ) {
				long mi_label;
				while( (is.peek()) == '#' ) {
					getline(is,line);
				}
				is >> mi_label;
				labels.push_back(mi_label);
			}
			if( verbose >= 3 ) {
				cout << "  labels: ";
				dai::operator << (cout, labels);
				cout << endl;
			}

			// counts
			vector<map<size_t, int> > counts;
			map<size_t, int> count;
			// the following caputres whitespace at the end of the line and the line feed
			getline(is,line);

			if ( (is.peek()) == '#') {
				// get #
				is.get();
				// counts are upcoming
				if (is.peek() == '!') {
					// get !
					is.get();
					char nextChar;
					size_t nextPos;
					int nextCount;
					if (is.peek() == '[') {
						for( size_t mi = 0; mi < nr_members; mi++ ) {
							count.clear();
							// next is [
							is >> nextChar;
							while (nextChar != ']') {
								is >> nextPos;
								// get : between pos and count
								is.get();
								is >> nextCount;
								count[nextPos] = nextCount;
								is >> nextChar;
							}
							counts.push_back(count);

							while (is.peek() == ' ') {
								is.get();
							}
						}
						getline(is,line);
					}
					if( verbose >= 3 ) {
						cout << "  counts: ";
						dai::operator << (cout, counts);
						cout << endl;
					}
				}
			}

			// read dimensions of variables
			vector<size_t> dims;
			size_t nrDims = max(nr_members, nrSuperVars);
			size_t nrStates = 1;
			for( size_t mi = 0; mi < nrDims; mi++ ) {
				size_t mi_dim;
				while( (is.peek()) == '#' ) {
					getline(is,line);
				}
				is >> mi_dim;
				nrStates *= mi_dim;
				dims.push_back(mi_dim);
			}
			if( verbose >= 3 ) {
				cout << "  dimensions: ";
				dai::operator <<(cout, dims);
				cout << endl;
			}

			// add the Factor
			VarSet I_vars;
			for( size_t mi = 0; mi < nr_members; mi++ ) {
				map<long,size_t>::iterator vdi = vardims.find( labels[mi] );
				if( vdi != vardims.end() ) {
					// check whether dimensions are consistent
					if( vdi->second != dims[mi] )
					DAI_THROW(INVALID_FACTORGRAPH_FILE);
				} else {
					vardims[labels[mi]] = dims[mi];
				}
				I_vars |= Var(labels[mi], dims[mi]);
			}
			facs->push_back( CFactor( I_vars, 0.0 ) );
			// add counts, in case no these were read before
			if (counts.size() > 0 ) {
				// resize factor table if necessary
				facs->back().resize(nrStates);
				// rearrange counts according to order of label
				for(VarSet::const_iterator j = I_vars.begin(); j != I_vars.end(); j++ ) {
					long search_for = j->label();
					vector<long>::iterator j_loc = find(labels.begin(),labels.end(),search_for);
					facs->back().counts()[j - I_vars.begin()] = counts[j_loc - labels.begin()];
				}
			} else {
				for(VarSet::const_iterator j = I_vars.begin(); j != I_vars.end(); j++ ) {
					count.clear();
					long search_for = j->label();
					vector<long>::iterator j_loc = find(labels.begin(),labels.end(),search_for);
					count[j_loc - labels.begin()] = 1;
					facs->back().counts()[j - I_vars.begin()] = count;
				}
			}

			if( verbose >= 3 ) {
				cout << "  counts: ";
				dai::operator <<(cout, facs->back().counts()) << endl;
			}

			// calculate permutation sigma (internally, members are sorted)
			vector<size_t> sigma(dims.size(),0);
			VarSet::iterator j = I_vars.begin();
			size_t varCount = 0;
			for( ; j != I_vars.end(); j++ ) {
				long search_for = j->label();
				vector<long>::iterator j_loc = find(labels.begin(),labels.end(),search_for);
				// sigmas
				if (nrSuperVars == 0 || nrSuperVars == nr_members) {
					sigma[varCount] = j_loc - labels.begin();
					varCount++;
				} else {
					for (map<size_t, int>::const_iterator pos = facs->back().counts()[j - I_vars.begin()].begin(); pos != facs->back().counts()[j - I_vars.begin()].end(); pos++) {
						sigma[varCount] = pos->first;
						varCount++;
					}
				}
			}
			if( verbose >= 3 ) {
				cout << "  sigma: ";
				dai::operator <<(cout, sigma);
				cout << endl;
			}

			// store sigma
			vector<size_t> undoSigma(sigma.size(),0);
			for (size_t mi=0; mi<sigma.size(); mi++) {
				undoSigma[sigma[mi]] = mi;
			}
			facs->back().sigma() = undoSigma;
			if( verbose >= 3 ) {
				cout << "  undoSigma: ";
				dai::operator <<(cout, facs->back().sigma());
				cout << endl;
			}

			// position
			// TODO should this be relative to size of supervars?
			map<size_t, size_t> position;
			for (size_t i=0; i<labels.size(); i++) {
				position[labels[i]] = i;
			}
			facs->back().position() = position;
			if( verbose >= 3 ) {
				cout << "  position: ";
				dai::operator <<(cout, facs->back().position());
				cout << endl;
			}

			// calculate multindices
			Permute permindex( dims, sigma );

			// read values
			size_t nr_nonzeros;
			while( (is.peek()) == '#' ) {
				getline(is,line);
			}
			is >> nr_nonzeros;

			if( verbose >= 3 ) {
				cout << "  nonzeroes: " << nr_nonzeros << endl;
			}
			for( size_t k = 0; k < nr_nonzeros; k++ ) {
				size_t li;
				double val;
				while( (is.peek()) == '#' ) {
					getline(is,line);
				}
				is >> li;
				while( (is.peek()) == '#' ) {
					getline(is,line);
				}
				is >> val;

				// store value, but permute indices first according
				// to internal representation
				facs->back()[permindex.convert_linear_index( li )] = val;
			}
		}

		if( verbose >= 3 ) {
			cout << "factors:";
			dai::operator <<(cout, *facs);
			cout << endl;
		}


	} catch (char *e) {
		cout << e << endl;
	}
}

/// Writes a FactorGraph to an output stream
// nrVars
// counts
// dimensions
// #nonZeros
// nonZeroEntry1
// ...
// nonZeroEntryN
ostream& operator << (ostream& os, const CFactorGraph& fg) {
	os << fg.nrFactors() << endl;

	for( size_t I = 0; I < fg.nrFactors(); I++ ) {
		os << endl;
		// # of vars in factor
		size_t nrGndVars = 0;
		for (size_t i=0; i<fg.factor(I).counts().size(); i++) {
			nrGndVars += fg.factor(I).counts()[i].size();
		}
		os << "#!" << nrGndVars << endl;

		// # of super vars in factor
		os << fg.factor(I).vars().size() << endl;

		// label of variables
		for (size_t i=0; i<fg.factor(I).vars().size(); i++) {
			if (i > 0) {
				os << " ";
			}
			size_t varIdx;
			if (nrGndVars == fg.factor(I).vars().size()) {
				// the variables can be ordered according to sigma
				varIdx = fg.nbF(I,fg.factor(I).sigma()[i]);
			} else {
				varIdx = fg.nbF(I,i);
			}
			os << fg.var(varIdx).label();
		}
		os << endl;


		// counts
		os << "#!";
		for (size_t i=0; i<fg.factor(I).vars().size(); i++) {
			if (i > 0) {
				os << " ";
			}

			os << "[";
			size_t varPos;
			if (nrGndVars == fg.factor(I).vars().size()) {
				// the variables can be ordered according to sigma
				varPos = fg.factor(I).sigma()[i];
			} else {
				varPos = i;
			}
			for (map<size_t, int>::const_iterator iter=fg.factor(I).counts()[varPos].begin(); iter!=fg.factor(I).counts()[varPos].end(); iter++) {
				if (iter != fg.factor(I).counts()[varPos].begin())
					os  << ",";
				os << (*iter).first << ":" << (*iter).second;
			}
			os << "]";
		}
		os << endl;

		// dimensions (for gnd variables)
		vector<size_t> dims;
		for (size_t i=0; i<fg.factor(I).vars().size(); i++) {
			if (i > 0) {
				os << " ";
			}

			size_t varPos;
			size_t varIdx;
			if (nrGndVars == fg.factor(I).vars().size()) {
				// the variables can be ordered according to sigma
				varPos = fg.factor(I).sigma()[i];
				varIdx = fg.nbF(I,fg.factor(I).sigma()[i]);
			} else {
				varPos = i;
				varIdx = fg.nbF(I,i);
			}

			for (size_t j=0; j<fg.factor(I).counts()[varPos].size(); j++) {
				if (j > 0) {
					os << " ";
				}
				os << fg.var(varIdx).states();
				dims.push_back(fg.var(fg.nbF(I,i)).states());
			}
		}
		os << endl;

		size_t nr_nonzeros = 0;
		for( size_t k = 0; k < fg.factor(I).states(); k++ ){
			if( fg.factor(I)[k] != 0.0 ) {
				nr_nonzeros++;
			}
		}
		os << nr_nonzeros << endl;

		vector<double> tmpVec(fg.factor(I).states());

		Permute perm;
		if (fg.factor(I).sigma().size() != dims.size()) {
			//TODO sigma should be read from object instead => sigma needs to be created when compressed
			vector<size_t> gndSigma(dims.size(),0);

			size_t gndVarPos = 0;
			for (size_t varPos=0; varPos<fg.factor(I).vars().size(); varPos++ ) {
				for (map<size_t, int>::const_iterator pos = fg.factor(I).counts()[varPos].begin(); pos != fg.factor(I).counts()[varPos].end(); pos++) {
					gndSigma[pos->first] = gndVarPos;
					gndVarPos++;
				}
			}

			perm = Permute(dims, gndSigma);
		} else {
			perm = Permute(dims, fg.factor(I).sigma());
		}

		for (size_t k=0; k<fg.factor(I).states(); k++) {
			tmpVec[perm.convert_linear_index(k)] = fg.factor(I)[k];
		}

		for (size_t k=0; k<fg.factor(I).states(); k++) {
			if (tmpVec[k] != 0.0) {
				char buf[20];
				sprintf(buf,"%18.14g", tmpVec[k]);
				os << k << " " << buf << endl;
			}
		}
	}

	return(os);
}

CFactorGraph CFactorGraph::createCleanedGraph(CFactorGraph& cfg, map<long, int>& fixedVariables) throw (const char*) {
	vector<CFactor> newFacs;
		map<long, int>::const_iterator fixedVarIter;
		for (size_t i=0; i<cfg.nrFactors(); i++) {
			bool remove = false;
			VarSet newVarSet(cfg.factor(i).vars());
			foreach (const dai::BipartiteGraph::Neighbor tmpVar, cfg.nbF(i)) {
				fixedVarIter = fixedVariables.find(cfg.var(tmpVar).label());
				if (fixedVarIter != fixedVariables.end()) {
					// factor contains the variable
					// determine zero state
					size_t zeroState = 0;
					for (size_t j=0; j<cfg.factor(i).states();j++) {
						if (cfg.factor(i)[j] <= 1) {
							zeroState = j;
							break;
						}
					}

					// check sign of variable in factor
					map< Var, size_t > states = cfg.factor(i).vars().calcStates(zeroState);
					if (states.find(cfg.var(tmpVar))->second != (size_t)fixedVarIter->second) {
						// variable has same sign in factor
						remove = true;
						break;
					} else {
						// variable has opposite sing in factor
						newVarSet /= cfg.var(tmpVar);
					}
				}
			}

			if (!remove) {
				if (newVarSet.size() == cfg.factor(i).vars().size()) {
					newFacs.push_back(cfg.factor(i));
				} else if (newVarSet.size() > 0) {
					// create new factor with reduced number of variables
					CFactor newFac(newVarSet, stream::FactorGraph<CFactor, Var>::NON_ZERO_VAL);
					State state(cfg.factor(i).vars());
					map< Var, size_t > newStates;
					for(size_t j=0; j<cfg.factor(i).states(); j++) {
						if (cfg.factor(i)[j] <= 1) {
							for (vector<Var>::const_iterator varIter=cfg.factor(i).vars().begin();varIter!=cfg.factor(i).vars().end();varIter++) {
								if (newVarSet.contains(*varIter)) {
									newStates[*varIter] = state(*varIter);
								}
							}
							break;
						}
						state++;
					}
					size_t zeroState = newVarSet.calcState(newStates);
					newFac[zeroState] = stream::FactorGraph<CFactor, Var>::ZERO_VAL;
					newFacs.push_back(newFac);
				} else {
					throw "ERROR: newVarSet doesn't contain any variables";
				}
			}
		}

		return CFactorGraph(newFacs);
}

} // end of namespace

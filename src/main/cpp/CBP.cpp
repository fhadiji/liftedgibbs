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
 * CBP.cpp
 *
 *  Created on: Jul 13, 2009
 *      Author: Fabian Hadiji
 *
 * Code originally copied from <src/bp.cpp> from libDAI
 */


#include <iostream>
#include <sstream>
#include <map>
#include <set>
#include <algorithm>
#include <STREAM/CBP/CBP.h>
#include <dai/util.h>
#include <dai/properties.h>

namespace stream {

using namespace std;
using namespace dai;

const char *CBP::Name = "CBP";


void CBP::setProperties( const PropertySet &opts ) {
    assert( opts.hasKey("tol") );
    assert( opts.hasKey("maxiter") );
    assert( opts.hasKey("logdomain") );
    assert( opts.hasKey("updates") );

    props.tol = opts.getStringAs<double>("tol");
    props.maxiter = opts.getStringAs<size_t>("maxiter");
    props.logdomain = opts.getStringAs<bool>("logdomain");
    props.updates = opts.getStringAs<Properties::UpdateType>("updates");

    if( opts.hasKey("verbose") )
        props.verbose = opts.getStringAs<size_t>("verbose");
    else
        props.verbose = 0;
    if( opts.hasKey("damping") )
        props.damping = opts.getStringAs<double>("damping");
    else
        props.damping = 0.0;
    if( opts.hasKey("inference") )
        props.inference = opts.getStringAs<Properties::InfType>("inference");
    else
        props.inference = Properties::InfType::SUMPROD;
    if( opts.hasKey("compress") ) {
    	props.compress = opts.getStringAs<bool>("compress");
    } else {
    	props.compress = true;
    }
    if( opts.hasKey("compressionAlg")) {
    	props.compressionAlg = opts.getStringAs<string>("compressionAlg");
    } else {
    	props.compressionAlg = string("Simple");
    }
}

void CBP::setProperty( const string &key, const string &value ){
    setProperties(getProperties().Set( key, value ));
}

PropertySet CBP::getProperties() const {
    PropertySet opts;
    opts.Set( "tol", props.tol );
    opts.Set( "maxiter", props.maxiter );
    opts.Set( "verbose", props.verbose );
    opts.Set( "logdomain", props.logdomain );
    opts.Set( "updates", props.updates );
    opts.Set( "damping", props.damping );
    opts.Set( "inference", props.inference );
    return opts;
}


string CBP::printProperties() const {
    stringstream s( stringstream::out );
    s << "[";
    s << "tol=" << props.tol << ",";
    s << "maxiter=" << props.maxiter << ",";
    s << "verbose=" << props.verbose << ",";
    s << "logdomain=" << props.logdomain << ",";
    s << "updates=" << props.updates << ",";
    s << "damping=" << props.damping << ",";
    s << "inference=" << props.inference << ",";
    s << "compress=" << props.compress << ",";
    s << "compressionAlg=" << props.compressionAlg << "]";
    return s.str();
}

void CBP::construct() {

	// prepare datastructures for compression
	for (size_t i=0; i<nrVars(); i++) {
		_gndVarToSuperVar[i] = i;
	}

	for (size_t i=0; i<nrFactors(); i++) {
		_gndFacToSuperFac[i] = i;
	}

	// create edge properties
	_edges.clear();
	_edges.reserve( nrVars() );
	for( size_t i = 0; i < nrVars(); ++i ) {
		_edges.push_back( vector<EdgeProp>() );
		_edges[i].reserve( nbV(i).size() );
		foreach( const Neighbor &I, nbV(i) ) {

			EdgeProp newEP;
			size_t edgeCount = factor(I).counts()[I.dual].size();
			newEP.message = vector<Prob>(edgeCount,Prob( var(i).states() ));
			newEP.newMessage = vector<Prob>(edgeCount,Prob( var(i).states() ));
			newEP.count = vector<int>(edgeCount);
			newEP.index = vector<ind_t>(edgeCount);
			newEP.nrPos = edgeCount;

			// simulate orginal varSet with possibly more variables, must ensure that the number of variables is equal to number of ground variables
			VarSet gndVarSet;
			size_t varCount = 0;
			foreach(const Neighbor &tmpVar, nbF(I)) {
				for(map<size_t, int>::const_iterator iter=factor(I).counts()[tmpVar.iter].begin(); iter!=factor(I).counts()[tmpVar.iter].end();iter++) {
					gndVarSet |= Var(varCount, var(tmpVar).states());
					varCount++;
				}
			}

			varCount = 0;
			foreach(const Neighbor &tmpVar, nbF(I)) {
				size_t pos=0;
				for(map<size_t, int>::const_iterator iter=factor(I).counts()[tmpVar.iter].begin(); iter!=factor(I).counts()[tmpVar.iter].end();iter++,pos++) {
					if (tmpVar == i) {
						// assumes that counts are iterated in increases order of positions
						size_t sortedPos = factor(I).sigma()[(*iter).first];
						newEP.count[pos] = (*iter).second;
						newEP.index[pos].reserve( factor(I).states() );
						for( IndexFor k( Var(sortedPos, var(i).states()), gndVarSet ); k >= 0; ++k ) {
							newEP.index[pos].push_back( k );
						}
					}
					varCount++;
				}
			}

			_edges[i].push_back( newEP );
		}
	}
}


void CBP::init() {
	if (props.compress) {
		if( props.verbose >= 3)
			cout << "Compression..." << endl;
		compress();
	}

	double c = props.logdomain ? 0.0 : 1.0;
	for( size_t i = 0; i < nrVars(); ++i ) {
		foreach( const Neighbor &I, nbV(i) ) {
			size_t edgeCount = factor(I).counts()[I.dual].size();
			for (size_t j=0; j<edgeCount; j++) {
				message( i, I.iter, j ).fill( c );
				newMessage( i, I.iter, j ).fill( c );
			}
		}
	}
}

void CBP::init(vector<size_t> &varColors) {
	_initialVarColors = varColors;
	init();
}


// message from Factor _I to Variable i
void CBP::calcNewMessage( size_t i, size_t _I, size_t &pos ) {
	// calculate updated message I->i
	size_t I = nbV(i,_I);
	_messages++;
	if (props.verbose >= 4) {
		cout << "m_f" << I << "->v" << i << var(i) << "(" << pos << ")" << endl;
	}

	/* OPTIMIZED VERSION */
	Prob prod( factor(I).p() );
	if( props.logdomain )
		prod.takeLog();

	// Calculate product of incoming messages and factor I
	size_t pos1, pos2;
	foreach( const Neighbor &j, nbF(I) ) {
		for (pos1=0; pos1 < nrPos(j, j.dual); pos1++) {
			if(!(j == i  && pos1 == pos)) {     // for all j in I \ i

				size_t newPos1 = pos1;
				// calculate new pos is case of AnyPosition
				if (count(j,j.dual,pos1) == 0) {
					for (size_t k=pos1 - 1 ; k>=0; k--) {
						if (count(j,j.dual,k) > 0) {
							//TODO this is not 100% correct!
							// one needs to find the correct and corresponding column, not any column with count > 0
							// FOR CNFS: use order in cnfs that positive vars come first, followed by negative vars
							newPos1 = k;
							break;
						}
					}
				}

				if (props.verbose >= 4) {
					cout << "\tm_v" << j << var(j) << "->f" << I << "(" << newPos1 << ") = " << endl;
				}

				// prod_j will be the product of messages coming into j
				Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0 );
				foreach( const Neighbor &J, nbV(j) ) {
					for (pos2 = 0; pos2<nrPos(j,J.iter); pos2++) {
						if (!(J == I  && pos2 == newPos1)) {
							if (props.verbose >= 4) {
								cout << "\t\tm_f" << J << "->v" << j << var(j) << "(" << pos2 << ") ^ "<< count(j,J.iter,pos2) << " = " << message( j, J.iter, pos2) << endl;
							}
							if (count(j,J.iter,pos2) == 1) {
								if( props.logdomain ) {
										prod_j += message( j, J.iter, pos2  );
								} else {
										prod_j *= message( j, J.iter, pos2  );
								}
							} else if (count(j,J.iter,pos2) > 1) {
								if( props.logdomain ) {
										prod_j += message( j, J.iter, pos2  ) * count(j,J.iter,pos2);
								} else {
										prod_j *= message( j, J.iter, pos2  ) ^ count(j,J.iter,pos2);
								}
							}

						} else {
							if (props.verbose >= 4) {
								cout << "\t\tm_f" << J << "->v" << j << var(j) << "(" << pos2 << ") ^ "<< count(j,J.iter,pos2) - 1 << " = " << message(j, J.iter, pos2 ) << endl;
							}
							if (count(j,J.iter,pos2) == 2) {
								if( props.logdomain ) {
									prod_j += message(j, J.iter, pos2);
								} else {
									prod_j *= message(j, J.iter, pos2);
								}
							} else if (count(j,J.iter,pos2) > 2) {
								if( props.logdomain ) {
									prod_j += message(j, J.iter, pos2) * (count(j,J.iter,pos2) - 1);
								} else {
									prod_j *= message(j, J.iter, pos2) ^ (count(j,J.iter,pos2) - 1);
								}
							}
						}
					}

				}
				if (props.verbose >= 4) {
					cout << "\t= " << prod_j << endl;
				}

				size_t _I = j.dual;
				// ind is the precalculated IndexFor(j,I) i.e. to x_I == k corresponds x_j == ind[k]
				const ind_t &ind = index(j, _I, pos1);

				if (props.verbose >= 4) {
					cout << "\t";
					dai::operator <<(cout, ind) << endl;
				}

				// multiply prod with prod_j
				for( size_t r = 0; r < prod.size(); ++r )
					if( props.logdomain )
						prod[r] += prod_j[ind[r]];
					else
						prod[r] *= prod_j[ind[r]];


				if (props.verbose >= 4) {
					cout << "\t= " << prod << endl;
				}

			}
		}
	}

	if( props.logdomain ) {
		prod -= prod.maxVal();
		prod.takeExp();
	}

	// Marginalize onto i
	Prob marg( var(i).states(), 0.0 );
	// ind is the precalculated IndexFor(i,I) i.e. to x_I == k corresponds x_i == ind[k]
	const ind_t ind = index(i,_I, pos);

	if (props.verbose >= 4) {
		dai::operator <<(cout, ind);
		cout << endl;
	}

	if( props.inference == Properties::InfType::SUMPROD ) {
		for( size_t r = 0; r < prod.size(); ++r ) {
			marg[ind[r]] += prod[r];
		}
	} else
		for( size_t r = 0; r < prod.size(); ++r )
			if( prod[r] > marg[ind[r]] )
				marg[ind[r]] = prod[r];
	marg.normalize();

	// Store result
	if( props.logdomain )
		newMessage(i,_I, pos) = marg.log();
	else
		newMessage(i,_I, pos) = marg;

	if (props.verbose >= 4) {
		cout << "= " << marg << endl;
	}
}


// CBP::run does not check for NANs for performance reasons
// Somehow NaNs do not often occur in BP...
double CBP::run() {
	if( props.verbose >= 1 )
        cout << "Starting " << identify() << "...";
    if( props.verbose >= 3)
       cout << endl;

    double tic = toc();

    Diffs diffs(nrVars(), 1.0);

    vector<Factor> old_beliefs;
    old_beliefs.reserve( nrVars() );
    for( size_t i = 0; i < nrVars(); ++i )
        old_beliefs.push_back( beliefV(i) );

    _messages = 0;
    size_t  pos;
    for( _iters=0; _iters < props.maxiter && diffs.maxDiff() > props.tol; ++_iters ) {
        if( props.updates == Properties::UpdateType::PARALL ) {
            // Parallel updates
            for( size_t i = 0; i < nrVars(); ++i ) {
                foreach( const Neighbor &I, nbV(i) ) {
                	for (pos=0; pos<nrPos(i, I.iter); pos++) {
                		if (count(i, I.iter,pos) > 0) {
                			calcNewMessage( i, I.iter, pos );
                		}
                	}
                }
            }

            for( size_t i = 0; i < nrVars(); ++i ) {
                foreach( const Neighbor &I, nbV(i) ) {
                	for (pos=0; pos<nrPos(i, I.iter); pos++) {
               			updateMessage( i, I.iter, pos );
                	}
                }
            }
        }

        // calculate new beliefs and compare with old ones
        for( size_t i = 0; i < nrVars(); ++i ) {
            Factor nb( beliefV(i) );
            diffs.push( dist( nb, old_beliefs[i], Prob::DISTLINF ) );
            old_beliefs[i] = nb;
        }

        if( props.verbose >= 3 )
            cout << Name << "::run:  maxdiff " << diffs.maxDiff() << " after " << _iters+1 << " passes (" << toc() - tic << " secs.) " << endl;
    }

    if( diffs.maxDiff() > _maxdiff )
        _maxdiff = diffs.maxDiff();

    if( props.verbose >= 1 ) {
        if( diffs.maxDiff() > props.tol ) {
            if( props.verbose == 1 )
                cout << endl;
                cout << Name << "::run:  WARNING: not converged within " << props.maxiter << " passes (" << toc() - tic << " seconds)...final maxdiff:" << diffs.maxDiff() << endl;
        } else {
            if( props.verbose >= 3 )
                cout << Name << "::run:  ";
                cout << "converged in " << _iters << " passes (" << toc() - tic << " seconds)." << endl;
        }
    }

    return diffs.maxDiff();
}


Factor CBP::beliefV( size_t i ) const {
    Prob prod( var(i).states(), props.logdomain ? 0.0 : 1.0 );
    foreach( const Neighbor &I, nbV(i) ) {
    	for (size_t pos=0; pos<nrPos(i, I.iter); pos++) {
    		if (count(i, I.iter, pos) == 1) {
    			if( props.logdomain ) {
    				prod += newMessage( i, I.iter, pos );
    			} else {
    				prod *= newMessage( i, I.iter, pos );
    			}
    		} else if (count(i, I.iter, pos) > 1) {
    			if( props.logdomain ) {
    				prod += newMessage( i, I.iter, pos ) * count(i, I.iter, pos);
    			} else {
    				prod *= newMessage( i, I.iter, pos ) ^ count(i, I.iter, pos);
    			}
    		}
    	}
    }

    if( props.logdomain ) {
        prod -= prod.maxVal();
        prod.takeExp();
    }

    prod.normalize();
    return( Factor( var(i), prod ) );
}

Factor CBP::beliefF( size_t I ) const {
	Prob p;
	calcBeliefF( I, p );

	if( props.logdomain ) {
		p -= p.max();
		p.takeExp();
	}
	p.normalize();

	VarSet gndVarSet;
	size_t nrLabel = 0;
	foreach( const Neighbor &tmpVar, nbF(I) ) {
		for (size_t tmpPos=0; tmpPos < nrPos(tmpVar, tmpVar.dual); tmpPos++) {
			gndVarSet |= Var(nrLabel, var(tmpVar).states());
			nrLabel++;
		}
	}
	return( Factor( gndVarSet, p ) );
}


void CBP::calcBeliefF( size_t I, Prob &p ) const {
	CFactor Fprod( factor( I ) );
	Prob &prod = Fprod.p();
	if( props.logdomain )
		prod.takeLog();

	size_t pos1, pos2;
	foreach( const Neighbor &j, nbF(I) ) {
		for (pos1=0; pos1 < nrPos(j, j.dual); pos1++) {

			size_t newPos1 = pos1;
			// calculate new pos is case of AnyPosition
			if (count(j,j.dual,pos1) == 0) {
				for (size_t k=pos1 - 1 ; k>=0; k--) {
					if (count(j,j.dual,k) > 0) {
						newPos1 = k;
						break;
					}
				}
			}

			Prob prod_j( var(j).states(), props.logdomain ? 0.0 : 1.0 );
			foreach( const Neighbor &J, nbV(j) ) {
				for (pos2 = 0; pos2<nrPos(j,J.iter); pos2++) {
					if (!(J == I  && pos2 == newPos1)) {
						if( props.logdomain )
							prod_j += newMessage( j, J.iter, pos2) * count(j,J.iter,pos2);
						else
							prod_j *= newMessage( j, J.iter, pos2) ^ count(j,J.iter,pos2);
					} else {
						if( props.logdomain )
							prod_j += newMessage( j, J.iter, pos2) * (count(j,J.iter,pos2) - 1);
						else
							prod_j *= newMessage( j, J.iter, pos2) ^ (count(j,J.iter,pos2) - 1);
					}
				}
			}
			size_t _I = j.dual;
			const ind_t & ind = index(j, _I, pos1);

			for( size_t r = 0; r < prod.size(); ++r ) {
				if( props.logdomain )
					prod[r] += prod_j[ind[r]];
				else
					prod[r] *= prod_j[ind[r]];
			}
		}
	}

	p = prod;
}

Factor CBP::belief (const Var &n) const {
    return( beliefV( findVar( n ) ) );
}



string CBP::identify() const {
    return string(Name) + printProperties();
}


void CBP::init( const VarSet &ns ) {
    for( VarSet::const_iterator n = ns.begin(); n != ns.end(); ++n ) {
        size_t ni = findVar( *n );
        foreach( const Neighbor &tmpFac, nbV( ni ) ) {
        	size_t edgeCount = factor(tmpFac).counts()[tmpFac.dual].size();
        	for (size_t pos=0; pos<edgeCount; pos++) {
        		message( ni, tmpFac.iter, pos).fill( props.logdomain ? 0.0 : 1.0 );
        	}
        }
    }
}

void CBP::compress() {

	double tic = toc();

	setCompressionAlg(props.compressionAlg);
	_compressAlg->setVerbose(props.verbose);
	CFactorGraph liftedFg;
	if (_initialVarColors.size() > 0) {
		liftedFg = _compressAlg->compress((CFactorGraph &)(*this), _initialVarColors);
	} else {
		liftedFg = _compressAlg->compress((CFactorGraph &)(*this));
	}

	_compressIters = _compressAlg->getIterations();

	if( props.verbose >= 3) {
		cout << "supernodes: " << liftedFg.nrVars() << " (" << nrVars() << ")"<< endl;
		cout << "superfactor: " << liftedFg.nrFactors() << " (" << nrFactors() << ")" << endl;
		cout << "time to compress " << toc() - tic << " seconds" << endl;
	}

	CFactorGraph::operator=( liftedFg );
	construct();
	createVarMapping(_compressAlg, liftedFg);
	createFacMapping(_compressAlg, liftedFg);

}

void CBP::createVarMapping(CompressInterface * compressAlg, CFactorGraph& liftedFg) {
	vector<size_t> varColorVec = compressAlg->getVarColorVec();
	map<size_t, Var> varRepr = compressAlg->getVarRepr();

	for (size_t i=0; i<varColorVec.size();i++) {
		_gndVarToSuperVar[i] = liftedFg.findVar(varRepr.find(varColorVec[i])->second);
	}
}

void CBP::createFacMapping(CompressInterface * compressAlg, CFactorGraph& liftedFg) {
	map<size_t, vector<size_t> > facClustering = compressAlg->getFacClustering();
	size_t cluster = 0;
	for (map<size_t, vector<size_t> >::const_iterator mapIter=facClustering.begin(); mapIter!=facClustering.end(); mapIter++) {
		for (size_t i=0; i<mapIter->second.size(); i++) {
			_gndFacToSuperFac[mapIter->second[i]] = cluster;
		}
		cluster++;
	}
}

vector<size_t> CBP::clusterV(size_t superVar) const {
	vector<size_t> cluster;
	for (map<size_t, size_t>::const_iterator iter=_gndVarToSuperVar.begin(); iter != _gndVarToSuperVar.end(); iter++) {
		if (iter->second == superVar) {
			cluster.push_back(iter->first);
		}
	}
	return cluster;
}

vector<size_t> CBP::clusterF(size_t superFac) const {
	vector<size_t> cluster;
	for (map<size_t, size_t>::const_iterator iter=_gndFacToSuperFac.begin(); iter != _gndFacToSuperFac.end(); iter++) {
		if (iter->second == superFac) {
			cluster.push_back(iter->first);
		}
	}
	return cluster;
}

} // end of namespace STREAM

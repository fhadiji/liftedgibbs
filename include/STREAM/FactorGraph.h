/*
 * FactorGraph.h
 *
 *  Created on: Jan 6, 2010
 *      Author: Fabian Hadiji
 */

#ifndef FACTORGRAPH_H_
#define FACTORGRAPH_H_

#include <dai/factorgraph.h>
#include <STREAM/CBP/CBP.h>
#include <bitset>

namespace stream {

struct Signature {
	int pos;
	int neg;
};

template<class F, class V>
class FactorGraph {

public:

	// only for CNF conversion
	static double ZERO_VAL;
	static double NON_ZERO_VAL;

	/// Stores the neighborhood structure
	dai::BipartiteGraph                    G;

	/// Shorthand for BipartiteGraph::Neighbor
	typedef dai::BipartiteGraph::Neighbor  Neighbor;

	/// Shorthand for BipartiteGraph::Neighbors
	typedef dai::BipartiteGraph::Neighbors Neighbors;

	/// Shorthand for BipartiteGraph::Edge
	typedef dai::BipartiteGraph::Edge      Edge;

	/// Iterator over factors
	typedef typename std::vector<F>::iterator iterator;

	/// Constant iterator over factors
	typedef typename std::vector<V>::const_iterator const_iterator;


private:
	/// Stores the variables
	std::vector<V>		_vars;
	/// Stores the factors
	std::vector<F>		_factors;

	size_t _verbose;

	/// only for lifted case
	// represents the counts originally introduced in CBP paper
	// _counts is the number of indenctical factors that send equal messages to a variable
	std::vector<std::vector<std::vector<int> > > _counts;
	// _edgeCounts is the number of appearances in a super factor for a specific position
	std::vector<std::vector<std::vector<int> > > _edgeCounts;
	size_t _compressIters;
	std::map<size_t, size_t> _gndVarToSuperVar;
	std::map<size_t, size_t> _gndFacToSuperFac;
	std::map<size_t, V> _varRepr;
	std::map<size_t, size_t> _facRepr;

public:
	/// Default constructor
	FactorGraph() : G(), _vars(), _factors() {}

	FactorGraph(const dai::FactorGraph & fg) : G(), _vars(), _factors() {
		std::vector<F> newFactors(fg.factors().size());
		for (size_t i=0; i<fg.factors().size(); i++) {
			newFactors[i] = F(fg.factor(i));
		}
		create(newFactors);
		initCounts();
	}

	FactorGraph(const CFactorGraph & cfg);


	FactorGraph(const std::vector<F> &factors) : G(), _vars(), _factors() {
		create(factors);
		initCounts();
	}

	// create FactorGraph from a, possibly compressed, CFactorGraph
	FactorGraph(const CBP & cbp) : G(), _vars(), _factors() {
		std::vector<F> newFactors(cbp.factors().size());
		for (size_t i=0; i<cbp.factors().size(); i++) {
			newFactors[i] = F(cbp.factor(i));
		}
		create(newFactors);
		initCounts();

		_gndVarToSuperVar = std::map<size_t, size_t>(cbp.allReprV());

		int posCount;
		int negCount;
		std::vector<int> signs;
		std::map<size_t,int> counts;
		std::map<size_t,int>::const_iterator iter;
		for (size_t i=0; i<nrFactors();i++) {
			signs = createSigns(cbp.factor(i).p());

			foreach(const Neighbor tmpVar, nbF(i)) {
				counts = cbp.factor(i).counts(tmpVar.iter);

				posCount = 0;
				negCount = 0;
				_counts[tmpVar][tmpVar.dual][0] = 0;
				_counts[tmpVar][tmpVar.dual][1] = 0;
				for(iter=counts.begin(); iter != counts.end(); ++iter) {
					if (signs[cbp.factor(i).sigma()[iter->first]] == -1) {
						if (iter->second > 0) {
							_counts[tmpVar][tmpVar.dual][1] = iter->second;
						}
						posCount++;
					} else {
						if (iter->second > 0) {
							_counts[tmpVar][tmpVar.dual][0] = iter->second;
						}
						negCount++;
					}
				}
				_edgeCounts[tmpVar][tmpVar.dual][0] = negCount;
				_edgeCounts[tmpVar][tmpVar.dual][1] = posCount;
			}
		}
	}

	/// Returns constant reference the \a i 'th variable
	const V & var(size_t i) const { return _vars[i]; }
	/// Returns constant reference to all variables
	const std::vector<V> & vars() const { return _vars; }

	/// Returns constant reference to \a I 'th factor
	const F & factor(size_t I) const { return _factors[I]; }
	/// Returns constant reference to all factors
	const std::vector<F> & factors() const { return _factors; }

	/// Returns number of variables
	size_t nrVars() const { return vars().size(); }
	/// Returns number of factors
	size_t nrFactors() const { return factors().size(); }
	/// Calculates number of edges
	size_t nrEdges() const { return G.nrEdges(); }

	/// Returns the index of a particular variable
	/** \note Time complexity: O(nrVars())
	 *  \throw OBJECT_NOT_FOUND if the variable is not part of this factor graph
	 */
	size_t findVar( const V &n ) const {
		size_t i = find( vars().begin(), vars().end(), n ) - vars().begin();
		if( i == nrVars() )
			DAI_THROW(OBJECT_NOT_FOUND);
		return i;
	}

	/// Returns constant reference to neighbors of the \a i 'th variable
	const Neighbors & nbV( size_t i ) const { return G.nb1(i); }
	/// Returns constant reference to neighbors of the \a I 'th factor
	const Neighbors & nbF( size_t I ) const { return G.nb2(I); }
	/// Returns constant reference to the \a _I 'th neighbor of the \a i 'th variable
	const Neighbor & nbV( size_t i, size_t _I ) const { return G.nb1(i)[_I]; }
	/// Returns constant reference to the \a _i 'th neighbor of the \a I 'th factor
	const Neighbor & nbF( size_t I, size_t _i ) const { return G.nb2(I)[_i]; }

	void setVerbose(size_t level) {
		_verbose = level;
	}

	void readFromFile( const char *filename );

	void compress();
	void initCounts();
	void initMappings() {
		// default initialization of var mapping
		for (size_t i=0; i<nrVars(); i++) {
			_gndVarToSuperVar[i] = i;
		}
		for (size_t i=0; i<nrFactors(); i++) {
			_gndFacToSuperFac[i] = i;
		}
	}

	void createVarMapping(CompressInterface* compressAlg, CFactorGraph& liftedFg) {
		for (size_t i=0; i<compressAlg->getVarColorVec().size();i++) {
			_gndVarToSuperVar[i] = findVar(compressAlg->getVarRepr().find(compressAlg->getVarColorVec()[i])->second);
		}
	}

	void createFactorMapping(CompressInterface* compressAlg, CFactorGraph& liftedFg) {
		size_t cluster = 0;
		for (std::map<size_t, std::vector<size_t> >::const_iterator mapIter=compressAlg->getFacClustering().begin(); mapIter!=compressAlg->getFacClustering().end(); mapIter++) {
			for (size_t i=0; i<mapIter->second.size(); i++) {
				_gndFacToSuperFac[mapIter->second[i]] = cluster;
			}
			cluster++;
		}
	}

	std::map<size_t, size_t> & gndVarToSuperVar () {
		return _gndVarToSuperVar;
	}

	std::map<size_t, size_t> & gndFacToSuperFac () {
		return _gndFacToSuperFac;
	}

	const size_t getCompressIterations() const {
		return _compressIters;
	}

	// only lifted case
	// each cluster is determined by the representative variable
	std::vector<size_t> clusterV(size_t superVar) const {
		std::vector<size_t> cluster;
		for (std::map<size_t, size_t>::const_iterator iter=_gndVarToSuperVar.begin(); iter != _gndVarToSuperVar.end(); iter++) {
			if (iter->second == superVar) {
				cluster.push_back(iter->first);
			}
		}
		return cluster;
	}

	std::vector<size_t> clusterF(size_t superFac) const {
		std::vector<size_t> cluster;
		for (std::map<size_t, size_t>::const_iterator iter=_gndFacToSuperFac.begin(); iter != _gndFacToSuperFac.end(); iter++) {
			if (iter->second == superFac) {
				cluster.push_back(iter->first);
			}
		}
		return cluster;
	}

	const size_t reprV(size_t gndIdx) const {
		return _gndVarToSuperVar.find(gndIdx)->second;
	}

	const size_t reprF(size_t gndIdx) const {
		return _gndFacToSuperFac.find(gndIdx)->second;
	}

	const std::vector<std::vector<std::vector<int> > > & counts() const { return _counts; }
	std::vector<std::vector<std::vector<int> > > & counts() { return _counts; }
	const std::vector<int>& count(size_t i, size_t _I) const { return _counts[i][_I]; }
	std::vector<int>& count(size_t i, size_t _I) { return _counts[i][_I]; }

	std::vector<std::vector<std::vector<int> > > & edgeCounts() { return _edgeCounts; }
	const std::vector<std::vector<std::vector<int> > > & edgeCounts() const { return _edgeCounts; }
	const std::vector<int>& edgeCount(size_t i, size_t _I) const { return _edgeCounts[i][_I]; }
	std::vector<int>& edgeCount(size_t i, size_t _I) { return _edgeCounts[i][_I]; }

	void fixVars(const std::map<dai::Var, int>& vars);
	void fixVarClusters(const std::map<dai::Var, int>& varClusters);

	static FactorGraph createCleanedGraph(const FactorGraph<F, V> & wfg, const std::map<long, int>& fixedVariables) throw (const char*);

	static std::vector<int> createSigns (const dai::TProb<dai::Real> &table) {
		size_t nrVars = (size_t)log2(table.size());
		std::vector<int> signs(nrVars, 0);

		std::bitset<sizeof(size_t) * 8 > bs;
		for (size_t j = 0; j < table.size(); j++) {
			if (table[j] <= 1) {
				bs = std::bitset<sizeof(size_t) * 8 >(j);
				break;
			}
		}

		for (size_t j=0; j<nrVars; j++) {
			if (bs[j] == 0) {
				signs[j] = -1;
			} else {
				signs[j] = 1;
			}
		}

		return signs;
	}

private:
	void create(const std::vector<F> &factors ) {
		// add factors, obtain variables
		std::set<V> varset;
		_factors.reserve( factors.size() );
		size_t nrEdges = 0;
		for(typename std::vector<F>::const_iterator p2 = factors.begin(); p2 != factors.end(); p2++ ) {
			_factors.push_back( *p2 );
			copy( p2->vars().begin(), p2->vars().end(), inserter( varset, varset.begin() ) );
			nrEdges += p2->vars().size();
		}

		// add vars
		_vars.reserve( varset.size() );
		for(typename std::set<V>::const_iterator p1 = varset.begin(); p1 != varset.end(); p1++ ) {
			_vars.push_back( *p1 );
		}
		// create graph structure
		constructGraph( nrEdges );
		initMappings();
	}

	/// Part of constructors (creates edges, neighbors and adjacency matrix)
	void constructGraph( size_t nrEdges ) {
		// create a mapping for indices
		dai::hash_map<size_t, size_t> hashmap;

		for( size_t i = 0; i < vars().size(); i++ ) {
			hashmap[var(i).label()] = i;
		}

		// create edge list
		std::vector<Edge> edges;
		edges.reserve( nrEdges );
		for( size_t i2 = 0; i2 < nrFactors(); i2++ ) {
			const dai::SmallSet<V>& ns = factor(i2).vars();
			for(typename dai::SmallSet<V>::const_iterator q = ns.begin(); q != ns.end(); q++ ) {
				edges.push_back( Edge(hashmap[q->label()], i2) );
			}
		}

		// create bipartite graph
		G.construct( nrVars(), nrFactors(), edges.begin(), edges.end() );
	}


	void initVarColors(std::vector<size_t> &varSigs, std::vector<std::vector<size_t> > &varInbox);
	void initFacColors(std::vector<size_t> &facSigs, std::vector<std::vector<size_t> > &facInbox);

	void iterate(
			std::vector<size_t> &varSigs,
			std::vector<size_t> &facSigs,
			std::vector<std::vector<size_t> > &varInbox,
			std::vector<std::vector<size_t> > &facInbox,
			std::map<size_t, size_t> &varClustering,
			std::map<size_t, size_t> &facClustering);

	size_t facColor(size_t i, std::vector<size_t> &varSigs, std::vector<size_t> &facSigs, std::vector<std::vector<size_t> > &facInbox);
	size_t varColor(size_t i, std::vector<size_t> &varSigs, std::vector<size_t> &facSigs, std::vector<std::vector<size_t> > &varInbox);

	bool hasConverged(std::map<size_t, size_t> &varClustering);

	void createFacColorVec(std::vector<size_t> &facColorVec, std::map<size_t,size_t> &facClustering,  std::vector<size_t> &facSigs);
	void createVarColorVec(std::vector<size_t> &varColorVec, std::map<size_t, size_t> &varClustering,  std::vector<size_t> &varSigs);
	void createCompressed(const std::map<size_t, size_t> &facClustering, const std::map<size_t, size_t> &varClustering, std::vector<size_t> &facColorVec, std::vector<size_t> &varColorVec);
	void createCounts();


};

template<class F, class V>
double FactorGraph<F,V>::ZERO_VAL = 0.000001;
template<class F, class V>
double FactorGraph<F,V>::NON_ZERO_VAL = 10000000;

template<class F, class V>
std::ostream& operator<< (std::ostream& os, const FactorGraph<F, V>& fg);

}

#endif /* FACTORGRAPH_H_ */

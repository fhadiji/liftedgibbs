/*
 * CFactorGraph.h
 *
 *  Created on: Aug 24, 2009
 *      Author: stream
 */

#ifndef CFACTORGRAPH_H_
#define CFACTORGRAPH_H_

#include <iostream>
#include <map>
#include <dai/bipgraph.h>
#include <dai/factorgraph.h>
#include <STREAM/CBP/CFactor.h>

namespace stream {

class CFactorGraph {
public:
	/// Stores the neighborhood structure
	dai::BipartiteGraph G;

	/// Shorthand for BipartiteGraph::Neighbor
	typedef dai::BipartiteGraph::Neighbor Neighbor;

	/// Shorthand for BipartiteGraph::Neighbors
	typedef dai::BipartiteGraph::Neighbors Neighbors;

	/// Shorthand for BipartiteGraph::Edge
	typedef dai::BipartiteGraph::Edge Edge;

private:
	std::vector<dai::Var> _vars;
	std::vector<CFactor> _factors;
	/// maps the labe of a var to a certain state
	std::map<long, size_t> _evidence;
	std::map<size_t,CFactor>  _backup;
	/// sets the verbosity level
	long _verbose;

public:
	/// Default constructor
	CFactorGraph() :
		G(), _vars(), _factors(), _evidence(), _backup(), _verbose(1) {
	}

	/// Copy constructor
	CFactorGraph(const CFactorGraph & x) :
		G(x.G), _vars(x._vars), _factors(x._factors), _evidence(x._evidence), _backup(x._backup), _verbose(1) {
	}

	/// Assignment operator
	CFactorGraph & operator=(const CFactorGraph & x) {
		if (this != &x) {
			G = x.G;
			_vars = x._vars;
			_factors = x._factors;
			_evidence = x._evidence;
			_backup = x._backup;
			_verbose = x._verbose;
		}
		return *this;
	}

	/// Constructs a FactorGraph from a vector of CFactors
	CFactorGraph(const std::vector<CFactor>& factors);
	CFactorGraph(const std::vector<dai::Factor>& factors);

	CFactorGraph(const dai::FactorGraph & fg) {
		std::vector<CFactor> cfactors;
		for (size_t i = 0; i < fg.nrFactors(); i++) {
			cfactors.push_back(CFactor(fg.factor(i)));
		}
		CFactorGraph::operator=( CFactorGraph(cfactors) );
	}

	/// Destructor
	virtual ~CFactorGraph() {
	}

	/// Create (virtual default constructor)
	virtual CFactorGraph* create() const {
		return new CFactorGraph();
	}

	/// Clone *this (virtual copy constructor)
	virtual CFactorGraph* clone() const {
		return new CFactorGraph(*this);
	}

	/// Returns const reference to i'th variable
	const dai::Var & var(size_t i) const {
		return _vars[i];
	}
	/// Returns const reference to all factors
	const std::vector<dai::Var> & vars() const {
		return _vars;
	}
	/// Returns reference to I'th factor
	CFactor & factor(size_t I) {
		return _factors[I];
	}
	/// Returns const reference to I'th factor
	const CFactor & factor(size_t I) const {
		return _factors[I];
	}
	/// Returns const reference to all factors
	const std::vector<CFactor> & factors() const {
		return _factors;
	}
	/// Returns number of variables
	size_t nrVars() const {
		return vars().size();
	}
	/// Returns number of factors
	size_t nrFactors() const {
		return factors().size();
	}
	/// Calculates number of edges
	size_t nrEdges() const {
		return G.nrEdges();
	}

	/// Provides read access to neighbors of variable
	const Neighbors & nbV(size_t i) const {
		return G.nb1(i);
	}
	/// Provides full access to neighbors of variable
	Neighbors & nbV(size_t i) {
		return G.nb1(i);
	}
	/// Provides read access to neighbors of factor
	const Neighbors & nbF(size_t I) const {
		return G.nb2(I);
	}
	/// Provides full access to neighbors of factor
	Neighbors & nbF(size_t I) {
		return G.nb2(I);
	}
	/// Provides read access to neighbor of variable
	const Neighbor & nbV(size_t i, size_t _I) const {
		return G.nb1(i)[_I];
	}
	/// Provides full access to neighbor of variable
	Neighbor & nbV(size_t i, size_t _I) {
		return G.nb1(i)[_I];
	}
	/// Provides read access to neighbor of factor
	const Neighbor & nbF(size_t I, size_t _i) const {
		return G.nb2(I)[_i];
	}
	/// Provides full access to neighbor of factor
	Neighbor & nbF(size_t I, size_t _i) {
		return G.nb2(I)[_i];
	}

    /// Returns the index of a particular variable
    size_t findVar( const dai::Var & n ) const {
        size_t i = find( vars().begin(), vars().end(), n ) - vars().begin();
        return i;
    }

	void setSigma(size_t facIdx, std::vector<size_t> sigma) {
		factor(facIdx).sigma() = sigma;
		factor(facIdx).position().clear();
		for (size_t i = 0; i < factor(facIdx).vars().size(); i++) {
			factor(facIdx).position()[nbF(facIdx, i)] = sigma[i];
		}
	}

	const std::vector<size_t> & getSigma (size_t facIdx) const {
		return factor(facIdx).sigma();
	}

	void setPosition(size_t facIdx, std::map<size_t, size_t> &position) {
		factor(facIdx).position() = position;
	}

	const std::map<size_t, size_t> & getPosition (size_t facIdx) const {
		return factor(facIdx).position();;
	}

	void resetPositions() {
		for (size_t i = 0; i < nrFactors(); i++) {
			for (std::map<size_t, size_t>::iterator iter = factor(i).position().begin(); iter
					!= factor(i).position().end(); ++iter) {
				factor(i).position()[iter->first] = 0;
			}
		}
	}

	const std::map<long, size_t> & evidence() {
		return _evidence;
	}

	void setVerbose (long value) {
		_verbose = value;
	}

    /// Set the content of the I'th factor and make a backup of its old content if backup == true
    virtual void setFactor( size_t I, const CFactor &newFactor, bool backup = false ) {
        assert( newFactor.vars() == factor(I).vars() );
        if( backup )
            backupFactor( I );
        _factors[I] = newFactor;
    }

    /// Set the contents of all factors as specified by facs and make a backup of the old contents if backup == true
    virtual void setFactors( const std::map<size_t, CFactor> & facs, bool backup = false ) {
        for( std::map<size_t, CFactor>::const_iterator fac = facs.begin(); fac != facs.end(); fac++ ) {
            if( backup )
                backupFactor( fac->first );
            setFactor( fac->first, fac->second );
        }
    }

	/// Clamp variable n to value i (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$);
	virtual void clamp(const dai::Var & n, size_t i, bool backup = false  );

	/// Clamp MLN Factor I to value state (i.e. multiply with a Kronecker delta \f$\delta_{x_n, i}\f$);
	virtual void clampMlnFactor(size_t I, bool clampState, float weight, bool backup = false);

	/// Set all factors interacting with the i'th variable 1
	virtual void makeCavity(unsigned i, bool backup );

	/// Restores the I'th Factor from the backup (it should be backed up first)
	void backupFactor( size_t I );

	/// Backup the factors specified by indices in facs
	virtual void backupFactors( const std::set<size_t> & facs );

	/// Makes a backup of all factors connected to a set of variables
	void backupFactors( const dai::VarSet &ns );


	void restoreFactor( size_t I );
	void restoreFactors( const dai::VarSet &ns );
	void restoreFactors();

	static void readFactors (std::vector<CFactor>* facs, std::istream& is, long verbose);

	/// Read the factors from a file
//	static std::vector<CFactor> readFactorsFromFile(const char *filename);

	/// Reads a FactorGraph from a file
	void readFromFile(const char *filename);

	/// Writes a FactorGraph to a file
	void writeToFile(const char *filename) const;

	/// Writes a FactorGraph to a GraphViz .dot file
	void printDot(std::ostream& os) const;

	static CFactorGraph createCleanedGraph (CFactorGraph& cfg, std::map<long, int>& fixedVariables) throw (const char*) ;

	// Friends
	friend std::ostream& operator <<(std::ostream& os, const CFactorGraph& fg);
	friend std::istream& operator >>(std::istream& is, CFactorGraph& fg);

private:
	/// Part of constructors (creates edges, neighbors and adjacency matrix)
	void constructGraph(size_t nrEdges);
};

} // end namespace

#endif /* CFACTORGRAPH_H_ */

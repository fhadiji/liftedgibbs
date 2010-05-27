/*
 * Compress.h
 *
 *  Created on: Jul 27, 2009
 *      Author: stream
 */

#ifndef COMPRESS_H_
#define COMPRESS_H_

#include <boost/functional/hash.hpp>
#include <sstream>
#include "STREAM/CBP/CFactorGraph.h"
#include <dai/daialg.h>
#include <dai/factorgraph.h>
#include <dai/util.h>
#include <dai/factor.h>
#include <dai/properties.h>

namespace stream {

class CompressInterface {
protected:
	// gives access to original graph
	CFactorGraph _cfg;
	// maps var idx to cluster idx
	std::vector<size_t> _varColorVec;
	// maps factor idx to cluster idx
	std::vector<size_t> _facColorVec;
	// contains the signatures of factors
	typedef size_t Signature;
	std::vector<Signature> _facSigs;
	// contains the signatures of a variables
	std::vector<Signature> _varSigs;
	// stores the incoming colors of the variables
	std::vector<std::vector<size_t> > _varInbox;
	// stores the incoming colors of the factors
	std::vector<std::vector<size_t> > _facInbox;
	// saves variable representatives (cluster idx => var idx)
	std::map<size_t, dai::Var> _varRepr;
	// saves factor representatives (cluster idx => fac idx)
	std::map<size_t, size_t> _facRepr;
	// verbose mode works similar to the one in CBP
	size_t _verbose;
	// for hashing color boxes
	boost::hash<std::vector<size_t> > hashVector;

	int _iters;
	size_t _varClusterCountOld;
	size_t _facClusterCountOld;
	std::map<size_t, std::vector<size_t> > _facClustering;
	std::map<size_t, std::vector<size_t> > _varClustering;
	std::vector<size_t> _initialVarColors;

public:
	CompressInterface() : _verbose(0) {
	}

	CFactorGraph compress(const CFactorGraph& cfg);
	CFactorGraph compress(const CFactorGraph& cfg, const std::vector<size_t> &varColors);
	virtual void execute();
	void init();
	bool hasConverged();
	void iterate();
	CFactorGraph createCFactorGraph();

	std::vector<size_t>& getVarColorVec() {
		return _varColorVec;
	}
	std::vector<size_t>& getFacColorVec() {
		return _facColorVec;
	}
	std::map<size_t, dai::Var> getVarRepr() {
		return _varRepr;
	}
	std::map<size_t, size_t> getFacRepr() {
		return _facRepr;
	}

	std::map<size_t, std::vector<size_t> >& getFacClustering () {
		return _facClustering;
	}

	void setVerbose(size_t verbose) {
		_verbose = verbose;
	}
	void setCfg(const CFactorGraph& cfg) {
		_cfg = cfg;
	}

	const int getIterations() const {
		return _iters;
	}

	static void rearrangeCnfFactor(stream::CFactor &cnfClause);

protected:
	virtual size_t facColor(size_t i);
	virtual size_t varColor(size_t i);
	virtual void initVarColors();
	virtual void initFacColors();
	virtual void createFacClustering();
	virtual void createVarClustering();
	virtual std::string name() {return std::string(""); }
	virtual std::vector<std::map<size_t, int> > createCounts(size_t &gndFactor, dai::VarSet &superVarSet);
};

class SimpleCompress: public CompressInterface {
protected:
	virtual std::string name() {
		return std::string("SimpleCompress");
	}
};

class PositionCompress: public CompressInterface {
protected:
	virtual std::string name() {
		return std::string("PositionCompress");
	}
	virtual void initFacColors();
	virtual size_t facColor(size_t i);
	virtual size_t varColor(size_t i);
};

class AdvancedCompress: public CompressInterface {
public:
	virtual void execute();
protected:
	virtual std::string name() {
		return std::string("AdvancedCompress");
	}
};

class AnyPositionCompress: public CompressInterface {
protected:
	virtual std::string name() {
		return std::string("AnyPositionCompress");
	}
	virtual void initFacColors();
	virtual size_t facColor(size_t i);
	virtual size_t varColor(size_t i);
	virtual std::vector<std::map<size_t, int> > createCounts(size_t &gndFactor, dai::VarSet &superVarSet);
private:
	void print(const std::vector<size_t> & v, const int size);
	void permute(std::vector<size_t> & v, const int start, const int n,
			std::vector<std::vector<size_t> > & perms);
};

class AnyPositionCnfCompress: public AnyPositionCompress {

private:
	// zero states in ordered factors (ordered = first pos vars then negative ones)
	std::vector<size_t> _zeroStates;

	virtual std::string name() {
		return std::string("AnyPositionCnfCompress");
	}
	virtual void initFacColors();
	virtual size_t facColor(size_t i);
	virtual std::vector<std::map<size_t, int> > createCounts(size_t &gndFactor, dai::VarSet &superVarSet);
};

}

#endif /* COMPRESS_H_ */

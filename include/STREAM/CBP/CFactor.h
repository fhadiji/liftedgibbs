/*
 * CFactor.h
 *
 *  Created on: Aug 24, 2009
 *      Author: stream
 */

#ifndef CFACTOR_H_
#define CFACTOR_H_

#include <iostream>
#include <cmath>
#include <dai/prob.h>
#include <dai/index.h>
#include <dai/varset.h>
#include <dai/factor.h>


namespace stream {

// predefine TCFactor<T> class
template<typename T> class TCFactor;

/// Represents a factor with probability entries represented as Real
typedef TCFactor<dai::Real> CFactor;

/// Represents a probability factor.
/** A \e factor is a function of the Cartesian product of the state
 *  spaces of some set of variables to the nonnegative real numbers.
 *  More formally, if \f$x_i \in X_i\f$ for all \f$i\f$, then a factor
 *  depending on the variables \f$\{x_i\}\f$ is a function defined
 *  on \f$\prod_i X_i\f$ with values in \f$[0,\infty)\f$.
 *
 *  A Factor has two components: a VarSet, defining the set of variables
 *  that the factor depends on, and a TProb<T>, containing the values of
 *  the factor for all possible joint states of the variables.
 *
 *  \tparam T Should be castable from and to double.
 */
template<typename T> class TCFactor {
private:
	dai::VarSet _vs;
	dai::TProb<T> _p;
	/// stores original ordering of the variables in the factor
	std::vector<size_t> _sigma;
	/// stores for the factor a map which contains the position of each var
	std::map<size_t, size_t> _position;
	// count
	std::vector<std::map<size_t,int> > _counts;


public:
	/// Construct Factor with empty VarSet
	TCFactor(dai::Real p = 1.0) :
		_vs(), _p(1, p), _sigma(), _position(), _counts() {
	}

	/// Construct Factor from VarSet
	TCFactor(const dai::VarSet& ns) :
		_vs(ns), _p(_vs.nrStates()) {
		defaults();
	}

	/// Construct Factor from VarSet and initial value
	TCFactor(const dai::VarSet& ns, dai::Real p) :
		_vs(ns), _p(_vs.nrStates(), p) {
		defaults();
	}

	/// Construct Factor from VarSet and initial array
	TCFactor(const dai::VarSet& ns, const dai::Real *p) :
		_vs(ns), _p(_vs.nrStates(), p) {
		defaults();
	}

	/// Construct Factor from VarSet and TProb<T>
	TCFactor(const dai::VarSet& ns, const dai::TProb<T>& p) :
		_vs(ns), _p(p) {
		defaults();
	}

	/// Construct Factor from Var
	TCFactor(const dai::Var& n) :
		_vs(n), _p(n.states()) {
		defaults();
	}

	TCFactor(const dai::Factor &factor) : _vs(factor.vars()), _p(factor.p()) {
		defaults();
	}

	/// Copy constructor
	TCFactor(const TCFactor<T> &x) :
		_vs(x._vs), _p(x._p), _sigma(x._sigma), _position(x._position), _counts(x._counts){
	}

	/// Assignment operator
	TCFactor<T> & operator=(const TCFactor<T> &x) {
		if (this != &x) {
			_vs = x._vs;
			_p = x._p;
			_sigma = x._sigma;
			_position = x._position;
			_counts = x._counts;
		}
		return *this;
	}

	/// Returns const reference to probability entries
	const dai::TProb<T> & p() const {
		return _p;
	}
	/// Returns reference to probability entries
	dai::TProb<T> & p() {
		return _p;
	}

	/// Returns const reference to variables
	const dai::VarSet & vars() const {
		return _vs;
	}

	/// Returns the number of possible joint states of the variables
	size_t states() const {
		return _p.size();
	}


	const std::vector<size_t> & sigma() const {
		return _sigma;
	}

	std::vector<size_t> & sigma() {
		return _sigma;
	}

	const std::map<size_t, size_t> & position() const {
		return _position;
	}

	std::map<size_t, size_t> & position() {
		return _position;
	}

	const std::vector<std::map<size_t,int> > & counts() const {
		return _counts;
	}

	std::vector<std::map<size_t, int> > & counts() {
		return _counts;
	}

	const std::map<size_t,int> & counts(size_t varPos) const {
		return _counts[varPos];
	}

	std::map<size_t,int> & counts(size_t varPos) {
		return _counts[varPos];
	}

	const int & count(size_t varPos, size_t pos) const {
		return _counts[varPos][pos];
	}

	int & count(size_t varPos, size_t pos) {
		return _counts[varPos][pos];
	}

	void resize(size_t n) {
		_p = dai::TProb<T>(n,0.0);
	}

	/// Returns a copy of the i'th probability value
	T operator[](size_t i) const {
		return _p[i];
	}

	/// Returns a reference to the i'th probability value
	T& operator[](size_t i) {
		return _p[i];
	}

	/// Sets all probability entries to p
	TCFactor<T> & fill(T p) {
		_p.fill(p);
		return (*this);
	}

	/// Fills all probability entries with random values
	TCFactor<T> & randomize() {
		_p.randomize();
		return (*this);
	}

	/// Returns product of *this with x
	TCFactor<T> operator*(T x) const {
		CFactor result = *this;
		result.p() *= x;
		return result;
	}

	/// Multiplies each probability entry with x
	TCFactor<T>& operator*=(T x) {
		_p *= x;
		return *this;
	}

	/// Returns quotient of *this with x
	TCFactor<T> operator/(T x) const {
		CFactor result = *this;
		result.p() /= x;
		return result;
	}

	/// Divides each probability entry by x
	TCFactor<T>& operator/=(T x) {
		_p /= x;
		return *this;
	}

	/// Returns product of *this with another Factor
	TCFactor<T> operator*(const TCFactor<T>& Q) const;

	/// Returns quotient of *this with another Factor
	TCFactor<T> operator/(const TCFactor<T>& Q) const;

	/// Multiplies *this with another Factor
	TCFactor<T>& operator*=(const TCFactor<T>& Q) {
		return (*this = (*this * Q));
	}

	/// Divides *this by another Factor
	TCFactor<T>& operator/=(const TCFactor<T>& Q) {
		return (*this = (*this / Q));
	}

	/// Returns sum of *this and another Factor (their vars() should be identical)
	TCFactor<T> operator+(const TCFactor<T>& Q) const {
#ifdef DAI_DEBUG
		assert( Q._vs == _vs );
#endif
		TCFactor<T> sum(*this);
		sum._p += Q._p;
		return sum;
	}

	/// Returns difference of *this and another Factor (their vars() should be identical)
	TCFactor<T> operator-(const TCFactor<T>& Q) const {
#ifdef DAI_DEBUG
		assert( Q._vs == _vs );
#endif
		TCFactor<T> sum(*this);
		sum._p -= Q._p;
		return sum;
	}

	/// Adds another Factor to *this (their vars() should be identical)
	TCFactor<T>& operator+=(const TCFactor<T>& Q) {
#ifdef DAI_DEBUG
		assert( Q._vs == _vs );
#endif
		_p += Q._p;
		return *this;
	}

	/// Subtracts another Factor from *this (their vars() should be identical)
	TCFactor<T>& operator-=(const TCFactor<T>& Q) {
#ifdef DAI_DEBUG
		assert( Q._vs == _vs );
#endif
		_p -= Q._p;
		return *this;
	}

	/// Adds scalar to *this
	TCFactor<T>& operator+=(T q) {
		_p += q;
		return *this;
	}

	/// Subtracts scalar from *this
	TCFactor<T>& operator-=(T q) {
		_p -= q;
		return *this;
	}

	/// Returns sum of *this and a scalar
	TCFactor<T> operator+(T q) const {
		TCFactor<T> result(*this);
		result._p += q;
		return result;
	}

	/// Returns difference of *this with a scalar
	TCFactor<T> operator-(T q) const {
		TCFactor<T> result(*this);
		result._p -= q;
		return result;
	}

	/// Returns *this raised to some power
	TCFactor<T> operator^(dai::Real a) const {
		TCFactor<T> x;
		x._vs = _vs;
		x._p = _p ^ a;
		return x;
	}

	/// Raises *this to some power
	TCFactor<T>& operator^=(dai::Real a) {
		_p ^= a;
		return *this;
	}

    /// Equal-to operator (compares varset and states)
    bool operator == ( const TCFactor<T>& fac ) const {
    	if (!(_vs == fac._vs)) {
    		return false;
    	}
    	if (_p.size() == fac._p.size() ) {
    		if (dai::dist(_p, fac._p, dai::Prob::DISTLINF) != 0) {
    			return false;
    		}
    	} else {
    		return false;
    	}
    	return true;
    }

	/// Returns exp of *this
	TCFactor<T> exp() const {
		TCFactor<T> e;
		e._vs = _vs;
		e._p = _p.exp();
		return e;
	}

	/// Returns logarithm of *this
	TCFactor<T> log() const {
		TCFactor<T> l;
		l._vs = _vs;
		l._p = _p.log();
		return l;
	}

	/// Normalizes *this Factor
	T normalize(typename dai::Prob::NormType norm = dai::Prob::NORMPROB) {
		return _p.normalize(norm);
	}

	/// Returns a normalized copy of *this
	TCFactor<T> normalized(typename dai::Prob::NormType norm = dai::Prob::NORMPROB) const {
		TCFactor<T> result;
		result._vs = _vs;
		result._p = _p.normalized(norm);
		return result;
	}

    /// Returns unnormalized marginal; ns should be a subset of vars()
    TCFactor<T> partSum(const dai::VarSet & ns) const;

	/// Returns (normalized by default) marginal; ns should be a subset of vars()
	TCFactor<T> marginal(const dai::VarSet & ns, bool normed = true) const {
		if (normed)
			return partSum(ns).normalized();
		else
			return partSum(ns);
	}



private:

	void defaults() {
		defaultSigma();
		defaultPosition();
		defaultCount();
	}

	void defaultSigma () {
		_sigma.reserve(_vs.size());
		for (size_t j = 0; j < _vs.size(); j++) {
			_sigma.push_back(j);
		}
	}

	void defaultPosition () {
		size_t pos = 0;
		for (std::vector<dai::Var>::const_iterator iter=_vs.begin(); iter<_vs.end(); iter++) {
			_position[iter->label()] = pos;
			pos++;
		}
	}

	void defaultCount () {
		// count needs to have size of _vs because of index generation in lifted case
		for (std::vector<dai::Var>::const_iterator iter=_vs.begin(); iter<_vs.end(); iter++) {
			std::map<size_t, int> count;
			count[iter - _vs.begin()] = 1;
			_counts.push_back(count);
		}
	}


};

template<typename T> TCFactor<T> TCFactor<T>::partSum(const dai::VarSet & ns) const {
#ifdef DAI_DEBUG
    assert( ns << _vs );
#endif

    TCFactor<T> res( ns, 0.0 );

    dai::IndexFor i_res( ns, _vs );
    for( size_t i = 0; i < _p.size(); i++, ++i_res )
        res._p[i_res] += _p[i];

    return res;
}

template<typename T> TCFactor<T> TCFactor<T>::operator*(const TCFactor<T>& Q) const {
	TCFactor<T> prod(_vs | Q._vs, 0.0);

	dai::IndexFor i1(_vs, prod._vs);
	dai::IndexFor i2(Q._vs, prod._vs);

	for (size_t i = 0; i < prod._p.size(); i++, ++i1, ++i2)
		prod._p[i] += _p[i1] * Q._p[i2];

	return prod;
}

template<typename T> TCFactor<T> TCFactor<T>::operator/(const TCFactor<T>& Q) const {
	TCFactor<T> quot(_vs + Q._vs, 0.0);

	dai::IndexFor i1(_vs, quot._vs);
	dai::IndexFor i2(Q._vs, quot._vs);

	for (size_t i = 0; i < quot._p.size(); i++, ++i1, ++i2)
		quot._p[i] += _p[i1] / Q._p[i2];

	return quot;
}

/// Writes a Factor to an output stream
template<typename T> std::ostream& operator<<(std::ostream& os,
		const TCFactor<T>& P) {
	os << "(" << P.vars() << " <";
	for (size_t i = 0; i < P.states(); i++)
		os << P[i] << " ";
	os << ">; S<";
	for (size_t i=0; i<P.sigma().size(); i++) {
		os << P.sigma()[i] << " ";
	}
	os << ">; P{";
    for(std::map<size_t, size_t>::const_iterator iter = P.position().begin(); iter != P.position().end(); ++iter ) {
    	os << iter->first << "->" << iter->second << " ";
    }
	os << "}; C";

	dai::operator <<(os, P.counts());

	os << "})";


	return os;
}

} /* end namespace */

#endif /* CFACTOR_H_ */

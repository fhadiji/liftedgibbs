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
 * CBP.h
 *
 *  Created on: Jul 13, 2009
 *      Author: stream
 */

#ifndef CBP_H_
#define CBP_H_

#include <STREAM/CBP/CFactorGraph.h>
#include <STREAM/CBP/Compress.h>
#include <string>
#include <dai/daialg.h>
#include <dai/properties.h>
#include <dai/enum.h>
#include <dai/factorgraph.h>

namespace stream {

class CBP : public CFactorGraph {

	friend class CompressInterface;

	private:
        typedef std::vector<size_t> ind_t;
        struct EdgeProp {
            std::vector<ind_t>  index;
            std::vector<dai::Prob>   message;
            std::vector<dai::Prob>   newMessage;
            std::vector<int> count;
            size_t nrPos;
        };
        std::vector<std::vector<EdgeProp> > _edges;
        /// Maximum difference encountered so far
        double _maxdiff;
        /// Number of iterations needed
        size_t _iters;

        // maps gnd var idxs to super var idxs
        std::map<size_t, size_t> _gndVarToSuperVar;

        // maps gnd factor idxs to super factor idxs
        std::map<size_t, size_t> _gndFacToSuperFac;

		// compression algorithm
		CompressInterface * _compressAlg;

		// initial coloring of the variables, needed for adapted clustering
		std::vector<size_t> _initialVarColors;

		/// Number of messages sent during message passing
		size_t _messages;

		/// number of iterations of the compression algorithm
		size_t _compressIters;

    public:
        /// Parameters of this inference algorithm
        struct Properties {
            /// Enumeration of possible update schedules
            DAI_ENUM(UpdateType,SEQFIX,SEQRND,SEQMAX,PARALL)

            /// Enumeration of inference variants
            DAI_ENUM(InfType,SUMPROD,MAXPROD)

            /// Verbosity
            size_t verbose;

            /// Maximum number of iterations
            size_t maxiter;

            /// Tolerance
            double tol;

            /// Do updates in logarithmic domain?
            bool logdomain;

            /// Damping constant
            double damping;

            /// Update schedule
            UpdateType updates;

            /// Type of inference: sum-product or max-product?
            InfType inference;

            // Do compression before BP?
            bool compress;

            // Name of the compression algorithm
            std::string compressionAlg;

        } props;

        /// Name of this inference algorithm
        static const char *Name;

    public:
        /// Default constructor
        CBP() : CFactorGraph(), _edges(), _maxdiff(0.0), _iters(0U), props() { }

        /// Copy constructor
        CBP( const CBP &x ) :  CFactorGraph(x), _edges(x._edges), _maxdiff(x._maxdiff), _iters(x._iters), _compressAlg(x._compressAlg) , props(x.props) { }

        /// Assignment operator
        CBP& operator=( const CBP &x ) {
            if( this != &x ) {
                CFactorGraph::operator=( x );
                _edges = x._edges;
                _maxdiff = x._maxdiff;
                _iters = x._iters;
                _compressAlg = x._compressAlg;
                props = x.props;
            }
            return *this;
        }

        // Construct from CFactorGraph cfg and PropertySet opts
        CBP( const CFactorGraph & cfg, const dai::PropertySet &opts ) : CFactorGraph(cfg), _edges(), _maxdiff(0.0), _iters(0U), props() {
            setProperties( opts );
            construct();
        }

        void compress ();

		template <typename T> bool veccomp (std::vector<T> &vec1, std::vector<T> &vec2) {
			if (vec1.size() != vec2.size())
				return false;
			for (size_t i=0; i<vec1.size(); i++) {
				if (vec1[i] != vec2[i])
					return false;
			}
			return true;
		}

		static size_t maxPos (const dai::Factor& fac) {
			size_t maxPos = -1;
			double maxVal = 0;
			for (size_t i=0; i<fac.states(); i++) {
				if (fac[i] > maxVal) {
					maxVal = fac[i];
					maxPos = i;
				}
			}
			return maxPos;
		}

		const dai::Prob & message(size_t i, size_t _I, size_t pos) const { return _edges[i][_I].message[pos]; }
        dai::Prob & message(size_t i, size_t _I, size_t pos) { return _edges[i][_I].message[pos]; }
    	const dai::Prob & newMessage(size_t i, size_t _I, size_t pos) const { return _edges[i][_I].newMessage[pos]; }
        dai::Prob & newMessage(size_t i, size_t _I, size_t pos) { return _edges[i][_I].newMessage[pos]; }

        virtual CBP* clone() const { return new CBP(*this); }
        virtual CBP* create() const { return new CBP(); }
        virtual std::string identify() const;
        virtual dai::Factor belief( const dai::Var &n ) const;
        virtual void init();
        virtual void init(std::vector<size_t> &varColors);
        virtual void init( const dai::VarSet &ns );
        virtual double run();
        virtual double maxDiff() const { return _maxdiff; }
        virtual size_t Iterations() const { return _iters; }
    	const size_t getMessages() const {
    		return _messages;
    	}
    	const size_t getCompressIterations() const {
    		return _compressIters;
    	}
        virtual dai::Factor beliefV( size_t i ) const;
        virtual dai::Factor beliefF( size_t I ) const;
        void createVarMapping(CompressInterface * compressAlg, CFactorGraph& liftedFg);
        void createFacMapping(CompressInterface * compressAlg, CFactorGraph& liftedFg);

    	virtual std::vector<size_t> clusterV(size_t superVar) const;
    	virtual std::vector<size_t> clusterF(size_t superFac) const;
        size_t reprV(size_t gndIdx) {
        	return _gndVarToSuperVar.find(gndIdx)->second;
        }
        size_t reprF(size_t gndIdx) {
        	return _gndFacToSuperFac.find(gndIdx)->second;
        }
       	const std::map<size_t, size_t> & allReprV() const {
			return _gndVarToSuperVar;
		}
        void setProperty( const std::string &key, const std::string &value );


    private:
		/// modified original methods to pass position
        ind_t & index(size_t i, size_t _I, size_t &pos) { return _edges[i][_I].index[pos]; }
        const ind_t & index(size_t i, size_t _I, size_t &pos) const { return _edges[i][_I].index[pos]; }
        /// addtional methods for CBP
        int & count(const size_t &var, const size_t &facIter, size_t &pos) { return _edges[var][facIter].count[pos]; }
        const int & count(const size_t &var, const size_t &facIter, size_t &pos) const { return _edges[var][facIter].count[pos]; }
        size_t & nrPos(const size_t &var, const size_t &facIter) { return _edges[var][facIter].nrPos; }
        const size_t & nrPos(const size_t &var, const size_t &facIter) const { return _edges[var][facIter].nrPos; }

        void calcNewMessage( size_t i, size_t _I, size_t &pos );
        void calcBeliefF( size_t I, dai::Prob &p ) const;
        void updateMessage( size_t i, size_t _I, size_t &pos ) {
            if( props.damping == 0.0 ) {
                message(i,_I, pos) = newMessage(i,_I, pos);
            } else {
                message(i,_I, pos) = (message(i,_I, pos) ^ props.damping) * (newMessage(i,_I, pos) ^ (1.0 - props.damping));
            }
        }


        void construct();
        /// Set Props according to the PropertySet opts, where the values can be stored as std::strings or as the type of the corresponding Props member
        void setProperties( const dai::PropertySet &opts );
        dai::PropertySet getProperties() const;
        std::string printProperties() const;

        // new methods
    	void createVarReprs ();
    	void createFacReprs ();
        void setCompressionAlg (std::string name) {
        	if (name.compare("Simple") == 0){
        		_compressAlg = new SimpleCompress;
        	} else if (name.compare("Position") == 0) {
        		_compressAlg = new PositionCompress;
        	}
#ifndef NO_FANCY_COMP_ALG
			else if (name.compare("AnyPosition") == 0) {
				_compressAlg = new AnyPositionCompress;
			} else if (name.compare("AnyPositionCnf") == 0) {
				_compressAlg = new AnyPositionCnfCompress;
			} else if (name.compare("Advanced") == 0) {
				_compressAlg = new AdvancedCompress;
			}
#endif
        	else {
        		std::cout << "ERROR: Wrong name ('" << name << "')for compression algorithm";
        	}
        }

        void setCompressionAlg (CompressInterface *alg) {
        	_compressAlg = alg;
        }
};

}

#endif /* CBP_H_ */

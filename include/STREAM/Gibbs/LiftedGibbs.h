/*
 * LiftedGibbs.h
 *
 *  Created on: Mar 26, 2010
 *      Author: Fabian Hadiji
 *
 *  based on gibbs.h from libDAI
 */

#ifndef LIFTEDGIBBS_H_
#define LIFTEDGIBBS_H_


#include <STREAM/FactorGraph.h>
#include <dai/daialg.h>
#include <dai/properties.h>


namespace stream {

class LiftedGibbs {

public:

	/// Parameters for Gibbs
	struct Properties {
		/// Total number of iterations
		size_t iters;

		/// Number of "burn-in" iterations
		size_t burnin;

		/// Verbosity (amount of output sent to stderr)
		size_t verbose;
	} props;

	/// Name of this inference algorithm
	static const char *Name;

	private:
		///
		FactorGraph<CFactor, dai::Var> _fg;
		/// Type used to store the counts of various states
		typedef std::vector<size_t> _count_t;
		/// Type used to store the joint state of all variables
		typedef std::vector<size_t> _state_t;
		/// Number of samples counted so far (excluding burn-in)
		size_t _sample_count;
		/// State counts for each variable
		std::vector<_count_t> _var_counts;
		/// Current joint state of all variables
		_state_t _states;

		/// schedule for calculating and drawing random vars
		std::vector<size_t> _schedule;

		/// Ground configurations for each super factor for each linear index
		std::vector<std::vector<std::vector<std::vector<size_t> > > > _configMapping;

		static std::string _name;

	public:
		/// Default constructor
		LiftedGibbs() : _fg(), _sample_count(0), _var_counts(), _states() {}

		/// Construct from FactorGraph \a fg and PropertySet \a opts
		/** \param opts Parameters @see Properties
		 */
		LiftedGibbs( const FactorGraph<CFactor, dai::Var> &fg, const dai::PropertySet &opts ) : _fg(fg), _sample_count(0), _var_counts(), _states() {
			setProperties( opts );
			construct();
		}

		virtual dai::Real run();
		virtual std::string identify() const { return _name + printProperties(); }
		virtual dai::Factor beliefV( size_t i ) const;
		virtual void setProperties( const dai::PropertySet &opts );
		virtual std::string printProperties() const;

		/// Draw the current joint state of all variables from a uniform random distribution
		void randomizeStates();

		void setSchedule(std::vector<size_t>& schedule) {_schedule = schedule; }
		const std::vector<size_t>& getSchedule() {return _schedule; }

		// uses internal factor graph
		void createConfigMapping();
		// uses original ground factor graph
		void createConfigMapping(CFactorGraph& cfg);

	private:
		/// Helper function for constructors
		void construct();
		/// Draw state of variable \a i randomly from its conditional distribution and update the current state
		void resampleVar( size_t i );
		/// Calculate conditional distribution of variable \a i, given the current state
		dai::Prob getVarDist( size_t i );
		/// Calculates linear index into factor \a I corresponding to the current state of a particular ground MB
		size_t getFactorEntry( const std::vector<size_t>& groundMB);
		/// Calculates the differences between linear indices into factor \a I corresponding with a state change of variable \a i
		size_t getFactorEntryDiff( const std::vector<size_t>& groundMB, size_t i );
		void updateCounts();
};

}

#endif /* LIFTEDGIBBS_H_ */

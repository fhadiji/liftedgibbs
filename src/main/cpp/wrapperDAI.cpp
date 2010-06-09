/*
    Copyright (C) 2010
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
#include <boost/python/call.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/docstring_options.hpp>
#include <boost/python/enum.hpp>
#include <boost/python/errors.hpp>
#include <boost/python/module.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/pure_virtual.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/suite/indexing/map_indexing_suite.hpp>
#include <boost/python/wrapper.hpp>

#include <iostream>
#include <string>
#include <fstream>

#include <dai/smallset.h>
#include <dai/varset.h>
#include <dai/var.h>
#include <dai/bp.h>
#include <dai/index.h>
#include <dai/jtree.h>
#include <dai/gibbs.h>
#include <dai/evidence.h>
#include <dai/emalg.h>

using namespace dai;
using namespace std;
using namespace boost::python;

void (BP::*BPInit1)() = &BP::init;

void (JTree::*JTinit)() = (&JTree::init);
dai::Factor (JTree::*JTbelief)(const Var&) const = (&JTree::belief);

const std::vector<size_t>& (Gibbs::*gibbs_state) () const = (&Gibbs::state);

Real (EMAlg::*em_iterate)() = (&EMAlg::iterate);

size_t (Var::*const_label)() const = &Var::label;
size_t (Var::*const_states)() const = &Var::states;

const VarSet& (dai::Factor::*const_vars)() const = &dai::Factor::vars;
const dai::Prob& (dai::Factor::*factor_p_const)() const = &dai::Factor::p;

const dai::Factor & (dai::FactorGraph::*const_factor)(size_t) const = &dai::FactorGraph::factor;
void (dai::FactorGraph::*clamp_idx)(size_t, size_t, bool) = &dai::FactorGraph::clamp;
const dai::FactorGraph::Neighbors & (dai::FactorGraph::*fgVarNeighbors)(size_t) const = &dai::FactorGraph::nbV;
const dai::FactorGraph::Neighbors & (dai::FactorGraph::*fgFactorNeighbors)(size_t) const = &dai::FactorGraph::nbF;
const dai::FactorGraph::Neighbor & (dai::FactorGraph::*fgFactorNeighbor)(size_t, size_t) const = &dai::FactorGraph::nbF;

void (dai::FactorGraph::*backupFactors)(const VarSet&) = (&dai::FactorGraph::backupFactors);
void (dai::FactorGraph::*restoreFactors)(const VarSet&) = (&dai::FactorGraph::restoreFactors);
void (dai::FactorGraph::*restoreAllFactors)() = (&dai::FactorGraph::restoreFactors);

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(prob_normalize_overloads, normalize, 0, 1);

PropertySet& propertySet_set (PropertySet &ps, string key, string val) {
	return ps.Set(key,val);
}

string propertySet_get(PropertySet &ps, string key) {
	return ps.GetAs<std::string>(key);
}

void propertySet_read (PropertySet &ps, const char* filename) {
	ifstream propertyFile(filename);
	propertyFile >> ps;
	propertyFile.close();
}

void multifor_next ( multifor& mfor ) {
	mfor++;
}

void factor_setitem(dai::Factor& factor, size_t index, double value)
{
	if (index >= 0 && index < factor.states()) {
		factor[index] = value;
	} else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
	}
}

double factor_getitem(dai::Factor &factor, size_t index)
{
	if (index >= 0 && index < factor.states()) {
		return factor[index];
	} else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
		return -1;
	}
}

void prob_setitem(dai::Prob& prob, size_t index, double value) {
	if (index >= 0 && index < prob.size()) {
		prob[index] = value;
	} else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
	}
}

double prob_getitem(dai::Prob& prob, size_t index) {
	if (index >= 0 && index < prob.size()) {
		return prob[index];
	} else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
		return -1;
	}
}


dai::Factor vector_factor_getitem(vector<dai::Factor> &factors, size_t index)
{
	if (index >= 0 && index < factors.size()) {
		return factors[index];
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
		return -1;
	}
}

void vector_factor_setitem(vector<dai::Factor> &factors, size_t index, dai::Factor value) {
	if (index >= 0 && index < factors.size()) {
		factors[index] = value;
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
	}
}

void setVarLabel(Var &var, size_t newLabel) {
	var.label() = newLabel;
}

size_t stateOfVar (State &state, const Var &var) {
	return state(var);
}

void incrState(State &state) {
	state++;
}

EMAlg createEM(Evidence &evi, JTree &infAlg, string filename ) {
	ifstream msteps_file( filename.c_str() );
	EMAlg em(evi, infAlg, msteps_file);
	msteps_file.close();
	return em;
}

string fgToDotString (dai::FactorGraph &fg) {
	stringstream out;
	fg.printDot(out);
	return out.str();
}

void export_dai() {
	def("rnd_seed", rnd_seed);
	def("createEM", createEM);

	class_<vector<map<Var, size_t> > >("vector_map_var_sizet")
		.def(init<size_t>())
		.def(vector_indexing_suite<vector<map<Var, size_t> > >())
		;

	class_<Var>("Var", init<long, size_t>())
		.def(str(self))
		.def("label", const_label)
		.def("states", const_states)
		.def("setLabel", setVarLabel)
		;

	class_<vector<Var> >("vector_var")
		.def(str(self))
		.def(vector_indexing_suite<std::vector<Var> >())
		.def("push_back", &vector<Var>::push_back)
		.def("size", &vector<Var>::size)
		.def("clear", &vector<Var>::clear)
		.def("reserve", &vector<Var>::reserve)
		;

	class_<SmallSet<Var> >("smallSet_var")
		;

	class_<VarSet, bases<SmallSet<Var> > >("VarSet")
		.def(init<SmallSet<Var> >())
		.def(self |=  Var())
		.def(self & self)
		.def(self | self)
		.def("__iter__", boost::python::iterator<VarSet, return_internal_reference<> >())
		.def(str(self))
		.def("size", &VarSet::size)
		.def("calcState", &VarSet::calcState)
		.def("calcStates", &VarSet::calcStates)
		;

	class_<dai::Factor>("Factor", init<VarSet>())
		.def(init<VarSet, double>())
		.def(init<const dai::VarSet&, const dai::Prob&>())
		.def("__getitem__", &factor_getitem)
		.def("__setitem__", &factor_setitem)
		.def("__len__", &dai::Factor::states)
		.def(str(self))
		.def("vars", const_vars, return_internal_reference<>())
		.def("states", &dai::Factor::states)
		.def("p", factor_p_const, return_internal_reference<>())
		;

	class_<vector<dai::Factor> >("vector_factor")
		.def(init<size_t>())
		.def(str(self))
		.def("__iter__", boost::python::iterator<vector<dai::Factor>, return_internal_reference<> >())
		.def("__getitem__", &vector_factor_getitem)
		.def("__setitem__", &vector_factor_setitem)
		.def("push_back", &vector<dai::Factor>::push_back)
		.def("size", &vector<dai::Factor>::size)
		;

	class_<dai::Prob>("Prob", init<size_t, double>())
		.def("draw", &dai::Prob::draw)
		.def("normalize", &dai::Prob::normalize, prob_normalize_overloads())
		.def("__getitem__", prob_getitem)
		.def("__setitem__", prob_setitem)
		.def(str(self))
		;

	class_<dai::FactorGraph>("FactorGraph")
		.def(init<vector<dai::Factor> >())
		.def(str(self))
		.def("ReadFromFile", &dai::FactorGraph::ReadFromFile)
		.def("nrEdges", &dai::FactorGraph::nrEdges)
		.def("nrFactors", &dai::FactorGraph::nrFactors)
		.def("nrVars", &dai::FactorGraph::nrVars)
		.def("factor", const_factor, return_internal_reference<>())
		.def("factors", &dai::FactorGraph::factors, return_internal_reference<>())
		.def("nbF", fgFactorNeighbor, return_internal_reference<>())
		.def("nbF", fgFactorNeighbors, return_internal_reference<>())
		.def("var", &dai::FactorGraph::var, return_internal_reference<>())
		.def("vars", &dai::FactorGraph::vars, return_internal_reference<>())
		.def("findVar", &dai::FactorGraph::findVar)
		.def("nbV", fgVarNeighbors, return_internal_reference<>())
		.def("clamp", clamp_idx)
		.def("fgToDotString",&fgToDotString)
		.def("backupFactors", backupFactors)
		.def("restoreFactors", restoreFactors)
		.def("restoreAllFactors", restoreAllFactors)
		.def("restoreFactor", &dai::FactorGraph::restoreFactor)
		;

	class_<BipartiteGraph::Neighbor>("Neighbor")
		.def_readonly("node", &BipartiteGraph::Neighbor::node)
		.def_readonly("dual", &BipartiteGraph::Neighbor::dual)
		.def_readonly("iter", &BipartiteGraph::Neighbor::iter)
		;

	class_<BipartiteGraph::Neighbors>("Neighbors")
		.def(vector_indexing_suite<BipartiteGraph::Neighbors>())
		;

	class_<BP, bases<dai::FactorGraph> >("BP", init<dai::FactorGraph, PropertySet>())
		.def("init", BPInit1)
		.def("run", &BP::run)
		.def("beliefV", &BP::beliefV)
		.def("beliefF", &BP::beliefF)
		.def("getIterations", &BP::Iterations)
		;

	class_<JTree, bases<dai::FactorGraph> >("JTree", init<dai::FactorGraph, PropertySet>())
		.def("init", JTinit)
		.def("run", &JTree::run)
		.def("belief", JTbelief)
		.def("beliefV", &JTree::beliefV)
		;

	class_<Gibbs>("Gibbs", init<dai::FactorGraph, PropertySet>())
		.def("run", &Gibbs::run)
		.def("state", gibbs_state, return_internal_reference<>())
		.def("beliefV", &Gibbs::beliefV)
		;

	class_<EMAlg>("EMAlg", no_init)
		.def("hasSatisfiedTermConditions", &EMAlg::hasSatisfiedTermConditions)
		.def("iterate", em_iterate)
		.def("Iterations", &EMAlg::Iterations)
		;

	class_<PropertySet>("PropertySet")
		.def(str(self))
		.def("Set", propertySet_set, return_internal_reference<>())
		.def("Get", propertySet_get)
		.def("read", propertySet_read)
		;

	class_<Evidence>("Evidence")
		.def(init<vector<map<Var, size_t> > & >())
		;

	class_<State>("State", init<VarSet>())
		.def("stateOfVar", stateOfVar)
		.def("incr", incrState)
		;

	class_<multifor>("Multifor", init<vector<size_t> >())
		.def("valid", &multifor::valid)
		.def("index", &multifor::operator size_t)
		.def("__getitem__", &multifor::operator[])
		.def("next", multifor_next)
		;

	class_<Permute>("Permute", init<vector<Var> >())
		.def("convertLinearIndex", &Permute::convertLinearIndex)
		;
}

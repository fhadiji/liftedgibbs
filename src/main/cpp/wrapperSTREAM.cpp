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
 * libSTREAMWrapper.cpp
 *
 *  Created on: Jul 13, 2009
 *      Author: Fabian Hadiji
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

#include <STREAM/CBP/CBP.h>
#include <STREAM/CBP/CFactorGraph.h>
#include <STREAM/Gibbs/LiftedGibbs.h>
#include <STREAM/FactorGraph.h>

using namespace dai;
using namespace std;
using namespace stream;
using namespace boost::python;

void export_dai();


// STREAM functions
void (CBP::*CBPInit1)() = (&CBP::init);
void (CBP::*CBPInit2)(vector<size_t>&) = (&CBP::init);

void (CFactorGraph::*backupCFactors)(const VarSet&) = (&CFactorGraph::backupFactors);
void (CFactorGraph::*restoreCFactors)(const VarSet&) = (&CFactorGraph::restoreFactors);
void (CFactorGraph::*restoreAllCFactors)() = (&CFactorGraph::restoreFactors);

const dai::Prob& (CFactor::*const_p)() const = &CFactor::p;

const CFactor & (CFactorGraph::*const_cfactor)(size_t) const = &CFactorGraph::factor;
const CFactorGraph::Neighbors & (CFactorGraph::*cfgVarNeighbors)(size_t) const = &CFactorGraph::nbV;
const CFactorGraph::Neighbors & (CFactorGraph::*cfgFactorNeighbor)(size_t) const = &CFactorGraph::nbF;
const CFactorGraph::Neighbor & (CFactorGraph::*cfgFactorNeighbors)(size_t, size_t) const = &CFactorGraph::nbF;

CFactorGraph (CompressInterface::*compressCfg)(const CFactorGraph&) = &CompressInterface::compress;

void (LiftedGibbs::*createConfigMapping)() = &LiftedGibbs::createConfigMapping;
void (LiftedGibbs::*createConfigMappingCfg)(CFactorGraph& cfg) = &LiftedGibbs::createConfigMapping;

vector<map<size_t, int> > & cfactor_getCounts(CFactor &fac) {
	return fac.counts();
}

void createVarMapping (stream::FactorGraph<CFactor, Var>& fg, CompressInterface& compressAlg, CFactorGraph& liftedFg) {
	fg.createVarMapping(&compressAlg, liftedFg);
}

void createFactorMapping(stream::FactorGraph<CFactor, Var>& fg, CompressInterface& compressAlg, CFactorGraph& liftedFg) {
	fg.createFactorMapping(&compressAlg, liftedFg);
}

void cfactor_setCounts(CFactor &fac, vector<map<size_t, int> > & counts) {
	fac.counts() = counts;
}

void cfactor_setitem(CFactor& factor, size_t index, double value)
{
	if (index >= 0 && index < factor.states()) {
		factor[index] = value;
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
	}
}

double cfactor_getitem(CFactor &factor, size_t index)
{
	if (index >= 0 && index < factor.states()) {
		return factor[index];
	}
	else {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		boost::python::throw_error_already_set();
		return -1;
	}
}


BOOST_PYTHON_MODULE(libSTREAMWrapper)
{
	docstring_options doc_options(DEMO_DOCSTRING_SHOW_ALL);

	def("srand", srand);

	/* std containers */
	class_<map<long, int> >("map_long_int")
		.def(map_indexing_suite<map<long, int> >())
		;

	class_<map<Var, size_t> >("map_var_sizet")
		.def(str(self))
		.def(map_indexing_suite<map<Var, size_t> >())
		.def("clear", &map<Var, size_t>::clear)
		;

	class_<map<size_t,int> >("map_sizet_int")
		.def(map_indexing_suite<map<size_t, int> >())
		;

	class_<vector<size_t> >("vector_sizet")
		.def(init<size_t>())
		.def(vector_indexing_suite<vector<size_t> >())
		.def("push_back", &vector<size_t>::push_back)
		.def("size", &vector<size_t>::size)
		;

	class_<vector<int> >("vector_int")
		.def(init<size_t>())
		.def(vector_indexing_suite<vector<int> >())
		;

	class_<vector<map<size_t,int> > >("vector_map_sizet_int")
		.def(init<size_t>())
		.def(vector_indexing_suite<vector<map<size_t, int> > >())
		;

	class_<vector<CFactor> >("vector_cfactor")
		.def(init<size_t>())
		.def(vector_indexing_suite<std::vector<CFactor> >())
		.def("push_back", &vector<CFactor>::push_back)
		;

	class_<CFactor>("CFactor", init<VarSet>())
		.def(init<VarSet, double>())
		.def(init<const dai::VarSet&, const dai::Prob&>())
		.def("__getitem__", &cfactor_getitem)
		.def("__setitem__", &cfactor_setitem)
		.def("__len__", &CFactor::states)
		.def(str(self))
		.def("vars", &CFactor::vars, return_internal_reference<>())
		.def("states", &CFactor::states)
		.def("getCounts", &cfactor_getCounts, return_internal_reference<>())
		.def("setCounts", &cfactor_setCounts)
		.def("p", const_p, return_internal_reference<>())
		;

	class_<CFactorGraph>("CFactorGraph")
		.def(init<vector<CFactor> >())
		.def(init<vector<dai::Factor> >())
		.def(init<dai::FactorGraph>())
		.def(init<CFactorGraph>())
		.def(str(self))
		.def("readFromFile", &CFactorGraph::readFromFile)
		.def("nrEdges", &CFactorGraph::nrEdges)
		.def("nrFactors", &CFactorGraph::nrFactors)
		.def("nrVars", &CFactorGraph::nrVars)
		.def("factor", const_cfactor, return_internal_reference<>())
		.def("factors", &CFactorGraph::factors, return_internal_reference<>())
		.def("nbF", cfgFactorNeighbors, return_internal_reference<>())
		.def("nbF", cfgFactorNeighbor, return_internal_reference<>())
		.def("var", &CFactorGraph::var, return_internal_reference<>())
		.def("vars", &CFactorGraph::vars, return_internal_reference<>())
		.def("nbV", cfgVarNeighbors, return_internal_reference<>())
		.def("findVar", &CFactorGraph::findVar)
		.def("clamp", &CFactorGraph::clamp)
		/* CFactorGraph specific functions */
		.def("setSigma", &CFactorGraph::setSigma)
		.def("getSigma", &CFactorGraph::getSigma, return_internal_reference<>())
		.def("createCleanedGraph", &CFactorGraph::createCleanedGraph)
		.staticmethod("createCleanedGraph")
		.def("backupFactors", backupCFactors)
		.def("restoreFactors", restoreCFactors)
		.def("restoreAllFactors", restoreAllCFactors)
		.def("restoreFactor", &CFactorGraph::restoreFactor)
		.def("clampMlnFactor", &CFactorGraph::clampMlnFactor)
		;

	class_<CBP, bases<CFactorGraph> >("CBP", init<CFactorGraph, PropertySet>())
		.def(init<dai::FactorGraph, PropertySet>())
		.def("init", CBPInit1)
		.def("initColored", CBPInit2)
		.def("run", &CBP::run)
		.def("maxDiff",&CBP::maxDiff)
		.def("beliefV", &CBP::beliefV)
		.def("beliefF", &CBP::beliefF)
		.def("getIterations", &CBP::Iterations)
		/* CFactorGraph specific functions*/
		.def("reprV", &CBP::reprV)
		.def("reprF", &CBP::reprF)
		.def("createVarMapping", &CBP::createVarMapping)
		.def("createFacMapping", &CBP::createFacMapping)
		.def("clusterV", &CBP::clusterV)
		.def("clusterF", &CBP::clusterF)
		.def("setProperty", &CBP::setProperty)
		.def("getMessages", &CBP::getMessages)
		.def("getCompressIterations", &CBP::getCompressIterations)
		;

	class_<CompressInterface>("CompressInterface", no_init)
		.def("compress", compressCfg)
		;

	class_<PositionCompress, bases<CompressInterface> >("PositionCompress")
		.def("setCfg", &PositionCompress::setCfg)
		.def("init", &PositionCompress::init)
		.def("iterate", &PositionCompress::iterate)
		.def("createCFactorGraph", &PositionCompress::createCFactorGraph)
		.def("hasConverged", &PositionCompress::hasConverged)
		;

	

	class_<LiftedGibbs>("LiftedGibbs", init<const stream::FactorGraph<CFactor, Var>&, const PropertySet& >())
		.def("createConfigMapping", createConfigMapping)
		.def("createConfigMapping", createConfigMappingCfg)
		.def("run", &LiftedGibbs::run)
		.def("beliefV", &LiftedGibbs::beliefV)
		;

	class_<stream::FactorGraph<CFactor, Var> >("FactorGraphProb")
		.def(init<CFactorGraph>())
		.def("createVarMapping", &stream::FactorGraph<CFactor, Var>::createVarMapping)
		.def("createFactorMapping", &stream::FactorGraph<CFactor, Var>::createFactorMapping)
		.def("readFromFile", &stream::FactorGraph<CFactor, Var>::readFromFile)
		.def("reprV", &stream::FactorGraph<CFactor, Var>::reprV)
		.def("reprF", &stream::FactorGraph<CFactor, Var>::reprF)
		;

	export_dai();
}

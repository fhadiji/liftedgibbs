#!/usr/bin/python
#===============================================================================
#    Copyright (C) 2010
#    Babak Ahmadi [babak dot ahmadi at iais dot fraunhofer dot de]
#    Fabian Hadiji [fabian dot hadiji at iais dot fraunhofer dot de]
#    Kristian Kersting (coordination) [kristian dot kersting at iais dot fraunhofer dot de]
#
#    STREAM Project at
#        Fraunhofer IAIS, Sankt Augustin, Germany, and 
#        KDML, Unversity of Bonn, Germany 
#
#    This file is part of libSTREAM.
#
#    libSTREAM is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    libSTREAM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
# 
#    You should have received a copy of the GNU General Public License 
#    along with this program; if not, see <http://www.gnu.org/licenses/>.
#===============================================================================
'''
Created on Mar 26, 2010

@author: Fabian Hadiji
'''
import sys;
import time;
from libSTREAMWrapper import *;

if sys.platform == 'win32':
    defaultTimer = time.clock
else:
    defaultTimer = time.time

def usage():
    print 'usage: compareLiftedSampling.py fgFilename'

if __name__ == '__main__':
    if len(sys.argv) != 2:
        usage();
        exit();
    else:
        fgFilename = sys.argv[1];
        
    print "Comparing Lifted Sampling on file", fgFilename
    
    rnd_seed(int(time.time()))
    
    fg = FactorGraph();
    fg.ReadFromFile(fgFilename)
    
    bpOpts = PropertySet();
    bpOpts.read('property_set_bp');
    startTime = defaultTimer();
    bp = BP(fg, bpOpts);
    bp.init();
    bp.run();
    bpTime  = defaultTimer() - startTime;
    
    jTreeOpts = PropertySet();
    jTreeOpts.read('property_set_jtree');
    startTime = defaultTimer();
    jTree = JTree(fg, jTreeOpts);
    jTree.init();
    jTree.run();
    jTreeTime  = defaultTimer() - startTime;
    
    gibbsOpts = PropertySet();
    gibbsOpts.read('property_set_gibbs');
    startTime = defaultTimer();
    gibbs = Gibbs(fg, gibbsOpts);
    gibbs.run();
    gibbsTime = defaultTimer() - startTime;
    
    groundFg = FactorGraphProb();
    groundFg.readFromFile(fgFilename);
    startTime = defaultTimer();
    groundGibbs = LiftedGibbs(groundFg, gibbsOpts)
    groundGibbs.createConfigMapping();
    groundGibbs.run();
    groundGibbsTime = defaultTimer() - startTime;
    
    cfg = CFactorGraph();
    cfg.readFromFile(fgFilename)
    startTime = defaultTimer();
    compressAlg = PositionCompress();
    liftedFg = FactorGraphProb(compressAlg.compress(cfg));
    liftedFg.createVarMapping(compressAlg, cfg);
    liftedFg.createFactorMapping(compressAlg, cfg);
    liftedGibbs = LiftedGibbs(liftedFg, gibbsOpts)
    liftedGibbs.createConfigMapping(cfg);
    liftedGibbs.run();
    liftedGibbsTime = defaultTimer() - startTime;
    
    bpMaxNorm = -1;
    gibbsMaxNorm = -1;
    liftedGibbsMaxNorm = -1;
    groundGibbsMaxNorm = -1;
    print
    print "%-7s  %-8s  %-8s  %-8s  %-8s  %-8s" % ("","JTree", "BP", "DAIGibbs", "GGibbs", "LGibbs")
    for i in range(fg.nrVars()):
        print '%-7d  %f  %f  %f  %f  %f' % (i, jTree.beliefV(i)[1], bp.beliefV(i)[1], gibbs.beliefV(i)[1], groundGibbs.beliefV(i)[1] ,liftedGibbs.beliefV(liftedFg.reprV(i))[1]);
        
        if abs(jTree.beliefV(i)[1] - bp.beliefV(i)[1]) > bpMaxNorm:
            bpMaxNorm = abs(jTree.beliefV(i)[1] - bp.beliefV(i)[1]);
        if abs(jTree.beliefV(i)[1] - gibbs.beliefV(i)[1]) > gibbsMaxNorm:
            gibbsMaxNorm = abs(jTree.beliefV(i)[1] - gibbs.beliefV(i)[1]);
        if abs(jTree.beliefV(i)[1] - liftedGibbs.beliefV(liftedFg.reprV(i))[1]) > liftedGibbsMaxNorm:
            liftedGibbsMaxNorm = abs(jTree.beliefV(i)[1] - liftedGibbs.beliefV(liftedFg.reprV(i))[1]);
        if abs(jTree.beliefV(i)[1] - groundGibbs.beliefV(i)[1]) > groundGibbsMaxNorm:
            groundGibbsMaxNorm = abs(jTree.beliefV(i)[1] - groundGibbs.beliefV(i)[1]); 
        
    print "%-7s  %-8s  %f  %f  %f  %f" % ("MAXNORM","", bpMaxNorm, gibbsMaxNorm, groundGibbsMaxNorm, liftedGibbsMaxNorm)
    print "%-7s  %f  %f  %f  %f  %f" % ("RUNTIME", jTreeTime, bpTime, gibbsTime, groundGibbsTime, liftedGibbsTime)
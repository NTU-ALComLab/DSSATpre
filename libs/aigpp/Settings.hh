/**************************************************************
 *
 *       AIGPP Package // Settings.hh
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 111 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#ifndef AIGPP_SETTINGS_HH
#define AIGPP_SETTINGS_HH

#include <string>

namespace aigpp {

class Settings
{
   public:
    Settings();

    void load(const std::string& fileName);

    /* Node replacement */
    enum NodeReplacement
    {
        NODEREPLACEMENT_NONE,
        NODEREPLACEMENT_CONESIZE
    };

   private:
    NodeReplacement _NodeReplacement;

   public:
    const NodeReplacement& getNodeReplacement() const { return _NodeReplacement; };
    void                   setNodeReplacement(const NodeReplacement& t) { _NodeReplacement = t; };

    /* Quantifier Scheduling */
    enum QuantifierScheduling
    {
        QUANTIFIERSCHEDULING_NONE,
        QUANTIFIERSCHEDULING_REMAINING_NODES,
        QUANTIFIERSCHEDULING_MIN_DEPTH
    };

   private:
    QuantifierScheduling _QuantifierScheduling;

   public:
    const QuantifierScheduling& getQuantifierScheduling() const { return _QuantifierScheduling; };
    void                        setQuantifierScheduling(const QuantifierScheduling& t) { _QuantifierScheduling = t; };

    /* BDD Sweeping */
    enum BDDSweeping
    {
        BDDSWEEPING_NONE,
        BDDSWEEPING_COFACTOR
    };

   private:
    BDDSweeping _BDDSweeping;
    bool        _BDDSweepingFinalReorder;
    bool        _BDDSweepingInitializeOrder;

   public:
    const BDDSweeping& getBDDSweeping() const { return _BDDSweeping; };
    void               setBDDSweeping(const BDDSweeping& t) { _BDDSweeping = t; };

    const bool& getBDDSweepingFinalReorder() const { return _BDDSweepingFinalReorder; };
    void        setBDDSweepingFinalReorder(const bool& t) { _BDDSweepingFinalReorder = t; };

    const bool& getBDDSweepingInitializeOrder() const { return _BDDSweepingInitializeOrder; };
    void        setBDDSweepingInitializeOrder(const bool& t) { _BDDSweepingInitializeOrder = t; };

    /* Redundancy detection/removal */
    enum RedundancyRemoval
    {
        REDUNDANCYREMOVAL_NONE,
        REDUNDANCYREMOVAL_QUANTIFICATION
    };

   private:
    RedundancyRemoval _RedundancyRemoval;

   public:
    const RedundancyRemoval& getRedundancyRemoval() const { return _RedundancyRemoval; };
    void                     setRedundancyRemoval(const RedundancyRemoval& t) { _RedundancyRemoval = t; };

   public:
    /* simulation vector merging */
    bool getSimVectorMerging() const;
    void setSimVectorMerging(bool b);

    /*!
     * Garbage collection is activated when |DEAD| > r * ( |DEAD| + |ALIVE| ),
     * e.g. the set of all nodes consists of at least (r*100)% dead nodes. Default
     * is r=0.95
     */
    double getGCRatio() const;
    void   setGCRatio(double d);

    /* default = 10 */
    int  getSATResetInterval() const;
    void setSATResetInterval(int i);

    bool getDumpSimulation() const { return _dumpSimulation; };
    void setDumpSimulation(bool v) { _dumpSimulation = v; };

    bool getDecomposedQuantification() const { return _decomposedQuantification; };
    void setDecomposedQuantification(bool v) { _decomposedQuantification = v; };

   private:
    bool   _simVectorMerging;
    double _gcRatio;
    int    _satResetInterval;
    bool   _dumpSimulation;
    bool   _decomposedQuantification;

   public:
    /***** quantification stuff *****/
    bool quantifyUseScheduling() const { return _quantifyUseScheduling; }
    void setQuantifyUseScheduling(bool t) { _quantifyUseScheduling = t; }
    bool _quantifyUseScheduling;

    bool quantifyUseInterpolation() const { return _quantifyUseInterpolation; }
    void setQuantifyUseInterpolation(bool t) { _quantifyUseInterpolation = t; }
    bool _quantifyUseInterpolation;

    bool quantifyUseDualInterpolation() const { return _quantifyUseDualInterpolation; }
    void setQuantifyUseDualInterpolation(bool t) { _quantifyUseDualInterpolation = t; }
    bool _quantifyUseDualInterpolation;

    double quantifyInterpolationTimeout() const { return _quantifyInterpolationTimeout; }
    void   setQuantifyInterpolationTimeout(double t) { _quantifyInterpolationTimeout = t; }
    double _quantifyInterpolationTimeout;
};

}  // namespace aigpp

#include "Settings.icc"

#endif /* AIGPP_SETTINGS_HH */

/**************************************************************
 *
 *       AIGPP Package // Settings.cc
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

#include "Settings.hh"

/* std */
#include <fstream>

/* LRABS utilities */
#include <lrabsutil/ConfigMap.hh>

aigpp::Settings::Settings() :
    _NodeReplacement(NODEREPLACEMENT_CONESIZE),
    _QuantifierScheduling(QUANTIFIERSCHEDULING_NONE),
    _BDDSweeping(BDDSWEEPING_NONE),
    _BDDSweepingFinalReorder(true),
    _BDDSweepingInitializeOrder(false),
    _RedundancyRemoval(REDUNDANCYREMOVAL_NONE),
    _simVectorMerging(false),
    _gcRatio(0.5),
    _satResetInterval(10),
    _dumpSimulation(false),
    _decomposedQuantification(false),
    _quantifyUseScheduling(true),
    _quantifyUseInterpolation(true),
    _quantifyUseDualInterpolation(true),
    _quantifyInterpolationTimeout(5.0)
{}

void
aigpp::Settings::load(const std::string& fileName)
{
    std::ifstream configFile(fileName.c_str());
    // file not found -> return silently
    if (!configFile) {
        std::cout << "config file not found: " << fileName << std::endl;
        throw;
    }
    // file found -> load config map
    else {
        std::cout << "loading config file \"" << fileName << "\"" << std::endl;

        lrabs::ConfigMap m;
        m.readMap(configFile);
        m.writeMap(std::cout);

        std::set<std::string> validKeys;
        validKeys.insert("node_replacement");
        validKeys.insert("quantifier_scheduling");
        validKeys.insert("bdd_sweeping");
        validKeys.insert("bdd_sweeping_final_reorder");
        validKeys.insert("bdd_sweeping_initialize_order");
        validKeys.insert("redundancy_removal");
        validKeys.insert("simvector_merging");
        validKeys.insert("gc_ratio");
        validKeys.insert("sat_reset");
        validKeys.insert("dump_simulation");
        validKeys.insert("decomposed_quantification");

        validKeys.insert("quantify_use_scheduling");
        validKeys.insert("quantify_use_interpolation");
        validKeys.insert("quantify_use_dual_interpolation");
        validKeys.insert("quantify_interpolation_timeout");

        // some keys are invalid -> abort
        if (!m.allKeysAreValid(validKeys)) {
            std::cout << "ERROR: the specified config file \"" << fileName << "\" contains invalid keys!" << std::endl;

            throw;
        }

        // all keys are valid -> set properties
        else {
            {
                std::map<std::string, NodeReplacement> mapping;
                mapping["none"]      = NODEREPLACEMENT_NONE;
                mapping["cone_size"] = NODEREPLACEMENT_CONESIZE;
                m.getValue("node_replacement", mapping, _NodeReplacement);
            }

            {
                std::map<std::string, QuantifierScheduling> mapping;
                mapping["none"]            = QUANTIFIERSCHEDULING_NONE;
                mapping["remaining_nodes"] = QUANTIFIERSCHEDULING_REMAINING_NODES;
                mapping["min_depth"]       = QUANTIFIERSCHEDULING_MIN_DEPTH;
                m.getValue("quantifier_scheduling", mapping, _QuantifierScheduling);
            }

            {
                std::map<std::string, BDDSweeping> mapping;
                mapping["none"]     = BDDSWEEPING_NONE;
                mapping["cofactor"] = BDDSWEEPING_COFACTOR;
                m.getValue("bdd_sweeping", mapping, _BDDSweeping);
            }

            m.getValue("bdd_sweeping_final_reorder", _BDDSweepingFinalReorder);

            m.getValue("bdd_sweeping_initialize_order", _BDDSweepingInitializeOrder);

            {
                std::map<std::string, RedundancyRemoval> mapping;
                mapping["none"]           = REDUNDANCYREMOVAL_NONE;
                mapping["quantification"] = REDUNDANCYREMOVAL_QUANTIFICATION;
                m.getValue("redundancy_removal", mapping, _RedundancyRemoval);
            }

            m.getValue("simvector_merging", _simVectorMerging);

            m.getValue("gc_ratio", _gcRatio);

            m.getValue("sat_reset", _satResetInterval);

            m.getValue("dump_simulation", _dumpSimulation);

            m.getValue("decomposed_quantification", _decomposedQuantification);

            m.getValue("quantify_use_scheduling", _quantifyUseScheduling);
            m.getValue("quantify_use_interpolation", _quantifyUseInterpolation);
            m.getValue("quantify_use_dual_interpolation", _quantifyUseDualInterpolation);
            m.getValue("quantify_interpolation_timeout", _quantifyInterpolationTimeout);
        }
    }
}

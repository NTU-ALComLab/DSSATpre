/**************************************************************
 *
 *       AIGPP // Manager_export.cc
 *
 *       Copyright (C) 2006 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: 717 $
 *         $Author: pigorsch $
 *
 ***************************************************************/

#include "Manager.hh"

#include <cassert>
#include <iomanip>
#include <sstream>

#include <lrabsutil/String.hh>

#define _unused(x) ((void)(x))

void
aigpp::Manager::exportAIGER(std::ostream& os, aigpp::Manager::AIGERMode mode) const
{
    assert(mode == AIGER_ASCIIMode);
    _unused(mode);

    std::map<unsigned int, std::size_t> indexmap;

    /* FORMAT: aag MAXINDEX INPUTS LATCHES OUTPUTS GATES */
    os << "aag " << nodeCount() << " " << variableCount() << " 0 0 " << nodeCount() - variableCount() << '\n';

    for (std::size_t i = 0; i < variableCount(); ++i) {
        os << (i + 1) * 2 << '\n';
        indexmap[_variables[i]->index()] = 2 * (i + 1);
    }

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (!(n->isVar()) && n->refCount() > 0) {
            std::size_t newindex = 2 * (indexmap.size() + 1);
            indexmap[n->index()] = newindex;
        }
    }

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (!(n->isVar()) && n->refCount() > 0) {
            std::size_t p1idx = indexmap[n->parent1().node()->index()];
            if (n->parent1().isInverted() != n->parent1().node()->isNAND()) {
                ++p1idx;
            }

            std::size_t p2idx = indexmap[n->parent2().node()->index()];
            if (n->parent2().isInverted() != n->parent2().node()->isNAND()) {
                ++p2idx;
            }

            os << indexmap[n->index()] << " " << p1idx << " " << p2idx << '\n';
        }
    }

    for (std::size_t i = 0; i < variableCount(); ++i) {
        os << "i" << i << " " << _variableNames[i] << '\n';
    }
}

void
aigpp::Manager::exportRootsAIGER(std::ostream& os) const
{
    std::vector<aigpp::Edge> outputs;
    std::vector<std::string> outputNames;

    int count = 0;
    for (const auto& e : _extRefTable->_targets) {
        if (e._refCount != 0) {
            outputs.push_back((Edge)(e._target));
            outputNames.push_back("root" + lrabs::toString(count));
            ++count;
        }
    }

    exportAIGER(os, outputs, outputNames);
}

void
aigpp::Manager::exportAIGER(std::ostream& os, const std::vector<aigpp::EdgeRef>& outputs,
                            const std::vector<std::string>& outputNames, AIGERMode mode) const
{
    std::vector<aigpp::Edge> outputs2;
    outputs2.reserve(outputs.size());

    for (const auto& p : outputs) {
        outputs2.push_back((Edge)(p.getInternal()));
    }

    if (mode == AIGER_BinaryMode) {
        exportAIGERBinary(os, outputs2, outputNames);
    } else {
        exportAIGER(os, outputs2, outputNames);
    }
}

void
aigpp::Manager::exportAIGER(std::ostream& os, const std::vector<aigpp::Edge>& outputs,
                            const std::vector<std::string>& outputNames) const
{
    assert(outputs.size() == outputNames.size());
    assert(outputs.size() >= 1);

    /* node* -> aiger node index */
    std::map<const Node*, std::size_t> aindex;
    /* init with constant0 node */
    aindex[nullptr] = 0;

    std::vector<std::size_t> input_aindex;

    /* get shared support -> inputs */
    std::vector<int> inputs = sharedSupport(outputs);

    /* init state bits */
    for (int v : inputs) {
        assert(isVarIndexValid(v));
        Node* n = _variables[v];
        assert(n->isVar());

        assert(aindex.find(n) == aindex.end());

        std::size_t ai = 2 * aindex.size();
        aindex[n]      = ai;

        input_aindex.push_back(ai);
    }

    /* collect functions */
    std::stack<const Node*> pending;

    for (const auto& v : outputs) {
        pending.push(v.node());
    }

    std::vector<std::vector<std::size_t>> ands;

    while (!pending.empty()) {
        const Node* n = pending.top();

        if (aindex.find(n) != aindex.end()) {
            pending.pop();
            continue;
        }

        assert(!(n->isVar()));

        const Node* n1  = n->parent1().node();
        const auto  ai1 = aindex.find(n1);
        if (ai1 == aindex.cend()) {
            pending.push(n1);
            continue;
        }

        const Node* n2  = n->parent2().node();
        const auto  ai2 = aindex.find(n2);
        if (ai2 == aindex.cend()) {
            pending.push(n2);
            continue;
        }

        std::size_t ai = 2 * aindex.size();
        aindex[n]      = ai;

        std::vector<std::size_t> andv(3);
        andv[0] = ai;
        andv[1] = ai1->second ^ ((n->parent1().isInverted() != n1->isNAND()) ? 1 : 0);
        andv[2] = ai2->second ^ ((n->parent2().isInverted() != n2->isNAND()) ? 1 : 0);
        ands.push_back(andv);
    }

    os << "aag " << (aindex.size() - 1) << " " << input_aindex.size() << " 0 " << outputs.size() << " " << ands.size()
       << '\n';

    /* print inputs (inputs + LCs) */
    for (const std::size_t v : input_aindex) {
        os << v << '\n';
    }

    /* print outputs */
    for (const Edge& e : outputs) {
        std::size_t t = aindex[e.node()];
        if (e.isInverted()) {
            t = t ^ 1u;
        }
        if (e.node() != nullptr) {
            if (e.node()->isNAND()) {
                t = t ^ 1u;
            }
        }

        os << t << '\n';
    }

    /* print ands */
    for (std::vector<std::vector<std::size_t>>::const_iterator p = ands.begin(); p != ands.end(); ++p) {
        assert(p->size() == 3);
        os << (*p)[0] << " " << (*p)[1] << " " << (*p)[2] << '\n';
    }

    /* print input labels */
    std::size_t index = 0;
    for (std::vector<int>::const_iterator v = inputs.begin(); v != inputs.end(); ++v) {
        os << "i" << index << " " << *_variables[*v] << '\n';
        ++index;
    }

    /* print output labels */
    index = 0;
    for (const auto& v : outputNames) {
        os << "o" << index << " " << v << '\n';
        ++index;
    }
}

/*
 * Export AIGER file in binary format:
 *
 * M: maximum variable index
 * I: number of inputs
 * L: number of latches (always 0)
 * O: number of outputs
 * A: number of AND gates
 *
 * Variable ordering:
 * inputs: I, And gates: I+A
 * unsigned literals: input: 2*I, AND gates: 2*(I+A)
 *
 * Binary encoded file:
 * 1. header: aig M I L O A
 * 2. list of outputs (ASCII format)
 * 3. AND gates (binary encoded)
 * 4. symbol table (ASCII format)
 * 5. comment section (ASCII format)
 */
void
aigpp::Manager::exportAIGERBinary(std::ostream& os, const std::vector<Edge>& outputs,
                                  const std::vector<std::string>& outputNames) const
{
    assert(outputs.size() == outputNames.size());
    assert(outputs.size() >= 1);

    /* node* -> aiger node index */
    std::map<const Node*, std::size_t> aindex;
    /* init with constant0 node */
    aindex[nullptr] = 0;

    std::vector<std::size_t> input_aindex;

    /* get shared support -> inputs */
    std::vector<int> inputs = sharedSupport(outputs);

    /* init state bits */
    for (const int v : inputs) {
        assert(isVarIndexValid(v));
        Node* n = _variables[v];
        assert(n->isVar());

        assert(aindex.find(n) == aindex.end());

        std::size_t ai = 2 * aindex.size();  // unsigned inputs: 2*1, ..., 2*I
        aindex[n]      = ai;

        input_aindex.push_back(ai);
    }

    /* collect functions */
    std::stack<const Node*> pending;

    for (const auto& v : outputs) {
        pending.push(v.node());
    }

    std::vector<std::vector<std::size_t>> ands;

    while (!pending.empty()) {
        const Node* n = pending.top();

        if (aindex.find(n) != aindex.end()) {
            pending.pop();
            continue;
        }

        assert(!(n->isVar()));

        const Node*                                        n1  = n->parent1().node();
        std::map<const Node*, std::size_t>::const_iterator ai1 = aindex.find(n1);
        if (ai1 == aindex.end()) {
            pending.push(n1);
            continue;
        }

        const Node*                                        n2  = n->parent2().node();
        std::map<const Node*, std::size_t>::const_iterator ai2 = aindex.find(n2);
        if (ai2 == aindex.end()) {
            pending.push(n2);
            continue;
        }

        std::size_t lhs = 2 * (input_aindex.size() + ands.size() + 1);  // unsigned AND literals: 2*(1+I), ..., 2*(I+A)
        aindex[n]       = lhs;

        // AND gate in ASCII: lhs (gateNum) rhs0 (input 1) rhs1 (input2)
        // AND gate in binary: delta0 (= lhs - rhs0) delta1 (= rhs0 - rhs1)

        size_t rhs0 = ai1->second ^ ((n->parent1().isInverted() != n1->isNAND()) ? 1 : 0);
        size_t rhs1 = ai2->second ^ ((n->parent2().isInverted() != n2->isNAND()) ? 1 : 0);

        if (rhs0 < rhs1) {
            const size_t temp = rhs0;
            rhs0              = rhs1;
            rhs1              = temp;
        }

        assert(lhs > rhs0);
        assert(rhs0 >= rhs1);

        std::vector<std::size_t> andv(2);
        andv[0] = lhs - rhs0;
        andv[1] = rhs0 - rhs1;
        ands.push_back(andv);
    }

    // header in ASCII
    os << "aig " << (aindex.size() - 1) << " " << input_aindex.size() << " 0 " << outputs.size() << " " << ands.size()
       << '\n';

    /* print outputs - ASCII*/
    for (std::size_t i = 0; i != outputs.size(); ++i) {
        Edge        e = outputs[i];
        std::size_t t = aindex[e.node()];
        if (e.isInverted()) {
            t = t ^ 1;
        }
        if (e.node() != nullptr) {
            if (e.node()->isNAND()) {
                t = t ^ 1;
            }
        }

        os << t << '\n';
    }

    /* print ands - binary*/
    for (const auto& delta : ands) {
        assert(delta.size() == 2);
        writeBinary(os, delta[0]);
        writeBinary(os, delta[1]);
    }

    /* print input labels - ASCII*/
    std::size_t index = 0;
    for (int i : inputs) {
        os << "i" << index << " " << *_variables[i] << '\n';
        ++index;
    }

    /* print output labels - ASCII */
    index = 0;
    for (const auto& n : outputNames) {
        os << "o" << index << " " << n << '\n';
        ++index;
    }
}

/*
 * \brief Writes value in binary format into ofstream .
 *
 * The binary encoding of size_t delta in AIGER is the sequence of i bites
 * b0,..., bi: 1w0, 1w1, ..., 1w(i-1), 0wi.
 */
void
aigpp::Manager::writeBinary(std::ostream& os, size_t delta) const
{
    unsigned char ch;

    while (delta & ~0x7F) {          // delta is not jet completely binary encoded - more
                                     // than 7 bits left
        ch = (delta & 0x7f) | 0x80;  // set bit 1 - 7 from delta and set MSB (bit 8) to 1
        os << ch;
        delta >>= 7;  // shift 7 bits
    }
    ch = (delta & ~0x80);  // less or equal than 7 bits left, MSB is set to 0
    os.write(((char*)&ch),
             sizeof(ch));  // write to file and moves file pointer ahead sizeof(ch) bytes
}

void
aigpp::Manager::exportSimulation(std::ostream& os) const
{
    assert(simCreation() == true);

    os << variableCount() << " " << (SimVector::BinCount * SimVector::BinSize) / 8 << '\n';

    for (std::size_t v = 0; v < variableCount(); ++v) {
        os << _variables[v]->sim() << '\n';
    }
}

void
aigpp::Manager::exportSeparationMatrix(std::ostream& os, int simvalues)
{
    assert(simCreation());

    /*
     * export a binary 2d-matrix: the rows correspond to node pairs, the columns
     * correspond to simulation values. cell(r,c)=1 iff  pair #r is separated by
     * simulation value #c
     *
     * an example of the output format for 5 pairs and 3 simulation values:
     * 5 3
     * 0 0 1
     * 1 1 0
     * 1 0 0
     * 1 0 1
     * 0 1 0
     */

    /* force insertion of pending new counter examples */
    updateSimTable(/*force=*/true);

    std::size_t pairs = (_nodeCount * (_nodeCount - 1)) / 2;
    if (simvalues <= 0 || simvalues > SIMVECSIZE) {
        simvalues = SIMVECSIZE;
    }

    os << pairs << " " << simvalues << std::endl;

    SimVector xorVec;
    for (Node* n1 = _nodes; n1 != nullptr; n1 = n1->next()) {
        for (Node* n2 = n1->next(); n2 != nullptr; n2 = n2->next()) {
            xorVec = SimVector::Xor(n1->sim(), n2->sim());

            for (int i = 0; i != simvalues; ++i) {
                if (xorVec.getBit(i)) {
                    os << "1 ";
                } else {
                    os << "0 ";
                }
            }

            os << std::endl;
        }
    }
}

void
aigpp::Manager::exportBLIF(std::ostream& os) const
{
    os << ".model aigdump" << std::endl;

    os << "\n# inputs" << std::endl;

    for (Node* n : _variables) {
        os << ".inputs " << *n << std::endl;
    }

    os << "\n# outputs" << std::endl;
    bool needConst    = false;
    int  currentIndex = 0;
    for (const auto& e : _extRefTable->_targets) {
        if (e._refCount != 0) {
            if (e._target.isConstant()) {
                needConst = true;
            }

            os << ".outputs "
               << "out" << currentIndex << std::endl;
            ++currentIndex;
        }
    }

    os << "\n# input mapping" << std::endl;

    for (Node* n : _variables) {
        os << ".names " << *n << " " << n << std::endl << "1 1" << std::endl;
    }

    if (needConst) {
        std::cout << "\n# constant" << std::endl;
        os << ".names const0" << std::endl;
    }

    os << "\n# internal signals" << std::endl;

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (n->refCount() == 0) {
            continue;
        }

        if (n->isVar()) {
            continue;
        }

        os << ".names " << n->parent1().node() << " " << n->parent2().node() << " " << n << std::endl;

        os << (n->parent1().isInverted() ? "0" : "1") << (n->parent2().isInverted() ? "0" : "1") << " "
           << (n->isNAND() ? "0" : "1") << std::endl;
    }

    os << "\n# output mapping" << std::endl;

    currentIndex = 0;
    for (const auto& e : _extRefTable->_targets) {
        if (e._refCount != 0) {
            if (e._target.isConstant()) {
                os << ".names const0 out" << currentIndex << std::endl;
            } else {
                os << ".names " << e._target.node() << " out" << currentIndex << std::endl;
            }

            if (e._target.isInverted()) {
                os << "0 1" << std::endl;
            } else {
                os << "1 1" << std::endl;
            }

            ++currentIndex;
        }
    }

    os << "\n.end" << std::endl;
}

void
aigpp::Manager::dumpGraphViz(std::ostream& os, bool displayNodeIndexes) const
{
    std::vector<Edge> roots;

    for (std::vector<ExtRefTable::Item>::const_iterator e = _extRefTable->_targets.begin();
         e != _extRefTable->_targets.end(); ++e) {
        if (e->_refCount != 0 && !(e->_target.isConstant())) {
            roots.push_back(e->_target);
        }
    }

    std::stack<Node*> pending;
    std::set<Node*>   processed;

    os << "digraph G {" << std::endl;

    for (Node* n = _nodes; n != nullptr; n = n->next()) {
        if (n->refCount() == 0) {
            continue;
        }

        pending.push(n);

        while (!pending.empty()) {
            Node* p = pending.top();
            pending.pop();

            if (processed.find(p) != processed.end()) {
                continue;
            } else if (p->isVar()) {
                if (!displayNodeIndexes) {
                    os << "n" << p << " [label=\"" << *p << "\" shape=box];" << std::endl;
                } else {
                    os << "n" << p << " [label=\"" << p->index() << ": " << *p << "\" shape=box];" << std::endl;
                }

                processed.insert(p);

                continue;
            } else {
                Node* p1 = p->parent1().node();
                Node* p2 = p->parent2().node();

                if (!displayNodeIndexes) {
                    os << "n" << p << " [label=\"\" "
                       << "shape=circle" << (p->isNAND() ? " style=dotted" : "") << "];" << std::endl
                       << "n" << p << " -> n" << p1 << (p->parent1().isInverted() ? "[style=dotted]" : "") << ";"
                       << std::endl
                       << "n" << p << " -> n" << p2 << (p->parent2().isInverted() ? "[style=dotted]" : "") << ";"
                       << std::endl;
                } else {
                    os << "n" << p << " [label=\"" << p->index() << "\" "
                       << "shape=circle" << (p->isNAND() ? " style=dotted" : "") << "];" << std::endl
                       << "n" << p << " -> n" << p1 << (p->parent1().isInverted() ? "[style=dotted]" : "") << ";"
                       << std::endl
                       << "n" << p << " -> n" << p2 << (p->parent2().isInverted() ? "[style=dotted]" : "") << ";"
                       << std::endl;
                }

                processed.insert(p);

                pending.push(p1);
                pending.push(p2);
            }
        }
    }

    os << "subgraph roots {" << std::endl;

    for (std::vector<aigpp::Edge>::const_iterator r = roots.begin(); r != roots.end(); ++r) {
        os << "root" << (r - roots.begin()) << " [shape=\"house\"];" << std::endl;
        os << "root" << (r - roots.begin()) << " -> n" << (r->isConstant() ? nullptr : r->node())
           << (r->isInverted() ? "[style=dotted]" : "") << ";" << std::endl;
    }

    os << "}" << std::endl;

    os << "}" << std::endl;
}

void
aigpp::Manager::dumpGraphViz(std::ostream& os, const std::vector<aigpp::EdgeRef>& roots, bool displayNodeIndexes) const
{
    std::vector<InternalEdgeRef> roots2 = EdgeRef::toInternal(roots);
    dumpGraphViz(os, roots2, displayNodeIndexes);
}

void
aigpp::Manager::dumpGraphViz(std::ostream& os, const std::vector<aigpp::InternalEdgeRef>& roots,
                             bool displayNodeIndexes) const
{
    std::stack<Node*> pending;
    std::set<Node*>   processed;

    os << "digraph G {" << std::endl;

    bool needConst = false;
    for (const auto& r : roots) {
        if (r.isConstant()) {
            needConst = true;
        } else {
            pending.push(r.node());
        }
    }

    if (needConst) {
        os << "n0 [label=\"false\" shape=box];" << std::endl;
    }

    while (!pending.empty()) {
        Node* p = pending.top();
        pending.pop();

        if (processed.find(p) != processed.end()) {
            continue;
        } else if (p->isVar()) {
            if (!displayNodeIndexes) {
                os << "n" << p << " [label=\"" << *p << "\" shape=box];" << std::endl;
            } else {
                os << "n" << p << " [label=\"" << p->index() << ": " << *p << "\" shape=box];" << std::endl;
            }
            processed.insert(p);

            continue;
        } else {
            Node* p1 = p->parent1().node();
            Node* p2 = p->parent2().node();

            if (!displayNodeIndexes) {
                os << "n" << p << " [label=\"\" "
                   << "shape=circle" << (p->isNAND() ? " style=dotted" : "") << "];" << std::endl
                   << "n" << p << " -> n" << p1 << (p->parent1().isInverted() ? "[style=dotted]" : "") << ";"
                   << std::endl
                   << "n" << p << " -> n" << p2 << (p->parent2().isInverted() ? "[style=dotted]" : "") << ";"
                   << std::endl;
            } else {
                os << "n" << p << " [label=\"" << p->index() << "\" "
                   << "shape=circle" << (p->isNAND() ? " style=dotted" : "") << "];" << std::endl
                   << "n" << p << " -> n" << p1 << (p->parent1().isInverted() ? "[style=dotted]" : "") << ";"
                   << std::endl
                   << "n" << p << " -> n" << p2 << (p->parent2().isInverted() ? "[style=dotted]" : "") << ";"
                   << std::endl;
            }

            processed.insert(p);

            pending.push(p1);
            pending.push(p2);
        }
    }

    os << "subgraph roots {" << std::endl;

    for (auto r = roots.begin(); r != roots.end(); ++r) {
        os << "root" << (r - roots.begin()) << " [shape=\"house\"];" << std::endl;
        os << "root" << (r - roots.begin()) << " -> n" << (r->isConstant() ? nullptr : r->node())
           << (r->isInverted() ? "[style=dotted]" : "") << ";" << std::endl;
    }

    os << "}" << std::endl;

    os << "}" << std::endl;
}

void
aigpp::Manager::dumpGraphViz(std::ostream& os, const aigpp::EdgeRef& root, bool displayNodeIndexes) const
{
    dumpGraphViz(os, root.getInternal(), displayNodeIndexes);
}

void
aigpp::Manager::dumpGraphViz(std::ostream& os, const aigpp::InternalEdgeRef& root, bool displayNodeIndexes) const
{
    std::vector<InternalEdgeRef> roots(1, root);
    dumpGraphViz(os, roots, displayNodeIndexes);
}

/*
<?xml version="1.0" encoding="UTF-8"?>
<graphml xmlns="http://graphml.graphdrawing.org/xmlns"
    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
    xsi:schemaLocation="http://graphml.graphdrawing.org/xmlns
     http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd">
  <graph id="G" edgedefault="undirected">
    <node id="n0"/>
    <node id="n1"/>
    <node id="n2"/>
    <node id="n3"/>
    <node id="n4"/>
    <node id="n5"/>
    <node id="n6"/>
    <node id="n7"/>
    <node id="n8"/>
    <node id="n9"/>
    <node id="n10"/>
    <edge source="n0" target="n2"/>
    <edge source="n1" target="n2"/>
    <edge source="n2" target="n3"/>
    <edge source="n3" target="n5"/>
    <edge source="n3" target="n4"/>
    <edge source="n4" target="n6"/>
    <edge source="n6" target="n5"/>
    <edge source="n5" target="n7"/>
    <edge source="n6" target="n8"/>
    <edge source="n8" target="n7"/>
    <edge source="n8" target="n9"/>
    <edge source="n8" target="n10"/>
  </graph>
</graphml>
*/

void
aigpp::Manager::exportGraphML(std::ostream& os, const std::vector<aigpp::EdgeRef>& roots) const
{
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
       << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n"
       << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
       << "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n"
       << "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n";

    os << "<graph id=\"AIG\" edgedefault=\"directed\">\n";

    std::stack<Node*> pending;
    std::set<Node*>   processed;

    bool needConst = false;
    for (const auto& r : roots) {
        if (r.isConstant()) {
            needConst = true;
        } else {
            pending.push(r.getInternal().node());
        }
    }

    if (needConst) {
        os << "<node id=\"n0\"/>\n";
    }

    while (!pending.empty()) {
        Node* p = pending.top();
        pending.pop();

        if (processed.find(p) != processed.end()) {
            continue;
        } else if (p->isVar()) {
            os << "<node id=\"n" << p << "\"/>\n";
            processed.insert(p);

            continue;
        } else {
            Node* p1 = p->parent1().node();
            Node* p2 = p->parent2().node();

            os << "<node id=\"n" << p << "\"/>\n";
            os << "<edge source=\"n" << p << "\" target=\"n" << p1 << "\"/>\n";
            os << "<edge source=\"n" << p << "\" target=\"n" << p2 << "\"/>\n";

            processed.insert(p);

            pending.push(p1);
            pending.push(p2);
        }
    }

    for (auto r = roots.begin(); r != roots.end(); ++r) {
        os << "<node id=\"root" << (r - roots.begin()) << "\"/>\n";
        os << "<edge source=\"root" << (r - roots.begin()) << "\" target=\"n"
           << (r->isConstant() ? nullptr : r->getInternal().node()) << "\"/>\n";
    }

    os << "</graph>\n"
       << "</graphml>\n";
}

std::string
interpolateColor(int v, int maxv)
{
    int minr = 0x0f;
    int ming = 0x97;
    int minb = 0x59;

    int maxr = 0x2c;
    int maxg = 0x75;
    int maxb = 0xbf;

    double ar = double(maxr - minr) / double(maxv - 1);
    double ag = double(maxg - ming) / double(maxv - 1);
    double ab = double(maxb - minb) / double(maxv - 1);

    double br = minr - ar;
    double bg = ming - ag;
    double bb = minb - ab;

    auto r = static_cast<int>(static_cast<double>(ar) * v + static_cast<double>(br));
    auto g = static_cast<int>(static_cast<double>(ag) * v + static_cast<double>(bg));
    auto b = static_cast<int>(static_cast<double>(ab) * v + static_cast<double>(bb));

    std::ostringstream oss;
    oss << "#" << std::hex << std::setfill('0') << std::setw(2) << std::uppercase << r << std::hex << std::setfill('0')
        << std::setw(2) << std::uppercase << g << std::hex << std::setfill('0') << std::setw(2) << std::uppercase << b;

    return oss.str();
}

void
aigpp::Manager::exportGML(std::ostream& os, const std::vector<aigpp::EdgeRef>& roots) const
{
    bool dark = false;

    os << "graph [\n"
       << "id 0\n"
       << "version 0\n"
       << "directed 1\n"
       << "hierarchical 1\n";

    std::stack<Node*> pending;
    std::set<Node*>   processed;

    /* compute refcount */
    std::map<Node*, int> refcount;
    for (const auto& r : roots) {
        if (!(r.isConstant())) {
            pending.push(r.getInternal().node());
        }
    }

    while (!pending.empty()) {
        Node* p = pending.top();

        if (refcount.find(p) != refcount.end()) {
            pending.pop();
        } else if (p->isVar()) {
            refcount[p] = 0;
            pending.pop();
        } else if (refcount.find(p->parent1().node()) == refcount.end()) {
            pending.push(p->parent1().node());
        } else if (refcount.find(p->parent2().node()) == refcount.end()) {
            pending.push(p->parent2().node());
        } else {
            Node* p1 = p->parent1().node();
            Node* p2 = p->parent2().node();

            ++refcount[p1];
            ++refcount[p2];

            refcount[p] = 0;

            pending.pop();
        }
    }

    for (const auto& r : roots) {
        ++refcount[r.getInternal().node()];
    }
    refcount.erase(nullptr);

    int maxrefcount = -1;
    for (const auto p : refcount) {
        if (p.second > maxrefcount) {
            maxrefcount = p.second;
        }
    }

    bool needConst = false;
    for (const auto& r : roots) {
        if (r.isConstant()) {
            needConst = true;
        } else {
            pending.push(r.getInternal().node());
        }
    }

    std::map<Node*, int> ids;
    int                  nextid = 0;

    double minSize = 30;
    double maxSize = 120;

    if (needConst) {
        int id       = nextid++;
        ids[nullptr] = id;
        os << "node [\n"
           << "id " << id << '\n'
           << "label \"F\"\n"
           << "graphics [\n"
           << "type \"oval\"\n"
           << "hasOutline   0\n";

        if (dark) {
            os << "outline \"#FFFFFF\"\n";
        }

        os << "]\n"
           << "]\n";
    }

    while (!pending.empty()) {
        Node* p = pending.top();

        if (processed.find(p) != processed.end()) {
            pending.pop();
        } else if (p->isVar()) {
            double s  = minSize + ((maxSize - minSize) * (refcount[p])) / maxrefcount;
            int    id = nextid++;
            ids[p]    = id;
            os << "node [\n"
               << "id " << id
               << '\n'
               //               << "label \"" << *p << "\"\n"
               << "graphics [\n"
               << "type \"oval\"\n"
               << "w " << s << '\n'
               << "h " << s << '\n'
               << "hasOutline   0\n"
               << "fill \"" << interpolateColor(refcount[p], maxrefcount) << "\"\n";

            if (dark) {
                os << "outline \"#FFFFFF\"\n";
            }

            os << "]\n"
               << "]\n";

            processed.insert(p);
            pending.pop();
        } else if (processed.find(p->parent1().node()) == processed.end()) {
            pending.push(p->parent1().node());
        } else if (processed.find(p->parent2().node()) == processed.end()) {
            pending.push(p->parent2().node());
        } else {
            Node* p1 = p->parent1().node();
            Node* p2 = p->parent2().node();

            assert(ids.find(p->parent1().node()) != ids.end());
            assert(ids.find(p->parent2().node()) != ids.end());

            double s = minSize + ((maxSize - minSize) * (refcount[p])) / maxrefcount;

            int id1 = ids[p1];
            int id2 = ids[p2];
            int id  = nextid++;
            ids[p]  = id;
            os << "node [\n"
               << "id " << id << '\n'
               << "graphics [\n"
               << "type \"oval\"\n"
               << "w " << s << '\n'
               << "h " << s << '\n'
               << "hasOutline   0\n"
               << "fill \"" << interpolateColor(refcount[p], maxrefcount) << "\"\n";

            if (dark) {
                os << "outline \"#FFFFFF\"\n";
            }

            os << "]\n"
               << "]\n";

            os << "edge [\n"
               << "source " << id << '\n'
               << "target " << id1 << '\n'
               << "graphics [\n"
               << "targetArrow \"standard\"\n"
               << "fill \"#FCD605\"\n";

            os << "]\n"
               << "]\n";

            os << "edge [\n"
               << "source " << id << '\n'
               << "target " << id2 << '\n'
               << "graphics [\n"
               << "targetArrow \"standard\"\n"
               << "fill \"#FCD605\"\n";

            os << "]\n"
               << "]\n";

            processed.insert(p);
            pending.pop();
        }
    }

    /*
    for( std::vector<aigpp::EdgeRef>::const_iterator r = roots.begin();
         r != roots.end(); ++r )
    {
        int id = nextid++;

        os << "node [\n"
           << "id " << id << '\n'
           << "label \"root" << ( r - roots.begin() ) << "\"\n"
           << "graphics [\n"
           << "type \"rectangle\"\n"
           << "fill \"#8888FF\"\n";

        if( dark )
        {
            os << "outline \"#FFFFFF\"\n";
        }

        os << "]\n"
           << "]\n";

        int id_other = ids[r->getInternal().node()];

        os << "edge [\n"
           << "source " << id << '\n'
           << "target " << id_other << '\n'
           << "graphics [\n"
           << "targetArrow \"standard\"\n";

        if( dark )
        {
            os << "fill \"#FFFFFF\"\n";
        }

       if( r->getInternal().isInverted() )
       {
           os << "style \"dashed_dotted\"\n";
       }

        os << "]\n"
           << "]\n";
    }
    */

    os << "]\n";
}

std::vector<std::vector<int>>
aigpp::Manager::createCNF(const aigpp::EdgeRef& root, int& rootLiteral) const
{
    using Clause = std::vector<int>;
    using MapN2T = std::map<Node*, int>;

    std::vector<Clause> clauses;

    if (root.isConstant()) {
        rootLiteral = 0;

        if (structurallyFalse(root)) {
            /* add an empty clause */
            clauses.emplace_back();
        }

        return clauses;
    }

    std::stack<Node*> pending;
    MapN2T            reversemapping;

    if (root.isVariable()) {
        /* root is (possibly negated) variable -> satisfiable! */
        rootLiteral = 0;

        /* return an empty clause set */
        return clauses;
    }

    pending.push(root.getInternal().node());
    MapN2T::const_iterator p1, p2;

    while (!pending.empty()) {
        Node* n = pending.top();

        /* node already processed */
        if (reversemapping.find(n) != reversemapping.end()) {
            pending.pop();
        }
        /* non-processed variable */
        else if (n->isVar()) {
            reversemapping[n] = 2 * (int)reversemapping.size();
            pending.pop();
        }
        /* parent1 not processed yet */
        else if ((p1 = reversemapping.find(n->parent1().node())) == reversemapping.end()) {
            pending.push(n->parent1().node());
        }
        /* parent2 not processed yet */
        else if ((p2 = reversemapping.find(n->parent2().node())) == reversemapping.end()) {
            pending.push(n->parent2().node());
        }
        /* both parents processed */
        else {
            pending.pop();

            /* get literals */
            int tseitinLit = 2 * (int)reversemapping.size();

            int lit1   = p1->second ^ (n->parent1().isInverted() ? 1 : 0);
            int lit2   = p2->second ^ (n->parent2().isInverted() ? 1 : 0);
            int litAnd = tseitinLit ^ (n->isNAND() ? 1 : 0);

            /* update mapping(s) */
            reversemapping[n] = tseitinLit;

            /* create clauses */

            /*  x = y * z <-> ( !x + y ) * ( !x + z ) * ( x + !y + !z ) */
            Clause clause2(2), clause3(3);

            clause2[0] = litAnd ^ 1;
            clause2[1] = lit1;
            clauses.push_back(clause2);

            /*clause2[0] = litAnd ^ 1;*/
            clause2[1] = lit2;
            clauses.push_back(clause2);

            clause3[0] = litAnd;
            clause3[1] = lit1 ^ 1;
            clause3[2] = lit2 ^ 1;
            clauses.push_back(clause3);
        }
    }

    MapN2T::const_iterator proot = reversemapping.find(root.getInternal().node());
    assert(proot != reversemapping.end());

    rootLiteral = proot->second ^ (root.isInverted() ? 1 : 0);

    return clauses;
}

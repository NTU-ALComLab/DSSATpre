#include "EquivalenceChecker.hh"

#include <cassert>

void
aigpp::EquivalenceChecker::setOriginalRoots(const std::vector<aigpp::EdgeRef>& roots)
{
    assert(roots_.empty());

    man_.toggleAutomaticFRAIGing(false);

    /* look for original manager */
    Manager* origManager = getManager(roots);

    /* create variables */
    copyVariables(origManager);

    /* copy roots */
    std::map<Node*, InternalEdgeRef> mapping;
    for (const auto& r : roots) {
        const InternalEdgeRef& e    = r.getInternal();
        InternalEdgeRef        newe = copyCone(e, mapping);
        roots_.push_back(man_.DebugGetExternal(newe));
    }
}

void
aigpp::EquivalenceChecker::setOriginalRoots(const std::vector<aigpp::InternalEdgeRef>& roots)
{
    assert(roots_.empty());

    man_.toggleAutomaticFRAIGing(false);

    /* look for original manager */
    Manager* origManager = getManager(roots);

    /* create variables */
    copyVariables(origManager);

    /* copy roots */
    std::map<Node*, InternalEdgeRef> mapping;
    for (const auto& r : roots) {
        InternalEdgeRef newe = copyCone(r, mapping);
        roots_.push_back(man_.DebugGetExternal(newe));
    }
}

bool
aigpp::EquivalenceChecker::checkNewRoots(const std::vector<aigpp::EdgeRef>& roots)
{
    assert(roots.size() == roots_.size());
    /* look for original manager */
    Manager* origManager = getManager(roots);

    /* create variables */
    copyVariables(origManager);

    /* copy roots */
    std::vector<EdgeRef>             newroots;
    std::map<Node*, InternalEdgeRef> mapping;
    for (const auto& r : roots) {
        const InternalEdgeRef& e    = r.getInternal();
        InternalEdgeRef        newe = copyCone(e, mapping);
        newroots.push_back(man_.DebugGetExternal(newe));
    }

    /* check roots */
    bool diff = false;
    for (unsigned int i = 0; i != roots_.size(); ++i) {
        if (!functionallyEquivalent(roots_[i], newroots[i])) {
            std::cout << "found non-equivalent nodes #" << i << std::endl;
            std::cout << roots_[i] << std::endl << newroots[i] << std::endl;

            diff = true;
        }
    }

    return !diff;
}

bool
aigpp::EquivalenceChecker::checkNewRoots(const std::vector<aigpp::InternalEdgeRef>& roots)
{
    assert(roots.size() == roots_.size());

    /* look for original manager */
    Manager* origManager = getManager(roots);

    /* create variables */
    copyVariables(origManager);

    /* copy roots */
    std::vector<EdgeRef>             newroots;
    std::map<Node*, InternalEdgeRef> mapping;
    for (const auto& r : roots) {
        InternalEdgeRef newe = copyCone(r, mapping);
        newroots.push_back(man_.DebugGetExternal(newe));
    }

    /* check roots */
    bool diff = false;
    for (unsigned int i = 0; i != roots_.size(); ++i) {
        if (!functionallyEquivalent(roots_[i], newroots[i])) {
            std::cout << "found non-equivalent nodes #" << i << std::endl;
            /*
            std::cout << roots_[i] << std::endl
                      << newroots[i] << std::endl;
            */
            diff = true;
        }
    }

    return !diff;
}

aigpp::InternalEdgeRef
aigpp::EquivalenceChecker::copyCone(const aigpp::InternalEdgeRef&            root,
                                    std::map<Node*, aigpp::InternalEdgeRef>& mapping)
{
    if (root.isConstant()) {
        if (root.functionallyTrue()) {
            return const1;
        } else {
            assert(root.functionallyFalse());
            return const0;
        }
    }

    std::map<Node*, InternalEdgeRef>::const_iterator p1, p2;

    static std::stack<Node*> pending;
    pending.push(root.node());
    while (!pending.empty()) {
        if (mapping.find(pending.top()) != mapping.end()) {
            pending.pop();
        } else if (pending.top()->isVar()) {
            int v = pending.top()->varIndex();
            assert(man_.isVarIndexValid(v));
            mapping[pending.top()] = man_.variable(v).getInternal();
            pending.pop();
        } else if ((p1 = mapping.find(pending.top()->parent1().node())) == mapping.end()) {
            pending.push(pending.top()->parent1().node());
        } else if ((p2 = mapping.find(pending.top()->parent2().node())) == mapping.end()) {
            pending.push(pending.top()->parent2().node());
        } else {
            InternalEdgeRef e1     = p1->second.notIf(pending.top()->parent1().isInverted());
            InternalEdgeRef e2     = p2->second.notIf(pending.top()->parent2().isInverted());
            InternalEdgeRef r      = e1 & e2;
            mapping[pending.top()] = r.notIf(pending.top()->isNAND());
            pending.pop();
        }
    }

    std::map<Node*, InternalEdgeRef>::const_iterator p = mapping.find(root.node());
    assert(p != mapping.end());

    return p->second.notIf(root.isInverted());
}

void
aigpp::EquivalenceChecker::copyVariables(aigpp::Manager* manager)
{
    if (manager == nullptr) {
        return;
    }
    if (manager->variableCount() <= man_.variableCount()) {
        return;
    }

    for (std::size_t v = man_.variableCount(); v != manager->variableCount(); ++v) {
        man_.addVariable();
    }
}

aigpp::Manager*
aigpp::EquivalenceChecker::getManager(const std::vector<aigpp::EdgeRef>& roots) const
{
    Manager* manager = nullptr;

    for (const auto& r : roots) {
        const InternalEdgeRef& e = r.getInternal();
        Node*                  n = e.node();

        if (n != nullptr) {
            if (manager == nullptr) {
                manager = n->manager();
            } else {
                assert(manager == n->manager());
            }
        }
    }

    return manager;
}

aigpp::Manager*
aigpp::EquivalenceChecker::getManager(const std::vector<aigpp::InternalEdgeRef>& roots) const
{
    Manager* manager = nullptr;

    for (const auto& r : roots) {
        Node* n = r.node();

        if (n != nullptr) {
            if (manager == nullptr) {
                manager = n->manager();
            } else {
                assert(manager == n->manager());
            }
        }
    }

    return manager;
}

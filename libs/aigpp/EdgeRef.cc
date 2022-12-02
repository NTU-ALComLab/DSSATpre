/**************************************************************
 *
 *       AIGPP Package // EdgeRef.cc
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

#include "EdgeRef.hh"

/* aigpp includes */
#include "ExtRefTable.hh"
#include "Manager.hh"
#include "Node.hh"
#include "VarSupport.hh"

aigpp::EdgeRef::EdgeRef() : _extRefTable(nullptr), _index(0)
{}

aigpp::EdgeRef::EdgeRef(const aigpp::EdgeRef& base) : _extRefTable(base._extRefTable), _index(base._index)
{
    if (isValid()) {
        _extRefTable->ref(_index);
    }
}

aigpp::EdgeRef::EdgeRef(aigpp::EdgeRef&& base) noexcept : _extRefTable(base._extRefTable), _index(base._index)
{
    base._extRefTable = nullptr;
    base._index       = 0;
}

aigpp::EdgeRef::~EdgeRef()
{
    if (isValid()) {
        _extRefTable->deref(_index);
    }
}

aigpp::EdgeRef&
aigpp::EdgeRef::operator=(const aigpp::EdgeRef& base)
{
    if (_extRefTable == base._extRefTable && _index == base._index) {
        return *this;
    }

    /* increment new target ref counter first */
    if (base.isValid()) {
        base._extRefTable->ref(base._index);
    }

    /* decrement old ref counter */
    if (isValid()) {
        _extRefTable->deref(_index);
    }

    /* update members */
    _extRefTable = base._extRefTable;
    _index       = base._index;

    return *this;
}

aigpp::EdgeRef&
aigpp::EdgeRef::operator=(aigpp::EdgeRef&& base)
{
    if (_extRefTable == base._extRefTable && _index == base._index) {
        return *this;
    }

    /* decrement old ref counter */
    if (isValid()) {
        _extRefTable->deref(_index);
    }

    /* update members */
    _extRefTable = base._extRefTable;
    _index       = base._index;

    base._extRefTable = nullptr;
    base._index       = 0;

    return *this;
}

void
aigpp::EdgeRef::invalidate()
{
    if (isValid()) {
        _extRefTable->deref(_index);
    }

    _extRefTable = nullptr;
    _index       = 0;
}

/* inline */
bool
aigpp::EdgeRef::isValid() const
{
    return _extRefTable != nullptr;
}

bool
aigpp::EdgeRef::isConstant() const
{
    return getInternal().isConstant();
}

bool
aigpp::EdgeRef::isVariable() const
{
    Node* n = getInternal().node();
    assert(n != nullptr);
    return n->isVar();
}

bool
aigpp::EdgeRef::isInverted() const
{
    return getInternal().isInverted();
}

int
aigpp::EdgeRef::variableIndex() const
{
    return getInternal().node()->varIndex();
}

const std::string&
aigpp::EdgeRef::variableName() const
{
    Node* n = getInternal().node();
    assert(n != nullptr);
    assert(n->isVar());

    return n->manager()->variableName(n->varIndex());
}

/* inline */
std::size_t
aigpp::EdgeRef::nodeCount() const
{
    return getInternal().nodeCount();
}

std::vector<int>
aigpp::EdgeRef::support() const
{
    Node* internal = getInternal().node();

    /* constant node */
    if (internal == nullptr) {
        return std::vector<int>();
    }

    Manager* m = internal->manager();

    VarSupport supp = internal->support();

    /* supp must not be empty since "internal" is a non-constant node! */
    assert(supp.vars() > 0);

    std::vector<int> result;
    result.reserve(supp.vars());

    for (std::size_t v = 0; v != m->variableCount(); ++v) {
        if (supp.hasVar(v)) {
            result.push_back((int)v);
        }
    }

    assert(result.size() == supp.vars());

    return result;
}

std::vector<aigpp::EdgeRef>
aigpp::EdgeRef::supportVariables() const
{
    Node* internal = getInternal().node();

    /* constant node */
    if (internal == nullptr) {
        return std::vector<EdgeRef>();
    }

    Manager* m = internal->manager();

    VarSupport supp = internal->support();

    /* supp must not be empty since "internal" is a non-constant node! */
    assert(supp.vars() > 0);

    std::vector<EdgeRef> result;
    result.reserve(supp.vars());

    for (std::size_t v = 0; v != m->variableCount(); ++v) {
        if (supp.hasVar(v)) {
            result.push_back(m->variable(v));
        }
    }

    assert(result.size() == supp.vars());

    return result;
}

std::pair<bool, aigpp::VariableAssignment>
aigpp::EdgeRef::satisfyingAssignment(bool minimize) const
{
    InternalEdgeRef internal = getInternal();

    /* check for constants */
    if (internal.isConstant()) {
        if (structurallyFalse(internal)) {
            /* const0 -> not satisfiable */
            return std::pair<bool, VariableAssignment>(false, VariableAssignment());
        } else {
            /* const1 -> satisfiable for all assignments */
            return std::pair<bool, VariableAssignment>(true, VariableAssignment());
        }
    }

    /* call backend function */
    Manager* m = internal.node()->manager();
    return m->satisfyingAssignment(internal, minimize);
}

aigpp::EdgeRef
aigpp::EdgeRef::notIf(bool i) const
{
    if (i) {
        return Not(*this);
    } else {
        return *this;
    }
}

aigpp::EdgeRef
aigpp::EdgeRef::cofactor(const aigpp::EdgeRef& variable) const
{
    return EdgeRef(_extRefTable, getInternal().cofactor(variable.getInternal()));
}

aigpp::EdgeRef
aigpp::EdgeRef::cofactor(const std::vector<aigpp::EdgeRef>& variables) const
{
    std::vector<InternalEdgeRef> variables2;
    variables2.reserve(variables.size());
    std::transform(variables.cbegin(), variables.cend(), std::back_inserter(variables2),
                   [](const auto& p) -> InternalEdgeRef { return p.getInternal(); });

    return EdgeRef(_extRefTable, getInternal().cofactor(variables2));
}

aigpp::EdgeRef
aigpp::EdgeRef::compose(const std::vector<aigpp::EdgeRef>& variables,
                        const std::vector<aigpp::EdgeRef>& replacements) const
{
    std::vector<InternalEdgeRef> variables2    = toInternal(variables);
    std::vector<InternalEdgeRef> replacements2 = toInternal(replacements);

    return EdgeRef(_extRefTable, getInternal().compose(variables2, replacements2));
}

aigpp::EdgeRef
aigpp::EdgeRef::existentialQ(const aigpp::EdgeRef& variable) const
{
    return EdgeRef(_extRefTable, getInternal().existentialQ(variable.getInternal()));
}

aigpp::EdgeRef
aigpp::EdgeRef::existentialQ(const std::vector<aigpp::EdgeRef>& variables) const
{
    std::vector<InternalEdgeRef> variables2;
    variables2.reserve(variables.size());
    std::transform(variables.cbegin(), variables.cend(), std::back_inserter(variables2),
                   [](const auto& p) -> InternalEdgeRef { return p.getInternal(); });

    return EdgeRef(_extRefTable, getInternal().existentialQ(variables2));
}

aigpp::EdgeRef
aigpp::EdgeRef::universalQ(const aigpp::EdgeRef& variable) const
{
    return EdgeRef(_extRefTable, getInternal().universalQ(variable.getInternal()));
}

aigpp::EdgeRef
aigpp::EdgeRef::universalQ(const std::vector<aigpp::EdgeRef>& variables) const
{
    std::vector<InternalEdgeRef> variables2;
    variables2.reserve(variables.size());

    std::transform(variables.cbegin(), variables.cend(), std::back_inserter(variables2),
                   [](const auto& p) -> InternalEdgeRef { return p.getInternal(); });

    return EdgeRef(_extRefTable, getInternal().universalQ(variables2));
}

/* inline */
const aigpp::InternalEdgeRef&
aigpp::EdgeRef::getInternal() const
{
    assert(isValid());
    return _extRefTable->lookup(_index);
}

aigpp::EdgeRef::EdgeRef(aigpp::ExtRefTable* extRefTable, const InternalEdgeRef& internal) :
    _extRefTable(extRefTable),
    _index(_extRefTable->insert(internal))
{}

std::vector<aigpp::InternalEdgeRef>
aigpp::EdgeRef::toInternal(const std::vector<EdgeRef>& v)
{
    std::vector<InternalEdgeRef> result;
    result.reserve(v.size());

    std::transform(v.cbegin(), v.cend(), std::back_inserter(result),
                   [](const auto& p) -> InternalEdgeRef { return p.getInternal(); });

    return result;
}

std::vector<aigpp::EdgeRef>
aigpp::EdgeRef::fromInternal(ExtRefTable* extRefTable, const std::vector<InternalEdgeRef>& v)
{
    std::vector<EdgeRef> result;
    result.reserve(v.size());

    std::transform(v.cbegin(), v.cend(), std::back_inserter(result),
                   [&extRefTable](const auto& p) -> EdgeRef { return EdgeRef(extRefTable, p); });

    return result;
}

std::ostream&
aigpp::operator<<(std::ostream& os, const EdgeRef& e)
{
    os << e.getInternal();
    return os;
}

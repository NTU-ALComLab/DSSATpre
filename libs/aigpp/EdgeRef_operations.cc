/**************************************************************
 *
 *       AIGPP Package // EdgeRef_operations.cc
 *
 *       Copyright (C) 2007 Florian Pigorsch
 *
 *       Author:
 *         Florian Pigorsch
 *         University of Freiburg
 *         pigorsch@informatik.uni-freiburg.de
 *
 *       Last revision:
 *         $Revision: $
 *         $Author: $
 *
 ***************************************************************/

#include "EdgeRef.hh"

/* aigpp */
#include "InternalEdgeRef.hh"
#include "Manager.hh"

aigpp::EdgeRef
aigpp::Not(const aigpp::EdgeRef& e)
{
    const InternalEdgeRef& i = e.getInternal();

    InternalEdgeRef result = !i;

    return EdgeRef(e._extRefTable, result);
}

aigpp::EdgeRef
aigpp::And(const aigpp::EdgeRef& e1, const aigpp::EdgeRef& e2)
{
    const InternalEdgeRef& i1 = e1.getInternal();
    const InternalEdgeRef& i2 = e2.getInternal();

    InternalEdgeRef result = i1 & i2;

    return EdgeRef(e1._extRefTable, result);
}

aigpp::EdgeRef
aigpp::Or(const aigpp::EdgeRef& e1, const aigpp::EdgeRef& e2)
{
    const InternalEdgeRef& i1 = e1.getInternal();
    const InternalEdgeRef& i2 = e2.getInternal();

    InternalEdgeRef result = i1 | i2;

    return EdgeRef(e1._extRefTable, result);
}

aigpp::EdgeRef
aigpp::Xor(const aigpp::EdgeRef& e1, const aigpp::EdgeRef& e2)
{
    const InternalEdgeRef& i1 = e1.getInternal();
    const InternalEdgeRef& i2 = e2.getInternal();

    InternalEdgeRef result = i1 ^ i2;

    return EdgeRef(e1._extRefTable, result);
}

aigpp::EdgeRef
aigpp::Equiv(const aigpp::EdgeRef& e1, const aigpp::EdgeRef& e2)
{
    const InternalEdgeRef& i1 = e1.getInternal();
    const InternalEdgeRef& i2 = e2.getInternal();

    InternalEdgeRef result = !(i1 ^ i2);

    return EdgeRef(e1._extRefTable, result);
}

aigpp::EdgeRef
aigpp::Mux(const aigpp::EdgeRef& select, const aigpp::EdgeRef& eThen, const aigpp::EdgeRef& eElse)
{
    const InternalEdgeRef& iSelect = select.getInternal();
    const InternalEdgeRef& iThen   = eThen.getInternal();
    const InternalEdgeRef& iElse   = eElse.getInternal();

    InternalEdgeRef result = iSelect.Mux(iThen, iElse);

    return EdgeRef(select._extRefTable, result);
}

aigpp::EdgeRef
aigpp::multiAnd(aigpp::Manager& manager, const std::vector<aigpp::EdgeRef>& e)
{
    std::vector<InternalEdgeRef> internal = EdgeRef::toInternal(e);

    InternalEdgeRef result = InternalEdgeRef::multiAnd(internal);

    return EdgeRef(manager.getExtRefTable(), result);
}

aigpp::EdgeRef
aigpp::multiOr(aigpp::Manager& manager, const std::vector<aigpp::EdgeRef>& e)
{
    std::vector<InternalEdgeRef> internal = EdgeRef::toInternal(e);

    InternalEdgeRef result = InternalEdgeRef::multiOr(internal);

    return EdgeRef(manager.getExtRefTable(), result);
}

bool
aigpp::structurallyEquivalent(const aigpp::EdgeRef& e1, const aigpp::EdgeRef& e2)
{
    const InternalEdgeRef& i1 = e1.getInternal();
    const InternalEdgeRef& i2 = e2.getInternal();

    return structurallyEquivalent(i1, i2);
}

bool
aigpp::structurallyAntivalent(const aigpp::EdgeRef& e1, const aigpp::EdgeRef& e2)
{
    const InternalEdgeRef& i1 = e1.getInternal();
    const InternalEdgeRef& i2 = e2.getInternal();

    return structurallyAntivalent(i1, i2);
}

bool
aigpp::structurallyTrue(const aigpp::EdgeRef& e)
{
    const InternalEdgeRef& i = e.getInternal();

    return structurallyTrue(i);
}

bool
aigpp::structurallyFalse(const aigpp::EdgeRef& e)
{
    const InternalEdgeRef& i = e.getInternal();

    return structurallyFalse(i);
}

bool
aigpp::functionallyEquivalent(const aigpp::EdgeRef& e1, const aigpp::EdgeRef& e2)
{
    const InternalEdgeRef& i1 = e1.getInternal();
    const InternalEdgeRef& i2 = e2.getInternal();

    return functionallyEquivalent(i1, i2);
}

bool
aigpp::functionallyAntivalent(const aigpp::EdgeRef& e1, const aigpp::EdgeRef& e2)
{
    const InternalEdgeRef& i1 = e1.getInternal();
    const InternalEdgeRef& i2 = e2.getInternal();

    return functionallyAntivalent(i1, i2);
}

bool
aigpp::functionallyTrue(const aigpp::EdgeRef& e)
{
    const InternalEdgeRef& i = e.getInternal();

    return functionallyTrue(i);
}

bool
aigpp::functionallyFalse(const aigpp::EdgeRef& e)
{
    const InternalEdgeRef& i = e.getInternal();

    return functionallyFalse(i);
}

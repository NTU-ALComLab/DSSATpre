#include "AigerParser.hh"

#include <stack>
#include <unordered_map>

#include <lrabsutil/Assert.hh>
#include <lrabsutil/String.hh>

#define _unused(x) ((void)(x))

aigpp::AigerModel::AigerModel() :
    maxvar_(0),
    numInputs_(0),
    numLatches_(0),
    numOutputs_(0),
    numAnds_(0),
    computedDependencies_(false)
{
    types_.emplace_back(CONSTANT0, -1);
}

aigpp::AigerModel::~AigerModel()
{
    clear();
}

void
aigpp::AigerModel::parse(std::istream& is)
{
    clear();

    std::string magic;
    is >> magic >> maxvar_ >> numInputs_ >> numLatches_ >> numOutputs_ >> numAnds_;
    assert(magic == "aag" || magic == "aig");
    assert(maxvar_ >= (numInputs_ + numLatches_ + numAnds_));

    assert(types_.size() == 1);
    types_.insert(types_.end(), maxvar_, TypeInfo(UNKNOWN, -1));

    assert(inputs_.empty());
    inputs_.assign(numInputs_, (InputNode*)nullptr);

    assert(latches_.empty());
    latches_.assign(numLatches_, (LatchNode*)nullptr);

    assert(outputs_.empty());
    outputs_.assign(numOutputs_, (OutputNode*)nullptr);

    assert(ands_.empty());
    ands_.assign(numAnds_, (AndNode*)nullptr);

    // AIGER is in ASCII format
    if (magic == "aag") {
        for (int i = 0; i != numInputs_; ++i) {
            int id;
            is >> id;
            assert((id & 1) == 0);
            int index = id / 2;
            assert(index <= maxvar_);

            TypeInfo& t = types_[index];
            assert(t.type_ == UNKNOWN);
            t.type_  = INPUT;
            t.index_ = i;

            auto n     = new InputNode();
            n->signal_ = id;
            assert(inputs_[i] == 0);
            inputs_[i] = n;
        }

        for (int i = 0; i != numLatches_; ++i) {
            int id1, id2;
            is >> id1 >> id2;
            assert((id1 & 1) == 0);
            int index1 = id1 / 2;
#ifndef NDEBUG
            assert(index1 <= maxvar_);
            int index2 = id2 / 2;
            assert(index2 <= maxvar_);
#endif
            TypeInfo& t = types_[index1];
            assert(t.type_ == UNKNOWN);
            t.type_  = LATCH;
            t.index_ = i;

            auto n     = new LatchNode();
            n->signal_ = id1;
            n->input_  = id2;
            assert(latches_[i] == 0);
            latches_[i] = n;
        }

        for (int i = 0; i != numOutputs_; ++i) {
            int id;
            is >> id;
#ifndef NDEBUG
            int index = id / 2;
            assert(index <= maxvar_);
#endif
            auto n     = new OutputNode();
            n->signal_ = id;
            assert(outputs_[i] == 0);
            outputs_[i] = n;
        }

        for (int i = 0; i != numAnds_; ++i) {
            int id, id1, id2;
            is >> id >> id1 >> id2;
            assert((id & 1) == 0);
            int index = id / 2;
#ifndef NDEBUG
            assert(index <= maxvar_);
            int index1 = id1 / 2;
            assert(index1 <= maxvar_);
            int index2 = id2 / 2;
            assert(index2 <= maxvar_);
#endif
            TypeInfo& t = types_[index];
            assert(t.type_ == UNKNOWN);
            t.type_  = AND;
            t.index_ = i;

            auto n     = new AndNode();
            n->signal_ = id;
            n->input1_ = id1;
            n->input2_ = id2;
            assert(ands_[i] == 0);
            ands_[i] = n;
        }
    } else {  // AIGER format is binary
        for (int i = 0; i != numInputs_; ++i) {
            const int id = 2 * (i + 1);
            assert((id & 1) == 0);
            const int index = id / 2;
            assert(index <= maxvar_);

            TypeInfo& t = types_[index];
            assert(t.type_ == UNKNOWN);
            t.type_  = INPUT;
            t.index_ = i;

            auto n     = new InputNode();
            n->signal_ = id;
            assert(inputs_[i] == 0);
            inputs_[i] = n;
        }

        assert(numLatches_ == 0);

        for (int i = 0; i != numOutputs_; ++i) {
            int id = 0;
            is >> id;
#ifndef NDEBUG
            const int index = id / 2;
            assert(index <= maxvar_);
#endif
            auto n     = new OutputNode();
            n->signal_ = id;
            assert(outputs_[i] == 0);
            outputs_[i] = n;
        }

        const char ch = is.get();
        _unused(ch);
        assert('\n' == ch);  // read line break before binary block

        for (int i = 0; i != numAnds_; ++i) {
            const unsigned int id = 2 * (numInputs_ + (i + 1));
            assert((id & 1) == 0);
            const int index = id / 2;
            assert(index <= maxvar_);

            const size_t delta0 = readBinary(is);
            const size_t rhs0   = id - delta0;
#ifndef NDEBUG
            const int index1 = rhs0 / 2;
            assert(index1 <= maxvar_);
#endif
            const size_t delta1 = readBinary(is);
            const size_t rhs1   = rhs0 - delta1;
#ifndef NDEBUG
            const int index2 = rhs1 / 2;
            assert(index2 <= maxvar_);
#endif
            TypeInfo& t = types_[index];
            assert(t.type_ == UNKNOWN);
            t.type_  = AND;
            t.index_ = i;

            auto n     = new AndNode();
            n->signal_ = id;
            n->input1_ = rhs0;
            n->input2_ = rhs1;
            assert(ands_[i] == 0);
            ands_[i] = n;
        }
    }

    /* read signal names */
    std::string s;
    while (is >> s) {
        if (s.empty()) {
            break;
        } else if (s == "c") {
            /* start of comments section */
            std::getline(is, s);
            while (std::getline(is, s)) {
                comments_.push_back(s);
            }
            break;
        }

        else if (s[0] == 'i' || s[0] == 'l' || s[0] == 'o') {
            assert(s.size() >= 2);
            int index = 0;
            for (unsigned int i = 1; i != s.size(); ++i) {
                assert(s[i] >= '0' && s[i] <= '9');
                index = index * 10 + (s[i] - '0');
            }

            std::string annotation;
            std::getline(is, annotation);
            lrabs::trim(annotation);

            if (s[0] == 'i') {
                assert(index < numInputs_);
                assert(inputs_[index] != 0);
                inputs_[index]->name_ = annotation;
            } else if (s[0] == 'l') {
                assert(index < numLatches_);
                assert(latches_[index] != 0);
                latches_[index]->name_ = annotation;
            } else if (s[0] == 'o') {
                assert(index < numOutputs_);
                assert(outputs_[index] != 0);
                outputs_[index]->name_ = annotation;
            }
        } else {
            assert(false);
        }
    }

    checkConsistency();
}

/*
 *\brief Read characters in binary format and converts to unsigned int.
 *
 * An unsigned int x is binary encoded in the sequence of i bites b0,..., bi:
 *1w0, 1w1, ..., 1w(i-1), 0wi. This function reads characters until the last
 *character (0wi), that defines x, is read and decodes it to unsigned int.
 */
std::size_t
aigpp::AigerModel::readBinary(std::istream& is) const
{
    std::size_t   x = 0, i = 0;
    unsigned char ch;

    while ((ch = getChar(is)) & 0x80u) {  // first bit is set to 1 --> input word is not ended yet
        x |= (ch & 0x7fu) << (7u * i++);  // set bits 1 - 7, bit 8 (MSB) is set to 0,
                                          // shift i times 7 bits and increment i
    }

    return x | (ch << (7u * i));  // first bit was 0 --> input word is ended, shift 7 bits
}

/*
 * \brief Functions reads one char ch from istream and returns ch, if ch is not
 * EOF.
 */
unsigned char
aigpp::AigerModel::getChar(std::istream& is) const
{
    const unsigned char ch = is.get();

#ifndef NDEBUG
    if (is.eof()) std::cerr << "Unexpected EOF.\n\n";
    assert(!is.eof());
    if (is.fail()) std::cerr << "Reading char failed.\n\n";
    assert(!is.fail());
#endif

    return ch;
}

aigpp::AigerModel::AIGRepresentation
aigpp::AigerModel::computeAIGRepresentation(aigpp::Manager& man) const
{
    AIGRepresentation repr;
    repr.inputs.reserve(numInputs_);
    repr.latchInputFunctions.reserve(numLatches_);
    repr.latches.reserve(numLatches_);
    repr.outputs.reserve(numOutputs_);

    for (const TypeInfo& p : types_) {
        p.computed_ = false;
    }

    using AIGMap = std::unordered_map<int, aigpp::EdgeRef>;

    AIGMap amap;
    amap[0]             = man.getConst0();
    types_[0].computed_ = true;

    for (InputNode* p : inputs_) {
        EdgeRef e = (man.hasVariableName(p->name_) ? man.lookupVariable(p->name_) : man.addVariable(p->name_));
        repr.inputs.push_back(e);
        amap[p->signal_ / 2]             = e;
        types_[p->signal_ / 2].computed_ = true;
    }

    for (LatchNode* p : latches_) {
        EdgeRef e = (man.hasVariableName(p->name_) ? man.lookupVariable(p->name_) : man.addVariable(p->name_));
        repr.latches.push_back(e);
        amap[p->signal_ / 2]             = e;
        types_[p->signal_ / 2].computed_ = true;
    }

    std::stack<int> pending;
    for (LatchNode* p : latches_) {
        pending.push(p->input_ / 2);
    }
    for (OutputNode* p : outputs_) {
        pending.push(p->signal_ / 2);
    }

    while (!pending.empty()) {
        const TypeInfo& t = types_[pending.top()];

        if (t.computed_) {
            pending.pop();
            continue;
        }

        assert(t.type_ == AND);

        AndNode* n = ands_[t.index_];
        if (!types_[n->input1_ / 2].computed_) {
            pending.push(n->input1_ / 2);
            continue;
        }
        if (!types_[n->input2_ / 2].computed_) {
            pending.push(n->input2_ / 2);
            continue;
        }

        EdgeRef p1 = amap[n->input1_ / 2].notIf((n->input1_ & 1u) == 1u);
        EdgeRef p2 = amap[n->input2_ / 2].notIf((n->input2_ & 1u) == 1u);

        EdgeRef e = p1 & p2;

        amap[pending.top()] = e;
        t.computed_         = true;
    }

    for (LatchNode* p : latches_) {
        assert(types_[p->input_ / 2].computed_);
        repr.latchInputFunctions.push_back(amap[p->input_ / 2].notIf((p->input_ & 1u) == 1u));
    }
    for (OutputNode* p : outputs_) {
        assert(types_[p->signal_ / 2].computed_);
        repr.outputs.push_back(amap[p->signal_ / 2].notIf((p->signal_ & 1u) == 1u));
    }

    return repr;
}

std::vector<std::vector<unsigned int>>
aigpp::AigerModel::outputs2CNF(const std::map<std::string, unsigned int>& inputMapping, unsigned int& nextfreevar) const
{
    std::vector<std::vector<unsigned int>> cnf;

    std::map<int, unsigned int> signal2cnflit;

    for (InputNode* i : inputs_) {
        auto p = inputMapping.find(i->name_);
        assert(p != inputMapping.end());

        signal2cnflit[i->signal_ / 2] = 2 * (p->second);
    }

    std::stack<int> pending;

    for (OutputNode* o : outputs_) {
        if (o->signal_ > 1) {
            pending.push(o->signal_ / 2);
        }
    }

    while (!pending.empty()) {
        int sig = pending.top();

        if (signal2cnflit.find(sig) != signal2cnflit.end()) {
            pending.pop();
            continue;
        }

        const TypeInfo& t = types_[sig];

        assert(t.type_ == AND);

        AndNode* a = ands_[t.index_];

        if (signal2cnflit.find(a->input1_ / 2) == signal2cnflit.end()) {
            pending.push(a->input1_ / 2);
            continue;
        } else if (signal2cnflit.find(a->input2_ / 2) == signal2cnflit.end()) {
            pending.push(a->input2_ / 2);
            continue;
        }

        unsigned int out = 2 * nextfreevar;
        ++nextfreevar;
        signal2cnflit[sig] = out;

        unsigned int a1 = signal2cnflit[a->input1_ / 2u] + (a->input1_ & 1u);
        unsigned int a2 = signal2cnflit[a->input2_ / 2u] + (a->input2_ & 1u);

        std::vector<unsigned int> clause3(3), clause2(2);

        clause3[0] = out;
        clause3[1] = a1 ^ 1u;
        clause3[2] = a2 ^ 1u;

        cnf.push_back(clause3);

        clause2[0] = out ^ 1u;
        clause2[1] = a1;

        cnf.push_back(clause2);

        clause2[0] = out ^ 1u;
        clause2[1] = a2;

        cnf.push_back(clause2);

        pending.pop();
    }

    for (OutputNode* o : outputs_) {
        if (o->signal_ > 1) {
            std::vector<unsigned int> unit(1);
            unit[0] = signal2cnflit[o->signal_ / 2u] + (o->signal_ & 1u);
            cnf.push_back(unit);
        } else {
            if (o->signal_ == 0) {
                std::cout << "WARNING: aiger output is 'FALSE' -> generating "
                             "unsatisfiable CNF\n";
                cnf.clear();
                cnf.emplace_back();
                break;
            } else {
                std::cout << "WARNING: aiger output is 'TRUE' -> omitting output int "
                             "CNF generation\n";
            }
        }
    }

    return cnf;
}

void
aigpp::AigerModel::computeDependencies()
{
    if (computedDependencies_) {
        return;
    }
    computedDependencies_ = true;

    /*
     * first vector: indices of inputs
     * second vector: indices of latches
     */
    using SetPair = std::pair<std::set<int>, std::set<int>>;
    using DepMap  = std::unordered_map<int, SetPair>;
    DepMap dmap;

    for (TypeInfo& p : types_) {
        p.computed_ = false;
    }

    types_[0].computed_ = true;

    for (InputNode* p : inputs_) {
        SetPair vp;
        // vp.first.insert( ( *p )->signal_/2 );
        vp.first.insert(types_[p->signal_ / 2].index_);
        dmap[p->signal_ / 2] = vp;

        types_[p->signal_ / 2].computed_ = true;
    }

    for (LatchNode* p : latches_) {
        SetPair vp;
        vp.second.insert(types_[p->signal_ / 2].index_);
        dmap[p->signal_ / 2] = vp;

        types_[p->signal_ / 2].computed_ = true;
    }

    std::stack<int> pending;
    for (LatchNode* p : latches_) {
        pending.push(p->input_ / 2);
    }
    for (OutputNode* p : outputs_) {
        pending.push(p->signal_ / 2);
    }

    while (!pending.empty()) {
        const TypeInfo& t = types_[pending.top()];

        if (t.computed_) {
            pending.pop();
            continue;
        }

        assert(t.type_ == AND);

        AndNode* n = ands_[t.index_];
        if (!types_[n->input1_ / 2].computed_) {
            pending.push(n->input1_ / 2);
            continue;
        }
        if (!types_[n->input2_ / 2].computed_) {
            pending.push(n->input2_ / 2);
            continue;
        }

        const SetPair vp1 = dmap[n->input1_ / 2];
        const SetPair vp2 = dmap[n->input2_ / 2];

        /* compute vp := union( vp1, vp2 ) */
        SetPair vp = vp1;
        for (int p : vp2.first) {
            vp.first.insert(p);
        }
        for (int p : vp2.second) {
            vp.second.insert(p);
        }

        dmap[pending.top()] = vp;
        t.computed_         = true;
    }

    for (LatchNode* p : latches_) {
        assert(types_[p->input_ / 2].computed_);
        p->dependsOnInputs_  = dmap[p->input_ / 2].first;
        p->dependsOnLatches_ = dmap[p->input_ / 2].second;
    }
    for (OutputNode* p : outputs_) {
        assert(types_[p->signal_ / 2].computed_);
        p->dependsOnInputs_  = dmap[p->signal_ / 2].first;
        p->dependsOnLatches_ = dmap[p->signal_ / 2].second;
    }
}

void
aigpp::AigerModel::printDependencies(std::ostream& os) const
{
    assert(computedDependencies_);

    for (unsigned int i = 0; i != latches_.size(); ++i) {
        const LatchNode* n = latches_[i];
        os << "latch[" << i << "/" << latches_.size() << "] name='" << n->name_ << "' index=" << n->signal_
           << " input=" << n->input_ << " depends on\n";

        os << "latch(es)";
        for (int p : n->dependsOnLatches_) {
            os << " " << p;
        }
        os << '\n';

        os << "input(s)";
        for (int p : n->dependsOnInputs_) {
            os << " " << p;
        }
        os << '\n';
    }

    for (unsigned int i = 0; i != outputs_.size(); ++i) {
        const OutputNode* n = outputs_[i];
        os << "output[" << i << "/" << outputs_.size() << "] name='" << n->name_ << "' index=" << n->signal_
           << " depends on\n";

        os << "latch(es)";
        for (int p : n->dependsOnLatches_) {
            os << " " << p;
        }
        os << '\n';

        os << "input(s)";
        for (int p : n->dependsOnInputs_) {
            os << " " << p;
        }
        os << '\n';
    }
}

void
aigpp::AigerModel::clear()
{
    maxvar_     = 0;
    numInputs_  = 0;
    numLatches_ = 0;
    numOutputs_ = 0;
    numAnds_    = 0;

    types_.clear();
    types_.emplace_back(CONSTANT0, -1);

    for (InputNode* p : inputs_) {
        delete p;
    }
    inputs_.clear();

    for (LatchNode* p : latches_) {
        delete p;
    }
    latches_.clear();

    for (OutputNode* p : outputs_) {
        delete p;
    }
    outputs_.clear();

    for (AndNode* p : ands_) {
        delete p;
    }
    ands_.clear();

    comments_.clear();
}

void
aigpp::AigerModel::checkConsistency() const
{
#ifndef NDEBUG
    /* check if latch inputs are connected to defined signals */
    for (LatchNode* p : latches_) {
        assert(types_[p->input_ / 2].type_ != UNKNOWN);
    }

    /* check if outputs are connected to defined signals */
    for (OutputNode* p : outputs_) {
        assert(types_[p->signal_ / 2].type_ != UNKNOWN);
    }

    /* check if ands are connected to defined signals */
    for (AndNode* p : ands_) {
        assert(types_[p->input1_ / 2].type_ != UNKNOWN);
        assert(types_[p->input2_ / 2].type_ != UNKNOWN);
    }
#endif
    /* check if there are no cycles */
    /* TODO */
}

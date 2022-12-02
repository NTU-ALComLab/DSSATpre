#ifndef AIGPP_AIGERPARSER_HH
#define AIGPP_AIGERPARSER_HH

#include <iostream>
#include <vector>

#include <lrabsutil/Exception.hh>

#include "EdgeRef.hh"
#include "Manager.hh"

/* sample aiger file

   aag 7 2 1 2 4               aag M(max var index) I(number of inputs) L(number
   of latches) O(number of outputs) A(number of and gates) 2 input 0 'enable' 4
   input 1         'reset' 6 8                         latch 0         Q next(Q)
   6                           output 0        Q
   7                           output 1        !Q
   8 4 10                      AND gate 0      reset & (enable ^ Q)
   10 13 15                    AND gate 1      enable ^ Q
   12 2 6                      AND gate 2      enable & Q
   14 3 7                      AND gate 3      !enable & !Q
   l0 latchQ
   i0 enable
   i1 reset
   o0 Q
   o1 notQ
   c
   flip flop

*/
namespace aigpp {

class AigerModel
{
   public:
    AigerModel();
    AigerModel(AigerModel&& other)      = default;
    AigerModel(const AigerModel& other) = delete;
    ~AigerModel();
    AigerModel& operator=(AigerModel&& other) = default;
    AigerModel& operator=(const AigerModel& other) = delete;

    void          parse(std::istream& is);
    std::size_t   readBinary(std::istream& is) const;
    unsigned char getChar(std::istream& is) const;

    struct AIGRepresentation
    {
        std::vector<EdgeRef> inputs;
        std::vector<EdgeRef> latchInputFunctions;
        std::vector<EdgeRef> latches;
        std::vector<EdgeRef> outputs;
    };

    AIGRepresentation computeAIGRepresentation(Manager& man) const;

    std::vector<std::vector<unsigned int>> outputs2CNF(const std::map<std::string, unsigned int>& inputMapping,
                                                       unsigned int&                              nextfreevar) const;

    void computeDependencies();
    void printDependencies(std::ostream& os = std::cout) const;

    void clear();
    void checkConsistency() const;

    int maxvar_;
    int numInputs_;
    int numLatches_;
    int numOutputs_;
    int numAnds_;

    enum Type
    {
        UNKNOWN,
        CONSTANT0,
        INPUT,
        AND,
        LATCH
    };
    struct TypeInfo
    {
        TypeInfo() : type_(UNKNOWN), index_(-1), computed_(false) {}

        TypeInfo(Type t, int i) : type_(t), index_(i), computed_(false) {}

        Type         type_;
        int          index_;
        mutable bool computed_;
    };

    std::vector<TypeInfo> types_;

    struct InputNode
    {
        int         signal_ = 0;
        std::string name_;
    };
    struct LatchNode
    {
        std::string name_;
        int         signal_ = 0;
        int         input_  = 0;

        std::set<int> dependsOnLatches_;
        std::set<int> dependsOnInputs_;
    };
    struct OutputNode
    {
        std::string name_;
        int         signal_ = 0;

        std::set<int> dependsOnLatches_;
        std::set<int> dependsOnInputs_;
    };
    struct AndNode
    {
        int signal_ = 0;
        int input1_ = 0;
        int input2_ = 0;
    };

    std::vector<InputNode*>  inputs_;
    std::vector<LatchNode*>  latches_;
    std::vector<OutputNode*> outputs_;
    std::vector<AndNode*>    ands_;

    std::vector<std::string> comments_;

    bool computedDependencies_;
};

}  // namespace aigpp

#endif /* AIGPP_AIGERPARSER_HH */

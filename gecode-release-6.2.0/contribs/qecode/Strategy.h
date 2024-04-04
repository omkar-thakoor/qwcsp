#pragma once
/****   , [ bobocheTree.hh ],
 Copyright (c) 2008 Universite d'Orleans - Jeremie Vautard

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *************************************************************************/


#include <vector>
#include <iostream>
#include <string>
#include "StrategyNode.h"
#include "qecode.h"

namespace qecode
{
    class Strategy;

    /** Internal represnetation of a Strategy node tree. Basically a shared pointer
     * TODO : reimplement using modern c++ features
     */
    class StrategyImp
    {
        friend class Strategy;

    private:
        unsigned int mReferencesNumber; ///< number of Strategy objects pointing to this
        StrategyNode mNodeInformation; ///< Node information
        std::vector <Strategy> mChildNodes; ///< children of this node
        unsigned int mTodos; ///< number of child nodes marked as 'to do' in this subtree
        StrategyImp *father; ///< father node (if exists)

        void todosUpdate(int i);

        StrategyImp();

    public:
        explicit StrategyImp(const StrategyNode& tag);

        ~StrategyImp();
    };

/** \brief Strategy of a QCSP+ problem
 *
 * This class represents a solution of a QCSP+. Basically it consists in the tree-representation of the winning strategy.
 * 3 spacial cases exists : the trivially true strategy, the trivially false strategy (used for the UNSAT answer),
 * and the "Dummy" node, used to link together each tree of a strategy which first player is universal (such a strategy is not a tree but a forest)
 */
    class QECODE_VTABLE_EXPORT Strategy
    {
        friend class StrategyImp;

    private:
        StrategyImp *imp;

        void todosUpdate(int i)
        { imp->todosUpdate(i); }

        explicit Strategy(StrategyImp *imp);

    public:
        QECODE_EXPORT Strategy(); ///< default constructor
        QECODE_EXPORT explicit Strategy(StrategyNode tag); ///< builds a strategy on a given strategy node (deprecated)

        /** \brief builds a one node (sub)strategy
         *
         * this method builds a one-node strategy that will typically be attached as child of another strategy. A strategy node embeds informations about quantification,
         * scope and values of variables that must be provided
         * @param qt quantifier of the scope this node represents
         * @param VMin index of the first variable of the scope this node represents
         * @param VMax index of the last variable of the scope this node represents
         * @param values values taken by the variables between VMin and VMax in this part of the strategy
         */
        QECODE_EXPORT Strategy(bool qt, int VMin, int VMax, int scope, std::vector<int> values);

        QECODE_EXPORT Strategy(const Strategy &tree); ///< copy constructor
        QECODE_EXPORT Strategy &operator=(const Strategy &rvalue);

        QECODE_EXPORT ~Strategy();

        /// DEPRECATED returns the StrategyNode object corresponding to the root of this strategy
        QECODE_EXPORT const StrategyNode& getTag() const;
        QECODE_EXPORT forceinline unsigned int degree() const
        { return imp->mChildNodes.size(); }///< returns this strategy's number of children
        QECODE_EXPORT Strategy getChild(int i) const; ///< returns the i-th child of this strategy
        QECODE_EXPORT Strategy getFather();

        QECODE_EXPORT bool hasFather();

        QECODE_EXPORT Strategy
        getSubStrategy(std::vector<int> position); ///< returns the substrategy at given position, dummy if does not exists
        QECODE_EXPORT void
        attach(const Strategy& child);///< attach the strategy given in parameter as a child of the current strategy
        QECODE_EXPORT void detach(unsigned int i); ///< detach the i-th child of this strategy (provided it exists)
        QECODE_EXPORT void detach(const Strategy& son);

        QECODE_EXPORT forceinline bool isFalse() const
        { return imp->mNodeInformation.isFalse(); } ///< returns wether this strategy represents the UNSAT answer
        QECODE_EXPORT forceinline bool isTrue() const
        { return imp->mNodeInformation.isTrue(); } ///< returns wether this strategy is trivially true
        QECODE_EXPORT forceinline bool isComplete() const
        { return ((imp->mTodos) == 0); } ///< returns wether this is a complete sub-strategy (without todo nodes)
        QECODE_EXPORT forceinline bool isDummy() const
        { return imp->mNodeInformation.isDummy(); } ///< returns wether this strategy is a set of
        QECODE_EXPORT forceinline bool isTodo() const
        { return imp->mNodeInformation.isTodo(); } ///< return wether this strategy is a "ToDo" marker or not

        QECODE_EXPORT static Strategy STrue(); ///< returns the trivially true strategy
        QECODE_EXPORT static Strategy SFalse(); ///< returns the trivially false strategy
        QECODE_EXPORT static Strategy Dummy(); ///< returns a "dummy" node
        QECODE_EXPORT static Strategy Stodo(); ///< returns a "todo" node
        QECODE_EXPORT forceinline bool quantifier() const
        { return imp->mNodeInformation.quantifier; } ///< returns the quantifier of the root (true for universal, false for existential)
        QECODE_EXPORT forceinline int VMin() const
        { return imp->mNodeInformation.Vmin; } ///< returns the index of the first variable of the scope of the root
        QECODE_EXPORT forceinline int VMax() const
        { return imp->mNodeInformation.Vmax; } ///< returns the index of the last variable of the scope of the root
        QECODE_EXPORT forceinline int scope() const
        { return imp->mNodeInformation.scope; } ///< returns the scope of the root
        QECODE_EXPORT forceinline const std::vector<int>& values() const
        { return imp->mNodeInformation.valeurs; } ///< returns the values taken by the variables of the scope of the root in this (sub)strategy
        QECODE_EXPORT forceinline int value(int var) const
        { return imp->mNodeInformation.valeurs[var]; }

        QECODE_EXPORT forceinline const void *id() const
        { return imp; } ///< Return an identifier for this strategy (this identifier is shared among multiples instances of this strategy)
        QECODE_EXPORT std::vector<int> getPosition();

        QECODE_EXPORT void getAllValues(std::vector<int> &vals);

        QECODE_EXPORT forceinline int getNbTodos() const
        { return imp->mTodos; }

        QECODE_EXPORT int checkIntegrity();
    };

    void printStrategy(std::ostream& os,const Strategy &s, int depth);
    
    std::ostream& operator<<(std::ostream& os,const Strategy& s);


}
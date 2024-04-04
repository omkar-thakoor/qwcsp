#pragma once
/****   , [ QCOPPlus.hh ],
 Copyright (c) 2009 Universite d'Orleans - Jeremie Vautard

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

#include "QSpace.h"
#include "vartype.h"
#include "OptVar.h"
#include <vector>

namespace qecode
{
    class QECODE_VTABLE_EXPORT QProblem
    {
    private:
        /// Optimization criteria
        struct Opts
        {
            std::vector<OptVar *> vars; // Several Optvars in universal scopes, one in existential ones
            int opt_type;        // Used in existential scopes : 0 for ANY, 1 for MIN, 2 for MAX.
        };
        int *nvar;
        int n;
        int nbSpaces;
        VarType *type_of_v;

        bool *Quantifiers;
        QSpace **rules;
        QSpace *goal;
        Opts *optim;
        int *nbVarBySpace;
        int *whichSpaceOwns;
        bool *varInitialised;

        int currentDeclareSpace;

        QProblem()=default;

    public:
        /// Total number of variables
        forceinline QECODE_EXPORT int getNVariables() const
        { return n; }

        /** \brief  Usual constructor
         *
         *  This is the first step for building a QCSP+/QCOP+ problem. The number of variables and the quantifier of
         *  each rqsets is declared here. However, each variable will have to be defined after that (using
         *  QIntVar or QBoolBar methods).
         *
         *  @param ns number of rqsets in the prefix.
         *  @param quant (Array) quantifier for each rqsets (true: forall, false: exists).
         *  @param nv Array of integer which contains the number of variables by rqset.
         */
        QECODE_EXPORT QProblem(int ns, bool *quant, int *nv);

        QECODE_EXPORT ~QProblem();

        /** \brief Define an integer variable in the quantified space
         *
         *  Defines an integer variable in the quantifies space using a fully declared domain.
         *  @param var Number of the variable to be  defined.
         *  @param dom The initial domain of the variable.
         */
        QECODE_EXPORT void QIntVar(int var, int min, int max);

        /** \brief Defines an integer variable in the quantified space
         *
         *  Defines an integer variable in the quantifies space using a fully declared domain.
         *  @param var Number of the variable to be  defined.
         *  @param dom The initial domain of the variable.
         */
        QECODE_EXPORT void QIntVar(int var, Gecode::IntSet dom);

        /** \brief Defines a boolean variable in the quantified space
         *
         *  Defines a boolean variable with a full initial domain in the quantified space.
         *  @param var Number of the variable to be defined.
         */
        QECODE_EXPORT void QBoolVar(int var);

        /** \brief Returns the current space for constraint declaration
         *
         *  Return the Gecode::space we are currently declaring constraints in.
         *  The nextScope() method should be used to move to the next space once the current one is fully defined.
         */
        QECODE_EXPORT QSpace *space();

        /** \brief Returns an integer variable from the space we are currently declaring
         *
         *  Returns an integer variable from the cpace we are currently declaring. Will abort if the variable is not integer.
         *  @param n The number of the variable to return.
         */
        QECODE_EXPORT Gecode::IntVar var(int n);

        /** \brief Returns a boolean variable from the space we are currently declaring
         *
         * Returns a boolean variable from the space we are currently declaring. Will abort if the variable is not boolean.
         *  @param n The number of the variable to return.
         */
        QECODE_EXPORT Gecode::BoolVar bvar(int n);

        /** \brief Switch to the next space for constraint declaration
         *
         *  Switches to the next space for constraint declaration. var, bvar and space methods will now refer to the
         *  CSP corresponding to the next rqset (or to the goal space if the end of the prefix has been reached).
         *  Returns the new space number, or -1 if it was called while there was no next space to declare constraints in.
         */
        QECODE_EXPORT int nextScope();

        /** \brief Readies the problem for solving
         *
         * This method prepares internal structures for solving. It must be invoked once every variable and constraints
         * have been declared, and before it is given to the solver.
         * Calling this method is not mandatory anymore, although it is recommended to make sure that everything has
         * been well declared (in particular, it will check that each variables are defined.
         * For the existential scopes on which an optimization condition have not been defined yet, this method will
         * post a "Any" optimization condition (i.e. do not optimize).
         */
        QECODE_EXPORT void makeStructure();

        /** \brief returns an aggregate of the problem
         *
         * Creates an aggregate at a given universal scope. This aggregate is an optimization variable that an existential scope can use for optimizing.
         * @param scope the scope where this aggregate will be defined. Must be an unversal scope
         * @param opt The optimization variable we want to aggregate. There must not be any universal scope between an agregate and thisoptimization variable.
         * @param agg The aggregator function that will be used for this aggregate.
         */
        QECODE_EXPORT OptVar *getAggregate(int scope, OptVar *opt, Aggregator *agg);

        /** returns an existential optimization variable
         *
         * Creates an optimization variable that will have the value of an existential variable defined in the problem.
         * @param var the number of the existential variable that must be considered
         */
        QECODE_EXPORT OptVar *getExistential(int var);

        /** set an optimization condition on an existential scope
         *
         * set an optimization condition on a given existential scope of the problem. An optimizaiton condition is composed of an optimization variable (aggregate or existential variable), and of an
         * optimization type.
         * @param scope the scope on which the optimization condition is posted. Must be existential.
         * @param optType the optimization type of the condision we post. 0 is for ANY, 1 is for MIN, 2 is for MAX.
         * @param var the optimization variable to be minimized/maximized
         */
        QECODE_EXPORT void optimize(int scope, int optType, OptVar *var);
        //    QECODE_EXPORT void print();

        QECODE_EXPORT forceinline bool getVariableQuentifier(int v) const
        { return Quantifiers[whichSpaceOwns[v]]; } ///< returns uantifier of variable 'v'
        QECODE_EXPORT forceinline bool getScopeQuantifier(int scope) const
        { return Quantifiers[scope]; } ///< returns quantifier of scope 'scope'
        QECODE_EXPORT int getNbSpaces() const; ///< returns the number of scopes of the problem
        QECODE_EXPORT forceinline int getScopeSize(int scope)
        { return nbVarBySpace[scope]; }///< returns the total number of variables reachable in scope 'scope'
        QECODE_EXPORT QSpace *getSpace(
                int scopeIndex);///< returns a copy of the Gecode::Space corresponding to the given restricted quantifier of the mProblem
        QECODE_EXPORT QSpace *getGoal();

        QECODE_EXPORT int
        getOptType(int scope); ///< returns the optimization type of the given scope of the problem
        QECODE_EXPORT OptVar *
        getOptVar(int scope);///< returns the optimization variable of the 'scope'-th scope of the problem
        QECODE_EXPORT QProblem *clone(); ///< makes a copy of the quantified problem
        QECODE_EXPORT forceinline int getOwnerScope(int var)
        { return whichSpaceOwns[var]; }
    };

}
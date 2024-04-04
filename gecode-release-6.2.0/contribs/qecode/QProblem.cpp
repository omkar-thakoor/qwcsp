/****   , [ QCOPPlus.cc ],
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

#include "QProblem.h"

using namespace Gecode;
using namespace std;
namespace qecode
{
    QProblem::QProblem(int ns, bool *quant, int *nv)
    {
        n = 0;
        this->nvar = nv;
        for (int i = 0; i < ns; i++)
        {
            n += nv[i];
        }
        nbSpaces = ns;
//	v=new void*[n];
        type_of_v = new VarType[n];
        Quantifiers = quant;
        whichSpaceOwns = new int[n];
        nbVarBySpace = new int[nbSpaces];
        rules = new QSpace *[ns];
        nbVarBySpace[0] = nv[0];
        rules[0] = new QSpace(nbVarBySpace[0]);
        for (int i = 1; i < nbSpaces; i++)
        {
            nbVarBySpace[i] = nbVarBySpace[i - 1] + nv[i];
            rules[i] = new QSpace(nbVarBySpace[i]);
        }
        goal = new QSpace(n);
        for (unsigned int i = 0; i < n; i++)
        {
            int lespace = 0;
            while (nbVarBySpace[lespace] <= i) lespace++;
            whichSpaceOwns[i] = lespace;
        }


        varInitialised = new bool[n];
        for (unsigned int i = 0; i < n; i++) varInitialised[i] = false;
        currentDeclareSpace = 0;

        optim = new Opts[nbSpaces];
        for (int i = 0; i < nbSpaces; i++)
            optim[i].opt_type = 0;
    }


    QProblem::~QProblem()
    {
        for (int i = 0; i < nbSpaces; i++)
        {
            delete rules[i];
        }

        delete goal;
        delete[] rules;
    }


    int QProblem::getNbSpaces() const
    {
        return nbSpaces;
    }


    void QProblem::QIntVar(int var, int min, int max)
    {
        if (varInitialised[var])
        {
            cout << "Variable " << var << "  Already created !!" << endl;
            abort();
        }

        for (int i = whichSpaceOwns[var]; i < nbSpaces; i++)
        {
            rules[i]->mVars[var] = new IntVar(*rules[i], min, max);
            rules[i]->mVarTypes[var] = VarType::VTYPE_INT;
        }
        goal->mVars[var] = new IntVar(*goal, min, max);
        goal->mVarTypes[var] = VarType::VTYPE_INT;
        varInitialised[var] = true;
        type_of_v[var] = VarType::VTYPE_INT;
    }


    void QProblem::QIntVar(int var, IntSet dom)
    {
        if (varInitialised[var])
        {
            cout << "Variable " << var << "  Already created !!" << endl;
            abort();
        }

        for (int i = whichSpaceOwns[var]; i < nbSpaces; i++)
        {
            rules[i]->mVars[var] = new IntVar(*rules[i], dom);
            rules[i]->mVarTypes[var] = VarType::VTYPE_INT;
        }
        goal->mVars[var] = new IntVar(*goal, dom);
        goal->mVarTypes[var] = VarType::VTYPE_INT;
        varInitialised[var] = true;
        type_of_v[var] = VarType::VTYPE_INT;
    }


    void QProblem::QBoolVar(int var)
    {
        if (varInitialised[var])
        {
            cout << "Variable " << var << " Already created !!" << endl;
            abort();
        }

        for (int i = whichSpaceOwns[var]; i < nbSpaces; i++)
        {
            rules[i]->mVars[var] = new BoolVar(*rules[i], 0, 1);
            rules[i]->mVarTypes[var] = VarType::VTYPE_BOOL;
        }
        goal->mVars[var] = new BoolVar(*goal, 0, 1);
        goal->mVarTypes[var] = VarType::VTYPE_BOOL;
        varInitialised[var] = true;
        type_of_v[var] = VarType::VTYPE_BOOL;
    }

    QSpace *QProblem::space()
    {
        if (currentDeclareSpace < nbSpaces)
            return rules[currentDeclareSpace];
        if (currentDeclareSpace == nbSpaces)
            return goal;
        return nullptr;
    }


    IntVar QProblem::var(int n)
    {
        if (!varInitialised[n])
        {
            cout << "Variable " << n << " not initialized !" << endl;
            abort();
        }
        if (type_of_v[n] != VarType::VTYPE_INT)
        {
            cout << "Variable " << n << " is not INT" << endl;
            abort();
        }
        return *(static_cast<IntVar *>(space()->mVars[n]));
    }


    BoolVar QProblem::bvar(int n)
    {
        if (!varInitialised[n])
        {
            cout << "Variable " << n << " not initialized !" << endl;
            abort();
        }
        if (type_of_v[n] != VarType::VTYPE_BOOL)
        {
            cout << "Variable " << n << " is not BOOL" << endl;
            abort();
        }
        return *(static_cast<BoolVar *>(space()->mVars[n]));
    }


    int QProblem::nextScope()
    {
        if (currentDeclareSpace == -1) return -1;
        currentDeclareSpace++;
        if (currentDeclareSpace > nbSpaces) return -1;
        return currentDeclareSpace;
    }

    void QProblem::makeStructure()
    {
        for (unsigned int i = 0; i < n; i++)
        {
            if (!varInitialised[i])
            {
                cout << "Can't make structure : variable " << i << " not initialised" << endl;
                abort();
            }
        }

        for (int i = 0; i < nbSpaces; i++)
        {
            if (!Quantifiers[i])
                if (optim[i].vars.empty())
                    optimize(i, 0, getExistential((i == 0) ? 0 : nbVarBySpace[i - 1]));
        }

        for (unsigned int i = 0; i < this->nbSpaces; i++)
        {
            if (rules[i]->status() == SS_FAILED)
            {
                cout << "MakeStructure : rule space " << i << " is already failed." << endl;
            }
        }
        if (goal->status() == SS_FAILED)
        {
            cout << "MakeStructure : goal space is already failed." << endl;
        }

    }


    OptVar *QProblem::getAggregate(int scope, OptVar *opt, Aggregator *agg)
    {
        if (!Quantifiers[scope])
        {
            cout << "Try to get aggregate on existential scope" << endl;
            abort();
        } // aggregateur sur existentiel
        if (opt->getScope() < scope)
        {
            cout << "aggregated variable out of aggregator scope" << endl;
            abort();
        } // Variable aggr�g�e avant aggregateur
        for (int i = scope + 1; i < opt->getScope(); i++) // Universelle entre aggregateur et variable aggr�g�e
            if (Quantifiers[i])
            {
                cout << "Universal scope between variable and aggregator" << endl;
                abort();
            }

        auto *zeopt = new UnivOptVar(scope, opt, agg);
        optim[scope].vars.push_back(zeopt);
        return zeopt;
    }


    OptVar *QProblem::getExistential(int var)
    {
        if (Quantifiers[whichSpaceOwns[var]])
        {
            cout << "can't create existOptVar : variable " << var << " is universal" << endl;
            abort();
        } // Var est universelle
        auto *opt = new ExistOptVar(var, whichSpaceOwns[var]);
        return opt;
    }


    void QProblem::optimize(int scope, int optType, OptVar *var)
    {
        if (var->getScope() < scope)
        {
            cout << "optvar out of optimizing scope" << endl;
            abort();
        } // Variable � optimiser avant d�cideur
        for (int i = scope; i < var->getScope(); i++)
        {
            if (Quantifiers[i])  // universelle entre variable et d�cideur
            {
                cout << "universal scope between variable and optimizing scope" << endl;
                abort();
            }
        }
        //    cout<<"I add an optVar* for scope "<<scope<<endl;
        optim[scope].vars.clear();
        optim[scope].vars.push_back(var);
        optim[scope].opt_type = optType;
    }


    QSpace *QProblem::getSpace(int scope)
    {
        if (scope < 0 || scope > nbSpaces)
        {
//            cout << "I return NULL coz of bad scope value (<0)" << endl;
            return nullptr;
        }
        if (scope == nbSpaces)
        {
            if (goal->status() == SS_FAILED)
            {
//                cout << "I return NULL coz goal is failed" << endl;
                return nullptr;
            }
            //        cout<<"I return the goal"<<endl;
            return static_cast<QSpace *>(goal->clone());
        }
        if (rules[scope]->status() == SS_FAILED)
        {
//            cout << "I return NULL coz scope " << scope << " is failed" << endl;
            return nullptr;
        }
        //    cout<<"I return the rule "<<scope<<endl;
        return (static_cast<QSpace *>(rules[scope]->clone()));
    }


    QSpace *QProblem::getGoal()
    {
        if (goal->status() == SS_FAILED) return nullptr;
        return static_cast<QSpace *>(goal->clone());
    }


    int QProblem::getOptType(int scope)
    {
        if (Quantifiers[scope])
        {
            cout << "Try to get OptType on universal scope" << endl;
            abort();
        }
        return optim[scope].opt_type;
    }


    OptVar *QProblem::getOptVar(int scope)
    {
        if (Quantifiers[scope]) abort();
        //    if (optim[scope].vars.size() == 0) cout<<"no optvar to return, ca va planter...";
        //    cout<<"getoptvar returns optvar* of scope "<<scope;
        return static_cast<OptVar *>(optim[scope].vars[0]);
    }

    QProblem *QProblem::clone()
    {
        bool *cloneQuantifiers = new bool[nbSpaces];
        int *cloneNvar = new int[nbSpaces];
        auto **cloneRules = new QSpace *[nbSpaces];
        Opts *cloneOptims = new Opts[nbSpaces];
        int *cloneNVarBySpace = new int[nbSpaces];
        for (unsigned int i = 0; i < nbSpaces; i++)
        {
            cloneQuantifiers[i] = this->Quantifiers[i];
            cloneNvar[i] = this->nvar[i];
            cloneRules[i] = static_cast<QSpace *>(this->rules[i]->clone());
            cloneOptims[i] = this->optim[i];
            cloneNVarBySpace[i] = this->nbVarBySpace[i];
        }

        //	void** v = new void*[this->n];
        auto *cloneTypeOfV = new VarType[n];
        int *cloneWichSpaceOwns = new int[n];
        for (unsigned int i = 0; i < n; i++)
        {
            cloneTypeOfV[i] = this->type_of_v[i];
            cloneWichSpaceOwns[i] = this->whichSpaceOwns[i];
        }

        QProblem *ret = new QProblem();
        ret->nvar = cloneNvar;
        ret->n = this->n;
        ret->nbSpaces = this->nbSpaces;
        ret->type_of_v = cloneTypeOfV;
        ret->Quantifiers = cloneQuantifiers;
        ret->rules = cloneRules;
        ret->goal = static_cast<QSpace *>(this->goal->clone());
        ret->optim = cloneOptims;
        ret->nbVarBySpace = cloneNVarBySpace;
        ret->whichSpaceOwns = cloneWichSpaceOwns;
        ret->varInitialised = this->varInitialised;
        ret->currentDeclareSpace = this->currentDeclareSpace;

        return ret;
    }
}
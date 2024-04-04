/****   , [ QCSPPlusUnblockable.cc ],
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

#include "QCSPPlusUnblockable.h"
namespace qecode
{

    QcspUnblockable::QcspUnblockable(int ns, bool *quant, int *nv)
    {
        n = 0;
        for (int i = 0; i < ns; i++)
        {
            n += nv[i];
        }
        nbSpaces = ns;
        //    cout<<"QCSPU construction : "<<n<<" vriables"<<endl;
        v = new void *[n];
        type_of_v = new VarType[n];
        Quantifiers = quant;
        whichSpaceOwns = new int[n];
        nbVarBySpace = new int[nbSpaces];
        //    cout<<"Before rule creation : "<<arul<<endl;
        arul = new QSpace(n);
        //    cout<<"after rule creation : "<<arul<<endl;
        nbVarBySpace[0] = nv[0];
        for (int i = 1; i < n; i++)
        {
            nbVarBySpace[i] = nbVarBySpace[i - 1] + nv[i];
        }

        arul = new QSpace(n);
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

        vars = new vector<int>[nbSpaces];
        bvars = new vector<int>[nbSpaces];
    }


    QcspUnblockable::~QcspUnblockable()
    {
        delete arul;
        delete goal;
    }

    int QcspUnblockable::spaces()
    {
        return nbSpaces;
    }

    void QcspUnblockable::QIntVar(int var, int min, int max)
    {
        if (varInitialised[var])
        {
            cout << "Variable " << var << "  Already created !!" << endl;
            abort();
        }
        //    cout<<"Qintvar : arul = "<<arul<<endl;

        arul->mVars[var] = new IntVar(*arul, min, max);
        arul->mVarTypes[var] = VarType::VTYPE_INT;

        goal->mVars[var] = new IntVar(*goal, min, max);
        goal->mVarTypes[var] = VarType::VTYPE_INT;
        varInitialised[var] = true;
        type_of_v[var] = VarType::VTYPE_INT;
    }

    void QcspUnblockable::QIntVar(int var, IntSet dom)
    {
        if (varInitialised[var])
        {
            cout << "Variable " << var << "  Already created !!" << endl;
            abort();
        }

        arul->mVars[var] = new IntVar(*arul, dom);
        arul->mVarTypes[var] = VarType::VTYPE_INT;
        goal->mVars[var] = new IntVar(*goal, dom);
        goal->mVarTypes[var] = VarType::VTYPE_INT;
        varInitialised[var] = true;
        type_of_v[var] = VarType::VTYPE_INT;
    }


    void QcspUnblockable::QBoolVar(int var)
    {
        if (varInitialised[var])
        {
            cout << "Variable " << var << " Already created !!" << endl;
            abort();
        }

        arul->mVars[var] = new BoolVar(*arul, 0, 1);
        arul->mVarTypes[var] = VarType::VTYPE_BOOL;
        goal->mVars[var] = new BoolVar(*goal, 0, 1);
        goal->mVarTypes[var] = VarType::VTYPE_BOOL;
        varInitialised[var] = true;
        type_of_v[var] = VarType::VTYPE_BOOL;
    }

    QSpace *QcspUnblockable::space()
    {
        if (currentDeclareSpace < nbSpaces)
        {
            //cout<<"Return space arul"<<endl; cout.flush();
            return arul;
        }
        if (currentDeclareSpace == nbSpaces)
        {
            //cout<<"Return space Goal"<<endl; cout.flush();
            return goal;
        }
        cout << "Return null in space()" << endl;
        return NULL;
    }


    IntVar QcspUnblockable::var(int n)
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

    BoolVar QcspUnblockable::bvar(int n)
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

    int QcspUnblockable::nextScope()
    {
        if (currentDeclareSpace == -1) return -1;
        currentDeclareSpace++;
        if (currentDeclareSpace > nbSpaces) return -1;
        return currentDeclareSpace;
    }

    void QcspUnblockable::makeStructure()
    {
        for (unsigned int i = 0; i < n; i++)
        {
            if (!varInitialised[i])
            {
                cout << "Can't make structure : variable " << i << " not initialised" << endl;
                abort();
            }
        }
        for (unsigned int i = 0; i < nbSpaces; i++)
        {
            unsigned int nbint = 0;
            unsigned int nbbool = 0;
            for (unsigned int j = 0; j < nbVarBySpace[i]; j++)
            {
                if (type_of_v[j] == VarType::VTYPE_INT)
                    nbint++;
                else
                    nbbool++;
            }
            nbint = 0;
            nbbool = 0;
            for (unsigned int j = 0; j < nbVarBySpace[i]; j++)
            {
                if (type_of_v[j] == VarType::VTYPE_INT)
                {
                    (vars[i]).push_back(j);
                } else
                {
                    (bvars[i]).push_back(j);
                }
            }
        }
    }

    forceinline bool QcspUnblockable::qt_of_var(int v)
    {
        return Quantifiers[whichSpaceOwns[v]];
    }

    QSpace *QcspUnblockable::getSpace(int scope)
    {
        if (scope < 0 || scope > nbSpaces)
        {
            cout << "I return NULL coz of bad scope value (<0)" << endl;
            return NULL;
        }
        if (scope == nbSpaces)
        {
            if (goal->status() == SS_FAILED)
            {
                cout << "I return NULL coz goal is failed" << endl;
                return NULL;
            }
            //        cout<<"I return the goal"<<endl;
            return static_cast<QSpace *>(goal->clone());
        }
        if (arul->status() == SS_FAILED)
        {
            cout << "I return NULL coz scope " << scope << " is failed" << endl;
            return NULL;
        }
        //    cout<<"I return the rule "<<scope<<endl;
        QSpace *ret = (static_cast<QSpace *>(arul->clone()));
        IntVarArgs iva(vars[scope].size());
        BoolVarArgs bva(bvars[scope].size());
        //    cout << "sizes : " <<iva.size() << " " << bva.size()<<endl;
        for (int i = 0; i < iva.size(); i++)
        {
            int idx = (vars[scope])[i];
            iva[i] = *(static_cast<IntVar *>(ret->mVars[idx]));
        }
        for (int i = 0; i < bva.size(); i++)
        {
            int idx = (bvars[scope])[i];
            bva[i] = *(static_cast<BoolVar *>(ret->mVars[idx]));
        }
        br->branch(ret, iva, bva);

        return ret;
    }


    QSpace *QcspUnblockable::getGoal()
    {
        if (goal->status() == SS_FAILED) return NULL;
        return static_cast<QSpace *>(goal->clone());
    }


    void QcspUnblockable::branch(UnblockableBranching *b)
    {
        br = b;
    }
}
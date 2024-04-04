/****   , [ QCSPSolver.cc ],
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

#include "QCSPSolver.h"
#include <climits>

using namespace Gecode;
using namespace std;
namespace qecode
{
    inline vector<int> getTheValues(QSpace* sol,int vmin,int vmax)
    {
        vector<int> zevalues;
        //	cout<<sol<<" "<<vmin<<" "<<vmax<<endl;
        if (vmax > (sol->getNVars())) cout<<"getTheValues mal appele"<<endl;
        for (int i = vmin; i<=vmax; i++) {
            //	cout<<i<<" ";
            //	cout.flush();
            switch (sol->mVarTypes[i]) {
            case VarType::VTYPE_INT :
                zevalues.push_back( (static_cast<IntVar*>(sol->mVars[i]))->val() );
                break;
            case VarType::VTYPE_BOOL :
                zevalues.push_back( (static_cast<BoolVar*>(sol->mVars[i]))->val() );
                break;
            default :
                cout<<"4Unknown variable type"<<endl;
                abort();
            }
        }
        //	cout<<endl;

        return zevalues;
    }

    QCSPSolver::QCSPSolver(QProblem* sp)
    {
        this->sp = sp;
        nbRanges=new int;
    }

    Strategy QCSPSolver::solve(unsigned long int& nodes, unsigned int pLimit, bool allStrategies)
    {
        this->limit=pLimit;
        QSpace* espace=sp->getSpace(0);
        Search::Options o;
        Search::Base<QSpace>* solutions = new DFS<QSpace>(espace,o);
        return rSolve(sp,0,solutions,nodes,allStrategies);
    }



    Strategy QCSPSolver::rSolve(QProblem* qs, int scope, Search::Base<QSpace> *L, unsigned long int& nodes, bool allStrategies)
    {
        nodes++;
        auto* sol = static_cast<QSpace*>(L->next());
        Strategy ret=Strategy::Dummy();
        bool LwasEmpty = true;
        bool atLeastOneExistential = false;
        while ((sol != nullptr) ) {
            LwasEmpty=false;
            vector<int> assignments = getTheValues(sol,0,sol->getNVars()-1);
            Strategy result;

            if (scope == (qs->getNbSpaces() - 1) ) { // last scope reached. Verify the goal...
                QSpace* g = qs->getGoal();
                for (int i=0; i<g->getNVars(); i++) {
                    switch (g->mVarTypes[i]) {
                    case VarType::VTYPE_INT :
                        rel(*g,*(static_cast<IntVar*>(g->mVars[i])) == assignments[i]);
                        break;
                    case VarType::VTYPE_BOOL :
                        rel(*g,*(static_cast<BoolVar*>(g->mVars[i])) == assignments[i]);
                        break;
                    default :
                        cout<<"Unknown variable type"<<endl;
                        abort();
                    }
                }
                Gecode::DFS<QSpace> solutions(g);
                QSpace* goalsol = solutions.next();
                if (goalsol == nullptr) {
                    delete g;
                    result = Strategy::SFalse();
                } else {
                    int vmin = ( (scope==0)? 0 : (qs->getScopeSize(scope - 1)) );
                    int vmax = (qs->getScopeSize(scope)) - 1;
                    vector<int> zevalues=getTheValues(sol,vmin,vmax);
                    result = Strategy::STrue();
                    //                result=Strategy(qs->getScopeQuantifier(scope),vmin,vmax,scope,zevalues);
                    //                result.attach(Strategy::STrue());
                    delete g;
                    //	delete sol;
                    delete goalsol;
                }
            }

            else { // This is not the last scope...
                QSpace* espace = qs->getSpace(scope+1);
                for (int i=0; i<assignments.size(); i++) {
                    switch (espace->mVarTypes[i]) {
                    case VarType::VTYPE_INT :
                        rel(*espace,*(static_cast<IntVar*>(espace->mVars[i])) == assignments[i]);
                        break;
                    case VarType::VTYPE_BOOL :
                        rel(*espace,*(static_cast<BoolVar*>(espace->mVars[i])) == assignments[i]);
                        break;
                    default :
                        cout<<"Unknown variable type"<<endl;
                        abort();
                    }
                }
                Search::Options o;
                Search::Base<QSpace>* solutions = new DFS<QSpace>(espace,o);
                delete espace;
                result=rSolve(qs,scope+1,solutions,nodes,allStrategies);
            }

            int vmin = ( (scope == 0) ? 0 : (qs->getScopeSize(scope - 1)) );
            int vmax = (qs->getScopeSize(scope)) - 1;
            vector<int> zevalues=getTheValues(sol,vmin,vmax);
            delete sol;
            if (qs->getScopeQuantifier(scope)) { // current scope is universal
                if (result.isFalse()) { // one branch fails
                    delete L;
                    return Strategy::SFalse();
                } else {
                    Strategy toAttach(true,vmin,vmax,scope,zevalues);
                    toAttach.attach(result);
                    ret.attach(toAttach);
                }
            } else { //current scope is existential
                if (!result.isFalse()) { // result is not the trivially false strategy...
                    atLeastOneExistential =true;
                    Strategy toAttach;
                    if (allStrategies) {
                        // We want to save every possible strategies. Each correct existential branch will be saved
                        if (scope >= limit) toAttach = Strategy::STrue();
                        else {
                            toAttach = Strategy(qs->getScopeQuantifier(scope), vmin, vmax, scope, zevalues);
                            toAttach.attach(result);
                        }
                        ret.attach(toAttach);
                    } else {
                        //We want only one possible strategy. We found an assignment which leads to a valid substrategy. So, we return it immediately
                        delete L;
                        if (scope >= limit) return Strategy::STrue();
                        ret = Strategy(qs->getScopeQuantifier(scope), vmin, vmax, scope, zevalues);
                        ret.attach(result);
                        return ret;
                    }
                }
            }
            sol = static_cast<QSpace*>(L->next());
        }
        delete L;
        if (scope>limit)
            ret = Strategy::STrue();
        if (qs->getScopeQuantifier(scope))  //universal scope
            return (LwasEmpty ? Strategy::STrue() : ret);
        else // existnetial Scope
            return (atLeastOneExistential ? ret : Strategy::SFalse());
    }
}
/****   , [ Worker.cc ],
 Copyright (c) 2010 Universite de Caen Basse Normandie- Jeremie Vautard

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

#include "Worker.h"

namespace qecode
{

    inline vector<int> getTheValues(QSpace *sol, int vmin, int vmax)
    {
        vector<int> zevalues;
//	cout<<sol<<" "<<vmin<<" "<<vmax<<endl;
        if (vmax > (sol->getNVars())) cout << "getTheValues mal appele" << endl;
        for (int i = vmin; i <= vmax; i++)
        {
            //	cout<<i<<" ";
            //	cout.flush();
            switch (sol->mVarTypes[i])
            {
                case VarType::VTYPE_INT :
                    zevalues.push_back((static_cast<IntVar *>(sol->mVars[i]))->val());
                    break;
                case VarType::VTYPE_BOOL :
                    zevalues.push_back((static_cast<BoolVar *>(sol->mVars[i]))->val());
                    break;
                default :
                    cout << "4Unknown variable type" << endl;
                    abort();
            }
        }
//	cout<<endl;

        return zevalues;
    }


    QWorker::~QWorker()=default;

    void QWorker::stopAndReturn()
    {
        access.acquire();
        if (stopandforget < 1) stopandforget = 1;
        access.release();
    }

    void QWorker::stopAndForget()
    {
        access.acquire();
        if (stopandforget < 2) stopandforget = 2;
        access.release();
    }

    bool QWorker::mustStop()
    {
        return (stopandforget > 0);
    }

    void QWorker::wake()
    {
        goToWork.signal();
    }

    vector<int> QWorker::workPosition()
    {
        return currentWork.root();
    }

    void QWorker::run()
    {

        while (true)
        {
            access.acquire();
            stopandforget = 0;
            finished = false;
            todo.clear();
            access.release();
            QWork cur = wm->getWork(this);
            if (cur.isStop())
            {
                return;
            } else if (cur.isWait())
            {
                goToWork.wait();
            }
            else
            {
                access.acquire();
                currentWork = cur;
                access.release();
                Strategy ret = solve();
                access.acquire();
                bool forget = (stopandforget == 2);
                access.release();
                if (!forget)
                {
                    wm->returnWork(this, ret, todo, cur.root());
                } else
                {
                    for (auto & i : todo)
                    {
                        i.clean();
                    }
                    todo.clear();
                }
            }
        }
    }

    Strategy QWorker::solve()
    {
        return rsolve(currentWork.getScope(), currentWork.getRemaining());
    }


    Strategy QWorker::rsolve(int scope,Search::Base<QSpace> *L)
    {
        access.acquire();
        bool forget = (stopandforget == 2);
        access.release();
        if (forget)
        {
            delete L;
            return Strategy::Dummy();
        }

        QSpace *sol = static_cast<QSpace *>(L->next());
        Strategy ret = Strategy::Dummy();
        bool LwasEmpty = true;


        while ((sol != nullptr))
        {
            LwasEmpty = false;
            vector<int> assignments = getTheValues(sol, 0, sol->getNVars() - 1);
            Strategy result;

            if (scope == (problem->getNbSpaces() - 1))
            { // last scope reached. Verify the goal...
                QSpace *g = problem->getGoal();
                for (int i = 0; i < g->getNVars(); i++)
                {
                    switch (g->mVarTypes[i])
                    {
                        case VarType::VTYPE_INT :
                            rel(*g, *(static_cast<IntVar *>(g->mVars[i])) == assignments[i]);
                            break;
                        case VarType::VTYPE_BOOL :
                            rel(*g, *(static_cast<BoolVar *>(g->mVars[i])) == assignments[i]);
                            break;
                        default :
                            cout << "Unknown variable type" << endl;
                            abort();
                    }
                }
                Gecode::DFS<QSpace> solutions(g);
                QSpace *goalsol = solutions.next();
                if (goalsol == nullptr)
                {
                    delete g;
                    delete L;
                    result = Strategy::SFalse();
                } else
                {
                    int vmin = ((scope == 0) ? 0 : (problem->getScopeSize(scope - 1)));
                    int vmax = (problem->getScopeSize(scope)) - 1;
                    vector<int> zevalues = getTheValues(sol, vmin, vmax);
                    result = Strategy(problem->getScopeQuantifier(scope), vmin, vmax, scope, zevalues);
                    result.attach(Strategy::STrue());
                    delete g;
                    //	delete sol;
                    delete goalsol;
                }
            } else
            { // This is not the last scope...
                QSpace *espace = problem->getSpace(scope + 1);
                for (int i = 0; i < assignments.size(); i++)
                {
                    switch (espace->mVarTypes[i])
                    {
                        case VarType::VTYPE_INT :
                            rel(*espace, *(static_cast<IntVar *>(espace->mVars[i])) == assignments[i]);
                            break;
                        case VarType::VTYPE_BOOL :
                            rel(*espace, *(static_cast<BoolVar *>(espace->mVars[i])) == assignments[i]);
                            break;
                        default :
                            cout << "Unknown variable type" << endl;
                            abort();
                    }
                }
                Options o;
                Search::Base<QSpace> *solutions = new DFS<QSpace>(espace,o);
                delete espace;

                access.acquire();
                forget = (stopandforget == 2);
                bool stop = (stopandforget == 1);
                access.release();
                if (forget)
                {
                    delete sol;
                    delete solutions;
                    return Strategy::Dummy();
                }
                if (stop)
                {
                    vector<int> root1;
                    if (scope != 0)
                    { root1 = getTheValues(sol, 0, problem->getScopeSize(scope - 1) - 1); }
                    QWork current(scope, root1, L);
                    vector<int> root2 = getTheValues(sol, 0, problem->getScopeSize(scope) - 1);
                    QWork onemore(scope + 1, root2, solutions);
                    todo.push_back(current);
                    todo.push_back(onemore);
                    int vmin = ((scope == 0) ? 0 : (problem->getScopeSize(scope - 1)));
                    int vmax = (problem->getScopeSize(scope)) - 1;
                    vector<int> zevalues = getTheValues(sol, vmin, vmax);
                    Strategy toAttach(problem->getScopeQuantifier(scope), vmin, vmax, scope, zevalues);
                    toAttach.attach(Strategy::Stodo());
                    ret.attach(toAttach);
                    ret.attach(Strategy::Stodo());
                    return ret;
                }

                result = rsolve(scope + 1, solutions);
            }

            int vmin = ((scope == 0) ? 0 : (problem->getScopeSize(scope - 1)));
            int vmax = (problem->getScopeSize(scope)) - 1;
            vector<int> zevalues = getTheValues(sol, vmin, vmax);
            delete sol;
            access.acquire();
            forget = (stopandforget == 2);
            access.release();
            if (forget)
            {
                delete L;
                return Strategy::Dummy();
            }
            if (problem->getScopeQuantifier(scope))
            { // current scope is universal
                if (result.isFalse()) // one branch fails
                {
                    delete L;
                    return Strategy::SFalse();
                } else
                {
                    Strategy toAttach(true, vmin, vmax, scope, zevalues);

                    toAttach.attach(result);
                    ret.attach(toAttach);
                }
            } else
            { //current scope is existential
                if (!result.isFalse())
                { // result not the truivilally false strategy...
                    if (result.isComplete())
                    { // ...and has no todo nodes. So...
                        ret = Strategy(problem->getScopeQuantifier(scope), vmin, vmax, scope, zevalues);
                        ret.attach(result);
                        delete L;
                        return ret; // we return it immediately
                    } else
                    { // If result is not false, but still has todo nodes...
                        Strategy toAttach(problem->getScopeQuantifier(scope), vmin, vmax, scope, zevalues);
                        toAttach.attach(result); // the node corresponding to the current scope iis added...
                        ret.attach(
                                toAttach); // and the corresponding strategy is attached to a Dummy node. Next loop in the outer while will need it...
                    }
                }
            }
            sol = static_cast<QSpace *>(L->next());
        }
        delete L;
        if (problem->getScopeQuantifier(scope))  //universal scope
            return (LwasEmpty ? Strategy::STrue() : ret);
        else
        {
            if (ret.isComplete())
                return Strategy::SFalse();
            else
                return ret;
        }
    }
}
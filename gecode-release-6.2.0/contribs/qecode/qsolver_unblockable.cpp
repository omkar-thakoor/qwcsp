/****   , [ QSolverUnblockable.cc ],
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

#include "qsolver_unblockable.h"
namespace qecode
{

    QSolverUnblockable::QSolverUnblockable(QcspUnblockable *sp)
    {
        this->sp = sp;
    }

    Strategy QSolverUnblockable::solve(unsigned long int &nodes)
    {
        vector<int> plop;
        plop.clear();
        return rSolve(sp, 0, plop, nodes);
    }

    Strategy
    QSolverUnblockable::rSolve(QcspUnblockable *qs, int scope, vector<int> assignments, unsigned long int &nodes)
    {
        nodes++;
        //cout<<"rSolve for scope "<<scope<<" with assignments ";
        //    for (int i=0;i<assignments.size();i++) cout<<assignments[i]<<" ";
        //    cout<<endl;
        //////////////////////////////////////////////////////////////////////////////////////////
        // First case : Unblockable QCSP+ : Goal checking.                                      //
        //  If the goal is failed, even before the end of the search, we fail.                  //
        //////////////////////////////////////////////////////////////////////////////////////////
        QSpace *g = qs->getGoal();
        if (g == nullptr)
        { return Strategy::SFalse(); }
        for (int i = 0; i < assignments.size(); i++)
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
                    cout << "1Unknown variable type" << endl;
                    abort();
            }
        }
        if (g->status() == SS_FAILED)
        {
            delete g;
            return Strategy::SFalse();
        }
        delete g;

        if (scope == qs->spaces())
        {
            return Strategy::STrue();
        }

            /////////////////////////////////////////////////////////////////////////////////////////
            // Second case : we are in the middle of the problem...                                //
            /////////////////////////////////////////////////////////////////////////////////////////
        else
        {
            QSpace *espace = qs->getSpace(scope);
            if (espace == nullptr) cout << "I caught a NULL for scope " << scope << ". I will crash..." << endl;
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
                        cout << "2Unknown variable type" << endl;
                        abort();
                }
            }

            // Second case, first subcase : current scope is universal
            /////////////////////////////////////////////////////////
            if (qs->quantification(scope))
            {
                if (espace->status() == SS_FAILED)
                {
                    delete espace;
                    return Strategy::STrue();
                }

                DFS<QSpace> solutions(espace);
                QSpace *sol = solutions.next();
                if (sol == nullptr)
                {
                    delete espace;
                    return Strategy::STrue();
                }

                Strategy retour = Strategy::Dummy();
                while (sol != nullptr)
                {
                    vector<int> assign;
                    for (int i = 0; i < sp->nbVarInScope(scope); i++)
                    {
                        switch (sol->mVarTypes[i])
                        {
                            case VarType::VTYPE_INT :
                                assign.push_back((static_cast<IntVar *>(sol->mVars[i]))->val());
                                break;
                            case VarType::VTYPE_BOOL :
                                assign.push_back((static_cast<BoolVar *>(sol->mVars[i]))->val());
                                break;
                            default :
                                cout << "3Unknown variable type" << endl;
                                abort();
                        }
                    }

                    int vmin = ((scope == 0) ? 0 : (qs->nbVarInScope(scope - 1)));
                    int vmax = (qs->nbVarInScope(scope)) - 1;
                    vector<int> zevalues;
                    for (int i = vmin; i <= vmax; i++)
                    {
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
                    Strategy toAttach(true, vmin, vmax, scope, zevalues);

                    Strategy son = rSolve(qs, scope + 1, assign, nodes);
                    if (son.isFalse())
                    {
                        delete sol;
                        delete espace;
                        return Strategy::SFalse();
                    }
                    toAttach.attach(son);
                    retour.attach(toAttach);
                    delete sol;
                    sol = solutions.next();
                } // end of while
                delete espace;
                return retour;
            } // end of if(universal)

                // Second case, second subcase : current scope is existential
                ////////////////////////////////////////////////////////////
            else
            {
                if ((espace->status()) == SS_FAILED)
                {
                    delete espace;
                    return Strategy::SFalse();
                }

                DFS<QSpace> solutions(espace);
                QSpace *sol = solutions.next();
                if (sol == nullptr)
                {
                    delete espace;
                    return Strategy::SFalse();
                }
                while (sol != nullptr)
                {
                    vector<int> assign;
                    for (int i = 0; i < sp->nbVarInScope(scope); i++)
                    {
                        //                    cout << "i = "<<i<<endl;
                        switch (sol->mVarTypes[i])
                        {
                            case VarType::VTYPE_INT :
                                assign.push_back((static_cast<IntVar *>(sol->mVars[i]))->val());
                                break;
                            case VarType::VTYPE_BOOL :
                                assign.push_back((static_cast<BoolVar *>(sol->mVars[i]))->val());
                                break;
                            default :
                                cout << "5Unknown variable type" << endl;
                                abort();
                        }
                    } // end for

                    int vmin = ((scope == 0) ? 0 : qs->nbVarInScope(scope - 1));
                    int vmax = qs->nbVarInScope(scope) - 1;
                    vector<int> zevalues;
                    for (int i = vmin; i <= vmax; i++)
                    {
                        switch (sol->mVarTypes[i])
                        {
                            case VarType::VTYPE_INT :
                                zevalues.push_back((static_cast<IntVar *>(sol->mVars[i]))->val());
                                break;
                            case VarType::VTYPE_BOOL :
                                zevalues.push_back((static_cast<BoolVar *>(sol->mVars[i]))->val());
                                break;
                            default :
                                cout << "6unknown Variable type" << endl;
                                abort();
                        }
                    }
                    Strategy candidate(false, vmin, vmax, scope, zevalues);

                    Strategy son_of_candidate = rSolve(qs, scope + 1, assign, nodes);
                    if (son_of_candidate.isFalse()) candidate = Strategy::SFalse();
                    else candidate.attach(son_of_candidate);

                    if (!candidate.isFalse())
                    {
                        delete sol;
                        delete espace;
                        return candidate;
                    }
                    delete sol;
                    sol = solutions.next();
                } // end while sol != null
                delete espace;
                return Strategy::SFalse();
            } // end if..else (existential)
        }// end "Second case"
    }


//////////////////////
//////////////////////
//////////////////////

    QSolverUnblockable2::QSolverUnblockable2(QProblem *sp)
    {
        this->sp = sp;
    }

    Strategy QSolverUnblockable2::solve(unsigned long int &nodes)
    {
        vector<int> plop;
        plop.clear();
        return rSolve(sp, 0, plop, nodes);
    }

    Strategy QSolverUnblockable2::rSolve(QProblem *qs, int scope, vector<int> assignments, unsigned long int &nodes)
    {
        nodes++;
        //cout<<"rSolve for scope "<<scope<<" with assignments ";
        //    for (int i=0;i<assignments.size();i++) cout<<assignments[i]<<" ";
        //    cout<<endl;
        //////////////////////////////////////////////////////////////////////////////////////////
        // First case : Unblockable QCSP+ : Goal checking.                                      //
        //  If the goal is failed, even before the end on f the search, we fail.                //
        //////////////////////////////////////////////////////////////////////////////////////////
        QSpace *g = qs->getGoal();
        if (g == nullptr)
        { return Strategy::SFalse(); }
        for (int i = 0; i < assignments.size(); i++)
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
                    cout << "1Unknown variable type" << endl;
                    abort();
            }
        }
        if (g->status() == SS_FAILED)
        {
            delete g;
            return Strategy::SFalse();
        }
        delete g;

        if (scope == qs->getNbSpaces())
        {
            return Strategy::STrue();
        }

            /////////////////////////////////////////////////////////////////////////////////////////
            // Second case : we are in the middle of the problem...                                //
            /////////////////////////////////////////////////////////////////////////////////////////
        else
        {
            QSpace *espace = qs->getSpace(scope);
            if (espace == nullptr) cout << "I caught a NULL for scope " << scope << ". I will crash..." << endl;
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
                        cout << "2Unknown variable type" << endl;
                        abort();
                }
            }

            // Second case, first subcase : current scope is universal
            /////////////////////////////////////////////////////////
            if (qs->getScopeQuantifier(scope))
            {
                if (espace->status() == SS_FAILED)
                {
                    delete espace;
                    return Strategy::STrue();
                }

                DFS<QSpace> solutions(espace);
                QSpace *sol = solutions.next();
                if (sol == nullptr)
                {
                    delete espace;
                    return Strategy::STrue();
                }

                Strategy retour = Strategy::Dummy();
                while (sol != nullptr)
                {
                    vector<int> assign;
                    for (int i = 0; i < sp->getScopeSize(scope); i++)
                    {
                        switch (sol->mVarTypes[i])
                        {
                            case VarType::VTYPE_INT :
                                assign.push_back((static_cast<IntVar *>(sol->mVars[i]))->val());
                                break;
                            case VarType::VTYPE_BOOL :
                                assign.push_back((static_cast<BoolVar *>(sol->mVars[i]))->val());
                                break;
                            default :
                                cout << "3Unknown variable type" << endl;
                                abort();
                        }
                    }

                    int vmin = ((scope == 0) ? 0 : (qs->getScopeSize(scope - 1)));
                    int vmax = (qs->getScopeSize(scope)) - 1;
                    vector<int> zevalues;
                    for (int i = vmin; i <= vmax; i++)
                    {
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
                    Strategy toAttach(true, vmin, vmax, scope, zevalues);

                    Strategy son = rSolve(qs, scope + 1, assign, nodes);
                    if (son.isFalse())
                    {
                        delete sol;
                        delete espace;
                        return Strategy::SFalse();
                    }
                    toAttach.attach(son);
                    retour.attach(toAttach);
                    delete sol;
                    sol = solutions.next();
                } // end of while
                delete espace;
                return retour;
            } // end of if(universal)

                // Second case, second subcase : current scope is existential
                ////////////////////////////////////////////////////////////
            else
            {
                if ((espace->status()) == SS_FAILED)
                {
                    delete espace;
                    return Strategy::SFalse();
                }

                DFS<QSpace> solutions(espace);
                QSpace *sol = solutions.next();
                if (sol == nullptr)
                {
                    delete espace;
                    return Strategy::SFalse();
                }
                while (sol != nullptr)
                {
                    vector<int> assign;
                    for (int i = 0; i < sp->getScopeSize(scope); i++)
                    {
                        //                    cout << "i = "<<i<<endl;
                        switch (sol->mVarTypes[i])
                        {
                            case VarType::VTYPE_INT :
                                assign.push_back((static_cast<IntVar *>(sol->mVars[i]))->val());
                                break;
                            case VarType::VTYPE_BOOL :
                                assign.push_back((static_cast<BoolVar *>(sol->mVars[i]))->val());
                                break;
                            default :
                                cout << "5Unknown variable type" << endl;
                                abort();
                        }
                    } // end for

                    int vmin = ((scope == 0) ? 0 : qs->getScopeSize(scope - 1));
                    int vmax = qs->getScopeSize(scope) - 1;
                    vector<int> zevalues;
                    for (int i = vmin; i <= vmax; i++)
                    {
                        switch (sol->mVarTypes[i])
                        {
                            case VarType::VTYPE_INT :
                                zevalues.push_back((static_cast<IntVar *>(sol->mVars[i]))->val());
                                break;
                            case VarType::VTYPE_BOOL :
                                zevalues.push_back((static_cast<BoolVar *>(sol->mVars[i]))->val());
                                break;
                            default :
                                cout << "6unknown Variable type" << endl;
                                abort();
                        }
                    }
                    Strategy candidate(false, vmin, vmax, scope, zevalues);

                    Strategy son_of_candidate = rSolve(qs, scope + 1, assign, nodes);
                    if (son_of_candidate.isFalse()) candidate = Strategy::SFalse();
                    else candidate.attach(son_of_candidate);

                    if (!candidate.isFalse())
                    {
                        delete sol;
                        delete espace;
                        return candidate;
                    }
                    delete sol;
                    sol = solutions.next();
                } // end while sol != null
                delete espace;
                return Strategy::SFalse();
            } // end if..else (existential)
        }// end "Second case"
    }
}
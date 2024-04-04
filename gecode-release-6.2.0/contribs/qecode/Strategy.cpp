/****   , [ bobocheTree.cc ],
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

#include "Strategy.h"
using namespace std;
namespace qecode
{
    StrategyImp::StrategyImp()
    {
        cout << "Default constructor of StrategyImp should not be called !" << endl;
        mReferencesNumber = 1;
        mNodeInformation = StrategyNode::Dummy();
        mTodos = 0;
        father = nullptr;
    }

    void StrategyImp::todosUpdate(int i)
    {
        mTodos += i;
        if (father != nullptr) father->todosUpdate(i);
    }

    StrategyImp::StrategyImp(const StrategyNode& tag)
    {
        //    cout<<"Strategy imp constructor"<<endl;
        mReferencesNumber = 1;
        mNodeInformation = tag;
        mTodos = 0;
        father = nullptr;
        //    cout<<"strategyimp constructor fini"<<endl;
    }


    StrategyImp::~StrategyImp()
    {
        for (auto & node : mChildNodes)
        {
            if ((node.imp->father) == this)
            {
                node.imp->father = nullptr;
            }
        }
    }


    Strategy::Strategy()
    {
        //    cout<<"strategy default"<<endl;
        StrategyNode tag = StrategyNode::Dummy();
        imp = new StrategyImp(tag);
    }

    Strategy::Strategy(StrategyNode tag)
    {
        //    cout<<"Strategy with tag"<<endl;
        imp = new StrategyImp(tag);
        //cout<<"passed imp creation. End of strategy creator"<<endl;
    }

    Strategy::Strategy(StrategyImp *imp)
    {
        this->imp = imp;
        this->imp->mReferencesNumber++;
    }

    Strategy::Strategy(bool qt, int VMin, int VMax, int scope, vector<int> values)
    {
        //    cout<<"strategy with values"<<endl;
        StrategyNode tag(2, qt, VMin, VMax, scope);
        tag.valeurs = values;
        imp = new StrategyImp(tag);

    }


    Strategy::Strategy(const Strategy &tree)
    {
        //    cout<<"Strategy copy"<<endl;
        imp = tree.imp;
        (imp->mReferencesNumber)++;
    }


    Strategy &Strategy::operator=(const Strategy &rvalue)
    {
        //    cout<<"Strategy = "<<endl;
        if (imp != NULL)
        {
            (imp->mReferencesNumber)--;
            if ((imp->mReferencesNumber) == 0)
            {
                //        cout<<"no more references for the imp. Delete"<<endl;
                delete imp;
            }
        }
        imp = rvalue.imp;
        (imp->mReferencesNumber)++;
        return *this;
    }


    Strategy::~Strategy()
    {
        //    cout<<"strategy destructor"<<endl;
        (imp->mReferencesNumber)--;
        if ((imp->mReferencesNumber) == 0)
        {
            //        cout<<"no more references for the imp. Delete"<<endl;
            delete imp;
        }
    }


    const StrategyNode &Strategy::getTag() const
    {
        return imp->mNodeInformation;
    }

    Strategy Strategy::getFather()
    {
        if (hasFather()) return Strategy(imp->father);
        return Dummy();
    }

    bool Strategy::hasFather()
    {
        if (imp->father != NULL)
        {
            for (int i = 0; i < imp->father->mChildNodes.size(); i++)
            {
                if ((imp->father->mChildNodes[i].imp) == imp)
                    return true;
            }
        }
        return false;
    }

    Strategy Strategy::getChild(int i) const
    {
        if (i < 0 || i >= degree())
        {
            cout << "Child " << i << " does not exist" << endl;
            abort();
        }
        return imp->mChildNodes[i];
    }

    Strategy Strategy::getSubStrategy(vector<int> position)
    {
        if (position.empty()) return *this;
        int deg = degree();
        if (deg == 0)
        {
            cout << "Did not find substrategy" << endl;
            return Strategy::Dummy();
        }
        for (int i = 0; i < deg; i++)
        {
            Strategy child = getChild(i);
            bool ok = true;
            if (child.values().size() == 0)
            {
                ok = false;
            }
            for (int j = 0; (j < child.values().size()) && ok; j++)
            {
                if (child.value(j) != position[j]) ok = false;
            }
            if (ok)
            {
                position.erase(position.begin(), position.begin() + (child.values().size()));
                return child.getSubStrategy(position);
            }
        }
        cout << "Did not find substrategy" << endl;
        return Strategy::Dummy();
    }

    void Strategy::attach(const Strategy& child)
    {
        if (child.isDummy())
        {
            int todosToAdd = 0;

            for (int i = 0; i < child.degree(); i++)
            {
                this->attach(child.getChild(i));
            }
        } else
        {
            imp->mChildNodes.push_back(child);
            todosUpdate(child.imp->mTodos);
            (child.imp)->father = this->imp;
        }
    }

    void Strategy::detach(const Strategy& son)
    {

        auto it = imp->mChildNodes.begin();
        while (it != (imp->mChildNodes.end()) && ((*it).id() != son.id()))
        {
            it++;
        }
        if (it != imp->mChildNodes.end())
        {
            todosUpdate(0 - ((*it).imp->mTodos));
            (*it).imp->father = nullptr;
            imp->mChildNodes.erase(it);
        }
    }


    void Strategy::detach(unsigned int i)
    {
        if (imp->mChildNodes.size() < i) return;

        auto it = imp->mChildNodes.begin() + i;
        todosUpdate(0 - ((*it).imp->mTodos));
        (*it).imp->father = nullptr;
        imp->mChildNodes.erase(it);
    }

    Strategy Strategy::STrue()
    {
        Strategy ret(StrategyNode::STrue());
        return ret;
    }

    Strategy Strategy::SFalse()
    {
        Strategy ret(StrategyNode::SFalse());
        return ret;
    }

    Strategy Strategy::Dummy()
    {
        Strategy ret(StrategyNode::Dummy());
        return ret;
    }

    Strategy Strategy::Stodo()
    {
        Strategy ret(StrategyNode::Todo());
        ret.imp->mTodos = 1;
        return ret;
    }

    vector<int> Strategy::getPosition()
    {
        vector<int> ret;
        Strategy asc = *this;
        while (!asc.isDummy())
        {
//		cout<<"GetPosition adding "<<asc.values().size()<<" elements to a vector of size "<<ret.size()<<endl;
            vector<int> ret2;
            ret2.reserve(ret.size() + asc.values().size() + 1);
            for (int i = 0; i < asc.values().size(); i++)
            {
                ret2.push_back(asc.values()[i]);
            }
            for (int i = 0; i < ret.size(); i++)
            {
                ret2.push_back(ret[i]);
            }
            ret = ret2;
            asc = asc.getFather();
        }
        return ret;
    }

    int Strategy::checkIntegrity()
    {
        int ret = 0;
        for (unsigned int i = 0; i < (this->degree()); i++)
        {
            if ((((imp->mChildNodes[i]).imp)->father) != imp)
            {
                ret++;
                cout << (((imp->mChildNodes[i]).imp)->father) << " should be " << imp << endl;
            }
            ret += (getChild(i).checkIntegrity());
        }
        return ret;
    }

    void Strategy::getAllValues(std::vector<int> &vals) {
        vals.insert(vals.end(), getTag().valeurs.begin(), getTag().valeurs.end());
        for (int i=0;i<degree();i++) getChild(i).getAllValues(vals);
    }


    void printStrategy(ostream& os,const Strategy &s, int depth)
    {
        StrategyNode node = s.getTag();
        os << std::string(depth, ' ');
        if (s.isTrue())
            os<<"TRUE"<<endl;
        else if (s.isFalse())
            os<<"FALSE"<<endl;
        else
            os << "type " << node.type
               << " qt " << node.quantifier
               << " vmin " << node.Vmin
               << " vmax " << node.Vmax
               << " scope " << node.scope
               << " - ";
        for (int valeur : s.getTag().valeurs)
            os<<valeur<<" ";
        os<<endl;
        os << std::string(depth, ' ');
        os<<s.degree()<<" child(ren)"<<endl;
        for (int i=0;i<s.degree();i++)
        {
            os << std::string(depth, ' ');
            os<<"Child "<<i<<" : "<<endl;
            printStrategy(os,s.getChild(i), depth + 1);
        }
    }
    
    ostream &operator<<(ostream &os, const Strategy &s)
    {
        printStrategy(os,s,0);
        return os;
    }
}
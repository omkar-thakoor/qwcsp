/*****************************************************************[myspace.cc]
Copyright (c) 2007, Universite d'Orleans - Jeremie Vautard, Marco Benedetti,
Arnaud Lallouet.

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
*****************************************************************************/
#include "QSpace.h"
#include <iostream>


using namespace std;
using namespace Gecode;

namespace qecode
{
    QSpace::QSpace(unsigned int nv)
    {
        //  cout <<"Space with "<<getNVariables<<" variables"<<endl;
        mNbVars = nv;
        mVars = new void *[nv];
        mVarTypes = new VarType[nv];
    }


    QSpace::~QSpace()
    {
        for (int i = 0; i < mNbVars; i++)
            switch (mVarTypes[i])
            {
                case VarType::VTYPE_INT:
                    delete static_cast<IntVar *>(mVars[i]);
                    break;
                case VarType::VTYPE_BOOL:
                    delete static_cast<BoolVar *>(mVars[i]);
                    break;
                default:
                    cout << "Unsupported variable type" << endl;
                    abort();
            }

        delete[] mVars;
        delete[] mVarTypes;
    }


    QSpace::QSpace(QSpace& ms) : Space( ms)
    {
        mNbVars = ms.mNbVars;
        mVars = new void *[mNbVars];
        mVarTypes = new VarType[mNbVars];
        for (int i = 0; i < mNbVars; i++)
        {
            mVarTypes[i] = ms.mVarTypes[i];
            switch (mVarTypes[i])
            {
                case VarType::VTYPE_INT:
                    mVars[i] = new IntVar(*(static_cast<IntVar *>(ms.mVars[i])));
                    (static_cast<IntVar *>(mVars[i]))->update(*this, *(static_cast<IntVar *>(ms.mVars[i])));
                    break;
                case VarType::VTYPE_BOOL:
                    mVars[i] = new BoolVar(*(static_cast<BoolVar *>(ms.mVars[i])));
                    (static_cast<BoolVar *>(mVars[i]))->update(*this, *(static_cast<BoolVar *>(ms.mVars[i])));
                    break;
                default:
                    cout << "Unsupported variable type" << endl;
                    abort();
            }
        }
    }


    QSpace* QSpace::copy() { return new QSpace( *this); }

    IntVarArgs QSpace::getIntVars(unsigned int idMax)
    {
        int cpt = 0;
        int i = 0;
        if (mNbVars < idMax) idMax = mNbVars;

        for (int i = 0; i < idMax; i++)
        {
            if (mVarTypes[i] == VarType::VTYPE_INT) cpt++;
        }
        IntVarArgs ret(cpt);
        cpt = 0;
        for (i = 0; i < idMax; i++)
        {
            if (mVarTypes[i] == VarType::VTYPE_INT)
            {
                ret[cpt] = *(static_cast<IntVar *>(mVars[i]));
                cpt++;
            }
        }

        return ret;
    }

    BoolVarArgs QSpace::getBoolVars(unsigned int idMax)
    {
        int cpt = 0;
        if (mNbVars < idMax) idMax = mNbVars;

        for (int i = 0; i < idMax; i++)
        {
            if (mVarTypes[i] == VarType::VTYPE_BOOL) cpt++;
        }
        BoolVarArgs ret(cpt);
        cpt = 0;
        for (int i = 0; i < idMax; i++)
        {
            if (mVarTypes[i] == VarType::VTYPE_BOOL)
            {
                ret[cpt] = *(static_cast<BoolVar *>(mVars[i]));
                cpt++;
            }
        }

        return ret;
    }
}

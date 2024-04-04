#pragma once
/************************************************************ Work.hh
Copyright (c) 2010 Universite de Caen Basse Normandie - Jeremie Vautard

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
#include<vector>
#include "Strategy.h"
#include "qecode.h"
#include "gecode/support.hh"
#include "gecode/search.hh"
#include "gecode/search/support.hh"
#include "QSpace.h"

using namespace std;
using namespace Gecode;
using namespace Gecode::Support;
namespace qecode
{
    class QECODE_VTABLE_EXPORT QWork
    {
    private:
        bool wait;
        bool stop;
        vector<int> theRoot;
        Gecode::Search::Base<QSpace>*remaining;
        int scope;
    public:
        QECODE_EXPORT QWork()
        {
            wait = stop = false;
            remaining = nullptr;
            scope = -1;
        } ///< creates a dummy work. Don't try to solve it.
        QECODE_EXPORT QWork(int scope, vector<int> root, Gecode::Search::Base<QSpace> *todo);

        QECODE_EXPORT ~QWork();

        QECODE_EXPORT static QWork Wait();

        QECODE_EXPORT static QWork Stop();

        QECODE_EXPORT forceinline bool isWait() const
        { return wait; }

        QECODE_EXPORT forceinline bool isStop() const
        { return stop; }

        QECODE_EXPORT forceinline Gecode::Search::Base<QSpace> *getRemaining()
        { return this->remaining; }

        QECODE_EXPORT const vector<int>& root() const
        {
            return theRoot;
        }

        QECODE_EXPORT forceinline int getScope() const
        { return this->scope; }

        QECODE_EXPORT forceinline void clean()
        { if (!wait && !stop) delete remaining; }
    };
}
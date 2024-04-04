#pragma once
/****   , [ Worker.hh ],
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

#include "qecode.h"
#include<list>
#include<vector>
#include <ctime>
#include "QProblem.h"
#include "QSpace.h"
#include "Strategy.h"
#include "AbstractWorker.h"
#include "gecode/support.hh"
#include "gecode/search.hh"
#include "gecode/search/support.hh"
#include "WorkManager.h"

using namespace Gecode::Support;
using namespace Gecode::Search;
using namespace Gecode::Search::Sequential;

namespace qecode
{

    class QECODE_VTABLE_EXPORT QWorker : public AQWorker
    {

    private:
        WorkManager *wm;
        Gecode::Support::Event goToWork;
        Gecode::Support::Mutex access;
        int stopandforget; // 0 : continue, 1 : stop and return, 2 : stop and forget.
        QProblem *problem;
        QWork currentWork;
        list<QWork> todo;
        bool finished;

        Strategy rsolve(int scope,Gecode::Search::Base<QSpace> *L);


    public:
        QECODE_EXPORT explicit QWorker(WorkManager *wm)
        {
            this->wm = wm;
            this->problem = wm->problem->clone();
        }

        QECODE_EXPORT ~QWorker() override;

        QECODE_EXPORT void run() override;

        QECODE_EXPORT void stopAndReturn() override;

        QECODE_EXPORT void stopAndForget() override;

        QECODE_EXPORT vector<int> workPosition() override;

        QECODE_EXPORT bool mustStop() override;

        QECODE_EXPORT void wake() override;

        QECODE_EXPORT    Strategy solve();
    };
}
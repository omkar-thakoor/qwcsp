/************************************************************ Work.cc
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

#include "Work.h"

#include <utility>
namespace qecode
{

    QWork::QWork(int scope, vector<int> root, Gecode::Search::Base<QSpace> *todo)
    {
        this->theRoot = std::move(root);
        this->remaining = todo;
        this->scope = scope;
        this->wait = this->stop = false;
    }

    QWork::~QWork()= default;

    QWork QWork::Wait()
    {
        QWork ret;
        ret.wait = true;
        ret.stop = false;
        ret.remaining = nullptr;
        return ret;
    }

    QWork QWork::Stop()
    {
        QWork ret;
        ret.wait = false;
        ret.stop = true;
        ret.remaining = nullptr;
        return ret;
    }
}
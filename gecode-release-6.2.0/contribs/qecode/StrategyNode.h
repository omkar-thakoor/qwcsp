#pragma once
/****   , [ StrategyNode.hh ],
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
 #include "qecode.h"
#include <vector>
namespace qecode
{
    class QECODE_VTABLE_EXPORT StrategyNode
    {
    public:
        int type; // 0 = dummy, 1 = empty, 2 = normal, 3 = Todo.
        bool quantifier;
        int Vmin;
        int Vmax;
        int scope;
        std::vector<int> valeurs;

        QECODE_EXPORT StrategyNode();

        QECODE_EXPORT StrategyNode(int type, bool qt, int Vmin, int Vmax, int scope);

        QECODE_EXPORT ~StrategyNode();

        QECODE_EXPORT static StrategyNode STrue()
        { return {1, true, -1, -1, -1}; }

        QECODE_EXPORT static StrategyNode SFalse()
        { return {1, false, -1, -1, -1}; }

        QECODE_EXPORT static StrategyNode Dummy()
        { return {0, true, -1, -1, -1}; }

        QECODE_EXPORT static StrategyNode Todo()
        { return {3, false, -1, -1, -1}; }

        QECODE_EXPORT forceinline bool isFalse() const
        { return ((type == 1) && !quantifier); }

        QECODE_EXPORT forceinline bool isTrue() const
        { return ((type == 1) && quantifier); }

        QECODE_EXPORT forceinline bool isDummy() const
        { return (type == 0); }

        QECODE_EXPORT forceinline bool isTodo() const
        { return (type == 3); }
    };
}

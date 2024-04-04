#pragma once
/*****************************************************************[myspace.hh]
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

#include "qecode.h"
#include "gecode/minimodel.hh"
#include "vartype.h"

/** \brief A simple extension of Gecode::Space class
 *
 *  A simple extension of the Space class from Gecode, in order to have access to the variables it contains.
 */

namespace qecode
{
    class QECODE_VTABLE_EXPORT QSpace : public Gecode::Space
    {
    protected :
        /// NUmber of variables (size of mNbVars array)
        unsigned int mNbVars;

    public:
        /// This space's variables
        void** mVars;
        /// Variables types (intvar or boolVar)
        VarType* mVarTypes;

        /** \brief Constructor of a space with a fixed number of variables
         *
         * Builds a space which will contain  getNVariables variables (the variables themselves are however not declared).
         *  @param nv the number of variable the space must contain.
         */
        QECODE_EXPORT explicit QSpace(unsigned int nv);

        QECODE_EXPORT unsigned int getNVars() const { return mNbVars; }
        QECODE_EXPORT QSpace(QSpace& ms);

        QECODE_EXPORT QSpace* copy() override;

        QECODE_EXPORT~QSpace() override;

        //	QECODE_EXPORT int getValue(unsigned int i); ///< returns the value of variable i. If boolean : 0 or 1 (false / true).

        /// \brief Get all integer variables with index < idMax
        QECODE_EXPORT Gecode::IntVarArgs getIntVars(unsigned int idMax);

        /// \brief Get all boolean variables with index < idMax
        QECODE_EXPORT Gecode::BoolVarArgs getBoolVars(unsigned int idMax);
    };
}

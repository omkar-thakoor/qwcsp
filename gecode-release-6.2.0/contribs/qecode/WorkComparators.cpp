#include "WorkComparators.h"

bool qecode::QuantifierThenDepthComparator::cmp(qecode::QWork a, qecode::QWork b)
{
    bool q1 = (mProblem->getVariableQuentifier(a.root().size()) != mExistsFirst);
    bool q2 = (mProblem->getVariableQuentifier(b.root().size()) != mExistsFirst);
    if (q1 && !q2) return true;
    if (!q1 && q2) return false;
    int d1 = a.root().size();
    int d2 = b.root().size();
    if ((d1 < d2) != mDeepestFirst) return true;
    return false;
}

qecode::QuantifierThenDepthComparator::QuantifierThenDepthComparator(qecode::QProblem *p, bool existsFirst,
                                                                     bool deepestFirst) :
        mProblem(p),mExistsFirst(existsFirst),mDeepestFirst(deepestFirst)
{}

bool qecode::DepthComparator::cmp(qecode::QWork a, qecode::QWork b)
{
    int d1 = a.root().size();
    int d2 = b.root().size();
    if ((d1 < d2) != mDeepestFirst) return true;
    return false;
}

bool qecode::DeepFirstComparator::cmp(qecode::QWork a, qecode::QWork b)
{
    if (a.root().size() > b.root().size()) return true;
    if (a.root().size() < b.root().size()) return false;
    return ((a.getRemaining()) < (b.getRemaining()));
}

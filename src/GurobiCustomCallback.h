#ifndef GUROBICUSTOMCALLBACK_H_
#define GUROBICUSTOMCALLBACK_H_

#include <gurobi_c++.h>
#include <limits>
#include <fstream>

#include "RunningTime.h"
#include "UpdateHistory.h"
#include "Recorder.h"

/** \brief This callback records the objective function value when a new solution is found by the MIP solver */
class GurobiCustomCallbackType0 : public GRBCallback
{
private:

    LPObjUpdateHistory *updateHistory;

    double objbst = std::numeric_limits<double>::max();
    double objbnd = 0;

public:
    GurobiCustomCallbackType0(LPObjUpdateHistory *updateHistory) {
        this->updateHistory = updateHistory;
    }

protected:
    void callback()
    {
        try
        {
            if (where == GRB_CB_MIPSOL)
            {
                // MIP solution callback
                double new_bnd = getDoubleInfo(GRB_CB_MIPSOL_OBJBND);
                if (new_bnd > objbnd)
                    objbnd = new_bnd;
                objbst = getDoubleInfo(GRB_CB_MIPSOL_OBJBST);
                auto ttime = RunningTime::GetInstance().getCurrentRunningTime();
                updateHistory->push_back_bounds(objbnd, objbst, ttime);
            }
            else if (where == GRB_CB_MIP)
            {
                // MIP callback
                double new_bnd = getDoubleInfo(GRB_CB_MIP_OBJBND);
                if (new_bnd > objbnd)
                {
                    auto ttime = RunningTime::GetInstance().getCurrentRunningTime();
                    objbnd = new_bnd;
                    objbst = getDoubleInfo(GRB_CB_MIP_OBJBST);
                    updateHistory->push_back_bounds(objbnd, objbst, ttime);
                }
            }
        }
        catch (GRBException e)
        {
            std::cerr << "Error number: " << e.getErrorCode() << std::endl;
            std::cerr << e.getMessage() << std::endl;
        }
        catch (...)
        {
            std::cerr << "Error during callback" << std::endl;
        }
    }
};

/** This callback records the assignments of the variable when a new solution is found by the MIP solver */
class GurobiCustomCallbackType1 : public GRBCallback
{
private:
    LPAssignmentsUpdateHistory *updateHistory;

    GRBVar *xvars; // internal Gurobi variables
    int len;

public:
    GurobiCustomCallbackType1(
        LPAssignmentsUpdateHistory *updateHistory,
        GRBVar *xvars,
        int len)
    {
        this->updateHistory = updateHistory;
        this->xvars = xvars;
        this->len = len;
    }

protected:
    void callback()
    {
        try
        {
            if (where == GRB_CB_MIPSOL)
            {
                // MIP solution callback
                auto ttime = RunningTime::GetInstance().getCurrentRunningTime();
                auto solution = getSolution(xvars, len);
                std::vector<double> assignments;
                for (unsigned i = 0; i < len; i++)
                    assignments.push_back(solution[i]);
                updateHistory->push_back_lp_assignments(assignments, ttime);
            }
        }
        catch (GRBException e)
        {
            std::cerr << "Error number: " << e.getErrorCode() << std::endl;
            std::cerr << e.getMessage() << std::endl;
        }
        catch (...)
        {
            std::cerr << "Error during callback" << std::endl;
        }
    }
};

#endif // GUROBICUSTOMCALLBACK_H_
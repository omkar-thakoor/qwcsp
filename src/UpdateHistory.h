#ifndef UPDATEHISTORY_H_
#define UPDATEHISTORY_H_

#include <vector>

enum update_t {
    numerical_bounds,
    assignments
};

class UpdateHistory
{
protected:
    int iter_idx = -1;
    std::vector<double> timeHistory;

public:
    virtual update_t update_type() = 0;

    bool has_next_update()
    {
        return iter_idx + 1 < timeHistory.size();
    }

    double next_update()
    {
        iter_idx++;
        return timeHistory[iter_idx];
    }
};

class LPObjUpdateHistory : public UpdateHistory
{
private:
    std::vector<double> lbHistory;
    std::vector<double> ubHistory;

public:
    update_t update_type() { return numerical_bounds; } 

    void push_back_bounds(double lb, double ub, double ttime)
    {
        lbHistory.push_back(lb);
        ubHistory.push_back(ub);
        timeHistory.push_back(ttime);
    }

    double getLB() { return lbHistory[iter_idx]; }

    double getUB() { return ubHistory[iter_idx]; }
};

class AssignmentsUpdateHistory : public UpdateHistory
{
public:
    update_t update_type() { return assignments; }
    virtual bool assignmentAt(int idx) = 0;
    virtual bool bestAssignmentAt(int idx) = 0;
};

class LPAssignmentsUpdateHistory : public AssignmentsUpdateHistory
{
private:
    std::vector<std::vector<double>> lpAssignmentsUpdateHistory;

public:
    void push_back_lp_assignments(std::vector<double> &assignments, double ttime)
    {
        lpAssignmentsUpdateHistory.push_back(assignments);
        timeHistory.push_back(ttime);
    }

    bool assignmentAt(int idx)
    {
        return lpAssignmentsUpdateHistory[iter_idx][idx] > 0.5;
    }

    bool bestAssignmentAt(int idx)
    {
        return lpAssignmentsUpdateHistory[lpAssignmentsUpdateHistory.size() - 1][idx] > 0.5;
    }
};

class VCUpdateHistory : public AssignmentsUpdateHistory
{
private:
    std::vector<std::vector<int>> vcHistory;

public:
    void emplace_back_vc(int *begin, int *end, double ttime)
    {
        vcHistory.emplace_back(begin, end);
        timeHistory.push_back(ttime);
    }

    bool assignmentAt(int idx)
    {
        return vcHistory[iter_idx][idx] == 1;
    }

    bool bestAssignmentAt(int idx)
    {
        return vcHistory[vcHistory.size() - 1][idx] == 1;
    }
};

#endif // UPDATEHISTORY_H_
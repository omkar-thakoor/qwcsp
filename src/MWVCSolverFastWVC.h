#ifndef MWVCSOLVERNEW_H_
#define MWVCSOLVERNEW_H_

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/maximum_weighted_matching.hpp>

#include "ConstraintCompositeGraph.h"
#include "MWVCSolver.h"
#include "FastWVCSolver.h"

template <class CCG = ConstraintCompositeGraph<>>
class MWVCSolverFastWVC : public MWVCSolver<CCG>
{
private:
    std::unique_ptr<FastWVCSolver<>> fs;

    double best_weight;
    double best_comp_time;

public:

    MWVCSolverFastWVC(FastWVCSolver<>* fs, VCUpdateHistory* updateHistory)
    {
        this->fs.reset(fs);
        this->updateHistory = updateHistory;
        fs->set_update_history(updateHistory);
    }

    bool is_in_best_vc(int idx)
    {
        return this->updateHistory->bestAssignmentAt(idx);
    }

    virtual double solve(const typename CCG::graph_t &g,
                         map<typename CCG::variable_id_t, bool> &out)
    {
        fs->reset(g, this->vertex_to_id);        
        fs->set_seed(time(NULL));        

        fs->constructVC();
        fs->localSearch();

        if (fs->checkSolution() == 1)
        {
            best_weight = fs->get_best_weight();
        
            cout << "MWVC best weight=" << best_weight << endl;

            auto vertex_id_map = get(boost::vertex_name_t::vertex_name, g);
            auto all_v = vertices(g);
            for (auto it = all_v.first; it != all_v.second; ++it)
            {
                auto id = vertex_id_map[*it];
                if (id >= 0)
                    out[id] = fs->is_in_best_vc(this->vertex_to_id[*it]);
            }
        }
        else
        {
            cout << ", the solution is wrong." << endl;
        }
        
        fs->freeMemory();

        return best_weight;
    }            
};

#endif // MWVCSOLVERNEW_H_
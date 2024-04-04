/** \file MWVCSolverFastWVC.h
 *
 * Define the MWVC solver that uses the Minimum Weighted Vertex Cover Problem solver FastWVC.
 */

#ifndef MWVCSOLVERFASTWVC_H_
#define MWVCSOLVERFASTWVC_H_

#include <time.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <string>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstring>
#include <float.h>
#include <boost/graph/adjacency_list.hpp>

#include "ConstraintCompositeGraph.h"
#include "MWVCSolver.h"
#include "fastWVC.h"

template <class CCG = ConstraintCompositeGraph<>>
class MWVCSolverFastWVC : public MWVCSolver<CCG>
{
private:
    /** Alias of \p CCG::graph_t. */
    typedef typename CCG::graph_t graph_t;
    /** Alias of \p CCG::vertex_t. */
    typedef typename CCG::vertex_t vertex_t;
    typename std::map<vertex_t, int> idx_map = {};

public:

    void init(const typename CCG::graph_t &g)
    {   
        auto u_v_num = num_vertices(g);
        auto u_e_num = num_edges(g);

        if (max(u_v_num, u_e_num) > INT_MAX)
            cout << "The result will likely be wrong. The number of vertices/edges is too high.";

        v_num = (int)u_v_num;
        e_num = (int)u_e_num;

        edge = new Edge[e_num];
        edge_weight = new double[e_num];
        uncov_stack = new int[e_num];
        index_in_uncov_stack = new int[e_num];
        dscore = new double[v_num + 1];
        time_stamp = new llong[v_num + 1];
        v_edges = new int *[v_num + 1];
        v_adj = new int *[v_num + 1];
        v_degree = new int[v_num + 1];
        v_weight = new double[v_num + 1];
        v_in_c = new int[v_num + 1];
        remove_cand = new int[v_num + 1];
        index_in_remove_cand = new int[v_num + 1];
        best_v_in_c = new int[v_num + 1];
        conf_change = new int[v_num + 1];
        tabu_list = new int[v_num + 1];

        fill_n(v_degree, v_num + 1, 0);
        fill_n(tabu_list, v_num + 1, 0);
        fill_n(v_in_c, v_num + 1, 0);
        fill_n(dscore, v_num + 1, 0.0);
        fill_n(conf_change, v_num + 1, 1);
        fill_n(time_stamp, v_num + 1, 0);
        fill_n(edge_weight, e_num, 1.0);

        auto vertex_weight_map = boost::get(boost::vertex_weight, g);

        // vertices
        auto all_v = vertices(g);
        int v = 1;
        for (auto it = all_v.first; it != all_v.second; ++it)
        {
            v_weight[v] = (double) vertex_weight_map[*it];
            idx_map[*it] = v;
            v++;
        }

        // edges
        auto all_e = edges(g);
        int e_idx = 0;
        int v1, v2;
        for (auto it = all_e.first; it != all_e.second; ++it)
        {
            v1 = idx_map[source(*it, g)];
            v2 = idx_map[target(*it, g)];

            v_degree[v1]++;
            v_degree[v2]++;

            edge[e_idx].v1 = v1;
            edge[e_idx].v2 = v2;
            e_idx++;
        }

        v_adj[0] = 0;
        v_edges[0] = 0;
        for (int v = 1; v < v_num + 1; v++)
        {
            v_adj[v] = new int[v_degree[v]];
            v_edges[v] = new int[v_degree[v]];
        }
        
        int *v_degree_tmp = new int[v_num + 1];
        fill_n(v_degree_tmp, v_num + 1, 0);

        for (int e = 0; e < e_num; e++)
        {
            v1 = edge[e].v1;
            v2 = edge[e].v2;

            v_edges[v1][v_degree_tmp[v1]] = e;
            v_edges[v2][v_degree_tmp[v2]] = e;

            v_adj[v1][v_degree_tmp[v1]] = v2;
            v_adj[v2][v_degree_tmp[v2]] = v1;

            v_degree_tmp[v1]++;
            v_degree_tmp[v2]++;
        }
        delete[] v_degree_tmp;
    }

    virtual double solve(const typename CCG::graph_t &g,
                         map<typename CCG::variable_id_t, bool> &out)
    {
        this->init(g);

        initial_time = RunningTime::GetInstance().getCurrentRunningTime();
        // update_listener.init(g, idx_map);

        // Init seed using the current time
        seed = time(NULL);
        srand(seed);

        fastWVC_mode = 0;

        start = chrono::steady_clock::now();

        ConstructVC();
        LocalSearch();

        if (CheckSolution() == 1)
        {
            cout << "MWVC best weight=" << best_weight << ", time=" << best_comp_time << endl;

            auto vertex_id_map = get(boost::vertex_name_t::vertex_name, g);
            auto all_v = vertices(g);
            for (auto it = all_v.first; it != all_v.second; ++it)
            {
                int idx = idx_map[*it];
                auto id = vertex_id_map[*it];
                if (id >= 0)
                    out[id] = best_v_in_c[idx] == 1;
            }
        }
        else
        {
            cout << ", the solution is wrong." << endl;
        }

        FreeMemory();

        return best_weight;
    }            
};

#endif // MWVCSOLVERFASTWVC_H_
#ifndef FASTWVCSOLVER_H_
#define FASTWVCSOLVER_H_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <random>
#include <cmath>
#include <cstring>
#include <boost/graph/adjacency_list.hpp>

#include "ConstraintCompositeGraph.h"
#include "RunningTime.h"

using namespace std;

#define spop(stack) stack[--stack ## _fill_pointer]
#define spush(item, stack) stack[stack ## _fill_pointer++] = item

#define TOLERANCE 1e-6

typedef long long llong;
typedef unsigned int uint;

struct Edge
{
    int v1;
    int v2;
};

template <class CCG = ConstraintCompositeGraph<>>
class FastWVCSolver
{
    private:

        typedef typename CCG::vertex_t vertex_t;
        typedef typename CCG::weight_t weight_t;

        llong max_steps;
        llong step;
        int try_step;
        uint seed;
        int cutoff_time;
        int mode;
    
        int v_num;
        int e_num;
  
        Edge *edge;
        weight_t *edge_weight;

        weight_t *dscore;
        llong *time_stamp;

        weight_t *v_weight;
        int **v_edges;
        int **v_adj;
        int *v_degree;

        int c_size;
        int *v_in_c;
        int *remove_cand;
        int *index_in_remove_cand;
        int remove_cand_size;
        weight_t now_weight;

        int best_c_size;
        int *best_v_in_c;
        double best_comp_time;
        llong best_step;
        weight_t best_weight;

        int *uncov_stack;
        int uncov_stack_fill_pointer;
        int *index_in_uncov_stack;

        int *conf_change;
        int *tabu_list;

        weight_t ave_weight;
        weight_t delta_total_weight;
        weight_t threshold;
        double p_scale;

        VCUpdateHistory *updateHistory;

    public:

        weight_t get_best_weight() { return best_weight; }
    
        bool is_in_best_vc(int idx) { return best_v_in_c[idx] == 1; }

        void set_seed(uint seed) 
        { 
            this->seed = seed; 
            srand(seed);
        }

        void set_update_history(VCUpdateHistory *updateHistory)
        {
            this->updateHistory = updateHistory;
        }

        void reset(const typename CCG::graph_t &g, typename std::map<vertex_t, int> &vertex_to_id)
        {
            auto u_v_num = num_vertices(g);
            auto u_e_num = num_edges(g);

            if (max(u_v_num, u_e_num) > INT_MAX)
                cout << "The result will likely be wrong. The number of vertices/edges is too high.";

            v_num = (int)u_v_num;
            e_num = (int)u_e_num;
            
            edge = new Edge[e_num];
            edge_weight = new weight_t[e_num];
            uncov_stack = new int[e_num];
            index_in_uncov_stack = new int[e_num];
            dscore = new weight_t[v_num + 1];
            time_stamp = new llong[v_num + 1];
            v_edges = new int *[v_num + 1];
            v_adj = new int *[v_num + 1];
            v_degree = new int[v_num + 1];
            v_weight = new weight_t[v_num + 1];
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
                v_weight[v] = (weight_t) vertex_weight_map[*it];
                vertex_to_id.insert(std::make_pair(*it, v));
                v++;
            }

            // edges
            auto all_e = edges(g);
            int e_idx = 0;
            int s, t;
            for (auto it = all_e.first; it != all_e.second; ++it)
            {
                s = vertex_to_id[source(*it, g)];
                t = vertex_to_id[target(*it, g)];

                v_degree[s]++;
                v_degree[t]++;

                edge[e_idx].v1 = s;
                edge[e_idx].v2 = t;
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
                s = edge[e].v1;
                t = edge[e].v2;

                v_edges[s][v_degree_tmp[s]] = e;
                v_edges[t][v_degree_tmp[t]] = e;

                v_adj[s][v_degree_tmp[s]] = t;
                v_adj[t][v_degree_tmp[t]] = s;

                v_degree_tmp[s]++;
                v_degree_tmp[t]++;
            }
            delete[] v_degree_tmp;
            
            mode = 0;
        }

        void constructVC()
        {
            int e;
            int s, t;
            double sdd, tdd;

            uncov_stack_fill_pointer = 0;
            c_size = 0;
            best_weight = (weight_t)(~0U >> 1);
            now_weight = 0;

            for (e = 0; e < e_num; e++)
            {
                s = edge[e].v1;
                t = edge[e].v2;

                if (v_in_c[s] == 0 && v_in_c[t] == 0)
                {
                    sdd = (double)v_degree[s] / (double)v_weight[s];
                    tdd = (double)v_degree[t] / (double)v_weight[t];
                    if (compare(sdd, tdd) > 0)
                    {
                        v_in_c[s] = 1;
                        now_weight += v_weight[s];
                    }
                    else
                    {
                        v_in_c[t] = 1;
                        now_weight += v_weight[t];
                    }
                    c_size++;
                }
            }

            int *save_v_in_c = new int[v_num + 1];
            memcpy(save_v_in_c, v_in_c, sizeof(int) * (v_num + 1));
            int save_c_size = c_size;
            weight_t save_weight = now_weight;

            int times = 50;
            vector<int> blocks(e_num / 1024 + 1);
            for (int i = 0; i < e_num / 1024 + 1; i++)
            {
                blocks[i] = i;
            }

            while (times-- > 0)
            {
                fill_n(v_in_c, v_num + 1, 0);
                c_size = 0;
                now_weight = 0;
                shuffle(blocks.begin(), blocks.end(), default_random_engine(seed));

                for (auto &block : blocks)
                {
                    auto begin = block * 1024;
                    auto end = block == e_num / 1024 ? e_num : begin + 1024;
                    int tmpsize = end - begin + 1;
                    vector<int> idx(tmpsize);
                    for (int i = begin; i < end; i++)
                    {
                        idx[i - begin] = i;
                    }
                    while (tmpsize > 0)
                    {
                        int i = rand() % tmpsize;
                        Edge e = edge[idx[i]];
                        s = e.v1;
                        t = e.v2;
                        swap(idx[i], idx[--tmpsize]);
                        if (v_in_c[s] == 0 && v_in_c[t] == 0)
                        {
                            sdd = (double)v_degree[s] / (double)v_weight[s];
                            tdd = (double)v_degree[t] / (double)v_weight[t];
                            if (compare(sdd, tdd) > 0)
                            {
                                v_in_c[s] = 1;
                                now_weight += v_weight[s];
                            }
                            else
                            {
                                v_in_c[t] = 1;
                                now_weight += v_weight[t];
                            }
                            c_size++;
                        }
                    }
                }
                if (compare(now_weight, save_weight) < 0)
                {
                    save_weight = now_weight;
                    save_c_size = c_size;
                    memcpy(save_v_in_c, v_in_c, sizeof(int) * (v_num + 1));
                }
            }

            now_weight = save_weight;
            c_size = save_c_size;
            memcpy(v_in_c, save_v_in_c, sizeof(int) * (v_num + 1));
            delete[] save_v_in_c;

            for (e = 0; e < e_num; e++)
            {
                s = edge[e].v1;
                t = edge[e].v2;

                if (v_in_c[s] == 1 && v_in_c[t] == 0)
                {
                    dscore[s] -= edge_weight[e];
                }
                else if (v_in_c[t] == 1 && v_in_c[s] == 0)
                {
                    dscore[t] -= edge_weight[e];
                }
            }

            resetRemoveCand();
            for (int v = 1; v < v_num + 1; v++)
            {
                if (v_in_c[v] == 1 && compare(dscore[v], (weight_t)0) == 0)
                {
                    remove(v);
                }
            }
            updateBestSolution();
        }
        
        void localSearch()
        {
            int add_v, remove_v, update_v = 0;
            step = 1;
            try_step = 100;

            ave_weight = 1;
            delta_total_weight = 0;
            p_scale = 0.3;
            threshold = (weight_t)std::floor(0.5 * v_num);

            while (true)
            {
                updateBestSolution();
                update_v = updateTargetSize();

                if (step % try_step == 0)
                {
                    if (RunningTime::GetInstance().isTimeOut())
                        return;
                }

                if (remove_cand_size == 0)
                    return;

                remove_v = chooseRemoveV();
                remove(remove_v);
                time_stamp[remove_v] = step;

                fill_n(tabu_list, v_num + 1, 0);

                while (uncov_stack_fill_pointer > 0)
                {
                    add_v = chooseAddV(remove_v, update_v);
                    add(add_v);
                    updateEdgeWeight();
                    tabu_list[add_v] = 1;
                    time_stamp[add_v] = step;
                }
                removeRedundant();
                step++;
                update_v = 0;
            }
        }

        int checkSolution()
        {
            int e;

            for (e = 0; e < e_num; ++e)
            {
                if (best_v_in_c[edge[e].v1] != 1 && best_v_in_c[edge[e].v2] != 1)
                {
                    cout << ", uncovered edge " << e;
                    return 0;
                }
            }
            return 1;
        }        

        void freeMemory()
        {
            int v;
            for (v = 0; v < v_num + 1; v++)
            {
                delete[] v_adj[v];
                delete[] v_edges[v];
            }

            delete[] conf_change;
            delete[] best_v_in_c;
            delete[] index_in_remove_cand;
            delete[] remove_cand;
            delete[] v_in_c;
            delete[] v_weight;
            delete[] v_degree;
            delete[] v_adj;
            delete[] v_edges;
            delete[] time_stamp;
            delete[] dscore;
            delete[] index_in_uncov_stack;
            delete[] uncov_stack;
            delete[] edge_weight;
            delete[] edge;
        }
    
    private:

        void resetRemoveCand()
        {
            int v;
            int j = 0;

            for (v = 1; v < v_num + 1; v++)
            {
                if (v_in_c[v] == 1)
                {
                    remove_cand[j] = v;
                    index_in_remove_cand[v] = j;
                    j++;
                }
                else
                {
                    index_in_remove_cand[v] = 0;
                }
            }

            remove_cand_size = j;
        }

        void uncover(int e)
        {
            index_in_uncov_stack[e] = uncov_stack_fill_pointer;
            spush(e, uncov_stack);
        }

        void cover(int e)
        {
            int index, last_uncov_edge;
            last_uncov_edge = spop(uncov_stack);
            index = index_in_uncov_stack[e];
            uncov_stack[index] = last_uncov_edge;
            index_in_uncov_stack[last_uncov_edge] = index;
        }

        void add(int v)
        {
            int i, e, n;
            int edge_count = v_degree[v];

            v_in_c[v] = 1;
            c_size++;
            dscore[v] = -dscore[v];
            now_weight += v_weight[v];

            remove_cand[remove_cand_size] = v;
            index_in_remove_cand[v] = remove_cand_size++;

            for (i = 0; i < edge_count; i++)
            {
                e = v_edges[v][i];
                n = v_adj[v][i];

                if (v_in_c[n] == 0)
                {
                    dscore[n] -= edge_weight[e];
                    conf_change[n] = 1;
                    cover(e);
                }
                else
                {
                    dscore[n] += edge_weight[e];
                }
            }
        }

        void remove(int v)
        {
            int i, e, n;
            int edge_count = v_degree[v];

            v_in_c[v] = 0;
            c_size--;
            dscore[v] = -dscore[v];
            conf_change[v] = 0;

            int last_remove_cand_v = remove_cand[--remove_cand_size];
            int index = index_in_remove_cand[v];
            remove_cand[index] = last_remove_cand_v;
            index_in_remove_cand[last_remove_cand_v] = index;
            index_in_remove_cand[v] = 0;

            now_weight -= v_weight[v];

            for (i = 0; i < edge_count; i++)
            {
                e = v_edges[v][i];
                n = v_adj[v][i];

                if (v_in_c[n] == 0)
                {
                    dscore[n] += edge_weight[e];
                    conf_change[n] = 1;
                    uncover(e);
                }
                else
                {
                    dscore[n] -= edge_weight[e];
                }
            }
        }

        int updateTargetSize()
        {
            int v;
            int best_remove_v;
            double best_dscore;
            double dscore_v;

            best_remove_v = remove_cand[0];
            best_dscore = (double)v_weight[best_remove_v] / (double)abs(dscore[best_remove_v]);

            if (compare(dscore[best_remove_v], (weight_t)0) != 0)
            {
                for (int i = 1; i < remove_cand_size; i++)
                {
                    v = remove_cand[i];
                    if (compare(dscore[v], (weight_t)0) == 0) break;
                    dscore_v = (double)v_weight[v] / (double)abs(dscore[v]);
                    if (compare(dscore_v, best_dscore) > 0)
                    {
                        best_dscore = dscore_v;
                        best_remove_v = v;
                    }
                }
            }

            remove(best_remove_v);

            return best_remove_v;
        }

        int chooseRemoveV()
        {
            int i, v;
            double dscore_v, dscore_remove_v;
            int remove_v = remove_cand[rand() % remove_cand_size];
            int to_try = 50;

            for (i = 1; i < to_try; i++)
            {
                v = remove_cand[rand() % remove_cand_size];
                dscore_v = (double)v_weight[v] / (double)abs(dscore[v]);
                dscore_remove_v = (double)v_weight[remove_v] / (double)abs(dscore[remove_v]);

                if (tabu_list[v] == 1)
                {
                    continue;
                }
                if (compare(dscore_v, dscore_remove_v) < 0)
                {
                    continue;
                }
                if (compare(dscore_v, dscore_remove_v) > 0)
                {
                    remove_v = v;
                }
                else if (time_stamp[v] < time_stamp[remove_v])
                {
                    remove_v = v;
                }
            }
            return remove_v;
        }

        int chooseAddFromV()
        {
            int v;
            int add_v = 0;
            double improvement = 0.0;
            double dscore_v;

            for (v = 1; v < v_num + 1; v++)
            {
                if (v_in_c[v] == 1)
                {
                    continue;
                }
                if (conf_change[v] == 0)
                {
                    continue;
                }
                dscore_v = (double)dscore[v] / (double)v_weight[v];
                if (compare(dscore_v, improvement) > 0)
                {
                    improvement = dscore_v;
                    add_v = v;
                }
                else if (compare(dscore_v, improvement) == 0)
                {
                    if (time_stamp[v] < time_stamp[add_v])
                    {
                        add_v = v;
                    }
                }
            }
            return add_v;
        }

        int chooseAddV(int remove_v, int update_v = 0)
        {
            int i, v;
            int add_v = 0;
            double improvement = 0.0;
            double dscore_v;

            int tmp_degree = v_degree[remove_v];

            int degree_sum;
            int *adjp = v_adj[remove_v];
            for (i = 0; i < tmp_degree; i++)
            {
                v = adjp[i];
                if (v_in_c[v] == 1)
                {
                    continue;
                }
                if (conf_change[v] == 0)
                {
                    continue;
                }
                dscore_v = (double)dscore[v] / (double)v_weight[v];
                if (compare(dscore_v, improvement) > 0)
                {
                    improvement = dscore_v;
                    add_v = v;
                }
                else if (compare(dscore_v, improvement) == 0)
                {
                    if (time_stamp[v] < time_stamp[add_v])
                    {
                        add_v = v;
                    }
                }
            }
            v = remove_v;
            if (conf_change[v] == 1 && v_in_c[v] == 0)
            {
                dscore_v = (double)dscore[v] / (double)v_weight[v];
                if (compare(dscore_v, improvement) > 0)
                {
                    improvement = dscore_v;
                    add_v = v;
                }
                else if (compare(dscore_v, improvement) == 0)
                {
                    if (time_stamp[v] < time_stamp[add_v])
                    {
                        add_v = v;
                    }
                }

            }

            if (update_v != 0)
            {
                tmp_degree = v_degree[update_v];
                adjp = v_adj[update_v];
                for (i = 0; i < tmp_degree; i++)
                {
                    v = adjp[i];
                    if (v_in_c[v] == 1)
                    {
                        continue;
                    }
                    if (conf_change[v] == 0)
                    {
                        continue;
                    }
                    dscore_v = (double)dscore[v] / (double)v_weight[v];
                    if (compare(dscore_v, improvement) > 0)
                    {
                        improvement = dscore_v;
                        add_v = v;
                    }
                    else if (compare(dscore_v, improvement) == 0)
                    {
                        if (time_stamp[v] < time_stamp[add_v])
                        {
                            add_v = v;
                        }
                    }
                }
                v = update_v;
                if (conf_change[v] == 1 && v_in_c[v] == 0)
                {
                    dscore_v = (double)dscore[v] / (double)v_weight[v];
                    if (compare(dscore_v, improvement) > 0)
                    {
                        improvement = dscore_v;
                        add_v = v;
                    }
                    else if (compare(dscore_v, improvement) == 0)
                    {
                        if (time_stamp[v] < time_stamp[add_v])
                        {
                            add_v = v;
                        }
                    }
                }

            }

            return add_v;
        }

        void updateBestSolution()
        {
            int v;

            if (compare(now_weight, best_weight) < 0)
            {
                for (v = 1; v < v_num + 1; v++)
                {
                    best_v_in_c[v] = v_in_c[v];
                }
                best_weight = now_weight;
                best_c_size = c_size;
                best_comp_time = RunningTime::GetInstance().getCurrentRunningTime();
                best_step = step; 
                updateHistory->emplace_back_vc(
                    best_v_in_c, 
                    best_v_in_c + v_num * sizeof(best_v_in_c[0]), 
                    best_comp_time
                );

            }    
        }

        void removeRedundant()
        {
            int v;
            for (int i = 0; i < remove_cand_size; i++)
            {
                v = remove_cand[i];
                if (v_in_c[v] == 1 && compare(dscore[v], (weight_t)0) == 0)
                {
                    remove(v);
                    i--;
                }
            }
        }

        void forgetEdgeWeights()
        {
            int v, e;
            int new_total_weitght = 0;

            for (v = 1; v < v_num + 1; v++)
            {
                dscore[v] = 0;
            }

            for (e = 0; e < e_num; e++)
            {
                edge_weight[e] = edge_weight[e] * p_scale;
                new_total_weitght += edge_weight[e];

                if (v_in_c[edge[e].v1] + v_in_c[edge[e].v2] == 0)
                {
                    dscore[edge[e].v1] += edge_weight[e];
                    dscore[edge[e].v2] += edge_weight[e];
                }
                else if (v_in_c[edge[e].v1] + v_in_c[edge[e].v2] == 1)
                {
                    if (v_in_c[edge[e].v1] == 1)
                    {
                        dscore[edge[e].v1] -= edge_weight[e];
                    }
                    else
                    {
                        dscore[edge[e].v2] -= edge_weight[e];
                    }
                }
            }
            ave_weight = new_total_weitght / e_num;
        }

        void updateEdgeWeight()
        {
            int i, e;

            for (i = 0; i < uncov_stack_fill_pointer; i++)
            {
                e = uncov_stack[i];
                edge_weight[e] += 1;
                dscore[edge[e].v1] += 1;
                dscore[edge[e].v2] += 1;
                if (mode % 2 == 1)
                {
                    conf_change[edge[e].v1] = 1;
                    conf_change[edge[e].v2] = 1;
                }
            }

            delta_total_weight += uncov_stack_fill_pointer;

            if (mode / 2 == 1)
            {
                if (delta_total_weight >= e_num)
                {
                    ave_weight += 1;
                    delta_total_weight -= e_num;
                }

                if (compare(ave_weight, threshold) >= 0)
                {
                    forgetEdgeWeights();
                }
            }
        }

        /** \brief Compare two numbers a and b. 
            \return +1 if a > b,
                     0 if a = b,
                    -1 if a < b 
        */
        int compare(int a, int b)
        {
            if (a == b) return 0;
            return a > b ? 1 : -1;
        }

        int compare(long a, long b)
        {
            if (a == b) return 0;
            return a > b ? 1 : -1;
        }

        int compare(llong a, llong b)
        {
            if (a == b) return 0;
            return a > b ? 1 : -1;
        }

        int compare(float a, float b)
        {   
            if (abs(a - b) <= TOLERANCE) return 0;
            return a > b ? 1 : -1;
        }

        int compare(double a, double b)
        {
            if (abs(a - b) <= TOLERANCE) return 0;
            return a > b ? 1 : -1;
        }
};

#endif // FASTWVCSOLVER_H_
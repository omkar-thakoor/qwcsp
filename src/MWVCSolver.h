/*
  Copyright (c) 2016-2017 Hong Xu

  This file is part of WCSPLift.

  WCSPLift is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  WCSPLift is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with WCSPLift.  If not, see <http://www.gnu.org/licenses/>.
*/

/** \file MWVCSolver.h
 *
 * Defines the base class for all solvers that solve the MWVC problem on a CCG.
 */
#ifndef MWVCSOLVER_H_
#define MWVCSOLVER_H_

#include "ConstraintCompositeGraph.h"
#include "UpdateHistory.h"

/** The base class for all solvers that solve the MWVC problem on a CCG.
 *
 * \tparam The type of the CCG. It defaults to \p ConstraintCompositeGraph<>.
 */
template <class CCG = ConstraintCompositeGraph<> >
class MWVCSolver
{
protected:
    typedef typename CCG::graph_t graph_t;
    typedef typename boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
    typedef typename boost::graph_traits<graph_t>::edge_descriptor edge_t;

    std::map<vertex_t, int> vertex_to_id;

    AssignmentsUpdateHistory* updateHistory;

public:
    /** \brief Solve the MWVC problem on a given CCG.
     *
     * \param[in] g The graph to solve the MWVC problem on.
     *
     * \param[out] out A map indicating whether each variable is in the MWVC. Note that if \p out is
     * not empty when it is passed in, it will be merged with the existing key-value pairs.
     *
     * \return The total weight of the MWVC.
     */
    virtual double solve(const typename CCG::graph_t& g,
                         typename std::map<typename CCG::variable_id_t, bool>& out) = 0;

    bool has_next_update() { return updateHistory->has_next_update(); }

    double next_update(const typename CCG::graph_t &g,
                       std::map<typename CCG::variable_id_t, bool> &out)
    {   
        auto ttime = updateHistory->next_update();

        auto vertex_id_map = get(boost::vertex_name_t::vertex_name, g);
        auto all_v = vertices(g);
        for (auto it = all_v.first; it != all_v.second; ++it)
        {
            auto id = vertex_id_map[*it];
            if (id >= 0)
                out[id] = updateHistory->assignmentAt(vertex_to_id[*it]);
        }
        return ttime;
    }
};

#endif // MWVCSOLVER_H_

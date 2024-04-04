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

/** \file WCSPInstance.h
 *
 * Define the definition of WCSP instances and the constraints in it. Algorithms that can be applied
 * directly on WCSP instances, i.e., the min-sum message passing algorithm \cite xkk17 and integer
 * linear programming \cite xkk17a, are also included in this file.
 */

#ifndef WCSPINSTANCE_H_
#define WCSPINSTANCE_H_

#include <array>
#include <cstdint>
#include <istream>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include <gecode/minimodel.hh>
#include "../gecode-release-6.2.0/contribs/qecode/QProblem.h"
#include "../gecode-release-6.2.0/contribs/qecode/qsolver_qcop.h"
#include <boost/dynamic_bitset.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <cblas.h>

#include "global.h"
#include "Recorder.h"
#include "RunningTime.h"
//#include "LinearProgramSolver.h"

using namespace Gecode;
using namespace Gecode::Int;
using namespace qecode;
using namespace std;

/** \brief Weighted constraint.
 *
 * This class describes a weighted constraint. It consists of the IDs of the variables as well as
 * how the weight corresponding to each assignment.
 *
 * \tparam VarIdType The type of variable IDs. It must be an integer type and defaults to \p
 * intmax_t.
 *
 * \tparam WeightType The type of weights. It must be a numeric type and defaults to \p double.
 *
 * \tparam ValueType The type of values of all variables in the constraint. It should usually be a
 * bitset type and defaults to \p boost::dynamic_bitset<>.
 */
template < class VarIdType = intmax_t,
            class WeightType = double,
            class ValueType = boost::dynamic_bitset<> > // The type of the values to represent each constraint>
class Constraint
{
public:
    /** Alias of \p VarIDType. */
    typedef VarIdType variable_id_t;

    /** Alias of \p ValueType. */
    typedef ValueType value_t;

    /** Alias of \p WeightType. */
    typedef WeightType weight_t;

private:
    // This class is only used for Polynomial key comparison. This function ensures keys with
    // a smaller number of variables must precede larger number of variables.
    class PolynomialKeyComparison
    {
    public:
        bool operator () (const std::set<variable_id_t>& a, const std::set<variable_id_t>& b) const
        {
            if (a.size() > b.size())
                return true;

            if (a.size() < b.size())
                return false;

            return a < b;
        }
    };
public:
    /** The polynomial form of constraints consisting of the coefficient of each term (each term is
     * represented by an assignment of values to variables). */
    typedef std::map<std::set<variable_id_t>, weight_t, PolynomialKeyComparison> Polynomial;

private:
    std::vector<variable_id_t> variables;
    std::map<value_t, weight_t> weights;
    weight_t maxWeight = 0;


public:
    /** Alias of \p Constraint<VarIdType, WeightType, ValueType>, i.e., the class type itself. */
    typedef Constraint<VarIdType, WeightType, ValueType> self_t;

    /** \brief Get the list of the IDs of variables in the constraint.
     *
     * \return A list of variable IDs.
     */
    inline const std::vector<variable_id_t>& getVariables() const noexcept
    {
        return variables;
    }

    /** \brief Get the ID of the <tt>i</tt>'th variable.
     *
     * \param[in] i The index of the variable of which to get the ID.
     *
     * \return The ID of the <tt>i</tt>'th variable.
     */
    inline variable_id_t getVariable(size_t i) const noexcept
    {
        return variables.at(i);
    }

    /** \brief Set the IDs of the variables in the constraint from a \p std::vector<variable_id_t>
     * object.
     *
     * \param[in] variables: A list of variable IDs to set.
     */
    inline void setVariables(const std::vector<variable_id_t>& variables) noexcept
    {
        this->variables = variables;
        this->weights.clear();
    }

    /** \brief Set the IDs of the variables in the constraint from a range of iterators.
     *
     * \param[in] b the beginning iterator.
     *
     * \param[in] e the ending iterator.
     */
    template<class Iter>
    inline void setVariables(Iter b, Iter e) noexcept
    {
        variables.clear();
        std::copy(b, e, std::back_inserter(variables));
    }

    /** \brief Set the weight of a given assignment of values to variables specified by \p v to \p
     * w.
     *
     * \param[in] v An assignment of values to variables.
     *
     * \param[in] w A weight.
     *
     * \return A reference to the class object itself.
     */
    inline self_t& setWeight(const value_t& v, weight_t w)
    {
        weights[v] = w;
        if (w > maxWeight) maxWeight = w;
        return *this;
    }

    /** \brief Get the weight of a given assignment of values to variables specified by \p v.
     *
     * \param[in] v The assignment of values to variables.
     *
     * \return The weight of a given assignment of values to variables specified by \p v.
     */
    inline weight_t getWeight(const value_t& v) const noexcept
    {
        auto it = weights.find(v);
        if (it == weights.end()) // not found, return 0
            return weight_t();

        return it->second;
    }

    /** \brief Get a map object that maps assignments of values to variables to weights.
     *
     * \return A map object that maps assignments of values to variables to weights.
     */
    inline const std::map<value_t, weight_t> getWeights() const noexcept
    {
        return weights;
    }

    inline weight_t getMaxWeight() const noexcept
    {
        return maxWeight;
    }

    /** \brief Represent the constraint using a \p boost::property_tree::ptree object.
     *
     * \param[out] t A \p boost::property_tree::ptree object that represents the constraint.
     */
    void toPropertyTree(boost::property_tree::ptree& t) const
    {
        using namespace boost::property_tree;

        t.clear();

        // put the variables
        ptree vs;
        for (auto& v : variables)
            vs.push_back(
                std::make_pair("", ptree(std::to_string(v))));

        t.put_child("variables", vs);

        // put the weights
        ptree ws;
        for (auto& w : weights)
        {
            ptree we;

            for (size_t i = 0; i < getVariables().size(); ++ i)
            {
                if (w.first[i])
                    we.push_back(std::make_pair("", ptree("1")));
                else
                    we.push_back(std::make_pair("", ptree("0")));
            }

            we.push_back(std::make_pair("", ptree(std::to_string(w.second))));

            ws.push_back(std::make_pair("", we));
        }

        t.put_child("weights", ws);
    }

    /** \brief Compute the coefficients of the polynomial converted from a constraint according to
     * \cite k08.
     *
     * \param[out] p A Constraint::Polynomial object corresponding to the polynomial representation
     * of the constraint.
     */
    void toPolynomial(Polynomial& p) const noexcept
    {
        size_t s = getVariables().size();

        // we basically solve A * x = b

        size_t max_int = ((size_t) 1) << s;

        // matrix a
        double * a = new double[max_int * max_int];
        memset(a, 0, max_int * max_int * sizeof(double));

/// \cond
#define ENTRY(a, x, y) a[(y) * max_int + (x)]
/// \endcond
        // set the first column to 1
        for (size_t i = 0; i < max_int; ++ i)
            ENTRY(a, i, 0) = 1.0;
        // Set the lower triangle. The row number corresponds to the assignments of variables, and
        // the col number corresponds to terms of the polynomial in the following order:
        // X_1, X_2, X_1 X_2, X_3, X_1 X_3, X_2 X_3, X_1 X_2 X_3...
        for (size_t i = 1; i < max_int; ++ i)
            for (size_t j = 1; j <= i; ++ j)
                // (i|j)==i : the 1-bits integer i includes all the 1-bits of j
                ENTRY(a, i, j) = ((i | j) == i) ? 1.0 : 0.0;
#undef ENTRY

        // assign the vector b
        double * x = new double[max_int];
        for (size_t i = 0; i < max_int; ++ i)
            x[i] = getWeight(value_t(s, i));

        cblas_dtrsv(CblasColMajor, CblasLower, CblasNoTrans, CblasUnit, max_int, a, max_int, x, 1);
        delete[] a;

        // Insert all the coefficients into the polynomial.
        for (size_t i = 0; i < max_int; ++ i)
        {
            typename Polynomial::key_type k;

            for (size_t j = 0; j < s; ++ j)
                if ((((size_t) 1u) << j) & i)
                    k.insert(getVariable(j));

            p[k] += x[i];
        }

        delete[] x;
    }
};

/** \brief A WCSP instance.
 *
 * This class describes a WCSP instance. It consists a set of constraints.
 *
 * \tparam VarIdType The type of variable IDs. It must be an integer type and defaults to \p
 * intmax_t.
 *
 * \tparam WeightType The type of weights. It must be a numeric type and defaults to \p double.
 *
 * \tparam ConstraintValueType The type of values of all variables in the constraints in the WCSP
 * instance. It should usually be a bitset type and defaults to \p boost::dynamic_bitset<>.
 */
template <class VarIdType = intmax_t, class WeightType = double,
          class ConstraintValueType = boost::dynamic_bitset<> >
class WCSPInstance
{
public:
    /** Alias of \p VarIdType. */
    typedef VarIdType variable_id_t;
    /** Alias of \p WeightType. */
    typedef WeightType weight_t;
    /** Alias of \p ConstraintValueType. */
    typedef ConstraintValueType constraint_value_t;

public:
    /** Alias of the Constraint type in the WCSP instance, i.e., \p
     * Constraint<variable_id_t,weight_t,constraint_value_t>.
     */
    typedef Constraint<variable_id_t, weight_t, constraint_value_t> constraint_t;

private:
    std::vector<constraint_t> constraints;
    size_t numVars;
    weight_t totMaxWeight = 0;
    std::vector<int> optType;

public:

    /** \brief Get the opt types.
     *
     * \return A list of \p constraint_t objects.
     */
    inline const std::vector<int> getOptType() const noexcept
    {
        return optType;
    }

    /** \brief Get the list of constraints.
     *
     * \return A list of \p constraint_t objects.
     */
    inline const std::vector<constraint_t>& getConstraints() const noexcept
    {
        return constraints;
    }

    /** \brief Get the <tt>i</tt>'th constraint.
     *
     * \param[in] i The index of the Constraint object to get.
     *
     * \return The <tt>i</tt>'th constraint.
     */
    inline const constraint_t& getConstraint(size_t i) const noexcept
    {
        return constraints.at(i);
    }

    /** \brief Compute the total weight corresponding to an assignment of values to all variables.
     *
     * \param[in] assignments The assignment of values to all variables.
     *
     * \return The computed total weight.
     */
    inline weight_t computeTotalWeight(
        const std::map<variable_id_t, bool>& assignments) const noexcept
    {
        weight_t tw = 0;

        for (const auto& c : constraints)
        {
            constraint_value_t val(c.getVariables().size());
            for (size_t i = 0; i < c.getVariables().size(); ++ i)
            {
                try
                {
                    val[i] = assignments.at(c.getVariables().at(i));
                } catch (std::out_of_range& e)
                {
                    // Here, it means the variable assignments does not hold c.getVariables().at(i)
                    // -- the CCG does not contain that variable. This may well be that the variable
                    // itself is a "dummy" variable -- the value of it does not affect the total
                    // weight. In this case, we always assign zero (or anything else) to it.
                    val[i] = 0;
                }
            }

            tw += c.getWeight(val);
        }

        return tw;
    }

public:
    /** \brief Load a problem instance from a stream.
     *
     * \param[in] f The input stream.
     */
    void load(std::istream& f);

    /** \brief Load the problem in DIMACS format. The format specification is at
     * http://graphmod.ics.uci.edu/group/WCSP_file_format
     *
     * \param[in] f The input stream.
     */
    void loadDimacs(std::istream& f);

    /** \brief Load the problem in UAI format. The format specification is at
     * http://www.hlt.utdallas.edu/~vgogate/uai14-competition/modelformat.html
     *
     * \param[in] f The input stream.
     */
    void loadUAI(std::istream& f);

    /** List of supported input file formats.
     */
    enum class Format
    {
        DIMACS,
        UAI
    };

    /** \brief Construct a WCSP instance from an input stream in a given format.
     *
     * \param[in] f The input stream.
     *
     * \param[in] format The format of \p f.
     */
    WCSPInstance(std::istream& f, Format format)
    {
        switch (format)
        {
        case Format::DIMACS:
            loadDimacs(f);
            break;
        case Format::UAI:
            loadUAI(f);
            break;
        }
    }

    /** \brief Convert this WCSP instance to a human-readable string.
     *
     * \return The human-readable string.
     */
    std::string toString() const noexcept
    {
        std::stringstream ss;

        ss << '{' << std::endl;
        for (auto& c : constraints)
        {
            ss << c;
        }

        return ss.str();
    }

    /** \brief Solve the WCSP instance using the min-sum message passing algorithm.
     *
     * \param[in] delta The threshold for convergence determination. That is, the algorithm
     * terminates iff all messages change within \p delta compared with the previous iteration.
     *
     * \return A solution to the WCSP instance.
     */
    std::map<variable_id_t, bool> solveUsingMessagePassing(weight_t delta) const noexcept
    {
        std::map<variable_id_t, std::set<const constraint_t*> > constraint_list_for_v;

        // Initialize all messages.
        std::map<std::pair<variable_id_t, const constraint_t*>, std::array<weight_t, 2> > msgs_v_to_c;
        std::map<std::pair<const constraint_t*, variable_id_t>, std::array<weight_t, 2> > msgs_c_to_v;
        for (const auto& c : constraints)
        {
            for (auto v : c.getVariables())
            {
                constraint_list_for_v[v].insert(&c);
                msgs_v_to_c[std::make_pair(v, &c)] = {0, 0};
                msgs_c_to_v[std::make_pair(&c, v)] = {0, 0};
            }
        }

        bool converged = false;
        uintmax_t num_iterations = 0;
        do
        {
            if (RunningTime::GetInstance().isTimeOut())  // time's up
                break;

            ++ num_iterations;

            converged = true;

            auto msgs_v_to_c0 = msgs_v_to_c;
            auto msgs_c_to_v0 = msgs_c_to_v;

            // Update all the messages from variables to constraints.
            for (auto& it : msgs_v_to_c)
            {
                auto v = it.first.first;
                auto c = it.first.second;

                std::array<weight_t, 2> m = {0, 0};

                for (auto nc : constraint_list_for_v.at(v))
                {
                    if (nc != c)
                    {
                        const auto m_to_v = msgs_c_to_v0.find(std::make_pair(nc, v))->second;

                        m[0] += m_to_v.at(0);
                        m[1] += m_to_v.at(1);
                    }
                }

                it.second = m;

                // is it now convergent?
                if (converged)
                {
                    auto old_m = msgs_v_to_c0.at(it.first);

                    if (std::abs(old_m.at(0) - m.at(0)) > delta ||
                        std::abs(old_m.at(1) - m.at(1)) > delta)
                        converged = false;
                }
            }

            // Update all the messages from constraints to variables.
            for (auto& it : msgs_c_to_v)
            {
                auto c = it.first.first;
                auto v = it.first.second;

                std::vector<std::array<weight_t, 2> > m_to_c;
                m_to_c.reserve(c->getVariables().size());
                std::vector<variable_id_t> v_ids;
                v_ids.reserve(c->getVariables().size());
                size_t self_index;   // index of the variable itself
                for (size_t i = 0; i < c->getVariables().size(); ++ i)
                {
                    auto nv = c->getVariable(i);

                    m_to_c.push_back(msgs_v_to_c0.at(std::make_pair(nv, c)));
                    v_ids.push_back(nv);

                    if (nv == v)
                        self_index = i;
                }

                std::array<weight_t, 2> m = {
                    std::numeric_limits<weight_t>::max(), std::numeric_limits<weight_t>::max()
                };

                // TODO: Change the iteration which does not have an upper limit of number of bits.
                for (unsigned long i = 0; i < (1ul << c->getVariables().size()); ++ i)
                {
                    constraint_value_t values(c->getVariables().size(), i);

                    weight_t sum = c->getWeight(values);

                    for (size_t j = 0; j < c->getVariables().size(); ++ j)
                        if (c->getVariable(j) != v)
                            sum += msgs_v_to_c0.at(
                                std::make_pair(c->getVariable(j), c))[values[j] ? 1 : 0];

                    size_t index_update = values[self_index] ? 1 : 0;

                    if (m[index_update] > sum)
                        m[index_update] = sum;
                }

                auto m_min = m[0] < m[1] ? m[0] : m[1];
                m[0] -= m_min;
                m[1] -= m_min;

                it.second = m;

                // is it now convergent?
                if (converged)
                {
                    auto old_m = msgs_c_to_v0.at(it.first);

                    if (std::abs(old_m.at(0) - m.at(0)) > delta ||
                        std::abs(old_m.at(1) - m.at(1)) > delta)
                        converged = false;
                }
            }
        } while (!converged);

        std::cout << "Number of iterations: " << num_iterations << std::endl;

        // compute the variable values
        std::map<variable_id_t, std::array<weight_t, 2> > assignments_w;

        converged = true;
        for (const auto& c : constraints)
            for (auto v : c.getVariables())
            {
                assignments_w[v][0] += msgs_c_to_v[std::make_pair(&c, v)][0];
                assignments_w[v][1] += msgs_c_to_v[std::make_pair(&c, v)][1];
                if (isinf(assignments_w[v][0]) || isinf(assignments_w[v][1]))
                    converged = false;
            }

        std::map<variable_id_t, bool> assignments;
        for (const auto& a : assignments_w)
            assignments[a.first] = a.second[0] > a.second[1] ? a.second[1] : a.second[0];

        Recorder::getInstance().set_mp_converged(converged);
        if (!converged)
            std::cout << "*** Message passing not converged! ***" << std::endl;

        return assignments;
    }

    /** \brief Solve the WCSP instance using linear programming.
     *
     * \param[in] lps The linear program solver to be used.
     *
     * \return A solution to the WCSP instance.
     
    std::map<variable_id_t, bool> solveUsingLinearProgramming(LinearProgramSolver& lps) const noexcept
    {
        lps.reset();
        lps.setTimeLimit(RunningTime::GetInstance().getTimeLimit().count());

        lps.setObjectiveType(LinearProgramSolver::ObjectiveType::MIN);

        std::vector<std::vector<LinearProgramSolver::variable_id_t> > lp_variables(
            this->getConstraints().size());

        // add LP variables and constraints
        for (size_t i = 0; i < this->getConstraints().size(); ++ i)
        {
            size_t num_values = 1 << this->getConstraint(i).getVariables().size();
            lp_variables[i].reserve(num_values);
            for (size_t j = 0; j < num_values; ++ j)  // add all variables
                lp_variables[i].push_back(
                    lps.addVariable(
                        this->getConstraint(i).getWeight(
                            constraint_value_t(this->getConstraint(i).getVariables().size(), j))));

            // Only one value in a WCSP constraint can take effect.
            lps.addConstraint(lp_variables[i], std::vector<double>(num_values, 1.0), 1.0,
                              LinearProgramSolver::ConstraintType::EQUAL);
        }

        // Add constraints that make sure LP variables that represent overlapped WCSP variables are
        // consistent.
        for (size_t i = 0; i < lp_variables.size(); ++ i)
        {
            for (size_t j = 0; j < i; ++ j)
            {
                // Get the overlap map.
                std::vector<std::pair<size_t, size_t> > overlapped_variables;

                // Whether a bit is overlapping?
                std::vector<bool> vs_overlap_p[2];
                vs_overlap_p[0].resize(this->getConstraint(i).getVariables().size());
                vs_overlap_p[1].resize(this->getConstraint(j).getVariables().size());

                for (size_t k0 = 0; k0 < this->getConstraint(i).getVariables().size(); ++ k0)
                    for (size_t k1 = 0; k1 < this->getConstraint(j).getVariables().size(); ++ k1)
                    {
                        if (this->getConstraint(i).getVariable(k0) ==
                            this->getConstraint(j).getVariable(k1))
                        {
                            overlapped_variables.push_back(std::make_pair(k0, k1));
                            vs_overlap_p[0][k0] = true;
                            vs_overlap_p[1][k1] = true;
                            break;
                        }
                    }

                constraint_value_t vs[2];
                vs[0].resize(this->getConstraint(i).getVariables().size());
                vs[1].resize(this->getConstraint(j).getVariables().size());

                // For each overlapped variable assignment, add a new LP constraint.
                constraint_value_t overlapped(overlapped_variables.size(), 0);
                do {
                    for (size_t k = 0; k < overlapped_variables.size(); ++ k)
                    {
                        vs[0][overlapped_variables[k].first] = overlapped[k];
                        vs[1][overlapped_variables[k].second] = overlapped[k];
                    }

                    // bits that represent non-overlapped variables
                    constraint_value_t vs_non_overlapped[2];
                    vs_non_overlapped[0].resize((vs[0].size() - overlapped_variables.size()));
                    vs_non_overlapped[1].resize((vs[1].size() - overlapped_variables.size()));

                    std::vector<LinearProgramSolver::variable_id_t> lp_vs;  // LP variables
                    std::vector<double> coefs;
                    size_t reserved_size = (2 << vs_non_overlapped[0].size()) + 2 << vs_non_overlapped[1].size();
                    lp_vs.reserve(reserved_size);
                    coefs.reserve(reserved_size);

                    for (size_t z = 0; z < 2; ++ z)
                    {
                        do {
                            for (size_t k = 0, l = 0; k < vs_non_overlapped[z].size(); ++ k, ++ l)
                            {
                                for (; vs_overlap_p[z][l]; ++ l)
                                    ;

                                vs[z][l] = vs_non_overlapped[z][k];
                            }

                            if (z == 0)
                            {
                                lp_vs.push_back(lp_variables[i][vs[z].to_ulong()]);
                                coefs.push_back(1.0);
                            }
                            else
                            {
                                lp_vs.push_back(lp_variables[j][vs[z].to_ulong()]);
                                coefs.push_back(-1.0);
                            }

                            bitset_increase_1(vs_non_overlapped[z]);
                        } while (!vs_non_overlapped[z].none());
                    }

                    lps.addConstraint(lp_vs, coefs, 0.0,
                                      LinearProgramSolver::ConstraintType::EQUAL);

                    bitset_increase_1(overlapped);
                } while (!overlapped.none());
            }
        }

        std::vector<double> assignments;
        lps.solve(assignments);

        // Compute the solutions from assignments
        std::map<variable_id_t, bool> solution;
        for (size_t i = 0; i < lp_variables.size(); ++ i)
        {
            for (size_t j = 0; j < lp_variables[i].size(); ++ j)
            {
                if (assignments[lp_variables[i][j]] > 0.99)   // the corresponding constraint holds
                {
                    constraint_value_t value(this->getConstraint(i).getVariables().size(), j);
                    for (size_t k = 0; k < this->getConstraint(i).getVariables().size(); ++ k)
                        solution[this->getConstraint(i).getVariable(k)] = value[k];
                    break;
                }
            }
        }

        return solution;
    }*/

    /** \brief Solve the QWCSP instance using QCOP (QeCode)
     *
     * \param[in] lps The quantifiers
     *
     * \return A solution to the WCSP instance.
     */
    std::map<variable_id_t, bool> solveUsingQCOP() const noexcept
    {
        int numLevels = numVars+1;
        bool* quants = new bool[numLevels];
        int* scopeVars = new int[numLevels];

        std::fill_n(quants, numLevels, false);
        std::fill_n(scopeVars, numVars, 1);
        scopeVars[numLevels-1] = constraints.size() + 1; //Add one for total cost
        QProblem problem(numLevels, quants, scopeVars);

        //cout<<"QP declared"<<endl;
        for (int varInd = 0; varInd < numVars; varInd++) {
            problem.QIntVar(varInd,0,1);
            IntVarArgs branchNext(varInd+1);
            for (int i=0;i<=varInd;i++)
                branchNext[i] = problem.var(i);
            branch(*(problem.space()),branchNext,INT_VAR_SIZE_MIN(),INT_VAL_MIN());
            problem.nextScope();
        }
        
        //cout<<"var scopes done"<<endl;
        int totMaxWeight = 0;
        for (size_t i = 0; i < constraints.size(); i++) {
            int mw = (int) constraints[i].getMaxWeight();
            problem.QIntVar(numVars + i,0,mw);
            totMaxWeight += mw;
            size_t num_constr_vars = constraints[i].getVariables().size();
            size_t num_values = 1 << num_constr_vars;
            
            //cout<<"constr "<<i<<" var added"<<endl;
            for (size_t j = 0; j < num_values; ++ j) {  // add all variables
                constraint_value_t values(num_constr_vars, j);
                int weight = (int) constraints[i].getWeight(values);
                BoolExpr valBoolPrefix[num_constr_vars];
                valBoolPrefix[0] = (problem.var(constraints[i].getVariable(0)) == values[0]);
                //cout<<"constr "<<i<<" val "<<values<<" obtained"<<endl;
                for (size_t k = 1; k < num_constr_vars; k++) {
                    valBoolPrefix[k] = (
                        valBoolPrefix[k-1] && (problem.var(constraints[i].getVariable(k)) == values[k]));
                }
                rel(*(problem.space()), valBoolPrefix[num_constr_vars-1] >> (problem.var(numVars + i) == weight));
            }
        }
        
        
        LinIntExpr costSumPrefix[constraints.size()];
        costSumPrefix[0] = problem.var(numVars);
        for (size_t i = 1; i < constraints.size(); i++) {
            costSumPrefix[i] = (costSumPrefix[i-1] + problem.var(numVars + i));
        }
        //cout<<"sum prefs evaled"<<endl;
        problem.QIntVar(numVars + constraints.size(),0,totMaxWeight); //total cost
        rel(*(problem.space()), costSumPrefix[constraints.size()-1] == problem.var(numVars + constraints.size()));
        //cout<<"total cost constr added"<<endl;
        
        IntVarArgs branchFin(numVars + scopeVars[numVars]);
        for (int i = 0; i< numVars + scopeVars[numVars]; i++)
            branchFin[i] = problem.var(i);
        branch(*(problem.space()),branchFin,INT_VAR_SIZE_MIN(),INT_VAL_MIN());
        
        OptVar* costopt = problem.getExistential(numVars + constraints.size());
        //problem.optimize(numVars, 1, costopt)
        for (int varInd = numVars - 1; varInd >= 0; varInd--) {
            problem.optimize(varInd, optType[varInd], costopt);
        }

        problem.makeStructure();
        qecode::QCOP_solver sol(&problem);
        
        unsigned long int nodes=0;

        Strategy strategy = sol.solve(nodes);
        vector<int> solVals;
        solVals.reserve(numVars + scopeVars[numVars]);
        strategy.getAllValues(solVals);
        cout<<nodes<<" Nodes visited."<<endl;
        for (int i : solVals) cout<<i<<" ";
        cout<<endl;
        //return 0;
        std::map<variable_id_t, bool> solution;
        return solution;
    }

    /** \brief Solve the QWCSP instance using QCOP (QeCode)
     *
     * \param[in] lps The quantifiers
     *
     * \return A solution to the WCSP instance.
     */
    std::map<variable_id_t, bool> solveWCSPUsingQCOP() const noexcept
    {
        int numLevels = 2;
        bool* quants = new bool[numLevels];
        int* scopeVars = new int[numLevels];

        std::fill_n(quants, numLevels, false);
        scopeVars[0] = numVars;
        scopeVars[1] = constraints.size() + 1; //Add one for total cost
        QProblem problem(numLevels, quants, scopeVars);

        //cout<<"QP declared"<<endl;
        for (int varInd = 0; varInd < numVars; varInd++) {
            problem.QIntVar(varInd,0,1);
        }
        IntVarArgs branchOne(numVars);
        for (int i=0; i<numVars; i++) {
            branchOne[i] = problem.var(i);
        }
        branch(*(problem.space()),branchOne,INT_VAR_SIZE_MIN(),INT_VAL_MIN());
        problem.nextScope();

        cout<<"vars scope done"<<endl;
        int totMaxWeight = 0;
        for (size_t i = 0; i < constraints.size(); i++) {
            int mw = (int) constraints[i].getMaxWeight();
            problem.QIntVar(numVars + i,0,mw);
            totMaxWeight += mw;
            size_t num_constr_vars = constraints[i].getVariables().size();
            size_t num_values = 1 << num_constr_vars;
            
            cout<<"constr "<<i<<" var added"<<endl;
            for (size_t j = 0; j < num_values; ++ j) {  // add all variables
                constraint_value_t values(num_constr_vars, j);
                int weight = (int) constraints[i].getWeight(values);
                BoolExpr valBoolPrefix[num_constr_vars];
                valBoolPrefix[0] = (problem.var(constraints[i].getVariable(0)) == values[0]);
                //cout<<"constr "<<i<<" val "<<values<<" obtained"<<endl;
                for (size_t k = 1; k < num_constr_vars; k++) {
                    valBoolPrefix[k] = (
                        valBoolPrefix[k-1] && (problem.var(constraints[i].getVariable(k)) == values[k]));
                }
                rel(*(problem.space()), valBoolPrefix[num_constr_vars-1] >> (problem.var(numVars + i) == weight));
            }
        }

        LinIntExpr costSumPrefix[constraints.size()];
        costSumPrefix[0] = problem.var(numVars);
        for (size_t i = 1; i < constraints.size(); i++) {
            costSumPrefix[i] = (costSumPrefix[i-1] + problem.var(numVars + i));
        }
        cout<<"sum prefs evaled"<<endl;
        problem.QIntVar(numVars + constraints.size(),0,totMaxWeight); //total cost
        rel(*(problem.space()), costSumPrefix[constraints.size()-1] == problem.var(numVars + constraints.size()));
        cout<<"total cost constr added"<<endl;
        
        IntVarArgs branchFin(numVars + scopeVars[1]);
        for (int i = 0; i< numVars + scopeVars[1]; i++)
            branchFin[i] = problem.var(i);
        branch(*(problem.space()),branchFin,INT_VAR_SIZE_MIN(),INT_VAL_MIN());
        cout<<"total cost constr added1"<<endl;
        
        OptVar* costopt = problem.getExistential(numVars + constraints.size());
        //problem.optimize(numVars, 1, costopt)
        problem.optimize(0, 1, costopt);
        cout<<"total cost constr added2"<<endl;

        problem.makeStructure();
        qecode::QCOP_solver sol(&problem);
        
        unsigned long int nodes=0;
        cout<<"total cost constr added3"<<endl;

        Strategy strategy = sol.solve(nodes);
        vector<int> solVals;
        solVals.reserve(numVars + scopeVars[1]);
        cout<<"total cost constr added4"<<endl;

        strategy.getAllValues(solVals);
        cout<<nodes<<" Nodes visited."<<endl;
        for (int i : solVals) cout<<i<<" ";
        cout<<endl;
        //return 0;
        std::map<variable_id_t, bool> solution;
        return solution;
    }

    void test() {
        int NCustomer = 9;
        int NArc = 5;
        int*c = new int[NCustomer*NArc]; // [NArc*i+j] : Initial cost for client i to reach way j
        int* d = new int[NCustomer]; // d[i] : amount of data client i wants to transmit
        int* u = new int[NCustomer]; // u[i] : maximum cost client i will accept
        int max = 19;

        c[0*NArc+0]=5;    c[0*NArc+1]=1;    c[0*NArc+2]=8;    c[0*NArc+3]=6;    c[0*NArc+4]=6;
        c[1*NArc+0]=2;    c[1*NArc+1]=8;    c[1*NArc+2]=6;    c[1*NArc+3]=3;    c[1*NArc+4]=1;
        c[2*NArc+0]=0;    c[2*NArc+1]=8;    c[2*NArc+2]=0;    c[2*NArc+3]=8;    c[2*NArc+4]=6;
        c[3*NArc+0]=2;    c[3*NArc+1]=4;    c[3*NArc+2]=7;    c[3*NArc+3]=5;    c[3*NArc+4]=8;
        c[4*NArc+0]=8;    c[4*NArc+1]=2;    c[4*NArc+2]=4;    c[4*NArc+3]=3;    c[4*NArc+4]=4;
        c[5*NArc+0]=5;    c[5*NArc+1]=4;    c[5*NArc+2]=6;    c[5*NArc+3]=4;    c[5*NArc+4]=1;
        c[6*NArc+0]=5;    c[6*NArc+1]=6;    c[6*NArc+2]=7;    c[6*NArc+3]=4;    c[6*NArc+4]=8;
        c[7*NArc+0]=5;    c[7*NArc+1]=8;    c[7*NArc+2]=1;    c[7*NArc+3]=3;    c[7*NArc+4]=6;
        c[8*NArc+0]=7;    c[8*NArc+1]=8;    c[8*NArc+2]=8;    c[8*NArc+3]=3;    c[8*NArc+4]=5;

        d[0]=7;    d[1]=16;    d[2]=16;    d[3]=19;    d[4]=18;    d[5]=13;    d[6]=8;    d[7]=11;    d[8]=18;
        u[0]=122;    u[1]=190;    u[2]=113;    u[3]=285;    u[4]=247;    u[5]=255;    u[6]=143;    u[7]=121;    u[8]=139;
        IntArgs carg(NCustomer*NArc,c); // Copy c in an IntArgs for further constraint posting
        IntArgs darg(NCustomer,d); // Copy d in an IntArgs for further constraint posting
        IntArgs uarg(NCustomer,u); // Copy u in an IntArgs for further constraint posting

        bool q[] = {false,true,false};
        int* nv = new int[3];
        nv[0]=NArc;
        nv[1]=1;
        nv[2]=9;

        QProblem problem(3,q,nv);
        int postarfs[] = {3,7,11,15,19};
        IntSet thePossibleTariffs(postarfs,5);
        for (int i=0;i<NArc;i++)
            problem.QIntVar(i,thePossibleTariffs); // tariff for way i
        IntVarArgs branch1(NArc);
        for (int i=0;i<NArc;i++)
            branch1[i] = problem.var(i);
        branch(*(problem.space()),branch1,INT_VAR_SIZE_MIN(),INT_VAL_MIN());
        problem.nextScope();

        problem.QIntVar(NArc,0,NCustomer-1); // k
        IntVarArgs branch2(NArc+1);
        for (int i=0;i<NArc+1;i++)
            branch2[i] = problem.var(i);
        branch(*(problem.space()),branch2,INT_VAR_SIZE_MIN(),INT_VAL_MIN());
        problem.nextScope();

        problem.QIntVar(NArc+1,0,NArc-1); // a
        problem.QIntVar(NArc+2,0,1000000); // cost
        problem.QIntVar(NArc+3,0,1000000); // Income
        IntVar a(problem.var(NArc+1));
        IntVar cost(problem.var(NArc+2));
        IntVar income(problem.var(NArc+3));
        problem.QIntVar(NArc+4,0,1000000);
        problem.QIntVar(NArc+5,0,1000000);
        problem.QIntVar(NArc+6,0,1000000);
        problem.QIntVar(NArc+7,0,1000000);
        problem.QIntVar(NArc+8,0,1000000);
        problem.QIntVar(NArc+9,0,1000000);
        IntVar aux1(problem.var(NArc+4)); // k* NArc + a
        IntVar aux2(problem.var(NArc+5)); // c[k*Narc+a]
        IntVar aux3(problem.var(NArc+6)); // t[a]
        IntVar aux4(problem.var(NArc+7)); // d[k]
        IntVar aux5(problem.var(NArc+8)); // c[]+t[]
        IntVar aux6(problem.var(NArc+9)); // u[k]
        IntVar k(problem.var(NArc));
        rel(*(problem.space()), aux1 == ( NArc * k + a) );
        element(*(problem.space()),carg,aux1,aux2);
        IntVarArgs t(NArc);
        for (int i=0;i<NArc;i++) t[i]=problem.var(i);
        element(*(problem.space()),t,a,aux3);
        element(*(problem.space()),darg,k,aux4);
        rel(*(problem.space()), aux5 == aux2 + aux3);
        mult(*(problem.space()),aux5,aux4,cost); // cost = aux5 * aux4
        mult(*(problem.space()),aux3,aux4,income);
        element(*(problem.space()),uarg,k,aux6);
        rel(*(problem.space()),cost <= aux6);

        IntVarArgs branch3(NArc+10);
        for (int i=0;i<NArc+10;i++)
            branch3[i] = problem.var(i);
        branch(*(problem.space()),branch3,INT_VAR_SIZE_MIN(),INT_VAL_MIN());

        OptVar* costopt = problem.getExistential(NArc+2);
        OptVar* incomeopt = problem.getExistential(NArc+3);
        problem.optimize(2,1,costopt); // at scope 2, we minimize (1) the variable cost
        AggregatorSum somme;
        OptVar* sumvar = problem.getAggregate(1,incomeopt,&somme);
        problem.optimize(0,2,sumvar);


        QCOP_solver sol(&problem);
        unsigned long int nodes=0;

        Strategy strategy = sol.solve(nodes);
        std::cout<<nodes<<" Nodes visited."<<std::endl;
        std::cout<<strategy<<std::endl;
        return;
    }
};

template <class VarIdType, class WeightType, class ConstraintValueType>
void WCSPInstance<VarIdType, WeightType, ConstraintValueType>::loadDimacs(std::istream& f)
{
    std::string line;
    std::string tmp;
    size_t tmpint, nc; // number of variables, number of constraints

    // first line: problem_name number_of_variables max_domain_size number_of_constraints
    // global_upper_bound
    std::getline(f, line);
    std::istringstream ss(line);
    ss >> tmp; // name
    ss >> numVars;
    ss >> tmp; // max domain size, must be 2
    if (tmp != "2")
        throw std::domain_error("Domain size must be 2");
    ss >> nc;
    ss >> tmp; // ignore global upper bound
    optType.reserve(numVars);

    // Opt types
    std::getline(f, line);
    ss.clear();
    ss.str(line);
    for (size_t i = 0; i < numVars; ++ i)
    {
        ss >> tmpint;
        optType.push_back(tmpint);
        if (optType[i] != 1 && optType[i] != 2)
            throw std::domain_error("Opt type must be 1 or 2"); // QCOP syntax for min/max resp.
    }

    // domain size of all variables
    std::getline(f, line);
    ss.clear();
    ss.str(line);
    for (size_t i = 0; i < numVars; ++ i)
    {
        int tmp;
        ss >> tmp;
        if (tmp != 2)
            throw std::domain_error("Domain size must be 2");
    }

    constraints.clear();
    constraints.resize(nc);

    // iterate over constraints
    for (size_t i = 0; i < nc; ++ i)
    {
        size_t arity;
        weight_t default_cost;
        size_t ntuples; // number of tuples not having the default cost
        std::getline(f, line);
        ss.clear();
        ss.str(line);

        ss >> arity;

        // load the variables in the constraint
        std::vector<variable_id_t> variables;
        variables.reserve(arity);
        for (size_t j = 0; j < arity; ++ j)
        {
            variable_id_t vid;
            ss >> vid;
            variables.push_back(vid);
        }
        constraints[i].setVariables(std::move(variables));
        ss >> default_cost;
        if (default_cost > std::abs(1e-6))   // non-zero default cost
        {
            constraint_value_t values;
            values.resize(constraints[i].getVariables().size());
            unsigned long te = (1u << constraints[i].getVariables().size());
            for (unsigned long l = 0; l < te; ++ l)
            {
                for (unsigned long k = 0; k < values.size(); ++ k)
                    values.set(k, (l & (1u << k)) ? 1 : 0);
                constraints[i].setWeight(values, default_cost);
            }
        }
        ss >> ntuples;

        // load the entries in the constraints
        for (size_t j = 0; j < ntuples; ++ j)
        {
            weight_t cost;
            constraint_value_t values;
            values.resize(constraints[i].getVariables().size());
            std::getline(f, line);
            ss.clear();
            ss.str(line);
            for (size_t k = 0; k < values.size(); ++ k)
            {
                int val;
                ss >> val;
                values.set(k, val ? 1 : 0);
            }
            ss >> cost;
            constraints[i].setWeight(std::move(values), cost);
        }
        totMaxWeight += constraints[i].getMaxWeight();
    }
}

template <class VarIdType, class WeightType, class ConstraintValueType>
void WCSPInstance<VarIdType, WeightType, ConstraintValueType>::loadUAI(std::istream& f)
{
    std::string line;
    std::string tmp;
    size_t nv, nc; // number of variables, number of constraints

    // first line: MARKOV, just ignore it
    std::getline(f, tmp);

    // second line: number of variables
    {
        std::getline(f, line);
        std::istringstream ss(line);
        ss >> nv; // number of variables
    }

    // third line: arity
    {
        std::getline(f, line);
        std::istringstream ss(line);
        for (size_t i = 0; i < nv; ++ i)
        {
            ss >> tmp; // max domain size, must be 2
            if (tmp != "2")
                throw std::domain_error("Domain size must be 2");
        }
    }

    // fourth line: number of cliques
    {
        std::getline(f, line);
        std::istringstream ss(line);
        ss >> nc;
        constraints.clear();
        constraints.resize(nc);
    }

    // iterate over constraints
    for (size_t i = 0; i < nc; ++ i)
    {
        size_t arity;
        std::getline(f, line);
        std::istringstream ss(line);

        ss >> arity;

        // Load the variables in the constraint. We load it reversely because UAI format arrange
        // their constraint value reversely as we do.
        std::vector<variable_id_t> variables(arity);

        for (size_t j = 0; j < arity; ++ j)
        {
            variable_id_t vid;
            ss >> vid;
            variables[arity-j-1] = vid;
        }

        constraints[i].setVariables(std::move(variables));
    }

    // iterate over constraints
    for (size_t i = 0; i < nc; ++ i)
    {
        size_t ntuples;
        f >> ntuples;  // number of entries

        // load the entries in the constraints
        std::vector<weight_t> costs(ntuples);
        for (size_t j = 0; j < ntuples; ++ j)
            f >> costs[j];

        // normalization constant
        weight_t sum(0);
        for (auto& c : costs)
            sum += c;

        for (size_t j = 0; j < ntuples; ++ j)
        {
            costs[j] /= sum;
            costs[j] = -std::log(costs[j]);
            if (!std::isfinite(costs[j]))
                costs[j] = 1e6;
            constraint_value_t values(constraints[i].getVariables().size(), j);
            constraints[i].setWeight(std::move(values), costs[j]);
        }
    }
}

template <class VarIdType, class WeightType, class ConstraintValueType>
void WCSPInstance<VarIdType, WeightType, ConstraintValueType>::load(std::istream& f)
{
    using namespace boost::property_tree;

    ptree pt;
    read_json(f, pt);

    auto& constraints_tree = pt.get_child("constraints");
    constraints.clear();
    constraints.reserve(constraints_tree.size());
    // iterate over all constraints
    for (auto& constraint : pt.get_child("constraints"))
    {
        auto& cons_tree = constraint.second;

        // append one more element to the constraint list
        constraints.resize(constraints.size() + 1);
        constraint_t& cons = constraints.back();

        std::vector<variable_id_t> variables;

        // iterate over variables
        for (auto& v : cons_tree.get_child("variables"))
            variables.push_back(v.second.get_value<variable_id_t>());

        cons.setVariables(std::move(variables));

        // iterate over weights
        for (auto& w : cons_tree.get_child("weights"))
        {
            constraint_value_t values;
            values.resize(cons.getVariables().size());
            size_t i = 0;

            for (auto it = w.second.begin(); it != w.second.end(); ++ it)
            {
                if (std::next(it) == w.second.end()) // the last one is the weight
                {
                    cons.setWeight(std::move(values), it->second.get_value<weight_t>());
                    break;
                }

                values.set(i++, it->second.get_value<bool>());
            }
        }
    }
}

/** \brief Write a Constraint object to a stream in a human-readable form.
 *
 * \param[in] o The stream to write to.
 *
 * \param[in] c The constraint object to be write to \p o.
 *
 * \return The stream \p o.
 */
template <class ...T>
std::ostream& operator << (std::ostream& o, const Constraint<T...>& c)
{
    boost::property_tree::ptree t;
    c.toPropertyTree(t);
    boost::property_tree::write_json(o, t, true);
    return o;
}

#endif /* WCSPINSTANCE_H_ */

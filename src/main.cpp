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

#include <iostream>
#include <string.h>
#include <regex>
#include <boost/algorithm/string.hpp>
#include <boost/program_options.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#if (BOOST_VERSION < 106500)
#define NO_LZMA
#endif

#ifndef NO_LZMA
#include <boost/iostreams/filter/lzma.hpp>
#endif

#include "global.h"
#include "Recorder.h"
#include "RunningTime.h"
#include "WCSPInstance.h"
#include "ConstraintCompositeGraph.h"
//#include "LinearProgramSolverGurobi.h"
//#include "MWVCSolverLinearProgramming.h"
//#include "MWVCSolverMessagePassing.h"
//#include "MWVCSolverFastWVC.h"
//#include "KernelizerLinearProgramming.h"

int main(int argc, const char **argv)
{
    // Parsing command line arguments and print the help message if needed
    namespace po = boost::program_options;
    po::variables_map vm;
    char mwvc_solver;
    double eps, delta, delta1;
    std::string ccg_out_file;
    std::string parameters; // parameters of some MWVC solvers
    char file_format;
    double initial_solution_time_limit, time_limit;
    double ttime;

    {
        po::options_description desc("Options");
        desc.add_options()
            ("help,h", "print help message")
            ("ccg,c", po::value<>(&ccg_out_file)->default_value(""),
             "specify the file that the output CCG should be written into (default is stdout)")
            ("file-format,f", po::value<>(&file_format)->default_value('d')->notifier(
                [](char x){
                    if (x != 'd' && x != 'u')
                        throw po::validation_error(
                            po::validation_error::invalid_option_value,
                            "file-format", std::string(1, x)); 
                }),
             "specify the input file format:\n"
             "d: DIMACS\n"
             "u: UAI\n")
            ("ccg-only,g", "print the CCG only without solving the MWVC problem on it")
            ("no-kernelization,k", "don't kernelize")
            ("kernelization-only,K", "exit after kernelization")
            ("message-passing,M", "do message passing on the original problem directly")
            ("linear-programming,L", "use linear programming to solve the original problem directly")
            ("mwvc-solver,m", po::value<>(&mwvc_solver)->default_value('l')->notifier(
                [](char x){
                    if (std::string("lamfh").find(x) == std::string::npos)
                        throw po::validation_error(
                            po::validation_error::invalid_option_value,
                            "mwvc-solver", std::string(1, x));            
                }),
             "MWVC solver:\n"
             "l: linear program solver\n"
             "m: message passing\n"
             "f: fastWVC\n"
             "h: hybrid solver where the initial solution is produced by fastWVC and then the linear program solver is used to find the optimum\n")
            ("parameters,p", po::value<>(&parameters)->default_value(""))
            ("initial-solution-time-limit,i",
             po::value<double>(&initial_solution_time_limit)->default_value(std::numeric_limits<double>::max()),
             "for the hybrid solvers only, specify the time limit for finding the initial solution")
            ("time-limit,t",
             po::value<double>(&time_limit)->default_value(std::numeric_limits<double>::max()),
             "specify the time limit in seconds")
            ("record-results,r", po::value<std::string>()->default_value("")->implicit_value(""), 
             "record the results on file");

        po::options_description desc_hidden("Hidden Options");
        desc_hidden.add_options()
            ("input-file", po::value<std::string>()->required(), "specify input file");

        po::positional_options_description desc_p;
        desc_p.add("input-file", 1);

        po::options_description all_options("Options");
        all_options.add(desc).add(desc_hidden);
        try
        {
            po::store(po::command_line_parser(argc, argv).options(all_options).positional(desc_p).run(), vm);

            if (vm.count("help")) // print help message
            {
                std::cerr << "Usage: " << argv[0] << " [options] [input-file]" << std::endl;
                std::cerr << desc << std::endl;
                return 1;
            }

            if (vm.count("input-file") != 1)
                throw po::error_with_no_option_name("An input file must be specified.");

            po::notify(vm);

            // validate p
            if (mwvc_solver == 'a')
            {
                std::vector<std::string> p;
                boost::split(p, parameters, boost::is_any_of(","));
                if (p.size() != 3)
                    throw po::validation_error(
                        po::validation_error::invalid_option_value, "parameters",
                        "3 parameters (epsilon, delta, delta1) separated with commas are expected.");
                eps = std::stod(p.at(0));
                delta = std::stod(p.at(1));
                delta1 = std::stod(p.at(2));
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
            std::cerr << "Run this program with the \"--help\" option to see help." << std::endl;
            return 1;
        }
    }

    std::string input_file = vm["input-file"].as<std::string>();

    // find compression format
    std::string cformat = "";
    if (std::regex_search(input_file, std::regex(".gz$")))
        cformat = "gz";
    else if (std::regex_search(input_file, std::regex(".xz$")))
        cformat = "xz";

    // decompress input if compressed
    std::ifstream ifs(input_file, (std::ios_base::in | std::ios_base::binary));

    if (!ifs)
    {
        std::cerr << "Unable to read file \"" << vm["input-file"].as<std::string>() << '"' << std::endl;
        return 2;
    }

    boost::iostreams::filtering_streambuf<boost::iostreams::input> dstream_buf;
    if (cformat == "gz")
    {
        dstream_buf.push(boost::iostreams::gzip_decompressor());
    }
    else if (cformat == "xz")
    {
#ifndef NO_LZMA
        dstream_buf.push(boost::iostreams::lzma_decompressor());
#else
        cerr << "Error: compiling with Boost version 1.65 or higher is needed to allow to read xz compressed wcsp format files." << endl;
#endif
    }
    dstream_buf.push(ifs);
    std::istream in(&dstream_buf);

    WCSPInstance<>::Format fformat;
    switch (file_format)
    {
    case 'd':
        fformat = WCSPInstance<>::Format::DIMACS;
        break;
    case 'u':
        fformat = WCSPInstance<>::Format::UAI;
        break;
    }

    RunningTime::GetInstance().setStartingTime(std::chrono::high_resolution_clock::now());

    WCSPInstance<> instance(in, fformat);

    if (time_limit != std::numeric_limits<double>::max())  // time limit is set
        RunningTime::GetInstance().setTimeLimit(std::chrono::duration<double>(time_limit));
    
    if (mwvc_solver == 'h' && initial_solution_time_limit != std::numeric_limits<double>::max())  // time limit is set
        RunningTime::GetInstance().setTimeLimit(std::chrono::duration<double>(initial_solution_time_limit));

    if (vm.count("record-results"))
    {
        std::string results_file_name = vm["record-results"].as<std::string>();
        Recorder::getInstance().set_of_path(input_file, results_file_name);
    }

    if (vm.count("message-passing"))
    {
        Recorder::getInstance().set_mode("M");
        auto assignments = instance.solveUsingMessagePassing(1e-6);
        ttime = RunningTime::GetInstance().getCurrentRunningTime();
        Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
        Recorder::getInstance().print_results(assignments);
        if (vm.count("no-logging") != 1)
            Recorder::getInstance().store_results();
        return 0;
    }
    /*
    if (vm.count("linear-programming"))
    {
        Recorder::getInstance().set_mode("L");
        LinearProgramSolverGurobi lps;
        LPObjUpdateHistory *updateHistory = new LPObjUpdateHistory();
        lps.setUpdateHistory(updateHistory);
        auto assignments = instance.solveUsingLinearProgramming(lps);
        //TODO: handle timeout exception
        while (updateHistory->has_next_update())
        {
            auto ttime = updateHistory->next_update();
            Recorder::getInstance().record_best_bounds(
                updateHistory->getLB(),
                updateHistory->getUB(),
                ttime);
        }
        Recorder::getInstance().set_opt_proven();
        Recorder::getInstance().print_results(assignments);
        if (vm.count("record-results"))
            Recorder::getInstance().store_results();
        return 0;
    }*/

    ConstraintCompositeGraph<> ccg;

    WCSPInstance<>::constraint_t::Polynomial p;
    for (const auto &c : instance.getConstraints())
        c.toPolynomial(p);
    // s is the remnant weight
    ConstraintCompositeGraph<>::weight_t s = ccg.addPolynomial(p);

    // The reason that we use a map instead of a set to represent the vertex cover is that after
    // each step, we can see what variables have been assigned and what have not.
    std::map<ConstraintCompositeGraph<>::variable_id_t, bool> assignments;

    //ccg.simplify(assignments);
    

    //std::cout << "Variables simplified out: " << assignments.size() << std::endl;
    std::cout << "==========================" << std::endl;
    auto stats = ccg.getStatistics();
   // std::cout << ccg << std::endl;
    //ccg.toDimacs(std::cout, true);
    
    std::cout << "Number of variables: " << stats[0] << std::endl;
    std::cout << "Number of type 1 auxiliary variables: " << stats[1] << std::endl;
    std::cout << "Number of type 2 auxiliary variables: " << stats[2] << std::endl;
    std::cout << "=====test qcop========" << std::endl;
    //instance.test();
    //int opt[] = {1,1,1};
    instance.solveWCSPUsingQCOP();
    std::cout << "=====test ccg qcop========" << std::endl;
    ccg.solveWCSPCCGUsingQCOP(instance.getOptType());
    if (!ccg_out_file.empty())
    {
        std::ofstream ofs(ccg_out_file);
        if (!ofs)
        {
            std::cerr << "Failed to open file \"" << ccg_out_file << '"' << std::endl;
            abort();
        }
        ofs.exceptions(std::ofstream::failbit | std::ofstream::badbit);
        ccg.toDimacs(ofs, true);
    }

    std::cout << std::setprecision(std::numeric_limits<decltype(s)>::digits10 + 1) << "s = " << s << std::endl;

    if (vm.count("ccg-only")) // don't solve the MWVC problem
        return 0;

    /*
    do
    {
        auto g = *ccg.getGraph();
        if (vm.count("no-kernelization"))
        {
            std::string mode = "-no-kern";
            mode.insert(0, 1, mwvc_solver);
            Recorder::getInstance().set_mode(mode);
            std::cout << "================================" << std::endl
                      << "|| No kernelization performed ||" << std::endl
                      << "================================" << std::endl;
        }
        else
        {
            std::string mode = "-kern";
            mode.insert(0, 1, mwvc_solver);
            Recorder::getInstance().set_mode(mode);
            size_t cur_assignment_size = -1;

            for (size_t i = 1; cur_assignment_size != assignments.size(); ++i)
            {
                cur_assignment_size = assignments.size();

                KernelizerLinearProgramming<> klp(new LinearProgramSolverGurobi());
                klp.kernelize(g, assignments);
                std::cout << "==========================" << std::endl;
                std::cout << "After the " << i << "th kernelization, number of variables resolved: " << assignments.size() << std::endl;
                std::cout << "After the " << i << "th kernelization, number of variables left: " << ccg.getNumberOfVariables() - assignments.size() << std::endl;
                std::cout << "==========================" << std::endl;

                // Don't continue if all variables have been factored out.
                if (assignments.size() >= ccg.getNumberOfVariables())
                    break;
            }

            if (assignments.size() >= ccg.getNumberOfVariables())
            {
                ttime = RunningTime::GetInstance().getCurrentRunningTime();
                Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
                break;
            }
        }

        if (vm.count("kernelization-only"))
            std::exit(0);

        switch (mwvc_solver)
        {
        case 'l':
        {
            LinearProgramSolverGurobi *lps = new LinearProgramSolverGurobi();
            MWVCSolverLinearProgramming<> mwvc_solver_lp(lps, new LPAssignmentsUpdateHistory());
            try
            {
                mwvc_solver_lp.solve(g, assignments);
            }
            catch (LinearProgramSolver::TimeOutException e)
            {
                while (mwvc_solver_lp.has_next_update())
                {
                    ttime = std::min(
                        mwvc_solver_lp.next_update(g, assignments),
                        RunningTime::GetInstance().getTimeLimit().count());
                    Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
                }
                if (vm.count("no-logging") != 1)
                    Recorder::getInstance().store_results();
                throw e;
            }
            Recorder::getInstance().set_opt_proven();
            while (mwvc_solver_lp.has_next_update())
            {
                ttime = mwvc_solver_lp.next_update(g, assignments);
                Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
            }
            break;
        }
        case 'm':
        {
            MWVCSolverMessagePassing<> mwvc_solver_mp(1e-6);
            mwvc_solver_mp.solve(g, assignments);
            ttime = RunningTime::GetInstance().getCurrentRunningTime();
            Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
            break;
        }
        case 'f':
        {
            MWVCSolverFastWVC<> mwvc_solver_fastwvc(new FastWVCSolver<>(), new VCUpdateHistory());
            mwvc_solver_fastwvc.solve(g, assignments);
            while (mwvc_solver_fastwvc.has_next_update())
            {
                ttime = mwvc_solver_fastwvc.next_update(g, assignments);
                Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
            }
            break;
        }
        case 'h':
        {
            bool initialSolutionFound = false;
            MWVCSolverFastWVC<> mwvc_solver_fastwvc(new FastWVCSolver<>(), new VCUpdateHistory());
            mwvc_solver_fastwvc.solve(g, assignments);
            while (mwvc_solver_fastwvc.has_next_update())
            {
                initialSolutionFound = true;
                ttime = mwvc_solver_fastwvc.next_update(g, assignments);
                Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
            }
            Recorder::getInstance().set_initial_sol_idx();
            RunningTime::GetInstance().setTimeLimit(std::chrono::duration<double>(time_limit));
            MWVCSolverLinearProgramming<> mwvc_solver_lp(new LinearProgramSolverGurobi(), new LPAssignmentsUpdateHistory());
            if (initialSolutionFound)
            {
                int v_num = num_vertices(g);
                double *initialSolution = new double[v_num];
                for (int v = 0; v < v_num; ++v)
                    initialSolution[v] = mwvc_solver_fastwvc.is_in_best_vc(v + 1);
                mwvc_solver_lp.setInitialSolution(initialSolution);
            }
            try
            {
                mwvc_solver_lp.solve(g, assignments);
            }
            catch (LinearProgramSolver::TimeOutException e)
            {
                while (mwvc_solver_lp.has_next_update())
                {
                    ttime = std::min(
                        mwvc_solver_lp.next_update(g, assignments),
                        RunningTime::GetInstance().getTimeLimit().count());
                    Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
                }
                if (vm.count("no-logging") != 1)
                    Recorder::getInstance().store_results();
                throw e;
            }
            Recorder::getInstance().set_opt_proven();
            while (mwvc_solver_lp.has_next_update())
            {
                ttime = mwvc_solver_lp.next_update(g, assignments);
                Recorder::getInstance().record_best_assignments(instance, assignments, ttime);
            }
            break;
        }
        }
    } while (0);
    */

    std::cout << "=================================================" << std::endl;

    Recorder::getInstance().print_results(assignments);
    if (vm.count("record-results"))
        Recorder::getInstance().store_results();

    return 0;
}

#ifndef RECORDER_H_
#define RECORDER_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <regex>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "RunningTime.h"

#define CONVERGED 1
#define NOT_CONVERGED 2
#define DOT "_dot_"

class Recorder
{

    private:
        int mp_outcome = 0;
        int initial_sol_idx = -1;
        bool opt_proven = false;
        std::string mode = "L";
        std::string test_rel_path = "tests/test.wcsp";
        std::string of_path = "results/results.json"; 
        std::vector<double> ub_updates;
        std::vector<double> lb_updates;
        std::vector<double> time_updates;
        Recorder() {}

    public:

        void set_mode(std::string mode) { this->mode = mode; }

        void set_of_path(std::string ipath, std::string of_path = "")
        {
            this->test_rel_path = remove_path_to_root(ipath);

            if (of_path.size() > 0)
                this->of_path = of_path;
            else 
            {
                std::string results_dir_path, results_log_file_name;
                if (test_rel_path.find_first_of("/") == std::string::npos ||
                    test_rel_path.find_first_of("/") == test_rel_path.find_last_of("/")) {
                    results_dir_path = "results/";
                    results_log_file_name = "results.json";
                } else {
                    results_dir_path = path_to_root(ipath) + "results/";
                    boost::filesystem::create_directories(results_dir_path);
                    unsigned long start_idx = test_rel_path.find_first_of("/") + 1;
                    unsigned long end_idx = test_rel_path.find_last_of("/");
                    results_log_file_name = test_rel_path.substr(start_idx, end_idx - start_idx) + ".json";
                    std::replace(results_log_file_name.begin(), results_log_file_name.end(), '/', '_');
                }

                this->of_path = results_dir_path + results_log_file_name;
            }
        }

        void set_mp_converged(bool converged) { this->mp_outcome = converged ? CONVERGED : NOT_CONVERGED; }

        void set_opt_proven() { opt_proven = true; }

        void set_initial_sol_idx() { this->initial_sol_idx = ub_updates.size() - 1; }

        template <class T2>
        void print_results(const T2& assignments)
        {
            if (RunningTime::GetInstance().isTimeOut())
                std::cout << "Timeout solution" << std::endl;
            if (ub_updates.size() > 0) 
            {
                std::cout << "Best assignments:" << std::endl;
                std::cout << "ID\tassignment" << std::endl;
                for (auto a : assignments)
                    std::cout << a.first << '\t' << a.second << std::endl;
                std::cout << "Optimal value: " << ub_updates.back() << std::endl;
            }
        }

        template <class T1, class T2>
        void record_best_assignments(const T1 &instance, const T2 &assignments, const double timestamp)
        {   
            if (lb_updates.size() > 0)
                lb_updates.push_back(lb_updates.back());
            auto newUb = instance.computeTotalWeight(assignments);
            ub_updates.push_back(newUb);
            time_updates.push_back(timestamp);
        }

        void record_best_bounds(const double newLb, const double newUb, const double timestamp)
        {
            lb_updates.push_back(newLb);
            ub_updates.push_back(newUb);
            time_updates.push_back(timestamp);
        }

        void store_results()
        {    
            // read existing data
            boost::property_tree::ptree pt;
            try {
                boost::property_tree::json_parser::read_json(of_path, pt);
            } catch (const boost::property_tree::json_parser_error& e1) {}

            std::string test_name = test_rel_path.substr(test_rel_path.find_last_of("/") + 1);
            std::string test_name_no_ext = test_name.substr(
                0,
                std::min(
                    test_name.rfind(".wcsp"),
                    test_name.rfind(".uai")
                )
            );
            test_name_no_ext = replace(test_name_no_ext, ".", DOT);
            auto prefix = test_name_no_ext + "." + mode; 

            pt.put(prefix + ".testRelPath", test_rel_path);

            std::string status;
            if (ub_updates.size() > 0)
            {
                pt.put(prefix + ".bestObj", ub_updates.back());
                pt.put(prefix + ".time", time_updates.back());
                status = opt_proven ? "SC" : "S";
                
            } else {  
                status = "UNK";
            }

            pt.put(prefix + ".status", status);

            if (mp_outcome == CONVERGED || mp_outcome == NOT_CONVERGED) 
                pt.put(prefix + ".converged", mp_outcome == CONVERGED);

            if (ub_updates.size() > 0)
            {
                if (initial_sol_idx >= 0)
                    pt.put(prefix + ".initalSolutionIdx", initial_sol_idx);

                boost::property_tree::ptree ptree_bounds, ptree_times;

                if (lb_updates.size() == 0) 
                {
                    for (auto w: ub_updates) 
                    {
                        boost::property_tree::ptree ptree_obj;
                        ptree_obj.put("", w);
                        ptree_bounds.push_back(std::make_pair("", ptree_obj));
                    }
                } 
                else 
                {
                    for(unsigned i = 0; i < ub_updates.size(); i++) {
                        boost::property_tree::ptree ptree_node_bounds;
                        boost::property_tree::ptree lb;
                        boost::property_tree::ptree ub;
                        lb.put("", lb_updates[i]);
                        ub.put("", ub_updates[i]);
                        ptree_node_bounds.push_back(std::make_pair("", lb));
                        ptree_node_bounds.push_back(std::make_pair("", ub));
                        ptree_bounds.push_back(std::make_pair("", ptree_node_bounds));
                    }
                }

                for (auto t: time_updates) 
                {
                    boost::property_tree::ptree ptree_time;
                    ptree_time.put("", t);
                    ptree_times.push_back(std::make_pair("", ptree_time));
                }

                pt.put_child(prefix + ".boundsHistory", ptree_bounds);
                pt.put_child(prefix + ".timeHistory", ptree_times);
            }

            write_json(of_path, pt);
        }

        static Recorder& getInstance()
        {
            static std::unique_ptr<Recorder> r;   // the only Recorder instance in the world

            if (!r)   // first time, initialize.
                r.reset(new Recorder());

            return *r;
        }
    
    private:
        std::string path_to_root(std::string path)
        {
            auto split_on = path.rfind("./") == std::string::npos ? 0 : path.rfind("./") + 2;
            return path.substr(0, split_on);
        }

        std::string remove_path_to_root(std::string path)
        {
            auto split_on = path.rfind("./") == std::string::npos ? 0 : path.rfind("./") + 2;
            return path.substr(split_on);
        }

        std::string replace(std::string s, std::string old_v, std::string new_v)
        {
            std::size_t begin = s.find(old_v);
            while (begin != std::string::npos)
            {
                s.replace(begin, old_v.size(), new_v);
                begin = s.find(old_v);
            }
            return s;
        }

        void write_json(const std::string & path, const boost::property_tree::ptree & ptree)
        {
            std::ostringstream oss;
            boost::property_tree::json_parser::write_json(oss, ptree);
            std::regex reg("\\\"([0-9]+\\.{0,1}[0-9]*)\\\"");
            std::string result = std::regex_replace(oss.str(), reg, "$1");  //replace string double with double
                                                                            //(boost does not provide full
                                                                            //compatibility with json, yet)
            result = replace(result, DOT, ".");
            std::ofstream file;
            file.open(path);
            file << result;
            file.close();
        }
    };

#endif //RECORDER_H_
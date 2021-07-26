//   Copyright 2020 Robert P. Rambo
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.


#include <iostream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <sys/stat.h>
#include <sastools/include/vector3.h>
#include "src/Aligner.h"

namespace po = boost::program_options;
namespace {
    const size_t ERROR_IN_COMMAND_LINE = 1;
    const size_t SUCCESS = 0;
    const size_t ERROR_UNHANDLED_EXCEPTION = 2;
}


bool fileExists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}



int main(int argc, char** argv) {


    std::string ref, tar, fileList;
    bool allowMirror = false; // po::bool_switch defaults to false
    bool useBackbone = false;
    bool remap = false;
    bool useCPD=false;
    bool useRandom = false;

    std::string descText =
            "\n  To align a set of PDB files : aligner --list files.inp\n";
    descText += "\n     To align against a model : aligner -r reference.pdb -t target.pdb ";
    descText += "\n                                                   ";
    descText += "\n     files.inp should be a list of pdb files to align";
    descText += "\n";
    descText += "\n   Program uses coherent point drift (CPD) algorithm for point set registration";
    descText += "\n   C++ library by Gadomski, P.J. (December 2016).";
    descText += "\n   CPD is default and can be switched to SUPCOMB like method with -c flag";
    descText += "\n   Mirror images are always checked by default, use --noMirror to switch off";
    descText += "\n   ";


    po::options_description desc(descText);

    desc.add_options()
            ("help,h", "Print help messages")
            ("list,l", po::value<std::string>(&fileList), "List of PDB files to align for map calculation")
            ("reference,r", po::value<std::string>(&ref), "reference PDB file to align to")
            ("target,t", po::value<std::string>(&tar), "PDB coordinate file that will be moving")
            ("noMirror,n", po::bool_switch(&allowMirror), "Reduce bead radius defining the lattice")
            ("noCPD,c", po::bool_switch(&useCPD), "Default is CPD, flag switches to alternate method")
            ("useBackbone,b", po::bool_switch(&useBackbone), "align using only backbone atoms")
            ("remap,m", po::bool_switch(&remap), "remap aligned set to minimum spacing grid")
            ("fast,f", po::bool_switch(&useRandom), "random sampling of initial poses instead of spherical grid search")
            ;
    po::positional_options_description positionalOptions;
    positionalOptions.add("dat", 1);
    po::variables_map vm;

    try {

        po::store(po::command_line_parser(argc, argv).options(desc).positional(positionalOptions).run(), vm);

        // checking options list
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n";
            return SUCCESS;
        }

        if (fileList.size() > 3) {
            Aligner alignIt(fileList, 0.11, !allowMirror, !useCPD, !useRandom);
            alignIt.superImpose();
            return 1;
        } else if (remap){

            if (!fileExists(ref) ){ // check file exists
                SASTOOLS_UTILS_H::logger("CANNOT READ FILE", ref);
                throw std::invalid_argument("** ERROR => file not adequate " + ref);
            }

            Aligner alignIt(ref);



        } else {
            // do pair, requires reference and target to be defined
            // validate files then

            if (!fileExists(ref) ){ // check file exists
                SASTOOLS_UTILS_H::logger("CANNOT READ FILE", ref);
                throw std::invalid_argument("** ERROR => file not adequate " + ref);
            }

            if (!fileExists(tar) ){ // check file exists
                SASTOOLS_UTILS_H::logger("CANNOT READ FILE", tar);
                throw std::invalid_argument("** ERROR => file not adequate " + tar);
            }

            Aligner(ref, tar, useBackbone, !allowMirror, !useCPD, !useRandom);
        }

    } catch (boost::program_options::required_option& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (boost::program_options::error& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        return ERROR_IN_COMMAND_LINE;
    } catch (const std::invalid_argument& e){
        std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
        std::cerr<<"Type "<<typeid(e).name()<<std::endl;
    }

    return 0;
}
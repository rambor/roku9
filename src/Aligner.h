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
#ifndef ROKU9_ALIGNER_H
#define ROKU9_ALIGNER_H
//
#include "Eigen/Dense"
#include "PointCloud.h"
#include <cpd/matrix.hpp>
#include <cpd/rigid.hpp>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <utility>
#include <boost/regex.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem/operations.hpp>
#include <utility>
#include <random>

struct ScoreSet {
    /// keep the score
    unsigned int ref;
    unsigned int target;
    bool mirror;
    float score;
    float nsd;

    bool operator<(const ScoreSet& a) const {
        return score < a.score;
    }
};


class Aligner {
    std::string filename;
    std::vector<std::string> pdbfiles;
    std::vector<PointCloud> pointClouds;
    float subselection = 0.17;
    float subselectionUpper = 0.37;
    float edge_radius = 0.0f;
    unsigned int totalTrials = 13157;
    unsigned int unsuccessfulCount = 2157;
    unsigned int theta_limit = 37;
    unsigned int phi_limit = 37;
    unsigned int psi_limit = 37;

    unsigned int totalFiles;
    float bead_radius;
    bool tryBothEnantiomorphs=true;
    bool useCPD = true;
    bool useRandom = false;


    bool checkFile();
    void extractCoordinates();
    unsigned long factorial(unsigned int n_value);

    void cpdAlgorithm(PointCloud &refModel, PointCloud &tarModel);

    //void heapAlgorithm(int maxindex);
    void heapAlgorithm(std::vector<unsigned int> & c_array, std::vector<unsigned int> & permutationIndices, unsigned int * i_index);
    bool nSelectK(std::vector<unsigned int> & permutationIndices,  unsigned int n_index);


    void printVector(const vector3 & vector);

public:
    Aligner();
    Aligner(std::string reffilename, std::string tarfilename, bool isAtomistic, bool tryBoth, bool useCPD, bool useRandom);
    Aligner(std::string filename, float subselection, bool tryBoth, bool useCPD, bool useRandom);
    Aligner(std::string filename);


    std::string getFileExt(const std::string &s);

    void superImpose();

    float calculateVarianceScore(cpd::Matrix &moved, PointCloud &refModel, PointCloud &tarModel);
    float calculateVarianceScore(Eigen::MatrixXf &moved, PointCloud &refModel, PointCloud &tarModel);
    float getCorrespondenceScore(Eigen::MatrixXf &moved, PointCloud &refModel, PointCloud &tarModel);



    void fillMatrixUsingPermutationIndices(std::vector<unsigned int> &permutationIndices, PointCloud &pointCloud,
                                           Eigen::MatrixXf &matrix);

    Eigen::MatrixXf getRotation(Eigen::MatrixXf &moving, Eigen::MatrixXf &fixed);


    float testAssignment(Eigen::MatrixXf &fixed, Eigen::MatrixXf &moving, PointCloud &refModel, PointCloud &tarModel);

    unsigned long nChoosek(unsigned int n_index, unsigned int k_index);


    float calculateAlignmentErrorCVX(PointCloud &refModel, PointCloud &tarModel);


    void
    fillTargetMartrixUsingSelectedSet(std::set<unsigned int> &workingSet, PointCloud &pointCloud, Eigen::MatrixXf &matrix);

    void createReMappedFilteredSet(std::vector<vector3> & vectorOfCoordinates);

    void fillCrossValidationMatrix(Eigen::MatrixXf &matrix, PointCloud &tarModel);

    void ICPCutoff(PointCloud &refModel, PointCloud &tarModel);

    Eigen::MatrixXf assignedCorrespondencesWithCutoff(float cutoff, PointCloud &refModel, PointCloud &tarModel, bool subsampling);

    float unNormalizedHausdorfScore(Eigen::MatrixXf &transformed, PointCloud &refModel, PointCloud &tarModel);
    float unNormalizedHausdorfScore(cpd::Matrix &transformed, PointCloud &refModel);

    Eigen::MatrixXf randomRotation(Eigen::MatrixXf matrix);
    Eigen::MatrixXf getRandomRotation();

    float NSDScore(cpd::Matrix &transformed, cpd::Matrix &fixed, float dminTar, float dminRef);

    void writeLogFile(std::string name, std::vector<ScoreSet> & scores, std::vector<ScoreSet> & paired);
    void writeKDEFileList(std::string basename, std::vector<std::string> & files);


    void align(ScoreSet &scored, PointCloud &referenceModel, PointCloud &tarModel, bool updateModel);

    Eigen::MatrixXf generateRotationMatrix(float theta, float phi, float psi);

    const PointCloud & getPointCloudAt(unsigned int index) const { return pointClouds[index];}

    // requires sorted vector
    inline float getMedian(std::vector<float> &values){
        std::sort(values.begin(), values.end());
        auto total = (unsigned int)values.size();
        if ((total & 1) == 0){ // even
            unsigned int mid = total/2;
            return (values[mid] + values[mid+1])*0.5f;
        } else {
            return values[(total+1)/2 - 1];
        }
    }

    void AddToCoordinateVector(std::vector<vector3> & fillme, const vector3 * pVecs, unsigned int total);
    void writeCoordinateVectorToFile(std::vector<vector3> & vectorOfCoordinates);

    float getEpsilon(int total);
};


#endif //ROKU9_ALIGNER_H

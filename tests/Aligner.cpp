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

#include <support.hpp>
#include "gtest/gtest.h"
#include "../src/Aligner.h"
#include "../src/PointCloud.h"

class AlignerTest : public ::testing::Test {


public:
    PointCloud base9_ref;
    PointCloud base9_tar;

    AlignerTest() : ::testing::Test(),
                    base9_ref(PointCloud(tests::fixture("rfd1_9.pdb"), false, 1.0)),
                    base9_tar(PointCloud(tests::fixture("rfd1_9r.pdb"), false, 1.0)) {
    }

};

// check the calculation of the hausdorfscore
TEST_F(AlignerTest, checkHausdorfScore) {
    //Aligner(std::string filename, float subselection, bool tryBoth, bool useCPD, bool useRandom);
    Aligner align_direct = Aligner(tests::loadInp("aligner.inp"), 0.11, false, false, false);

    auto refModel = align_direct.getPointCloudAt(0);
    auto tarModel = align_direct.getPointCloudAt(1);

    Eigen::MatrixXf untransformedTarget(tarModel.getTotalPointsInCloud(), 3);

    auto coordinates = tarModel.getCenteredCoordinates();
    for(size_t i=0; i<tarModel.getTotalPointsInCloud(); ++i){
        const vector3 * vec = &(coordinates[i]);
        untransformedTarget(i, 0) = vec->x;
        untransformedTarget(i, 1) = vec->y;
        untransformedTarget(i, 2) = vec->z;
    }

    float tmpscore = align_direct.unNormalizedHausdorfScore(untransformedTarget, refModel, tarModel);

    std::cout << " Hausdorf score unaligned :: " << tmpscore << std::endl;
    EXPECT_GT(tmpscore, 1);

    // perform calculation with perfect correspondance
    tmpscore = align_direct.unNormalizedHausdorfScore(untransformedTarget, base9_tar, tarModel);
    EXPECT_FLOAT_EQ(tmpscore, 0.0f);
}


// check the calculation of the hausdorfscore
TEST_F(AlignerTest, checkCalculateVarianceScore) {

    Aligner align_direct = Aligner(tests::loadInp("aligner.inp"), 0.11, false, false, false);

    auto refModel = align_direct.getPointCloudAt(0);
    auto tarModel = align_direct.getPointCloudAt(1); //rfd1_9r

    cpd::Matrix untransformedTarget(tarModel.getTotalPointsInCloud(), 3);

    auto coordinates = tarModel.getCenteredCoordinates();
    for(size_t i=0; i<tarModel.getTotalPointsInCloud(); ++i){
        const vector3 * vec = &(coordinates[i]);
        untransformedTarget(i, 0) = (double)vec->x;
        untransformedTarget(i, 1) = (double)vec->y;
        untransformedTarget(i, 2) = (double)vec->z;
    }

    float tmpscore = align_direct.calculateVarianceScore(untransformedTarget, refModel, tarModel);
    std::cout << " SVD score aligned 1 :: " << tmpscore << std::endl;
    EXPECT_LT(tmpscore, 0.9) << "CVX hull points should not match due to floating point error " << tmpscore;

    // perform calculation with perfect correspondance
    tmpscore = align_direct.calculateVarianceScore(untransformedTarget, base9_tar, tarModel);
    std::cout << " SVD score aligned 2 :: " << tmpscore << std::endl;
    EXPECT_NEAR(tmpscore, 0.0f, FLT_EPSILON);
}

TEST_F(AlignerTest, checkConstructor){
    Aligner align_direct = Aligner(tests::loadInp("aligner.inp"), 0.11, false, false, false);


}
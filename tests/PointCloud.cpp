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
#include "../src/PointCloud.h"


class PointCloudTest : public ::testing::Test {

public:
    PointCloud base9_n;
    PointCloud base9_t;

    PointCloudTest() : ::testing::Test(),
                       base9_n(PointCloud(tests::fixture("rfd1_9.pdb"), false, 1.0)),
                       base9_t(PointCloud(tests::fixture("rfd1_9r.pdb"), false, 1.0)) {
    }
};


TEST_F(PointCloudTest, IsTheSameCoordinates) {
    EXPECT_EQ(base9_n.getTotalPointsInCloud(),base9_t.getTotalPointsInCloud());
}

TEST_F(PointCloudTest, checkClosestDistance) {
    EXPECT_NEAR(base9_n.getClosestDistance(), base9_t.getClosestDistance(), 0.0015);
    EXPECT_NEAR(base9_n.getClosestDistance(), 2*4.350, 0.01); // twice bead radius
}

//tests centered coordinates are not same as original
TEST_F(PointCloudTest, checkCoordinates) {

    unsigned int total = base9_n.getTotalPointsInCloud();

    auto coordsN = base9_n.getCoordinates();
    auto coordsC = base9_n.getCenteredCoordinates();

    for(unsigned int n=0; n<total; n++){
        const vector3 &pR = *(coordsN + n);
        const vector3 &pC = *(coordsC + n);

        EXPECT_NE(pR.x, pC.x);
        EXPECT_NE(pR.y, pC.y);
        EXPECT_NE(pR.z, pC.z);
    }


    for(unsigned int n=0; n<total; n++){
        const vector3 &pR = *(coordsN + n);
        const vector3 &pC = *(coordsC + n);

        EXPECT_NEAR(pR.x, pC.x + base9_n.getCentering_vec().x, 0.001);
        EXPECT_NEAR(pR.y, pC.y + base9_n.getCentering_vec().y, 0.001);
        EXPECT_NEAR(pR.z, pC.z + base9_n.getCentering_vec().z, 0.001);
    }
}
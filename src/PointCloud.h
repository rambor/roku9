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

#ifndef ROKU9_POINTCLOUD_H
#define ROKU9_POINTCLOUD_H

#include "Eigen/Dense"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <utility>
#include <sastools/include/vector3.h>
#include <sastools/include/PDBModel.h>
#include <random>
#include <boost/regex.hpp>
#include <cpd/matrix.hpp>

#ifdef __cplusplus
extern "C" {
#endif
#include <libqhull/qhull_a.h>
#ifdef __cplusplus
}
#endif


class PointCloud {

    std::string filename;
    bool isAtomistic = false, isMirror = false,  useKept=false;
    float bead_radius=0.0f, cvx_volume=0.0f, surface_area=0.0f, edge_radius=0.0f;

    std::vector<vector3> coordinates;
    std::set<unsigned int> cvx_hull;
    std::set<unsigned int> non_cvx_hull;
    std::set<unsigned int> cvx_hull_subset;
    std::set<unsigned int> cross_validation_set;
    std::set<unsigned int> all_indices;
    Eigen::MatrixXf rotation;

    float closestDistance= FLT_MAX;
    float avgDistance = 0.0f;
    float smax;


    vector3 centering_vec, translation_vec;
    std::vector<vector3> centered_coordinates;
    std::vector<vector3> transformed_coordinates;
    std::vector<unsigned int> shell_assignments;

    unsigned int totalCoordinates;

    float score=0.0f, nsd=0.0f;

    void extractCoordinates();

public:

    PointCloud(std::string pdbfile, bool isAtomistic, float bead_radius);

    PointCloud(const PointCloud & cloud2);

    // move assignment operator
    PointCloud & operator=(const PointCloud & cloudToMove){
        if (&cloudToMove == this){
            return *this;
        }

        PointCloud tmp(cloudToMove);
        tmp.swap(*this);
        return *this;
    }

    // move constructor
    PointCloud(PointCloud && cloudToMove) noexcept  {
        cloudToMove.swap(*this);
    }

    ~PointCloud() = default;

    void swap(PointCloud & other) noexcept {
        using std::swap;
        other.filename = filename;
        other.isAtomistic = isAtomistic;
        other.isMirror = isMirror;
        other.useKept = useKept;
        other.bead_radius = bead_radius;
        other.cvx_volume = cvx_volume;
        other.surface_area = surface_area;
        other.avgDistance = avgDistance;
        other.smax = smax;
        other.closestDistance = closestDistance;
        other.score = score;
        other.nsd = nsd;
        other.totalCoordinates = totalCoordinates;
        other.edge_radius = edge_radius;

        std::swap(cvx_hull, other.cvx_hull);
        std::swap(non_cvx_hull, other.non_cvx_hull);
        std::swap(cvx_hull_subset, other.cvx_hull_subset);
        std::swap(cross_validation_set, other.cross_validation_set);
        std::swap(all_indices, other.all_indices);

        std::swap(coordinates, other.coordinates);
        std::swap(centered_coordinates, other.centered_coordinates);
        std::swap(transformed_coordinates, other.transformed_coordinates);
        std::swap(shell_assignments, other.shell_assignments);

        other.centering_vec = vector3(centering_vec);
        other.translation_vec = vector3(translation_vec);
        other.rotation = rotation;
    }

    size_t getTotalPointsInCloud(){ return coordinates.size(); }

    const vector3 * getCenteredCoordinates() const { return centered_coordinates.data();}

    std::vector<vector3> & getCenteredCoordinatesVector() { return centered_coordinates;}

    const vector3 * getCoordinates() const { return coordinates.data();}

    const vector3 * getTransformedCoordinates() const { return transformed_coordinates.data();}

    void setTransformedCoordinatesEigen(Eigen::MatrixXf & matrix);
    void setTransformedCoordinates(cpd::Matrix& matrix);

    const std::string &getFilename() const;

    const std::set<unsigned int> & getCvx_hull() const;
    const std::set<unsigned int> & getCvx_hull_Subset() const;
    const std::set<unsigned int> & getNon_Cvx_hull_Subset() const;
    const std::set<unsigned int> & getCrossValidationSet() const;
    const std::set<unsigned int> & getAllIndices() const;

    void setComplexHull();

    float getClosestDistance() const;
    float getAvgDistance() const;

    bool getKept(){ return useKept;}
    void setKept(bool val){ useKept = val;}
    void setRandomCVXPoints();

    void setScore(float value){ score = value; }
    void setNSD(float value){ nsd = value; }
    float getScore(){ return score;}
    float getNSD(){ return nsd;}

    std::string writeToFile();
    void  writeToFile(std::string name);
    std::string writeCenteredCoordinates(std::string name);
    std::string writeCenteredCoordinatesAtomistic(std::string name);
    std::string writeOriginalCoordinates();

    void writeToCSVFile();

    void setCrossValidationSet(float percent);

    void writeRandomSubsetTofile(float percent);

    void makeMirror();

    void resetEnantiomer();

    void replaceTransformedCoordinates(const vector3 * newVectors);

    void translateTransformedCoordinates(const vector3 &vec);

    const vector3 &getCentering_vec() const;
    void calculateClosestDistance();

    const Eigen::MatrixXf & getRotationMatrix() const { return rotation;}

    void setRotationMatrix(cpd::Matrix & rotation);
    void setRotationMatrix(Eigen::MatrixXf & rotation);

    void updateRotationMatrix(Eigen::MatrixXf & rotation);

    float getEdgeRadius() { return edge_radius;}

    float getSmax() const;

    float assignPointsToShell(float shell_thickness);
};


#endif //ROKU9_POINTCLOUD_H

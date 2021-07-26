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
//
//   Created by Robert P Rambo on 13/07/2018.
//
#include "PointCloud.h"


PointCloud::PointCloud(std::string pdbfile, bool isAtomistic, float br) : filename(std::move(pdbfile)), isAtomistic(isAtomistic), bead_radius(br) {
    /*
     *  assume pdbfile was already validated
     */

    this->extractCoordinates();
    this->calculateClosestDistance();
    this->setRandomCVXPoints();
    this->setCrossValidationSet(0.17);
}

/*
 * copy constructor
 */
PointCloud::PointCloud(const PointCloud & cloud2) : filename(cloud2.filename), isAtomistic(cloud2.isAtomistic), bead_radius(cloud2.bead_radius) {

    isMirror = false;

    for(unsigned int i=0; i<cloud2.totalCoordinates; i++){
        coordinates.emplace_back(vector3(cloud2.coordinates[i]));
        centered_coordinates.emplace_back(vector3(cloud2.centered_coordinates[i]));
    }

    if (cloud2.transformed_coordinates.size() > 0){
        this->transformed_coordinates.resize(cloud2.transformed_coordinates.size());
        for(unsigned int i=0; i<cloud2.totalCoordinates; i++){
            this->transformed_coordinates[i] = vector3(cloud2.transformed_coordinates[i]);
        }
    }

    centering_vec = vector3(cloud2.centering_vec);
    translation_vec = vector3(cloud2.translation_vec);

    closestDistance = cloud2.closestDistance;
    rotation = cloud2.rotation;

    totalCoordinates = cloud2.totalCoordinates;
    surface_area = cloud2.surface_area;
    cvx_volume = cloud2.cvx_volume;
    edge_radius = cloud2.edge_radius;
    avgDistance = cloud2.avgDistance;
    score = cloud2.score;
    nsd = cloud2.nsd;

    cvx_hull = std::set<unsigned int> (cloud2.cvx_hull);
    non_cvx_hull = std::set<unsigned int> (cloud2.non_cvx_hull);
    cvx_hull_subset = std::set<unsigned int> (cloud2.cvx_hull_subset);
    cross_validation_set = std::set<unsigned int> (cloud2.cross_validation_set);
    all_indices = std::set<unsigned int> (cloud2.all_indices);;

    isMirror = cloud2.isMirror;
    useKept=cloud2.useKept;

}


/*
 * Assume first file in list is the reference
 */
void PointCloud::extractCoordinates(){

    PDBModel baseModel(filename, false, true);

    smax = baseModel.getSMax();

    isAtomistic = (baseModel.getTotalUniqueAtoms() == 1 && isAtomistic) ? false : isAtomistic;

    unsigned int totalCoords = baseModel.getTotalCoordinates();

    const float * pX = baseModel.getX();
    const float * pY = baseModel.getY();
    const float * pZ = baseModel.getZ();

    if (isAtomistic){ // only take backbone
        for(unsigned int i = 0; i<totalCoords; i++){
            if (baseModel.belongsToResidue(i) && baseModel.isBackbone(i)){
                coordinates.emplace_back(vector3(pX[i],pY[i],pZ[i]));
                centered_coordinates.emplace_back(vector3(baseModel.getCenteredXVec()[i], baseModel.getCenteredYVec()[i], baseModel.getCenteredZVec()[i]));
            }
        }
    } else {
        for(unsigned int i = 0; i<totalCoords; i++){
            coordinates.emplace_back(vector3(pX[i],pY[i],pZ[i]));
            centered_coordinates.emplace_back(vector3(baseModel.getCenteredXVec()[i], baseModel.getCenteredYVec()[i], baseModel.getCenteredZVec()[i]));
        }
    }

    centering_vec = vector3(*baseModel.getCenteringVector());

    SASTOOLS_UTILS_H::logger("Total Coordinates", std::to_string(coordinates.size()));

    totalCoordinates = coordinates.size();
    transformed_coordinates.resize(totalCoordinates);
    rotation = Eigen::MatrixXf(3,3);

    // assuming points are not overlapping?  Each point belongs in spherical shell
    float testdis, edgesum=0.0f;
    float count = 0.0f;
    const auto * pVecs = coordinates.data();
    for(unsigned int i=0; i<(totalCoords-1); i++){
        const auto * pVec1 = &pVecs[i];
        float edge = FLT_MAX;
        for(unsigned int j=(i+1); j<totalCoords; j++){
            testdis = (*pVec1 - (pVecs[j])).length();
            if (testdis < edge){
                edge = testdis;
            }
        }
        edgesum += edge;
        count += 1.0f;
    }
    edge_radius = edgesum/count*0.5f;
    setComplexHull();
    logger("","FINISHED EXTRACTING COORDINATES");
}

float PointCloud::getSmax() const {
    return smax;
}


/*
 * for each point in set, calculate closest neighboring distance
 */
void PointCloud::calculateClosestDistance(){

    float avgDis=0.0f;
    float dis;

    vector3 * const pCoords = coordinates.data();
    float tempClosest = FLT_MAX;

    for (unsigned int n=0; n < coordinates.size(); n++) {

        const vector3 &pRef = pCoords[n];

        for(unsigned int j=n+1; j<coordinates.size(); j++){
            dis = (pCoords[j] - pRef).length();
            avgDis += dis;

            if (dis < tempClosest){
                tempClosest = dis;
            }
        }
        // calculate average closest distance
    }
    closestDistance = tempClosest;
    avgDistance = avgDis/(float)(coordinates.size()*(coordinates.size()-1)/2);
}



const std::string &PointCloud::getFilename() const {
    return filename;
}

void PointCloud::setComplexHull(){
    logger("","SETTING CVX HULL");
    char flags[] = "qhull FA";
    unsigned int numpoints = 3*totalCoordinates;
    coordT points[numpoints];

    const vector3 * ptr = coordinates.data();

    std::vector<unsigned int> active_indices(totalCoordinates);

    for (unsigned int i=0; i<totalCoordinates; i++){
        const vector3 & vec = ptr[i];
        points[3*i] = vec.x;
        points[3*i+1] = vec.y;
        points[3*i+2] = vec.z;
        active_indices[i] = i;
        all_indices.insert(i);
    }

    qh_new_qhull (3, totalCoordinates, points, 0, flags, nullptr, nullptr);

    vertexT * vertices = qh vertex_list;
    auto totalV = (unsigned int)qh num_vertices;

//    for(unsigned int i = workingLimit; i < num; i++) {
//        unsigned int value = ptr[i];
//        beadToPoint(testPoint, pModel->getBead(value));
//        // exclude HULL points, for each bead, determine if outside HULL
//        qh_findbestfacet (testPoint, qh_ALL, &bestdist, &isoutside);
//
//        if (!isoutside){
//            hull.insert(value);
//        }
//    }

    // only move CVX hull points

    cvx_hull.clear();
    for (unsigned int v = 0; v < totalV; v++) { //
        cvx_hull.insert( active_indices[qh_pointid( vertices->point)]);
        vertices = vertices->next;
    }

    for(unsigned int i=0; i<totalCoordinates; i++){
        if (cvx_hull.find(i) == cvx_hull.end()){
            non_cvx_hull.insert(i);
        }
    }

    surface_area = (float)(qh totarea);
    cvx_volume = (float)(qh totvol);
    qh_freeqhull(true);
logger("","FINISHED CVX HULL");
}

/**
 * pick 5 random points on the Convex Hull
 */
void PointCloud::setRandomCVXPoints(){

    // calculate all pairwise distances, exclude the pair that constitute longest distance
    cvx_hull_subset.clear();
    std::vector<unsigned int> indices(cvx_hull.size());
    std::copy(cvx_hull.begin(), cvx_hull.end(), indices.begin());

    float dis = FLT_MIN;
    unsigned int first=0, second=1;

    // find the longest distance pair
    for(unsigned int i=0; i<cvx_hull.size(); i++){
        const vector3 & vec1 = coordinates[indices[i]];
        unsigned int next = i+1;
        for(unsigned int j=next; j<cvx_hull.size(); j++){
            float temp = (vec1 - coordinates[indices[j]]).length();
            if (temp > dis){
                dis = temp;
                first = i;
                second = j;
            }
        }
    }

    auto limit  = (unsigned int)(indices.size()-2);
    auto iter = std::find(indices.begin(), indices.end(), indices[first]);

    if (iter == indices.end()){
        std::cout << " EXITING " << std::endl;
        throw std::invalid_argument("** ERROR => CORRUPTED INPUT FILE : NON SENSE CVX HULL");
    }

    std::iter_swap(iter, indices.begin() + (indices.size()-1));

    iter = std::find(indices.begin(), indices.end(), indices[second]);
    std::iter_swap(iter, indices.begin() + limit);

    //std::sort(indices.begin(), indices.begin() + limit);
    // pick random point, find the furthest point
    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(indices.begin(), indices.begin() + limit, gen);

    dis = FLT_MIN;

    // store first point
    cvx_hull_subset.insert(indices[0]);
    const vector3 & vec1 = coordinates[indices[0]];

    limit--;
    std::iter_swap(indices.begin() , indices.begin() + limit);

    unsigned int keepit = 0;
    for(unsigned int j=0; j<limit; j++){
        float tmp = (vec1 - coordinates[indices[j]]).length();
        if (tmp > dis){
            dis = tmp;
            keepit = j;
        }
    }

    // store second point
    cvx_hull_subset.insert(indices[keepit]);
    const vector3 & vec2 = coordinates[indices[keepit]];

    limit--;
    std::iter_swap(indices.begin() + keepit, indices.begin() + limit);

    // find next point that is furthest from both by area
    dis = FLT_MIN;
    for(unsigned int j=0; j<limit; j++){
        const vector3 & vec3 = coordinates[indices[j]];
        float tmp = 0.5f*((vec1 - vec3).length() + (vec2 - vec3).length() + (vec1 - vec2).length());
        if (tmp > dis){
            dis = tmp;
            keepit = j;
        }
    }
    // store third point
    cvx_hull_subset.insert(indices[keepit]);
    const vector3 & vec3 = coordinates[indices[keepit]];
    // move to limit so not in use
    limit--;
    std::iter_swap(indices.begin() + keepit, indices.begin() + limit);

    // find fourth point that maximizes area of triangular pyramid
    dis = FLT_MIN;
    for(unsigned int j=0; j<limit; j++){
        const vector3 & vec = coordinates[indices[j]];
        float tmp1 = 0.5f*((vec1 - vec).length() + (vec2 - vec).length() + (vec1 - vec2).length());
        float tmp2 = 0.5f*((vec3 - vec).length() + (vec2 - vec).length() + (vec3 - vec2).length());
        float tmp3 = 0.5f*((vec1 - vec).length() + (vec3 - vec).length() + (vec1 - vec3).length());
        float tmp = tmp1 + tmp2 + tmp3;
        if (tmp > dis){
            dis = tmp;
            keepit = j;
        }
    }

    // store fourth point
    cvx_hull_subset.insert(indices[keepit]);
    const vector3 & vec4 = coordinates[indices[keepit]];
    limit--;
    std::iter_swap(indices.begin() + keepit, indices.begin() + limit);

    // find last point, largest area and furthest from vec4
    dis = FLT_MIN;
    for(unsigned int j=0; j<limit; j++){
        const vector3 & vec = coordinates[indices[j]];
        float tmp1 = 0.5f*((vec1 - vec).length() + (vec2 - vec).length() + (vec1 - vec2).length());
        float tmp2 = 0.5f*((vec3 - vec).length() + (vec2 - vec).length() + (vec3 - vec2).length());
        float tmp3 = 0.5f*((vec1 - vec).length() + (vec3 - vec).length() + (vec1 - vec3).length());
        float tmp = tmp1 + tmp2 + tmp3 + (vec - vec4).length();
        if (tmp > dis){
            dis = tmp;
            keepit = j;
        }
    }

    cvx_hull_subset.insert(indices[keepit]);

    // add subset to non
    for(auto ind : cvx_hull_subset){
        non_cvx_hull.insert(ind);
    }

}


void PointCloud::setCrossValidationSet(float percent){

    std::vector<unsigned int> indices(this->getTotalPointsInCloud());
    unsigned int * const ptr = (this->getTotalPointsInCloud() != 0) ? indices.data() : nullptr;
    for(unsigned int i = 0; i < this->getTotalPointsInCloud(); i++) {
        ptr[i] = i;
    }

    // partition indices by moving select cvx points to back
    auto limit  = (unsigned int)(this->getTotalPointsInCloud()-1);
    for(auto ind : cvx_hull_subset){
        auto iter = std::find(indices.begin(), indices.end(), ind);
        if (iter != indices.end()){
            std::iter_swap(iter, indices.begin() + limit);
            limit--;
        } else { // throw error

        }
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    std::shuffle(indices.begin(), indices.begin() + limit, gen);

    auto grabLimit = (unsigned int)(this->getTotalPointsInCloud()*percent);

    cross_validation_set = std::set<unsigned int>(cvx_hull_subset);
    for(unsigned int i=0; i<grabLimit; i++){
        cross_validation_set.insert(indices[i]);
    }
}


const std::set<unsigned int> & PointCloud::getCvx_hull() const {
    return cvx_hull;
}

const std::set<unsigned int> & PointCloud::getCvx_hull_Subset() const {
    return cvx_hull_subset;
}

const std::set<unsigned int> & PointCloud::getNon_Cvx_hull_Subset() const {
    return non_cvx_hull;
}

const std::set<unsigned int> & PointCloud::getAllIndices() const {
    return all_indices;
}

float PointCloud::getClosestDistance() const {
    return closestDistance;
}

float PointCloud::getAvgDistance() const {
    return avgDistance;
}

/**
 * matrix is [3 x N], N is the number of cols
 * need this as default to keep right handed rotation
 * @param matrix
 */
void PointCloud::setTransformedCoordinatesEigen(Eigen::MatrixXf &matrix) {

//    transformedMatrix = std::move(matrix);
    for(long i=0; i < matrix.cols(); i++){
//        vector3 & vec = transformed_coordinates[i];
        auto colInUse = matrix.col(i);
        transformed_coordinates[i] = vector3(colInUse(0), colInUse(1), colInUse(2));
//        vec.x = colInUse(0);
//        vec.y = colInUse(1);
//        vec.z = colInUse(2);

//        vec.x = matrix(0,i);
//        vec.y = matrix(1,i);
//        vec.z = matrix(2,i);
//        transformed_coordinates.emplace_back(vector3(matrix(i,0), matrix(i,1), matrix(i,2)));
    }
}

/*
 * matrix is double precision
 */
void PointCloud::setTransformedCoordinates(cpd::Matrix &matrix) {
    transformed_coordinates.clear();
    for (cpd::Matrix::Index row = 0; row < matrix.rows(); ++row) {
        transformed_coordinates.emplace_back(vector3((float)matrix(row,0), (float)matrix(row,1), (float)matrix(row,2)));
    }
}


/*
 * matrix is double precision
 */
void PointCloud::translateTransformedCoordinates(const vector3 & vec) {
    translation_vec = vector3(vec);

    for(auto & tvec : transformed_coordinates){
        tvec += vec;
    }
}

void PointCloud::replaceTransformedCoordinates(const vector3 * newVectors) {
    transformed_coordinates.clear();
    for(unsigned int i=0; i<totalCoordinates; i++){
        transformed_coordinates.emplace_back(vector3(newVectors[i]));
    }

}

std::string PointCloud::writeCenteredCoordinatesAtomistic(std::string subname) {

    FILE * pFile;

    boost::filesystem::path p1(filename);
    std::string basename = p1.filename().c_str();
    std::string name = basename.substr(0, basename.size()-4) + "_" + subname +".pdb";


    pFile = fopen(name.c_str(), "w");

    fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", this->filename.c_str());
    for (unsigned int n=0; n < coordinates.size(); n++) {
        const vector3 & vec = centered_coordinates[n];
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, "CA ", "ALA", "A", (n+1), vec.x, vec.y, vec.z );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
    return name;
}


std::string PointCloud::writeCenteredCoordinates(std::string subname) {

    FILE * pFile;

    boost::filesystem::path p1(filename);
    std::string basename = p1.filename().c_str();
    std::string name = basename.substr(0, basename.size()-4) + "_" + subname +".pdb";

    pFile = fopen(name.c_str(), "w");

    fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", this->filename.c_str());
    for (unsigned int n=0; n < coordinates.size(); n++) {
        const vector3 & vec = centered_coordinates[n];
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, "CA ", "ALA", "A", (n+1), vec.x, vec.y, vec.z );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
    return name;
}

std::string PointCloud::writeOriginalCoordinates() {

    FILE * pFile;

    boost::filesystem::path p1(filename);
    std::string basename = p1.filename().c_str();
    std::string name = basename.substr(0, basename.size()-4) + "_ref.pdb";

    pFile = fopen(name.c_str(), "w");

    fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", this->filename.c_str());
    for (unsigned int n=0; n < coordinates.size(); n++) {
        const vector3 & vec = coordinates[n];
        fprintf(pFile, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, "CA ", "ALA", "A", (n+1), vec.x, vec.y, vec.z );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
    return name;
}


std::string PointCloud::writeToFile() {

    FILE * pFile;
    std::string ending = (isMirror) ? "_m_r.pdb" : "_r.pdb";

    boost::filesystem::path p1(filename);

    SASTOOLS_UTILS_H::logger("FILENAME", p1.filename().c_str());
    SASTOOLS_UTILS_H::logger("PARENT PATH", p1.parent_path().c_str());
    SASTOOLS_UTILS_H::logger("ROOT PATH", p1.root_path().c_str());

    filename = p1.filename().c_str();

    std::string name = filename.substr(0, filename.size()-4) + ending;
    SASTOOLS_UTILS_H::logger("WRITING FILE", name);

    pFile = fopen(name.c_str(), "w");

    fprintf(pFile,"REMARK                       : ROTATION MATRIX\n");
    fprintf(pFile,"REMARK                 ROW 1 : %.3f %.3f %.3f\n", rotation(0,0), rotation(0,1), rotation(0,2));
    fprintf(pFile,"REMARK                 ROW 2 : %.3f %.3f %.3f\n", rotation(1,0), rotation(1,1), rotation(1,2));
    fprintf(pFile,"REMARK                 ROW 3 : %.3f %.3f %.3f\n", rotation(2,0), rotation(2,1), rotation(2,2));
    fprintf(pFile,"REMARK                       : \n");
    fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", this->filename.c_str());
    fprintf(pFile,"REMARK      CENTERING VECTOR : %.3f %.3f %.3f\n", centering_vec.x, centering_vec.y, centering_vec.z);
    fprintf(pFile,"REMARK                       : \n");
    fprintf(pFile,"REMARK                       : \n");
    fprintf(pFile,"REMARK                       : \n");
    fprintf(pFile,"REMARK                       : \n");


    if (isAtomistic){ // only take backbone
        PDBModel baseModel(filename, false, true);

        unsigned int totalCoords = baseModel.getTotalCoordinates();

        const float * pX = baseModel.getCenteredXVec();
        const float * pY = baseModel.getCenteredYVec();
        const float * pZ = baseModel.getCenteredZVec();
        Eigen::MatrixXf coorVec(3,1);
        Eigen::MatrixXf transVec(3,1);
        transVec(0,0) = translation_vec.x;
        transVec(1,0) = translation_vec.y;
        transVec(2,0) = translation_vec.z;

        for(unsigned int i = 0; i<totalCoords; i++){
            // transform coordinates and write to file
            coorVec(0,0) = (isMirror) ? -pX[i] : pX[i];
            coorVec(1,0) = pY[i];
            coorVec(2,0) = pZ[i];
            // center coordinates
            // translate
            Eigen::MatrixXf rotated = rotation*coorVec + transVec;
            unsigned int indexn = i+1;
            unsigned int resid = baseModel.getResid()[i];
            std::string resi = baseModel.getResidueAt(i);

            fprintf(pFile, "ATOM  %5i %.4s %.3s A%4i    %8.3f%8.3f%8.3f  1.00100.00\n", indexn, baseModel.getAtomTypeByIndex(i).c_str(), resi.c_str(), resid, rotated(0,0), rotated(1,0), rotated(2,0) );
        }

    } else {

        for (unsigned int n=0; n < transformed_coordinates.size(); n++) {
            vector3 & vec = transformed_coordinates[n];
            auto vecx = (double)vec.x;
            auto vecy = (double)vec.y;
            auto vecz = (double)vec.z;
            unsigned int indexn = n+1;
            fprintf(pFile, "ATOM  %5i  CA  ALA A%4i    %8.3f%8.3f%8.3f  1.00100.00\n", indexn, indexn, vecx, vecy, vecz );
        }
    }


    fprintf(pFile,"END\n");
    fclose(pFile);
    return name;
}

void PointCloud::writeToFile(std::string newname) {

    FILE * pFile;
    std::string ending = newname + ".pdb";

    boost::filesystem::path p1(filename);

//    std::cout << " -> " << p1.filename().c_str() << std::endl;
//    std::cout << " -> " << p1.parent_path() << std::endl;
//    std::cout << " -> " << p1.root_path() << std::endl;

    filename = p1.filename().c_str();

    std::cout << " Writing file :: " << newname << std::endl;

    pFile = fopen(ending.c_str(), "w");

    fprintf(pFile,"REMARK                       : ROTATION MATRIX\n");
    fprintf(pFile,"REMARK                 ROW 1 : %.3f %.3f %.3f\n", rotation(0,0), rotation(0,1), rotation(0,2));
    fprintf(pFile,"REMARK                 ROW 2 : %.3f %.3f %.3f\n", rotation(1,0), rotation(1,1), rotation(1,2));
    fprintf(pFile,"REMARK                 ROW 3 : %.3f %.3f %.3f\n", rotation(2,0), rotation(2,1), rotation(2,2));
    fprintf(pFile,"REMARK                       : \n");
    for (size_t n=0; n < transformed_coordinates.size(); n++) {
        vector3 & vec = transformed_coordinates[n];
        double vecx = (double)vec.x;
        double vecy = (double)vec.y;
        double vecz = (double)vec.z;
        unsigned int indexn = n+1;
        fprintf(pFile, "ATOM  %5i  CA  ALA A%4i    %8.3f%8.3f%8.3f  1.00100.00\n", indexn, indexn, vecx, vecy, vecz );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}


void PointCloud::writeRandomSubsetTofile(float percent){
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> randomDis(0.0, 1);
    float theta, phi, z, sqrtz, sint, cost;
    auto limit = (unsigned int)(percent*coordinates.size());
    Eigen::MatrixXf reflection_matrix = Eigen::MatrixXf::Identity(3, 3);
    Eigen::MatrixXf rot(3,3);

    for (int r=0; r<10; r++){

        std::shuffle(coordinates.begin(), coordinates.end(), gen);
        // generate random matrix
        theta = (float)(2.0*randomDis(gen)*M_PI);
        phi = (float)(2.0*M_PI*randomDis(gen));
        z = (float)(2.0*randomDis(gen));

        sqrtz = sqrt(z);
        sint = sin(theta);
        cost = cos(theta);
        Eigen::Vector3f vector(sin(phi)*sqrtz, cos(phi)*sqrtz, sqrt(2.0f-z));

        rot(0,0) = cost;
        rot(0,1) = sint;
        rot(0,2) = 0.0f;
        rot(1,0) = -sint;
        rot(1,1) = cost;
        rot(1,2) = 0;
        rot(2,0) = 0;
        rot(2,1) = 0;
        rot(2,2) = 1;

        Eigen::MatrixXf matrix = (vector*vector.transpose()-reflection_matrix)*rot;

        FILE * pFile;
        std::string ending = "_"+std::to_string(r)+".pdb";
        std::string name = filename.substr(0, filename.size()-4) + ending;

        pFile = fopen(name.c_str(), "w");

        fprintf(pFile,"REMARK  CENTERED COORDINATES : %s\n", this->filename.c_str());

        for (unsigned int n=0; n < limit; n++) {
            const vector3 & vec = coordinates[n];
            Eigen::Vector3f coorvec (vec.x, vec.y, vec.z);
            Eigen::Matrix3Xf rotvec = matrix*coorvec;
            fprintf(pFile, "%-6s%5i %4s %3s %1s%4i    %8.3f%8.3f%8.3f  1.00100.00\n", "ATOM", n+1, "CA ", "ALA", "A", n, rotvec(0), rotvec(1), rotvec(2) );
        }
        fprintf(pFile,"END\n");
        fclose(pFile);
    }

}

void PointCloud::writeToCSVFile() {

    FILE * pFile;
    std::string ending = (isMirror) ? ".csv" : ".csv";
    std::string name = filename.substr(0, filename.size()-4) + ending;

    pFile = fopen(name.c_str(), "w");

    for (unsigned int n=0; n < coordinates.size(); n++) {
        const vector3 & vec = coordinates[n];
        fprintf(pFile, "%.3f,%.3f,%.3f\n", vec.x, vec.y, vec.z );
    }
    fclose(pFile);
}

const std::set<unsigned int> &PointCloud::getCrossValidationSet() const {
    return cross_validation_set;
}


void PointCloud::resetEnantiomer(){

    for(unsigned int i = 0; i<this->totalCoordinates; i++){
        coordinates[i].x *= -1.0f;
        centered_coordinates[i].x *= -1.0f;
    }

    if(transformed_coordinates.size() > 0){
        for(unsigned int i = 0; i<this->totalCoordinates; i++){
            transformed_coordinates[i].x *= -1.0f;
        }
    }

    centering_vec.x *= -1.0f;
    isMirror = false;
    useKept = false;
}


void PointCloud::makeMirror(){

    for(unsigned int i = 0; i<this->totalCoordinates; i++){
        coordinates[i].x *= -1.0f;
        centered_coordinates[i].x *= -1.0f;
    }

    if(!transformed_coordinates.empty()){
        for(unsigned int i = 0; i<this->totalCoordinates; i++){
            transformed_coordinates[i].x *= -1.0f;
        }
    }

    centering_vec.x *= -1.0f;
    isMirror = true;
}

const vector3 &PointCloud::getCentering_vec() const {
    return centering_vec;
}

void PointCloud::setRotationMatrix(cpd::Matrix & rot) {

    for(unsigned int i=0; i<3; i++){
        rotation.row(i) << rot(i,0), rot(i,1), rot(i,2);
    }
}

void PointCloud::setRotationMatrix(Eigen::MatrixXf & rot) {

    for(unsigned int i=0; i<3; i++){
        rotation.row(i) = rot.row(i);
    }
}


void PointCloud::updateRotationMatrix(Eigen::MatrixXf & rot) {
    rotation = rot*rotation;
}

float PointCloud::assignPointsToShell(float shell_thickness) {
    unsigned int index = 0;
    for (auto & vec : centered_coordinates){
        shell_assignments[index] = (unsigned int)std::floor(vec.length()/shell_thickness);
        index++;
    }

    return 0;
}

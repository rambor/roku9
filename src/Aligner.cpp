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
//   Created by Robert P Rambo on 24/01/2019.
//

#include "Aligner.h"

Aligner::Aligner(){

}
/*
 * need a list of files to align
 *
 * for each file, use as reference and align both images of target
 *
 * keep best match evaluated using Hausdorf metric
 *
 */
Aligner::Aligner(std::string filname, float subselection, bool tryBoth, bool useCPD, bool useRandom) :
        filename(std::move(filname)),
        subselection(subselection),
        tryBothEnantiomorphs(tryBoth),
        useCPD(useCPD),
        useRandom(useRandom){

    // validate each file in the list
    std::string ext = this->getFileExt(filename);
    std::cout << "** CHECKING FILE " << ext << std::endl;

    /*
     * compile list of coordinates into single array of pointClouds
     * validate file with checkFile and that each PDB file is available
     */
    if (ext == "inp"){
        if (this->checkFile()){
            // extract coordinates
            this->extractCoordinates();
        } else {
            throw std::invalid_argument("** ERROR => CORRUPTED INPUT FILE : " + filename);
        }
    } else {
        throw std::invalid_argument("** ERROR => INCORRECT INPUT FILE MUST be .inp : " + filename);
    }
}


Aligner::Aligner(std::string filename) : filename(std::move(filename)){

    pointClouds.emplace_back(PointCloud(this->filename, false, 1.0f));

    this->createReMappedFilteredSet(pointClouds[0].getCenteredCoordinatesVector() );

}

/*
 *
 */
Aligner::Aligner(std::string reffilename, std::string tarfilename, bool isAtomistic, bool tryBoth, bool useCPD, bool useRandom) :
        filename(std::move(reffilename)),
        tryBothEnantiomorphs(tryBoth),
        useCPD(useCPD),
        useRandom(useRandom) {

    PointCloud referenceModel(filename, isAtomistic, 1.0f);
    PointCloud target(tarfilename, isAtomistic, 1.0f);

    float maxSmax = std::max(referenceModel.getSmax(), target.getSmax());
    float delta = 3.4;

    referenceModel.writeCenteredCoordinates("centered");

    ScoreSet score ({0,1,false,0.0f,0.0f});

    align(score, referenceModel, target, true); // transformed coordinates are saved to target

    target.setKept(true);
    // translate model to reference position
    std::cout << " => TRANSLATING TO REFERENCE" << std::endl;
    target.translateTransformedCoordinates(referenceModel.getCentering_vec());

    std::cout << " => FINISHED ALIGNING, WRITING FILE" << std::endl;
    target.writeToFile();
}

/*
 * should check of list of files that exists
 *
 */
bool Aligner::checkFile() {

    bool returnMe = false;
    // read in file
    std::ifstream data (filename, std::ifstream::in);

    if (data.is_open()) {

        boost::regex format("pdb");
        std::string line;

        while(!data.eof()) //
        {
            getline(data, line); //this function grabs a line and moves to next line in file
            /*
             * require at least two columns (1: q, 2: I(q), 3: sigma)
             */
            if (isspace(line[0])){
                line.erase(line.begin(), std::find_if(line.begin(), line.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            }

            std::vector<std::string> tempLine;
            boost::split(tempLine, line, boost::is_any_of("\t  "), boost::token_compress_on);

            if(boost::regex_search(tempLine[0], format)){
                SASTOOLS_UTILS_H::logger("CHECKING PDB COORDINATE FILE", tempLine[0]);

                if (!boost::filesystem::exists(tempLine[0])){
                    throw std::invalid_argument("** ERROR => FILE NOT FOUND : " + tempLine[0]);
                }

                FileClass fileclass(tempLine[0]);
                //std::cout << ":: " << fileclass.getFullPath() << std::endl;
                pdbfiles.push_back(fileclass.getFullPath());

            } else {
                std::cout << "NO FILE FOUND FOR LINE  =>  " << line <<std::endl;
            }
        }
        returnMe = true;
    }

    totalFiles = (unsigned int)pdbfiles.size(); // totalfiles to align
    data.close();
    return returnMe;
}


std::string Aligner::getFileExt(const std::string& s) {

    size_t i = s.rfind('.', s.length());
    if (i != std::string::npos) {
        return(s.substr(i+1, s.length() - i));
    }

    return("");
}

/*
 * for the validated pdb coordinate files, create PointCloud objects
 */
void Aligner::extractCoordinates(){

    for(const auto & file : pdbfiles){
        pointClouds.emplace_back(PointCloud(file, false, 1.0f));
    }
}


void Aligner::superImpose(){

    // random
//    std::random_device rd;
//    std::mt19937 gen(rd());
    //std::shuffle(pointClouds.begin(), pointClouds.end(), gen);
    /*
     * find the best aligned pair
     */
//    PointCloud & refModel = pointClouds[0];
//    refModel.writeRandomSubsetTofile(0.13);
//    exit(0);

    std::vector<ScoreSet> scores;

    // create cross correlation table
    Eigen::MatrixXf ccTable(totalFiles, totalFiles);

    unsigned permutation = 0;
    unsigned int totalPermutations = (totalFiles *(totalFiles-1))/2;
    for(unsigned int i=0; i<totalFiles; i++){ // all pairwise comparisons

        PointCloud & refModel = pointClouds[i];
        refModel.setKept(false);

        for(unsigned int j=i+1; j<totalFiles; j++){
            std::cout << "  ALIGNMENT ROUND => " << (permutation + 1) << " of " << totalPermutations << std::endl;
            PointCloud & tarModel = pointClouds[j];

            scores.emplace_back( ScoreSet {i,j,false,0.0f,0.0f} );
            align(scores[permutation], refModel, tarModel, false);

            //reset tarModel
//            if (scores[permutation].mirror == true){
//                //tarModel.resetEnantiomer();
//                pointClouds[i] = PointCloud(tarModel.getFilename(), false, 1.0f);
//            }
            //tarModel.setKept(false);
            permutation++;
        }
    }

    // find the reference model that provides lowest variance
    float invN = 1.0f/(float)(totalFiles-1);
    float median=0, bestVar = FLT_MAX;
    unsigned int keptRefIndex=0;

    std::vector<float>tempScores(totalFiles-1);

    for(unsigned int i=0; i<totalFiles; i++){
        // down the column
        float xsquared=0.0f, sumx=0.0f;

        for(unsigned int row=0; row<i; row++){
            int index = row*totalFiles - (row*(row+1)/2) + (i - row) - 1;
            float value = scores[index].score;
            tempScores[row] = value;
            xsquared += value*value;
            sumx += value;
        }

        //across the row
        for(unsigned int col=(i+1); col<totalFiles; col++){
            int index = i*totalFiles - i*(i+1)/2 - i + col - 1;
            float value = scores[index].score;
            tempScores[col-1] = value;
            xsquared += value*value;
            sumx += value;
        }

        float tempvar = invN*xsquared - invN*invN*sumx*sumx;

        if (tempvar < bestVar){
            std::printf("FOUND BETTER REFERENCE %s :: %.3E < %.3E \n", pointClouds[i].getFilename().c_str(), tempvar, bestVar);
            bestVar = tempvar;
            keptRefIndex = i;
            median = getMedian(tempScores); // calculate median absolute deviation
        }
    }

    // if I reject a model, do median absolute deviation?
    for(unsigned int row=0; row<keptRefIndex; row++){
        int index = row*totalFiles - (row*(row+1)/2) + (keptRefIndex - row) - 1;
        float value = scores[index].score;
        tempScores[row] = std::abs(value-median);
    }

    //across the row
    for(unsigned int col=(keptRefIndex+1); col<totalFiles; col++){
        int index = keptRefIndex*totalFiles - keptRefIndex*(keptRefIndex+1)/2 - keptRefIndex + col - 1;
        float value = scores[index].score;
        tempScores[col-1] = std::abs(value-median);
    }

    float mad = 1.4826f*getMedian(tempScores);

//    std::sort(scores.begin(), scores.end());
//    for(auto & scr : scores){
//        std::cout << scr.ref << " :: " << scr.target << " " << scr.mirror << " " << scr.score << " " << scr.nsd <<  std::endl;
//    }

    /*
     * assemble the models
     * first element in scores is the reference
     * first pair in scores should be entry in pairedList
     */

    std::vector<std::string> filesList;
    std::vector<ScoreSet> pairedList;

    scores.clear();
    logger("SETTING COMMON REFERENCE", pointClouds[keptRefIndex].getFilename());

    PointCloud * referenceModel = &(pointClouds[keptRefIndex]);
    referenceModel->setKept(false); // insures we used centered coordinates
    filesList.push_back(referenceModel->writeOriginalCoordinates());

    std::vector<vector3> finalcoordinates;
    AddToCoordinateVector(finalcoordinates, referenceModel->getCoordinates(), referenceModel->getTotalPointsInCloud());

    // all models should be aligned to reference model in its original, non-centered position
    edge_radius=referenceModel->getEdgeRadius();
    float counter = 1.0f;

    for(unsigned int i=0; i<keptRefIndex; i++){
        PointCloud target(pointClouds[i].getFilename(), false, 1.0f);

        logger("TARGET", target.getFilename());
        scores.emplace_back( ScoreSet {keptRefIndex,i,false,0.0f,0.0f} );

        ScoreSet & scored = scores[scores.size()-1];
        align(scored, *referenceModel, target, true);

        target.setKept(true);
        // translate model to reference position
        logger("", "TRANSLATING TO REFERENCE");
        target.translateTransformedCoordinates(referenceModel->getCentering_vec());

        if ( (scored.score-median)/mad < 2.5){ // keep
            AddToCoordinateVector(finalcoordinates, target.getTransformedCoordinates(), target.getTotalPointsInCloud());
            edge_radius += target.getEdgeRadius();
            counter += 1.0f;
        }

        filesList.push_back(target.writeToFile());

        pairedList.emplace_back(ScoreSet {scored.ref,scored.target,scored.mirror,scored.score,scored.nsd});
    }

    for(unsigned int i=(keptRefIndex+1); i<pointClouds.size(); i++){
        PointCloud target(pointClouds[i].getFilename(), false, 1.0f);

        logger("TARGET", target.getFilename());
        scores.emplace_back( ScoreSet {keptRefIndex,i,false,0.0f,0.0f} );

        ScoreSet & scored = scores[scores.size()-1];
        align(scored, *referenceModel, target, true);

        // translate model to reference position
        logger("", "TRANSLATING TO REFERENCE");
        target.translateTransformedCoordinates(referenceModel->getCentering_vec());

        if ( (scored.score-median)/mad < 2.5){ // keep
            AddToCoordinateVector(finalcoordinates, target.getTransformedCoordinates(), target.getTotalPointsInCloud());
            edge_radius += target.getEdgeRadius();
            counter += 1.0f;
        }

        filesList.push_back(target.writeToFile());

        pairedList.emplace_back(ScoreSet {scored.ref,scored.target,scored.mirror,scored.score,scored.nsd});
    }

    // write consolidated list to file

    writeCoordinateVectorToFile(finalcoordinates);
    writeLogFile("roku9.log", scores, pairedList);
    writeKDEFileList("list", filesList);
}

void Aligner::createReMappedFilteredSet(std::vector<vector3> & vectorOfCoordinates){
    DminType dmin_values = SASTOOLS_UTILS_H::getDminValues(vectorOfCoordinates);

    float dmin_supremum = dmin_values.dmin_supremum;
//    float dmin_infimum = dmin_values.dmin_infimum;
//    float average_dmin = dmin_values.average_dmin;
//    float stdev_dmin = dmin_values.stdev_dmin;

    logger("DMIN SUPREMUM", formatNumber(dmin_values.dmin_supremum, 2));
    logger("DMIN INFIMUM", formatNumber(dmin_values.dmin_infimum, 2));
    logger("DMIN AVERAGE", formatNumber(dmin_values.average_dmin, 2));
    logger("DMIN STDEV", formatNumber(dmin_values.stdev_dmin, 2));

    //formula fwhm = 2.355*sigma which should account for 76% area in distribution
    // dmin supremum should correspond to actual grid spacing used in ab initio modeling

    float grid_spacing = std::max(dmin_supremum*0.666666f, dmin_values.stdev_dmin*2.355f);
    //float grid_spacing = dmin_values.average_dmin;

    // make grid based on average_dmin
    const vector3 * pVec = vectorOfCoordinates.data();

    int total = vectorOfCoordinates.size();

    float minX= pVec->x;
    float maxX= pVec->x;
    float minY= pVec->y;
    float maxY= pVec->y;
    float minZ= pVec->z;
    float maxZ= pVec->z;

    for(int i = 1; i<total; i++){
        const vector3 &pVec1 = *(pVec + i);
        if (pVec1.x < minX){
            minX = pVec1.x;
        }
        if (pVec1.x > maxX){
            maxX = pVec1.x;
        }
        if (pVec1.y < minY){
            minY = pVec1.y;
        }
        if (pVec1.y > maxY){
            maxY = pVec1.y;
        }
        if (pVec1.z < minZ){
            minZ = pVec1.z;
        }
        if (pVec1.z > maxZ){
            maxZ = pVec1.z;
        }
    }

    vector3 x_max_y_min_z_min(maxX, minY, minZ);
    vector3 x_max_y_max_z_min(maxX, maxY, minZ);
    vector3 x_min_y_max_z_min(minX, maxY, minZ);
    vector3 x_max_y_max_z_max(maxX, maxY, maxZ);

    float length_x = 2*(maxX - minX) + grid_spacing;
    float length_y = 2*(maxY - minY) + grid_spacing;
    float length_z = 2*(maxZ - minZ) + grid_spacing;

    // make hexagonal grid
    logger("=> CREATING GRID","BOXED SEARCH SPACE");
    float br = grid_spacing*0.5;
    float inv_bead_radius = 1.0f/br;

    const auto x_prime_index = (int) std::ceil(length_x*inv_bead_radius*0.5);
    const auto z_prime_index = (int) std::ceil(length_z*inv_bead_radius*0.5);
    const auto y_prime_index = (int) std::ceil(length_y*inv_bead_radius*0.5);

    const vector3 a1 = vector3(2*br, 0, 0);
    const vector3 a2 = vector3(-br, br*sqrt(3.0f), 0);
    const vector3 a3 = vector3(0, 0, (float)(4.0f*sqrt(6)/3.0f*br)); // seems fine, spreads out

    // basis 1 is (0,0,0)
    const vector3 basis2 = a1*(2.0f/3.0f) + a2*(1.0f/3.0f) + a3/2.0f;

    float xlimit = length_x*0.5f;
    float ylimit = length_y*0.5f;
    float zlimit = length_z*0.5f;

    std::vector<vector3> points;

    for (int n3 = -z_prime_index; n3<z_prime_index; n3++){
        for (int n2 = -y_prime_index; n2 <= y_prime_index; n2++) {
            for (int n1 = -x_prime_index; n1 <= x_prime_index; n1++) {
                const vector3 temp = (a1 * n1 + a2 * n2 + a3 * n3);
                if (std::abs(temp.x) < xlimit && std::abs(temp.y) < ylimit && std::abs(temp.z) < zlimit){
                    points.emplace_back(vector3(temp.x, temp.y, temp.z));
                }

                vector3 temp2 = temp + basis2;
                if (std::abs(temp2.x) < xlimit && std::abs(temp2.y) < ylimit && std::abs(temp2.z) < zlimit) {
                    points.emplace_back(vector3(temp2.x, temp2.y, temp2.z));
                }
            }
        }
    }



    std::string residue_index;

    std::string name = "remapped_unfiltered.pdb";
    FILE * pFile;
    pFile = fopen(name.c_str(), "w");
    const unsigned int totalAtoms = vectorOfCoordinates.size();
    fprintf(pFile,"REMARK  ACCEPTED COORDINATES : %s\n", name.c_str());
    fprintf(pFile,"REMARK 265                   EDGE RADIUS : %.3f\n", grid_spacing);

    std::vector<bool> useIt(totalAtoms);
    for(unsigned int i=0; i<totalAtoms; i++){
        useIt[i] = true;
    }

    const char *colour[6] = { "A", "B", "C", "D", "E", "F" };

    int indexer = 1;
    int ind = 0;
    for(auto & pt : points){
        for (unsigned int n=0; n < totalAtoms; n++) {
            if (useIt[n]){
                pVec = &vectorOfCoordinates[n];
                if ((pt - *pVec).length() <= dmin_values.average_dmin){ // don't want the output model to be inflated by the upscale of the grid
                    residue_index = boost::lexical_cast<std::string>(indexer+1);
                    fprintf(pFile, "ATOM  %5i  CA  ALA %s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", indexer, colour[ind], residue_index.c_str(), pt.x, pt.y, pt.z);
                    useIt[n] = false;
                    indexer++;
                    break;
                }
            }
        }

        if (indexer%9999 == 0){
            indexer = 1;
            ind += 1;
        }
    }
    fprintf(pFile,"END\n");
    fclose(pFile);


    name = "centered_ref.pdb";
    pFile = fopen(name.c_str(), "w");
    fprintf(pFile,"REMARK  ACCEPTED COORDINATES : %s\n", name.c_str());
    fprintf(pFile,"REMARK 265                   EDGE RADIUS : %.3f\n", grid_spacing);

    indexer = 1;
    ind = 0;

    for (unsigned int n=0; n < totalAtoms; n++) {
        pVec = &vectorOfCoordinates[n];

        residue_index = boost::lexical_cast<std::string>(indexer+1);
        fprintf(pFile, "ATOM  %5i  CA  ALA %s%4s    %8.3f%8.3f%8.3f  1.00100.00\n", indexer, colour[ind], residue_index.c_str(), pVec->x, pVec->y, pVec->z);
        indexer++;

        if (indexer%9999 == 0){
            indexer = 0;
            ind += 1;
        }
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}

void Aligner::ICPCutoff(PointCloud & refModel, PointCloud & tarModel){

    Eigen::MatrixXf untransformedTarget(3, tarModel.getTotalPointsInCloud());

    auto coordinates = tarModel.getCenteredCoordinates();
    for(size_t i=0; i<tarModel.getTotalPointsInCloud(); ++i){
        const vector3 * vec = &(coordinates[i]);
        untransformedTarget(0,i) = vec->x;
        untransformedTarget(1,i) = vec->y;
        untransformedTarget(2,i) = vec->z;
    }


    Eigen::MatrixXf untransformedTest(3, 2);

    untransformedTest(0,0) = 0;
    untransformedTest(1,0) = 0.5;
    untransformedTest(2,0) = 0.9;
    untransformedTest(0,1) = 0;
    untransformedTest(1,1) = -0.5;
    untransformedTest(2,1) = -0.9;


    // correspondences are only assigned using a cut-off (< 1 Angstrom?)
//    tarModel.setCrossValidationSet(0.17);
//    Eigen::MatrixXf cross_validation_matrix (tarModel.getCrossValidationSet().size(), 3);
//    this->fillCrossValidationMatrix(cross_validation_matrix, tarModel);
    /*
     *
     * start with a random rotations
     * to test varying correspondences
     *
     */
    float score = FLT_MAX;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "          SEARCHING INITIAL POSES       " << std::endl;
    std::cout << "----------------------------------------" << std::endl;

    std::clock_t startTime = std::clock();
    //std::clock_t profileTime = std::clock();

    double runtimeavg=0.0d;
    Eigen::MatrixXf keepIt, keptrotation;
    Eigen::MatrixXf rotation(3,3);

    int successIt=0;
    if (useRandom){
        std::cout << "        SEARCH :: SPHERICAL GRID"  << std::endl;
        std::cout << "  total trials :: " << (theta_limit*phi_limit*psi_limit) << std::endl;

        int total_sphere_grid_points = 100;
        std::cout << "  total trials :: " << (total_sphere_grid_points*psi_limit) << std::endl;
        // grid search
        // use Fibonacci grid

        auto inv_golden_ratio = (float)((2.0f*M_PI)/((1.0f + sqrtf(5.0))*0.5f));
        float epsilon = getEpsilon(total_sphere_grid_points);
        auto delta_psi = (float)(2.0*M_PI/(double)psi_limit);
        int totalinGrid = 0;

        for(int j=0; j < total_sphere_grid_points; j++){

            double angle_theta = j*inv_golden_ratio; // rotation around X
            double angle_phi = acos(1.0d - 2.0d*(j+epsilon)/(total_sphere_grid_points - 1 + 2.0d*epsilon)); // rotation around z

            for (unsigned int psi_n=0; psi_n<psi_limit; psi_n++){ // range is from 0 to 360
                // want to rotate about the vector

                Eigen::MatrixXf rotationMatrix = generateRotationMatrix(angle_theta, angle_phi, psi_n*delta_psi);
                Eigen::MatrixXf test_pose = rotationMatrix*untransformedTarget; // [3 x N] matrix

                Eigen::MatrixXf test_pose2 = rotationMatrix*untransformedTest; // [3 x N] matrix

                //std::cout << test_pose2(0,1) << " " << test_pose2(1,1)<< " " << test_pose2(2,1) << std::endl;

                tarModel.setTransformedCoordinatesEigen(test_pose); // generate a pose
                tarModel.setKept(true);

                // generate rotation matrix from rotated pose
                Eigen::MatrixXf transformedPAll = assignedCorrespondencesWithCutoff((tarModel.getClosestDistance()*1.1189990f), refModel, tarModel, false)*test_pose;

                float tmpscore = calculateVarianceScore(transformedPAll, refModel, tarModel);

                if (tmpscore < score){
                    //std::printf("=> RND (%6i) :: %.4f < %.4E\n", totalinGrid, tmpscore, score );

                    score = tmpscore;
                    keepIt = std::move(transformedPAll); // keep the best pose for starting refinement
                    keptrotation = std::move(rotationMatrix);
                    successIt++;
                }
                runtimeavg += (std::clock() - startTime)/(double) CLOCKS_PER_SEC;

                totalinGrid++;
            }

        }
//        exit(0);

//        auto delta_theta = (float)(2.0*M_PI/(double)theta_limit);
//        auto delta_phi = (float)(2.0*M_PI/(double)phi_limit);
//        auto delta_psi = (float)(2.0*M_PI/(double)psi_limit);
//
//        int totalinGrid = 0;
//        for (unsigned int theta_n=0; theta_n<theta_limit; theta_n++){ // range is from 0 to 180
//            float angle_theta = theta_n*delta_theta;
//            for (unsigned int phi_n=0; phi_n<phi_limit; phi_n++){ // range is from 0 to 360
//                float angle_phi = phi_n*delta_phi;
//                for (unsigned int psi_n=0; psi_n<psi_limit; psi_n++){ // range is from 0 to 360
//
//
//                    Eigen::MatrixXf rotationMatrix = generateRotationMatrix(angle_theta, angle_phi, psi_n*delta_psi);
//                    Eigen::MatrixXf test_pose = rotationMatrix*untransformedTarget; // [3 x N] matrix
//
//                    tarModel.setTransformedCoordinatesEigen(test_pose); // generate a pose
//                    tarModel.setKept(true);
//
////                    profileTime = std::clock();
//
//                    // generate rotation matrix from rotated pose
//                    Eigen::MatrixXf transformedPAll = assignedCorrespondencesWithCutoff((tarModel.getClosestDistance()*1.1189990f), refModel, tarModel, false)*test_pose;
//
////                    std::cout << " assigning " << ((std::clock() - profileTime)/(double) CLOCKS_PER_SEC) << std::endl;
////                    profileTime = std::clock();
//
//                    float tmpscore = calculateVarianceScore(transformedPAll, refModel, tarModel);
////                    std::cout << " variance score " << ((std::clock() - profileTime)/(double) CLOCKS_PER_SEC) << std::endl;
//
//                    if (tmpscore < score){
//                        std::printf("=> RND (%6i) :: %.4f < %.4E\n", totalinGrid, tmpscore, score );
////                        tarModel.writeToFile("model_"+std::to_string(successIt));
//                        score = tmpscore;
//                        keepIt = std::move(transformedPAll); // keep the best pose for starting refinement
//                        keptrotation = std::move(rotationMatrix);
//                        successIt++;
//                    }
//                    runtimeavg += (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
//
//                    totalinGrid++;
//                }
//            }
//        }
//        tarModel.writeCenteredCoordinates("unrefined");
    } else { // fast mode
        std::cout << "        SEARCH :: RANDOM ROTATION SAMPLING"  << std::endl;
        std::cout << "  total trials :: " << totalTrials << std::endl;
        Eigen::MatrixXf randomRotMatrix(3,3);
        Eigen::MatrixXf randomlyRotated(3, untransformedTarget.cols()), transformedPAll(3, untransformedTarget.cols());

        for(unsigned int i=0; i< totalTrials; i++){
            /*
             * steps:
             * 1. randomly rotate movingPWhole centered at center-of-mass
             * 2. then perform SVD ICP rotation to minimize difference to reference
             * 3. score it, keep best transformed set
             * lack of memory, each trial is fully independent of prior
             */
            startTime = std::clock();
            //Eigen::MatrixXf randomlyRotated = randomRotation(untransformedTarget);
            randomRotMatrix = getRandomRotation();
            randomlyRotated = randomRotMatrix*untransformedTarget;

            tarModel.setTransformedCoordinatesEigen(randomlyRotated);
            tarModel.setKept(true);

            rotation = assignedCorrespondencesWithCutoff((tarModel.getClosestDistance()*0.610f), refModel, tarModel, false);
            transformedPAll = rotation*randomlyRotated;

            float tmpscore = calculateVarianceScore(transformedPAll, refModel, tarModel);
            if (tmpscore < score){
                printf("=> RND (%5i) :: %.4f < %.4E\n", i, tmpscore, score );
                score = tmpscore;
                keepIt = std::move(transformedPAll); //keep the best pose for starting refinement
                keptrotation = rotation*randomRotMatrix;
            }
            runtimeavg += (std::clock() - startTime)/(double) CLOCKS_PER_SEC;
        }
    }

    std::cout << "   AVG  runtime :: " << runtimeavg/(double)totalTrials << std::endl;

    float tmpscore;
    tarModel.setRotationMatrix(keptrotation);
    tarModel.setTransformedCoordinatesEigen(keepIt);
    tarModel.setKept(true);

    std::cout << "----------------------------------------" << std::endl;
    std::cout << "          REFINING SELECTED POSE        " << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    printf(" STARTING SCORE :: %.4f\n", score);
    // try a bunch of random poses, lack of memory between each try, but keep the best
    unsigned int count = 0, counter=0;
    Eigen::MatrixXf tempTransformedAll;

    while (count < unsuccessfulCount){ // should really scale with complexity of model
        // probability of grabbing correct correspondance
        rotation = assignedCorrespondencesWithCutoff((tarModel.getClosestDistance()*0.735f), refModel, tarModel, true);

        tempTransformedAll = rotation*keepIt; // should be [3 x N]

        tmpscore = calculateVarianceScore(tempTransformedAll, refModel, tarModel);
        if (tmpscore < score){
            printf("=> RFN (%5i) :: %.4f < %.4f\n", counter, tmpscore, score );
            score = tmpscore;
            keepIt = std::move(tempTransformedAll); // should be [N x 3] matrix
            tarModel.setTransformedCoordinatesEigen(keepIt);

            keptrotation = std::move(rotation);
            tarModel.updateRotationMatrix(keptrotation);

            count = 0;
            successIt++;
        }
        count++;
        counter++;
    }

    logger("iteration", formatNumber(counter));
    logger("reference", refModel.getFilename());
    logger("target", tarModel.getFilename());

    if (count == unsuccessfulCount){
        logger("MAX ITERATIONS REACHED", formatNumber(unsuccessfulCount));
    }

    logger("best score", formatNumber(score, 4));

    //tarModel.setRotationMatrix(keptrotation);
    //tarModel.updateRotationMatrix(keptrotation);
    // refine using all available points within cutoff
    rotation = assignedCorrespondencesWithCutoff((tarModel.getClosestDistance()*0.637f), refModel, tarModel, false);
    tempTransformedAll = rotation*keepIt;

    tmpscore = calculateVarianceScore(tempTransformedAll, refModel, tarModel);

    logger("","FINAL ASSIGNMENT");
    logger("FINAL", formatNumber(tmpscore, 5));
    logger("FINAL", formatNumber(score, 5));

    if (tmpscore < score){
        logger("ACCEPTED", formatNumber(tmpscore, 5));
        score = tmpscore;
        keepIt = std::move(tempTransformedAll);
//        tarModel.setRotationMatrix(rotation);
        tarModel.updateRotationMatrix(rotation);
    }

    tarModel.setScore(score);
    tarModel.setNSD(unNormalizedHausdorfScore(keepIt, refModel, tarModel));
    tarModel.setTransformedCoordinatesEigen(keepIt);
}



Eigen::MatrixXf Aligner::getRotation(Eigen::MatrixXf & moving, Eigen::MatrixXf & fixed){
    Eigen::MatrixXf cov = moving.transpose() * fixed;
    Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::NoQRPreconditioner> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXf s_matrix = svd.singularValues().asDiagonal();
    Eigen::MatrixXf reflection_matrix = Eigen::MatrixXf::Identity(3, 3);
    float d_value = (svd.matrixV() * svd.matrixU().transpose()).determinant();
    reflection_matrix(2, 2) = d_value;
    return (svd.matrixU() * reflection_matrix * svd.matrixV().transpose());
}




/*
 * Should I only score CVX points?
 *
 */
float Aligner::testAssignment(Eigen::MatrixXf & fixed, Eigen::MatrixXf & moving, PointCloud & refModel, PointCloud & tarModel){

    Eigen::MatrixXf rotation = getRotation(moving, fixed);

    Eigen::MatrixXf movingPWhole(tarModel.getTotalPointsInCloud(), 3);

    auto coordinates = tarModel.getCenteredCoordinates();
    for(size_t i=0; i<tarModel.getTotalPointsInCloud(); ++i){
        const vector3 * vec = &(coordinates[i]);
        movingPWhole(i, 0) = vec->x;
        movingPWhole(i, 1) = vec->y;
        movingPWhole(i, 2) = vec->z;
    }

    movingPWhole = movingPWhole*rotation;

    return std::abs(getCorrespondenceScore(movingPWhole, refModel, tarModel));
}




void Aligner::fillMatrixUsingPermutationIndices(std::vector<unsigned int> &permutationIndices, PointCloud &pointCloud,
                                                Eigen::MatrixXf &matrix){

    auto cvx_indices = pointCloud.getCvx_hull();
    std::vector<unsigned int> points(cvx_indices.begin(), cvx_indices.end());

    int rowCount = 0;
    auto coordinates = pointCloud.getKept() ? pointCloud.getTransformedCoordinates() : pointCloud.getCenteredCoordinates();
    for(auto & ind : permutationIndices){
        const vector3 * vec = &(coordinates[points[ind - 1]]);
        matrix(rowCount,0) = vec->x;
        matrix(rowCount,1) = vec->y;
        matrix(rowCount,2) = vec->z;
        rowCount++;
    }
}


/**
 * Uses Coherent Point Drift algorithm
 * will update transformed coordinates in target Model (tarModel)
 * @param refModel
 * @param tarModel
 */
void Aligner::cpdAlgorithm(PointCloud &refModel, PointCloud &tarModel){
    /*
     * both matrices are double precision
     */
    cpd::Matrix fixed(refModel.getTotalPointsInCloud(), 3);
    cpd::Matrix moving(tarModel.getTotalPointsInCloud(), 3);

    const vector3 * const coordinates = refModel.getKept() ? refModel.getTransformedCoordinates() : refModel.getCenteredCoordinates();
    const vector3 * const tarCoordinates = tarModel.getKept() ? tarModel.getTransformedCoordinates() : tarModel.getCenteredCoordinates();

    std::cout << "              :: ASSEMBLING MATRICES " << refModel.getTotalPointsInCloud() << " " << tarModel.getTotalPointsInCloud() << std::endl;
    std::cout << "    REF MODEL :: " << refModel.getTotalPointsInCloud() << std::endl;
    std::cout << "    TAR MODEL :: " << tarModel.getTotalPointsInCloud() << std::endl;
    std::cout << "              :: " << std::endl;

    for(size_t i=0; i < refModel.getTotalPointsInCloud(); ++i){
        const vector3 &vec = coordinates[i];
        fixed(i, 0) = (double)vec.x;
        fixed(i, 1) = (double)vec.y;
        fixed(i, 2) = (double)vec.z;
    }

    for(size_t i=0; i < tarModel.getTotalPointsInCloud(); ++i){
        const vector3 & vec = tarCoordinates[i];
        moving(i, 0) = (double)vec.x;
        moving(i, 1) = (double)vec.y;
        moving(i, 2) = (double)vec.z;
    }

    std::cout << "              :: STARTING CPD Algorithm" << std::endl;

    cpd::Rigid rigid;
    rigid.outliers(0.2);
//    rigid.tolerance(1e-5);
    rigid.reflections(false);
    rigid.scale(false);

    std::clock_t startTime = std::clock();

    cpd::RigidResult result = rigid.run(fixed, moving);

    cpd::Matrix transformed = result.points; // transformed with respect to center of fixed

    std::cout << "    iteration :: " << result.iterations << std::endl;
    std::cout << "    reference :: " << refModel.getFilename() << std::endl;
    std::cout << "       target :: " << tarModel.getFilename() << std::endl;
    std::cout << "      runtime :: " << (std::clock() - startTime)/(double) CLOCKS_PER_SEC << std::endl;

//    float nsdscore = NSDScore(transformed, fixed, tarModel.getClosestDistance(), refModel.getClosestDistance());
    float nsdscore = unNormalizedHausdorfScore(transformed, refModel);
    float initial_score = calculateVarianceScore(transformed, refModel, tarModel);
//    float initial_score = result.sigma2;
    std::cout << "       scores :: " << initial_score << " ( NSD " << nsdscore << " ) " << std::endl;

    /*
     * At this point, we have
     */
    tarModel.setScore(initial_score);
    tarModel.setNSD(nsdscore);
    tarModel.setTransformedCoordinates(transformed);
    tarModel.setRotationMatrix(result.rotation);
}

/*
 * return rotation matrix
 */
Eigen::MatrixXf Aligner::assignedCorrespondencesWithCutoff(float cutoff, PointCloud & refModel, PointCloud & tarModel, bool subsampling){

    const vector3 * refcoordinates = refModel.getKept() ? refModel.getTransformedCoordinates() : refModel.getCenteredCoordinates();
    const vector3 * tarCoordinates = tarModel.getKept() ? tarModel.getTransformedCoordinates() : tarModel.getCenteredCoordinates();

    auto totalCoords = (unsigned int)tarModel.getTotalPointsInCloud();
    auto totalRefCoords = (unsigned int)refModel.getTotalPointsInCloud();
    float temp, dis;
    unsigned int keep;
    bool accepted;

    std::map<unsigned int, unsigned int> tar_to_ref;

    unsigned int count = 0;
    Eigen::MatrixXf cov;

    if (subsampling){
        // correspondences should only be calculated off selected indices to speed up
        std::random_device rd;
        std::mt19937 gen(rd());
        /*
         * for selected indices - find corresponding index in reference
         *  grab random set of coordinates
         */
        std::uniform_real_distribution<> randomDis(subselection,0.57);
        std::uniform_int_distribution<unsigned int> randomInt(0, (totalCoords-1));

        auto percent = (float)randomDis(gen);
        auto subTotal = (unsigned int)((totalCoords*percent) > 5 ? (totalCoords*percent) : 5); // need at least 5 points in space to align

        std::set<unsigned int> useMe;
        unsigned int countIt = 0;
        while(countIt < subTotal){
            unsigned int tryMe = randomInt(gen);
            if (useMe.find(tryMe) == useMe.end()){
                useMe.insert(tryMe);
                countIt++;
            }
        }

        for(auto & usethis : useMe){

            const auto & target = tarCoordinates[usethis];
            dis = FLT_MAX;
            keep = 0;
            accepted = false;
            for(unsigned int r=0; r<totalRefCoords; r++){
                temp = (target-refcoordinates[r]).length();
                if (temp < dis && temp < cutoff){ // find closest point within cutoff
                    dis = temp;
                    keep = r;
                    accepted = true;
                }
            }
            // record the assignment
            if (accepted){
                tar_to_ref.insert(std::pair<unsigned int, unsigned int>(usethis, keep));
            }
        }


        if (tar_to_ref.size() < 5){
            while(tar_to_ref.size() < 5){
                tar_to_ref.clear();
                // increase cutoff and reselect
                cutoff += 0.1*cutoff;
                for(auto & usethis : useMe){

                    const auto & target = tarCoordinates[usethis];
                    dis = FLT_MAX;
                    keep = 0;
                    accepted = false;
                    for(unsigned int r=0; r<totalRefCoords; r++){
                        temp = (target-refcoordinates[r]).length();
                        if (temp < dis && temp < cutoff){ // find closest point within cutoff
                            dis = temp;
                            keep = r;
                            accepted = true;
                        }
                    }
                    // record the assignment
                    if (accepted){
                        tar_to_ref.insert(std::pair<unsigned int, unsigned int>(usethis, keep));
                        if (tar_to_ref.size() == 5){
                            break;
                        }
                    }
                }
            }
        }

        Eigen::MatrixXf fixed(tar_to_ref.size(),3);
        Eigen::MatrixXf moving(tar_to_ref.size(),3);

        int row=0;
        for (auto & element : tar_to_ref){
            const vector3 & tarVec = tarCoordinates[element.first];
            moving.row(row) << tarVec.x, tarVec.y, tarVec.z;

//            moving(row,0) = tarVec.x;
//            moving(row,1) = tarVec.y;
//            moving(row,2) = tarVec.z;

            const vector3 & refVec = refcoordinates[element.second];
//            fixed(row,0) = refVec.x;
//            fixed(row,1) = refVec.y;
//            fixed(row,2) = refVec.z;
            fixed.row(row) << refVec.x, refVec.y, refVec.z;
            row++;
        }

        cov = fixed.transpose() * moving;

    } else {
        // create correspondences
        // find closest neighboring distance
        for(unsigned int i=0; i<totalCoords; i++){
            const auto & target = tarCoordinates[i];
            dis = FLT_MAX;
            keep = 0;
            accepted = false;
            for(unsigned int r=0; r<totalRefCoords; r++){
                temp = (target-refcoordinates[r]).length();
                if (temp < dis && temp < cutoff){
                    dis = temp;
                    keep = r;
                    accepted = true;
                }
            }
            // record the assignment
            if (accepted){
                tar_to_ref.insert(std::pair<unsigned int, unsigned int>(i, keep));
            }
        }

        // assemble matrices
        const auto totalRows = (unsigned int)tar_to_ref.size();

        Eigen::MatrixXf fixed(totalRows,3);
        Eigen::MatrixXf moving(totalRows,3);
        for(auto & kv : tar_to_ref){
            const vector3 & tarVec = tarCoordinates[kv.first];
            moving.row(count) << tarVec.x, tarVec.y, tarVec.z;
//            moving(count,0) = tarVec.x;
//            moving(count,1) = tarVec.y;
//            moving(count,2) = tarVec.z;

            const vector3 & refVec = refcoordinates[kv.second];
//            fixed(count,0) = refVec.x;
//            fixed(count,1) = refVec.y;
//            fixed(count,2) = refVec.z;
            fixed.row(count) << refVec.x, refVec.y, refVec.z;
            count++;
        }

        cov = fixed.transpose() * moving;
    }

    // perform the SVD
    Eigen::JacobiSVD<Eigen::MatrixXf, Eigen::NoQRPreconditioner> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::MatrixXf s_matrix = svd.singularValues().asDiagonal();

    Eigen::MatrixXf reflection_matrix = Eigen::MatrixXf::Identity(3, 3);
    reflection_matrix(2, 2) = (svd.matrixV() * svd.matrixU().transpose()).determinant(); // reflection correction

    return (svd.matrixU() * reflection_matrix * svd.matrixV().transpose());

}


/*
 * Calculate variance of CVX hull points using SVD
 * then calculate distribution around each CVX hull point in each respective point set, doesn't matter if mirror is not used
 * Matrix is double precision
 *
 * moved matrix should cover entire molecule
 */
float Aligner::calculateVarianceScore(cpd::Matrix &moved, PointCloud &refModel, PointCloud &tarModel){

    const auto & tar_indices = tarModel.getAllIndices();

    auto totalRows = (unsigned int)tar_indices.size();
    std::vector<vector3> transVec(moved.rows());
    vector3 * const pTVec = transVec.data();
    for(long r=0; r<moved.rows(); r++){
        pTVec[r] = vector3((float)moved(r,0), (float)moved(r,1), (float)moved(r,2));
    }

    const auto & cvx_hull_ref = refModel.getAllIndices();

    const vector3 * refcoordinates = refModel.getKept() ? refModel.getTransformedCoordinates() : refModel.getCenteredCoordinates();

//    float minSum=0.0f;
    unsigned int rowCount=0;

    float sumAvg = 0.0f;
    float sumSqAvg = 0.0f;

    for(const auto & target : transVec){

        float dis = FLT_MAX;

        for(auto & rind : cvx_hull_ref){ // find closest CVX point in reference
            float temp = (refcoordinates[rind] - target).length();
            if (temp<dis){
                dis = temp;
            }
        }

        sumAvg += dis;
        sumSqAvg += dis*dis;

//        nearest_fixed.push_back(index);
//        nearest_dis.push_back(dis);
//        minSum += dis*dis;
        rowCount++;
    }

//    minSum *= 1.0f/(float)totalRows;

    sumSqAvg *= 1.0f/(float)totalRows;
    sumAvg *= 1.0f/(float)totalRows;

//    float minSum2 = 0.0f;
    float sumAvg2 = 0.0f;
    float sumSqAvg2 = 0.0f;

    const auto & ref_set = refModel.getAllIndices();
    for(auto & rind : ref_set){

        const vector3 & refvec = refcoordinates[rind];

        float dis = FLT_MAX;

        for(long i=0; i<moved.rows(); i++){
            float temp = (refvec - pTVec[i]).length();
            if (temp<dis){
                dis = temp;
            }
        }
        sumAvg2 += dis;
        sumSqAvg2 += dis*dis;
//        minSum2 += dis*dis;
    }

    sumSqAvg2 *= 1.0f/(float)ref_set.size();
    sumAvg2 *= 1.0f/(float)ref_set.size();

    //return (std::sqrt(0.5f*(minSum+minSum2/(float)ref_set.size())));
    return 0.5f*((sumSqAvg - sumAvg*sumAvg) + (sumSqAvg2 - sumAvg2*sumAvg2))/std::max(tarModel.getClosestDistance(), refModel.getClosestDistance());
}


/*
 * Calculate score using points that exclude 5 extreme CVX hull point
 * moved should be [3 x N] matrix
 */
float Aligner::calculateVarianceScore(Eigen::MatrixXf &moved, PointCloud &refModel, PointCloud &tarModel){

    const auto & tar_indices = tarModel.getAllIndices(); // remove the extreme points
    const auto totalPoints = (unsigned int)tar_indices.size();
    auto total_N = (unsigned int)moved.cols();

//    std::vector<vector3> transVec(total_N);
//    vector3 * const pTVec = transVec.data();
//
//    for(long r=0; r<total_N; r++){ // column vectors for each coordinate
//        auto colInUse = moved.col(r); // reference to the row
//        pTVec[r] = vector3(colInUse(0), colInUse(1), colInUse(2));
//    }

    const auto & ref_set = refModel.getAllIndices();
    const vector3 * refcoordinates = refModel.getKept() ? refModel.getTransformedCoordinates() : refModel.getCenteredCoordinates();

    auto total_N_ref = ref_set.size();
    std::vector<float> refDistances(total_N_ref);
    float * pRefDis = refDistances.data();
    float * pDis;

    std::fill(refDistances.begin(), refDistances.end(), FLT_MAX);
//    float minSum=0.0f;
    float sumAvg = 0.0f;
    float sumSqAvg = 0.0f;
    float maxT = FLT_MIN;
    float maxF = FLT_MIN;
//    int outcount1 = 0;

    // take half the average distance between the reference and target, two points are in the same space if within this metric
//    float metric = 0.5f*(0.5f*(refModel.getClosestDistance() + tarModel.getClosestDistance()));

    vector3 target;
    for(long r=0; r<total_N; r++){ // column vectors for each coordinate
        auto colInUse = moved.col(r); // reference to the row
        target.x = colInUse(0);
        target.y = colInUse(1);
        target.z = colInUse(2);
        //_mm_set_ps(0, z, y, x);

        float dis = FLT_MAX;

        for(auto & rind : ref_set){ // calculate squared length from single target to all in ref
            float temp = (refcoordinates[rind] - target).sqlength();
            if (temp<dis){ // finding the smallest difference between r_t and all in reference
                dis = temp;
            }

            pDis = &pRefDis[rind]; // for each point in Ref, what is the smallest distance to target
            if (temp < *pDis){
                *pDis = temp;
            }
        }

        sumSqAvg += dis;
        dis = sqrtf(dis);

        if (dis > maxT){ // largest difference between target and reference
            maxT = dis;
        }

        sumAvg += dis;
        //sumSqAvg += dis*dis;
    }


//    for(auto & target : transVec){
//
//        float dis = FLT_MAX;
//
//        for(auto & rind : ref_set){ // find closest CVX point in reference
//            float temp = (refcoordinates[rind] - target).sqlength();
//            if (temp<dis){
//                dis = temp;
//            }
//
//            if (temp < pRefDis[rind]){
//                pRefDis[rind] = temp;
//            }
//        }
//
//        if (dis > maxT){
//            maxT = dis;
//        }
//
////        auto met = (int)std::floor(dis/metric);
////        if (met > 0){
////            outcount1 += (met);
////        }
//
//        sumSqAvg += dis;
//        dis = sqrtf(dis);
//        if (dis > maxT){
//            maxT = dis;
//        }
//
//        sumAvg += dis;
////        sumSqAvg += dis*dis;
////        minSum += dis*dis;
//    }

    //minSum *= 1.0f/(float)totalPoints;
    sumSqAvg *= 1.0f/(float)totalPoints;
    sumAvg *= 1.0f/(float)totalPoints;


//    float minSum2 = 0.0f;
    float sumAvg2 = 0.0f;
    float sumSqAvg2 = 0.0f;
//    int outcount2 = 0;

//    for(auto & rind : ref_set){
//
//        const vector3 & refvec = refcoordinates[rind];
//        float dis = FLT_MAX;
//
//        for(auto & vec : transVec){
//            float temp = (refvec - vec).length();
//            if (temp<dis){
//                dis = temp;
//            }
//        }
//
//        if (dis > maxF){ // find longest shortest distance
//            maxF = dis;
//        }
//
////        auto met = (int)std::floor(dis/metric);
////        if (met > 0){ // if shortest distance is outside half lattice distance, penalize
////            outcount2 += met;
////        }
//
//        sumAvg2 += dis;
//        sumSqAvg2 += dis*dis;
////        minSum2 += dis*dis;
//    }

    float dis;
    for(unsigned int i=0; i<total_N_ref; i++){
        dis = pRefDis[i];
        sumSqAvg2 += dis;
        dis = sqrtf(dis);
        if (dis > maxF){
            maxF = dis;
        }
        sumAvg2 += dis;
//        sumSqAvg2 += dis*dis;
    }


    // want to penalize large distances
    //minSum2 *= 1.0f/(float)ref_set.size();
    sumSqAvg2 *= 1.0f/(float)total_N_ref;
    sumAvg2 *= 1.0f/(float)total_N_ref;

//    float minmax = std::max(tarModel.getTotalPointsInCloud(), refModel.getTotalPointsInCloud());
    float minmaxdis = std::max(tarModel.getClosestDistance(), refModel.getClosestDistance());

    /*
     * minimize number of points that are more than 1 unit away from neighboring set
     */

    //return (std::sqrt(0.5f*(minSum+minSum2/(float)ref_set.size())));
    // average variance normalized to max lattice distance
//    return 0.5f*(((sumSqAvg - sumAvg*sumAvg) + (sumSqAvg2 - sumAvg2*sumAvg2)))/minmaxdis;

    // minimize on NSD/average variance
//    return 0.3f*(std::sqrt(0.5f*(minSum/((float)tarModel.getTotalPointsInCloud()*tarModel.getClosestDistance()*tarModel.getClosestDistance()) + minSum2/((float)refModel.getTotalPointsInCloud()*refModel.getClosestDistance()*refModel.getClosestDistance())))) +
//            (0.5f*(((sumSqAvg - sumAvg*sumAvg) + (sumSqAvg2 - sumAvg2*sumAvg2))))/minmaxdis;

//    std::cout << (outcount1 + outcount2)/minmax << " " << ((totalRows*maxT+ref_set.size()*maxF)/((float)(ref_set.size()+totalRows)))/minmaxdis << std::endl;

    // minimize magnitude of isolated points and Hausdorff betweeen each set
    // first term should be < 1
    // std::cout << " " << 0.2f*(outcount1 + outcount2)/minmax << " " <<  ((total_N*maxT+ref_set.size()*maxF)/((float)(ref_set.size()+total_N)))/minmaxdis << std::endl;
    // return 0.2f*(outcount1 + outcount2)/minmax + ((total_N*maxT+ref_set.size()*maxF)/((float)(ref_set.size()+total_N)))/minmaxdis;

    // minimize magnitude of isolated points and NSD
    //std::cout << (outcount1 + outcount2)/minmax << " " << (std::sqrt(0.5*(minSum/((float)tarModel.getTotalPointsInCloud()*tarModel.getClosestDistance()*tarModel.getClosestDistance()) + minSum2/((float)refModel.getTotalPointsInCloud()*refModel.getClosestDistance()*refModel.getClosestDistance())))) << std::endl;
    //return 0.2f*(outcount1 + outcount2)/minmax + (std::sqrt(0.5f*(minSum/((float)tarModel.getTotalPointsInCloud()*tarModel.getClosestDistance()*tarModel.getClosestDistance()) + minSum2/((float)refModel.getTotalPointsInCloud()*refModel.getClosestDistance()*refModel.getClosestDistance()))));

    // variance and Hausdorrf
//       std::cout << (total_N*maxT+ref_set.size()*maxF)/((float)(ref_set.size()+total_N))/minmaxdis << "  "<< (((sumSqAvg - sumAvg*sumAvg) + (sumSqAvg2 - sumAvg2*sumAvg2)))/minmaxdis << std::endl;

    // NSD plus Hausdorff metric for each set
//    return (totalRows*maxT+ref_set.size()*maxF)/((float)(ref_set.size()+totalRows))/minmaxdis + (std::sqrt(0.5*(minSum/((float)tarModel.getTotalPointsInCloud()*tarModel.getClosestDistance()*tarModel.getClosestDistance()) + minSum2/((float)refModel.getTotalPointsInCloud()*refModel.getClosestDistance()*refModel.getClosestDistance())))) ;
//    return 0.1f*(total_N*maxT+ref_set.size()*maxF)/((float)(ref_set.size()+total_N))/minmaxdis + 0.5f*(((sumSqAvg - sumAvg*sumAvg) + (sumSqAvg2 - sumAvg2*sumAvg2)))/minmaxdis;
    return (total_N*maxT+total_N_ref*maxF)/((float)(total_N_ref+total_N)*minmaxdis) + 0.1f*(0.5f*(((sumSqAvg - sumAvg*sumAvg) + (sumSqAvg2 - sumAvg2*sumAvg2))))/minmaxdis;
}


/*
 * use the CVX hull of moved coordinates to determine the correspondence in the fixed coordinates
 * then calculate distribution around each CVX hull point in each respective point set, doesn't matter if mirror is not used
 */
float Aligner::getCorrespondenceScore(Eigen::MatrixXf &moved, PointCloud &refModel, PointCloud &tarModel){

    const auto & cvx_indices = tarModel.getCvx_hull();

    std::vector<unsigned int> nearest_fixed;
    std::vector<float> nearest_dis;

    size_t totalInRef = refModel.getTotalPointsInCloud();

    std::cout << " rows " << moved.rows() << " :: " << moved.cols() << std::endl;

    for(auto & cvx : cvx_indices){ // for each CVX point in target find closest in reference
        auto row = moved.row(cvx);
        vector3 target (row(0), row(1), row(2));

        // find closest point in fixed;
        float dis = FLT_MAX;
        float temp;
        auto cvx_hull_ref = refModel.getCvx_hull();
        auto coordinates = refModel.getKept() ? refModel.getTransformedCoordinates() : refModel.getCenteredCoordinates();
        unsigned int index;

        for(auto & cvxr : cvx_hull_ref){ // find closest CVX point in reference
            temp = (coordinates[cvxr] - target).length();
            if (temp<dis){
                dis = temp; // keep the distance
                index = cvxr;
            }
        }

        nearest_fixed.push_back(index);
        nearest_dis.push_back(dis);
    }


    const unsigned int upperIndex = 3;
    std::vector<float> moved_distribution(upperIndex);
    std::vector<float> fixed_distribution(upperIndex);

    float shortestDis = 1.0f/tarModel.getClosestDistance();
    float shortestDisfixed = 1.0f/refModel.getClosestDistance();

    float within = 0.5f*(tarModel.getClosestDistance()+refModel.getClosestDistance());

    auto coordinates = tarModel.getCoordinates();
    auto fixcoordinates = refModel.getCoordinates();

    float dkl = 0.0f;
    int count=0;
    for(auto & cvx : cvx_indices){ // for each CVX point in target, calculate its distribution of points locally
        const auto & target = coordinates[cvx];

        std::fill(moved_distribution.begin(), moved_distribution.end(), 0.0f);
        std::fill(fixed_distribution.begin(), fixed_distribution.end(), 0.0f);

        // find closest point in fixed;
        for(unsigned int i=0; i<cvx; i++){
            auto temp = (unsigned int)std::floor((coordinates[i] - target).length()*shortestDis);
            if (temp<upperIndex){
                moved_distribution[temp] += 1.0f;
            }
        }

        for(unsigned int i=(cvx+1); i<totalInRef; i++){
            auto temp = (unsigned int)std::floor((coordinates[i] - target).length()*shortestDis);
            if (temp<upperIndex){
                moved_distribution[temp] += 1.0f;
            }
        }

        // do the paired index in reference (fixed) model
        unsigned int fixedindex = nearest_fixed[count];
        const auto & fixedref = fixcoordinates[fixedindex];

        for(unsigned int i=0; i<fixedindex; i++){
            auto temp = (unsigned int)std::floor((fixcoordinates[i] - fixedref).length()*shortestDisfixed);
            if (temp<upperIndex){
                fixed_distribution[temp] += 1.0f;
            }
        }

        for(unsigned int i=(fixedindex+1); i<totalInRef; i++){
            auto temp = (unsigned int)std::floor((fixcoordinates[i] - fixedref).length()*shortestDisfixed);
            if (temp<upperIndex){
                fixed_distribution[temp] += 1.0f;
            }
        }

        for(unsigned int i=0;i<upperIndex; i++){
            if (moved_distribution[i] < FLT_EPSILON){
                moved_distribution[i] = FLT_EPSILON;
            }
        }
        // calculate metric
        float movSum=0;
        float fixSum=0;
        for(unsigned int i=0;i<upperIndex; i++){
            movSum += moved_distribution[i];
            fixSum += fixed_distribution[i];
        }


        float invFix = 1.0f/fixSum;
        for(unsigned int i=0;i<upperIndex; i++){
            float fixit = fixed_distribution[i]*invFix;

            if (fixit > FLT_EPSILON){
                dkl += fixit*log(fixit/moved_distribution[i]*movSum);
            }
        }

        count++;
    }


    float sumit = 0.0f;
    for(unsigned int i=0; i<nearest_dis.size(); i++){
        float diff = (nearest_dis[i]/within);
        sumit += diff*diff;
    }

    dkl *= sumit/cvx_indices.size();
    return dkl;
}




//void Aligner::heapAlgorithm(int maxindex){
//
//    std::vector<int> c_array(maxindex);
//    std::vector<int> permutationIndices(maxindex);
//    for(int i=0; i<maxindex; i++){
//        c_array[i] = 0;
//        permutationIndices[i] = i;
//    }
//
//
//    int i_index=0;
//    while(i_index < maxindex){
//        if ( c_array[i_index] < i_index){
//
//            if (i_index % 2 == 0){
//                std::iter_swap(permutationIndices.begin(), permutationIndices.begin()+i_index);
//            } else {
//                std::iter_swap(permutationIndices.begin() + c_array[i_index], permutationIndices.begin()+i_index);
//            }
//            // output A - copy to matrix
//
//            c_array[i_index] += 1;
//            i_index = 0;
//
//        } else {
//            c_array[i_index] = 0;
//            i_index += 1;
//        }
//    }
//}

// pass in permutation array, c_rray and i_index
// get back

void Aligner::heapAlgorithm(std::vector<unsigned int> & c_array, std::vector<unsigned int> & permutationIndices, unsigned int * i_index){

    if ( c_array[*i_index] < *i_index){

        if ((*i_index) % 2 == 0){
            std::iter_swap(permutationIndices.begin(), permutationIndices.begin()+(*i_index));
        } else {
            std::iter_swap(permutationIndices.begin() + c_array[*i_index], permutationIndices.begin()+(*i_index));
        }
        // output A - copy to matrix

        c_array[(*i_index)] += 1;
        *i_index = 0;

    } else {
        c_array[(*i_index)] = 0;
        *i_index += 1;
    }
}

unsigned long Aligner::factorial(unsigned int n_value){
    unsigned long startit = 1;
    for(unsigned long i = 1; i <=n_value; ++i) {
        startit *= i;
    }
    return startit;
}


bool Aligner::nSelectK(std::vector<unsigned int> & permutationIndices, unsigned int n_index){

    auto k_index = (unsigned int)permutationIndices.size();

    for(long i = k_index-1; i >=0; i--){
        if(permutationIndices[i] < (n_index - k_index + i + 1)){

            permutationIndices[i]++;

            for(long j = (i+1); j < k_index; j++){
                permutationIndices[j] = permutationIndices[j-1] + 1;
            }
            return true;
        }
    }

    return false;
}


unsigned long Aligner::nChoosek(unsigned int n_index, unsigned int k_index){
    unsigned long value = 1;

    for(unsigned int i=1; i <= k_index; i++){
        value *= (n_index + 1 - i)/i;
    }
    return value;
}





/*
 * Hausdorf is for each point in S_1, calculate minimum distance
 * Un-normalized calculation
 * We use un-normalized as the minmization doesn't depend on number of beads in each set
 * i.e., derivative with respect to constant is zero
 *
 * Take the minimum of the set.
 *
 */
float Aligner::unNormalizedHausdorfScore(Eigen::MatrixXf &transformed, PointCloud &refModel, PointCloud &tarModel){

    auto refCoordinates = refModel.getKept() ? refModel.getTransformedCoordinates() : refModel.getCenteredCoordinates();

    float temp, dis, sum = 0.0f, sum2 = 0.0f;

    auto totalRows = (unsigned int)transformed.rows();
    std::vector<vector3> transVec(totalRows);
    vector3 * const pTVec = transVec.data();

    // min between target to reference
    float haus1 = FLT_MIN, haus2 = FLT_MIN;

    for(unsigned int r=0; r<totalRows; r++){
        pTVec[r] = vector3((float)transformed(r,0), (float)transformed(r,1), (float)transformed(r,2));
        const vector3 & vec = pTVec[r];
        dis = FLT_MAX;
        for(unsigned int i = 0; i < refModel.getTotalPointsInCloud(); i++) {
            temp = (refCoordinates[i]-vec).length();
            if (temp < dis){
                dis = temp;
            }
        }
        sum += dis*dis; // sum of all atoms

        if (dis > haus1){
            haus1 = dis;
        }
    }

    sum *= 1.0f/(float)totalRows;

    totalRows = (unsigned int)refModel.getTotalPointsInCloud();
    // min between target to reference
    for(unsigned int r=0; r<totalRows; r++){
        const vector3 & refVec = refCoordinates[r];

        dis = FLT_MAX;

        for(unsigned int i = 0; i < tarModel.getTotalPointsInCloud(); i++) {
            temp = (pTVec[i]-refVec).length();
            if (temp < dis){
                dis = temp;
            }
        }
        sum2 += dis*dis;

        if (dis > haus2){
            haus2 = dis;
        }
    }
    //std::cout << " runtime " << (std::clock() - startTime)/(double) CLOCKS_PER_SEC << std::endl;
    sum2 *= 1.0f/(float)totalRows;
    return std::sqrt(0.5f*(sum + sum2));
}

/*
 * Hausdorf is for each point in S_1, calculate minimum distance
 * Un-normalized calculation
 * We use un-normalized as the minmization doesn't depend on number of beads in each set
 * i.e., derivative with respect to constant is zero
 *
 * Take the minimum of the set.
 *
 */
float Aligner::unNormalizedHausdorfScore(cpd::Matrix &transformed, PointCloud &refModel){

    auto refCoordinates = refModel.getKept() ? refModel.getTransformedCoordinates() : refModel.getCenteredCoordinates();

    float temp, dis, sum = 0.0f, sum2 = 0.0f;

    auto totalRows = (unsigned int)transformed.rows();
    auto totalRowsInRef = (unsigned int)refModel.getTotalPointsInCloud();
    // min between target to reference
    float haus1 = FLT_MIN, haus2 = FLT_MIN;

    std::vector<vector3> transVec(totalRows);
    vector3 * const pTVec = transVec.data();
    for(unsigned int r=0; r<totalRows; r++){
        pTVec[r] = vector3((float)transformed(r,0), (float)transformed(r,1), (float)transformed(r,2));
    }


    for(unsigned int r=0; r<totalRows; r++){

        const vector3 & tvec = pTVec[r];

        dis = FLT_MAX;
        for(unsigned int i = 0; i < totalRowsInRef; i++) {
            temp = (refCoordinates[i]-tvec).length();
            if (temp < dis){
                dis = temp;
            }
        }
        sum += dis*dis; // sum of all atoms

        if (dis > haus1){
            haus1 = dis;
        }
    }

    sum *= 1.0f/(float)totalRows;


    // min between target to reference
    for(unsigned int r=0; r<totalRowsInRef; r++){
        const vector3 & refVec = refCoordinates[r];

        dis = FLT_MAX;

        for(unsigned int i = 0; i < totalRows; i++) {
            temp = (pTVec[i]-refVec).length();
            if (temp < dis){
                dis = temp;
            }
        }
        sum2 += dis*dis;

        if (dis > haus2){
            haus2 = dis;
        }
    }
    sum2 *= 1.0f/(float)totalRowsInRef;
    std::cout << "HAUSDORFF " << haus1 << " " << haus2 << std::endl;
    return std::sqrt(0.5f*(sum + sum2));
}

/*
 * Hausdorf is for each point in S_1, calculate minimum distance
 * Take the minimum of the set.
 */
float Aligner::NSDScore(cpd::Matrix &transformed, cpd::Matrix &fixed, float dminTar, float dminRef){

    const auto totalInFixed = (unsigned int)fixed.rows();

    float temp, dis, sum = 0.0f, sum2=0.0f;

    const auto totalRows = (unsigned int)transformed.rows();

    // min between target to reference
    float haus1 = FLT_MIN, haus2 = FLT_MIN;

    std::clock_t startTime = std::clock();
    for(unsigned int r=0; r<totalRows; r++){
//        auto row = transformed.row(r);
//        const vector3 tarvec= vector3(row(0), row(1), row(2));
        const vector3 tarvec = vector3((float)transformed(r,0), (float)transformed(r,1), (float)transformed(r,2));

        dis = FLT_MAX;

        for (cpd::Matrix::Index i = 0; i < fixed.rows(); ++i) {
//            row = fixed.row(i);
//            temp = (vector3((float)row(0), (float)row(1), (float)row(2))-tarvec).length();
            temp = (vector3(fixed(i,0), fixed(i,1), fixed(i,2))-tarvec).length();
            if (temp < dis){
                dis = temp;
            }
        }

        sum += dis*dis;

        if (dis > haus1){
            haus1 = dis;
        }
    }


    sum *= 1.0f/((float)totalRows*dminTar*dminTar);

    // min between reference to target
    for(unsigned int r=0; r<totalInFixed; r++){
//        auto row = fixed.row(r);
//        const vector3 refVec= vector3(row(0), row(1), row(2));
        const vector3 & refVec = vector3(fixed(r,0), fixed(r,1), fixed(r,2));

        dis = FLT_MAX;

        for(long i = 0; i < totalRows; i++) {
//            row = transformed.row(i);
            temp = (vector3((float)transformed(i,0), (float)transformed(i,1), (float)transformed(i,2))-refVec).length();
//            temp = (vector3((float)row(0), (float)row(1), (float)row(2))-refVec).length();
            if (temp < dis){
                dis = temp;
            }
        }
        sum2 += dis*dis;

        if (dis > haus2){
            haus2 = dis;
        }
    }
    sum2 *= 1.0f/((float)totalInFixed*dminRef*dminRef);
    std::cout << " runtime " << (std::clock() - startTime)/(double) CLOCKS_PER_SEC << std::endl;
    return std::sqrt(0.5f*(sum+sum2));
}





/**
 * calculate correspondences excluding CVX hull points
 *
 * @param refModel
 * @param tarModel
 * @return
 */
float Aligner::calculateAlignmentErrorCVX(PointCloud &refModel, PointCloud &tarModel) {

    auto nonCvxPointsTarget = tarModel.getNon_Cvx_hull_Subset();

    auto targetCoordinates = tarModel.getKept() ? tarModel.getTransformedCoordinates() : tarModel.getCenteredCoordinates();
    auto refCoordinates = refModel.getKept() ? refModel.getTransformedCoordinates() : refModel.getCenteredCoordinates();

    float score = 0;

    for(auto & pt : nonCvxPointsTarget){

        // find closest point in fixed;
        float dis = FLT_MAX;
        float temp;
        const auto & target = targetCoordinates[pt];

        for(unsigned int i = 0; i < refModel.getTotalPointsInCloud(); i++) {
            temp = (refCoordinates[i]-target).length();
            if (temp < dis){
                dis = temp;
            }
        }

        score += dis*dis;
    }

    return score;
}


void Aligner::fillCrossValidationMatrix(Eigen::MatrixXf & matrix, PointCloud & tarModel){

    auto & cvset = tarModel.getCrossValidationSet();
    auto coords = (tarModel.getKept()) ? tarModel.getTransformedCoordinates() : tarModel.getCenteredCoordinates();

    int rowCount = 0;
    for(auto ind : cvset){
        const vector3 & vec = coords[ind];
        matrix(rowCount,0) = vec.x;
        matrix(rowCount,1) = vec.y;
        matrix(rowCount,2) = vec.z;

        rowCount++;
    }
}



Eigen::MatrixXf Aligner::generateRotationMatrix(float theta, float phi, float psi){

    //float twoPI = (float)(2.0*M_PI);
    float costheta = cos(theta);
    float sintheta = sin(theta);
    float cosphi = cos(phi);
    float sinphi = sin(phi);
    float cospsi = cos(psi);
    float sinpsi = sin(psi);

    float sinphi_costheta= sinphi*costheta;
    float sinphi_sintheta= sinphi*sintheta;

//    Eigen::Matrix3f rx; //roll
//    rx << 1.0f, 0.0f, 0.0f, 0.0f, costheta, -sintheta, 0.0f, sintheta, costheta;
//
//    Eigen::Matrix3f ry; //pitch
//    ry << cosphi, 0.0f, sinphi, 0.0f, 1.0f, 0.0f, -sinphi, 0.0f, cosphi;
//
//    Eigen::Matrix3f rz; //yaw
//    rz << cospsi, -sinpsi, 0.0f, sinpsi, cospsi, 0.0f, 0.0f, 0.0f, 1.0f;
//
    Eigen::Matrix3f rxyz;
    rxyz << cospsi*cosphi, cospsi*sinphi_sintheta-sinpsi*costheta, cospsi*sinphi_costheta + sinpsi*sintheta,
            sinpsi*cosphi, sinpsi*sinphi_sintheta+cospsi*costheta, sinpsi*sinphi_costheta - cospsi*sintheta,
            -sinphi, cosphi*sintheta, cosphi*costheta;

    //return rx*ry*rz;
    return rxyz;
}

Eigen::MatrixXf Aligner::getRandomRotation() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> randomDis(0.0, 1);

    Eigen::MatrixXf reflection_matrix = Eigen::MatrixXf::Identity(3, 3);
    Eigen::MatrixXf rot(3,3);

    auto  theta = (float)(2.0*randomDis(gen)*M_PI);
    auto phi = (float)(2.0*M_PI*randomDis(gen));
    auto zed = (float)(2.0*randomDis(gen));

    float sqrtz = sqrt(zed);
    float sint = sin(theta);
    float cost = cos(theta);

    Eigen::Vector3f vector(sin(phi)*sqrtz, cos(phi)*sqrtz, sqrt(2.0f-zed));

    rot(0,0) = cost;
    rot(0,1) = sint;
    rot(0,2) = 0.0f;
    rot(1,0) = -sint;
    rot(1,1) = cost;
    rot(1,2) = 0.0f;
    rot(2,0) = 0.0f;
    rot(2,1) = 0.0f;
    rot(2,2) = 1.0f;

    return ((vector*vector.transpose()-reflection_matrix)*rot);
}


Eigen::MatrixXf Aligner::randomRotation(Eigen::MatrixXf toBeTransformed) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> randomDis(0.0, 1);

    Eigen::MatrixXf reflection_matrix = Eigen::MatrixXf::Identity(3, 3);
    Eigen::MatrixXf rot(3,3);

    auto  theta = (float)(2.0*randomDis(gen)*M_PI);
    auto phi = (float)(2.0*M_PI*randomDis(gen));
    auto zed = (float)(2.0*randomDis(gen));

    float sqrtz = sqrt(zed);
    float sint = sin(theta);
    float cost = cos(theta);

    Eigen::Vector3f vector(sin(phi)*sqrtz, cos(phi)*sqrtz, sqrt(2.0f-zed));

    rot(0,0) = cost;
    rot(0,1) = sint;
    rot(0,2) = 0.0f;
    rot(1,0) = -sint;
    rot(1,1) = cost;
    rot(1,2) = 0.0f;
    rot(2,0) = 0.0f;
    rot(2,1) = 0.0f;
    rot(2,2) = 1.0f;

    return ((vector*vector.transpose()-reflection_matrix)*rot)*toBeTransformed;
}


void Aligner::writeLogFile(std::string name, std::vector<ScoreSet> & scores, std::vector<ScoreSet> & paired){

    FILE * pFile;

    pFile = fopen(name.c_str(), "w");

    fprintf(pFile,"REMARK  AlignIt ALIGNMENT \n");
    fprintf(pFile,"REMARK             REFERENCE : %s\n", pointClouds[scores[0].ref].getFilename().c_str());
    fprintf(pFile, "Index :: Input file \n");

    for (size_t n=0; n < pointClouds.size(); n++) {
        fprintf(pFile, "%5i :: %s \n", (unsigned int)(n+1), pointClouds[n].getFilename().c_str());
    }

    fprintf(pFile,"REMARK  \n");
    if (useCPD){
        fprintf(pFile,"REMARK  USING CPD\n");
    } else {
        fprintf(pFile,"REMARK  USING ICP\n");
    }
    fprintf(pFile,"REMARK  PAIRWISE SCORES \n");
    fprintf(pFile,"REMARK  \n");
    int index = 1;

    for(auto & scr : scores){
        std::string flag = (scr.mirror > 0) ? " TRUE" : "FALSE";
        fprintf(pFile, "%-4i :: %4i - %-4i %s %.5f %.5f \n", index, (scr.ref+1), (scr.target+1), flag.c_str(), scr.score, scr.nsd);
        index++;
    }

    fprintf(pFile,"REMARK  \n");
    fprintf(pFile,"REMARK  FINAL PAIRINGS \n");
    fprintf(pFile,"REMARK  \n");

    index = 1;
    for(auto & scr : paired){
        std::string flag = (scr.mirror > 0) ? " TRUE" : "FALSE";
        fprintf(pFile, "%-4i :: %4i - %-4i %s %.4f %.3f %s (%i) - %s (%i) \n",
                index,
                (scr.ref+1),
                (scr.target+1),
                flag.c_str(),
                scr.score,
                scr.nsd,
                pointClouds[scr.ref].getFilename().c_str(),
                (scr.ref+1),
                pointClouds[scr.target].getFilename().c_str(),
                (scr.target+1)
        );
        index++;
    }

    fclose(pFile);
}


void Aligner::writeKDEFileList(std::string basename, std::vector<std::string> & files){

    FILE * pFile;
    std::string name = "kde_"+basename+".inp";
    pFile = fopen(name.c_str(), "w");

    for(auto & file : files){

        fprintf(pFile, "%s \n", file.c_str());
    }

    fclose(pFile);
}


void Aligner::align(ScoreSet & scored, PointCloud & referenceModel, PointCloud & tarModel, bool updateTarModel){

    if (useCPD){
        cpdAlgorithm(referenceModel, tarModel);
    } else {
        ICPCutoff(referenceModel, tarModel);
    }

    scored.score = tarModel.getScore();
    scored.nsd = tarModel.getNSD();

    if (tryBothEnantiomorphs){
        PointCloud mirrorImage(tarModel);
        mirrorImage.makeMirror();
        mirrorImage.setKept(false); // use centered coordinates
        std::cout << " **********                 ********** "  << std::endl;
        std::cout << "             TESTING MIRROR "  << std::endl;
        std::cout << " **********                 ********** "  << std::endl;

        if (useCPD){
            cpdAlgorithm(referenceModel, mirrorImage);
        } else {
            ICPCutoff(referenceModel, mirrorImage);
        }

        if (mirrorImage.getScore() < tarModel.getScore()){
//        if (pc.getNSD() < tarModel.getNSD()){
            std::cout << "             :: CHOOSING MIRROR "  << std::endl;

//            if (scored.nsd < mirrorImage.getNSD()){
//                std::cout << "          :: NSD "  << scored.nsd << " < " << mirrorImage.getNSD() << std::endl;
//            }

            if (updateTarModel){
                tarModel.makeMirror();
                tarModel.replaceTransformedCoordinates(mirrorImage.getTransformedCoordinates());

                tarModel.setRotationMatrix(const_cast<Eigen::MatrixXf &>(mirrorImage.getRotationMatrix()));

                tarModel.setScore(mirrorImage.getScore());
                tarModel.setNSD(mirrorImage.getNSD());
            }

            scored.mirror = true;
            scored.score = mirrorImage.getScore();
            scored.nsd = mirrorImage.getNSD();
        } else {
            scored.mirror = false;
        }
    }
}






void Aligner::printVector(const vector3 & vec){
    std::cout << "VEC :: " << vec.x << " " << vec.y << " " << vec.z << std::endl;
}

void Aligner::AddToCoordinateVector(std::vector<vector3> & fillme, const vector3 * pVecs, unsigned int total){

    for(unsigned int i=0; i<total; i++){
        fillme.emplace_back(pVecs[i]);
    }
}


void Aligner::writeCoordinateVectorToFile(std::vector<vector3> & vectorOfCoordinates) {

    std::string residue_index;

    std::string name = "aligned_set.pdb";
    FILE * pFile;
    pFile = fopen(name.c_str(), "w");
    unsigned int totalAtoms = vectorOfCoordinates.size();
    vector3 * pVec;

    fprintf(pFile,"REMARK  ACCEPTED COORDINATES : %s\n", name.c_str());
    fprintf(pFile,"REMARK 265                   EDGE RADIUS : %.3f\n", this->edge_radius);

    for (unsigned int n=0; n < totalAtoms; n++) {
        residue_index = boost::lexical_cast<std::string>(n+1);
        pVec = &vectorOfCoordinates[n];
        fprintf(pFile, "ATOM  %5i  CA  ALA A%4s    %8.3f%8.3f%8.3f  1.00100.00\n", n+1,residue_index.c_str(), pVec->x,pVec->y,pVec->z );
    }
    fprintf(pFile,"END\n");
    fclose(pFile);
}

float Aligner::getEpsilon(int total){

    if (total >= 600000){
        return 214;
    } else if (total >= 400000){
        return 75;
    } else if (total >= 11000){
        return 27;
    } else if (total >= 890){
        return 10;
    } else if (total >= 177){
        return 3.33;
    } else if (total >= 24){
        return 1.33;
    }

    return 0.33;
}

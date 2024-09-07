//
// Authror: Basit M. Imran.
// Date : 2024-04-17
// Copyright (c) Hybrid Dynamic Systems and Robot Locomotion Lab, Virginia Tech
//

#include "raisim/OgreVis.hpp"
#include "randyImguiPanel.hpp"
#include "raisimKeyboardCallback.hpp"
#include "helper.hpp"
#include "helper2.hpp"
#include "nmpc.hpp"
#include "OtherUtils.hpp"
#include "timer.h"
#include <cstdlib> // For rand() and srand()
#include <ctime>   // For time()
#include <fstream>
#include <iostream>

using std::cout;
using std::cin;

bool first_time0 = 1;
bool sec_time = 1;

void setupCallback() {
    raisim::OgreVis *vis = raisim::OgreVis::get();

    /// light
    vis->getLight()->setDiffuseColour(1, 1, 1);
    vis->getLight()->setCastShadows(false);
    Ogre::Vector3 lightdir(-3,3,-0.5); // Light shines on ROBOTS top/front/right side
    // Ogre::Vector3 lightdir(-3,-3,-0.5); // Light shines on ROBOTS top/front/left side
    lightdir.normalise();
    vis->getLightNode()->setDirection({lightdir});
    vis->setCameraSpeed(300);

    vis->addResourceDirectory(raisim::loadResource("material"));
    vis->loadMaterialFile("myMaterials.material");

    vis->addResourceDirectory(vis->getResourceDir() + "/material/skybox/violentdays");
    vis->loadMaterialFile("violentdays.material");

    /// shdow setting
    vis->getSceneManager()->setShadowTechnique(Ogre::SHADOWTYPE_TEXTURE_ADDITIVE);
    vis->getSceneManager()->setShadowTextureSettings(2048, 3);

    /// scale related settings!! Please adapt it depending on your map size
    // beyond this distance, shadow disappears
    vis->getSceneManager()->setShadowFarDistance(10);
    // size of contact points and contact forces
    vis->setContactVisObjectSize(0.03, 0.6);
    // speed of camera motion in freelook mode
    vis->getCameraMan()->setTopSpeed(5);
}

void plotGRFs(std::map<std::string, raisim::VisualObject>& objList, const std::vector<double>& GRFs, const std::vector<double>& pf, const std::vector<double>& contact) {
    // Ensure the vectors are of the correct size
    if (GRFs.size() < 12 || pf.size() < 12 || contact.size() < 4) {
        std::cerr << "Error: Input vectors are of incorrect size." << std::endl;
        return;
    }

    for (int i = 0; i < 4; ++i) {
        raisim::Vec<3> dir;
        for (int j = 0; j < 3; ++j) {
            dir[j] = GRFs[3 * i + j];
        }

        // Normalize the direction vector if it is not a zero vector
        double norm = dir.norm();
        if (norm > 1e-6) {  // Check if the vector is non-zero to avoid division by zero
            dir /= norm;
        } else {
            dir.setZero();
        }

        // Convert direction vector to a rotation matrix that aligns the z-axis with the direction vector
        raisim::Mat<3, 3> rot;
        if (contact[i] == 1 && norm > 1e-6) {
            raisim::zaxisToRotMat(dir, rot);
        } else {
            rot.setIdentity();  // Set rotation to identity if no contact or zero norm
        }

        // Visual object key
        std::string objKey = "GRF" + std::to_string(i + 1);
        objList[objKey].offset = {pf[3 * i], pf[3 * i + 1], pf[3 * i + 2]};
        objList[objKey].scale = {0.2, 0.2, 0.005 * norm};  // Scaling based on the norm of the GRF vector

        // Set the rotation matrix
        objList[objKey].rotationOffset = rot;
    }
}

int main(int argc, char *argv[]) {


    NMPC* nmpc;
    nmpc = new NMPC();
    // nmpc->code_gen();  // For now I have a separate subfolder for code generation. Any codegen changes should be made there.
    nmpc->init();

    const int NUMBER_OF_SIMS = 1;
    // ============================================================ //
    // =================== SETUP RAISIM/VISUALS =================== //
    // ============================================================ //
    /// create raisim world

    raisim::World::setActivationKey(raisim::loadResource("activation.raisim"));
    raisim::World world;
    world.setTimeStep(simfreq_raisim);

    raisim::OgreVis *vis = raisim::OgreVis::get();

    /// these method must be called before initApp
    vis->setWorld(&world);
    vis->setWindowSize(1792, 1200); // Should be evenly divisible by 16!!
    vis->setImguiSetupCallback(imguiSetupCallback); // These 2 lines make the interactable gui visible
    vis->setImguiRenderCallback(imguiRenderCallBack);
    vis->setKeyboardCallback(raisimKeyboardCallback);
    vis->setSetUpCallback(setupCallback);
    vis->setAntiAliasing(2);

    /// starts visualizer thread
    vis->initApp();
    /// create raisim objects
    raisim::TerrainProperties terrainProperties;
    terrainProperties.frequency = 0.0;
    terrainProperties.zScale = 0.0;
    terrainProperties.xSize = 300.0;
    terrainProperties.ySize = 300.0;
    terrainProperties.xSamples = 50;
    terrainProperties.ySamples = 50;
    terrainProperties.fractalOctaves = 0;
    terrainProperties.fractalLacunarity = 0.0;
    terrainProperties.fractalGain = 0.0;

    raisim::HeightMap *ground = world.addHeightMap(0.0, 0.0, terrainProperties);
    vis->createGraphicalObject(ground, "terrain", "checkerboard_blue");
    world.setDefaultMaterial(0.8, 0.0, 0.0); //surface friction could be 0.8 or 1.0
    vis->addVisualObject("extForceArrow", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    // Foot force visualization arrows
    vis->addVisualObject("GRF1", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF2", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF3", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);
    vis->addVisualObject("GRF4", "arrowMesh", "red", {0.0, 0.0, 0.0}, false, raisim::OgreVis::RAISIM_OBJECT_GROUP);

    auto& list = vis->getVisualObjectList();

    // ============================================================ //
    // ======================= SETUP Robot ======================== //
    // ============================================================ //
    std::vector<raisim::ArticulatedSystem*> A1;
    A1.push_back(world.addArticulatedSystem(raisim::loadResource("A1/A1_trunk.urdf")));

    vis->createGraphicalObject(A1.back(), "A1");
    A1.back()->setName("A1_Robot");
    
    // Episodic setup. Runs for Episodes = NUMBER_OF_SIMS
    for (size_t sim = 0; sim < NUMBER_OF_SIMS; ++sim)
    {
        Eigen::Matrix<double, 2*NUMBER_OF_AGENTS, 1> Pstart; 
    
        
        bool checkifrunMPC0 = 0;
        first_time0 = 1;
        sec_time = 1;

        // ============================================================ //
        // =================== SETUP ENVIRONMENT ====================== //
        // ============================================================ //

        // Agent coordinates for full-order URDF of the robot
        // {Pstart(0), Pstart(1), 0.12, 1, 0, 0, 0,0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6, 0.0, Pi/3, -2.6}

        // We can change the start position of the robot here
        // A1.back()->setGeneralizedCoordinate({Pstart(0), Pstart(1), 0.12, 1, 0, 0, 0});

        // For now start position is set to 0,0,0.8
        A1.back()->setGeneralizedCoordinate({0, 0, 0.8, 1, 0, 0, 0});
        A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);

        // jointTorqueFF = Eigen::Map< Eigen::Matrix<double,18,1> >(tau,18);
        // jointTorqueFF.block(0,0,6,1).setZero();

        // // Set the desired torques
        // A1.back()->setControlMode(raisim::ControlMode::FORCE_AND_TORQUE);

        // A1.back()->setGeneralizedForce(jointTorqueFF);



        // ============================================================ //
        // ================= VIEW AND RECORDING OPTIONS =============== //
        // ============================================================ //
        raisim::gui::showContacts = false;
        raisim::gui::showForces = false;
        raisim::gui::showCollision = false;
        raisim::gui::showBodies = true;

        std::string cameraview = "iso";
        bool panX = true;                // Pan view with robot during walking (X direction)
        bool panY = true;                // Pan view with robot during walking (Y direction)
        bool record = false;             // Record?
        double startTime = 0*ctrlHz;    // RecoPstartrding start time
        double simlength = nmpc->params.sim_length/nmpc->params.dt; //60*ctrlHz;   // Sim end time
        double fps = 30;            
        std::string directory = "/Sim_Recordings/";
        // std::string filename = "Payload_Inplace";
        std::string filename = "SRB_NMPC";
        // std::string filename = "inplace_sim";

        // ============================================================ //
        // ========================= VIEW SETUP ======================= //
        // ============================================================ //
        if(cameraview == "iso"){
            vis->getCameraMan()->getCamera()->setPosition(-1, -3, 0.5);
            vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(4.5*Pi/6-Pi/2));
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
        }else if(cameraview == "isoside"){
            vis->getCameraMan()->getCamera()->setPosition(1.1, -2, 0.5);
            vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(4*Pi/6-Pi/2));
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
        }else if(cameraview == "side"){
            vis->getCameraMan()->getCamera()->setPosition(0, -2, 0.5);
            vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(0));
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
        }else if(cameraview == "front"){
            vis->getCameraMan()->getCamera()->setPosition(2, 0, 0.5);
            vis->getCameraMan()->getCamera()->yaw(Ogre::Radian(Pi/2));
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(Pi/2));
        }else if(cameraview == "top"){
            vis->getCameraMan()->getCamera()->setPosition(2, -1, 6);
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(0));
        }else{
            vis->getCameraMan()->getCamera()->setPosition(1, -3, 2.5);
            vis->getCameraMan()->getCamera()->pitch(Ogre::Radian(1.0));
        }
        unsigned long mask = 0;
        if(raisim::gui::showBodies) mask |= raisim::OgreVis::RAISIM_OBJECT_GROUP;
        if(raisim::gui::showCollision) mask |= raisim::OgreVis::RAISIM_COLLISION_BODY_GROUP;
        if(raisim::gui::showContacts) mask |= raisim::OgreVis::RAISIM_CONTACT_POINT_GROUP;
        if(raisim::gui::showForces) mask |= raisim::OgreVis::RAISIM_CONTACT_FORCE_GROUP;
        vis->setVisibilityMask(mask);

        if(panX) raisim::gui::panViewX = panX;
        if(panY) raisim::gui::panViewY = panY;
        
        // ============================================================ //
        // ========================== RUN SIM ========================= //
        // ============================================================ //
        const std::string name = directory+filename+"_"+cameraview+".mp4";
        vis->setDesiredFPS(fps);
        long simcounter = 0;
        static bool added = false;

        while (!vis->getRoot()->endRenderingQueued() && simcounter < simlength-10){

            nmpc->run_NMPC();

            std::vector<double> com_pos = nmpc->get_com_pos();
            std::vector<double> GRF = nmpc->get_GRF();
            std::vector<double> contacts = nmpc->get_contacts();
            std::vector<double> feet_vec = nmpc->get_Pf();
            A1.back()->setGeneralizedCoordinate({com_pos[0], com_pos[1], com_pos[2], 1, 0, 0, 0});

            plotGRFs(list, GRF, feet_vec, contacts);

            world.integrate();        
            
            if (simcounter%3 == 0)
                vis->renderOneFrame();
            
            if (!vis->isRecording() & record & simcounter>=startTime)
                vis->startRecordingVideo(name);
            
            auto currentPos = vis->getCameraMan()->getCamera()->getPosition();
            if (raisim::gui::panViewX){
                Eigen::VectorXd jointPosTotal(18 + 1);
                Eigen::VectorXd jointVelTotal(18);
                jointPosTotal.setZero();
                jointVelTotal.setZero();
                A1.back()->getState(jointPosTotal, jointVelTotal);
                if (cameraview=="front"){
                    currentPos[0] = jointPosTotal(0)+2;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="side"){
                    currentPos[0] = jointPosTotal(0);
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="iso"){
                    currentPos[0] = jointPosTotal(0)+2;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="isoside"){
                    currentPos[0] = jointPosTotal(0)+1.1;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                }
            }
            if (raisim::gui::panViewY){
                Eigen::VectorXd jointPosTotal(18 + 1);
                Eigen::VectorXd jointVelTotal(18);
                jointPosTotal.setZero();
                jointVelTotal.setZero();
                A1.back()->getState(jointPosTotal, jointVelTotal);
                if (cameraview=="front"){
                    currentPos[1] = jointPosTotal(1);
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="side"){
                    currentPos[1] = jointPosTotal(1)-2;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="iso"){
                    currentPos[1] = jointPosTotal(1)-1;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                } else if(cameraview=="isoside"){
                    currentPos[1] = jointPosTotal(1)-2;
                    vis->getCameraMan()->getCamera()->setPosition(currentPos);
                }
            }
            simcounter++;
        }

        nmpc->logData();
        
    }

    // End recording if still recording
    if (vis->isRecording())
        vis->stopRecordingVideoAndSave();
    
    /// terminate the app
    vis->closeApp();
    return 0;
}
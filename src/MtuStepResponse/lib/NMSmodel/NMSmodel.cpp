/* -------------------------------------------------------------------------- *
 * CEINMS is a standalone toolbox for neuromusculoskeletal modelling and      *
 * simulation. CEINMS can also be used as a plugin for OpenSim either         *
 * through the OpenSim GUI or API. See https://simtk.org/home/ceinms and the  *
 * NOTICE file for more information. CEINMS development was coordinated       *
 * through Griffith University and supported by the Australian National       *
 * Health and Medical Research Council (NHMRC), the US National Institutes of *
 * Health (NIH), and the European Union Framework Programme 7 (EU FP7). Also  *
 * see the PROJECTS file for more information about the funding projects.     *
 *                                                                            *
 * Copyright (c) 2010-2015 Griffith University and the Contributors           *
 *                                                                            *
 * CEINMS Contributors: C. Pizzolato, M. Reggiani, M. Sartori,                *
 *                      E. Ceseracciu, and D.G. Lloyd                         *
 *                                                                            *
 * Author(s): C. Pizzolato, M. Reggiani, D.G. Lloyd, M. Sartori               *
 *                                                                            *
 * CEINMS is licensed under the Apache License, Version 2.0 (the "License").  *
 * You may not use this file except in compliance with the License. You may   *
 * obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0.*
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include <iostream>
using std::cout;
using std::endl;
#include <vector>
using std::vector;
#include <map>
using std::map;
#include "ceinms/NMSmodel.h"
#include <string>
using std::string;
#include "ceinms/DoF.h"
#include "ceinms/MuscleParameters.h"

#include <algorithm>

//#define DEBUG
#define CHECKS

namespace ceinms {



    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    NMSmodel<Activation, Tendon, mode>::NMSmodel() :
        hasOpenSimModel_(false)
    {}

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    NMSmodel<Activation, Tendon, mode>::NMSmodel(const NMSmodel<Activation, Tendon, mode>& orig) :
        muscles_(orig.muscles_), 
        muscleNames_(orig.muscleNames_), 
        dofNames_(orig.dofNames_), 
        hasOpenSimModel_(orig.hasOpenSimModel_),
        osimModelFilename_(orig.osimModelFilename_)
    {

        std::map<std::string, std::vector<std::string>> dofsToMuscles;
        orig.getMuscleNamesOnDofs(dofsToMuscles);
        for (auto& d : dofNames_) {
            DoFtype dof(d);
            for (auto& m : dofsToMuscles[d]) {
                auto i(std::distance(muscleNames_.begin(), std::find(muscleNames_.begin(), muscleNames_.end(), m)));
                dof.addNewMuscle(muscles_.begin() + i);
            }
            dofs_.emplace_back(dof);
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::addMuscle(const MTUtype& muscle) {

        muscles_.push_back(muscle);
        muscleNames_.push_back(muscle.getMuscleName());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::addDoF(const DoFtype& dof) {

        dofs_.push_back(dof);
        dofNames_.push_back(dof.getName());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    bool NMSmodel<Activation, Tendon, mode>::haveThisMuscle(const string& currentMuscle,
        vectorMTUitr& found) {

        found = muscles_.begin();
        while ((found != muscles_.end()) && (!found->compareMusclesId(currentMuscle)))
            found++;

        if (found == muscles_.end())
            return false;
        else
            return true;
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    bool NMSmodel<Activation, Tendon, mode>::haveTheseMuscles(const vector<string>& musclesNames) {

        vectorMTUitr found;
        for (vector<string>::const_iterator mIt = musclesNames.begin(); mIt != musclesNames.end(); ++mIt)
            if (!haveThisMuscle(*mIt, found))
                return false;
        return true;
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    bool NMSmodel<Activation, Tendon, mode>::compareMusclesNames(const vector<string>& muscleNames) const {

        vector<string>::const_iterator m1 = muscleNames.begin();
        for (vectorMTUconstItr m2 = muscles_.begin();
            m2 != muscles_.end();) {

            if (!(*m2).compareMusclesId(*m1))
                return false;
            m1++; m2++;
        }

        return true;
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    bool NMSmodel<Activation, Tendon, mode>::compareMusclesNamesOnDoF(const vector<string>& muscleNames,
        int whichDoF) const {

        return dofs_.at(whichDoF).compareMusclesNames(muscleNames);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setTime(const double& time) {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt; muscleIt != muscles_.end(); ++muscleIt)
            muscleIt->setTime(time);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setEmgs(const std::vector<double>& currentEmgData) {

#ifdef CHECKS
        if (currentEmgData.size() != muscles_.size()) {
            cout << "Error: the number of muscle excitations is different from the number of MTUs in the model\n";
            cout << "Found " << currentEmgData.size() << " excitations and " << muscles_.size() << " MTUs\n";
            exit(EXIT_FAILURE);
        }
#endif

        vector<double>::const_iterator emgIt;
        vectorMTUitr muscleIt = muscles_.begin();
        for (emgIt = currentEmgData.begin(); emgIt < currentEmgData.end(); ++emgIt) {
            muscleIt->setEmg(*emgIt);
            ++muscleIt;
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setEmgsSelective(const std::vector<double>& currentEmgData,
        const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<double>::const_iterator emgIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            emgIt = currentEmgData.begin() + *it;
            muscleIt->setEmg(*emgIt);
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setTime_emgs_updateActivations_pushState_selective(double time,
        const vector<double>& currentEmgData,
        const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<double>::const_iterator emgIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            emgIt = currentEmgData.begin() + *it;
            muscleIt->setTime(time);
            muscleIt->setEmg(*emgIt);
            muscleIt->updateActivation();
            muscleIt->pushState();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setMuscleTendonLengths(const vector<double>& currentLmtData) {

#ifdef CHECKS
        if (currentLmtData.size() != muscles_.size()) {
            cout << "Error: the number of musculotendon lengths is different from the number of MTUs in the model\n";
            cout << "Found " << currentLmtData.size() << " musculotendon lengths and " << muscles_.size() << " MTUs\n";
            exit(EXIT_FAILURE);
        }
#endif
        vector<double>::const_iterator lmtIt;
        vectorMTUitr muscleIt = muscles_.begin();
        for (lmtIt = currentLmtData.begin(); lmtIt < currentLmtData.end(); ++lmtIt) {
            muscleIt->setMuscleTendonLength(*lmtIt);
            ++muscleIt;
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setMuscleTendonLengthsSelective(const std::vector<double>& currentLmtData,
        const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<double>::const_iterator lmtIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            lmtIt = currentLmtData.begin() + *it;
            muscleIt->setMuscleTendonLength(*lmtIt);
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setMomentArms(const vector<double>& currentMaData, unsigned whichDof) {

        dofs_.at(whichDof).setMomentArms(currentMaData);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateState() {

        updateActivations();
        updateFibreLengthsAndVelocities();
        updateMuscleForces();
        updateTorques();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateState(const vector<unsigned>& selectedmusclesIndex) {

        updateActivations(selectedmusclesIndex);
        updateFibreLengthsAndVelocities(selectedmusclesIndex);
        updateMuscleForces(selectedmusclesIndex);
        updateTorques();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateState_OFFLINE() {

        updateActivations();
        updateFibreLengthsAndVelocities_OFFLINE();
        updateMuscleForces();
        updateTorques();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateState_OFFLINE(const vector<unsigned>& selectedmusclesIndex) {

        updateActivations(selectedmusclesIndex);
        updateFibreLengthsAndVelocities_OFFLINE(selectedmusclesIndex);
        updateMuscleForces(selectedmusclesIndex);
        updateTorques();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateState_HYBRID() {

        updateActivations();
        updateFibreLengthsAndVelocities_HYBRID();
        updateMuscleForces();
        updateTorques();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateState_HYBRID(const vector<unsigned>& selectedmusclesIndex) {

        updateActivations(selectedmusclesIndex);
        updateFibreLengthsAndVelocities_HYBRID(selectedmusclesIndex);
        updateMuscleForces(selectedmusclesIndex);
        updateTorques();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateState_END_OF_HYBRID_MINIMIZATION() {

        for (auto& mIt : muscles_)
            mIt.updateFibreLengthsAndVelocities_END_OF_HYBRID_MINIMIZATION();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateActivations() {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt; muscleIt != muscles_.end(); ++muscleIt)
            muscleIt->updateActivation();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateActivations(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->updateActivation();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengthsAndVelocities() {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt; muscleIt != muscles_.end(); ++muscleIt)
            muscleIt->updateFibreLengthAndVelocity();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengthsAndVelocities(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->updateFibreLengthAndVelocity();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengthsAndVelocities_OFFLINE() {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt; muscleIt != muscles_.end(); ++muscleIt)
            muscleIt->updateFibreLengthAndVelocity_OFFLINE();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengthsAndVelocities_OFFLINE(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->updateFibreLengthAndVelocity_OFFLINE();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengthsAndVelocities_HYBRID() {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt; muscleIt != muscles_.end(); ++muscleIt)
            muscleIt->updateFibreLengthAndVelocity_HYBRID();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengthsAndVelocities_HYBRID(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->updateFibreLengthAndVelocity_HYBRID();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengths_OFFLINEPREP() {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt; muscleIt != muscles_.end(); ++muscleIt)
            muscleIt->updateFibreLength_OFFLINEPREP();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengths_OFFLINEPREP(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->updateFibreLength_OFFLINEPREP();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateMuscleForces() {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt; muscleIt != muscles_.end(); ++muscleIt)
            muscleIt->updateMuscleForce();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateMuscleForces(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->updateMuscleForce();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateTorques() {

        vectorDoFitr dofIt = dofs_.begin();
        for (dofIt; dofIt != dofs_.end(); ++dofIt)
            dofIt->updateTorque();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::pushState() {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt; muscleIt != muscles_.end(); ++muscleIt)
            muscleIt->pushState();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::pushState(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->pushState();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setMuscleForces(const std::vector<double>& muscleForces) {

        vector<double>::const_iterator mtfIt;
        vectorMTUitr muscleIt = muscles_.begin();
        for (mtfIt = muscleForces.begin(); mtfIt < muscleForces.end(); ++mtfIt) {
            muscleIt->setMuscleForce(*mtfIt);
            ++muscleIt;
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setActivations(const std::vector<double>& activations) {

        vector<double>::const_iterator actIt;
        vectorMTUitr muscleIt = muscles_.begin();
        for (actIt = activations.begin(); actIt < activations.end(); ++actIt) {
            muscleIt->setActivation(*actIt);
            ++muscleIt;
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengthTraces() {

        // usata
        vectorMTUitr muscleIt;
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            muscleIt->updateFibreLengthTrace();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::updateFibreLengthTraces(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->updateFibreLengthTrace();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getActivations(vector<double>& activations) const {

        activations.clear();
        activations.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            activations.emplace_back(muscleIt->getActivation());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getEmgs(vector<double>& currentEMGData) const {

        currentEMGData.clear();
        currentEMGData.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt != muscles_.end(); ++muscleIt)
            currentEMGData.emplace_back(muscleIt->getEmg());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getPastEmgs(vector<double>& pastEMGData) const {

        pastEMGData.clear();
        pastEMGData.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt != muscles_.end(); ++muscleIt)
            pastEMGData.emplace_back(muscleIt->getPastEmg());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getNeuralActivations(std::vector<double>& neuralActivations) const {

        neuralActivations.clear();
        neuralActivations.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt != muscles_.end(); ++muscleIt)
            neuralActivations.emplace_back(muscleIt->getNeuralActivation());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getFiberLengths(vector<double>& fibreLengths) const {

        fibreLengths.clear();
        fibreLengths.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            fibreLengths.emplace_back(muscleIt->getFiberLength());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getNormFiberLengths(std::vector<double>& normFiberLengths) const {

        normFiberLengths.clear();
        normFiberLengths.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            normFiberLengths.emplace_back(muscleIt->getNormFiberLength());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getFiberVelocities(vector<double>& fiberVelocities) const {

        fiberVelocities.clear();
        fiberVelocities.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            fiberVelocities.emplace_back(muscleIt->getFiberVelocity());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getNormFiberVelocities(vector<double>& normFiberVelocities) const {

        normFiberVelocities.clear();
        normFiberVelocities.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            normFiberVelocities.emplace_back(muscleIt->getNormFiberVelocity());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getPennationAnglesAtT(vector<double>& pennationAnglesAtT) const {

        pennationAnglesAtT.clear();
        pennationAnglesAtT.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            pennationAnglesAtT.emplace_back(muscleIt->getPennationAngleAtT());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMuscleForces(vector<double>& muscleForces) const {

        muscleForces.clear();
        muscleForces.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt != muscles_.end(); ++muscleIt)
            muscleForces.emplace_back(muscleIt->getMuscleForce());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getTorques(vector<double>& torques) const {

        torques.clear();
        torques.reserve(dofs_.size());
        vectorDoFconstItr dofIt = dofs_.begin();
        for (dofIt = dofs_.begin(); dofIt < dofs_.end(); ++dofIt)
            torques.emplace_back(dofIt->getTorque());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMusclesContribution(std::vector<double>& musclesContribution, unsigned whichDof) const {

        musclesContribution = dofs_.at(whichDof).getMusclesContribution();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getStrengthCoefficients(vector<double>& strengthCoefficients) const {

        strengthCoefficients.clear();
        strengthCoefficients.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            strengthCoefficients.emplace_back(muscleIt->getStrengthCoefficient());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getTendonSlackLengths(vector<double>& tendonSlackLengths) const {

        tendonSlackLengths.clear();
        tendonSlackLengths.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            tendonSlackLengths.emplace_back(muscleIt->getTendonSlackLength());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getOptimalFibreLengths(vector<double>& optimalFibreLengths) const {

        optimalFibreLengths.clear();
        optimalFibreLengths.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            optimalFibreLengths.emplace_back(muscleIt->getOptimalFibreLength());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getShapeFactors(vector<double>& shapeFactors) const {

        shapeFactors.clear();
        shapeFactors.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            shapeFactors.emplace_back(muscleIt->getShapeFactor());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getActivationScales(vector<double>& activationScales) const {

        activationScales.clear();
        activationScales.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            activationScales.emplace_back(muscleIt->getActivationScale());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getC1Coefficients(vector<double>& c1Coefficients) const {

        c1Coefficients.clear();
        c1Coefficients.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            c1Coefficients.emplace_back(muscleIt->getC1());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getC2Coefficients(vector<double>& c2Coefficients) const {

        c2Coefficients.clear();
        c2Coefficients.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            c2Coefficients.emplace_back(muscleIt->getC2());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMaxContractionVelocities(vector<double>& maxContractionVelocities) const {

        maxContractionVelocities.clear();
        maxContractionVelocities.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt != muscles_.end(); ++muscleIt)
            maxContractionVelocities.emplace_back(muscleIt->getMaxContractionVelocity());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getDampings(vector<double>& dampings) const {

        dampings.clear();
        dampings.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt != muscles_.end(); ++muscleIt)
            dampings.emplace_back(muscleIt->getDamping());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMusclesIndexOnDof(vector<unsigned>& musclesIndex, unsigned whichDof) const {

        musclesIndex.clear();
        musclesIndex.reserve(muscles_.size());
        vector<string> muscleNamesOnDof;
        dofs_.at(whichDof).getMusclesNamesOnDof(muscleNamesOnDof);
        vector<string>::const_iterator mnIt = muscleNamesOnDof.begin(), pos;
        for (mnIt; mnIt != muscleNamesOnDof.end(); ++mnIt) {
            pos = std::find(muscleNames_.begin(), muscleNames_.end(), *mnIt);
            musclesIndex.emplace_back(std::distance(muscleNames_.begin(), pos));
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMuscleNamesOnDofs(vector<vector<string> >& muscleNamesOnDofs) const {

        muscleNamesOnDofs.clear();
        for (vectorDoFconstItr it(dofs_.begin()); it != dofs_.end(); ++it) {
            vector<string> muscleNamesOnDof;
            it->getMusclesNamesOnDof(muscleNamesOnDof);
            muscleNamesOnDofs.emplace_back(muscleNamesOnDof);
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMuscleNamesOnDofs(map <string, vector<string> >& muscleNamesOnDofs) const {

        muscleNamesOnDofs.clear();
        for (vectorDoFconstItr it(dofs_.begin()); it != dofs_.end(); ++it) {
            vector<string> muscleNamesOnDof;
            it->getMusclesNamesOnDof(muscleNamesOnDof);
            muscleNamesOnDofs[it->getName()] = muscleNamesOnDof;
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMusclesIndexFromLastDof(vector<unsigned int>& musclesIndexList,
        const vector<string>& whichDofs) {

        vector<string> lastDof;
        lastDof.push_back(whichDofs.back());
        getMusclesIndexFromDofs(musclesIndexList, lastDof);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMusclesIndexFromDofs(vector<unsigned>& musclesIndexList,
        const vector<string>& whichDofs) {

        musclesIndexList.clear();
        vector<string> currentMusclesList;
        vector<string> musclesList;
        bool found;

        for (unsigned i = 0; i < whichDofs.size(); ++i) {
            found = false;
            for (unsigned j = 0; j < dofs_.size(); ++j) {
                if (whichDofs.at(i) == dofs_.at(j).getName()) {
                    dofs_.at(j).getMusclesNamesOnDof(currentMusclesList);
                    for (unsigned k = 0; k < currentMusclesList.size(); ++k) {
                        musclesList.push_back(currentMusclesList.at(k));
                    }
                    found = true;
                }
            }
            if (!found) {
                cout << whichDofs.at(i) << " not found on current subject, impossible to calibrate\n";
                exit(EXIT_FAILURE);
            }
        }
        std::sort(musclesList.begin(), musclesList.end()); //sort all muscles names in the vector
        std::vector<std::string>::iterator last = std::unique(musclesList.begin(), musclesList.end()); //removes duplicated muscles

        for (unsigned i = 0; i < muscles_.size(); ++i) {
            bool muscFound = false;
            std::vector<std::string>::iterator musclesOfInterest(musclesList.begin());
            while (musclesOfInterest != last && !muscFound)
            {
                if (muscles_.at(i).compareMusclesId(*musclesOfInterest)) {
                    musclesIndexList.push_back(i);
                    muscFound = true;
                }
                musclesOfInterest++;
            }
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMusclesIndexFromMusclesList(std::vector<unsigned>& muscleIndexList, const vector<string>& musclesList) {

        //  cout << "you should test this function before using it\n";
        muscleIndexList.clear();
        muscleIndexList.reserve(musclesList.size());
        vectorMTUconstItr foundIt;
        unsigned j;

        for (unsigned i = 0; i < musclesList.size(); ++i) {
            foundIt = muscles_.begin(); j = 0;
            while ((foundIt != muscles_.end()) && (!foundIt->compareMusclesId(musclesList.at(i)))) {
                foundIt++; ++j;
            }
            if (foundIt != muscles_.end())
                muscleIndexList.push_back(j);
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getGroupMusclesBasedOnStrengthCoefficients(vector<double>& values,
        vector< vector<int> >& muscleGroups) {

        values.clear();
        muscleGroups.clear();

        unsigned int currentMuscle = 0;
        for (currentMuscle = 0; currentMuscle < muscles_.size(); ++currentMuscle) {
            double currentStrengthCoefficient = muscles_.at(currentMuscle).getStrengthCoefficient();
            bool found = false;
            vector<double>::const_iterator vIt;
            vector< vector< int > >::iterator mgIt;
            for (vIt = values.begin(), mgIt = muscleGroups.begin(); vIt < values.end(); ++vIt, ++mgIt)
                if (*vIt == currentStrengthCoefficient) {
                    found = true;
                    mgIt->push_back(currentMuscle);
                }

            if (!found) {
                values.push_back(currentStrengthCoefficient);
                muscleGroups.push_back(vector<int>());
                muscleGroups.at(muscleGroups.size() - 1).push_back(currentMuscle);
            }
        }
    }

    //TODO: brrr.. questa e' da sistemare
    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getGroupMusclesBasedOnStrengthCoefficientsFilteredByMusclesIndexList(vector<double>& values,
        vector< vector<int> >& muscleGroups,
        const vector<unsigned>& musclesIndexList) {
        getGroupMusclesBasedOnStrengthCoefficients(values, muscleGroups);
        vector< vector <int> > refinedMuscleGroups;
        vector<double> refinedValues;
        vector<bool> checkedMuscle(musclesIndexList.size(), false);

        for (int i = 0; i < musclesIndexList.size(); ++i) {
            if (!checkedMuscle.at(i)) {
                for (int mgIdx = 0; mgIdx < muscleGroups.size(); ++mgIdx) {
                    vector<int>::iterator it;
                    it = std::find(muscleGroups.at(mgIdx).begin(), muscleGroups.at(mgIdx).end(), musclesIndexList.at(i));
                    if (it != muscleGroups.at(mgIdx).end()) {
                        refinedMuscleGroups.push_back(muscleGroups.at(mgIdx));
                        refinedValues.push_back(values.at(mgIdx));
                        //then check as used the other muscles in the group
                        for (int muscleIdx = 0; muscleIdx < muscleGroups.at(mgIdx).size(); ++muscleIdx) {
                            for (int j = 0; j < musclesIndexList.size(); ++j)
                                if (muscleGroups.at(mgIdx).at(muscleIdx) == musclesIndexList.at(j))
                                    checkedMuscle.at(j) = true;
                        }
                    }
                }
            }
        }
        //copy refinedMuscleGroups to muscleGroups
        muscleGroups.clear();
        muscleGroups = refinedMuscleGroups;
        values.clear();
        values = refinedValues;

#ifdef DEBUG
        cout << "Muscle Groups:\n";
        for(int i = 0; i < muscleGroups.size(); ++i)
        {
            for(int j = 0; j < muscleGroups.at(i).size(); ++j)
                cout << muscleGroups.at(i).at(j) << " ";
            cout << endl;
        }
#endif
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMuscle(MTUtype& muscle,
        const string& muscleName) {

        for (int currentMuscle = 0; currentMuscle < muscles_.size(); ++currentMuscle) {
            if (muscles_.at(currentMuscle).compareMusclesId(muscleName)) {
                muscle = muscles_.at(currentMuscle);
                return;
            }
        }

        cout << "NMSmodel::getMuscle: did not found: " << muscleName << endl;
        exit(EXIT_FAILURE);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setStrengthCoefficientsBasedOnGroups(const vector<double>& values,
        const vector< vector<int> >& muscleGroups) {

        for (unsigned int i = 0; i < values.size(); ++i)
            for (unsigned int j = 0; j < muscleGroups.at(i).size(); ++j)
                muscles_.at(muscleGroups.at(i).at(j)).setStrengthCoefficient(values.at(i));
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setStrengthCoefficients(const vector<double>& strengthCoefficients){

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator strengthCoefficientsIt = strengthCoefficients.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt, ++strengthCoefficientsIt)
            muscleIt->setStrengthCoefficient(*strengthCoefficientsIt);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setShapeFactor(double shapeFactor) {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            muscleIt->setShapeFactor(shapeFactor);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setShapeFactors(const vector<double>& shapeFactors){

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator shapeFactorsIt = shapeFactors.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt, ++shapeFactorsIt)
            muscleIt->setShapeFactor(*shapeFactorsIt);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setActivationScales(const vector<double>& activationScales) {

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator activationScalesIt = activationScales.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt, ++activationScalesIt)
            muscleIt->setActivationScale(*activationScalesIt);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setC1(double C1) {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            muscleIt->setC1(C1);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setC1Coefficients(const vector<double>& c1Coefficients){

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator c1CoefficientsIt = c1Coefficients.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt, ++c1CoefficientsIt)
            muscleIt->setC1(*c1CoefficientsIt);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setC2(double C2){

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            muscleIt->setC2(C2);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setC2Coefficients(const vector<double>& c2Coefficients){

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator c2CoefficientsIt = c2Coefficients.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt, ++c2CoefficientsIt)
            muscleIt->setC2(*c2CoefficientsIt);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setTendonSlackLengths(const vector<double>& tendonSlackLengths){

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator tendonSlackLengthIt = tendonSlackLengths.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt, ++tendonSlackLengthIt)
            muscleIt->setTendonSlackLength(*tendonSlackLengthIt);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setOptimalFibreLengths(const vector<double>& optimalFibreLengths) {

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator optimalFibreLengthIt = optimalFibreLengths.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt, ++optimalFibreLengthIt)
            muscleIt->setOptimalFibreLength(*optimalFibreLengthIt);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setMaxContractionVelocities(const vector<double>& maxContractionVelocities) {

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator it = maxContractionVelocities.begin();
        for (muscleIt = muscles_.begin(); muscleIt != muscles_.end(); ++muscleIt, ++it)
            muscleIt->setMaxContractionVelocity(*it);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setDampings(const vector<double>& dampings) {

        vectorMTUitr muscleIt = muscles_.begin();
        vector<double>::const_iterator it = dampings.begin();
        for (muscleIt = muscles_.begin(); muscleIt != muscles_.end(); ++muscleIt, ++it)
            muscleIt->setDamping(*it);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setTendonTolerance(double tolerance) {
        for (auto &muscleIt : muscles_)
            muscleIt.setTendonTolerance(tolerance);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::resetFibreLengthTraces() {

        vectorMTUitr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            muscleIt->resetFibreLengthTrace();
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::resetFibreLengthTraces(const vector<unsigned>& selectedMusclesIndex) {

        vectorMTUitr muscleIt;
        vector<unsigned>::const_iterator it = selectedMusclesIndex.begin();
        for (it; it != selectedMusclesIndex.end(); ++it) {
            muscleIt = muscles_.begin() + *it;
            muscleIt->resetFibreLengthTrace();
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMomentArmsOnDof(vector<double>& momentArms, unsigned whichDof) const {

        dofs_.at(whichDof).getMomentArms(momentArms);
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    double NMSmodel<Activation, Tendon, mode>::getMusclesPenalty() const {

        double penalty = 0.;
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            penalty += muscleIt->getPenalty();
        return penalty;
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    double NMSmodel<Activation, Tendon, mode>::getMusclesPenalty(vector<unsigned>& selectedMusclesIndex) const {

        double penalty = 0.;
        unsigned int mILi = 0;
        for (unsigned int i = 0; i < muscles_.size() && mILi < selectedMusclesIndex.size(); ++i)
            if (i == selectedMusclesIndex.at(mILi)) {
                penalty += muscles_.at(i).getPenalty();
                ++mILi;
            }
        return penalty;
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMusclesPenaltyVector(vector<double>& penalties) const {

        penalties.clear();
        penalties.reserve(muscles_.size());
        vectorMTUconstItr muscleIt = muscles_.begin();
        for (muscleIt = muscles_.begin(); muscleIt < muscles_.end(); ++muscleIt)
            penalties.emplace_back(muscleIt->getPenalty());
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    double NMSmodel<Activation, Tendon, mode>::getGlobalEmDelay() const {

        return muscles_.begin()->getEmDelay();
    }


    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setOpenSimModelFilename(const std::string& openSimModelFilename) {

        osimModelFilename_ = openSimModelFilename;
        hasOpenSimModel_ = true;

    }


    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    std::string NMSmodel<Activation, Tendon, mode>::getOpenSimModelFilename() const {

        return osimModelFilename_;
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::getMusclesParameters(vector<MuscleParameters>& parameters) const {
        //this stuff with muscle parameters should be a call to the mtu instead..
        parameters.clear();
        parameters.resize(muscles_.size());

        for (unsigned int i = 0; i < muscles_.size(); ++i) {
            parameters.at(i).setC1(muscles_.at(i).getC1());
            parameters.at(i).setC2(muscles_.at(i).getC2());
            parameters.at(i).setShapeFactor(muscles_.at(i).getShapeFactor());
            parameters.at(i).setActivationScale(muscles_.at(i).getActivationScale());
            parameters.at(i).setOptimalFiberLength(muscles_.at(i).getOptimalFibreLength());
            parameters.at(i).setPennationAngle(muscles_.at(i).getPennationAngle());
            parameters.at(i).setTendonSlackLength(muscles_.at(i).getTendonSlackLength());
            parameters.at(i).setMaxIsometricForce(muscles_.at(i).getMaxIsometricForce());
            parameters.at(i).setStrengthCoefficient(muscles_.at(i).getStrengthCoefficient());
            parameters.at(i).setEmDelay(muscles_.at(i).getEmDelay());
            parameters.at(i).setMaxContractionVelocity(muscles_.at(i).getMaxContractionVelocity());
            parameters.at(i).setDamping(muscles_.at(i).getDamping());
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    void NMSmodel<Activation, Tendon, mode>::setMusclesParameters(const vector<MuscleParameters>& parameters)  {
        //carefull with this.. there's no control over the order of muscles..

        auto mIt(muscles_.begin());
        auto pIt(parameters.cbegin());
        for (pIt; pIt != parameters.cend(); ++pIt, ++mIt){
            mIt->setC1(pIt->getC1());
            mIt->setC2(pIt->getC2());
            mIt->setShapeFactor(pIt->getShapeFactor());
            mIt->setActivationScale(pIt->getActivationScale());
            mIt->setOptimalFibreLength(pIt->getOptimalFiberLength());
            mIt->setPennationAngle(pIt->getPennationAngle());
            mIt->setTendonSlackLength(pIt->getTendonSlackLength());
            mIt->setMaxIsometricForce(pIt->getMaxIsometricForce());
            mIt->setStrengthCoefficient(pIt->getStrengthCoefficient());
            mIt->setEmDelay(pIt->getEmDelay());
            mIt->setMaxContractionVelocity(pIt->getMaxContractionVelocity());
            mIt->setDamping(pIt->getDamping());
        }
    }

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    std::ostream& operator<< (std::ostream& output, const NMSmodel<Activation, Tendon, mode>& b) {

        output << "Current NMSmodel has " << b.muscles_.size() << " muscles:\n";
        for (typename std::vector< MTU<Activation, Tendon, mode> >::const_iterator m = b.muscles_.begin();
            m != b.muscles_.end(); ++m) {
            output << "--------- MUSCLE -----------\n";
            output << *m << std::endl;
        }
        output << "and " << b.dofs_.size() << " DoF:\n";
        for (typename vector< DoF<Activation, Tendon, mode> >::const_iterator dof = b.dofs_.begin();
            dof != b.dofs_.end(); ++dof) {
            output << "--------- DOF -----------\n";
            output << *dof << std::endl;
        }
        return output;
    }
}

#include "ceinms/NMSmodelPolicyTemplates.h"
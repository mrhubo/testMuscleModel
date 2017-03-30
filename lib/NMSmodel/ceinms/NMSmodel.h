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

#ifndef ceinms_NMSmodel_h
#define ceinms_NMSmodel_h

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <stdlib.h>
#include "MTU.h"
#include "DoF.h"
#include "MuscleParameters.h"

namespace ceinms {
    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    class NMSmodel;

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    std::ostream& operator<< (std::ostream& output, const NMSmodel<Activation, Tendon, mode>& b);

    template <typename Activation, typename Tendon, CurveMode::Mode mode>
    class NMSmodel
    {
    public:

        typedef Activation ActivationT;
        typedef Tendon TendonT;
        typedef MTU<Activation, Tendon, mode> MTUtype;
        typedef DoF<Activation, Tendon, mode> DoFtype;
        typedef typename std::vector<MTUtype>::iterator vectorMTUitr;
        typedef typename std::vector<DoFtype>::iterator vectorDoFitr;
        typedef typename std::vector<MTUtype>::const_iterator vectorMTUconstItr;
        typedef typename std::vector<DoFtype>::const_iterator vectorDoFconstItr;

        NMSmodel();
        NMSmodel(const NMSmodel& orig);// { std::cout << "YOU SHOULD NOT BE HERE!"; exit(EXIT_FAILURE); }
        virtual ~NMSmodel() {}
        bool haveThisMuscle(const std::string& currentMuscle, vectorMTUitr& found);
        bool haveTheseMuscles(const std::vector<std::string>& musclesNames);
        bool compareMusclesNames(const std::vector<std::string>& muscleNames) const;
        bool compareMusclesNamesOnDoF(const std::vector<std::string>& muscleNames, int whichDoF) const;

        void setTime(const double& time);
        void setEmgs(const std::vector<double>& currentEMGData);
        void setEmgsSelective(const std::vector<double>& currentEMGData, const std::vector<unsigned>& selectedMusclesIndex);
        void setTime_emgs_updateActivations_pushState_selective(double time, const std::vector<double>& currentEMGData, const std::vector<unsigned>& selectedMusclesIndex);
        void setMuscleTendonLengths(const std::vector<double>& currentLmtData);
        void setMuscleTendonLengthsSelective(const std::vector<double>& currentLmtData, const std::vector<unsigned>& selectedMusclesIndex);
        void setMomentArms(const std::vector<double>& currentMaData, unsigned whichDof);
        void updateState();
        void updateState(const std::vector<unsigned>& selectedMusclesIndex);
        void updateState_OFFLINE();
        void updateState_OFFLINE(const std::vector<unsigned>& selectedMusclesIndex);
        void updateState_HYBRID();
        void updateState_HYBRID(const std::vector<unsigned>& selectedMusclesIndex);
        void updateState_END_OF_HYBRID_MINIMIZATION();
        void updateActivations();
        void updateActivations(const std::vector<unsigned>& selectedMusclesIndex);
        void updateFibreLengthsAndVelocities();
        void updateFibreLengthsAndVelocities(const std::vector<unsigned>& selectedMusclesIndex);
        void updateFibreLengthsAndVelocities_OFFLINE();
        void updateFibreLengthsAndVelocities_OFFLINE(const std::vector<unsigned>& selectedMusclesIndex);
        void updateFibreLengthsAndVelocities_HYBRID();
        void updateFibreLengthsAndVelocities_HYBRID(const std::vector<unsigned>& selectedMusclesIndex);
        void updateFibreLengths_OFFLINEPREP();
        void updateFibreLengths_OFFLINEPREP(const std::vector<unsigned>& selectedMusclesIndex);
        void updateMuscleForces();
        void updateMuscleForces(const std::vector<unsigned>& selectedMusclesIndex);
        void updateTorques();
        void pushState();
        void pushState(const std::vector<unsigned>& selectedMusclesIndex);
        void updateFibreLengthTraces();
        void updateFibreLengthTraces(const std::vector<unsigned>& selectedMusclesIndex);
        void resetFibreLengthTraces();
        void resetFibreLengthTraces(const std::vector<unsigned>& selectedMusclesIndex);

        // called on model data initialization
        void addMuscle(const MTUtype& muscle);
        void addDoF(const DoFtype& dof);
        void setOpenSimModelFilename(const std::string& openSimModelFilename);
        std::string getOpenSimModelFilename() const;

        // called on model data initialization or in model calibration
        void setStrengthCoefficientsBasedOnGroups(
            const std::vector<double>& values,
            const std::vector< std::vector< int > >& muscleGroups);
        void setStrengthCoefficients(const std::vector<double>& strengthCoefficients);
        void setShapeFactor(double shapeFactor);
        void setShapeFactors(const std::vector<double>& shapeFactors);
        void setActivationScales(const std::vector<double>& activationScales);
        void setC1(double C1);
        void setC1Coefficients(const std::vector<double>& c1Coefficients);
        void setC2(double C2);
        void setC2Coefficients(const std::vector<double>& c2Coefficients);
        void setTendonSlackLengths(const std::vector<double>& tendonSlackLengths);
        void setActivations(const std::vector<double>& activations);
        void setOptimalFibreLengths(const std::vector<double>&  optimalFibreLengths);
        void setTendonTolerance(double tolerance = 0.0);
        void setMaxContractionVelocities(const std::vector<double>& maxContractionVelocities);
        void setDampings(const vector<double>& dampings);
        // called	for control purpose
        std::vector<std::string> getMuscleNames() const { return muscleNames_; }
        void getMuscleNames(std::vector<std::string>& muscleNames) const { muscleNames = muscleNames_; }
        void getMuscleNamesOnDofs(std::vector<std::vector<std::string> >& muscleNamesOnDofs) const;
        void getMuscleNamesOnDofs(std::map< std::string, std::vector<std::string> >& muscleNamesOnDofs) const;
        int  getNoMuscles() const { return muscles_.size(); }
        void getDoFNames(std::vector<std::string>& dofNames) const { dofNames = dofNames_; }
        std::vector<std::string> getDoFNames() const { return dofNames_; };
        int  getNoDoF() const { return dofs_.size(); }
        bool hasOpenSimModel() { return hasOpenSimModel_; }
        // called for get computed data
        void getEmgs(std::vector<double>& currentEMGData) const;
        void getPastEmgs(std::vector<double>& pastEMGData) const;
        void getActivations(std::vector<double>& activations) const;
        void getNeuralActivations(std::vector<double>& neuralActivations) const;
        void getFiberLengths(std::vector<double>& fiberLengths) const;
        void getNormFiberLengths(std::vector<double>& fiberLengths) const;
        void getFiberVelocities(std::vector<double>& fiberVelocities) const;
        void getNormFiberVelocities(std::vector<double>& normFiberVelocities) const;
        void getPennationAnglesAtT(std::vector<double>& pennationAnglesAtT) const; //instant pennationangle
        void getMuscleForces(std::vector<double>& muscleForces) const;
        void getTorques(std::vector<double>& torques) const;
        void getMusclesContribution(std::vector<double>& musclesContribution, unsigned whichDof) const;

        // used during calibration
        void getGroupMusclesBasedOnStrengthCoefficients(std::vector<double>& values,
            std::vector< std::vector< int > >& muscleGroups);
        void getGroupMusclesBasedOnStrengthCoefficientsFilteredByMusclesIndexList(std::vector<double>& values,
            std::vector< std::vector< int > >& muscleGroups, const std::vector<unsigned int>& musclesIndexList);
        void getMusclesIndexOnDof(std::vector<unsigned>& muscleIndex, unsigned whichDof) const;
        void getMusclesIndexFromDofs(std::vector<unsigned int>& musclesIndexList, const std::vector<std::string>& whichDofs);
        void getMusclesIndexFromLastDof(std::vector<unsigned int>& musclesIndexList, const std::vector<std::string>& whichDofs);
        void getMusclesIndexFromMusclesList(std::vector<unsigned>& muscleIndexList, const std::vector<std::string>& musclesList);
        double getShapeFactor() const { return muscles_.at(0).getShapeFactor(); }
        void getShapeFactors(std::vector<double>& shapeFactors) const;
        void getActivationScales(std::vector<double>& activationScales) const;
        double getC1() const { return muscles_.at(0).getC1(); }
        void getC1Coefficients(std::vector<double>& c1Coefficients) const;
        double getC2() const { return muscles_.at(0).getC2(); }
        void getC2Coefficients(std::vector<double>& c2Coefficients) const;
        void getStrengthCoefficients(std::vector<double>& strengthCoefficients) const;
        void getTendonSlackLengths(std::vector<double>& tendonSlackLengths) const;
        void getMuscle(MTUtype& muscle, const std::string& muscleName);
        double getMusclesPenalty() const;
        double getMusclesPenalty(std::vector<unsigned>& selectedMusclesIndex) const;
        void getMusclesPenaltyVector(std::vector<double>& penalties) const;
        void getMusclesParameters(std::vector< MuscleParameters >& parameters) const;
        void setMusclesParameters(const std::vector< MuscleParameters >& parameters);
        void setMuscleForces(const std::vector<double>& muscleTendonForces);
        void getMomentArmsOnDof(std::vector<double>& momentArms, unsigned whichDof) const;
        void getOptimalFibreLengths(std::vector<double>& optimalFibreLengths) const;
        void getMaxContractionVelocities(std::vector<double>& maxContractionVelocities) const;
        void getDampings(std::vector<double>& dampings) const;

        double getGlobalEmDelay() const; //NOTE: this will be replaced by singular muscles emdelay

        friend std::ostream& operator<< <> (std::ostream& output, const NMSmodel& b);

    private:
        std::vector<MTUtype> muscles_;
        std::vector<DoFtype> dofs_;
        std::vector<std::string> muscleNames_;
        std::vector<std::string> dofNames_;
        std::string osimModelFilename_;
        bool hasOpenSimModel_;
    };
}

#endif

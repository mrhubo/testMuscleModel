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
 * Author(s): C. Pizzolato, M. Reggiani                                       *
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

#include "ceinms/Curve.h"
#include <iostream>
using std::cout;
using std::endl;
#include <stdlib.h>
#include <cmath>
using std::cos;
#include "MTUutils.h"

namespace ceinms {
    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    MTU<Activation, Tendon, mode>::MTU()
        :id_(""), emDelay_(.0), c1_(0.), c2_(0.), shapeFactor_(0.), activation_(0.),
        optimalFibreLength_(0.), pennationAngle_(0.), tendonSlackLength_(0.),
        percentageChange_(0.), fibreLength_(0.), pennationAngleAtT_(0.),
        fibreVelocity_(0.), damping_(0.), maxIsometricForce_(0.), time_(.0), timeScale_(0.),
        strengthCoefficient_(0.), muscleForce_(0.), maxContractionVelocity_(0.), activationScale_(1.)  { }

    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    MTU<Activation, Tendon, mode>::MTU(std::string id)
        :id_(id), emDelay_(.0), c1_(0.), c2_(0.), shapeFactor_(0.), activation_(0.),
        optimalFibreLength_(0.), pennationAngle_(0.), tendonSlackLength_(0.),
        percentageChange_(0.), fibreLength_(0.), pennationAngleAtT_(0.),
        fibreVelocity_(0.), damping_(0.), maxIsometricForce_(0.), time_(.0), timeScale_(0.),
        strengthCoefficient_(0.), muscleForce_(0.), maxContractionVelocity_(0.), activationScale_(1.), tendonDynamic_(id)  { }

    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    MTU<Activation, Tendon, mode>::MTU(const MTU<Activation, Tendon, mode>& orig) {

        id_ = orig.id_;
        emDelay_ = orig.emDelay_;
        activationDynamic_ = orig.activationDynamic_;
        activation_ = orig.activation_;
        c1_ = orig.c1_;
        c2_ = orig.c2_;
        shapeFactor_ = orig.shapeFactor_;
        activationScale_ = orig.activationScale_;

        tendonDynamic_ = orig.tendonDynamic_;
        fibreVelocity_ = orig.fibreVelocity_;
        normFibreVelocity_ = orig.normFibreVelocity_;
        fibreLength_ = orig.fibreLength_;
        fibreLengthTrace_ = orig.fibreLengthTrace_;
        muscleForce_ = orig.muscleForce_;
        pennationAngleAtT_ = orig.pennationAngleAtT_;

        optimalFibreLength_ = orig.optimalFibreLength_;
        pennationAngle_ = orig.pennationAngle_;
        tendonSlackLength_ = orig.tendonSlackLength_;
        percentageChange_ = orig.percentageChange_;
        damping_ = orig.damping_;
        maxIsometricForce_ = orig.maxIsometricForce_;
        strengthCoefficient_ = orig.strengthCoefficient_;
        maxContractionVelocity_ = orig.maxContractionVelocity_;
        forceVelocityCurve_ = orig.forceVelocityCurve_;
        activeForceLengthCurve_ = orig.activeForceLengthCurve_;
        passiveForceLengthCurve_ = orig.passiveForceLengthCurve_;
        tendonForceStrainCurve_ = orig.tendonForceStrainCurve_;

        timeScale_ = orig.timeScale_;
        time_ = orig.time_;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    MTU<Activation, Tendon, mode>& MTU<Activation, Tendon, mode>::operator=(const MTU<Activation, Tendon, mode>& orig) {

        id_ = orig.id_;
        emDelay_ = orig.emDelay_;
        activationDynamic_ = orig.activationDynamic_;
        activation_ = orig.activation_;
        c1_ = orig.c1_;
        c2_ = orig.c2_;
        shapeFactor_ = orig.shapeFactor_;
        activationScale_ = orig.activationScale_;

        tendonDynamic_ = orig.tendonDynamic_;
        fibreVelocity_ = orig.fibreVelocity_;
        normFibreVelocity_ = orig.normFibreVelocity_;
        fibreLength_ = orig.fibreLength_;
        fibreLengthTrace_ = orig.fibreLengthTrace_;
        muscleForce_ = orig.muscleForce_;
        pennationAngleAtT_ = orig.pennationAngleAtT_;

        optimalFibreLength_ = orig.optimalFibreLength_;
        pennationAngle_ = orig.pennationAngle_;
        tendonSlackLength_ = orig.tendonSlackLength_;
        percentageChange_ = orig.percentageChange_;
        damping_ = orig.damping_;
        maxIsometricForce_ = orig.maxIsometricForce_;
        strengthCoefficient_ = orig.strengthCoefficient_;
        maxContractionVelocity_ = orig.maxContractionVelocity_;
        forceVelocityCurve_ = orig.forceVelocityCurve_;
        activeForceLengthCurve_ = orig.activeForceLengthCurve_;
        passiveForceLengthCurve_ = orig.passiveForceLengthCurve_;
        tendonForceStrainCurve_ = orig.tendonForceStrainCurve_;

        timeScale_ = orig.timeScale_;
        time_ = orig.time_;
        return *this;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setParametersToComputeForces(double optimalFibreLength,
        double pennationAngle,
        double tendonSlackLength,
        double percentageChange,
        double damping,
        double maxIsometricForce,
        double strengthCoefficient,
        double maxContractionVelocity) {
        optimalFibreLength_ = optimalFibreLength;
        pennationAngle_ = pennationAngle;
        tendonSlackLength_ = tendonSlackLength;
        percentageChange_ = percentageChange;
        damping_ = damping;
        maxIsometricForce_ = maxIsometricForce;
        strengthCoefficient_ = strengthCoefficient;
        maxContractionVelocity_ = maxContractionVelocity;

        tendonDynamic_.setParametersToComputeForces(optimalFibreLength,
            pennationAngle,
            tendonSlackLength,
            percentageChange,
            damping,
            maxIsometricForce,
            strengthCoefficient,
            maxContractionVelocity);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setParametersToComputeActivation(double c1, double c2, double shapeFactor, double activationScale) {

        c1_ = c1;
        c2_ = c2;
        shapeFactor_ = shapeFactor;
        activationScale_ = activationScale;
        activationDynamic_.setFilterParameters(c1_, c2_);
        activationDynamic_.setShapeFactor(shapeFactor_);
        activationDynamic_.setActivationScale(activationScale_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::resetState() {

        activationDynamic_.resetState();
        tendonDynamic_.resetState();
        fibreLengthTrace_.reset();
        activation_ = .0;
        time_ = 0;
        timeScale_ = 0.005;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setTime(const double& time) {

        timeScale_ = time - time_;
        time_ = time;
        tendonDynamic_.setTime(time);
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setEmg(double emg) {

        if (!activationDynamic_.setEmg(emg))
            std::cout << "WARNING! Wrong excitation data provided for muscle " << getMuscleName() << ", it should be between 0 and 1" << std::endl;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setMuscleTendonLength(double muscleTendonLength) {

        tendonDynamic_.setMuscleTendonLength(muscleTendonLength);
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateActivation() {

        activationDynamic_.updateActivation();
        activation_ = activationDynamic_.getActivation();
        tendonDynamic_.setActivation(activation_);
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreLengthAndVelocity() {

        updateFibreLength();
        //add extra point to control the derivative
     //   fibreLengthTrace_.addPointNoUpdate(time_ + timeScale_, fibreLength_);
        updateFibreVelocity();
     //   fibreLengthTrace_.removeLastPointNoUpdate();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreLength() {

        tendonDynamic_.updateFibreLength();
        fibreLength_ = tendonDynamic_.getFibreLength();
  //      cout << fibreLengthTrace_.getNoElements() << std::endl;
        fibreLengthTrace_.addPointNoUpdate(time_, fibreLength_);
        //helps control the derivative
    //    cout << fibreLengthTrace_.getNoElements() << std::endl;

    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreVelocity() {

        fibreLengthTrace_.refresh();
        fibreVelocity_ = fibreLengthTrace_.getFirstDerivative(time_);
        cout << fibreVelocity_ << endl;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreLength_OFFLINEPREP() {

        updateFibreLength();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreLengthAndVelocity_OFFLINE() {

        updateFibreLength_OFFLINE();
        updateFibreVelocity_OFFLINE();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreLength_OFFLINE() {

        tendonDynamic_.updateFibreLength();
        fibreLength_ = tendonDynamic_.getFibreLength();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreVelocity_OFFLINE() {

        fibreVelocity_ = fibreLengthTrace_.getFirstDerivative(time_);
  //      std::cout << id_ << " fibre velocity: " << fibreVelocity_ << std::endl;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreLengthAndVelocity_HYBRID() {

        updateFibreLength_OFFLINE();
        fibreLengthTrace_.addPoint(time_, fibreLength_);
//        fibreLengthTrace_.refresh();
        updateFibreVelocity_OFFLINE();
        fibreLengthTrace_.removeLastPointNoUpdate();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreLengthsAndVelocities_END_OF_HYBRID_MINIMIZATION() {

        fibreLengthTrace_.removeLastPointNoUpdate();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateMuscleForce() {

        double optimalFiberLengthAtT = optimalFibreLength_*(percentageChange_*(1.0 - activation_) + 1);
        double normFiberLengthAtT = fibreLength_ / optimalFiberLengthAtT;
        double normFiberLength = fibreLength_ / optimalFibreLength_;
        double fiberVel = fibreVelocity_ / optimalFibreLength_;
        if (fiberVel > maxContractionVelocity_)
            fiberVel = maxContractionVelocity_;
        if (fiberVel < -maxContractionVelocity_)
            fiberVel = -maxContractionVelocity_;
        double normFiberVelocity = fiberVel / maxContractionVelocity_;

        double fv = forceVelocityCurve_.getValue(normFiberVelocity);
        double fp = passiveForceLengthCurve_.getValue(normFiberLength);
        double fa = activeForceLengthCurve_.getValue(normFiberLengthAtT);
        double pennationAngleAtT = computePennationAngle(optimalFibreLength_);
        normFibreVelocity_ = normFiberVelocity;
        pennationAngleAtT_ = pennationAngleAtT; //for logging
        muscleForce_ = maxIsometricForce_*strengthCoefficient_*
            (fa*fv*activation_ + fp + damping_*normFiberVelocity)*
            cos(pennationAngleAtT);
        //clamp muscle force
        if (muscleForce_ < 0)
            muscleForce_ = 0;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::pushState() {

        activationDynamic_.pushState();
        tendonDynamic_.pushState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::resetFibreLengthTrace() {

        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::updateFibreLengthTrace() {

        fibreLengthTrace_.refresh();
        activationDynamic_.resetState();
        tendonDynamic_.resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setCurves(const CurveOffline& activeForceLengthCurve,
        const CurveOffline& passiveForceLengthCurve,
        const CurveOffline& forceVelocityCurve,
        const CurveOffline& tendonForceStrainCurve) {

        activeForceLengthCurve_ = activeForceLengthCurve;
        passiveForceLengthCurve_ = passiveForceLengthCurve;
        forceVelocityCurve_ = forceVelocityCurve;
        tendonForceStrainCurve_ = tendonForceStrainCurve;

        tendonDynamic_.setCurves(activeForceLengthCurve, passiveForceLengthCurve, forceVelocityCurve, tendonForceStrainCurve);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setActivation(double activation) {

        activation_ = activation;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setShapeFactor(double shapeFactor) {

        shapeFactor_ = shapeFactor;
        activationDynamic_.setShapeFactor(shapeFactor_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setActivationScale(double activationScale) {

        activationScale_ = activationScale;
        activationDynamic_.setActivationScale(activationScale_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setC1(double c1) {

        c1_ = c1;
        activationDynamic_.setFilterParameters(c1_, c2_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setC2(double c2) {

        c2_ = c2;
        activationDynamic_.setFilterParameters(c1_, c2_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setStrengthCoefficient(double strengthCoefficient) {

        strengthCoefficient_ = strengthCoefficient;
        tendonDynamic_.setStrengthCoefficient(strengthCoefficient_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setTendonSlackLength(double tendonSlackLength) {

        tendonSlackLength_ = tendonSlackLength;
        tendonDynamic_.setTendonSlackLength(tendonSlackLength_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setOptimalFibreLength(double optimalFiberLength) {

        optimalFibreLength_ = optimalFiberLength;
        tendonDynamic_.setOptimalFibreLength(optimalFibreLength_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setPennationAngle(double pennationAngle) {

        pennationAngle_ = pennationAngle;
        tendonDynamic_.setPennationAngle(pennationAngle_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setMaxIsometricForce(double maxIsometricForce) {

        maxIsometricForce_ = maxIsometricForce;
        tendonDynamic_.setMaxIsometricForce(maxIsometricForce_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setMaxContractionVelocity(double maxContractionVelocity) {

        maxContractionVelocity_ = maxContractionVelocity;
        tendonDynamic_.setMaxContractionVelocity(maxContractionVelocity_);
        resetState();
    }

    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setDamping(double damping) {
        damping_ = damping;
        tendonDynamic_.setDamping(damping_);
        resetState();
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    double MTU<Activation, Tendon, mode>::getPenalty() const {

        double penalty(tendonDynamic_.getPenalty());
        double const normalisedFibreRatio(fabs(fibreLength_ / optimalFibreLength_ - 1.0));
        if (normalisedFibreRatio > 0.5)
            penalty += normalisedFibreRatio*normalisedFibreRatio * 100;
        return penalty;
    }


    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    inline double MTU<Activation, Tendon, mode>::computePennationAngle(double optimalFiberLength) {

        return PennationAngle::compute(fibreLength_, optimalFiberLength, pennationAngle_);

    }

    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setTendonTolerance(double tolerance) {

        tendonDynamic_.setTolerance(tolerance);
    }

    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    std::ostream& operator<< (std::ostream& output, const MTU<Activation, Tendon, mode>& m)
    {
        output << "Name: " << m.id_ << endl;
        output << "emDelay: " << m.emDelay_ << endl;
        output << "C1: " << m.c1_ << " C2: " << m.c2_ << endl;
        output << "Shape Factor: " << m.shapeFactor_ << endl;
        output << "Activation Scale: " << m.activationScale_ << endl;
        output << "activeForceLength" << endl << m.activeForceLengthCurve_ << endl;
        output << "passiveForceLength" << endl << m.passiveForceLengthCurve_ << endl;
        output << "forceVelocity" << endl << m.forceVelocityCurve_ << endl;
        output << "tendonForceStrain " << endl << m.tendonForceStrainCurve_ << endl;
        output << "optimalFibreLength: " << m.optimalFibreLength_ << endl;
        output << "pennationAngle: " << m.pennationAngle_ << endl;
        output << "tendonSlackLength: " << m.tendonSlackLength_ << endl;
        output << "percentageChange: " << m.percentageChange_ << endl;
        output << "damping: " << m.damping_ << endl;
        output << "maxIsometricForce: " << m.maxIsometricForce_ << endl;
        output << "strengthCoefficient: " << m.strengthCoefficient_ << endl;
        // :TODO: valli a mettere anche nel costruttore di copia

        return output;
    }
    
    // set %FT
    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    void MTU<Activation, Tendon, mode>::setPercentFastTwitch(double percentFastTwitch) {
        percentFastTwitch_ = percentFastTwitch;
		tendonDynamic_.setPercentFastTwitch(percentFastTwitch_);
		resetState();
    }
    
	// set muscle mass
	template<typename Activation, typename Tendon, CurveMode::Mode mode>
	void MTU<Activation, Tendon, mode>::setMuscleMass(double muscleMass) {
		muscleMass_ = muscleMass;
		tendonDynamic_.setMuscleMass(muscleMass_);
		resetState();
	}

    // compute hAM
    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    inline double MTU<Activation, Tendon, mode>::computeAct_MainHeatRate() {
		tendonDynamic_.computeAct_MainHeatRate();
		return act_mainHeatRate_ = 1.28 * percentFastTwitch_ + 25;		
    }
    
    // compute shortening heat coefficients
    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    inline double MTU<Activation, Tendon, mode>::computeShorteningHeatCoeFT() {
		tendonDynamic_.computeShorteningHeatCoeFT();
        return shorteningHeatCoeFT_ = 4 * 25 / maxContractionVelocity_;
    }
    
    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    inline double MTU<Activation, Tendon, mode>::computeShorteningHeatCoeST() {
		tendonDynamic_.computeShorteningHeatCoeST();
        return shorteningHeatCoeST_ = 4 * 25 / (maxContractionVelocity_ / 2.5);
    }
    
    // compute shortening heat rate
    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    inline double MTU<Activation, Tendon, mode>::computeShorteningHeatRate() {
		tendonDynamic_.computeShorteningHeatRate();
        shorteningHeatRate_ = - shorteningHeatCoeST_ * normFibreVelocity_ * (1 - percentFastTwitch_ / 100) 
			- shorteningHeatCoeFT_ * normFibreVelocity_ * (percentFastTwitch_ / 100);
        return shorteningHeatRate_;
    }

	// compute lengthening heat rate
	template<typename Activation, typename Tendon, CurveMode::Mode mode>
	inline double MTU<Activation, Tendon, mode>::computeLengtheningHeatRate() {
		tendonDynamic_.computeLengtheningHeatRate();
		lengtheningHeatRate_ = 4 * shorteningHeatCoeST_ * normFibreVelocity_;
		return lengtheningHeatRate_;
	}

    // compute mechanical work rate
    template<typename Activation, typename Tendon, CurveMode::Mode mode>
    inline double MTU<Activation, Tendon, mode>::computeWorkRate() {
		tendonDynamic_.computeWorkRate();
        return workRate_ = muscleForce_ * fibreVelocity_ / muscleMass_;
    }

	// set parameters to calculate heat rate
	template<typename Activation, typename Tendon, CurveMode::Mode mode>
	inline double MTU<Activation, Tendon, mode>::setParametersToComputeHeatRate(double percentFastTwitch,
		double maxContractionVelocity,
		double optimalFibreLength,
		double activation,
		double activationScale) {
		percentFastTwitch_ = percentFastTwitch;
		maxContractionVelocity_ = maxContractionVelocity;
		optimalFibreLength_ = optimalFibreLength;
		activation_ = activation;
		activationScale_ = activationScale;

		tendonDynamic_.setParametersToComputeHeatRate(percentFastTwitch,
			maxContractionVelocity,
			optimalFibreLength,
			activation,
			activationScale);
		resetState();
	}

	// compute heat rate ( shortening Lce < Lce(opt) & ~Vce < 0)
	template<typename Activation, typename Tendon, CurveMode::Mode mode>
	inline double MTU<Activation, Tendon, mode>::computeHeatRate(){
		// update hAM, hSL, workRate
		//
		if (fibreLength_ <= optimalFibreLength_)
		{
			if (normFibreVelocity_ <= 0)
			{
				heatRate_ = act_mainHeatRate_ * pow(activation_, 0.6) * activationScale_
					+ shorteningHeatRate_ * pow(activation_, 2.0) * activationScale_
					- workRate_;
			}
			else
			{
				heatRate_ = act_mainHeatRate_ * pow(activation_, 0.6) * activationScale_
					+ lengtheningHeatRate_ * normFiberVelocity_ * activation_ * activationScale_
					- workRate_;
			}
		else
		{
			if (normFibreVelocity_ <= 0)
			{
				//heatRate = (0.4 + 0.6 * imsometricForce)
			}
			else
			{

			}
		}
		}
		

		return heatRate_;
	}
}

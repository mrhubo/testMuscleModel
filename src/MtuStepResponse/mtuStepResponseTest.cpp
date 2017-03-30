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
 * Author(s): C. Pizzolato, L. Tagliapietra                                   *
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
using namespace ceinms;
#include <vector>
#include <fstream>
#include <string>
#include <cmath>
#include "ceinms/NMSmodel.h"
#include "ceinms/Activation/ExponentialActivation.h"
#include "ceinms/Tendon/ElasticTendon_BiSec.h"
#include "ceinms/Tendon/StiffTendon.h"

using namespace std;
using MtuCurve = Curve<CurveMode::Mode::Offline, CurveMode::Interpolation::Cubic>;


MtuCurve getActiveForceLengthCurve() {

    MtuCurve activeForceLength;
    std::vector<double> afX{ 0.4035, 0.52725, 0.62875, 0.71875, 0.86125, 1.045, 1.2175, 1.43875, 1.61875 };
    std::vector<double> afY{ 0, 0.226667, 0.636667, 0.856667, 0.95, 0.993333, 0.77, 0.246667, 0, };
    activeForceLength.resetPointsWith(afX, afY);
    return activeForceLength;
}

MtuCurve getPassiveForceLengthCurve() {

    MtuCurve passiveForceLength;
    std::vector<double> pfX{ 1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6};
    std::vector<double> pfY{ 0, 0.035, 0.12, 0.26, 0.55, 1.17, 2};
    passiveForceLength.resetPointsWith(pfX, pfY);
    return passiveForceLength;
}


MtuCurve getForceVelocityCurve() {
    
    MtuCurve forceVelocity;
    std::vector<double> fvX{ -1, -0.6, -0.3, -0.1, 0, 0.1, 0.3, 0.6, 0.8 };
    std::vector<double> fvY{ 0, 0.08, 0.2, 0.55, 1, 1.4, 1.6, 1.7, 1.75 };
    forceVelocity.resetPointsWith(fvX, fvY);
    return forceVelocity;
}

MtuCurve getTendonForceLengthCurve() {

    MtuCurve tendonForceLength;
    std::vector<double> tfX{ 0, 0.00131, 0.00281, 0.00431, 0.00581, 0.00731, 0.00881, 0.0103, 0.0118, 0.0123, 9.2};
    std::vector<double> tfY{0, 0.0108, 0.0257, 0.0435, 0.0652, 0.0915, 0.123, 0.161, 0.208, 0.227, 345};
    tendonForceLength.resetPointsWith(tfX, tfY);
    return tendonForceLength;
}

std::vector<double> getRange(double min, double max, unsigned nPoints) {

    std::vector<double> x;
    double dx = (max - min) / static_cast<double>(nPoints - 1);
    for (unsigned i(0); i < nPoints; ++i)
        x.push_back(min + static_cast<double>(i)*dx);
    return x;
}


void printVectors(const vector<vector<double>>& data, const vector<std::string>& headers, const std::string& filename) {

    std::ofstream oF(filename);
    for (auto h : headers)
        oF << h << ",";
    oF << std::endl;

    int nCols(data.size());
    int nRows(data.front().size());
    for (int iRow(0); iRow < nRows; ++iRow){
        for (int iCol(0); iCol < nCols; ++iCol) 
            oF << data.at(iCol).at(iRow) << ",";
        oF << std::endl;
    }
    oF.close();

}


template<typename T>
void testHeatRate(T mtu) {

	std::vector<double> timeColumn, fl, fv, mf, act, hr;
	double time(0.);
	for (auto mtl : getRange(0, 0, 2500)) {
		mtu.setMuscleTendonLength(mtl);
		mtu.setPercentFastTwitch(20);
		mtu.setMaxContractionVelocity(12);
		mtu.setOptimalFibreLength(8.9);
		mtu.setTime(time += 0.001);
		mtu.setActivation(0.9);
		mtu.setActivationScale(1.0);
		mtu.updateFibreLengthAndVelocity();
		mtu.updateMuscleForce();
		mtu.pushState();
		mtu.computeHeatRate();
		timeColumn.emplace_back(time);
		fl.emplace_back(mtu.getFiberLength());
		fv.emplace_back(mtu.getFiberVelocity());
		mf.emplace_back(mtu.getMuscleForce());
		act.emplace_back(mtu.getActivation());
		hr.emplace_back(mtu.computeHeatRate());

	}
	printVectors({ timeColumn, fl, fv, mf, act, hr},
	{ "time","fl", "fv", "mf", "act", "hr"}, "heatRate.csv");

}

int main() {


    MTU<ExponentialActivation, ElasticTendon_BiSec, CurveMode::Mode::Online> mtu("muscle");
    mtu.setCurves(
        getActiveForceLengthCurve(),
        getPassiveForceLengthCurve(),
        getForceVelocityCurve(),
        getTendonForceLengthCurve());

	// mtu.setParametersToComputeForces(0.17, 0.13, 0.83, 0.15, 0.2, 1, 1);
    // mtu.setParametersToComputeHeatRate(....);
	testHeatRate(mtu);
    return 0;

}

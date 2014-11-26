#ifndef OPENSIM_ROCINFORCECONTROLLER_H_
#define OPENSIM_ROCINFORCECONTROLLER_H_

/*********************************************************/
/*                    ROCIN Controller                   */
/*********************************************************/

#include <Simbody.h>
#include <OpenSim/Simulation/Control/Controller.h>
#include <OpenSim/Simulation/MarkersReference.h>
#include <OpenSim/Simulation/Model/ExternalLoads.h>
#include <OpenSim/Simulation/MomentArmSolver.h>
#include <OpenSim/Simulation/Manager/Manager.h>
#include "MPC.h"
#include <OpenSim/Common/LinearFunction.h>
#include <OpenSim/Simulation/Control/ControlSet.h>

class ConstantController;
class ROCINActuatorSystem;
class ROCINActuatorSubsystem;


namespace OpenSim {

class ROCINForceController : public Controller {
OpenSim_DECLARE_CONCRETE_OBJECT(ROCINForceController, Controller)

public:
	ROCINForceController() {}
	ROCINForceController(const Model& aModel, const std::string& dataFileName, bool track_coords = false, double lookahead_window = 0.01, int lookahead_number = 1);
	~ROCINForceController();

    // initialize the controller
	void initController(const Model& aModel, const std::string& dataFileName, bool track_coords, double lookahead_window, int lookahead_number);
    // set up the internal model
	void setUpInternalModel(const Model& aModel);
    // add additional actuators
	void addResidualAndReserveActuators(Model& aModel, const std::string& actuatorFiles);
    // create external load file
	void createExternalLoads(const std::string& aExternalLoadsFileName, Model& aModel);
    // override the computeControls method for the controller
	void computeControls(const SimTK::State& s, SimTK::Vector& controls) const;
    // core function that computes the control signals
	void doControlPrecomputation(const SimTK::State& s);
    // initialize the MPC solver
	void initMPCSolver(int numExtraObservations = 0);
    // set the initial state from a file
	SimTK::State& setInitStateFromFile(Model& aModel, const std::string& aStateFileName);
    // set the initial state without a file
	SimTK::State& setInitState(Model& aModel, double tiReal);
    //  compute the optimal initial state
    void computeOptimalInitialState(SimTK::State& initState, double tiReal);
    // initialize the system
	SimTK::State& initSystem(Model& aModel);
    // set the penalty
	void setCoordTrackingPenalty(double s) { _coord_tracking_penalty = s; }
	void setSpeedTrackingPenalty(double s) { _speed_tracking_penalty = s; }
	void setActivationPenalty(double s) { _activation_penalty = s; }
	void setPDPenalty(double s) { _PD_penalty = s; }
    // set the other parameters
	void setUseImplicitIntegForMPC(bool s) { _mpcSolver->setUseImplicit(s); }
	void setUseTaylorExpansion(bool s) { _use_taylor_expansion = s; }
	void setLowPassCutOffFrequency(double s) { _lowpass_freq = s; }
	void setPenalizeControlDerivative(bool s) { _mpcSolver->setPenalizeUdot(s); }
	void setPenalizeOutputDerivate(bool s) { _mpcSolver->setPenalizeYdot(s); }
	void setUsingVaryingDynamics(bool s) { _mpcSolver->setUsingVaryingDynamics(s); }
	void setModelNaturalFrequency(double s) { _mpcSolver->setNaturalFrequency(s); }

    //evaluate the coord tracking errors
    static void evalModelCoordTrackingErrs(const Model& aModel, const std::string& coordFileName1, const std::string& coordFileName2, const std::string& outFileName);
	void evalCoordTrackingErrs(const std::string& coordFileName1, const std::string& coordFileName2, const std::string& outFileName);

    // add event handler (calling doControlPrecomputation) to the system		
	void addEventHandlerToSystem();

    //set and get whether to check the target time in the event handler
	void setCheckTargetTime(bool aTrueFalse) { _checkTargetTime = aTrueFalse; }
	bool getCheckTargetTime() const { return _checkTargetTime; }
    //set and get target time that will trigger the event handler
	void setTargetTime(double aTargetTime) { _tf = aTargetTime; }
	double getTargetTime() const { return _tf; }
    //set and get the time step for the event handler
	void setTargetDT(double aDT) { _targetDT = aDT; }
	double getTargetDT() const { return _targetDT; }


	double getInitTime() { return _mpcSolver->getInitTime(); }
	double getFinalTime() { return _mpcSolver->getFinalTime(); }
	
    // set and get the state activations
	static void setStateActivation(const Model& aModel, const Vector& acts,SimTK::State& s);
    static void getStateActivation(const Model& aModel, const SimTK::State& s, int numMuscles, Vector& acts);
    // set and get the state fiber lengths
	static void setStateFiberLength(const Model& aModel, const Vector& fiberLens, SimTK::State& s);
	static void getStateFiberLength(const Model& aModel, const SimTK::State& s, int numMuscles,Vector& fiberLens);
    // get the muscle forces
	static void getMuscleForces(const Model& aModel, const SimTK::State& s, int numMuscles, Vector& muscleForces);
    // root solver
	static void rootSolve(void (*evaluate)(const Model&, const SimTK::State&,const double*, const Vector&,  Vector&), const Model& aModel, const SimTK::State& s, const Vector& y, const Vector& ax, const Vector& bx, const Vector& tol, Vector& x, const double* para = NULL);
    // evaluate the muscle forces based on the fiber lengths
	static void evaluateMuscleForceBasedOnFiberLength(const Model& aModel,const SimTK::State& s, const double* para, const Vector& fiberLength, Vector& muscleForces);
    // evaluate the muscle forces based on the activations
	static void evaluateMuscleForceBasedOnActivation(const Model& aModel, const SimTK::State& s, const double* para, const Vector& activations, Vector& muscleForces);
	
protected:
	void connectToModel(Model& model) OVERRIDE_11;
private:
	void setNull();
    
	void doPrecomputation_Varying_Dynamics(const SimTK::State& s);
    void addVirtualActuators(Model& aModel, const Array<std::string>& coords_virtual);
    void setInternalTrajectoryCorrections(const SimTK::State& s);
    void evaluateInternalModelForces(const SimTK::State& s, double* targetForces, double* controls, double* forces);
    void getInternalModelActivations(const SimTK::State& s, double* acts);
    void rootSolveActivations(const SimTK::State& s, const Vector& muscleForces, const Vector& lm, const Vector& lm_dot, const Vector& acts_min, const Vector& acts_max, const Vector& tol_acts, Vector& activations) const;
    void rootSolveMuscleFiberLength(const SimTK::State& s, const Vector& muscleForces, const Vector& lm_min, const Vector& lm_max, const Vector& tol, Vector& muscleFiberLengthVec);
    void setStateActivation(const Vector& acts, SimTK::State& s) const;
    void setStateFiberLength(const Vector& fiberLens, SimTK::State& s) const;
    void getStateActivation(const SimTK::State& s, Vector& acts) const;
    void getStateFiberLength(const SimTK::State& s, Vector& fiberLens) const;

    // combine the multibody dynamics and muscle dynamics
    void combineSystemDynamics(const Matrix& A_mbs, const Matrix& B_mbs, const Vector& C_mbs, const Vector& A_fm_u, const Vector& B_fm_u, Matrix& A, Matrix& B, Vector& C);
    // compute the approximated linear relationship for muscle forces and excitations
    void computeMuscleRelation(const SimTK::State& s, double dt, int n_dynamics, Array<Vector>& A_fm_u_array, Array<Vector>& B_fm_u_array, Array<Vector>& A_act_u_array, Array<Vector>& B_act_u_array);

    // compute N that map generaltized speeds to qdot
    void computeNqdot(const SimTK::State& s, Matrix& Nqdot);
    // compute multibody dynamics linear equation that use actuator forces as inputs
    void computeMBSDynamics(const SimTK::State& s, Matrix& B, Vector& A);
    void computeMBSDynamicsNumerically(const SimTK::State& s, Matrix& B, Vector& A);
    // compute the muscle Jacobian
    void computeMuscleJacobian(const SimTK::State& s, Matrix& C_muscles);
    // compute multibody dynamics linear equation that use generalized forces as inputs
    void computeMBSGeneralizedForceDynamics(const SimTK::State& s, Matrix& B, Vector& A);
    void computeMBSGeneralizedForceDynamicsNumerically(const SimTK::State& s, Matrix& B, Vector& A);

    // get the observation variables from the state
    void getCurObservation(const SimTK::State& s, Vector& y, bool includeMuscleForce = false, bool includeMuscleFiberLength = false, bool includeActivation = false);
    // compute multibody dynamics that use actuator forces as inputs
    void computeMBSStateDynamics(const SimTK::State& s, const Vector& actuatorForces, bool useTaylorExpansion, Matrix& A, Matrix& B, Vector& C);
    // compute multibody dynamics that use generalized forces as inputs
    void computeMBSStateDynamics(const SimTK::State& s, const Vector& genForces, const Matrix& muscleJacobian, bool useTaylorExpansion, Matrix& A, Matrix& B, Vector& C);
    // compute partial derivatives
    void computeDynamicsPartialsFwdDiff(const SimTK::State& s, const Vector& actuatorForces, Matrix& dUDotdQ, Matrix& dUDotdU);
    void computeDynamicsPartialsCentralDiff(const SimTK::State& s, const Vector& actuatorForces, Matrix& dUDotdQ, Matrix& dUDotdU);

    void computeGeneralizedForceDynamicsPartialsFwdDiff(const SimTK::State& s, const Vector& genForces, Matrix& dUDotdQ, Matrix& dUDotdU);
    void computeGeneralizedForceDynamicsPartialsCentralDiff(const SimTK::State& s, const Vector& genForces, Matrix& dUDotdQ, Matrix& dUDotdU);


	MarkerData* _markerData;
	Storage* _coordData;
	MomentArmSolver* _momentArmSolver;

	GCVSplineSet* _qSet;
	GCVSplineSet* _uSet;
	GCVSplineSet* _uDotSet;


	Model* _internalModel;
	SimTK::RungeKuttaMersonIntegrator* _internalIntegrator;
	Manager* _internalManager;
	ConstantController* _internalController;


	ROCINActuatorSystem* _internalSys;
	ROCINActuatorSubsystem* _internalSubsys;

	double _lowpass_freq;

	double _coord_tracking_penalty;
	double _speed_tracking_penalty;
	double _activation_penalty;
	double _PD_penalty;

	double _lookahead_window;
	int _lookahead_number;

	bool _track_coords;
	int _numCoords;

	int _n_q;
	int _n_u;
	int _numMuscles;
	int _numMarkers;

	int _numNonMuscleActuators;
	Matrix _mapping_nonMuscleActuators;
	Matrix _CoordsSelection;


	SimTK::Array_<SimTK::MobilizedBodyIndex> _mobilizedBodyIndices;
	SimTK::Array_<SimTK::Vec3> _stationPositionsInBodies;

	ExternalLoads* _externalLoads;

	Vector _control_lowerbounds;
	Vector _control_upperbounds;

	Array<double> _actuator_forces;
	Array<double> _actuator_controls;
	Array<double> _actuator_controls_new;

	bool _use_taylor_expansion;

	Vector _muscle_min_control_array;
	Vector _muscle_max_control_array;

	Array<std::string> _virtual_actuator_coord_array;
	Array<double> _virtual_actuator_weight_array;
	Vector _virtual_actuator_min_control_array;
	Vector _virtual_actuator_max_control_array;

	Array<Vector> _gen_forces_array;
	Array<Vector> _actuator_forces_array;

	ControlSet _controlSet;


	double _tf;
	double _targetDT;

	bool _checkTargetTime;

	MPC* _mpcSolver;

	
};
}	//namespace

#endif
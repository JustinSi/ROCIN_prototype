#ifndef OPENSIM_STATEINITIALIZER_H_
#define OPENSIM_STATEINITIALIZER_H_

#include <SimTKmath.h>
#include "ControllerHelper.h"
#include <OpenSim/Simulation/Model/Model.h>

namespace OpenSim {

class StateInitializer;

class StateInitializerQP : public SimTK::OptimizerSystem{
public:
	StateInitializerQP() { _stateInitializer = NULL;}
	StateInitializerQP(StateInitializer* stateInitializer) { _stateInitializer = stateInitializer; }
	void setStateInitializer(StateInitializer* stateInitializer) { _stateInitializer = stateInitializer; }

	int objectiveFunc(const Vector& coefficients, bool new_coefficients, SimTK::Real& f) const;
	int constraintFunc( const Vector& coefficients, bool new_coefficients, Vector &constraints) const;
	int gradientFunc( const Vector& coefficients, bool new_coefficients, Vector &gradient ) const;
	int constraintJacobian( const Vector& coefficients, bool new_coefficients, Matrix& jac)  const;

	void testQP();

private:
	StateInitializer* _stateInitializer;
};

class StateInitializer {
	friend class StateInitializerQP;
public:
	StateInitializer(Model* model, ConstantController* controller, const SimTK::State& initState);
	~StateInitializer() {}
	void setInitQ(const Vector& Q) { _initState.updQ() = Q; }
	void setInitU(const Vector& U) { _initState.updU() = U; }
	void setInitZ(const Vector& Z) { _initState.updZ() = Z; }
	void setUDotRef(const Vector& UDot) {_UDotRef = UDot; }

	void setControlWeights(const Vector& w_controls) { _w_controls = w_controls; }
	void setActuatorForceWeights(const Vector& w_actuator_forces) { _w_actuator_forces = w_actuator_forces; }
	void setActsWeights(double w_acts) { _w_acts = w_acts; }
	void setAccWeights(double w_acc) { _w_acc = w_acc; }
	void setControlBounds(const Vector& lower_bd, const Vector& upper_bd) { _lower_bd_controls = lower_bd; _upper_bd_controls = upper_bd; }
	void setZBounds(const Vector& lower_bd, const Vector& upper_bd) { _lower_bd_Z = lower_bd; _upper_bd_Z = upper_bd; }
	void setActuatorForceBounds(const Vector& lower_bd, const Vector& upper_bd) { _lower_bd_actuatorforces = lower_bd; _upper_bd_actuatorforces = upper_bd; }
	void setFiberLengthBounds(const Vector& lower_bd, const Vector& upper_bd) { _lm_min = lower_bd; _lm_max = upper_bd; }
	void setActivationBounds(const Vector& lower_bd, const Vector& upper_bd) { _acts_min = lower_bd; _acts_max = upper_bd; }
	void setMBSDynamicsAndOptForces(const Matrix& A_mbs, const Vector& B_mbs, const Vector& Opt_f);

	void rootSolveActivations(const SimTK::State& s, const Vector& muscleForces, const Vector& lm, const Vector& lm_dot, const Vector& acts_min, const Vector& acts_max, const Vector& tol_acts, Vector& activations) const;
	void rootSolveMuscleFiberLength(const SimTK::State& s, const Vector& muscleForces, const Vector& lm_min, const Vector& lm_max, const Vector& tol, Vector& muscleFiberLengthVec) const;


	SimTK::State getOptimizedInitState();

private:
	SimTK::State _initState;
	Vector _UDotRef;
	Model* _model;
	ConstantController* _controller;
	Vector _w_controls;
	double _w_acts;
	double _w_acc;

	Vector _w_actuator_forces;

	int _n_controls;
	int _n_Z;
	int _n_muscles;
	Vector _lower_bd_controls;
	Vector _upper_bd_controls;
	Vector _lower_bd_Z;
	Vector _upper_bd_Z;

	Vector _opt_actuator_forces;

	Vector _lower_bd_actuatorforces;
	Vector _upper_bd_actuatorforces;

	Vector _lm_min;
	Vector _lm_max;
	Vector _acts_min;
	Vector _acts_max;

	Matrix _A_f;
	Vector _B_f;

};

}
#endif
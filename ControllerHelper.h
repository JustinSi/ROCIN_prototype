#ifndef OPENSIM_CONTROLLERHELPER_H_
#define OPENSIM_CONTROLLERHELPER_H_

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Common/FunctionSet.h>
#include "SimTKmath.h"
#include <OpenSim/Common/Signal.h>
#include <OpenSim/Simulation/Control/Controller.h>
#include <OpenSim/Simulation/Model/ControllerSet.h>
#include <OpenSim/Simulation/Control/ControlConstant.h>
#include "SimTKsimbody.h"
#include "SimTKcommon/internal/SubsystemGuts.h"
#include "SimTKcommon/internal/SystemGuts.h"
#include <OpenSim/Simulation/Control/ControlLinear.h>
#include <OpenSim/Simulation/Control/ControlSet.h>

using SimTK::Matrix;
using SimTK::Vector;
using SimTK::Vec3;
using namespace OpenSim;

class ROCINActuatorSystemRep : public SimTK::System::Guts{
    public:
    ROCINActuatorSystemRep() : SimTK::System::Guts( "CMCActuatorSystem", "2.0") {}
    
    // Required by the System::Guts interface.
    ROCINActuatorSystemRep* cloneImpl() const 
    {   return new ROCINActuatorSystemRep(*this); }

    // This system doesn't have constraints, prescribed motion, or events so
    // we don't need to implement any of the related System::Guts methods.

    SimTK_DOWNCAST( ROCINActuatorSystemRep, SimTK::System::Guts );

};

class ROCINActuatorSystem : public SimTK::System {
   public:
   explicit ROCINActuatorSystem() {
       adoptSystemGuts( new ROCINActuatorSystemRep() );
       SimTK::DefaultSystemSubsystem defsub(*this);
   }

   SimTK_PIMPL_DOWNCAST( ROCINActuatorSystem, SimTK::System );
};

class ROCINActuatorSubsystemRep : public SimTK::Subsystem::Guts {

  public:
  ROCINActuatorSubsystemRep( Model* model);

  ROCINActuatorSubsystemRep* cloneImpl() const;
  ~ROCINActuatorSubsystemRep();

  void setCompleteState(const SimTK::State& s);
  const SimTK::State& getCompleteState() const;

  int realizeSubsystemTopologyImpl(SimTK::State& s) const;
  int realizeSubsystemDynamicsImpl(const SimTK::State& s) const;
  void setSpeedCorrections(const double corrections[] );
  void setCoordinateCorrections(const double corrections[] );
//  void setDesiredCoords(const double desiredCoords[]);
//  void setDesiredSpeeds(const double desiredCoords[]);

  void setTransitTargetTime(double t, double interval) { _targetTime_transit = t; _interval_transit = interval; }
  void setIntegTargetTime(double t) { _targetTime_integ = t; }

  double getIntegTargetTime() { return _targetTime_integ; }

  void assignQU(double t, SimTK::Vector& q, SimTK::Vector& u) const;
  void getEstimatedUDot(double t, SimTK::Vector& uDot) const;
  void assignTimeAndQU(double t, SimTK::State& s) const;


  void setSpeedTrajectories(FunctionSet *aSet);
  void setCoordinateTrajectories(FunctionSet *aSet);
  FunctionSet* getCoordinateTrajectories() const ;
  FunctionSet* getSpeedTrajectories() const ;
  Model* getModel() const ;
  
  void holdCoordinatesConstant( double t );
  void releaseCoordinates();

  SimTK::State  _completeState;
  Model*        _model;
  bool          _holdCoordinatesConstant;
  double        _holdTime;
  Array<double> _qCorrections;
  Array<double> _uCorrections;
//  Array<double> _qDesired;
//  Array<double> _uDesired;

  double _targetTime_integ;


  double _targetTime_transit;
  double _interval_transit;

  mutable Array<double> _qWork;
  mutable Array<double> _uWork;
  
  FunctionSet *_qSet;  
  FunctionSet *_uSet;

  SimTK_DOWNCAST( ROCINActuatorSubsystemRep, SimTK::Subsystem::Guts );
};


class ROCINActuatorSubsystem : public SimTK::Subsystem {
    public:
    explicit ROCINActuatorSubsystem( ROCINActuatorSystem& sys, Model* model);

    ROCINActuatorSubsystemRep *rep;

	void setCompleteState(const SimTK::State& s)  {
        rep->setCompleteState( s );
    }

    const SimTK::State& getCompleteState() const {
        return( rep->getCompleteState());
    }

	void setSpeedTrajectories(FunctionSet* aSet) {
		rep->setSpeedTrajectories(aSet);
	}

     void setCoordinateTrajectories(FunctionSet *aSet)  {
        rep->setCoordinateTrajectories( aSet );
    }

	 FunctionSet* getSpeedTrajectories() const {
		 return ( rep->getSpeedTrajectories());
	 }

     FunctionSet* getCoordinateTrajectories() const {
        return( rep->getCoordinateTrajectories() );
    }
    void setSpeedCorrections(const double *corrections ) {
            rep->setSpeedCorrections( corrections );
    }
    void setCoordinateCorrections(const double *corrections ) {
            rep->setCoordinateCorrections( corrections );
    }

//	void setDesiredCoords(const double* desiredCoords){
//		rep->setDesiredCoords(desiredCoords);
//	}

//	void setDesiredSpeeds(const double* desiredSpeeds){
//		rep->setDesiredSpeeds(desiredSpeeds);
//	}

	void setTransitTargetTime(double t, double interval){
		rep->setTransitTargetTime(t,interval);
	}

	void setIntegTargetTime(double t){
		rep->setIntegTargetTime(t);
	}

	double getIntegTargetTime(){
		return rep->getIntegTargetTime();
	}

    Model* getModel() const {
        return( rep->getModel() );
    }

	void getEstimatedUDot(double t, SimTK::Vector& uDot) const {
		rep->getEstimatedUDot(t,uDot);
	}

	void assignTimeAndQU(double t, SimTK::State& s) const {
		rep->assignTimeAndQU(t,s);
	}


	void holdCoordinatesConstant(double t) {
		rep->holdCoordinatesConstant(t);
	}

	void releaseCoordinates() {
		rep->releaseCoordinates();
	}

    SimTK_PIMPL_DOWNCAST(ROCINActuatorSubsystem, SimTK::Subsystem);
};

class ConstantController : public Controller{
OpenSim_DECLARE_CONCRETE_OBJECT(ConstantController, Controller)
public:
	ConstantController(int N);
	~ConstantController() {}
	//void setControls(const SimTK::Vector& controls, double t);
	void setControls(const double* controls, double t);
	void computeControls(const SimTK::State& s, SimTK::Vector& controls) const;
	ControlSet& getControlSet() { return _controlSet; }
private:
	Vector _controls;
	ControlSet _controlSet;
};


#endif
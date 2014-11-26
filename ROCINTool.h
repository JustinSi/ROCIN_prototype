#ifndef OPENSIM_ROCINTOOL_H_
#define OPENSIM_ROCINTOOL_H_

/***********************************************/
/*       A Wrapper for ROCIN Controller        */
/***********************************************/

#include <OpenSim/Tools/osimToolsDLL.h>
#include "ROCINForceController.h"


namespace OpenSim
{

class __declspec(dllexport) ROCINTool : public ModelComponent {
	OpenSim_DECLARE_CONCRETE_OBJECT(ROCINTool, ModelComponent);
public:
    //model file to apply ROCIN controler to 
	OpenSim_DECLARE_PROPERTY(model_file,std::string,"model file");
    //file that specifies the additional virtual (residual actuators and reserve actuators) actuators
	OpenSim_DECLARE_PROPERTY(additional_actuator_file,std::string,"additional actuators");
    //file directory to write ouput 
	OpenSim_DECLARE_PROPERTY(results_directory,std::string,"result directory");
    //whether use an existing state file to initialize the state
	OpenSim_DECLARE_PROPERTY(use_init_state_file, bool, "whether or not to use initial state file");
    //the state file for initializing the state (only effect when use_init_state_file is set to be true)
	OpenSim_DECLARE_PROPERTY(init_state_file,std::string, "initial state file");
    //initial time of simulation
	OpenSim_DECLARE_PROPERTY(initial_time,double,"initial time");
    //final time of simulation
	OpenSim_DECLARE_PROPERTY(final_time,double,"final time");
    //maximum integration time step
	OpenSim_DECLARE_PROPERTY(maximum_integrator_step_size,double,"maximum integrator step size");
    //error tolerance for the integrator
	OpenSim_DECLARE_PROPERTY(integrator_error_tolerance,double,"integrator error tolerance");
    //external loads file, e.g. Ground Reaction Force File
	OpenSim_DECLARE_PROPERTY(external_loads_file,std::string,"external loads file");
    //desired kinematics file
	OpenSim_DECLARE_PROPERTY(desired_kinematics_file,std::string,"desired kinematics file");
    //lowpass cutoff frequency for filtering the desired_kinematics_file
	OpenSim_DECLARE_PROPERTY(lowpass_cutoff_frequency,double,"lowpass cutoff frequency");
    //whether to use implicit formulation for MPC (using implicit formulation might slow down the controller, but is necessary for multiple prediction step)
	OpenSim_DECLARE_PROPERTY(use_implicit_mpc_integ, bool, "whether or not to use implicit integration formulation for MPC");
    //whether to add PD (acceleration) penalty, i.e. use PD rule to compute a desired acceleration and put the acceleration error in the objective function
	OpenSim_DECLARE_PROPERTY(add_PD_penalty, bool, "whether or not to add PD feedback penalty on output derivative");
    //whether to use Taylor expansion when linearizing the multibody dynamics
	OpenSim_DECLARE_PROPERTY(use_taylor_expansion, bool, "whether or not to use taylor expansion to approximate the multibody dynamics");
    //control update horizon for MPC
	OpenSim_DECLARE_PROPERTY(update_horizon,double,"update horizon");
    //number of steps to predict, prediction_horizon = prediction_step*update_horizon;
	OpenSim_DECLARE_PROPERTY(prediction_step,int,"prediction step");
    //joint coordinate tracking penaly
	OpenSim_DECLARE_PROPERTY(coord_tracking_penalty,double,"coordinate tracking penalty");
    //joint velocity tracking penalty
	OpenSim_DECLARE_PROPERTY(speed_tracking_penalty,double,"speed tracking penalty");
    //activation penalty
	OpenSim_DECLARE_PROPERTY(activation_penalty,double,"activation penalty");
    //joint acceleration tracking penalty
	OpenSim_DECLARE_PROPERTY(PD_penalty, double, "PD controlled acceleration penalty");
    //natural frequency of the model (used to specify the critically damped PD gains)
	OpenSim_DECLARE_PROPERTY(natural_frequency, double, "Natural Frequency of the model");

	ROCINTool(const std::string& setUpFile);
	~ROCINTool();

	bool run();
	static void evalCoordTrackingErrs(const std::string& modelFile, const std::string& desired_kinematics_file, const std::string& state_file, const std::string& err_file);	

private:
	void constructProperties();
	ROCINForceController* _controller;
};
};
#endif

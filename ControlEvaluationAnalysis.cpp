#include "ControlEvaluationAnalysis.h"
#include <OpenSim/Simulation/Model/PathActuator.h>
#include <OpenSim/Actuators/CoordinateActuator.h>

namespace OpenSim {
ControlEvaluationAnalysis ::ControlEvaluationAnalysis(Model *aModel):
Analysis(aModel), _analysisStorage(1000)
{
	setName("ControlAnalysis");
	_momentArmSolver = NULL;

	setupStorage();
}

ControlEvaluationAnalysis::~ControlEvaluationAnalysis()
{
	if(_momentArmSolver!=NULL)
		delete _momentArmSolver;
}

void ControlEvaluationAnalysis::initMomentArmSolver()
{
	if(_momentArmSolver!=NULL)
		delete _momentArmSolver;

	_momentArmSolver = new MomentArmSolver(*_model);

	const SimTK::State& s = _model->getWorkingState();
	const SimTK::SimbodyMatterSubsystem& smss = _model->getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = _model->getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Velocity);

	const Set<Actuator>& actuatorSet = _model->getActuators();

	int numMuscles = 0;
	for(int i=0;i<actuatorSet.getSize();i++)
	{
		Actuator& act = actuatorSet.get(i);
		PathActuator* m = dynamic_cast<PathActuator*>(&act);
		
		if(m!=NULL)
		{
			numMuscles++;			
		}
	}

	_momentArmMask.resize(numMuscles);

	int n_u = _model->getNumCoordinates();

	
	int idx_musc = 0;
	for(int i=0;i<actuatorSet.getSize();i++)
	{
		Actuator& act = actuatorSet.get(i);
		PathActuator* m = dynamic_cast<PathActuator*>(&act);
		
		if(m!=NULL)
		{
			_momentArmMask[idx_musc].resize(n_u);
			for(int j=0;j<n_u;j++)
			{
				Coordinate& coord = _model->getCoordinateSet().get(j);
				if(!coord.isDependent(s))
				{
					if(_momentArmSolver->solve(s,coord,m->getGeometryPath()) != 0)
						_momentArmMask[idx_musc][j] = true;
					else
						_momentArmMask[idx_musc][j] = false;
				}
				else
					_momentArmMask[idx_musc][j] = false;
			}

			idx_musc++;
		}
	}
}

void ControlEvaluationAnalysis::setupStorage()
{
	_storageList.append(&_analysisStorage);
	_storageList.setMemoryOwner(false);
}

void ControlEvaluationAnalysis::constructColumnLabels()
{
	

	if(_model)
	{
		Array<std::string> labels;
		labels.append("time");

		int num_coords = _model->getNumCoordinates();

		for(int i=0;i<num_coords;i++)
		{
			Coordinate& coord = _model->getCoordinateSet().get(i);
			labels.append(coord.getName()+"_muscleContribution");
		}

		for(int i=0;i<num_coords;i++)
		{
			Coordinate& coord = _model->getCoordinateSet().get(i);
			labels.append(coord.getName()+"_resContribution");
		}

		_analysisStorage.setColumnLabels(labels);
	}
}

int ControlEvaluationAnalysis::begin(SimTK::State& s)
{
	if(!proceed()) return(0);

	constructColumnLabels();

	_analysisStorage.reset(s.getTime());

	int status = 0;
	if(_analysisStorage.getSize()<=0)
		status = record(s);

	return (status);
}

int ControlEvaluationAnalysis::step(const SimTK::State& s, int stepNumber)
{
	if(!proceed(stepNumber)) return(0);
	record(s);
	return(0);
}

int ControlEvaluationAnalysis::end(SimTK::State& s)
{
	if(!proceed()) return 0;
	record(s);
	return(0);
}

int ControlEvaluationAnalysis::
record(const SimTK::State& s)
{
	if(_model == NULL) return (-1);
	
	int n_u = _model->getNumCoordinates();
	SimTK::Vector result(n_u*2);
	result.setToZero();

	const SimTK::SimbodyMatterSubsystem& smss = _model->getMatterSubsystem();		
	const SimTK::MultibodySystem& mbs = _model->getMultibodySystem();
	mbs.realize(s,SimTK::Stage::Dynamics);

	const Set<Actuator>& actuatorSet = _model->getActuators();

	int idx_musc = 0;

	for(int i=0;i<actuatorSet.getSize();i++)
	{
		//Actuator& act = actuatorSet.get(i);
		//double actForce = act.getForce(s);

        ScalarActuator& act = (ScalarActuator&)actuatorSet.get(i);
        double actForce = act.getActuation(s);


		PathActuator* m = dynamic_cast<PathActuator*>(&act);
		
		if(m!=NULL)
		{
			
			for(int j=0;j<n_u;j++)
			{
				Coordinate& coord = _model->getCoordinateSet().get(j);
				if(!coord.isDependent(s) && _momentArmMask[idx_musc][j])
				{
					result[j] += _momentArmSolver->solve(s,coord,m->getGeometryPath())*actForce;
				}
			}

			idx_musc++;
		}
		else
		{
			CoordinateActuator* coordAct = dynamic_cast<CoordinateActuator*>(&act);			
			if(coordAct != NULL)
			{
				int coord_idx = _model->getCoordinateSet().getIndex(coordAct->getCoordinate(),0);
				result[n_u+coord_idx] = actForce;
			}
		}
	}

	_analysisStorage.append(s.getTime(),result);

	return(0);
}

int ControlEvaluationAnalysis::
printResults(const std::string &aBaseName,const std::string &aDir,
		double aDT,const std::string &aExtension)
{
	if(!getOn())
	{
		printf("ControlEvaluationAnalysis.printResults: Off - not printing.\n");
		return(0);
	}

	std::string prefix = aBaseName + "_"+getName()+"_";

	Storage::printResult(&_analysisStorage,prefix+"Result",aDir,aDT,aExtension);

	return(0);
}


}

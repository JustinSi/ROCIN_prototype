#ifndef OPENSIM_CONTROLEVALUATIONANALYSIS_H_
#define OPENSIM_CONTROLEVALUATIONANALYSIS_H_

#include <OpenSim/Simulation/Model/Model.h>
#include <OpenSim/Simulation/Model/Analysis.h>
//#include <OpenSim/Simulation/Model/Probe.h>
#include <OpenSim/Simulation/MomentArmSolver.h>

namespace OpenSim {

class __declspec(dllexport) ControlEvaluationAnalysis : public Analysis {
	OpenSim_DECLARE_CONCRETE_OBJECT(ControlEvaluationAnalysis, Analysis);

public:
	ControlEvaluationAnalysis(Model *aModel=0);
	~ControlEvaluationAnalysis();

	void initMomentArmSolver();

	virtual int printResults(const std::string &aBaseName,const std::string &aDir="",
		double aDT=-1.0,const std::string &aExtension=".sto");

	virtual int begin(SimTK::State& s);
	virtual int step(const SimTK::State& s, int stepNumber);
	virtual int end(SimTK::State& s);

	const Storage& getAnalysisStorage() const
	{
		return _analysisStorage;
	};
	Storage& updAnalysisStorage()
	{
		return _analysisStorage;
	}

protected:
    virtual int record(const SimTK::State& s );
private:
	void constructColumnLabels();
	void setupStorage();

private:

	Storage _analysisStorage;
	MomentArmSolver* _momentArmSolver;
	std::vector<std::vector<bool>> _momentArmMask;
};

}

#endif
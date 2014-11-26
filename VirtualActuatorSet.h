#ifndef OPENSIM_VIRTUALACTUATORSET_H_
#define OPENSIM_VIRTUALACTUATORSET_H_

#include<OpenSim/Simulation/Model/ComponentSet.h>

namespace OpenSim {

class VirtualActuator : public Object {
OpenSim_DECLARE_CONCRETE_OBJECT(VirtualActuator, Object);
public:
	OpenSim_DECLARE_PROPERTY(coordinate, std::string, "coordinate that the virtual actuator applies to");
	OpenSim_DECLARE_PROPERTY(weight, double, "penalty weight of the virtual actuator");
	OpenSim_DECLARE_PROPERTY(min_control, double, "minimum control value of the virtual actuator");
	OpenSim_DECLARE_PROPERTY(max_control, double, "maximum control value of the virtual actuator");

	VirtualActuator();
	VirtualActuator(const std::string& coord, double weight=1.0, double minVal = -200.0, double maxVal = 200.0);

private:
	void constructProperties();
};

class VirtualActuatorSet : public Set<VirtualActuator> {
OpenSim_DECLARE_CONCRETE_OBJECT(VirtualActuatorSet, Set<VirtualActuator>);
public:
	VirtualActuatorSet();
	VirtualActuatorSet(const std::string& aFileName);
	virtual ~VirtualActuatorSet();

	void getCoordinateArray(Array<std::string>& coord_array) const;
	void getWeightArray(Array<double>& weight_array) const;
	void getCoordAndWeightArray(Array<std::string>& coord_array, Array<double>& weight_array) const;
	void getMinAndMaxControlArray(SimTK::Vector& min_control_array, SimTK::Vector& max_control_array) const;
};

};

#endif
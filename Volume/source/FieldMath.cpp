#include "FieldMath.hpp"
using namespace lux;

VField lux::operator * (VField vf, SField sf) { return std::make_shared<VFMultiply>(vf, sf); }
VField lux::operator + (VField vf1, VField vf2) { return std::make_shared<VFAdd>(vf1, vf2); }
VField lux::operator - (VField vf1, VField vf2) { return std::make_shared<VFSubtract>(vf1, vf2); }

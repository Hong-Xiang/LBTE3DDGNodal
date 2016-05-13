#include "utilities.h"

dgn::quadrant dgn::quadrant_of_angle(num_t mu, num_t xi, num_t eta)
{
	if (mu >= 0.0 && xi >= 0.0 && eta >= 0.0)
		return quadrant::quadrant1;
	if (mu < 0.0 && xi >= 0.0 && eta >= 0.0)
		return quadrant::quadrant2;
	if (mu < 0.0 && xi < 0.0 && eta >= 0.0)
		return quadrant::quadrant3;
	if (mu >= 0.0 && xi < 0.0 && eta >= 0.0)
		return quadrant::quadrant4;
	if (mu >= 0.0 && xi >= 0.0 && eta < 0.0)
		return quadrant::quadrant5;
	if (mu < 0.0 && xi >= 0.0 && eta < 0.0)
		return quadrant::quadrant6;
	if (mu < 0.0 && xi < 0.0 && eta < 0.0)
		return quadrant::quadrant7;
	if (mu >= 0.0 && xi < 0.0 && eta < 0.0)
		return quadrant::quadrant8;
}

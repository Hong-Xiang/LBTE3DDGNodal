#include "utilities.h"

dgn::quadrant dgn::utilities::quadrant_of_angle(num_t mu, num_t xi, num_t eta)
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
	//if (mu >= 0.0 && xi < 0.0 && eta < 0.0)
	//	return quadrant::quadrant8;
	return quadrant::quadrant8;
}

bool dgn::utilities::quadrant_x_positive_flag(quadrant quad)
{
	return (quad == quadrant::quadrant1 || 
			quad == quadrant::quadrant4 || 
			quad == quadrant::quadrant5 || 
			quad == quadrant::quadrant8);
}
bool dgn::utilities::quadrant_y_positive_flag(quadrant quad)
{
	return (quad == quadrant::quadrant1 ||
			quad == quadrant::quadrant2 ||
			quad == quadrant::quadrant5 ||
			quad == quadrant::quadrant6);
}
bool dgn::utilities::quadrant_z_positive_flag(quadrant quad)
{
	return (quad == quadrant::quadrant1 ||
			quad == quadrant::quadrant2 ||
			quad == quadrant::quadrant3 ||
			quad == quadrant::quadrant4);
}

bool dgn::utilities::is_boundary(interface_direction idir, size_t np, size_t ix, size_t iy, size_t iz)
{
	bool ans = false;
	if (idir == interface_direction::B && ix == 0)
		ans = true;
	if (idir == interface_direction::F && ix == np-1)
		ans = true;
	if (idir == interface_direction::L && iy == 0)
		ans = true;
	if (idir == interface_direction::R && iy == np-1)
		ans = true;
	if (idir == interface_direction::D && iz == 0)
		ans = true;
	if (idir == interface_direction::U && iz == np-1)
		ans = true;
	return ans;
}

bool dgn::utilities::is_boundary(interface_direction idir, size_t np, size_t ix, size_t iy)
{
	bool ans = false;
	if (idir == interface_direction::B && ix == 0)
		ans = true;
	if (idir == interface_direction::F && ix == np - 1)
		ans = true;
	if (idir == interface_direction::L && iy == 0)
		ans = true;
	if (idir == interface_direction::R && iy == np - 1)
		ans = true;		
	return ans;
}

bool dgn::utilities::is_boundary(interface_direction idir, size_t np, size_t ix)
{
	bool ans = false;
	if (idir == interface_direction::B && ix == 0)
		ans = true;
	if (idir == interface_direction::F && ix == np - 1)
		ans = true;
	return ans;
}
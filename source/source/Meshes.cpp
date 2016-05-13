#include "Meshes.h"
#include "Elements.h"
#include "Interfaces.h"

DG3DNodal::Mesh::Mesh(size_t nx, size_t ny, size_t nz, num_t h, num_t xc, num_t yc, num_t zc)
{
	nx_ = nx;
	ny_ = ny;
	nz_ = nz;
	h_ = h;
	xc_ = xc; yc_ = yc; zc_ = zc;
}

void DG3DNodal::Mesh::generate_mesh()
{
	size_t total_elements = nx_*ny_*nz_;
	size_t total_interfaces = total_elements * 3 + 3 * nx_*ny_;
	il.reserve(total_interfaces);
	index3.resize(boost::extents[nx_][ny_][nz_]);
	m_ix.resize(nx_);
	m_iy.resize(ny_);
	m_iz.resize(nz_);
	//Generate Elements
	for(size_t i = 0; i < total_elements; ++i)
	{
		size_t ix, iy, iz;
		DG3DNodal::ind2sub3(nx_, i, ix, iy, iz);
		num_t xct, yct, zct;
		xct = ix*h_ - nx_*h_ / 2.0 + xc_;
		yct = iy*h_ - ny_*h_ / 2.0 + yc_;
		zct = iz*h_ - nz_*h_ / 2.0 + zc_;
		el.push_back(Element3D(i, 0, ix, iy, iz, xct, yct, zct));
		index3[ix][iy][iz] = i;
	}	
	size_t idi;
	for (size_t i = 0; i < el.size(); ++i)
	{
		size_t ix, iy, iz;		
		DG3DNodal::ind2sub3(nx_, i, ix, iy, iz);				
		idi = add_interface(0);		
		el.at(i).addInterface(InterfaceDirectionBack, il.at(idi).getHandle());		
		idi = add_interface(0);		
		el.at(i).addInterface(InterfaceDirectionLeft, il.at(idi).getHandle());		
		idi = add_interface(0);
		el.at(i).addInterface(InterfaceDirectionUnder, il.at(idi).getHandle());		
	}
	for (size_t ix = 0; ix < nx_; ix++)
	{
		for (size_t iy = 0; iy < ny_; iy++)
		{
			for (size_t iz = 0; iz < nz_; iz++)
			{
				size_t idx = index3[ix][iy][iz];
				if (ix == nx_-1)
				{
					size_t tid = add_interface(0);
					el.at(idx).addInterface(InterfaceDirectionFront, il.at(tid).getHandle());
				}
				else
				{
					size_t idf = index3[ix + 1][iy][iz];
					el.at(idx).addInterface(InterfaceDirectionFront, el.at(idf).get_interface(InterfaceDirectionBack));
				}
				if (iy == ny_-1)
				{
					size_t tid = add_interface(0);
					el.at(idx).addInterface(InterfaceDirectionRight, il.at(tid).getHandle());
				}
				else
				{
					size_t idf = index3[ix][iy+1][iz];
					el.at(idx).addInterface(InterfaceDirectionRight, el.at(idf).get_interface(InterfaceDirectionLeft));
				}
				if (iz == nz_-1)
				{
					size_t tid = add_interface(0);
					el.at(idx).addInterface(InterfaceDirectionTop, il.at(tid).getHandle());
				}
				else
				{
					size_t idf = index3[ix][iy][iz+1];
					el.at(idx).addInterface(InterfaceDirectionTop, el.at(idf).get_interface(InterfaceDirectionUnder));
				}
			}
		}
	}
	for (size_t ix = 0; ix < nx_; ix++)
	{
		for (size_t iy = 0; iy < ny_; iy++)
		{
			el.at(index3[ix][iy][0]).get_interface(InterfaceDirectionUnder)->markAsBoundary();
			el.at(index3[ix][iy][nz_-1]).get_interface(InterfaceDirectionTop)->markAsBoundary();
		}
	}
	for (size_t ix = 0; ix < nx_; ix++)
	{
		for (size_t iz = 0; iz < nz_; iz++)
		{
			el.at(index3[ix][0][iz]).get_interface(InterfaceDirectionLeft)->markAsBoundary();
			el.at(index3[ix][ny_-1][iz]).get_interface(InterfaceDirectionRight)->markAsBoundary();
		}
	}
	for (size_t iz = 0; iz < nz_; iz++)
	{
		for (size_t iy = 0; iy < ny_; iy++)
		{
			el.at(index3[0][iy][iz]).get_interface(InterfaceDirectionBack)->markAsBoundary();
			el.at(index3[nx_-1][iy][iz]).get_interface(InterfaceDirectionFront)->markAsBoundary();
		}
	}
#ifdef DEBUG
	for (size_t i = 0; i < il.size(); i++)
	{
		std::cout << il.at(i).id() << std::endl;
	}
#endif 
}

DG3DNodal::size_matrix_t DG3DNodal::Mesh::get_element_interface_matrix() const
{
	size_matrix_t ans;
	size_t totalElement = nx_ * ny_ * nz_;
	ans.resize(totalElement);
	size_t iid;
	for (size_t ie = 0; ie < totalElement; ++ie)
	{
		ans.at(ie).resize(6);
		iid = el.at(ie).get_interface(InterfaceDirectionBack)->id();
		ans.at(ie).at(InterfaceDirectionBack) = iid;
		iid = el.at(ie).get_interface(InterfaceDirectionFront)->id();
		ans.at(ie).at(InterfaceDirectionFront) = iid;
		iid = el.at(ie).get_interface(InterfaceDirectionLeft)->id();
		ans.at(ie).at(InterfaceDirectionLeft) = iid;
		ans.at(ie).at(InterfaceDirectionRight) = iid = el.at(ie).get_interface(InterfaceDirectionRight)->id();
		ans.at(ie).at(InterfaceDirectionUnder) = iid = el.at(ie).get_interface(InterfaceDirectionUnder)->id();
		ans.at(ie).at(InterfaceDirectionTop) = iid = el.at(ie).get_interface(InterfaceDirectionTop)->id();
	}
	return ans;
}

DG3DNodal::size_array_t DG3DNodal::Mesh::get_sweep_order(Quadrant quad) const
{
	size_array_t sweep_order;
	sweep_order.resize(total_element());
	int x0, x1, y0, y1, z0, z1;
	int step_x, step_y, step_z;
	if (quad == Quadrant1 || quad == Quadrant2 || quad == Quadrant5 || quad == Quadrant6)
	{
		x0 = 0; x1 = nx_ - 1; step_x = 1;
	}
	else
	{
		x0 = nx_ - 1; x1 = 0; step_x = -1;
	}
	if (quad == Quadrant1 || quad == Quadrant4 || quad == Quadrant5 || quad == Quadrant8)
	{
		y0 = 0; y1 = ny_ - 1; step_y = 1;
	}
	else
	{
		y0 = ny_ - 1; y1 = 0; step_y = -1;
	}
	if (quad == Quadrant1 || quad == Quadrant2 || quad == Quadrant3 || quad == Quadrant4)
	{
		z0 = 0; z1 = nz_ - 1; step_z = 1;
	}
	else
	{
		z0 = nz_ - 1; z1 = 0; step_z = -1;
	}
	size_t current_element = 0;
	for (int ix = x0; ix != (x1 + step_x); ix += step_x)
		for (int iy = y0; iy != (y1 + step_y); iy += step_y)
			for (int iz = z0; iz != (z1 + step_z); iz += step_z)
			{
				sweep_order.at(current_element) = index3[ix][iy][iz];
				current_element++;
			}
	return sweep_order;
}

DG3DNodal::size_array_t DG3DNodal::Mesh::get_refine_level() const
{
	size_array_t ans;
	size_t totalElement = nx_ * ny_ * nz_;
	ans.resize(totalElement);
	for (size_t ie = 0; ie < totalElement; ++ie)
	{
		ans.at(ie) = el.at(ie).get_refine_level();
	}
	return ans;
}

std::vector<DG3DNodal::vector3> DG3DNodal::Mesh::elements_center() const
{
	std::vector<vector3> ans;
	ans.resize(total_element());
	for (size_t i = 0; i < total_element(); i++)
	{
		ans.at(i).x = el.at(i).getXc();
		ans.at(i).y = el.at(i).getYc();
		ans.at(i).z = el.at(i).getZc();
	}
	return ans;
}

size_t DG3DNodal::Mesh::add_interface(mesh_refine_lvl_t refine_level)
{
	il.push_back(Interface(il.size(), refine_level));
	return il.size()-1;
}

void DG3DNodal::Mesh::print_index_of_elements(std::ostream& os)
{
	os << "id:\tix:\tiy:\tiz:" << std::endl;
	for (size_t i = 0; i < el.size(); i++)
	{
		os << i << "\t" << el.at(i).getIX() << "\t" << el.at(i).getIY() << "\t" << el.at(i).getIZ() << std::endl;
	}
}

std::vector<bool> DG3DNodal::Mesh::get_boundary_mark() const
{
	std::vector<bool> ans(il.size());
	for (size_t i = 0; i < il.size(); i++)
	{
		ans.at(i) = il.at(i).isBoundary();
	}
	return ans;
}


#include "Mesh.h"

Mesh::Mesh(size_t nx, size_t ny, size_t nz, num_t h, num_t xc, num_t yc, num_t zc)
	: nx_(nx), ny_(ny), nz_(nz), h_(h), xc_(xc), yc_(yc), zc_(zc)
{
	element_list_.clear();
	surface_list_.clear();
	sweep_order_cache_.clear();
	sweep_order_cache_.resize(kQuadrantTotal.at(3));
}

void Mesh::generate()
{
	generate_coarse();
}

#include <stack>
void Mesh::cache_sweep_order()
{

	for (size_t i = 0; i < kQuadrantTotal.at(3); i++)
	{
		quadrant quad = quadrant_list.at(i);

		size_t iquad = static_cast<size_t>(quad);

		std::vector<InterfaceDirection> search_directions;
		search_directions.clear();
		std::vector<bool> revert_sweep;
		revert_sweep.clear();

		if (quad == quadrant::quadrant1 || quad == quadrant::quadrant4 || quad == quadrant::quadrant5 || quad == quadrant::quadrant8)
		{
			search_directions.push_back(InterfaceDirection::B);
			revert_sweep.push_back(false);
		}
		else {
			search_directions.push_back(InterfaceDirection::F);
			revert_sweep.push_back(true);
		}
		if (quad == quadrant::quadrant1 || quad == quadrant::quadrant2 || quad == quadrant::quadrant5 || quad == quadrant::quadrant6)
		{
			search_directions.push_back(InterfaceDirection::L);
			revert_sweep.push_back(false);
		}
		else {
			search_directions.push_back(InterfaceDirection::R);
			revert_sweep.push_back(true);
		}
		if (quad == quadrant::quadrant1 || quad == quadrant::quadrant2 || quad == quadrant::quadrant3 || quad == quadrant::quadrant4)
		{
			search_directions.push_back(InterfaceDirection::D);
			revert_sweep.push_back(false);
		}
		else {
			search_directions.push_back(InterfaceDirection::U);
			revert_sweep.push_back(true);
		}



		std::vector<size_t> ans;
		ans.clear();

		std::vector<bool> added;
		added.resize(elements_total());
		for (size_t i = 0; i < elements_total(); i++)
		{
			added.at(i) = false;
		}
		std::stack<size_t> stk;
		while (!stk.empty())
			stk.pop();

		for (size_t i = 0; i < elements_total(); i++)
		{
			if (!added.at(i))
				stk.push(i);

			while (!stk.empty())
			{
				size_t e_now = stk.top();
				bool is_bd_el = true; //is boundary Element
				for (size_t j = 0; j < 3; j++)
				{
					size_t surface_search = element_list_.at(e_now).surface_adjacent(search_directions.at(j));
					const Surface& s_ref = surface_list_.at(surface_search);
					if (s_ref.is_boundary())
						continue;
					if (!revert_sweep.at(j)) {
						for (size_t k = 0; k < s_ref.pre_element_total(); k++)
						{
							size_t element_search = s_ref.element_pre(k);
							if (added.at(element_search))
								continue;
							is_bd_el = false;
							stk.push(element_search);
						}
					}
					else {
						for (size_t k = 0; k < s_ref.inc_element_total(); k++)
						{
							size_t element_search = s_ref.element_inc(k);
							if (added.at(element_search))
								continue;
							is_bd_el = false;
							stk.push(element_search);
						}
					}
				}

				if (is_bd_el)
				{
					ans.push_back(e_now);
					added.at(e_now) = true;
					stk.pop();
				}

			}
		}

		sweep_order_cache_.at(iquad) = ans;
	}
}

std::vector<size_t> Mesh::sweep_order(quadrant quad) const
{
	return sweep_order_cache_.at(static_cast<size_t>(quad));
}

boost::multi_array<size_t, 2> Mesh::element_interface_matrix() const
{
	boost::multi_array<size_t, 2> ans;
	ans.resize(boost::extents[element_list_.size()][kInterfaceTotal.at(3)]);
	for each (const Element& el in element_list_)
	{
		for (size_t idr = 0; idr < kInterfaceTotal.at(3); idr++)
		{
			InterfaceDirection dire = InterfaceDirection_list.at(idr);
			ans[el.id()][idr] = el.surface_adjacent(dire);
		}
	}
	return ans;
}

std::vector<size_t> Mesh::refine_level() const
{
	std::vector<size_t> ans;
	ans.clear();
	ans.resize(element_list_.size());
	for (size_t i = 0; i < ans.size(); i++)
	{
		ans.at(i) = element_list_.at(i).refine_level();
	}
	return ans;
}

num_t Mesh::h_element(size_t id) const
{
	return element_list_.at(id).h();
}

num_t Mesh::h_surface(size_t id) const
{
	return surface_list_.at(id).h();
}

std::vector<vector3> Mesh::center() const
{
	std::vector<vector3> ans;
	ans.clear();
	ans.resize(element_list_.size());
	for (size_t i = 0; i < ans.size(); i++)
	{
		ans.at(i) = vector3(element_list_.at(i).xc(), element_list_.at(i).yc(), element_list_.at(i).zc());
	}
	return ans;
}

std::vector<vector3> Mesh::center_surface() const
{
	std::vector<vector3> ans;
	ans.clear();
	ans.resize(surface_list_.size());
	for (size_t i = 0; i < ans.size(); i++)
	{
		ans.at(i) = vector3(surface_list_.at(i).xc(), surface_list_.at(i).yc(), surface_list_.at(i).zc());
	}
	return ans;
}

std::vector<size_t> Mesh::memory_element(size_t offset) const
{
	std::vector<size_t> ans;
	ans.clear();
	ans.resize(elements_total());
	size_t c_pos = 0;
	for (size_t i = 0; i < elements_total(); i++)
	{
		ans.at(i) = c_pos;
		if (element_list_.at(i).is_active())
		{
			c_pos += offset;
		}
	}
	return ans;
}

std::vector<size_t> Mesh::memory_surface(size_t offset) const
{
	std::vector<size_t> ans;
	ans.clear();
	ans.resize(surfaces_total());
	size_t c_pos = 0;
	for (size_t i = 0; i < surfaces_total(); i++)
	{
		ans.at(i) = c_pos;				
		c_pos += offset;	
	}
	return ans;
}

size_t Mesh::memory_element_total(size_t offset) const
{
	size_t ans = 0;
	for (size_t i = 0; i < element_list_.size(); i++)
	{
		if (element_list_.at(i).is_active())
		{
			ans += offset;
		}
	}
	return ans;
}

size_t Mesh::memory_surface_total(size_t offset) const
{
	size_t ans = 0;
	for (size_t i = 0; i < surface_list_.size(); i++)
	{		
			ans += offset;	
	}
	return ans;
}

size_t Mesh::add_surface(size_t refine_level, SurfaceDirection dir, num_t xc, num_t yc, num_t zc, num_t h)
{
	size_t id = surface_list_.size();
	surface_list_.push_back(Surface(id, refine_level, dir, xc, yc, zc, h));
//#ifdef _DEBUG
//	std::cerr << "add Surface called" << std::endl;
//#endif
	return id;
}

size_t Mesh::add_element(size_t refine_level, num_t xc, num_t yc, num_t zc, num_t h)
{
	size_t id = element_list_.size();
	element_list_.push_back(Element(id, refine_level, xc, yc, zc, h));
//#ifdef _DEBUG
//	std::cerr << "add Element called" << std::endl;
//#endif
	return id;
}

std::vector<bool> Mesh::boundary_mark(size_t offset) const
{
	std::vector<bool> ans;
	ans.clear();
	ans.resize(memory_surface_total(offset));
	for (size_t i = 0; i < surface_list_.size(); i++)
	{
		for (size_t j = 0; j < offset; j++)
		{
			ans.at(i*offset+j) = surface_list_.at(i).is_boundary();
		}
	}
	return ans;
}

SurfaceDirection Mesh::get_SurfaceDirection(size_t id) const
{
	return surface_list_.at(id).direction();
}

void Mesh::generate_coarse()
{
	size_t total_elements = nx_*ny_*nz_;
	size_t total_interfaces = total_elements * 3 + 3 * nx_*ny_;
	surface_list_.reserve(total_interfaces);
	boost::multi_array<size_t, 3> index3;
	index3.resize(boost::extents[nx_][ny_][nz_]);

	//Generate Elements	

	for (size_t i = 0; i < total_elements; i++)
	{
		size_t ix, iy, iz;
		xlib::ind2sub(nx_, ny_, nz_, i, ix, iy, iz);
		num_t xct, yct, zct;
		xct = ix*h_ - (nx_ - 1) * h_ / 2.0 + xc_;
		yct = iy*h_ - (ny_ - 1) * h_ / 2.0 + yc_;
		zct = iz*h_ - (nz_ - 1) * h_ / 2.0 + zc_;
		index3[ix][iy][iz] = add_element(0, xct, yct, zct, h_);
//#ifdef _DEBUG
//		std::cout << "ix= " << ix << "\tiy= " << iy << "\tiz= " << iz << "\tid3= " << index3[ix][iy][iz] << std::endl;
//#endif
	}



	size_t ids;
	for (size_t iz = 0; iz < nz_; iz++)
		for (size_t iy = 0; iy < ny_; iy++)
			for (size_t ix = 0; ix < nx_; ix++)
			{
				size_t ide = index3[ix][iy][iz];
				num_t xce = element_list_.at(ide).xc();
				num_t yce = element_list_.at(ide).yc();
				num_t zce = element_list_.at(ide).zc();
				num_t he = element_list_.at(ide).h();
				ids = add_surface(0, SurfaceDirection::X, xce-he/2.0, yce, zce, he);
				element_list_.at(ide).link_surface(InterfaceDirection::B, surface_list_.at(ids).id());
				surface_list_.at(ids).link_element(ElementDirection::inc, element_list_.at(ide).id());
				ids = add_surface(0, SurfaceDirection::Y, xce, yce-he/2.0, zce, he);
				element_list_.at(ide).link_surface(InterfaceDirection::L, surface_list_.at(ids).id());
				surface_list_.at(ids).link_element(ElementDirection::inc, element_list_.at(ide).id());
				ids = add_surface(0, SurfaceDirection::Z, xce, yce, zce-he/2.0, he);
				element_list_.at(ide).link_surface(InterfaceDirection::D, surface_list_.at(ids).id());
				surface_list_.at(ids).link_element(ElementDirection::inc, element_list_.at(ide).id());
			}

	size_t ide = 0;
	size_t idf = 0;
	for (size_t ix = 0; ix < nx_; ix++)
		for (size_t iy = 0; iy < ny_; iy++)
			for (size_t iz = 0; iz < nz_; iz++)
			{
				ide = index3[ix][iy][iz];
				if (ix == nx_ - 1)
				{
					num_t xce = element_list_.at(ide).xc();
					num_t yce = element_list_.at(ide).yc();
					num_t zce = element_list_.at(ide).zc();
					num_t he = element_list_.at(ide).h();
					size_t tid = add_surface(0, SurfaceDirection::X, xce+he/2.0, yce, zce, he);
					element_list_.at(ide).link_surface(InterfaceDirection::F, surface_list_.at(tid).id());
					surface_list_.at(tid).link_element(ElementDirection::pre, element_list_.at(ide).id());
				}
				else
				{
					idf = index3[ix + 1][iy][iz];
					element_list_.at(ide).link_surface(InterfaceDirection::F, element_list_.at(idf).surface_adjacent(InterfaceDirection::B));
					surface_list_.at(element_list_.at(idf).surface_adjacent(InterfaceDirection::B)).link_element(ElementDirection::pre, element_list_.at(ide).id());
				}
				if (iy == ny_ - 1)
				{
					num_t xce = element_list_.at(ide).xc();
					num_t yce = element_list_.at(ide).yc();
					num_t zce = element_list_.at(ide).zc();
					num_t he = element_list_.at(ide).h();
					size_t tid = add_surface(0, SurfaceDirection::Y, xce, yce+he/2.0, zce, he);
					element_list_.at(ide).link_surface(InterfaceDirection::R, surface_list_.at(tid).id());
					surface_list_.at(tid).link_element(ElementDirection::pre, element_list_.at(ide).id());
				}
				else
				{
					idf = index3[ix][iy + 1][iz];
					element_list_.at(ide).link_surface(InterfaceDirection::R, element_list_.at(idf).surface_adjacent(InterfaceDirection::L));
					surface_list_.at(element_list_.at(idf).surface_adjacent(InterfaceDirection::L)).link_element(ElementDirection::pre, element_list_.at(ide).id());
				}
				if (iz == (nz_ - 1))
				{
					num_t xce = element_list_.at(ide).xc();
					num_t yce = element_list_.at(ide).yc();
					num_t zce = element_list_.at(ide).zc();
					num_t he = element_list_.at(ide).h();
					size_t tid = add_surface(0, SurfaceDirection::Z, xce, yce, zce+he/2.0, he);
					element_list_.at(ide).link_surface(InterfaceDirection::U, surface_list_.at(tid).id());
					surface_list_.at(tid).link_element(ElementDirection::pre, element_list_.at(ide).id());
				}
				else
				{
					idf = index3[ix][iy][iz + 1];
					element_list_.at(ide).link_surface(InterfaceDirection::U, element_list_.at(idf).surface_adjacent(InterfaceDirection::D));					
					surface_list_.at(element_list_.at(idf).surface_adjacent(InterfaceDirection::D)).link_element(ElementDirection::pre, element_list_.at(ide).id());
				}
			}

	for (size_t ix = 0; ix < nx_; ix++)
	{
		for (size_t iy = 0; iy < ny_; iy++)
		{
			surface_list_.at(element_list_.at(index3[ix][iy][0]).surface_adjacent(InterfaceDirection::D)).mark_as_boundary();
			surface_list_.at(element_list_.at(index3[ix][iy][nz_ - 1]).surface_adjacent(InterfaceDirection::U)).mark_as_boundary();
		}
	}
	for (size_t ix = 0; ix < nx_; ix++)
	{
		for (size_t iz = 0; iz < nz_; iz++)
		{
			surface_list_.at(element_list_.at(index3[ix][0][iz]).surface_adjacent(InterfaceDirection::L)).mark_as_boundary();
			surface_list_.at(element_list_.at(index3[ix][ny_ - 1][iz]).surface_adjacent(InterfaceDirection::R)).mark_as_boundary();
		}
	}
	for (size_t iz = 0; iz < nz_; iz++)
	{
		for (size_t iy = 0; iy < ny_; iy++)
		{
			surface_list_.at(element_list_.at(index3[0][iy][iz]).surface_adjacent(InterfaceDirection::B)).mark_as_boundary();
			surface_list_.at(element_list_.at(index3[nx_ - 1][iy][iz]).surface_adjacent(InterfaceDirection::F)).mark_as_boundary();
		}
	}
}

Element::Element(size_t eId, size_t refine_level, num_t xc, num_t yc, num_t zc, num_t h) 
	: id_(eId), refine_level_(refine_level), xc_(xc), yc_(yc), zc_(zc), h_(h)
{
	surface_list_.resize(kInterfaceTotal.at(3));
	for each (auto sid in surface_list_)
		sid = -1;
	child_elements.clear();
	surface_local_id.resize(kInterfaceTotal.at(3));
	for each (auto lid in surface_local_id)
		lid = 0;
	refine_mark_ = false;
}

Element::~Element()
{
//#ifdef _DEBUG
//	std::cerr << "Element deconstruction called!" << std::endl;
//#endif
}

void Element::link_surface(const InterfaceDirection direction, const size_t sid, const size_t local_id /*= 0*/)
{
	surface_list_.at(static_cast<size_t>(direction)) = sid;
}

const size_t Element::surface_adjacent(InterfaceDirection dir) const
{
	return surface_list_.at(static_cast<size_t>(dir));
}

Surface::Surface(size_t id, size_t refine_level, SurfaceDirection direction, num_t xc, num_t yc, num_t zc, num_t h) 
	: id_(id), refine_level_(refine_level), direction_(direction), is_boundary_(false), 
	  xc_(xc), yc_(yc), zc_(zc), h_(h)
{
}

Surface::~Surface()
{
//#ifdef _DEBUG
//	std::cerr << "Surface deconstruction called!" << std::endl;
//#endif
}

void Surface::link_element(const ElementDirection ed, const size_t eid, const size_t local_id /*= 0*/)
{
	std::vector<size_t>& el = (ed == ElementDirection::pre ? element_pre_list_ : element_inc_list_);
	if (el.size() < local_id + 1)
		el.resize(local_id + 1);
	el.at(local_id) = eid;
}

#include <iomanip>

std::ostream& operator<<(std::ostream& os, const Mesh& m)
{
	os << std::setprecision(5);
	os << "information of Mesh" << std::endl;
	os << "total elements = " << m.elements_total() << std::endl;
	os << "total interface = " << m.surfaces_total() << std::endl;
	os << "========   elements info:  ================" << std::endl;
	os << "id\txc\tyc\tzc" << std::endl;
	std::vector<vector3> ec = m.center();
	for (size_t i = 0; i < m.elements_total(); i++)
	{
		os << i << "\t" << ec.at(i).x() << "\t" << ec.at(i).y() << "\t" << ec.at(i).z() << std::endl;
	}
	os << "===========================================" << std::endl;

	std::vector<vector3> sc = m.center_surface();
	os << "========   Surface info:   ================" << std::endl;
	os << "id\tbd_flag\tdiret\txc\tyc\tzc" << std::endl;
	std::vector<bool> bdf = m.boundary_mark(1);
	for (size_t i = 0; i < m.surfaces_total(); i++)
	{
		os << i << "\t" << bdf.at(i) << "\t" << static_cast<size_t>(m.surface_ref(i).direction()) << "\t" << sc.at(i).x() << "\t" << sc.at(i).y() << "\t" << sc.at(i).z() << std::endl;
	}
	os << "===========================================" << std::endl;

	os << "========   Element Surface link:   ========" << std::endl;
	boost::multi_array<size_t, 2> esam = m.element_interface_matrix();
	os << "id\tB\tF\tL\tR\tD\tU" << std::endl;
	for (size_t i = 0; i < m.elements_total(); i++)
	{
		os << i;
		for (size_t j = 0; j < kInterfaceTotal.at(3); j++)
		{
			os << "\t" << esam[i][j];
		}
		os << std::endl;
	}
	os << "===========================================" << std::endl;

	os << "========== sweep order ====================" << std::endl;
	for (size_t i = 0; i < kQuadrantTotal.at(3); i++)
	{
		quadrant quad = quadrant_list.at(i);
		std::vector<size_t> order = m.sweep_order(quad);
		os << "quadrant " << i << std::endl;
		for (size_t i = 0; i < order.size(); i++)
		{
			os << order.at(i) << "\t";
		}
		os << std::endl;
	}
	os << "===========================================" << std::endl;
	
	return os;
}

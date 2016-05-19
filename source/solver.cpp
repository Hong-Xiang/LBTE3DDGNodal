#include "solver.h"
//========	source_generator ================================
dgn::source_generator::source_generator(const system_matrix_angle& s, const mesh& m)
	: s_(&s), m_(&m)
{
	size_t basis_total = s.basis_element();
	size_t node_total = m.memory_element_total(basis_total);
	xp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
	yp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
	zp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
	x_.bind(node_total, xp_);
	y_.bind(node_total, yp_);
	z_.bind(node_total, zp_);
	calculate_coordinates();
}

dgn::source_generator::~source_generator()
{
	if (xp_ != nullptr)
		mkl_free(xp_);
	if (yp_ != nullptr)
		mkl_free(yp_);
	if (zp_ != nullptr)
		mkl_free(zp_);
}

void dgn::source_generator::calculate(vector_t& data) const
{
	
	if (data.size() != x_.size())
		throw(invalid_dimension());	
	for (size_t i = 0; i < data.size(); i++)
	{
		data(i) = analytical_solution::source(x_(i), y_(i), z_(i), s_->mu(), s_->xi(), s_->eta(), s_->sigma());
	}
}

void dgn::source_generator::calculate_coordinates()
{
	size_t basis_total = s_->basis_element();
	size_t node_total = m_->memory_element_total(basis_total);
	std::vector<size_t> mempos = m_->memory_element(basis_total);
	std::vector<vector3> center = m_->center();
	vector_t x_std = s_->x();
	vector_t y_std = s_->y();
	vector_t z_std = s_->z();
	vector_t x_loc, y_loc, z_loc;
	for (size_t i = 0; i < mempos.size(); i++)
	{
		x_loc.bind(basis_total, xp_ + mempos.at(i));
		y_loc.bind(basis_total, yp_ + mempos.at(i));
		z_loc.bind(basis_total, zp_ + mempos.at(i));
		x_loc.scal_add(1.0, x_std, center.at(i).x());
		y_loc.scal_add(1.0, y_std, center.at(i).y());
		z_loc.scal_add(1.0, z_std, center.at(i).z());
	}
}
//============================================================


//========	boundary_generator ===============================

dgn::boundary_generator::boundary_generator(const system_matrix_angle& s, const mesh& m)
	: s_(&s), m_(&m)
{
	size_t szb = s.basis_surface();
	size_t node_total = m.memory_surface_total(szb);
	xp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
	yp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
	zp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
	x_.bind(node_total, xp_);
	y_.bind(node_total, yp_);
	z_.bind(node_total, zp_);
	calculate_coordinates();
}

dgn::boundary_generator::~boundary_generator()
{
	if (xp_ != nullptr)
		mkl_free(xp_);
	if (yp_ != nullptr)
		mkl_free(yp_);
	if (zp_ != nullptr)
		mkl_free(zp_);
}

void dgn::boundary_generator::calculate(vector_t data) const
{
	if (data.size() != x_.size())
		throw(invalid_dimension());
	for (size_t i = 0; i < data.size(); i++)
	{
		data(i) = analytical_solution::solution(x_(i), y_(i), z_(i), s_->mu(), s_->xi(), s_->eta(), s_->sigma());
	}
}

void dgn::boundary_generator::calculate_coordinates()
{
	size_t basis_total = s_->basis_surface();
	size_t node_total = m_->memory_surface_total(basis_total);
	std::vector<size_t> mempos = m_->memory_surface(basis_total);
	std::vector<vector3> center = m_->center_surface();
	vector_t x_std = s_->xs();
	vector_t y_std = s_->ys();	
	vector_t x_loc, y_loc, z_loc;
	for (size_t i = 0; i < mempos.size(); i++)
	{
		x_loc.bind(basis_total, xp_ + mempos.at(i));
		y_loc.bind(basis_total, yp_ + mempos.at(i));
		z_loc.bind(basis_total, zp_ + mempos.at(i));

		//case X
		if (m_->get_surface_direction(i) == surface_direction::X)
		{
			vector3 c = center.at(i);
			for (size_t j = 0; j < x_loc.size(); j++)
			{
				x_loc(j) = c.x();
			}			
			y_loc.scal_add(1.0, x_std, center.at(i).y());
			z_loc.scal_add(1.0, y_std, center.at(i).z());
		}

		//case Y
		if (m_->get_surface_direction(i) == surface_direction::Y)
		{
			vector3 c = center.at(i);
			for (size_t j = 0; j < x_loc.size(); j++)
			{
				y_loc(j) = c.y();
			}
			x_loc.scal_add(1.0, x_std, center.at(i).x());
			z_loc.scal_add(1.0, y_std, center.at(i).z());
		}

		//case Z
		if (m_->get_surface_direction(i) == surface_direction::Z)
		{
			vector3 c = center.at(i);
			for (size_t j = 0; j < x_loc.size(); j++)
			{
				z_loc(j) = c.z();
			}
			x_loc.scal_add(1.0, x_std, center.at(i).x());
			y_loc.scal_add(1.0, y_std, center.at(i).y());
		}
	}
}

//============================================================

//========	solver ===========================================

dgn::solver::solver(const system_matrix_angle& s, const mesh& m) : s_(&s), m_(&m), sg_(s, m), bg_(s, m)
{

	source_p_ = xlib::mkl_ext::xcalloc<num_t>(m.memory_element_total(s.basis_element()));
	boundary_p_ = xlib::mkl_ext::xcalloc<num_t>(m.memory_surface_total(s.basis_surface()));
	solution_p_ = xlib::mkl_ext::xcalloc<num_t>(m.memory_element_total(s.basis_element()));
	source_.bind(m.memory_element_total(s.basis_element()), source_p_);
	boundary_.bind(m.memory_surface_total(s.basis_surface()), boundary_p_);
	solution_.bind(m.memory_element_total(s.basis_element()), solution_p_);

	//boundary mask	
	boundary_mask_p_ = xlib::mkl_ext::xcalloc<num_t>(m.memory_surface_total(s.basis_surface()));
	boundary_mask_.bind(m.memory_surface_total(s.basis_surface()), boundary_mask_p_);

	std::vector<bool> bdm = m_->boundary_mark(s_->basis_surface());
	if (bdm.size() != boundary_mask_.size())
		throw(invalid_dimension());
	for (size_t i = 0; i < boundary_.size(); i++)
	{
		if (bdm.at(i))
			boundary_mask_(i) = 0.0;
		else
			boundary_mask_(i) = 1.0;
	}
}

dgn::solver::~solver()
{
	if (source_p_ != nullptr)
		mkl_free(source_p_);
	if (solution_p_ != nullptr)
		mkl_free(solution_p_);
	if (boundary_p_ != nullptr)
		mkl_free(boundary_p_);
	if (boundary_mask_p_ != nullptr)
		mkl_free(boundary_mask_p_);
}

void dgn::solver::solve()
{
	std::vector<size_t> sweep_order;
	boost::multi_array<size_t, 2> eim;
	num_t mu, xi, eta;
	mu = s_->mu(); xi = s_->xi(); eta = s_->eta();
	sweep_order = m_->sweep_order(utilities::quadrant_of_angle(mu, xi, eta));

	eim.resize(boost::extents[m_->elements_total()][6]);	
	eim = m_->element_interface_matrix();

	
	size_t n_node_element = s_->basis_element();
	size_t n_node_surface = s_->basis_surface();

	

	
	
	std::vector<size_t> memory_element = m_->memory_element(n_node_element);
	std::vector<size_t> memory_surface = m_->memory_surface(n_node_surface);

	vector_t s_e;
	vector_t sor_e;
	vector_t bd_xpre, bd_xinc;
	vector_t bd_ypre, bd_yinc;
	vector_t bd_zpre, bd_zinc;
	vector_t bdmask_xinc, bdmask_yinc, bdmask_zinc;
	size_t ib_xpre, ib_xinc, ib_ypre, ib_yinc, ib_zpre, ib_zinc;
	interface_direction d_xpre, d_xinc, d_ypre, d_yinc, d_zpre, d_zinc;
	num_p bp = xlib::mkl_ext::xcalloc<num_t>(n_node_element);
	vector_t b(s_->basis_element(), bp);

	num_p bdtmpp = xlib::mkl_ext::xcalloc<num_t>(n_node_surface);
	vector_t bdtmp(n_node_surface, bdtmpp);

	quadrant quad = utilities::quadrant_of_angle(mu, xi, eta);
	bool xincf = utilities::quadrant_x_positive_flag(quad);
	bool yincf = utilities::quadrant_y_positive_flag(quad);
	bool zincf = utilities::quadrant_z_positive_flag(quad);

	for (size_t ids = 0; ids < sweep_order.size(); ids++)
	{
		// Binding vectors		
		size_t i = sweep_order.at(ids);
		d_xpre = xincf	? interface_direction::B
						: interface_direction::F;
		d_xinc = xincf	? interface_direction::F
						: interface_direction::B;
		d_ypre = yincf	? interface_direction::L
						: interface_direction::R;
		d_yinc = yincf	? interface_direction::R
						: interface_direction::L;
		d_zpre = zincf	? interface_direction::D
						: interface_direction::U;
		d_zinc = zincf	? interface_direction::U
						: interface_direction::D;

		ib_xpre = eim[i][static_cast<size_t>(d_xpre)];
		ib_xinc = eim[i][static_cast<size_t>(d_xinc)];
		ib_ypre = eim[i][static_cast<size_t>(d_ypre)];
		ib_yinc = eim[i][static_cast<size_t>(d_yinc)];
		ib_zpre = eim[i][static_cast<size_t>(d_zpre)];
		ib_zinc = eim[i][static_cast<size_t>(d_zinc)];

		bd_xpre.bind(n_node_surface, boundary_p_ + memory_surface.at(ib_xpre));		
		bd_xinc.bind(n_node_surface, boundary_p_ + memory_surface.at(ib_xinc));
		bd_ypre.bind(n_node_surface, boundary_p_ + memory_surface.at(ib_ypre));
		bd_yinc.bind(n_node_surface, boundary_p_ + memory_surface.at(ib_yinc));
		bd_zpre.bind(n_node_surface, boundary_p_ + memory_surface.at(ib_zpre));
		bd_zinc.bind(n_node_surface, boundary_p_ + memory_surface.at(ib_zinc));

		bdmask_xinc.bind(n_node_surface, boundary_mask_p_ + memory_surface.at(ib_xinc));
		bdmask_yinc.bind(n_node_surface, boundary_mask_p_ + memory_surface.at(ib_yinc));
		bdmask_zinc.bind(n_node_surface, boundary_mask_p_ + memory_surface.at(ib_zinc));

		sor_e.bind(n_node_element, source_.ptr() + memory_element.at(i));
		s_e.bind(n_node_element, solution_p_ + memory_element.at(i));
		

		// Calculate b

	
		
			b.gemx(CblasNoTrans, 1.0, s_->mass_matrix(), sor_e, 0.0);


			b.gemx(CblasNoTrans, 1.0, s_->lift_matrix(d_xpre), bd_xpre, 1.0);

	

			b.gemx(CblasNoTrans, 1.0, s_->lift_matrix(d_ypre), bd_ypre, 1.0);


			b.gemx(CblasNoTrans, 1.0, s_->lift_matrix(d_zpre), bd_zpre, 1.0);						


		// Solve matrix
			s_e.copy(b);
			s_e.gbtrs('N', s_->system_matrixlu());
			//s_e.gemx(CblasNoTrans, 1.0, s_->system_matrixiv(), b, 0.0);
		// Update boundary
			bdtmp.gemx(CblasTrans, 1.0, s_->flux_matrix(d_xinc), s_e, 0.0);		
			bdtmp.point_mul(bdtmp, bdmask_xinc);

			
			bd_xinc.add(bdtmp);
			

			bdtmp.gemx(CblasTrans, 1.0, s_->flux_matrix(d_yinc), s_e, 0.0);
			bdtmp.point_mul(bdtmp, bdmask_yinc);
			
			bd_yinc.add(bdtmp);
			
			bdtmp.gemx(CblasTrans, 1.0, s_->flux_matrix(d_zinc), s_e, 0.0);
			bdtmp.point_mul(bdtmp, bdmask_zinc);
			bd_zinc.add(bdtmp);

			//bd_xinc.gemx(CblasTrans, 1.0, s_->flux_matrix(d_xinc), s_e, 1.0);
			//bd_yinc.gemx(CblasTrans, 1.0, s_->flux_matrix(d_yinc), s_e, 1.0);
			//bd_zinc.gemx(CblasTrans, 1.0, s_->flux_matrix(d_zinc), s_e, 1.0);
	}
	if (bp != nullptr)
		mkl_free(bp);
	if (bdtmpp != nullptr)
		mkl_free(bdtmpp);
}

#include "globals.h"

void dgn::solver::solve_test(size_t id)
{
	analytical_solution::set_test_id(id);
	sg_.calculate(source_);
	bg_.calculate(boundary_);
	std::vector<bool> bdm = m_->boundary_mark(s_->basis_surface());
	if (bdm.size() != boundary_.size())
		throw(invalid_dimension());
	for (size_t i = 0; i < boundary_.size(); i++)
	{
		if (bdm.at(i) == false)
			boundary_(i) = 0.0;
	}

	solve();

}

void dgn::solver::solve_analytical(size_t id)
{
	analytical_solution::set_test_id(id);		
	vector_t x = sg_.x();
	vector_t y = sg_.y();
	vector_t z = sg_.z();
	for (size_t i = 0; i < solution_.size(); i++)
	{
		solution_(i) = analytical_solution::solution(x(i), y(i), z(i), s_->mu(), s_->xi(), s_->eta(), s_->sigma());
	}
}

void dgn::solver::solution_mean(vector_t data) const
{
	size_t n_element = m_->elements_total();
	vector_t v_e;
	std::vector<size_t> mempos = m_->memory_element(s_->basis_element());
	for (size_t i = 0; i < n_element; i++)
	{
		v_e.bind(s_->basis_element(), solution_p_ + mempos.at(i));
		data(i) = v_e.sum()/s_->basis_element();
	}
}

//============================================================
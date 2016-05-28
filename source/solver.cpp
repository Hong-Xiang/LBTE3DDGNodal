#include "solver.h"
#include <ctime>

////========	source_generator ================================
//source_generator::source_generator(const system_matrix_angle& s, const Mesh& m)
//	: s_(&s), m_(&m)
//{
//	size_t basis_total = s.basis_element();
//	size_t node_total = m.memory_element_total(basis_total);
//	xp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
//	yp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
//	zp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
//	x_.Bind(node_total, xp_);
//	y_.Bind(node_total, yp_);
//	z_.Bind(node_total, zp_);
//	calculate_coordinates();
//}
//
//source_generator::~source_generator()
//{
//	if (xp_ != nullptr)
//		mkl_free(xp_);
//	if (yp_ != nullptr)
//		mkl_free(yp_);
//	if (zp_ != nullptr)
//		mkl_free(zp_);
//}
//
//void source_generator::calculate(Vector& data) const
//{
//
//	if (data.size() != x_.size())
//		throw(invalid_dimension());
//	for (size_t ide = 0; ide < data.size(); ide++)
//	{
//		data(ide) = analytical_solution::source(x_(ide), y_(ide), z_(ide), s_->mu(), s_->xi(), s_->eta(), s_->sigma());
//	}
//}
//
//void source_generator::calculate_coordinates()
//{
//	size_t basis_total = s_->basis_element();
//	size_t node_total = m_->memory_element_total(basis_total);
//	std::vector<size_t> mempos = m_->memory_element(basis_total);
//	std::vector<vector3> center = m_->center();
//	Vector x_std = s_->x();
//	Vector y_std = s_->y();
//	Vector z_std = s_->z();
//	Vector x_loc, y_loc, z_loc;
//	for (size_t ide = 0; ide < mempos.size(); ide++)
//	{
//		x_loc.Bind(basis_total, xp_ + mempos.at(ide));
//		y_loc.Bind(basis_total, yp_ + mempos.at(ide));
//		z_loc.Bind(basis_total, zp_ + mempos.at(ide));
//		x_loc.scal_add(1.0, x_std, center.at(ide).x());
//		y_loc.scal_add(1.0, y_std, center.at(ide).y());
//		z_loc.scal_add(1.0, z_std, center.at(ide).z());
//	}
//}
////============================================================
//
//
////========	boundary_generator ===============================
//
//boundary_generator::boundary_generator(const system_matrix_angle& s, const Mesh& m)
//	: s_(&s), m_(&m)
//{
//	size_t szb = s.basis_surface();
//	size_t node_total = m.memory_surface_total(szb);
//	xp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
//	yp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
//	zp_ = xlib::mkl_ext::xcalloc<num_t>(node_total);
//	x_.Bind(node_total, xp_);
//	y_.Bind(node_total, yp_);
//	z_.Bind(node_total, zp_);
//	calculate_coordinates();
//}
//
//boundary_generator::~boundary_generator()
//{
//	if (xp_ != nullptr)
//		mkl_free(xp_);
//	if (yp_ != nullptr)
//		mkl_free(yp_);
//	if (zp_ != nullptr)
//		mkl_free(zp_);
//}
//
//void boundary_generator::calculate(Vector data) const
//{
//	if (data.size() != x_.size())
//		throw(invalid_dimension());
//	for (size_t ide = 0; ide < data.size(); ide++)
//	{
//		data(ide) = analytical_solution::solution(x_(ide), y_(ide), z_(ide), s_->mu(), s_->xi(), s_->eta(), s_->sigma());
//	}
//}
//
//void boundary_generator::calculate_coordinates()
//{
//	size_t basis_total = s_->basis_surface();
//	size_t node_total = m_->memory_surface_total(basis_total);
//	std::vector<size_t> mempos = m_->memory_surface(basis_total);
//	std::vector<vector3> center = m_->center_surface();
//	Vector x_std = s_->xs();
//	Vector y_std = s_->ys();
//	Vector x_loc, y_loc, z_loc;
//	for (size_t ide = 0; ide < mempos.size(); ide++)
//	{
//		x_loc.Bind(basis_total, xp_ + mempos.at(ide));
//		y_loc.Bind(basis_total, yp_ + mempos.at(ide));
//		z_loc.Bind(basis_total, zp_ + mempos.at(ide));
//
//		//case X
//		if (m_->get_SurfaceDirection(ide) == SurfaceDirection::X)
//		{
//			vector3 c = center.at(ide);
//			for (size_t j = 0; j < x_loc.size(); j++)
//			{
//				x_loc(j) = c.x();
//			}
//			y_loc.scal_add(1.0, x_std, center.at(ide).y());
//			z_loc.scal_add(1.0, y_std, center.at(ide).z());
//		}
//
//		//case Y
//		if (m_->get_SurfaceDirection(ide) == SurfaceDirection::Y)
//		{
//			vector3 c = center.at(ide);
//			for (size_t j = 0; j < x_loc.size(); j++)
//			{
//				y_loc(j) = c.y();
//			}
//			x_loc.scal_add(1.0, x_std, center.at(ide).x());
//			z_loc.scal_add(1.0, y_std, center.at(ide).z());
//		}
//
//		//case Z
//		if (m_->get_SurfaceDirection(ide) == SurfaceDirection::Z)
//		{
//			vector3 c = center.at(ide);
//			for (size_t j = 0; j < x_loc.size(); j++)
//			{
//				z_loc(j) = c.z();
//			}
//			x_loc.scal_add(1.0, x_std, center.at(ide).x());
//			y_loc.scal_add(1.0, y_std, center.at(ide).y());
//		}
//	}
//}
//
////============================================================
//
//const size_t solver::total_node_element(const Mesh& m, size_t np)
//{
//	return m.memory_element_total(system_matrix_angle::basis_total(np));
//}
//
//
//
//const size_t solver::total_node_surface() const
//{
//	return total_node_surface(*m_, s_->np());
//}
//
//const size_t solver::total_node_surface(const Mesh& m, size_t np)
//{
//	return m.memory_surface_total(system_matrix_angle::basis_lower_dim_total(np));
//}
//
//const size_t solver::total_node_element() const
//{
//	return total_node_element(*m_, s_->np());
//}
//
////========	solver ===========================================
//
//solver::solver(const system_matrix_angle& s, const Mesh& m) : s_(&s), m_(&m), sg_(s, m), bg_(s, m)
//{
//
//	source_p_ = xlib::mkl_ext::xcalloc<num_t>(m.memory_element_total(s.basis_element()));
//	boundary_p_ = xlib::mkl_ext::xcalloc<num_t>(m.memory_surface_total(s.basis_surface()));
//	solution_p_ = xlib::mkl_ext::xcalloc<num_t>(m.memory_element_total(s.basis_element()));
//	source_.Bind(m.memory_element_total(s.basis_element()), source_p_);
//	boundary_.Bind(m.memory_surface_total(s.basis_surface()), boundary_p_);
//	solution_.Bind(m.memory_element_total(s.basis_element()), solution_p_);
//
//	//boundary mask	
//	boundary_mask_p = xlib::mkl_ext::xcalloc<num_t>(m.memory_surface_total(s.basis_surface()));
//	boundary_mask_.Bind(m.memory_surface_total(s.basis_surface()), boundary_mask_p);
//
//	std::vector<bool> bdm = m_->boundary_mark(s_->basis_surface());
//	if (bdm.size() != boundary_mask_.size())
//		throw(invalid_dimension());
//	for (size_t ide = 0; ide < boundary_.size(); ide++)
//	{
//		if (bdm.at(ide))
//			boundary_mask_(ide) = 0.0;
//		else
//			boundary_mask_(ide) = 1.0;
//	}
//}
//
//solver::~solver()
//{
//	if (source_p_ != nullptr)
//		mkl_free(source_p_);
//	if (solution_p_ != nullptr)
//		mkl_free(solution_p_);
//	if (boundary_p_ != nullptr)
//		mkl_free(boundary_p_);
//	if (boundary_mask_p != nullptr)
//		mkl_free(boundary_mask_p);
//}
//
//void solver::solve()
//{
//	std::vector<bool> bdm = m_->boundary_mark(s_->basis_surface());
//	if (bdm.size() != boundary_.size())
//		throw(invalid_dimension());
//	for (size_t ide = 0; ide < boundary_.size(); ide++)
//	{
//		if (bdm.at(ide) == false)
//			boundary_(ide) = 0.0;
//	}
//
//	std::vector<size_t> sweep_order;
//	boost::multi_array<size_t, 2> esm;
//	num_t mu, xi, eta;
//	mu = s_->mu(); xi = s_->xi(); eta = s_->eta();
//	sweep_order = m_->sweep_order(Utilities::quadrant_of_angle(mu, xi, eta));
//
//	esm.resize(boost::extents[m_->elements_total()][6]);
//	esm = m_->element_interface_matrix();
//
//
//	size_t n_node_element = s_->basis_element();
//	size_t n_node_surface = s_->basis_surface();
//
//
//
//
//
//	std::vector<size_t> memory_element = m_->memory_element(n_node_element);
//	std::vector<size_t> memory_surface = m_->memory_surface(n_node_surface);
//
//	Vector s_e;
//	Vector sor_e;
//	Vector bd_xpre, bd_xinc;
//	Vector bd_ypre, bd_yinc;
//	Vector bd_zpre, bd_zinc;
//	Vector bdmask_xinc, bdmask_yinc, bdmask_zinc;
//	size_t ib_xpre, ib_xinc, ib_ypre, ib_yinc, ib_zpre, ib_zinc;
//	SurfaceDirection d_xpre, d_xinc, d_ypre, d_yinc, d_zpre, d_zinc;
//	num_p bp = xlib::mkl_ext::xcalloc<num_t>(n_node_element);
//	Vector b(s_->basis_element(), bp);
//	num_p btmpp = xlib::mkl_ext::xcalloc<num_t>(n_node_element);
//	Vector btmp(s_->basis_element(), btmpp);
//
//	num_p bdtmpp = xlib::mkl_ext::xcalloc<num_t>(n_node_surface);
//	Vector bdtmp(n_node_surface, bdtmpp);
//	num_p bdtmpp2 = xlib::mkl_ext::xcalloc<num_t>(n_node_surface);
//	Vector bdtmp2(n_node_surface, bdtmpp);
//
//
//	Quadrant quad = Utilities::quadrant_of_angle(mu, xi, eta);
//	bool xincf = Utilities::quadrant_x_positive_flag(quad);
//	bool yincf = Utilities::quadrant_y_positive_flag(quad);
//	bool zincf = Utilities::quadrant_z_positive_flag(quad);
//
////#ifdef _DEBUG
////	num_p solution_ana_p = xlib::mkl_ext::xcalloc<num_t>(m_->memory_element_total(s_->basis_element()));
////	num_p boundary_ana_p = xlib::mkl_ext::xcalloc<num_t>(m_->memory_surface_total(s_->basis_surface()));
////	Vector solution_ana(m_->memory_element_total(s_->basis_element()), solution_ana_p);
////	Vector boundary_ana(m_->memory_element_total(s_->basis_surface()), boundary_ana_p);
////	Vector xe = sg_.x(), ye = sg_.y(), ze = sg_.z();
////	for (size_t ide = 0; ide < m_->memory_element_total(s_->basis_element()); ide++)
////	{
////		solution_ana(ide) = analytical_solution::solution(xe(ide), ye(ide), ze(ide), mu, xi, eta, s_->sigma());
////	}
////	Vector xb = bg_.x(), yb = bg_.y(), zb = bg_.z();
////	for (size_t ide = 0; ide < m_->memory_surface_total(s_->basis_surface()); ide++)
////	{
////		boundary_ana(ide) = analytical_solution::solution(xb(ide), yb(ide), zb(ide), mu, xi, eta, s_->sigma());
////	}
////
////	Vector bda_xpre, bda_ypre, bda_zpre, bda_xinc, bda_yinc, bda_zinc, sa;
////#endif
//
//
//	//Reorder.
//	std::vector<size_t> pos;
//	pos.clear();
//	pos.resize(n_node_element);
//	for (size_t ide = 0; ide < n_node_element; ide++)
//	{
//		size_t ix, iy, iz, np;
//		np = s_->np();
//		xlib::ind2sub(np, np, np, ide, ix, iy, iz);
//		if (mu < 0)
//			ix = np - 1 - ix;
//		if (xi < 0)
//			iy = np - 1 - iy;
//		if (eta < 0)
//			iz = np - 1 - iz;
//		size_t nid;
//		xlib::sub2ind(np, np, np, nid, ix, iy, iz);
//		pos.at(ide) = nid;
//	}
//
//
//	num_p system_matrix_inverse_ordered_p = xlib::mkl_ext::xcalloc<num_t>(n_node_element*n_node_element);
//	matrix_t system_matrix_inverse_ordered(n_node_element, n_node_element, system_matrix_inverse_ordered_p);
//
//	for (size_t ide = 0; ide < n_node_element; ide++)
//	{
//		for (size_t j = 0; j < n_node_element; j++)
//		{
//			//system_matrix_inverse_ordered(ide, j) = s_->system_matrixiv()(ide, pos.at(j));
//			system_matrix_inverse_ordered(ide, j) = s_->system_matrixiv()(ide, j);
//		}
//	}
//
//	//xlib::imatlab::load_MAT_data(system_matrix_inverse_ordered_p, "invmat.mat", "invSysM");
//	//system_matrix_inverse_ordered.gemm(s_->system_matrix(), s_->system_matrixiv(), 1.0, 0.0);
//	//std::cout << system_matrix_inverse_ordered;
//	//num_p order_matrix_p = xlib::mkl_ext::xcalloc<num_t>(n_node_element*n_node_element);
//	//matrix_t order_matrix(n_node_element, n_node_element, order_matrix_p);
//	for (size_t ids = 0; ids < sweep_order.size(); ids++)
//	{
//		// Binding vectors		
//		size_t ide = sweep_order.at(ids);
//		d_xpre = xincf ? SurfaceDirection::B
//			: SurfaceDirection::F;
//		d_xinc = xincf ? SurfaceDirection::F
//			: SurfaceDirection::B;
//		d_ypre = yincf ? SurfaceDirection::L
//			: SurfaceDirection::R;
//		d_yinc = yincf ? SurfaceDirection::R
//			: SurfaceDirection::L;
//		d_zpre = zincf ? SurfaceDirection::D
//			: SurfaceDirection::U;
//		d_zinc = zincf ? SurfaceDirection::U
//			: SurfaceDirection::D;
//
//		ib_xpre = esm[ide][Utilities::id(d_xpre)];
//		ib_xinc = esm[ide][Utilities::id(d_xinc)];
//		ib_ypre = esm[ide][Utilities::id(d_ypre)];
//		ib_yinc = esm[ide][Utilities::id(d_yinc)];
//		ib_zpre = esm[ide][Utilities::id(d_zpre)];
//		ib_zinc = esm[ide][Utilities::id(d_zinc)];
//
//		bd_xpre.Bind(n_node_surface, boundary_p_ +ib_xpre));
//		bd_xinc.Bind(n_node_surface, boundary_p_ +ib_xinc));
//		bd_ypre.Bind(n_node_surface, boundary_p_ +ib_ypre));
//		bd_yinc.Bind(n_node_surface, boundary_p_ +ib_yinc));
//		bd_zpre.Bind(n_node_surface, boundary_p_ +ib_zpre));
//		bd_zinc.Bind(n_node_surface, boundary_p_ +ib_zinc));
//
//		bdmask_xinc.Bind(n_node_surface, boundary_mask_p +ib_xinc));
//		bdmask_yinc.Bind(n_node_surface, boundary_mask_p +ib_yinc));
//		bdmask_zinc.Bind(n_node_surface, boundary_mask_p +ib_zinc));
//
//		sor_e.Bind(n_node_element, source_.ptr() + memory_element.at(ide));
//		s_e.Bind(n_node_element, solution_p_ + memory_element.at(ide));
//
////#ifdef _DEBUG
////		bda_xpre.Bind(n_node_surface, boundary_ana_p +ib_xpre));
////		bda_ypre.Bind(n_node_surface, boundary_ana_p +ib_ypre));
////		bda_zpre.Bind(n_node_surface, boundary_ana_p +ib_zpre));
////		bda_xinc.Bind(n_node_surface, boundary_ana_p +ib_xinc));
////		bda_yinc.Bind(n_node_surface, boundary_ana_p +ib_yinc));
////		bda_zinc.Bind(n_node_surface, boundary_ana_p +ib_zinc));
////		sa.Bind(n_node_element, solution_ana_p + memory_element.at(ide));
////#endif
//		// Calculate b
//
//
//
//		b.gemx(CblasNoTrans, 1.0, s_->mass_matrix(), sor_e, 0.0);
//		b.gemx(CblasNoTrans, 1.0, s_->lift_matrix(d_xpre), bd_xpre, 1.0);
//		b.gemx(CblasNoTrans, 1.0, s_->lift_matrix(d_ypre), bd_ypre, 1.0);
//		b.gemx(CblasNoTrans, 1.0, s_->lift_matrix(d_zpre), bd_zpre, 1.0);
//
//
//		// Solve matrix
//		s_e.copy(b);
//
//		//b reorder
//		//std::vector<num_t> tmpb(n_node_element);
//		//for (size_t j = 0; j < n_node_element; j++)
//		//{
//		//	tmpb.at(pos.at(j)) = b(j);
//		//}
//		//for (size_t j = 0; j < n_node_element; j++)
//		//{
//		//	b(j) = tmpb.at(j);
//		//}
//		btmp.copy(b);
//		b.gemx(CblasNoTrans, 1.0, s_->reorder(), btmp, 0.0);
//		
//		//s_e.gemx(CblasNoTrans, 1.0, system_matrix_inverse_ordered, b);
//		//s_e.gbtrs('N', s_->system_matrixlu());
//
//
////s_e with reorder
//		s_e.gemx(CblasNoTrans, 1.0, system_matrix_inverse_ordered, b, 0.0);
//		//#ifdef _DEBUG
//		//		std::cout << system_matrix_inverse_ordered;
//		//		std::cout << b;
//		//		std::cout << s_e;
//		//#endif
//				//s_e reorder
//		//for (size_t j = 0; j < n_node_element; j++)
//		//{
//		//	tmpb.at(pos.at(j)) = s_e(j);
//		//}
//		//for (size_t j = 0; j < n_node_element; j++)
//		//{
//		//	s_e(j) = tmpb.at(j);
//		//}
//		btmp.copy(s_e);
//		s_e.gemx(CblasNoTrans, 1.0, s_->back_order(), btmp, 0.0);
//
//		//s_e.gemx(CblasNoTrans, 1.0, s_->system_matrixiv(), b, 0.0);
//		// Update boundary
//		bdtmp.gemx(CblasTrans, 1.0, s_->flux_matrix(d_xinc), s_e, 0.0);
//		bdtmp2.point_mul(bdtmp, bdmask_xinc);
//		bd_xinc.add(bdtmp2);
//
//		bdtmp.gemx(CblasTrans, 1.0, s_->flux_matrix(d_yinc), s_e, 0.0);
//		bdtmp2.point_mul(bdtmp, bdmask_yinc);
//		bd_yinc.add(bdtmp2);
//
//		bdtmp.gemx(CblasTrans, 1.0, s_->flux_matrix(d_zinc), s_e, 0.0);
//		bdtmp2.point_mul(bdtmp, bdmask_zinc);
//		bd_zinc.add(bdtmp2);
//
//#ifdef _DEBUG
//
//			//std::cout << "b" << std::endl << b;
//
//			//std::cout << "bda xpre" << std::endl << bda_xpre;
//			//std::cout << "bd xpre" << std::endl << bd_xpre;
//
//			//std::cout << "bda ypre" << std::endl << bda_ypre;
//			//std::cout << "bd ypre" << std::endl << bd_ypre;
//
//			//std::cout << "bda zpre" << std::endl << bda_zpre;
//			//std::cout << "bd zpre" << std::endl << bd_zpre;
//
//			//std::cout << "bda xinc" << std::endl << bda_xinc;
//			//std::cout << "bd xinc" << std::endl << bd_xinc;
//
//			//std::cout << "bda yinc" << std::endl << bda_yinc;
//			//std::cout << "bd yinc" << std::endl << bd_yinc;
//
//			//std::cout << "bda zinc" << std::endl << bda_zinc;
//			//std::cout << "bd zinc" << std::endl << bd_zinc;
//
//			//std::cout << "sa" << std::endl << sa;
//			//std::cout << "sn" << std::endl << s_e;
//		
//#endif
//		//bd_xinc.gemx(CblasTrans, 1.0, s_->flux_matrix(d_xinc), s_e, 0.0);
//		//bd_yinc.gemx(CblasTrans, 1.0, s_->flux_matrix(d_yinc), s_e, 0.0);
//		//bd_zinc.gemx(CblasTrans, 1.0, s_->flux_matrix(d_zinc), s_e, 0.0);
//	}
//	if (bp != nullptr)
//		mkl_free(bp);
//	if (bdtmpp != nullptr)
//		mkl_free(bdtmpp);
//	if (system_matrix_inverse_ordered_p != nullptr)
//		mkl_free(system_matrix_inverse_ordered_p);
//	if (btmpp != nullptr)
//		mkl_free(btmpp);
//}
//
//#include "globals.h"
//
//void solver::solve_test(size_t id)
//{
//	analytical_solution::set_test_id(id);
//	sg_.calculate(source_);
//	bg_.calculate(boundary_);
//	
//	solve();
//}
//
//void solver::solve_analytical(size_t id)
//{
//	analytical_solution::set_test_id(id);
//	Vector x = sg_.x();
//	Vector y = sg_.y();
//	Vector z = sg_.z();
//	for (size_t ide = 0; ide < solution_.size(); ide++)
//	{
//		solution_(ide) = analytical_solution::solution(x(ide), y(ide), z(ide), s_->mu(), s_->xi(), s_->eta(), s_->sigma());
//	}
//}
//
//void solver::solution_mean(Vector data) const
//{
//	size_t n_element = m_->elements_total();
//	Vector v_e;
//	std::vector<size_t> mempos = m_->memory_element(s_->basis_element());
//	for (size_t ide = 0; ide < n_element; ide++)
//	{
//		v_e.Bind(s_->basis_element(), solution_p_ + mempos.at(ide));
//		data(ide) = v_e.sum() / s_->basis_element();
//	}
//}
//
////============================================================


namespace discontinues_galerkin_nodal_solver {
  SolverBase::SolverBase() 
    : status_(SolverStatus::PreInit)
  {
  }
  void SolverBase::status_set(SolverStatus ss)
  {
    status_ = ss;
  }


  SolverSingleAngle::SolverSingleAngle(const MeshWithCoordinate & mesh_input, 
      std::string input_file_name)
    : mesh_(mesh_input)
  {
    problem_definition_.Load(input_file_name);
    is_test_ = !(problem_definition_.test_id()==kInfSize);
  }
  SolverSingleAngle::SolverSingleAngle(const MeshWithCoordinate & mesh_input, const Vector & source_ext, const Vector & boundary_ext, std::string input_file_name)
    : mesh_(mesh_input)
  {
    problem_definition_.Load(input_file_name);
    is_test_ = false;
    source_ext_ = source_ext;
    boundary_ext_ = boundary_ext;
  }
  void SolverSingleAngle::Procceed()
  {
    if (status() == SolverStatus::PreInit)
    {
      status_set(SolverStatus::PostInit);
      return;
    }
    if (status() == SolverStatus::PostInit)
    {
      status_set(SolverStatus::Closed);
      return;
    }
    if (status() == SolverStatus::Closed)
    {
      status_set(SolverStatus::OutputReady);
      return;
    }
    if (status() == SolverStatus::OutputReady)
    {
      status_set(SolverStatus::Closed);
      return;
    }
  }
  void SolverSingleAngle::Initialization()
  {
    ProblemDefinitionSingleAngle& pd = problem_definition_;
    system_ = std::make_shared<SystemMatrix>( pd.h(), pd.np(),
      pd.mu(), pd.xi(), pd.eta(), pd.sigma() );
    if (is_test_)
    {
      source_ = std::make_shared<AnalyticalSourceMeshCoordinates>(
          mesh_, *system_, pd.test_id());
      boundary_ = std::make_shared<AnalyticalBoundaryMeshCoordinate>(
          mesh_, *system_, pd.test_id());
    }else{
      source_ = std::make_shared<Source>(mesh_, *system_);
      boundary_ = std::make_shared<Boundary>(mesh_, *system_);
    }
    solution_ = std::make_shared<Solution>(mesh_, *system_);
  }

  void SolverSingleAngle::SetCondition()
  {
    if (is_test_)
    {
      std::static_pointer_cast<AnalyticalSourceMeshCoordinates>(source_)->CalculateValue();
      std::static_pointer_cast<AnalyticalBoundaryMeshCoordinate>(boundary_)->CalculateValue();
    }
    else {
      source_->Value().Copy(source_ext_);
      boundary_->Value().Copy(boundary_ext_);
    }
  }

  void SolverSingleAngle::Solve()
  {
    time_t t0, t1;
    t0 = std::clock();
    Quadrant quad = Utilities::quadrant_of_angle(mu(), xi(), eta());
    std::vector<size_t> sweep_order = mesh_.SweepOrder(quad);
    boost::multi_array<size_t, 2> esm = mesh_.ElementSurfaceMatrix();

    MKL_INT n_node_element = system_->n_basis_element();
    MKL_INT n_node_surface = system_->n_basis_surface();

    Vector sol_e;
    Vector sor_e;
    Vector bd_xpre, bd_xinc;
    Vector bd_ypre, bd_yinc;
    Vector bd_zpre, bd_zinc;
    Vector bdmask_xinc, bdmask_yinc, bdmask_zinc;
    Vector bdmask_xiv{ n_node_surface };
    Vector bdmask_yiv{ n_node_surface };
    Vector bdmask_ziv{ n_node_surface };
    
    size_t ib_xpre, ib_xinc, ib_ypre, ib_yinc, ib_zpre, ib_zinc;
    SurfaceElementRelation d_xpre, d_xinc, d_ypre, d_yinc, d_zpre, d_zinc;
    Vector b{ n_node_element };
    Vector btmp{ n_node_element };
    

    Vector bdtmp{ n_node_surface };
    	
    Vector bdtmp2{ n_node_surface };
            	
    bool xincf = Utilities::quadrant_x_positive_flag(quad);
    bool yincf = Utilities::quadrant_y_positive_flag(quad);
    bool zincf = Utilities::quadrant_z_positive_flag(quad);
    num_p boundary_p = boundary_->Value().p();
    num_p source_p = source_->Value().p();
    num_p solution_p = solution_->Value().p();
    num_p boundary_mask_p = mesh_.BoundaryMaskNodal().p();
    num_p boundary_mask_iv_p = mesh_.BoundaryMaskIvNodal().p();


    t1 = std::clock();
#ifdef MEATIME
    std::cout << "solve pre" << t1 - t0 << std::endl;
#endif
    time_t ta = 0, tb= 0, tc = 0, td = 0, te = 0, tf = 0;
    for (size_t ids = 0; ids < sweep_order.size(); ids++)
    {
      // Binding vectors		
      t0 = std::clock();
      size_t ide = sweep_order.at(ids);
      d_xpre = Utilities::upwind_source_direction_x(quad);
      d_xinc = Utilities::upwind_target_direction_x(quad);
      d_ypre = Utilities::upwind_source_direction_y(quad);
      d_yinc = Utilities::upwind_target_direction_y(quad);
      d_zpre = Utilities::upwind_source_direction_z(quad);
      d_zinc = Utilities::upwind_target_direction_z(quad);

      ib_xpre = esm[ide][Utilities::id(d_xpre)];
      ib_xinc = esm[ide][Utilities::id(d_xinc)];
      ib_ypre = esm[ide][Utilities::id(d_ypre)];
      ib_yinc = esm[ide][Utilities::id(d_yinc)];
      ib_zpre = esm[ide][Utilities::id(d_zpre)];
      ib_zinc = esm[ide][Utilities::id(d_zinc)];

      bd_xpre.Bind(n_node_surface, boundary_p + ib_xpre * n_node_surface);
      bd_xinc.Bind(n_node_surface, boundary_p + ib_xinc * n_node_surface);
      bd_ypre.Bind(n_node_surface, boundary_p + ib_ypre * n_node_surface);
      bd_yinc.Bind(n_node_surface, boundary_p + ib_yinc * n_node_surface);
      bd_zpre.Bind(n_node_surface, boundary_p + ib_zpre * n_node_surface);
      bd_zinc.Bind(n_node_surface, boundary_p + ib_zinc * n_node_surface);

      bdmask_xinc.Bind(n_node_surface, boundary_mask_p + ib_xinc * n_node_surface);
      bdmask_yinc.Bind(n_node_surface, boundary_mask_p + ib_yinc* n_node_surface);
      bdmask_zinc.Bind(n_node_surface, boundary_mask_p + ib_zinc* n_node_surface);
      bdmask_xiv.Bind(n_node_surface, boundary_mask_iv_p + ib_xinc * n_node_surface);
      bdmask_yiv.Bind(n_node_surface, boundary_mask_iv_p + ib_yinc* n_node_surface);
      bdmask_ziv.Bind(n_node_surface, boundary_mask_iv_p + ib_zinc* n_node_surface);      
      //bdmask_xiv.ScalarMulVectorAddScalar(-1.0, bdmask_xinc, 1.0);
      //bdmask_yiv.ScalarMulVectorAddScalar(-1.0, bdmask_yinc, 1.0);
      //bdmask_ziv.ScalarMulVectorAddScalar(-1.0, bdmask_zinc, 1.0);
      sor_e.Bind(n_node_element, source_p + ide * n_node_element);
      sol_e.Bind(n_node_element, solution_p + ide * n_node_element);


      t1 = std::clock();
      ta += t1 - t0;
      t0 = std::clock();

      b.AddScalarMulMatrixAddMulScalar(1.0, system_->mass_matrix(),
          Transport::N, sor_e, 0.0);
      b.AddScalarMulMatrixAddMulScalar(1.0, system_->lift_matrix(d_xpre),
          Transport::N, bd_xpre, 1.0);
      b.AddScalarMulMatrixAddMulScalar(1.0, system_->lift_matrix(d_ypre),
          Transport::N, bd_ypre, 1.0);
      b.AddScalarMulMatrixAddMulScalar(1.0, system_->lift_matrix(d_zpre),
          Transport::N, bd_zpre, 1.0);

      t1 = std::clock();
      tb += t1 - t0;
      t0 = std::clock();
      // Solve matrix
      sol_e.Copy(b);
      sol_e.Solve(system_->system_matrixlu());
      t1 = std::clock();
      tc += t1 - t0;
      t0 = std::clock();

      //s_e.gemx(CblasNoTrans, 1.0, s_->system_matrixiv(), b, 0.0);
      // Update boundary    		
      bdtmp.AddScalarMulMatrixAddMulScalar(1.0, system_->flux_matrix(d_xinc),
          Transport::N, sol_e, 0.0);
      bdtmp.PointMulVector(bdtmp, bdmask_xiv);

      bd_xinc.AddVector(1.0, bdtmp);

      bdtmp.AddScalarMulMatrixAddMulScalar(1.0, system_->flux_matrix(d_yinc),
          Transport::N, sol_e, 0.0);
      bdtmp.PointMulVector(bdtmp, bdmask_yiv);

      bd_yinc.AddVector(1.0, bdtmp);

      bdtmp.AddScalarMulMatrixAddMulScalar(1.0, system_->flux_matrix(d_zinc),
          Transport::N, sol_e, 0.0);
      bdtmp.PointMulVector(bdtmp, bdmask_ziv);

      bd_zinc.AddVector(1.0, bdtmp);
      t1 = std::clock();
      td += t1 - t0;
#ifdef DGNSDEBUG

      std::cout << "b" << std::endl << b;

      //std::cout << "bda xpre" << std::endl << bda_xpre;
      std::cout << "bd xpre" << std::endl << bd_xpre;

      //std::cout << "bda ypre" << std::endl << bda_ypre;
      std::cout << "bd ypre" << std::endl << bd_ypre;

      //std::cout << "bda zpre" << std::endl << bda_zpre;
      std::cout << "bd zpre" << std::endl << bd_zpre;

      //std::cout << "bda xinc" << std::endl << bda_xinc;
      std::cout << "bd xinc" << std::endl << bd_xinc;

      //std::cout << "bda yinc" << std::endl << bda_yinc;
      std::cout << "bd yinc" << std::endl << bd_yinc;

      //std::cout << "bda zinc" << std::endl << bda_zinc;
      std::cout << "bd zinc" << std::endl << bd_zinc;

      //std::cout << "sa" << std::endl << sa;
      std::cout << "sn" << std::endl << sol_e;

#endif
    }

#ifdef MEATIME
    std::cout << "bind " << ta << std::endl;
    std::cout << "calculate b " << tb << std::endl;
    std::cout << "solve " << tc << std::endl;
    std::cout << "update db " << td << std::endl;
#endif
  }
}



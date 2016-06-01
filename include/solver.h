#pragma once
#include <memory>
#include "global.h"
#include "system_matrices.h"
#include "mesh.h"
#include "source.h"
#include "boundary.h"
#include "solution.h"

namespace discontinues_galerkin_nodal_solver {	
  enum class SolverStatus {
    PreInit = 0,
    PostInit = 1,
    Closed = 2,
    OutputReady = 3
  };
  class SolverBase {
  public:
    SolverBase();
    SolverStatus status() const { return status_; }
    void status_set(SolverStatus ss);
    virtual void Procceed() = 0;
    virtual void Initialization() = 0;
    virtual void SetCondition() = 0;
    virtual void Solve() = 0;
  private:
    SolverStatus status_;
  };

  class SolverSingleAngle : public SolverBase{
  public:
    SolverSingleAngle(const MeshWithCoordinate& mesh, 
          std::string input_file_name = "input.txt");  
    SolverSingleAngle(const MeshWithCoordinate& mesh,
          const Vector& source_ext, const Vector& boundary_ext,
          std::string input_file_name = "input.txt");
    virtual void Procceed() override;
    virtual void Initialization() override;
    virtual void SetCondition() override;
    virtual void Solve() override;
    Vector solution() const { return solution_->Value(); }
    Vector xe() const { return mesh_.xe(); }
    Vector ye() const { return mesh_.ye(); }
    Vector ze() const { return mesh_.ze(); }
    Vector xs() const { return mesh_.xs(); }
    Vector ys() const { return mesh_.ys(); }
    Vector zs() const { return mesh_.zs(); }
    num_t mu() const { return system_->mu(); }
    num_t xi() const { return system_->xi(); }
    num_t eta() const { return system_->eta(); }
    num_t sigma() const { return system_->sigma(); }
  private:
    const MeshWithCoordinate& mesh_;
    ProblemDefinitionSingleAngle problem_definition_;
    std::shared_ptr<SystemMatrix> system_;
    std::shared_ptr<Source> source_;
    std::shared_ptr<Boundary> boundary_;
    std::shared_ptr<Solution> solution_; 
    Vector source_ext_, boundary_ext_;
    bool is_test_;
  };

	//class solver {
	//public:
	//	static const size_t total_node_element(const Mesh& m, size_t np);
	//	static const size_t total_node_surface(const Mesh& m, size_t np);

	//	const size_t total_node_element() const;
	//	const size_t total_node_surface() const;
	//public:
	//	solver(const system_matrix_angle& s, const Mesh& m);		
	//	~solver();
	//	//void source_set(vector_t s);
	//	//void boudnary_set(vector_t b);
	//	//vector_t source_get() const;
	//	//vector_t boudnary_get() const;

	//	void solve();

	//	void solve_test(size_t id);
	//	
	//	void solve_analytical(size_t id);

	//	num_p boundary_ptr() const { return boundary_p_; }
	//	num_p source_ptr() const { return source_p_; }
	//	num_p solution_ptr() const { return solution_p_; }

	//	vector_t x() const { return sg_.x(); }
	//	vector_t y() const { return sg_.y(); }
	//	vector_t z() const { return sg_.z(); }
	//	vector_t solution() const { return solution_; }
	//	vector_t x_b() const { return bg_.x(); }
	//	vector_t y_b() const { return bg_.y(); }
	//	vector_t z_b() const { return bg_.z(); }
	//	void solution_mean(vector_t data) const;

	//	size_t memory_total_solution() const {
	//		return m_->memory_element_total(s_->basis_element());
	//	}
	//private:
	//	source_generator sg_;
	//	boundary_generator bg_;
	//	num_p source_p_;
	//	num_p boundary_p_;
	//	num_p solution_p_;
	//	num_p boundary_mask_p_;

	//	vector_t source_;
	//	vector_t boundary_;
	//	vector_t solution_;
	//	vector_t boundary_mask_;
	//	const system_matrix_angle* s_;
	//	const Mesh* m_;

	//};

}
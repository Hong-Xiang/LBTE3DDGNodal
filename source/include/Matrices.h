#pragma once
#include <exception>
#include <vector>
#include "Global.h"

namespace DG3DNodal {

	class SystemMatrixStandard1D {
	public:
		static SystemMatrixStandard1D& getHandle();
		const vector_t getX(size_t Np) const;
		const matrix_t getM(size_t Np) const;
		const matrix_t getS(size_t Np) const;
		
	private:
		SystemMatrixStandard1D();
		void ReadMatFiles();
	private:
		static SystemMatrixStandard1D* handle_;
		std::vector<num_p> x;
		std::vector<num_p> s;
		std::vector<num_p> m;
		
	private:
		const size_t MINNP = 2;
		const size_t MAXNP = 20;
	};


	class SystemMatrixStandard2D {
	public:
		static SystemMatrixStandard2D& getHandle();
		const vector_t getX(size_t Np) const;
		const vector_t getY(size_t Np) const;
		const matrix_t getM(size_t Np) const;
		const matrix_t getSx(size_t Np) const;
		const matrix_t getSy(size_t Np) const;
		const size_t getTotalBasis(size_t Np) const { return Np*Np; }
	private:
		SystemMatrixStandard2D();
		void calculateMatrices();
	private:
		static SystemMatrixStandard2D* handle_;
		std::vector<num_p> x;
		std::vector<num_p> y;
		std::vector<num_p> sx;
		std::vector<num_p> sy;
		std::vector<num_p> m;
	private:
		const size_t MINNP = 2;
		const size_t MAXNP = 5;
	};

	class SystemMatrixStandard3D {
	public:
		static SystemMatrixStandard3D& getHandle();
		const vector_t getX(size_t Np) const;
		const vector_t getY(size_t Np) const;
		const vector_t getZ(size_t Np) const;
		const matrix_t getM(size_t Np) const;
		const matrix_t getSx(size_t Np) const;
		const matrix_t getSy(size_t Np) const;
		const matrix_t getSz(size_t Np) const;
		const matrix_t getLb(size_t Np) const;
		const matrix_t getLf(size_t Np) const;
		const matrix_t getLl(size_t Np) const;
		const matrix_t getLr(size_t Np) const;
		const matrix_t getLu(size_t Np) const;
		const matrix_t getLt(size_t Np) const;
		const matrix_t getBb(size_t Np) const;
		const matrix_t getBf(size_t Np) const;
		const matrix_t getBl(size_t Np) const;
		const matrix_t getBr(size_t Np) const;
		const matrix_t getBu(size_t Np) const;
		const matrix_t getBt(size_t Np) const;
		const size_t getTotalBasis(size_t Np) const { return Np*Np*Np; }
		const size_t getTotalBasis2d(size_t Np) const { return Np*Np; }
	private:
		SystemMatrixStandard3D();	
		static SystemMatrixStandard3D* handle_;
		void calculateMatrices();

	private:
		std::vector<num_p> x;
		std::vector<num_p> y;
		std::vector<num_p> z;
		std::vector<num_p> m;
		std::vector<num_p> sx;
		std::vector<num_p> sy;
		std::vector<num_p> sz;

		//Lift Flux (Np^2, Np^3)
		std::vector<num_p> lb;
		std::vector<num_p> lf;
		std::vector<num_p> ll;
		std::vector<num_p> lr;
		std::vector<num_p> lu;
		std::vector<num_p> lt;

		//In flux (Np^3, Np^2)
		std::vector<num_p> bb;
		std::vector<num_p> bf;
		std::vector<num_p> bl;
		std::vector<num_p> br;
		std::vector<num_p> bu;
		std::vector<num_p> bt;

	private:
		const size_t MINNP = 2;
		const size_t MAXNP = 4;
	};


	class SystemMatrixStandardOutOfOrder : public std::exception {};
}
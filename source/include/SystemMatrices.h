#pragma once
#include "Global.h"

namespace DG3DNodal {
	class SystemMatrix3D {
	public:
		SystemMatrix3D(num_t h, num_t mu, num_t xi, num_t eta, num_t sigma, size_t Np);
		~SystemMatrix3D();
		const SystemMatrix3D& getHandle() const;
	public:
		const vector_t getX() const;
		const vector_t getY() const;
		const vector_t getZ() const;
		const matrix_t getSysM() const;
		const matrix_t getMx() const;
		const matrix_t getMy() const;
		const matrix_t getMz() const;		
		const matrix_t getFx() const;
		const matrix_t getFy() const;
		const matrix_t getFz() const;
		size_t getRefineLevel() const;
		size_t getTotalBasis() const { return np_*np_*np_; }

		const matrix_t getSysMLU() const { return sysml; }
		const matrix_t getSysMI() const { return sysmi; }
		const lapack_int* getipiv() const { return ipiv; }		
	private:
		void calculate();
	private:
		num_t h_, mu_, xi_, eta_, sigma_;
		size_t np_;
	private:
		num_p xp, yp, zp;
		num_p sysmp, mxp, myp, mzp, fxp, fyp, fzp;
		vector_t x, y, z;
		matrix_t sysm, mx, my, mz, fx, fy, fz;
		size_t refine_level;
		num_p sysmlp, sysmup, sysmip;
		matrix_t sysml, sysmu, sysmi;
		lapack_int* ipiv;
	};
}
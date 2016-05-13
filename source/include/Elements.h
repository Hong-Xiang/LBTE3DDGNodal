#pragma once
#include "Global.h"
#include "Matrices.h"
#include "Interfaces.h"

namespace DG3DNodal {
	
	class Element3D {
	public:
		Element3D(size_t eId, mesh_refine_lvl_t refine_lvl, size_t ix, size_t iy, size_t iz, num_t xc, num_t yc, num_t zc) {
			m_eId = eId;
			m_ix = ix;
			m_iy = iy;
			m_iz = iz;
			m_xc = xc;
			m_yc = yc;
			m_zc = zc;
			m_refine_lvl = refine_lvl;
			interfaces.resize(6);
			for (size_t i = 0; i < 6; i++)
			{
				interfaces.at(i) = nullptr;
			}
		}
	public:
		num_t getXc() const { return m_xc; }
		num_t getYc() const { return m_yc; }
		num_t getZc() const { return m_zc; }
		size_t getIX() const { return m_ix; }
		size_t getIY() const { return m_iy; }
		size_t getIZ() const { return m_iz; }
		size_t getID() const { return m_eId; }
		void addChildElements(Quadrant quad) {//TODO;
		}
		void addInterface(InterfaceDirection direction, Interface* pI) {
			interfaces.at(direction) = pI;
		}
		Interface* get_interface(InterfaceDirection direction) const {
			return interfaces.at(direction);
		}
		const size_t get_refine_level() const { return m_refine_lvl; }
		const Element3D* handle() {
			return this;
		}
	public:		
		size_t m_eId;
		mesh_refine_lvl_t m_refine_lvl;
		num_t m_xc, m_yc, m_zc;
		size_t m_ix, m_iy, m_iz;
		std::vector<Interface*> interfaces;
		std::vector<Element3D*> child_elements;
	};
	
}
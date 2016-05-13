#pragma once
#include "Global.h"

namespace DG3DNodal {

	class Interface {
	public:
		Interface(size_t iid, mesh_refine_lvl_t refine_lvl) {
			m_start = 0;
			m_iid = iid;
			m_refine_lvl = refine_lvl;
			m_isBoundary = false;
		}
		Interface* getHandle() {
			return this;
		}

		size_t id() const {
			return m_iid; 
		}
		void markAsBoundary() {
			m_isBoundary = true;
		}
		bool isBoundary() const {
			return m_isBoundary;
		}
		size_t getStartPosition() const {
			return m_start;
		}
		void setStartPosition(size_t pos) {
			m_start = pos;
		}

	private:		
		size_t m_iid;
		mesh_refine_lvl_t m_refine_lvl;
		bool m_isBoundary;
		size_t m_start;
	};
}
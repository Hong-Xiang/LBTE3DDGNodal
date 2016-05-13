#pragma once
#include <mkl.h>
#include <iostream>
//Some definitions of helper types to enhance use of MKL
namespace MKLEXT
{

	template<typename Float_t>
	class Matrix {
	public:
		Matrix(size_t n_row, size_t n_column, Float_t* data, size_t inc_row = 1, size_t inc_column = 1)
			: data_(data), n_row_e(n_row), n_column_e(n_column), n_row_r(n_row*inc_row), n_column_r(n_column*inc_column), default_inc_row(inc_row), default_inc_column(inc_column)
		{
			isNull = false;
		}

		Matrix() : data_(nullptr){
			isNull = true;
		}

		void bind(size_t n_row, size_t n_column, Float_t* data, size_t inc_row = 1, size_t inc_column = 1) {
			data_ = data;
			n_row_e = n_row;
			n_column_e = n_column;
			n_row_r = n_row*inc_row;
			n_column_r = n_column*inc_column;
			default_inc_row = inc_row;
			default_inc_column = inc_column;
			isNull = false;
		}

		
		Float_t* getPtr() { return data_; }
		const Float_t* getPtr() const { return data_; }
		bool null_check() const { return isNull; }

		//size_t getRowIncrement() const { return default_inc_row; }
		//size_t getColumnIncrement() const { return default_inc_column; }
		//void setIncrement(size_t row_inc, size_t column_inc) {
		//	default_inc_row = row_inc;
		//	default_inc_column = column_inc;
		//}
		Float_t& operator()(size_t i, size_t j) { return data_[i*default_inc_row + j*default_inc_column*n_row_r]; }
		Float_t& operator()(size_t i, size_t j) const { return data_[i*default_inc_row + j*default_inc_column*n_row_r]; }
		void print(std::ostream& os = std::cout) const {
			for (size_t i = 0; i < n_row_e; i++)
			{
				for (size_t j = 0; j < n_column_e; j++)
				{
					os << this->operator()(i,j) << "\t";
				}
				os << std::endl;
			}
		}
	private:
		Float_t* data_;
		size_t n_row_e, n_column_e;	//number of effective row and column
		size_t n_row_r, n_column_r;	//number of real row and column
		size_t default_inc_row, default_inc_column;
		bool isNull;
	};

	template<typename Float_t>
	class Vector {
	public:
		Vector(size_t n_size, Float_t* data, size_t inc = 1)
			: data_(data), n_size_e(n_size), n_size_r(n_size*inc), default_inc(inc)
		{
			isNull = false;
		}

		Vector() 
			: data_(nullptr)
		{
			isNull = true;
		}
		
		void bind(size_t n_size, Float_t* data, size_t inc = 1)
		{
			data_ = data;
			n_size_e = n_size;
			n_size_r = n_size*inc;
			default_inc = inc;
			isNull = false;
		}

		Float_t* getPtr() { return data_; }
		const Float_t* getPtr() const { return data_; }

		bool null_check() const { return isNull; }

		Float_t& operator()(size_t i) { return data_[i*default_inc]; }
		Float_t& operator()(size_t i) const { return data_[i*default_inc]; }
		void print(std::ostream& os = std::cout) const {
			for (size_t i = 0; i < n_size_e; i++)
			{
				os << this->operator()(i) << "\t";
			}
			os << std::endl;
		}
	private:
		Float_t* data_;
		size_t n_size_e;
		size_t n_size_r;
		size_t default_inc;
		bool isNull;
	};
}
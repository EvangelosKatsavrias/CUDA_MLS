//   CUDA_MLS Framework
//
//   Copyright 2017-2018 Evangelos D. Katsavrias, Luxembourg
//
//   This file is part of the CUDA_MLS Framework.
//
//   CUDA_MLS Framework is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License version 3 as published by
//   the Free Software Foundation.
//
//   CUDA_MLS Framework is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CUDA_MLS Framework.  If not, see <https://www.gnu.org/licenses/>.
//
//   Contact Info:
//   Evangelos D. Katsavrias
//   email/skype: vageng@gmail.com
// -----------------------------------------------------------------------

#ifndef TRANSFORMATIONMATRIXHEADER
#define TRANSFORMATIONMATRIXHEADER


#include<stdexcept>
#include<exception>
#include<iostream>
#include<string>
#include<fstream>
#include<sstream>


template<class T>
class TransformationMatrix
{

	protected:
		T*	m_matrix = 0;
		size_t 	m_numOfComponents;
		size_t 	m_numOfDomainComponents;

	public:
		TransformationMatrix( size_t numOfDomainComponents = 3, size_t numOfComponents = 0 ) 
		{ 
			set_numOfDomainComponents(numOfDomainComponents); 
			if (numOfComponents) set_numOfComponents(numOfComponents); 
			set_identity( ); 
		}
		TransformationMatrix( T* matrix, size_t numOfDomainComponents = 3, size_t numOfComponents = 0 ) { set_matrix(matrix, numOfDomainComponents, numOfComponents ); }
		TransformationMatrix( std::istream& streamIn, size_t numOfDomainComponents = 3, size_t numOfComponents = 0 ) { set_matrix(streamIn, numOfDomainComponents, numOfComponents); }
		TransformationMatrix( TransformationMatrix<T>& matrix ) { set_matrix( matrix ); }
		TransformationMatrix( const TransformationMatrix<T>& matrix ) 
		{ 
			size_t domComp = matrix.get_numOfDomainComponents(); 
			size_t comp = matrix.get_numOfComponents();
			T* mat = matrix.get_matrix();
			set_matrix( mat, domComp, comp );
		}
		~TransformationMatrix() { delete[] m_matrix; }

		
		int set_component ( size_t element_row, size_t element_column, T value ) { m_matrix[ element_row*(m_numOfDomainComponents+1) +element_column ] = value; return 1; }
		int set_numOfDomainComponents ( size_t numOfComponents ) { m_numOfDomainComponents = numOfComponents; set_numOfComponents( (numOfComponents+1)*(numOfComponents+1) ); }
		int set_numOfComponents ( size_t numOfComponents ) { m_numOfComponents = numOfComponents; if (m_matrix) delete[] m_matrix; m_matrix = new T[m_numOfComponents]; }
		int set_identity( ) 
		{
			std::fill( m_matrix, m_matrix+m_numOfComponents, T(0) );
			for (size_t i = 0; i < m_numOfDomainComponents+1; i++ ) set_component(i, i, T(1));
		}
		int set_matrix( T* a )
		{ 
			for (size_t i = 0; i < m_numOfComponents; i++) m_matrix[i] = a[i]; 

		}
		int set_matrix( T* matrix, size_t numOfDomainComponents, size_t numOfComponents = 0 ) 
		{
			set_numOfDomainComponents(numOfDomainComponents);
			if ( numOfComponents ) set_numOfComponents( numOfComponents );
			set_matrix(matrix);
		}
		int set_matrix( std::istream& streamIn, size_t numOfDomainComponents, size_t numOfComponents = 0 ) 
		{ 
			set_numOfDomainComponents(numOfDomainComponents);
			if ( numOfComponents ) set_numOfComponents( numOfComponents );
			T temp; for ( size_t comp = 0; comp < numOfComponents; comp++ ) {streamIn >> temp; m_matrix[comp] = temp;} 
		}
		int set_matrix( TransformationMatrix<T>& matrix ) 
		{ 
			set_matrix( matrix.get_matrix(), matrix.get_numOfDomainComponents() ); 
		}


		typedef	T* iterator;
		typedef const T* const_iterator;
		T*	begin() { return m_matrix;}
		T*	end() 	{ return m_matrix+m_numOfComponents;}
		T	front() { return m_matrix[0]; }
		T	back() 	{ return *(end()-1); }
		size_t	size() 	{ return m_numOfComponents; }
		T	at(size_t index) 	{ return m_matrix[index]; }

		T* 	data 		( )						{ return m_matrix; }
		T& 	get_component 	( size_t element_row, size_t element_column ) 	{ return m_matrix[ element_row*(m_numOfDomainComponents+1) +element_column ]; }
		T* 	get_matrix 	( )			const			{ return m_matrix; }
		size_t 	get_numOfComponents( ) 			const			{ return m_numOfComponents; }
		size_t 	get_numOfDomainComponents( ) 		const			{ return m_numOfDomainComponents; }
		T& 	operator() 	( size_t element_row, size_t element_column ) 	{ return get_component( element_row, element_column ); }
		T& 	operator[] 	( size_t element ) 				{ return m_matrix[element]; }
		TransformationMatrix&	operator= ( const TransformationMatrix& matrix )	
		{ 
			size_t domComp = matrix.get_numOfDomainComponents(); 
			size_t comp = matrix.get_numOfComponents();
			T* mat = matrix.get_matrix();
			set_matrix( mat, domComp, comp );
			return *this;
		}
		TransformationMatrix&	operator* ( TransformationMatrix& follow) 
		{
			TransformationMatrix res = *this;
			TransformationMatrix& curr = *this;

			for ( size_t res_row = 0; res_row < m_numOfDomainComponents+1; res_row++ ) {
				for ( size_t res_col = 0; res_col < m_numOfDomainComponents+1; res_col++ ) {
					res(res_row, res_col) = 0;
					for ( size_t comp = 0; comp < m_numOfDomainComponents+1; comp++ ) {
					//	res(res_row, res_col) += follow( res_row, comp )*curr( comp, res_col );
						res(res_row, res_col) += curr(res_row, comp)*follow(comp, res_col );
					}
				}
			}

			set_matrix ( res );
			return *this;
		}
		int clear() { set_identity( ); }



		// =========================== special transformations
		int set_translation ( T* dx ) {
			set_identity();
			for ( size_t i = 0; i < m_numOfDomainComponents; i++) set_component(i, m_numOfDomainComponents, dx[i] );
		}

		TransformationMatrix& follow_translation(T* dx) {
			TransformationMatrix newTr(m_numOfDomainComponents); newTr.set_translation(dx);
			*this = newTr*(*this);
			return *this;
		}

		int set_rotation_x ( T theta_x ) {
			set_identity(); T c(cos(theta_x)), s(sin(theta_x));
			set_component(1, 1, c);
			set_component(2, 2, c);
			set_component(1, 2, -s);
			set_component(2, 1, s);
		}

		TransformationMatrix& follow_rotation_x(T theta_x) {
			TransformationMatrix newTr(m_numOfDomainComponents); newTr.set_rotation_x(theta_x);
			*this = newTr*(*this);
			return *this;
		}

		int set_rotation_y ( T theta_y ) {
			set_identity(); T c(cos(theta_y)), s(sin(theta_y));
			set_component(0, 0, c);
			set_component(2, 2, c);
			set_component(0, 2, s);
			set_component(2, 0, -s);
		}

		TransformationMatrix& follow_rotation_y(T theta_y) {
			TransformationMatrix newTr(m_numOfDomainComponents); newTr.set_rotation_y(theta_y);
			*this = newTr*(*this);
			return *this;
		}

		int set_rotation_z ( T theta_z ) {
			set_identity(); T c(cos(theta_z)), s(sin(theta_z));
			set_component(0, 0, c);
			set_component(1, 1, c);
			set_component(0, 1, -s);
			set_component(1, 0, s);
		}

		TransformationMatrix& follow_rotation_z(T theta_z) {
			TransformationMatrix newTr(m_numOfDomainComponents); 
			newTr.set_rotation_z(theta_z);
			*this = newTr*(*this);
			return *this;
		}

		int set_rotation ( T theta_z, T theta_y = 0, T theta_x = 0) {
			set_identity();
			if ( m_numOfDomainComponents == 2 ) follow_rotation_z( theta_z );
			else follow_rotation(theta_z, theta_y, theta_x);
		}

		TransformationMatrix& follow_rotation ( T theta_z, T theta_y = 0, T theta_x = 0 ) {
			if ( m_numOfDomainComponents > 2 ) { follow_rotation_x(theta_x); follow_rotation_y(theta_y); }	
			follow_rotation_z(theta_z);
			return *this;
		}


		int set_rotation (T u_x, T u_y, T u_z, T theta) {
			TransformationMatrix newTr(m_numOfDomainComponents); T c(cos(theta)), s(sin(theta));
			newTr(0, 0) = c + u_x*u_x*(1-c);
			newTr(0, 1) = u_x*u_y*(1-c) -u_z*s;
			newTr(0, 2) = u_x*u_z*(1-c) +u_y*s;
			newTr(1, 0) = u_x*u_y*(1-c) +u_z*s;
			newTr(1, 1) = c + u_y*u_y*(1-c);
			newTr(1, 2) = u_y*u_z*(1-c) -u_x*s;
			newTr(2, 0) = u_x*u_z*(1-c) -u_y*s;
			newTr(2, 1) = u_z*u_y*(1-c) +u_x*s;
			newTr(2, 2) = c + u_z*u_z*(1-c);
			*this = newTr;
		}


		TransformationMatrix& follow_rotation (T u_x, T u_y, T u_z, T theta)
		{
			TransformationMatrix newTr(m_numOfDomainComponents); newTr.set_rotation(u_x, u_y, u_z, theta);
			*this = newTr*(*this);
			return *this;
		}


		int set_rotation_quaternion (T w, T x, T y, T z) {
			T len = sqrt( x*x +y*y +z*z), u_x = x/len, u_y = y/len, u_z = z/len, theta = 2*atan2(len, w);
			set_rotation(u_x, u_y, u_z, theta);
		}

		TransformationMatrix& follow_rotation_quaternion (T w, T x, T y, T z)
		{
			TransformationMatrix newTr(m_numOfDomainComponents); newTr.set_rotation_quaternion(w, x, y, z);
			*this = newTr*(*this);
			return *this;
		}

		int set_rotation ( T x_0, T y_0, T z_0, T theta_z, T theta_y, T theta_x) {
			T x[3]{-x_0, -y_0, -z_0};
			set_translation( x );
			follow_rotation(theta_z, theta_y, theta_x);
			x[0] = x_0; x[1] = y_0; x[2] = z_0;
			follow_translation( x );
		}

		TransformationMatrix& follow_rotation ( T x_0, T y_0, T z_0, T theta_z, T theta_y, T theta_x ) {
			T x[3]{x_0, y_0, z_0};
			follow_translation( x );
			follow_rotation(theta_z, theta_y, theta_x);
			x[0] = -x_0; x[1] = -y_0; x[2] = -z_0;
			follow_translation( x );
			return *this;
		}
};


template<class T>
std::ostream& operator << ( std::ostream& out, TransformationMatrix<T>& matrix)
{
	for ( auto& element: matrix ) out << element << std::endl;
	return out;
}


template<class T>
std::istream& operator >> ( std::istream& in, TransformationMatrix<T>& matrix)
{
	for ( auto& element: matrix ) in >> element;
	return in;
}


/*
template<class T>
class TransformationMatrix2D: public TransformationMatrix<T>
{

	public:
		TransformationMatrix2D( T a11 = 1, T a12 = 0, T a13 = 0, 
					T a21 = 0, T a22 = 1, T a23 = 0, 
					T a31 = 0, T a32 = 0, T a33 = 1 )
		{ set_matrix( a11, a12, a13, a21, a22,  a23, a31,  a32,  a33 ); set_numOfComponents(9); }
		TransformationMatrix2D( T* matrix, size_t numOfComponents = 9 ) { set_matrix(matrix, numOfComponents); set_numOfComponents(9); }
		TransformationMatrix2D( std::istream& streamIn, size_t numOfComponents = 9 ) { set_matrix(streamIn, numOfComponents); set_numOfComponents(9); }
//		TransformationMatrix2D( TransformationMatrix<T>& matrix ) { set_matrix( matrix ); }
		~TransformationMatrix2D() {}


		int set_component ( size_t element_row, size_t element_column, T value ) { return this->m_matrix[ element_row*3 +element_column ] = value; }
		int set_numOfComponents ( size_t numOfComponents ) { this->m_numOfComponents = numOfComponents; }
		int set_identity( ) { set_matrix(); }
		int set_matrix( T a11 = 1, T a12 = 0, T a13 = 0, 
				T a21 = 0, T a22 = 1, T a23 = 0, 
				T a31 = 0, T a32 = 0, T a33 = 1)
		{
			this->m_matrix[0] = a11;  this->m_matrix[1] = a12;  this->m_matrix[2] = a13;
			this->m_matrix[3] = a21;  this->m_matrix[4] = a22;  this->m_matrix[5] = a23;
			this->m_matrix[6] = a31;  this->m_matrix[7] = a32;  this->m_matrix[8] = a33;
		}
		int set_matrix( T* matrix, size_t numOfComponents = 9 ) { TransformationMatrix<T>::set_matrix(matrix, numOfComponents); }
		int set_matrix( TransformationMatrix<T>& matrix ) { set_matrix( matrix.get_matrix(), matrix.get_numOfComponents() ); }
		int set_matrix( std::istream& streamIn, size_t numOfComponents = 9 ) { TransformationMatrix<T>::set_matrix(streamIn, numOfComponents); }


		typedef	T* iterator;
		typedef const T* const_iterator;
		T*	begin() { return m_matrix;}
		T*	end() 	{ return m_matrix+m_numOfComponents;}
		T	front() { return m_matrix[0]; }
		T	back() 	{ return *(end()-1); }
		size_t	size() 	{ return m_numOfComponents; }
		T	at(size_t index) 	{ return m_matrix[index]; }

		T* 	data 		( )						{ return m_matrix; }
		T& 	get_component 	( size_t element_row, size_t element_column ) 	{ return m_matrix[ element_row*3 +element_column ]; }
		T* 	get_matrix 	( )						{ return m_matrix; }
		size_t 	get_numOfComponents( ) 						{ return m_numOfComponents; }
		T& 	operator() 	( size_t element_row, size_t element_column ) 	{ return get_component( element_row, element_column ); }
		T& 	operator[] 	( size_t element ) 				{ return m_matrix[element]; }
//		TransformationMatrix&	operator= ( const TransformationMatrix& other )	{ set_matrix( other ); return *this; }
		TransformationMatrix&	operator* ( TransformationMatrix& follow) 
		{
			TransformationMatrix res(*this);
			TransformationMatrix& curr = *this;

			for ( size_t res_row = 0; res_row < 3; res_row++ ) {
				res(res_row, 0) = follow(res_row, 0)*curr(0, 0) +follow(res_row, 1)*curr(1, 0) +follow(res_row, 2)*curr(2, 0);
				res(res_row, 1) = follow(res_row, 0)*curr(0, 1) +follow(res_row, 1)*curr(1, 1) +follow(res_row, 2)*curr(2, 1);
				res(res_row, 2) = follow(res_row, 0)*curr(0, 2) +follow(res_row, 1)*curr(1, 2) +follow(res_row, 2)*curr(2, 2);
			}

			set_matrix ( res );
			return *this;
		}
		int clear() { set_identity( ); }



		// =========================== special transformations
		int set_translation ( T dx, T dy = 0 ) {
			set_identity();
			size_t colIndex(2);
			set_component(0, colIndex, dx);
			set_component(1, colIndex, dy);
		}

		TransformationMatrix& follow_translation(T dx, T dy = 0) {
			TransformationMatrix newTr; newTr.set_translation(dx, dy);
			*this = newTr*(*this);
			return *this;
		}


		int set_rotation_z ( T theta_z ) {
			set_identity(); T c(cos(theta_z)), s(sin(theta_z));
			set_component(0, 0, c);
			set_component(1, 1, c);
			set_component(0, 1, -s);
			set_component(1, 0, s);
		}


		TransformationMatrix& follow_rotation_z(T theta_z) {
			TransformationMatrix newTr; newTr.set_rotation_z(theta_z);
			*this = newTr*(*this);
			return *this;
		}


		int set_rotation ( T theta_z) {
			set_identity(); follow_rotation(theta_z);
		}

		TransformationMatrix& follow_rotation ( T theta_z) {
			follow_rotation_z(theta_z);
			return *this;
		}


		int set_rotation ( T x_0, T y_0, T theta_z) {
			set_translation( x_0, y_0 );
			follow_rotation(theta_z);
			follow_translation( -x_0, -y_0 );
		}


		TransformationMatrix& follow_rotation ( T x_0, T y_0, T theta_z ) {
			follow_translation( x_0, y_0 );
			follow_rotation(theta_z);
			follow_translation( -x_0, -y_0 );
			return *this;
		}


};
*/


#endif

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

#include"LinearSolver.h"
#include"gaussElimination.h"
#include"conjugateGradient.h"
//#include""
#include<iostream>


template<class T>
LinearSolverBase<T>::LinearSolverBase( int numOfColumns)
{
	set_numberOfColumns(numOfColumns);
	set_solverType( gaussEliminationWithBackSubst );
}

template<class T>
void LinearSolverBase<T>::set_numberOfColumns(int numOfCols)
{m_numCols = numOfCols;}


template<class T>
void LinearSolverBase<T>::set_solverFunction( void(*solverFunction)(T* matrix, T* rhs, T* sol, int numOfCols) )
{ m_solverFunction = solverFunction; m_linearSolverType = NonPredefinedSolver;}


template<class T>
void LUSolver<T>::set_solverFunction( void(*solverFunction)(T* lowerMatrix, T* upperMatrix, T* rhs, T* sol, int numOfCols) )
{ m_solverFunction = solverFunction; this->m_linearSolverType = NonPredefinedSolver; }


template<class T>
void LDUSolver<T>::set_solverFunction( void(*solverFunction)(T* lowerMatrix, T* upperMatrix, T* D, T* rhs, T* sol, int numOfCols) )
{ m_solverFunction = solverFunction; this->m_linearSolverType = NonPredefinedSolver; }


template<class T>
void LDLSolver<T>::set_solverFunction( void(*solverFunction)(T* lowerMatrix, T* D, T* rhs, T* sol, int numOfCols) )
{ m_solverFunction = solverFunction; this->m_linearSolverType = NonPredefinedSolver; }


template<class T>
void CholeskyBandedSolver<T>::set_solverFunction( void(*solverFunction)(T* lowerMatrix, T* rhs, T* sol, int numOfCols, int numOfBands) )
{ m_solverFunction = solverFunction; this->m_linearSolverType = NonPredefinedSolver; }


template<class T>
void CholeskySkylineSolver<T>::set_solverFunction( void(*solverFunction)(T* lowerMatrix, T* rhs, T* sol, int numOfCols, int* columnPointers) )
{ m_solverFunction = solverFunction; this->m_linearSolverType = NonPredefinedSolver; }


template<class T>
void LUSolver<T>::set_decompositionFunction( void(*decompositionFunction)(T* matrix, T* lowerMatrix, T* upperMatrix, int numOfCols) )
{ m_decompositionFunction = decompositionFunction; }


template<class T>
void LDUSolver<T>::set_decompositionFunction( void(*decompositionFunction)(T* matrix, T* lowerMatrix, T* upperMatrix, T* D, int numOfCols) )
{ m_decompositionFunction = decompositionFunction; }


template<class T>
void LDLSolver<T>::set_decompositionFunction( void(*decompositionFunction)(T* matrix, T* lowerMatrix, T* D, int numOfCols) )
{ m_decompositionFunction = decompositionFunction; }


template<class T>
void CholeskySolver<T>::set_decompositionFunction( void(*decompositionFunction)(T* matrix, T* lowerMatrix, int numOfCols) )
{ m_decompositionFunction = decompositionFunction; }


template<class T>
void CholeskyBandedSolver<T>::set_decompositionFunction( void(*decompositionFunction)(T* matrix, T* lowerMatrix, int numOfCols, int numOfBands) )
{ m_decompositionFunction = decompositionFunction; }


template<class T>
void CholeskySkylineSolver<T>::set_decompositionFunction( void(*decompositionFunction)(T* matrix, int numOfCols, int* columnPointers) )
{ m_decompositionFunction = decompositionFunction; }


template<class T>
void LinearSolverBase<T>::set_solverType(LinearSolverType type)
{
	switch (type)
	{
		case gaussEliminationWithBackSubst:  	{set_solverFunction( gaussElimination_backSubst ); break;}
		case ConjugateGradient:  		{set_solverFunction( CG ); break;}
		case LU_Doolittle:			{set_solverFunction( LU_Doolittle_solver ); set_decompositionFunction(LU_Doolittle_decomposer); break;}
		case LU_Crout:				{set_solverFunction( LU_Crout_solver ); set_decompositionFunction(LU_Crout_decomposer ); break;}
		case LDU:				{set_solverFunction( LDU_solver ); set_decompositionFunction(LDU_decomposer); break; }
		case LDL:				{set_solverFunction( LDL_solver ); set_decompositionFunction(LDL_decomposer); break; }
		case Cholesky: 				{set_solverFunction( Cholesky_solver ); set_decompositionFunction( Cholesky_decomposer ); break; }
		case CholeskyBanded: 			{set_solverFunction( Cholesky_solver_banded ); set_decompositionFunction(Cholesky_decomposer_banded ); break; }
		case CholeskySkyline: 			{set_solverFunction( Cholesky_solver_skyline ); set_decompositionFunction(Cholesky_decomposer_skylineInPlace ); break; }

	};
	m_linearSolverType = type;
}


template<class T>
void LinearSolverBase<T>::operator () (T* matrix, T* rhs, T* sol)
{ if (!sol) sol = rhs;
m_solverFunction(matrix, rhs, sol, m_numCols); }


template<class T>
void LUSolver<T>::operator () (T* matrix, T* rhs, T* sol)
{
if (!sol) sol = rhs; 

if ( this->get_givenMatrix() != matrix )
{
	this->set_givenMatrix( matrix);
	m_internalLowerMatrix.resize( this->m_numCols*this->m_numCols );
	m_internalUpperMatrix.resize( this->m_numCols*this->m_numCols );
	m_decompositionFunction( matrix, m_internalLowerMatrix.data(), m_internalUpperMatrix.data(), this->m_numCols );
}

	m_solverFunction( m_internalLowerMatrix.data(), m_internalUpperMatrix.data(), rhs, sol, this->m_numCols);

}


template<class T>
void LUSolver<T>::operator () (T* lowerMatrix, T* upperMatrix, T* rhs, T* sol)
{

if (!sol) sol = rhs; 
this->m_solverFunction( lowerMatrix, upperMatrix, rhs, sol, this->m_numCols);

}


template<class T>
void LDUSolver<T>::operator () (T* matrix, T* rhs, T* sol)
{

if (!sol) sol = rhs; 

if ( this->get_givenMatrix() != matrix )
{
	this->set_givenMatrix ( matrix);
	m_internalLowerMatrix.resize( this->m_numCols*this->m_numCols );
	m_internalUpperMatrix.resize( this->m_numCols*this->m_numCols );
	m_internalDiagonal.resize( this->m_numCols );
	m_decompositionFunction( matrix, m_internalLowerMatrix.data(), m_internalUpperMatrix.data(), m_internalDiagonal.data(), this->m_numCols );
}

	m_solverFunction( m_internalLowerMatrix.data(), m_internalUpperMatrix.data(), m_internalDiagonal.data(), rhs, sol, this->m_numCols);

}


template<class T>
void LDUSolver<T>::operator () (T* lowerMatrix, T* upperMatrix, T* D, T* rhs, T* sol)
{

if (!sol) sol = rhs; 
this->m_solverFunction( lowerMatrix, upperMatrix, D, rhs, sol, this->m_numCols);

}


template<class T>
void LDLSolver<T>::operator () (T* matrix, T* rhs, T* sol)
{

if (!sol) sol = rhs; 

if ( this->get_givenMatrix() != matrix )
{
	this->set_givenMatrix (matrix);
	m_internalLowerMatrix.resize( this->m_numCols*this->m_numCols );
	m_internalDiagonal.resize( this->m_numCols );
	m_decompositionFunction( matrix, m_internalLowerMatrix.data(), m_internalDiagonal.data(), this->m_numCols );
}

	m_solverFunction( m_internalLowerMatrix.data(), m_internalDiagonal.data(), rhs, sol, this->m_numCols);

}


template<class T>
void LDLSolver<T>::operator () (T* lowerMatrix, T* D, T* rhs, T* sol)
{

if (!sol) sol = rhs; 
this->m_solverFunction( lowerMatrix, D, rhs, sol, this->m_numCols);

}



template<class T>
void CholeskySolver<T>::operator () (T* matrix, T* rhs, T* sol)
{

if (!sol) sol = rhs; 

if ( this->get_givenMatrix () != matrix )
{
	this->set_givenMatrix (matrix);
	if ( !m_decomposedInputFlag )
	{
		if ( m_decomposeInPlaceFlag )	m_decompositionFunction( matrix, matrix, this->m_numCols );
		else {
			m_decomposedMatrix.resize( this->m_numCols*this->m_numCols );			
			m_decompositionFunction( matrix, m_decomposedMatrix.data(), this->m_numCols ); 
		}
	}
}

if ( m_decomposeInPlaceFlag ) this->m_solverFunction(this->get_givenMatrix(), rhs, sol, this->m_numCols);
else this->m_solverFunction( m_decomposedMatrix.data(), rhs, sol, this->m_numCols);

}


template<class T>
void CholeskyBandedSolver<T>::operator () (T* matrix, T* rhs, T* sol, int numOfColumns, int numOfBands)
{

if ( !sol ) sol = rhs; if ( !numOfBands ) set_numOfBands(numOfBands);

if ( this->get_givenMatrix() != matrix || numOfColumns || numOfBands )
{
	this->set_givenMatrix ( matrix);

	if ( m_bandedInputFlag ) {

		if ( numOfColumns ) set_numberOfColumns( numOfColumns);
		if ( numOfBands ) set_numOfBands(numOfBands);
		m_bandedMatrix->set_bandedMatrix(matrix);

		if ( !this->m_decomposedInputFlag ) {
			if ( this->m_decomposeInPlaceFlag ) m_decompositionFunction( this->get_givenMatrix(), this->get_givenMatrix(), get_numberOfColumns(), get_numOfBands() );
			else {
				m_bandedMatrix->set_content( );
				m_decompositionFunction( this->get_givenMatrix(), m_bandedMatrix->get_bandedMatrix(), get_numberOfColumns(), get_numOfBands() );
			}
		}


	}
	else {

		m_bandedMatrix->set_fromFullMatrix(matrix, numOfColumns, numOfBands);
		m_decompositionFunction( m_bandedMatrix->get_bandedMatrix(), m_bandedMatrix->get_bandedMatrix(), m_bandedMatrix->get_numOfColumns(), m_bandedMatrix->get_numOfBands() );

	}
}


this->m_solverFunction( m_bandedMatrix->get_bandedMatrix() , rhs, sol, get_numberOfColumns(), get_numOfBands() );

}


template<class T>
void CholeskySkylineSolver<T>::operator () (T* matrix, T* rhs, T* sol, int numOfColumns, int* columnPointers )
{

if ( !sol ) sol = rhs; 
if ( columnPointers ) set_columnPointers(columnPointers);
if ( numOfColumns ) set_numberOfColumns(numOfColumns);


if ( this->get_givenMatrix() != matrix )
{
	this->set_givenMatrix ( matrix); 

	if ( m_skylineInputFlag ) {
		m_skylineMatrix->set_skylineMatrix(matrix);
		if ( !this->m_decomposedInputFlag ) {
			if ( this->m_decomposeInPlaceFlag )
				m_decompositionFunction ( this->get_givenMatrix(), get_numberOfColumns(), get_columnPointers() );
			else {
				m_skylineMatrix->set_content(); std::copy(this->get_givenMatrix(), this->get_givenMatrix() +m_skylineMatrix->size(), m_skylineMatrix->data() );
				m_decompositionFunction( m_skylineMatrix->data(), get_numberOfColumns(), get_columnPointers() ); 
			}
		}
	}
	else {
		m_skylineMatrix->set_fromFullMatrix(matrix);
		m_decompositionFunction( m_skylineMatrix->get_skylineMatrix(), get_numberOfColumns(), get_columnPointers() );
	}
}

this->m_solverFunction( m_skylineMatrix->get_skylineMatrix() , rhs, sol, get_numberOfColumns(), get_columnPointers() );

}


std::istream& operator >> (std::istream& stream, LinearSolverType& var )
{
	int temp; stream >> temp; 
	var = static_cast<LinearSolverType> ( temp );
	return stream;
}



std::string get_LinearSolverType_name( LinearSolverType& type )
{
	std::string name;
	switch (type) {
		case NonPredefinedSolver: { name = "Non predefined"; break;  }
		case gaussEliminationWithBackSubst: { name = "Gauss elimination"; break;  }
		case ConjugateGradient: { name = "Conjugate gradient"; break;  }
		case LU_Doolittle: { name = "LU Doolittle"; break;  }
		case LU_Crout: { name = "LU Crout"; break;  }
		case LDU: { name = "LDU"; break;  }
		case LDL: { name = "LDL"; break;  }
		case Cholesky: { name = "Cholesky"; break;  }
		case CholeskyBanded: { name ="Cholesky banded"; break;  }
		case CholeskySkyline:{ name ="Cholesky skyline"; break;  } 
	}
	return name;
}


std::ostream& operator << ( std::ostream& out, LinearSolverType& type )
{
	out << get_LinearSolverType_name( type );
	return out;
}



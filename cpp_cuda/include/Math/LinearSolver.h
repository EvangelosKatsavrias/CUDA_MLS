#ifndef LINEARSOLVERHEADER
#define LINEARSOLVERHEADER

#include"gaussElimination.h"
#include"conjugateGradient.h"
#include"SystemConditionEvaluator.h"
#include"LUDecomposers.h"
#include"LUSolvers.h"
#include"matrix_memoryStorage.h"
#include<vector>
#include<iostream>


enum LinearSolverType{ NonPredefinedSolver, gaussEliminationWithBackSubst, ConjugateGradient, LU_Doolittle, LU_Crout, LDU, LDL, Cholesky, CholeskyBanded, CholeskySkyline };
std::ostream& operator << (std::ostream& out, LinearSolverType& type);
std::istream& operator >> (std::istream& in, LinearSolverType& type);
std::string get_LinearSolverType_name( LinearSolverType& type );


template<class T>
class LinearSolverBase
{
protected:
	int 				m_numCols;
	SystemConditionEvaluator<T>* 	m_systemConditionEvaluator;
	T				m_conditionNumber;
	LinearSolverType		m_linearSolverType;
	T* 				m_givenMatrix{0};

	void(*m_solverFunction)(T* matrix, T* rhs, T* sol, int n);

public:
	LinearSolverBase(int numOfColumns=2);

	void virtual 	set_numberOfColumns	(int numOfCols);
	void virtual	set_solverFunction 	(void(*)(T* A,T* b,T* x,int n) );
	void virtual	set_solverFunction 	(void(*)(T* L, T* U, T* b,T* x,int n) ) {}
	void virtual	set_solverFunction 	(void(*)(T* L,T* D,T* U,T* b,T* x,int n) ) {}
	void virtual	set_solverFunction 	(void(*)(T* L,T* b,T* x,int n, int n_b) ) {}
	void virtual	set_solverFunction 	(void(*)(T* L,T* b,T* x,int n, int* columnPointers) ) {}
	void		set_solverType		(LinearSolverType);
	int  virtual	get_numberOfColumns	() {return m_numCols;}

	void virtual	set_decompositionFunction  ( void(*)(T*, T*, int) ) {}
	void virtual	set_decompositionFunction  ( void(*)(T*, T*, T*, int) ) {}
	void virtual	set_decompositionFunction  ( void(*)(T*, T*, T*, T*, int) ) {}
	void virtual	set_decompositionFunction  ( void(*)(T*, T*, int, int) ) {}
	void virtual	set_decompositionFunction  ( void(*)(T*, int, int*) ) {}
	void virtual 	write_solverInfo(std::ostream& out ){ out << m_linearSolverType << std::endl; }
	int  virtual	reset () {}

	int 		set_givenMatrix(T* matrix) { m_givenMatrix = matrix; }
	T* 		get_givenMatrix() { return m_givenMatrix; }
	
	void ( *get_solverFunction() ) 	(T*,T*,T*,int) {return m_solverFunction;}
	void virtual operator () (T* matrix, T* rhs, T* sol=0);

	
//	static LinearSolver<T>* createLinearSolver (LinearSolverType solverType, int numOfColumns);

};

template class LinearSolverBase<float>;
template class LinearSolverBase<double>;


template<class T>
class GaussEliminationSolver: public LinearSolverBase<T>
{

public:
	GaussEliminationSolver( int numOfColumns=2 ): LinearSolverBase<T>(numOfColumns) { this->set_solverType( gaussEliminationWithBackSubst ); }

	void write_solverInfo(std::ostream& out ){ out << "Gauss elimination solver\n"; }
};
template class GaussEliminationSolver<float>;
template class GaussEliminationSolver<double>;


template<class T>
class ConjugateGradientSolver: public LinearSolverBase<T>
{

protected:

public:
	ConjugateGradientSolver( int numOfColumns=2 ): LinearSolverBase<T>(numOfColumns){ this->set_solverType( ConjugateGradient ); }

	void write_solverInfo(std::ostream& out ){ out << "Conjugate gradient solver\n"; }
};
template class ConjugateGradientSolver <float>;
template class ConjugateGradientSolver <double>;


template<class T>
class LUSolver: public LinearSolverBase<T>
{

protected:
	std::vector<T> 	m_internalLowerMatrix;
	std::vector<T> 	m_internalUpperMatrix;

	void (*m_decompositionFunction) (T* matrix, T* L, T* U, int n);
	void(*m_solverFunction)(T* lowerMatrix, T* upperMatrix, T* rhs, T* sol, int n);

public:
	LUSolver ( int numOfColumns=2): LinearSolverBase<T>(numOfColumns) { this->set_solverType( LU_Doolittle ); }

	void set_decompositionFunction  ( void(*)(T*, T*, T*, int) );
	void ( *get_decompositionFunction() ) ( T*, T*, T*, int ) { return m_decompositionFunction; }
	void set_solverFunction 	( void(*)(T*, T*, T*, T*, int) );
	void ( *get_solverFunction() ) 	(T*, T*, T*, T*, int) { return m_solverFunction; }
	T* get_lowerMatrix() { return m_internalLowerMatrix.data(); }
	T* get_upperMatrix() { return m_internalUpperMatrix.data(); }

	void write_solverInfo(std::ostream& out ){ out << "Solver with LU factorization\n"; }
	int reset() { this->set_givenMatrix(0); }

	void operator () (T* matrix, T* rhs, T* sol=0);
	void operator () (T* lowerMatrix, T* upperMatrix, T* rhs, T* sol);
};
template class LUSolver <float>;
template class LUSolver <double>;


template<class T>
class LUCroutSolver: public LUSolver<T>
{

public:
	using LUSolver<T>::set_decompositionFunction;
	using LUSolver<T>::set_solverFunction;
	LUCroutSolver ( int numOfColumns=2): LUSolver<T>(numOfColumns) { this->set_solverType( LU_Crout ); }

	int reset() { this->set_givenMatrix(0); }
	void write_solverInfo(std::ostream& out ){ out << "Solver with LU Crout factorization\n"; }
};
template class LUCroutSolver <float>;
template class LUCroutSolver <double>;


template<class T>
class LDUSolver: public LinearSolverBase<T>
{

protected:
	std::vector<T> m_internalLowerMatrix;
	std::vector<T> m_internalUpperMatrix;
	std::vector<T> m_internalDiagonal;

	void (*m_decompositionFunction) (T* matrix, T* L, T* U, T* D, int n);
	void (*m_solverFunction)(T* lowerMatrix, T* upperMatrix, T* D, T* rhs, T* sol, int n);

public:
	LDUSolver ( int numOfColumns=2): LinearSolverBase<T>(numOfColumns){ this->set_solverType( LDU ); }

	void set_solverFunction 	( void(*)(T*, T*, T*, T*, T*, int) );
	void set_decompositionFunction  ( void(*)(T*, T*, T*, T*, int) );
	void ( *get_decompositionFunction() ) ( T*, T*, T*, T*, int ) { return m_decompositionFunction; }
	void ( *get_solverFunction() ) 	(T*, T*, T*, T*, T*, int) { return m_solverFunction; }

	int reset() { this->set_givenMatrix(0); }
	void write_solverInfo(std::ostream& out ){ out << "Solver with LDU factorization\n"; }
	void operator () (T* matrix, T* rhs, T* sol=0);
	void operator () (T* lowerMatrix, T* upperMatrix, T* D, T* rhs, T* sol);
};
template class LDUSolver <float>;
template class LDUSolver <double>;


template<class T>
class LDLSolver: public LinearSolverBase<T>
{

protected:
	std::vector<T> m_internalLowerMatrix;
	std::vector<T> m_internalDiagonal;

	void (*m_decompositionFunction) (T* matrix, T* L, T* D, int n);
	void (*m_solverFunction)(T* lowerMatrix, T* D, T* rhs, T* sol, int n);

public:
	LDLSolver ( int numOfColumns=2): LinearSolverBase<T>(numOfColumns){ this->set_solverType( LDL ); }

	void set_solverFunction 	( void(*)(T*, T*, T*, T*, int) );
	void set_decompositionFunction  ( void(*)(T*, T*, T*, int) );
	void ( *get_decompositionFunction() ) ( T*, T*, T*, int ) { return m_decompositionFunction; }
	void ( *get_solverFunction() ) 	(T*, T*, T*, T*, int) { return m_solverFunction; }

	int reset() { this->set_givenMatrix(0); }
	void write_solverInfo(std::ostream& out ){ out << "Solver with LDL factorization\n"; }

	void operator () (T* matrix, T* rhs, T* sol=0);
	void operator () (T* lowerMatrix, T* D, T* rhs, T* sol);
};
template class LDLSolver <float>;
template class LDLSolver <double>;


template<class T>
class CholeskySolver: public LinearSolverBase<T>
{

protected:
	bool 	m_decomposedInputFlag;
	bool	m_decomposeInPlaceFlag;
	std::vector<T> m_decomposedMatrix;

	void (*m_decompositionFunction) (T* matrix, T* L, int n);

public:
	CholeskySolver ( int numOfColumns=2, bool decomposedInputFlag = 0, bool decomposeInPlaceFlag = 0 ): LinearSolverBase<T>(numOfColumns), m_decomposedInputFlag(decomposedInputFlag), m_decomposeInPlaceFlag(decomposeInPlaceFlag)
	{ this->set_solverType(Cholesky ); }

	void set_decompositionFunction  ( void(*)(T*, T*, int) );
	void ( *get_decompositionFunction() ) ( T*, T*, int ) { return m_decompositionFunction; }
	
	int reset() { this->set_givenMatrix(0); }
	void write_solverInfo(std::ostream& out ){ out << "Solver with Cholesky factorization\n"; }
	void operator () (T* matrix, T* rhs, T* sol=0);

};
template class CholeskySolver <float>;
template class CholeskySolver <double>;


template<class T>
class CholeskyBandedSolver: public virtual CholeskySolver<T>
{

protected:
	BandedMatrix<T>*	m_bandedMatrix;
	bool			m_bandedInputFlag;
	void (*m_decompositionFunction) (T* matrix, T* L, int n, int n_b);
	void (*m_solverFunction)(T* lowerMatrix, T* rhs, T* sol, int n, int n_b);

public:
	CholeskyBandedSolver ( int numOfColumns=2, bool bandedInputFlag=0, int numOfBands=1, bool decomposedInputFlag = 0, bool decomposeInPlaceFlag = 1 ): CholeskySolver<T>(numOfColumns, decomposedInputFlag, decomposeInPlaceFlag), m_bandedMatrix( new BandedMatrix<T>(numOfColumns, numOfBands) ), m_bandedInputFlag(bandedInputFlag)
	{ this->set_solverType( CholeskyBanded );  }

	~CholeskyBandedSolver ( ) { delete m_bandedMatrix; }


	void set_solverFunction 	( void(*)(T*, T*, T*, int, int) );
	void set_decompositionFunction  ( void(*)(T*, T*, int, int) );
	void ( *get_decompositionFunction() ) ( T*, T*, int, int ) { return m_decompositionFunction; }
	void ( *get_solverFunction() ) 	(T*, T*, T*, int, int) { return m_solverFunction; }

	void write_solverInfo(std::ostream& out ){ out << "Solver with Cholesky factorization in banded matrix\n"; }

	void operator () (T* matrix, T* rhs, T* sol=0) { operator()(matrix, rhs, sol, 0, 0); }
	void operator () (T* matrix, T* rhs, T* sol, int numOfColumns, int numOfBands=0);

	int reset() { this->set_givenMatrix(0); }

	void set_numberOfColumns(int numOfColumns ) { m_bandedMatrix->set_numOfColumns(numOfColumns); this->m_numCols = numOfColumns; }
	int get_numberOfColumns() { return m_bandedMatrix->get_numOfColumns(); }
	int set_bandedInputFlag(bool flag) { m_bandedInputFlag = flag; }
	bool get_bandedInputFlag( ) { return m_bandedInputFlag; }
	int get_numOfBands( ) { return m_bandedMatrix->get_numOfBands(); }
	int set_numOfBands(int numOfBands ) { m_bandedMatrix->set_numOfBands( numOfBands );  }
//	int set_bandedMatrix(BandedMatrix<T>* bmatrix ) { m_bandedMatrix = bmatrix; }
//	BandedMatrix<T>* get_bandedMatrix() { return m_bandedMatrix; }
	
};
template class CholeskyBandedSolver <float>;
template class CholeskyBandedSolver <double>;


template<class T>
class CholeskySkylineSolver: public virtual CholeskySolver<T>
{

protected:
	SkylineMatrix<T>*	m_skylineMatrix;
	bool			m_skylineInputFlag;
	void (*m_decompositionFunction) (T* matrix, int n, int* columnPointers);
	void (*m_solverFunction)(T* lowerMatrix, T* rhs, T* sol, int n, int* columnPointers);

public:
	CholeskySkylineSolver ( int numOfColumns=2, bool skylineInputFlag=0, int* columnPointers=0, bool decomposedInputFlag = 0, bool decomposeInPlaceFlag = 1 ): CholeskySolver<T>(numOfColumns, decomposedInputFlag, decomposeInPlaceFlag), m_skylineMatrix( new SkylineMatrix<T>(numOfColumns, columnPointers) ), m_skylineInputFlag(skylineInputFlag)
	{ if ( !decomposeInPlaceFlag ) this->m_decomposedMatrix.resize( this->m_numCols ); this->set_solverType( CholeskySkyline );  }

	~CholeskySkylineSolver () { delete m_skylineMatrix; }

	void set_solverFunction 	( void(*)(T*, T*, T*, int, int*) );
	void set_decompositionFunction  ( void(*)(T*, int, int*) );
	void ( *get_decompositionFunction() ) ( T*, int, int* ) { return m_decompositionFunction; }
	void ( *get_solverFunction() ) 	(T*, T*, T*, int, int*) { return m_solverFunction; }

	void write_solverInfo(std::ostream& out ){ out << "Solver with Cholesky factorization in skyline matrix\n"; }

	void operator () (T* matrix, T* rhs, T* sol=0 ) { operator()(matrix, rhs, sol, 0, 0); }
	void operator () (T* matrix, T* rhs, T* sol, int numberOfColumns, int* columnPointers );

	int reset() { this->set_givenMatrix(0); }
	void set_numberOfColumns(int numOfColumns ) { m_skylineMatrix->set_numOfColumns(numOfColumns); this->m_numCols = numOfColumns; }
	int get_numberOfColumns() { return m_skylineMatrix->get_numOfColumns(); }
	int set_skylineInputFlag(bool flag) { m_skylineInputFlag = flag; }
	bool get_skylineInputFlag( ) { return m_skylineInputFlag; }
	int* get_columnPointers ( ) { return m_skylineMatrix->get_columnPointers(); }
	int set_columnPointers (int* columnPointers ) { m_skylineMatrix->set_columnPointers( columnPointers );  }
//	int set_skylineMatrix(SkylineMatrix<T>* skymatrix ) { m_skylineMatrix = skymatrix; }
//	SkylineMatrix<T>* get_skylineMatrix() { return m_skylineMatrix; }
	
};
template class CholeskySkylineSolver <float>;
template class CholeskySkylineSolver <double>;



template<class T>
class LinearSolver
{

private:
	LinearSolverBase<T>* m_solver;


public:
	LinearSolver( LinearSolverType solverType = gaussEliminationWithBackSubst, int numOfColumns = 2 )
	{ m_solver = createLinearSolver( solverType, numOfColumns ); }


	void operator () (T* matrix, T* rhs, T* sol=0)
	{ m_solver->operator()( matrix, rhs, sol); }


	void write_solverInfo( std::ostream& out ){ m_solver->write_solverInfo(out); }


	void set_numberOfColumns (int numOfCols) { m_solver->set_numberOfColumns(numOfCols); }
	int reset () { m_solver->reset(); }
	int set_givenMatrix(T* matrix) { m_solver->set_givenMatrix(matrix); }
	T*  get_givenMatrix() { return m_solver->get_givenMatrix(); }
	

	static LinearSolverBase<T>* createLinearSolver (LinearSolverType solverType, int numOfColumns)
	{
		switch ( solverType ) {
			case gaussEliminationWithBackSubst: 	{ return new GaussEliminationSolver<T>( numOfColumns ); break;}
			case ConjugateGradient: 		{ return new ConjugateGradientSolver<T>( numOfColumns ); break;}
			case LU_Doolittle: 			{ return new LUSolver<T> ( numOfColumns ); break;}
			case LU_Crout: 				{ return new LUCroutSolver<T> ( numOfColumns ); break;}
			case LDU: 				{ return new LDUSolver<T> ( numOfColumns ); break;}
			case LDL: 				{ return new LDLSolver<T> ( numOfColumns ); break;}
			case Cholesky: 				{ return new CholeskySolver<T> ( numOfColumns); break; }
			case CholeskyBanded: 			{ return new CholeskyBandedSolver<T> (numOfColumns); break;}
			case CholeskySkyline: 			{ return new CholeskySkylineSolver<T>(numOfColumns); break;}
	}
}

};


#endif

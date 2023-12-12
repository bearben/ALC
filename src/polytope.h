#include <iostream>
#include <algorithm>
#include "armadillo"
#include "math.h"
#include "time.h"
#include "memory.h"
#include "boost/random.hpp"

#ifndef POLYVOL_H
#define POLYVOL_H

namespace polyvest{

class polytope{
public:
	polytope(int rows, int cols);
	~polytope();
	double 	matA(double val, int i, int j){ return (A(i, j) = val); }
	double 	matA(int i, int j){ return A(i, j); }
	double 	vecb(double val, int i){ return (b(i) = val); }
	double 	vecb(int i){ return b(i); }
	
	void	GetBounds(bool isIntBound);	// results in xmax and xmin
	void	GetBounds(int index, bool isIntBound);
	bool	Enlarge(double bound);
	bool 	AffineTrans();
	double 	EstimateVol(double epsilon, double delta, double coef);
	double	ExactVol();		// call vinci
	double	ExactCount();	// call barvinok
	double	Optimization(arma::vec &object);
	double	Optimization(arma::vec &object, bool isILP);

	void	PrepForWalk();
	void	Walk();
	arma::vec GetInvPoint(arma::vec point);
	arma::vec GetTransPoint(arma::vec point);
	bool	isInside(arma::vec point);
	bool	CUonBoundary(arma::vec point);
	
	double 	Volume() const { return vol; }
	double 	LatticeCount() const { return latcount; }
	void 	Print();
	polytope* Clone();

	bool 	msg_off;

	//reciprocal of beta, beta-cut
	double beta_r;
	
	//polytope denoted by: Ax<=b, A is an (m x n) matrix.
	int 	m, n;
	arma::mat A;
	arma::vec b;
	int 	*var_flag;	// 0 = Real; 1 = Integer

	// random walk point
	arma::vec x;
	
	// upper & lower bounds for each xi
	double 	*xmax;
	double 	*xmin;
	double 	*hpoffset;
	
	double 	vol;		// volume and det(A)
	double	latcount;	// count of integer points
	
private:
	double 	walk(int k);
	void 	genInitE(double &R2, arma::vec &Ori);

	double 	randd(double u){ return rand() * u / RAND_MAX; }
	int 	randi(int u){ return rand() % u; }

	arma::mat invT;		// to compute invert points
	arma::vec invOri;
	arma::mat applyT;	// to compute trans points
	arma::vec applyOri;
	
	double	determinant;
	int 	l;
	double 	*r2;

	arma::vec *B;
	arma::mat *Ai;
};

inline polytope::polytope(int rows, int cols) :
	msg_off(true),
	m(rows),
	n(cols),
	A(rows, cols),
	b(rows),
	x(cols),
	vol(0),
	invT(cols, cols),
	invOri(cols),
	applyT(cols, cols),
	applyOri(cols),
	determinant(0)
{
	beta_r = 2 * n;

	l = (int)(n * log((double)beta_r) / log((double)2)) + 2;
	r2 = new double[l];
	for (int i = 0; i < l; i++) 
		r2[i] = pow((double)2, (double)(2 * i) / n);

	var_flag = new int[n] ();
	for (int i = 0; i < n; i++) var_flag[i] = 1;

	xmax = new double[n] ();
	xmin = new double[n] ();
	hpoffset = new double[m] ();

	b.zeros();
	A.zeros();
	x.zeros();

	B = new arma::vec[n];
	Ai = new arma::mat[n];
}

inline polytope::~polytope(){
	delete []r2;
	delete []var_flag;
	delete []xmax;
	delete []xmin;
	delete []hpoffset;
	delete []B;
	delete []Ai;
}

}

#endif

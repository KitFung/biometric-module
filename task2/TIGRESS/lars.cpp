#include <string>
#include <vector>
#include <armadillo>
#include "lars.h"

using namespace std;
using namespace arma;

Lars::Lars(mat X, mat Y, int L){
	Beta = lars(X, Y, L);
}


mat Lars::lars(mat X, mat Y, int L, bool lasso){
	double eps = 1e-10;

	int n = X.n_rows;
	int p = X.n_cols;
	X = normalize(X);
	Y = Y - (mat(n, 1, fill::ones) * mean(Y));
	mat t;
	t << numeric_limits<double>::infinity();
	int T = t.size();

	mat beta(1, p, fill::zeros);
	mat mu(n, 1, fill::zeros);
	mat mu_old(n, 1, fill::zeros);

	mat gamma(L, 1, fill::zeros); 
	
	urowvec A;
	urowvec Ac(p);
	for (int a = 0; a < p; a++){
		Ac(a) = a;
	}
	int nVars = 0;
	int signOk = 1;
	uword i = 0;
	int t_prev = 0;

	mat beta_t = zeros(T,p);

	int ii = 0;
	mat tt = t;

	uvec addVar;
	while (nVars < L){
		mat c = trans(X) * (Y - mu);
		double C = as_scalar(max(abs(c)));
		
		if (C < eps || t.is_empty())
			break;

		if (i == 0)
			addVar = find(C == abs(c));
		if (signOk){
			nVars++;
			A.insert_cols(A.size(), addVar);
		}
		mat s_A = sign(c(A));
		Ac = setDiff(Ac, A);
		int nZeros = Ac.n_cols;

		mat X_A = X.cols(A);
		mat G_A = trans(X_A) * X_A;
		mat invG_A = inv(G_A);
		double L_A = as_scalar (1 / sqrt(trans(s_A) * invG_A * s_A));
		mat w_A = L_A * invG_A * s_A;
		mat u_A = X_A * w_A;
		mat a = trans(X) * u_A;
		mat beta_tmp(p, 1, fill::zeros);
		mat gammaTest(nZeros, 2, fill::zeros);

		if (nVars == L){
			gamma(i) = C / L_A;
		}
		else {
			for (uword j = 0; j < nZeros; j++){
				double jj = Ac(j);
				gammaTest(j, 0) = (C - c(jj)) / (L_A - a(jj));
				gammaTest(j, 1) = (C + c(jj)) / (L_A + a(jj));
			}
			uvec min_i; 
			urowvec min_j;

			gamma(i) = minplus(gammaTest, min_i, min_j);
			addVar = unique(trans(Ac.cols(min_i)));

		}
		for (uword o = 0; o < A.n_cols; o++){
			beta_tmp(A(o)) = beta(i, A(o)) + (gamma(i) * w_A(o));
		}
		
		/*if (lasso){
			signOk = 1;
			urowvec a(1);
			a.fill(i);
			gammaTest = - trans(beta.submat(a, A)) / w_A;
			
			urowvec min_i;
			uvec min_j;

			double gamma2 = minplus(gammaTest, min_i, min_j);

			if (gamma2 < gamma(i)){
				gamma(i) = gamma2;
				for (uword o = 0; o < A.n_cols; o++){
					beta_tmp(A(o)) = beta(i, A(o)) + (gamma(i) * w_A(o));
				}
				beta_tmp(A.cols(unique(min_i))) = 0;
				A.shed_col(minIndex(0));
				nVars--;
				signOk = 0;
			}
		}*/
		
		if (t(0) != numeric_limits<double>::infinity()){
			double t_now = norm(beta_tmp(A), 1);

			if (t_prev < t(0) && t_now >= t(0)){
				for (int o = 0; 0 < A.n_cols; o++){
					beta_t(ii, A(o)) = beta(i, A(o)) + L_A * (t(0) - t_prev) * w_A(o);
				}
				ii++;
			}
			t_prev = t_now;
		}

		mu = mu_old + gamma(i) * u_A;
		mu_old = mu;
		beta.insert_rows(beta.n_rows, trans(beta_tmp));

		i++;
	}
	return beta;
}



urowvec Lars::setDiff(urowvec a, urowvec b){
	urowvec maxList;
	maxList << a.max() << b.max();
	urowvec sortMat(maxList.max() + 1, fill::zeros);
	urowvec result;
	
	for (int i = 0; i < a.n_cols; i++){
		if (sortMat(a(i)) != 0)
			continue;
		sortMat(a(i)) = 1;
	}

	for (int i = 0; i < b.n_cols; i++){
		if (sortMat(b(i)) == 0)
			continue;
		sortMat(b(i)) = 0;
	}

	result = trans(find(sortMat));
	return result;
}

double Lars::minplus(mat X, uvec &i, urowvec &j){

	for (int col = 0; col < X.n_cols; col++){
		uvec a = find(imag(X.col(col)) != 0 );
		if (a.is_empty()) continue;
		for (int row = 0; row < a.n_rows; row++){
			X(a(row), col) = numeric_limits<double>::infinity();
		}
	}
	for (int col = 0; col < X.n_cols; col++){
		uvec a = find(X.col(col) <= 0);
		if (a.is_empty()) continue;
		for (int row = 0; row < a.n_rows; row++){
			X(a(row), col) = numeric_limits<double>::infinity();
		}
	}
	vector<int> tmpi, tmpj;
	double minValue = min(min(X));
	for (uword row = 0; row < X.n_rows; row++){
		for (uword col = 0; col < X.n_cols; col++){
			if (X(row, col) == minValue){
				tmpi.push_back(row);
				tmpj.push_back(col);
			}
		}
	}
	i = conv_to<uvec>::from(tmpi);
	j = conv_to<urowvec>::from(tmpj);

	return minValue;
}

mat Lars::normalize(mat X){
	int n = X.n_rows;
	X = (eye<mat>(n, n) - ((1 / (n * 1.0)) * ones<mat>(n, n))) * X;
	mat n22 = sqrt(diagvec(trans(X) * X));
	X = X / (ones<mat>(n, 1) * trans(n22));
	return X;
}

mat Lars::getBeta(){
	return Beta;
}
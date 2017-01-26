#ifndef TIGRESS_H
#define TIGRESS_H
#include <armadillo>

class Tigress{
	private:
		static arma::urowvec removeTG(int, arma::urowvec);
		
		static arma::mat stabilitySelection(arma::mat, arma::mat, int, int, double alpha, arma::uvec, int);
	public:
		static arma::cube tigress(int, int, double, arma::mat, arma::urowvec);
		static arma::mat normalize_data(arma::mat);

};
#endif

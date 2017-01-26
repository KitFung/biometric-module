#ifndef LARS_H
#define LARS_H
#include <armadillo>

class Lars{
	private:
		arma::mat Beta;
	public:
		Lars(arma::mat, arma::mat, int);
		arma::mat lars(arma::mat, arma::mat, int, bool=false);
		arma::urowvec setDiff(arma::urowvec, arma::urowvec);
		double minplus(arma::mat, arma::uvec &, arma::urowvec &);
		arma::mat normalize(arma::mat);
		arma::mat getBeta();

};

#endif
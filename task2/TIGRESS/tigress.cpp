#include <string>
#include <vector>
#include <armadillo>
#include "tigress.h"
#include "lars.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"

using namespace std;
using namespace arma;
using namespace tbb;



cube Tigress::tigress(int R, int L, double alpha, mat expdata, urowvec tfIndex){
	expdata = normalize_data(expdata);
	int ntf = tfIndex.n_cols;
	cube freq(ntf, L, expdata.n_cols);
	
	for (uword i = 0; i < expdata.n_cols;i++){
	//parallel_for(uword(0), expdata.n_cols, [&](uword i){
		urowvec TFofTG = removeTG(i, tfIndex);
		mat x = expdata.cols(TFofTG);	
		mat y = expdata.col(i);

		mat tmp = stabilitySelection(x, y, R, L, alpha, find(TFofTG >= 0), ntf); 
		freq.slice(i) = trans(tmp);
	}//);
	return freq;
}

mat Tigress::normalize_data(mat expdata){
	mat dataMean = mean(expdata);
	mat dataStd = stddev(expdata);

	return (expdata - (mat(expdata.n_rows, 1, fill::ones) * dataMean)) * diagmat(1 / dataStd);
	

}

urowvec Tigress::removeTG(int index, urowvec tfIndex){
	for(uword i = 0; i < tfIndex.size(); i++){
		if(tfIndex(i) != index) continue;
		tfIndex.shed_col(i);
	}
	return tfIndex;
}

mat Tigress::stabilitySelection(mat X, mat Y, int R, int L, double alpha, uvec predictorTF, int ntf){
	mat freq(L, ntf, fill::zeros);
	arma_rng::set_seed_random();
	for(int i = 0; i < floor(R/2); i++){
		
		mat reweightedData = X % repmat(alpha + (1 - alpha) * randu(1, X.n_cols), X.n_rows, 1);

		uvec listExperiment(X.n_rows);
		for(int a = 0; a < X.n_rows; a++ ){
			listExperiment(a) = a;
		}
		srand(NULL);
		random_shuffle(listExperiment.begin(), listExperiment.end());	

		int half = floor(X.n_rows / 2);
		Lars l1(reweightedData.rows(listExperiment.subvec(0, half - 1)), Y.rows(listExperiment.subvec(0, half - 1)), L);
		Lars l2(reweightedData.rows(listExperiment.subvec(half, listExperiment.n_rows - 1)), Y.rows(listExperiment.subvec(half, listExperiment.n_rows - 1)), L);

		mat result1 = l1.getBeta();
		mat result2 = l2.getBeta();

		freq.cols(predictorTF) = freq.cols(predictorTF) + abs(sign(result1(span(1, L), span())));
		freq.cols(predictorTF) = freq.cols(predictorTF) + abs(sign(result2(span(1, L), span())));

	}
	return freq / (R * 1.0);
}

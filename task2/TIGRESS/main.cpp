#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include "tbb/task_scheduler_init.h"
#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "tigress.h"

using namespace std;
using namespace arma;
using namespace tbb;

mat score_edges(cube, int);
mat rank_edges(mat, urowvec, int = 0);

int main(){
	ofstream outputDependency("Output-Dependency.txt");

	cout << "--------------- TIGRESS Method ---------------" << endl << endl;
	cout << "------------- 1. Retrieving Data -------------" << endl;

	string dataset, tfset;

	//Request input file name from user and store it into variable inputTxt
	cout << "Input Dataset: ";
	cin >> dataset;

	//Request input file name from user and store it into variable inputTxt
	cout << "Input TF: ";
	cin >> tfset;


	//Declare an input stream for reading data from inputTxt
	ifstream input(dataset);

	//Declare an input stream for reading data from inputTxt
	ifstream input2(tfset);

	//This variable is used to keep all genes' name
	vector<string> geneName;

	vector<string> tfList;

	//This variable is used to keep all genes' data
	mat genesData;


	string tfName;
	while (getline(input2,tfName)){
		tfList.push_back(tfName);
	}

	//This portion of code is used to move the pointer to next line in order to read a correct row of data
	string geneNameStr;
	getline(input, geneNameStr);

	geneNameStr.erase(geneNameStr.begin(), geneNameStr.begin() + geneNameStr.find(',') + 1);

	while (true){
		int curP = geneNameStr.find(',');
		if (curP == -1){
			geneName.push_back(geneNameStr);
			break;
		}
		geneName.push_back(geneNameStr.substr(0, curP));
		geneNameStr.erase(geneNameStr.begin(), geneNameStr.begin() + curP + 1);
	}



	//This variable is used to store input stream string
	string line;

	//Read the data line by line
	//This portion of code is used to decode the format of the data from .csv extension
	//First, get the gene's name and store into a vector
	//Then, get the data of particular gene and store into a vector
	mat exp;
	while (getline(input, line)){

		int curP = line.find(',');

		line.erase(line.begin(), line.begin() + curP + 1);

		vector<double> geneData;
		while (true){
			int curP = line.find(',');
			if (curP == -1){
				geneData.push_back(stod(line));
				break;
			}
			geneData.push_back(stod(line.substr(0, curP)));
			line.erase(line.begin(), line.begin() + curP + 1);
		}
		exp.insert_rows(exp.n_rows, conv_to <rowvec>::from(geneData));
	}
    cout << "------------- 2. Setting Parameters ----------" << endl;


	int R = 500;
	int L = 3;
	double alpha = 0.3;

	cout << "R\t= " << R << endl;
	cout << "L\t= " << L << endl;
	cout << "Alpha\t= " << alpha << endl << endl << endl;

	cout << "------------- 3. Inferring Network -----------" << endl;
	urowvec tfIndex;
	tfIndex <<0 << 1 << 2 << 3;

	cube freq = Tigress::tigress(R, L, alpha, exp, tfIndex);
	freq.print();
	cout << "------------- 4. Edges' Score ---------------" << endl;
	mat scores = score_edges(freq, 2);

	cout << "------------- 5. Edges' Ranking -------------" << endl;
	mat sortedEdges = rank_edges(scores, tfIndex, 0);

	cout << "------------- 6. Output Result ---------------" << endl;
	//output ranked important regulators
	for (int i = 0; i < sortedEdges.n_rows; i++){
		outputDependency << geneName.at(sortedEdges(i, 0)) << "\t" << geneName.at(sortedEdges(i, 1)) << "\t" << sortedEdges(i, 2)<<endl;
	}

	outputDependency.close();
	input.close();
	input2.close();
	return 0;
}

mat score_edges(cube freq, int method){
	int ntf = freq.n_rows;
	int L = freq.n_cols;
	int ngene = freq.n_slices;
	mat scores;
	switch (method){
	case 1:
		scores = reshape(max(freq, 1), ntf, ngene, 1);
		break;

	case 2:
		scores = reshape((zeros(ntf, 1, ngene) + freq(span(), span(0,0), span())) / 2 +
					sum((freq(span(), span(1, L - 1), span()) + freq(span(), span(0, L - 2), span())) / 2, 1), ntf, ngene, 1);
		scores = 1 / (L - 0.5) * scores;
		break;
	}
	return scores;
}

mat rank_edges(mat scores, urowvec tfIndex, int cutoff){
	int ntf = tfIndex.n_cols;
	int ngene = scores.n_cols;
	int k = 0;
	mat edges((ngene * ntf) - ntf, 3);
	mat sortedEdges((ngene * ntf) - ntf, 3);

	for (int i = 0; i < ntf; i++){
		for (int j = 0; j < ngene; j++){
			if (tfIndex(i) != j){
				edges(k, 0) = tfIndex(i);
				edges(k, 1) = j;
				edges(k, 2) = scores(tfIndex(i), j);
				k++;
			}
		}
	}
	uvec indexList = sort_index(edges.col(2), "descend");
	for (int i = 0; i < indexList.n_rows; i++){
		sortedEdges.row(i) = edges.row(indexList(i));
	}

	sortedEdges = sortedEdges.rows(find(sortedEdges.col(2) > 0));
	int nedges = sortedEdges.n_rows;

	if (cutoff != 0)
		sortedEdges = sortedEdges.head_rows(cutoff);


	return sortedEdges;
}

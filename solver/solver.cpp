#include <cfloat>
#include <cmath>
#include <vector>
#include <cstdio>

double eps = 100; //accuracy
double UPB, LOB, deltaL, glob = DBL_MAX; //obtained upped bound, lower bound, mean of delta*L, and real global minimizer

int nodes = 3, dim = 2; //amount of nodes per dimension, amount of dimensions in task

double(*compute)(double x1, double x2) = nullptr; //pointer to function that computes y at x[]
//double(*compute)(double *x) = nullptr; 

struct part {
	double *a;
	double *b;
	double LocUP, LocLO, deltaL;
	part(const double* toa, const double* tob, int ldim) {
		a = (double*)malloc(ldim * sizeof(double));
		b = (double*)malloc(ldim * sizeof(double));
		for (int i = 0; i < ldim; i++) {
			a[i] = toa[i];
			b[i] = tob[i];
		}
	}
};

std::vector<part> Partitions;

double getR(double delta) {
	return exp(delta);
}

void UpdateRecords(double LU, double LL) {
	if (LU < UPB) {
		UPB = LU;
		LOB = LL;
	}
}

int ChooseDim(double *a, double *b) { //выбор координаты для деления гиперинтервала
	double max = DBL_MIN, cr;
	int maxI = 0;
	for (int i = 0; i < dim; i++) {
		cr = fabs(b[i] - a[i]);
		if (cr > max) {
			max = cr;
			maxI = i;
		}
	}
	return maxI;
}

void GridEvaluator(double *a, double *b, double *Frp, double *LBp, double *dL) {
	double Fr = DBL_MAX, L = DBL_MIN, R, Fx, Fx1, LB;
	double curR, curdelta;
	int edge, change, flag;
	//code here
	double *step, *x, *x1;
	step = (double*)malloc(dim * sizeof(double));
	x = (double*)malloc(dim * sizeof(double));
	x1 = (double*)malloc(dim * sizeof(double));
	for (int i = 0; i < dim; i++) {
		step[i] = fabs(b[i] - a[i]) / (nodes - 1);
	}
	curR = getR(step[0]);
	curdelta = step[0];
	for (int i = 0; i < dim; i++) { //i - размерность, по которой считаем соседние точки
		for (int k = 0; k < dim; k++) {
			x[k] = a[k];
			x1[k] = a[k];
		}
		R = getR(step[i]);
		edge = i == 0 ? 1 : 0;
		change = i == (dim - 1) ? (dim - 2) : (dim - 1);
		int prechange = change; double loc;
		while (!(x[edge] > b[edge])) { //проход по всем точкам
			change = prechange;
			flag = 0;
			for (int j = 1; j < nodes; j++) {// from here
				x1[i] += step[i];
				Fx = compute(x[0], x[1]);
				Fx1 = compute(x1[0], x1[1]);
				Fr = Fx < Fr ? Fx : Fr;
				loc = fabs(Fx - Fx1) / step[i];
				if (loc > L) {
					L = loc;
					curR = R;
					curdelta = step[i];
				}
				x[i] += step[i];
			}//to here вычисления функции
			Fr = Fx1 < Fr ? Fx1 : Fr;//учитываем последнее вычисление
									 //from here
			while (!flag) {
				x[change] += step[change];
				if (change == edge) {
					flag = 1;
					continue;
				}
				if (x[change] > b[change]) {
					x[change] = a[change];
					change--;
					continue;
				}
				flag = 1;
			}//to here генерация новой точки
		}
	}
	LB = curR * L * curdelta;
	//code here
	*dL = LB;
	LB = Fr - LB;
	free(step); free(x); free(x1);
	*Frp = Fr;
	*LBp = LB;
}

int do_test(double *a, double* b) {
	UPB = DBL_MAX;
	if (glob == DBL_MAX) {
		return -1;
	}
	if (compute == nullptr) {
		return -2;
	}

	Partitions.emplace_back(a, b, dim);

	double *a1, *b1;

	a1 = (double*)malloc(dim * sizeof(double));
	b1 = (double*)malloc(dim * sizeof(double));

	if ((a1 == nullptr) || (b1 == nullptr)) {
		return -3;
	}

	while (Partitions.size() != 0) {
		unsigned int parts = Partitions.size();

		for (int i = 0; i < static_cast<int>(parts); i++) {
			double lUPB, lLOB, ldeltaL;
			GridEvaluator(Partitions[i].a, Partitions[i].b, &lUPB, &lLOB, &ldeltaL);
			Partitions[i].LocLO = lLOB;
			Partitions[i].LocUP = lUPB;
			Partitions[i].deltaL = ldeltaL;
		}

		for (unsigned int i = 0; i < parts; i++) {
			UpdateRecords(Partitions[i].LocUP, Partitions[i].LocLO);
		}

		for (unsigned int i = 0; i < parts; i++) {
			if ((Partitions[i].LocLO < (UPB - eps)) || (Partitions[i].deltaL > eps)) {
				int choosen = ChooseDim(Partitions[i].a, Partitions[i].b);

				for (int j = 0; j < dim; j++) { //make new edges for 2 new hyperintervals:
					if (j != choosen) {						//[a .. b1] [a1 .. b]
						a1[j] = Partitions[i].a[j];			//where a1 = [a[1], a[2], .. ,a[choosen] + b[choosen]/2, .. , a[dim] ]
						b1[j] = Partitions[i].b[j];			//and b1 = [b[1], b[2], .. ,a[choosen] + b[choosen]/2, .. , b[dim] ]
					}
					else {
						a1[j] = Partitions[i].a[j] + fabs(Partitions[i].b[j] - Partitions[i].a[j]) / 2.0;
						b1[j] = a1[j];
					}
				}
				try {
					Partitions.emplace_back(Partitions[i].a, b1, dim);
					Partitions.emplace_back(a1, Partitions[i].b, dim);
				}
				catch (std::bad_alloc e) {
					free(a1); free(b1);
					return -3;
				}
			}
		}

		for (unsigned int i = 0; i < parts; i++)
			Partitions.erase(Partitions.begin());
		Partitions.shrink_to_fit(); //to economy memory
	}
	free(a1); free(b1);
	Partitions.clear();
	return 0;
}

void checkErrors(int code) {
	switch (code) {
	case -1:
		fprintf(stderr, "glob value have not been initialized!\n");
		break;
	case -2: 
		fprintf(stderr, "pointer to computing function have not been initialized!\n");
		break;
	case -3:
		fprintf(stderr, "Sorry, amount of RAM on your device insufficient to solve this task, please upgrade :)\n");
		break;
	default:
		break;
	}
}
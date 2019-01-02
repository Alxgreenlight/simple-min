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
	int size;
	double *a = nullptr;
	double *b = nullptr;
	double LocUP, LocLO, deltaL;
};


struct Partitions {
	int size;
	int cur_alloc;
	struct part* base;
	Partitions() {
		size = 0;
		cur_alloc = 0;
		base = nullptr;
	}
	~Partitions() {
		for (int i = 0; i < size; i++) {
			free(base[i].a);
			free(base[i].b);
		}
		size = 0;
		if (cur_alloc != 0)
			free(base);
		cur_alloc = 0;;
	}
	void erase() {
		for (int i = 0; i < size; i++) {
			free(base[i].a);
			free(base[i].b);
		}
		size = 0;
		if (cur_alloc != 0)
			free(base);
		cur_alloc = 0;
	}
	void add(double* toa, double *tob, int dim) {
		if (cur_alloc == 0) {
			base = (struct part*)malloc(16 * sizeof(struct part));
			if (base) {
				cur_alloc = 16;
			}
			else {
				fprintf(stderr,"Error alloc\n");
				erase();
				return;
			}
		}
		if (size == cur_alloc) {
			base = (struct part*)realloc(base,(cur_alloc + 16) * sizeof(struct part));
			if (base) {
				cur_alloc += 16;
			}
			else {
				fprintf(stderr, "Error alloc\n");
				erase();
				return;
			}
		}
		base[size].size = dim;
		base[size].a = (double*)malloc(dim * sizeof(double));
		base[size].b = (double*)malloc(dim * sizeof(double));
		if ((!base[size].a)||(!base[size].b)) {
			fprintf(stderr, "Error alloc\n");
			erase();
			return;
		}
		for (int i = 0; i < dim; i++) {
			base[size].a[i] = toa[i];
			base[size].b[i] = tob[i];
		}
		size++;
	}
	part & operator[](int n) {
		if (n > size - 1) {
			fprintf(stderr, "Out of size\n");
			n = 0;
		}
		return base[n];
	}
	Partitions & operator=(Partitions& P) {
		erase();
		for (int i = 0; i < P.size; i++) {
			add(P[i].a, P[i].b, P[i].size);
		}
		return (*this);
	}
};


Partitions P, P1;

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

	P.add(a, b, dim);

	double *a1, *b1;

	a1 = (double*)malloc(dim * sizeof(double));
	b1 = (double*)malloc(dim * sizeof(double));

	if ((a1 == nullptr) || (b1 == nullptr)) {
		return -3;
	}

	while (P.size != 0) {
		unsigned int parts = P.size;

		for (int i = 0; i < static_cast<int>(parts); i++) {
			double lUPB, lLOB, ldeltaL;
			GridEvaluator(P[i].a, P[i].b, &lUPB, &lLOB, &ldeltaL);
			P[i].LocLO = lLOB;
			P[i].LocUP = lUPB;
			P[i].deltaL = ldeltaL;
		}

		for (unsigned int i = 0; i < parts; i++) {
			UpdateRecords(P[i].LocUP, P[i].LocLO);
		}

		for (unsigned int i = 0; i < parts; i++) {
			if ((P[i].LocLO < (UPB - eps)) || (P[i].deltaL > eps)) {
				int choosen = ChooseDim(P[i].a, P[i].b);

				for (int j = 0; j < dim; j++) { //make new edges for 2 new hyperintervals:
					if (j != choosen) {						//[a .. b1] [a1 .. b]
						a1[j] = P[i].a[j];			//where a1 = [a[1], a[2], .. ,a[choosen] + b[choosen]/2, .. , a[dim] ]
						b1[j] = P[i].b[j];			//and b1 = [b[1], b[2], .. ,a[choosen] + b[choosen]/2, .. , b[dim] ]
					}
					else {
						a1[j] = P[i].a[j] + fabs(P[i].b[j] - P[i].a[j]) / 2.0;
						b1[j] = a1[j];
					}
				}
					P1.add(P[i].a, b1, dim);
					P1.add(a1, P[i].b, dim);
			}
		}

		P = P1;
		P1.erase();
	}
	free(a1); free(b1);
	P.erase();
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
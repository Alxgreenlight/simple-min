#include <cfloat>
#include <cmath>
#include <vector>
#include <cstdio>


double eps = 100; //accuracy (100 is default value)
double UPB, LOB, deltaL, glob = DBL_MAX; //obtained upped bound, lower bound, value of delta*L, and real global minimizer
double sqrt2 = sqrt(2.0);

int nodes = 3, dim = 2; //amount of nodes per dimension, amount of dimensions in task

//double(*compute)(double x1, double x2) = nullptr; //pointer to function that computes y at x[]
double(*compute)(double *x) = nullptr; // this pointer should be uncommented and previous commented in case of function for evaluation
//takes an array x in parameter
//ALSO LOOK AT LINES 163 AND 165!!!

struct part {	//description of hyperinterval
	int size;
	double *a = nullptr;
	double *b = nullptr;
	double LocUP, LocLO, deltaL;
};


struct Partitions {		//all hyperintervals obtained on current step of algorithm
	int size;									//like std::vector, but in C language with only needed methods
	int cur_alloc; //memory already allocated for array of struct part
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
	int add(double* toa, double *tob, int dim) {
		if (cur_alloc == 0) {
			base = (struct part*)malloc(16 * sizeof(struct part));
			if (base) {
				cur_alloc = 16;
			}
			else {
				fprintf(stderr,"Error alloc\n");
				erase();
				return -1;
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
				return -1;
			}
		}
		base[size].size = dim;
		base[size].a = (double*)malloc(dim * sizeof(double));
		base[size].b = (double*)malloc(dim * sizeof(double));
		if ((!base[size].a)||(!base[size].b)) {
			fprintf(stderr, "Error alloc\n");
			erase();
			return -1;
		}
		for (int i = 0; i < dim; i++) {
			base[size].a[i] = toa[i];
			base[size].b[i] = tob[i];
		}
		size++;
		return 0;
	}
	part & operator[](int n) {
		if (n > size - 1) {
			fprintf(stderr, "Out of size\n");
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


Partitions P, P1; //P - hyperintervals on current step, P1 - hyperintervals on next step

double getR(double delta) { //r value for formula for etimating Lipschitz constant
	return 2.0*exp(delta);
}

void UpdateRecords(double LU, double LL) { //UPB - global upper bound, LU - local value of upper bound on current hyperinterval
	if (LU < UPB) {
		UPB = LU;
		LOB = LL;
	}
}

int ChooseDim(double *a, double *b) { //selection of coordinate for hyperinterval division
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
	double *step, *x, *x1;
	step = (double*)malloc(dim * sizeof(double)); //step of grid in every dimension
	x = (double*)malloc(dim * sizeof(double)); //algorithm now needs to evaluate function in two adjacent point at the same time
	x1 = (double*)malloc(dim * sizeof(double));
	for (int i = 0; i < dim; i++) {
		step[i] = fabs(b[i] - a[i]) / (nodes - 1);
	}
	curR = getR(step[0]);	//value of r for this dimension
	curdelta = step[0];		//delta value (step of grid) for this dimension
	for (int i = 0; i < dim; i++) { //i - размерность, по которой считаем соседние точки
		for (int k = 0; k < dim; k++) {
			x[k] = a[k];
			x1[k] = a[k];
		}
		R = getR(step[i]);
		edge = i == 0 ? 1 : 0;		//[0..n], there may be [i(fixed),1,..,n] or [0,..,n-1,i(fixed)]
		change = i == (dim - 1) ? (dim - 2) : (dim - 1);	//so we need edge and change values(change is lower digit for iterating over points)
		int prechange = change; double loc;
		while (!(x[edge] > b[edge])) { //all pount pass
			change = prechange;
			flag = 0;
			x[i] = a[i];
			x1[i] = a[i];
			for (int j = 1; j < nodes; j++) {// from here
				x1[i] += step[i];
				//Fx = compute(x[0],x[1]);
				Fx = compute(x);
				//Fx1 = compute(x1[0],x1[1]);
				Fx1 = compute(x1);
				Fr = Fx < Fr ? Fx : Fr;
				loc = fabs(Fx - Fx1) / step[i];
				if (loc > L) {
					L = loc;
					curR = R;
					curdelta = step[i];
				}
				x[i] += step[i];
			}//to here function evaluation
			//from here
			while (!flag) {
				if (change == i) {
					change--;
					continue;
				}
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
			}//to here new point generation
		}
	}
	//final calculation
	LB = curR * L * curdelta;
	*dL = LB;
	LB = Fr - LB;
	free(step); free(x); free(x1);
	*Frp = Fr;
	*LBp = LB;
}

int min_search(double *a, double* b) {
	UPB = DBL_MAX;	//upper bound
	if (glob == DBL_MAX) {	//check if uninitialized values
		return -1;
	}
	if (compute == nullptr) {
		return -2;
	}

	P.add(a, b, dim);	//add first hyperinterval

	double *a1, *b1;
	//if hyperinterval divides, 2 new hyperintervals with bounds [a..b1] [a1..b] creates
	a1 = (double*)malloc(dim * sizeof(double));
	b1 = (double*)malloc(dim * sizeof(double));

	if ((a1 == nullptr) || (b1 == nullptr)) {
		return -3;
	}

	while (P.size != 0) {	//each hyperinterval can be subdivided or pruned (if non-promisiable or fits accuracy)
		unsigned int parts = P.size;

		for (unsigned int i = 0; i < parts; i++) {
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
			if (!((P[i].LocLO > (UPB - eps)) || (P[i].deltaL < eps))) {
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
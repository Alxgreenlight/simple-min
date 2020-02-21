#include <iostream>
#include <fstream>
#include <limits>
#include <chrono>
#include <iomanip>
#include <iterator>
#include "mathexplib/testfuncs/benchmarks.hpp"
#include "solver/R_Optim_Pure_Parallel.hpp"
#include "util/helper.hpp"


using BM = Benchmark<double>;
BM *ptr;
std::chrono::steady_clock sc;
double maxdiff = std::numeric_limits<double>::min();
std::chrono::milliseconds atime_span = std::chrono::duration_values<std::chrono::milliseconds>::zero();
int summary;
//int sumef = 0;
unsigned long long int ae = 0;
unsigned long int ai = 0;
std::ofstream fp, adif;
double calc(const double *x)
{
	std::vector<double> xc;
	xc.assign(x, x + ptr->getDim());
	double v = ptr->calcFunc(xc);
	return v;
}
void findMin(const BM& bm, int nodes, double eps) {
	int dim = bm.getDim();
	double* a = new double[dim];
	double* b = new double[dim];
	double* x = new double[dim];
	for (int i = 0; i < dim; i++) {
		a[i] = bm.getBounds()[i].first;
		b[i] = bm.getBounds()[i].second;
	}
	PPrOptimizer<double> rOpt;
	double UPB, L;
	try {
		rOpt.init(dim, eps);
		auto start = sc.now();
		UPB = rOpt.search(dim, x, a, b, calc);
		auto end = sc.now();

		atime_span += std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
		auto time_span = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	
	double diff = 1.0*fabs(bm.getGlobMinY() - UPB);
	maxdiff = diff > maxdiff ? diff : maxdiff;
	unsigned long long int fevals;
	unsigned long int iters;
	rOpt.getInfo(fevals, iters, L);
	ae += fevals;
	ai += iters;
	fp << "Search complete:" << std::endl << "Upper bound: " << std::setprecision(4) << UPB << std::endl << \
		"At [";
	std::copy(x, x + dim, std::ostream_iterator<double>(fp, " "));
	fp << "]" << std::endl << "Real global minimum: " << std::setprecision(4) << bm.getGlobMinY() << std::endl << std::endl << "Diff: " << std::setprecision(4) << diff << std::endl;
	fp << "Function evaluations: " << fevals << " in " << iters << " iterations" << std::endl;
	fp << "With L = " << L << std::endl;
	fp << "Evaluation time: " << time_span.count() << "ms" << std::endl;
	if (diff > eps) {
		summary++;
	}
	//gs.printLinfo(adif);
	fp << std::endl << std::endl;
	fp << std::endl << std::endl;
	}
	catch (std::exception &e) {
		std::cerr << "Something wrong: " << e.what() << std::endl;
		fp << e.what() << std::endl;
		delete[]a; delete[]b; delete[]x;
		return;
	}
	delete[]a; delete[]b; delete[]x;
	return;
}

void testBench(const BM& bm, int nodes, double eps) {
	fp << "*************Testing benchmark**********" << std::endl;
	fp << bm;
	//adif << bm.getDesc() << ';';
	findMin(bm, nodes, eps);
	fp << "****************************************" << std::endl << std::endl;
	return;
}

int main(int argc, char* argv[]) {
	fp.open("mathexp_L_R_study.txt", std::ios::app);
	if (!fp.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}
	/*adif.open("Just_info.csv", std::ios::app);
	if (!adif.is_open()) {
		std::cerr << "Problem with file..." << std::endl;
		std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		std::cin.get();
		return -1;
	}*/

	//int nodes;
	double eps;

	/*std::cout << "Set number of nodes per dimension" << std::endl << \
		"It can significantly affect on performance!" << std::endl;
	std::cin >> nodes;
	while (std::cin.fail()) {
		std::cerr << "Please, repeat input" << std::endl;
		std::cin >> nodes;
	}*/

	if (argc < 2)
    {
        helper::help(benchlib);
        return 0;
    }

	eps = atof(argv[1]);
	if (eps < std::numeric_limits<double>::min()){
		std::cerr << "Accuracy defined incorrect, exit" << std::endl;
		return -1;
	}
	//adif << "Search with " << nodes << " nodes; eps = " << eps << ';' << std::endl;
	fp << "Search with " /*<< nodes */ << "fix nodes; eps = " << eps << ';' << std::endl;

	const int sum = 37;
	/*Benchmarks<double> tests;
	for (auto bm : tests) {
		ptr = &(*bm);
		std::cout << bm->getDesc();
		testBench(*bm, nodes, eps);
		std::cout << '|';
		sum++;
	}
	std::cout << std::endl;*/

	/*Ackley2Benchmark<double> ack2(3);
	ptr = &ack2;
	std::cout << ack2.getDesc();
	testBench(ack2, 0, eps);
	std::cout << "|";

	Alpine2Benchmark<double> Alp2;
	ptr = &Alp2;
	std::cout << Alp2.getDesc();
	testBench(Alp2, 0, eps);
	std::cout << "|";

	BartelsConnBenchmark<double> Brc;
	ptr = &Brc;
	std::cout << Brc.getDesc();
	testBench(Brc, 0, eps);
	std::cout << "|";

	BealeBenchmark<double> Bel;
	ptr = &Bel;
	std::cout << Bel.getDesc();
	testBench(Bel, 0, eps);
	std::cout << "|";

	BiggsEXP2Benchmark<double> Bx2;
	ptr = &Bx2;
	std::cout << Bx2.getDesc();
	testBench(Bx2, 0, eps);
	std::cout << "|";

	BiggsEXP3Benchmark<double> Bx3;
	ptr = &Bx3;
	std::cout << Bx3.getDesc();
	testBench(Bx3, 0, eps);
	std::cout << "|";

	BirdBenchmark<double> Brd;
	ptr = &Brd;
	std::cout << Brd.getDesc();
	testBench(Brd, 0, eps);
	std::cout << "|";

	Bohachevsky1Benchmark<double> Bh1;
	ptr = &Bh1;
	std::cout << Bh1.getDesc();
	testBench(Bh1, 0, eps);
	std::cout << "|";

	Bohachevsky2Benchmark<double> Bh2;
	ptr = &Bh2;
	std::cout << Bh2.getDesc();
	testBench(Bh2, 0, eps);
	std::cout << "|";

	Bohachevsky3Benchmark<double> Bh3;
	ptr = &Bh3;
	std::cout << Bh3.getDesc();
	testBench(Bh3, 0, eps);
	std::cout << "|";

	BoothBenchmark<double> Bth;
	ptr = &Bth;
	std::cout << Bth.getDesc();
	testBench(Bth, 0, eps);
	std::cout << "|";

	BoxBettsQuadraticSumBenchmark<double> Bbq;
	ptr = &Bbq;
	std::cout << Bbq.getDesc();
	testBench(Bbq, 0, eps);
	std::cout << "|";

	BraninRCOSBenchmark<double> brc;
	ptr = &brc;
	std::cout << brc.getDesc();
	testBench(brc, 0, eps);
	std::cout << "|";

	BraninRCOS2Benchmark<double> Brc2;
	ptr = &Brc2;
	std::cout << Brc2.getDesc();
	testBench(Brc2, 0, eps);
	std::cout << "|";
	
	BrentBenchmark<double> Brnt;
	ptr = &Brnt;
	std::cout << Brnt.getDesc();
	testBench(Brnt, 0, eps);
	std::cout << "|";

	Bukin2Benchmark<double> Bk2;
	ptr = &Bk2;
	std::cout << Bk2.getDesc();
	testBench(Bk2, 0, eps);
	std::cout << "|";

	Bukin4Benchmark<double> Bk4;
	ptr = &Bk4;
	std::cout << Bk4.getDesc();
	testBench(Bk4, 0, eps);
	std::cout << "|";

	CamelSixHumpBenchmark<double> Cm6;
	ptr = &Cm6;
	std::cout << Cm6.getDesc();
	testBench(Cm6, 0, eps);
	std::cout << "|";

	CamelThreeHumpBenchmark<double> Cm3;
	ptr = &Cm3;
	std::cout << Cm3.getDesc();
	testBench(Cm3, 0, eps);
	std::cout << "|";

	ComplexBenchmark<double> Comp;
	ptr = &Comp;
	std::cout << Comp.getDesc();
	testBench(Comp, 0, eps);
	std::cout << "|";

	CosineMixtureBenchmark<double> Csm;
	ptr = &Csm;
	std::cout << Csm.getDesc();
	testBench(Csm, 0, eps);
	std::cout << "|";

	CubeBenchmark<double> Cub;
	ptr = &Cub;
	std::cout << Cub.getDesc();
	testBench(Cub, 0, eps);
	std::cout << "|";

	DeckkersAartsBenchmark<double> Dka;
	ptr = &Dka;
	std::cout << Dka.getDesc();
	testBench(Dka, 0, eps);
	std::cout << "|";

	DixonPriceBenchmark<double> Dxp;
	ptr = &Dxp;
	std::cout << Dxp.getDesc();
	testBench(Dxp, 0, eps);
	std::cout << "|";

	DropWaveBenchmark<double> Dpw;
	ptr = &Dpw;
	std::cout << Dpw.getDesc();
	testBench(Dpw, 0, eps);
	std::cout << "|";

	EggCrateBenchmark<double> egg;
	ptr = &egg;
	std::cout << egg.getDesc();
	testBench(egg, 0, eps);
	std::cout << "|";

	ElAttarVidyasagarDuttBenchmark<double> elav;
	ptr = &elav;
	std::cout << elav.getDesc();
	testBench(elav, 0, eps);
	std::cout << "|";

	EngvallBenchmark<double> engv;
	ptr = &engv;
	std::cout << engv.getDesc();
	testBench(engv, 0, eps);
	std::cout << "|";

	Exp2Benchmark<double> exp2;
	ptr = &exp2;
	std::cout << exp2.getDesc();
	testBench(exp2, 0, eps);
	std::cout << "|";

	ExponentialBenchmark<double> Ex(3);
	ptr = &Ex;
	std::cout << Ex.getDesc();
	testBench(Ex, 0, eps);
	std::cout << "|";

	GoldsteinPriceBenchmark<double> Gdp;
	ptr = &Gdp;
	std::cout << Gdp.getDesc();
	testBench(Gdp, 0, eps);
	std::cout << "|";

	GramacyLee2Benchmark<double> Gm2;
	ptr = &Gm2;
	std::cout << Gm2.getDesc();
	testBench(Gm2, 0, eps);
	std::cout << "|";

	GramacyLee3Benchmark<double> Gm3;
	ptr = &Gm3;
	std::cout << Gm3.getDesc();
	testBench(Gm3, 0, eps);
	std::cout << "|";

	HansenBenchmark<double> han;
	ptr = &han;
	std::cout << han.getDesc();
	testBench(han, 0, eps);
	std::cout << "|";

	Hartman3Benchmark<double> Ht3;
	ptr = &Ht3;
	std::cout << Ht3.getDesc();
	testBench(Ht3, 0, eps);
	std::cout << "|";*/

	HelicalValleyBenchmark<double> Hv;
	ptr = &Hv;
	std::cout << Hv.getDesc();
	testBench(Hv, 0, eps);
	std::cout << "|"; 

	/*HimmelblauBenchmark<double> Hmb;
	ptr = &Hmb;
	std::cout << Hmb.getDesc();
	testBench(Hmb, 0, eps);
	std::cout << "|";

	HosakiBenchmark<double> Hsk;
	ptr = &Hsk;
	std::cout << Hsk.getDesc();
	testBench(Hsk, 0, eps);
	std::cout << "|";

	JennrichSampsonBenchmark<double> jsb;
	ptr = &jsb;
	std::cout << jsb.getDesc();
	testBench(jsb, 0, eps);
	std::cout << "|";

	LeonBenchmark<double> leon;
	ptr = &leon;
	std::cout << leon.getDesc();
	testBench(leon, 0, eps);
	std::cout << "|";

	MatyasBenchmark<double> mat;
	ptr = &mat;
	std::cout << mat.getDesc();
	testBench(mat, 0, eps);
	std::cout << "|";

	McCormickBenchmark<double> mcc;
	ptr = &mcc;
	std::cout << mcc.getDesc();
	testBench(mcc, 0, eps);
	std::cout << "|";

	Mishra5Benchmark<double> msh5;
	ptr = &msh5;
	std::cout << msh5.getDesc();
	testBench(msh5, 0, eps);
	std::cout << "|";

	ParsopoulosBenchmark<double> psp;
	ptr = &psp;
	std::cout << psp.getDesc();
	testBench(psp, 0, eps);
	std::cout << "|";

	PeriodicBenchmark<double> per;
	ptr = &per;
	std::cout << per.getDesc();
	testBench(per, 0, eps);
	std::cout << "|";

	Price1Benchmark<double> prc1;
	ptr = &prc1;
	std::cout << prc1.getDesc();
	testBench(prc1, 0, eps);
	std::cout << "|";

	Price2Benchmark<double> prc2;
	ptr = &prc2;
	std::cout << prc2.getDesc();
	testBench(prc2, 0, eps);
	std::cout << "|";

	Price3Benchmark<double> prc3;
	ptr = &prc3;
	std::cout << prc3.getDesc();
	testBench(prc3, 0, eps);
	std::cout << "|";

	Problem02Benchmark<double> pr2;
	ptr = &pr2;
	std::cout << pr2.getDesc();
	testBench(pr2, 0, eps);
	std::cout << "|";

	Problem04Benchmark<double> pr4;
	ptr = &pr4;
	std::cout << pr4.getDesc();
	testBench(pr4, 0, eps);
	std::cout << "|";

	Problem05Benchmark<double> pr5;
	ptr = &pr5;
	std::cout << pr5.getDesc();
	testBench(pr5, 0, eps);
	std::cout << "|";

	Problem06Benchmark<double> pr6;
	ptr = &pr6;
	std::cout << pr6.getDesc();
	testBench(pr6, 0, eps);
	std::cout << "|";

	QingBenchmark<double> qi;
	ptr = &qi;
	std::cout << qi.getDesc();
	testBench(qi, 0, eps);
	std::cout << "|";

	QuadraticBenchmark<double> qd;
	ptr = &qd;
	std::cout << qd.getDesc();
	testBench(qd, 0, eps);
	std::cout << "|";

	QuinticBenchmark<double> qnt(3);
	ptr = &qnt;
	std::cout << qnt.getDesc();
	testBench(qnt, 0, eps);
	std::cout << "|"; 

	RosenbrockBenchmark<double> rsb(3);
	ptr = &rsb;
	std::cout << rsb.getDesc();
	testBench(rsb, 0, eps);
	std::cout << "|";

	RosenbrockModifiedBenchmark<double> rsbm;
	ptr = &rsbm;
	std::cout << rsbm.getDesc();
	testBench(rsbm, 0, eps);
	std::cout << "|";

	RotatedEllipseBenchmark<double> re;
	ptr = &re;
	std::cout << re.getDesc();
	testBench(re, 0, eps);
	std::cout << "|";

	RotatedEllipse2Benchmark<double> re2;
	ptr = &re2;
	std::cout << re2.getDesc();
	testBench(re2, 0, eps);
	std::cout << "|";

	SchwefelBenchmark<double> sw(3);
	ptr = &sw;
	std::cout << sw.getDesc();
	testBench(sw, 0, eps);
	std::cout << "|";

	Schwefel1_2Benchmark<double> sw12(3);
	ptr = &sw12;
	std::cout << sw12.getDesc();
	testBench(sw12, 0, eps);
	std::cout << "|";

	ShubertBenchmark<double> sbb;
	ptr = &sbb;
	std::cout << sbb.getDesc();
	testBench(sbb, 0, eps);
	std::cout << "|";

	Shubert2Benchmark<double> sbb2;
	ptr = &sbb2;
	std::cout << sbb2.getDesc();
	testBench(sbb2, 0, eps);
	std::cout << "|";

	Shubert3Benchmark<double> sbb3;
	ptr = &sbb3;
	std::cout << sbb3.getDesc();
	testBench(sbb3, 0, eps);
	std::cout << "|";

	SphereBenchmark<double> sf(3);
	ptr = &sf;
	std::cout << sf.getDesc();
	testBench(sf, 0, eps);
	std::cout << "|";

	StyblinskiTangBenchmark<double> stb;
	ptr = &stb;
	std::cout << stb.getDesc();
	testBench(stb, 0, eps);
	std::cout << "|";

	SumSquaresBenchmark<double> ssq(3);
	ptr = &ssq;
	std::cout << ssq.getDesc();
	testBench(ssq, 0, eps);
	std::cout << "|";

	TrecanniBenchmark<double> trc;
	ptr = &trc;
	std::cout << trc.getDesc();
	testBench(trc, 0, eps);
	std::cout << "|";

	Trigonometric1Benchmark<double> trg1(3);
	ptr = &trg1;
	std::cout << trg1.getDesc();
	testBench(trg1, 0, eps);
	std::cout << "|";

	Ursem1Benchmark<double> urs;
	ptr = &urs;
	std::cout << urs.getDesc();
	testBench(urs, 0, eps);
	std::cout << "|";

	WWavyBenchmark<double> ww(3);
	ptr = &ww;
	std::cout << ww.getDesc();
	testBench(ww, 0, eps);
	std::cout << "|";

	WayburnSeader1Benchmark<double> wbs;
	ptr = &wbs;
	std::cout << wbs.getDesc();
	testBench(wbs, 0, eps);
	std::cout << "|";

	XinSheYang4Benchmark<double> xs4(3);
	ptr = &xs4;
	std::cout << xs4.getDesc();
	testBench(xs4, 0, eps);
	std::cout << "|";

	ZakharovBenchmark<double> zh(3);
	ptr = &zh;
	std::cout << zh.getDesc();
	testBench(zh, 0, eps);
	std::cout << "|";

	ZettlBenchmark<double> zt;
	ptr = &zt;
	std::cout << zt.getDesc();
	testBench(zt, 0, eps);
	std::cout << "|";

	ZirilliBenchmark<double> zi;
	ptr = &zi;
	std::cout << zi.getDesc();
	testBench(zi, 0, eps);
	std::cout << "|";*/




	fp << "Summary: " << summary << " of " << sum << std::endl;
	fp << "Total evaluation time: " << atime_span.count() << "ms" << std::endl;
	fp << "Maximum difference in this set: " << maxdiff << std::endl;
	fp << "Avg. evaluations: " << ae * 1.0 / sum << std::endl;
	fp << "Avg. iterations: " << ai * 1.0 / sum << std::endl;
	std::cout << "Total evaluation time: " << atime_span.count() << "ms" << std::endl;
	//fp << "Not so effectively: " << sumef << std::endl;

	fp.close();
	//adif.close();
	std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
	std::cin.get();
}

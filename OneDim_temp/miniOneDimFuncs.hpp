#ifndef MINIONEDIMFUNCS_HPP
#define MINIONEDIMFUNCS_HPP


#include <cmath>
#include <vector>
#include <string>
#include <memory>

# define MY_PI          3.141592653589793238462643383279502884L /* pi */

/*

All functions from http://infinity77.net/global_optimization/test_functions_1d.html

*/

namespace onedimopt {
	template <class T>
	class BasicFunc {
	protected:
		T RealX, RealY;
		T LB;
		T RB;
		std::string desc;
	public:
		BasicFunc(T l = -1.0, T r = 1.0) :LB(l), RB(r) {}
		virtual T calculate(const T coord) = 0;
		virtual std::string detDesc() const {
			return desc;
		}
		T RealminX() {
			return RealX;
		}
		T RealminY() {
			return RealY;
		}
		T LBound() {
			return LB;
		}
		T RBound() {
			return RB;
		}
	};
	template <class T>
	class Func1 : public BasicFunc<T> {
	public:
		Func1(const T l = 2.7, const T r = 7.5) :BasicFunc<T>(l, r) {
			this->RealX = 5.145735;
			this->RealY = -1.899599;
			this->desc = "Problem 02";
		}
		virtual T calculate(const T coord) {
			return sin(coord) + sin(10.0 / 3.0*coord);
		}
	};
	template <class T>
	class Func2 : public BasicFunc<T> {
	public:
		Func2(const T l = 1.9, const T r = 3.9) :BasicFunc<T>(l, r) {
			this->RealX = 2.868034;
			this->RealY = -3.85045;
			this->desc = "Problem 04";
		}
		virtual T calculate(const T coord) {
			return -1.0*exp(-1.0*coord)*(16.0*coord*coord - 24 * coord + 5);
		}
	};
	template <class T>
	class Func3 : public BasicFunc<T> {
	public:
		Func3(const T l = 0.0, const T r = 1.2) :BasicFunc<T>(l, r) {
			this->RealX = 0.96609;
			this->RealY = -1.48907;
			this->desc = "Problem 05";
		}
		virtual T calculate(const T coord) {
			return -1.0*sin(18.0*coord)*(1.4 - 3 * coord);
		}
	};

	template <class T>
	class Func4 : public BasicFunc<T> {
	public:
		Func4(const T l = 2.7, const T r = 7.5) :BasicFunc<T>(l, r) {
			this->RealX = 5.19978;
			this->RealY = -1.6013;
			this->desc = "Problem 07";
		}
		virtual T calculate(const T coord) {
			return sin(coord) + sin(10.0 / 3.0*coord) + log(coord) - 0.84*coord + 3;
		}
	};

	template <class T>
	class Func5 : public BasicFunc<T> {
	public:
		Func5(const T l = 0.0, const T r = 4.0) :BasicFunc<T>(l, r) {
			this->RealX = 0.224885;
			this->RealY = -0.788685;
			this->desc = "Problem 14";
		}
		virtual T calculate(const T coord) {
			return -1.0*exp(-1.0*coord)*sin(2 * MY_PI*coord);
		}
	};

	template <class T>
	class BMsOneDim {
		std::vector<std::shared_ptr<BasicFunc<T>>> v;
		void add(const std::shared_ptr<BasicFunc<T>> &ptr)
		{
			v.push_back(ptr);
		}
	public:
		BMsOneDim() {
			add(std::make_shared<Func1<double>>());
			add(std::make_shared<Func2<double>>());
			add(std::make_shared<Func3<double>>());
			add(std::make_shared<Func4<double>>());
			add(std::make_shared<Func5<double>>());

		}
		~BMsOneDim() {
			v.clear();
		}
		typename std::vector<std::shared_ptr<BasicFunc<T>>>::iterator begin()
		{
			return v.begin();
		}
		typename std::vector<std::shared_ptr<BasicFunc<T>>>::iterator end()
		{
			return v.end();
		}
		typename std::vector<std::shared_ptr<BasicFunc<T>>>::const_iterator begin() const
		{
			return v.begin();
		}
		typename std::vector<std::shared_ptr<BasicFunc<T>>>::const_iterator end() const
		{
			return v.end();
		}
	};


}

#endif

// naumann@stce.rwth-aachen.de

#pragma once

#include "Eigen/Dense"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

namespace ad {

/// tangent type
template<typename T, int n_t=1>
struct tangent_t { 
  /// primal function value
  T v; 
  /// directional derivative
  Eigen::Matrix<T,n_t,1> d; 
  /// access value
  T& value() { return v; }
  /// access tangent
  T& tangent(int i=0) { return d(i); }
  /// standard constructor
  tangent_t() : v(0), d(Eigen::Matrix<T,n_t,1>::Zero()) {}
  /// passive value constructor
  tangent_t(float v) : v(v), d(Eigen::Matrix<T,n_t,1>::Zero()) {}
};

template <typename T>
concept tangentT = requires(T x) {
  x.v; 
  x.d; 
  x.value(); 
  x.tangent(); 
  x.tangent(42); 
};

template<tangentT T> 
T operator+(const T& arg1, const T& arg2) {
  T res;
  res.v=arg1.v+arg2.v;
  res.d=arg2.d+arg1.d;
  return res;
}

template<tangentT T>
T& operator+=(T& lhs, const T& rhs) {
  lhs.v=lhs.v+rhs.v;
  lhs.d=lhs.d+rhs.d;
  return lhs;
}

template<tangentT T> 
T operator-(const T& arg1, const T& arg2) {
  T res;
  res.v=arg1.v-arg2.v;
  res.d=arg1.d-arg2.d;
  return res;
}

template<tangentT T>
T& operator-=(T& lhs, const T& rhs) {
  lhs.v=lhs.v-rhs.v;
  lhs.d=lhs.d-rhs.d;
  return lhs;
}

template<tangentT T>
T operator*(const T& arg1, const T& arg2) {
  T res;
  res.v=arg1.v*arg2.v;
  res.d=arg1.v*arg2.d+arg2.v*arg1.d;
  return res;
}

template<tangentT T>
T& operator*=(T& lhs, const T& rhs) {
  lhs.v=lhs.v*rhs.v;
  lhs.d=lhs.d*rhs.v+lhs.v*rhs.d;
  return lhs;
}

template<tangentT T> 
T operator/(const T& arg1, const T& arg2) {
  T res;
  res.v=arg1.v/arg2.v;
  res.d=arg1.d/arg2.v-arg1.v*arg2.d/(arg2.v*arg2.v);
  return res;
}

template<tangentT T>
T& operator/=(T& lhs, const T& rhs) {
  lhs.v=lhs.v/rhs.v;
  lhs.d=lhs.d/rhs.v-lhs.v*rhs.d/(rhs.v*rhs.v);
  return lhs;
}

template<tangentT T> 
T operator-(const T& arg) {
  T res;
  res.v=-arg.v;
  res.d=-arg.d;
  return res;
}

template<tangentT T>
std::ostream& operator<<(std::ostream& s, const T& arg) {
  s << arg.v;
  return s;
}

/// tape recording by overloading for custom data type
template<typename T, int n_a=1>
struct adjoint_t { 
  /// primal function value
  T v; 
  /// position in dependencies
  int vaddr; 
  /// register independent variable
  inline void register_input();
  /// record dependence and derivative of the result of an operation wrt. an argument 
  inline void record_arg(T) const;
  /// record the result of an operation and the number of its arguments 
  inline void record_res(int);
  /// access value
  T& value() { return v; };
  /// access adjoint
  T& adjoint(int=0);
  /// standard constructor
  adjoint_t() : v(0), vaddr(0) {};
  /// passive v constructor
  adjoint_t(float v) : v(v), vaddr(0) {};

  struct tape {
    static std::vector<int> dependencies;
    static std::vector<T> derivatives;
    static std::vector<Eigen::Matrix<T,n_a,1>> adjoints;
    static int n, adjoints_ctr; // for internal use

    /// reset tape for repeated use (e.g. gradient descent)
    static void reset();
    /// allocate and initialize vector of adjoints
    static void init_adjoints();
    /// initialize vector of adjoints for reuse
    static void reuse_adjoints();
    /// interpret tape
    static void interpret();
    /// prepare tape for visualization by dot
    static void todot(const std::string& filename);
    /// print dependencies, derivatives, and adjoints
    static void print();
  };

};

template <typename T, int n_a>
void adjoint_t<T,n_a>::register_input() {
	assert(tape::adjoints_ctr==tape::n+1);
	vaddr=tape::dependencies.size();
	tape::dependencies.push_back(tape::adjoints_ctr++);
	tape::n++;
}

template <typename T, int n_a>
void adjoint_t<T,n_a>::record_arg(T d) const {
	tape::dependencies.push_back(tape::dependencies[vaddr]);
	tape::derivatives.push_back(d);
}

template <typename T, int n_a>
void adjoint_t<T,n_a>::record_res(int num_args) {
	tape::dependencies.push_back(num_args);
	vaddr=tape::dependencies.size();
	tape::dependencies.push_back(tape::adjoints_ctr++);
}

template <typename T, int n_a>
T& adjoint_t<T,n_a>::adjoint(int i) {
	return tape::adjoints[tape::dependencies[vaddr]](i);
}

/// characteristics of adjoint type for static type deduction
template <typename T>
concept adjointT = requires(T x) {
	x.v;
	x.vaddr;
	x.register_input();
	x.record_arg(42);
	x.record_res(42);
	x.value();
	x.adjoint();
	x.adjoint(42);
};

template<adjointT T>
T operator+(const T& arg1, const T& arg2) {
	T res;
	arg1.record_arg(1);
	arg2.record_arg(1);
	res.record_res(2);
	res.v=arg1.v+arg2.v;
	return res;
}

template<adjointT T>
T operator-(const T& arg1, const T& arg2) {
	T res;
	arg1.record_arg(1);
	arg2.record_arg(-1);
	res.record_res(2);
	res.v=arg1.v-arg2.v;
	return res;
}

template<adjointT T>
T operator-(const T& arg) {
	T res;
	arg.record_arg(-1);
	res.record_res(1);
	res.v=-arg.v;
	return res;
}

template<adjointT T>
T operator*(const T& arg1, const T& arg2) {
	T res;
	arg1.record_arg(arg2.v);
	arg2.record_arg(arg1.v);
	res.record_res(2);
	res.v=arg1.v*arg2.v;
	return res;
}

template<adjointT T>
T operator/(const T& arg1, const T& arg2) {
	T res;
	arg1.record_arg(static_cast<decltype(arg1.v)>(1)/arg2.v);
	arg2.record_arg(-arg1.v/(arg2.v*arg2.v));
	res.record_res(2);
	res.v=arg1.v/arg2.v;
	return res;
}

template<adjointT T>
std::ostream& operator<<(std::ostream& s, const T& arg) {
	s << arg.v;
	return s;
}

template<typename T>
bool operator<(const T& arg1, const T& arg2) {
  return arg1.v<arg2.v;
}

template<typename T>
bool operator<=(const T& arg1, const T& arg2) {
  return arg1.v<=arg2.v;
}

template<typename T>
bool operator>=(const T& arg1, const T& arg2) {
  return arg1.v>=arg2.v;
}

template<typename T>
bool operator>(const T& arg1, const T& arg2) {
  return arg1.v>arg2.v;
}

template<typename T>
bool operator==(const T& arg1, const T& arg2) {
  return arg1.v==arg2.v;
}

template<typename T>
bool operator!=(const T& arg1, const T& arg2) {
  return arg1.v!=arg2.v;
}

template<typename T>
T sin(const T& arg) {
	T res;
	if constexpr (tangentT<T>) {
		res.v=sin(arg.v);
		res.d=cos(arg.v)*arg.d;
	}
	else if constexpr (adjointT<T>) {
		arg.record_arg(cos(arg.v));
		res.record_res(1);
		res.v=sin(arg.v);
	} 
	else {
		res=std::sin(arg);
	}
	return res;
}

template<typename T>
T cos(const T& arg) {
	T res;
	if constexpr (tangentT<T>) {
		res.v=cos(arg.v);
		res.d=-sin(arg.v)*arg.d;
	}
	else if constexpr (adjointT<T>) {
		arg.record_arg(-sin(arg.v));
		res.record_res(1);
		res.v=cos(arg.v);
	} 
	else {
		res=std::cos(arg);
	}
	return res;
}

template<typename T>
T tan(const T& arg) {
	T res;
	if constexpr (tangentT<T>) {
		res.v=tan(arg.v);
		res.d=arg.d/pow(cos(arg.v),2);
	}
	else if constexpr (adjointT<T>) {
	        decltype(arg.v) one=1;
		res.v=tan(arg.v);
		arg.record_arg(one/pow(cos(arg.v),2)); 
		res.record_res(1);
	} 
	else {
		res=std::tan(arg);
	}
	return res;
}

template<typename T>
T exp(const T& arg) {
	T res;
	if constexpr (tangentT<T>) {
		res.v=exp(arg.v);
		res.d=res.v*arg.d;
	}
	else if constexpr (adjointT<T>) {
		res.v=exp(arg.v);
		arg.record_arg(res.v);
		res.record_res(1);
	} 
	else {
		res=std::exp(arg);
	}
	return res;
}

template<typename T>
T sqrt(const T& arg) {
	T res;
	if constexpr (tangentT<T>) {
	        decltype(arg.v) two=2;
		res.v=sqrt(arg.v);
		res.d=arg.d/(two*res.v);
	}
	else if constexpr (adjointT<T>) {
		decltype(arg.v) one=1, two=2;
		res.v=sqrt(arg.v);
		arg.record_arg(one/(two*res.v)); 
		res.record_res(1);
	} 
	else {
		res=std::sqrt(arg);
	}
	return res;
}

template<typename T>
T pow(const T& arg, int n) {
	T res;
	if constexpr (tangentT<T>) {
		res.v=pow(arg.v,n);
		decltype(arg.v) n_=n;
		res.d=n_*pow(arg.v,n-1)*arg.d;
	}
	else if constexpr (adjointT<T>) {
		res.v=pow(arg.v,n);
		arg.record_arg(static_cast<decltype(arg.v)>(n)*pow(arg.v,n-1)); 
		res.record_res(1);
	} 
	else {
		res=std::pow(arg,n);
	}
	return res;
}

template<typename T>
T log10(const T& arg) {
	T res;
	if constexpr (tangentT<T>) {
		res.v=log10(arg.v);
		res.d=arg.d/res.v;
	}
	else if constexpr (adjointT<T>) {
		res.v=log10(arg.v);
		arg.record_arg(1./res.v);
		res.record_res(1);
	} 
	else {
		res=std::log10(arg);
	}
	return res;
}

template<typename T>
T fabs(const T& arg) {
	T res;
	if constexpr (tangentT<T>) {
		res.v=fabs(arg.v);
                if (arg.v>0) {
		  res.d=arg.d;
                } else {
		  res.d=-arg.d;
                } 
	}
	else if constexpr (adjointT<T>) {
		res.v=fabs(arg.v);
                if (arg.v>0) {
		  arg.record_arg(1.);
                } else {
		  arg.record_arg(-1.);
                }
		res.record_res(1);
	} 
	else {
		res=std::fabs(arg);
	}
	return res;
}

/// position zero reserved for passive variables in dependencies
template <typename T, int n_a>
std::vector<int> adjoint_t<T,n_a>::tape::dependencies{std::vector<int>(1,0)};

/// edges are enumerated starting from zero
template <typename T, int n_a>
std::vector<T> adjoint_t<T,n_a>::tape::derivatives; 

template <typename T, int n_a>
std::vector<Eigen::Matrix<T,n_a,1>> adjoint_t<T,n_a>::tape::adjoints;

/// entries 1 to n reserved for independent variables in vector of adjoints
template <typename T, int n_a>
int adjoint_t<T,n_a>::tape::n=0;

/// position zero reserved for passive variables in vector of adjoints
template <typename T, int n_a>
int adjoint_t<T,n_a>::tape::adjoints_ctr=1;

/// reset tape for repeated use (e.g. gradient descent)
template <typename T, int n_a>
void adjoint_t<T,n_a>::tape::reset() {
	dependencies.clear();
	dependencies.push_back(0);
	derivatives.clear();
	n=0;
	adjoints_ctr=1;
}

/// vector of adjoints allocated and intialized after recording for given value of adjoints_ctr
template <typename T, int n_a>
void adjoint_t<T,n_a>::tape::init_adjoints() {
  adjoints=std::vector<Eigen::Matrix<T,n_a,1>>(adjoints_ctr,Eigen::Matrix<T,n_a,1>::Zero());
}

/// vector of adjoints allocated and intialized after recording for given value of adjoints_ctr
template <typename T, int n_a>
void adjoint_t<T,n_a>::tape::reuse_adjoints() {
  for (int i=adjoints.size()-1;i>=0;--i) { 
    adjoints[i]=Eigen::Matrix<T,n_a,1>::Zero();
  }
}

template <typename T, int n_a>
void adjoint_t<T,n_a>::tape::interpret() {
  // virtual address (postion in dependencies; starting from last)
  int dep_ctr=dependencies.size()-1;
  // postion in derivatives (starting from last)
  int der_ctr=derivatives.size()-1;
  // interpret only operations recorded
  while (dep_ctr>n) {
    // virtual address of the result of the operation
    int result_vaddr=dep_ctr--;
    // number of arguments of the operation
    int num_args=dependencies[dep_ctr--];
    // adjoint of result from vector of adjoints; type depends on order of differentiation
    auto result_adj=adjoints[dependencies[result_vaddr]];
    for (int i=0;i<num_args;++i) { // increment adjoints of all arguments
      adjoints[dependencies[dep_ctr]]=
        adjoints[dependencies[dep_ctr]]+derivatives[der_ctr--]*result_adj; // chain rule
      dep_ctr--; // next argument
    }
  }
}

template <typename T, int n_a>
void adjoint_t<T,n_a>::tape::todot(const std::string& filename) {
  using namespace std;
  ofstream out(filename);
  out << "digraph {\nrankdir=LR\n";
  int i=dependencies.size()-1;
  int k=derivatives.size()-1;
  while (i>n) {
    int j=1;
    for (;j<=dependencies[i-1];++j) {
      out << dependencies[i-1-j] << "->" << dependencies[i]
          << " [label=\"" << derivatives[k--] << '\"';
      if (dependencies[i-1-j]==0) out << ",style=dotted";
      out << "];\n";
    }
    i-=j+1;
  }
  out << '}' << endl;
}

template <typename T, int n_a>
void adjoint_t<T,n_a>::tape::print() {
  using namespace std;
  cerr << "dependencies:\n";
  for (const auto& e:dependencies) cerr << e << ' ';
  cerr << "\nderivatives:\n";
  for (const auto& e:derivatives) cerr << e << ' ';
  cerr << "\nadjoints:\n";
  for (const auto& e:adjoints) {
    cerr << '(';
    for (int i=0;i<e.size()-1;++i) {
      cerr << e(i) << ' ';
    }
    cerr << e(e.size()-1) << ") ";
  }
  cerr << endl;
}

}

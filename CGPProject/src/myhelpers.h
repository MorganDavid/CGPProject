#ifndef MYHELPERS
#define MYHELPERS
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <complex>
#include "fftw3.h"
#include "mymatrix.h"

namespace myhelpers {
	template<class M>
	inline std::vector<M> getColFromMatrix(myMatrix<M> mat, const size_t x) {
		std::vector<M> out(mat.height);
		for (int i = 0; i < mat.height; ++i) {
			out[i] = (mat.data[i][x]);
		}
		return out;
	};

	template<class M>
	inline std::vector<M> getColFromMatrix(M** mat, const size_t height, const size_t x) {
		std::vector<M> out(height);
		for (int i = 0; i < height; ++i) {
			out[i] = (mat[i][x]);
		}
		return out;
	};
	/// <summary>
	/// Converts from fftw_complex to std::complex. O(n) time. 
	/// </summary>
	/// <param name="in"></param>
	/// <param name="length"></param>
	/// <returns></returns>
	inline std::vector<std::complex<double>> fftw_complex2std_complex(const fftw_complex* in, const size_t length) {
		std::vector<std::complex<double>> out(length);
		for (int i = 0; i < length; i++) {
			out[i] = std::complex<double>(in[i][0], in[i][1]);
		}
		return out;
	};

	/// <summary>
	/// Return index of k largest numbers in arr. O(n log n)ish
	/// </summary>
	/// <typeparam name="T">Must have defined < and ==, = op. </typeparam>
	/// <param name="arr"></param>
	/// <param name="k"></param>
	/// <returns>array of indicies into arr where the k largest numbrs are stored. descending order</returns>
	template<class T>
	inline std::vector<int> maxk(std::vector<T> arr, int k) {
		std::vector<int> inds(arr.size());
		int i = 0;
		std::iota(inds.begin(), inds.end(), i++); // = range (0:arr.end()).
		std::sort(inds.begin(), inds.end(), [&](int x, int y) {return arr[x] < arr[y]; }); // compare elements of arr, not inds. returns indicies instead.
		std::vector<int> sliced(inds.end() - k, inds.end());
		std::reverse(sliced.begin(), sliced.end()); //  make it descending order
		return sliced;
	};



}
#endif
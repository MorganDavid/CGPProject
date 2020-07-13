#include <iostream>
#include <string>
#include <limits>
#include <vector>
#include <numeric>
#include <algorithm>

namespace myhelpers {
	const char* error_list[] = { "Bad File","error2" };

	template<class M>
	std::vector<M> getColFromMatrix(myMatrix<M> mat, const size_t x) { 
		std::vector<M> out(mat.height);
		for (int i = 0; i < mat.height; ++i) {
			out[i]=(mat.data[i][x]);
		}
		return out;
	}

	/// <summary>
	/// Converts from fftw_complex to std::complex. O(n) time. 
	/// </summary>
	/// <param name="in"></param>
	/// <param name="length"></param>
	/// <returns></returns>
	std::vector<std::complex<double>> fftw_complex2std_complex(const fftw_complex* in, const size_t length) {
		std::vector<std::complex<double>> out(length);
		for (int i = 0; i < length; i++) {
			out[i] = std::complex<double>(in[i][0], in[i][1]);
		}
		return out;
	}

	/// <summary>
	/// Return index of k largest numbers in arr. O(n log n)ish
	/// </summary>
	/// <typeparam name="T">Must have defined < and ==, = op. </typeparam>
	/// <param name="arr"></param>
	/// <param name="k"></param>
	/// <returns>array of indicies into arr where the k largest numbrs are stored. descending order</returns>
	template<class T>
	std::vector<int> maxk(std::vector<T> arr, int k) {
		std::vector<int> inds(arr.size());
		int i = 0;
		std::iota(inds.begin(), inds.end(), i++); // = range (0:arr.end()).
		std::sort(inds.begin(), inds.end(), [&] (int x, int y) {return arr[x] < arr[y]; }); // compare elements of arr, not inds. returns indicies instead.
		std::vector<int> sliced(inds.end() - k, inds.end() );
		std::reverse(sliced.begin(), sliced.end()); //  make it descending order
		return sliced;
	}

	/// <summary>
	/// Handles error e and returns error in error_list.
	/// </summary>
	/// <param name="e"></param>
	/// <param name="extra"></param>
	/// TODO: setup proper error handling, make an std::format equivalent and pass params to this function. 
	void handleError(int e, std::string extra) {
		if (extra == "") extra = "none";
		try
		{
			throw e;
		}
		catch (int e)
		{
			std::cout << error_list[e] << " extra info: " + extra + " " << std::endl;;
		}
	}

	void copyCGPdataSet(struct dataSet from, struct dataSet to) {
		//Allocate new arrays
		
	}
}
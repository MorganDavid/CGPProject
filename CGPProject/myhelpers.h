#include <iostream>
#include <string>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;
namespace myhelpers {
	const char* error_list[] = { "Bad File","error2" };

	template<class M>
	M* getColFromMatrix(matrix<M> data,int x) {
		M* out = new M[data.size1()];
		for (int i = 0; i < data.size1(); ++i) {
			if (x > data.size1()) { throw "ERROR: x columns do not exist in m."; }
			out[i] = data(i,x);
		}
		return out;
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
}
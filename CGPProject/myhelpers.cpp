#include <iostream>

namespace myhelpers {

	// Gets the xth column from matrix m.
	// m must be row majour order, WARNING: this is really inefficient. 
	template <class M>
	M* getColFromMatrix(M** m, int x, int num_rows) {
		M* out = new M[num_rows];
		for (int i = 0; i < num_rows; ++i) {
			if (x>num_rows) { throw "ERROR: x columns do not exist in m."; }
			out[i] = m[i][x];
		}
		return out;
	}
}
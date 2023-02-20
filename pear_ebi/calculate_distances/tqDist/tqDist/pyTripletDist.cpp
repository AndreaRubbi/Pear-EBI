#include "TripletDistanceCalculator.h"

#include <Python.h>
#include <vector>

#ifdef _WIN32
#include <windows.h>

BOOL APIENTRY DllMain(HANDLE hModule, DWORD dwReason, LPVOID lpReserved) {
	return TRUE;
}
#endif

extern "C" {

#ifdef _WIN32
	__declspec(dllexport)
#endif
  unsigned long tripletDistance(const char *filename1, const char *filename2) {
    TripletDistanceCalculator tripletCalc;
    return tripletCalc.calculateTripletDistance(filename1, filename2);
  }

#ifdef _WIN32
	__declspec(dllexport)
#endif
  PyObject *allPairsTripletDistance(const char *filename) {
    TripletDistanceCalculator tripletCalc;
    const std::vector<std::vector<INTTYPE_REST> > &resultVector = tripletCalc.calculateAllPairsTripletDistance(filename);

    Py_Initialize();

    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();

    PyObject *result = PyList_New(0);
    for(std::vector<std::vector<INTTYPE_REST> >::const_iterator it1 = resultVector.begin(); it1 != resultVector.end(); ++it1) {
      PyObject *currentList = PyList_New(0);
      for(std::vector<INTTYPE_REST>::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2) {
	PyObject *i = PyLong_FromLong(*it2);
	PyList_Append(currentList, i);
	Py_DECREF(i);
      }
      PyList_Append(result, currentList);
      Py_DECREF(currentList);
    }

    PyGILState_Release(gstate);

    return result;
  }

#ifdef _WIN32
	__declspec(dllexport)
#endif
  PyObject *pairsTripletDistance(const char *filename1, const char *filename2) {
    TripletDistanceCalculator tripletCalc;
    const std::vector<INTTYPE_REST> &resultVector = tripletCalc.pairs_triplet_distance(filename1, filename2);

    Py_Initialize();

    PyGILState_STATE gstate;
    gstate = PyGILState_Ensure();

    PyObject *result = PyList_New(0);
    for(std::vector<INTTYPE_REST>::const_iterator it = resultVector.begin(); it != resultVector.end(); ++it) {
      PyObject *i = PyLong_FromLong(*it);
      PyList_Append(result, i);
      Py_DECREF(i);
    }

    PyGILState_Release(gstate);

    return result;
  }

}

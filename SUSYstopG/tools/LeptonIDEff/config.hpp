//rcp like config file
//Chad Jarvis

#ifndef config_hpp
#define config_hpp

#include <vector>
#include <string>
#include <iostream>

using namespace std;

class config {
public: 
  vector <vector <string> > parselist;
  int listsize;
  config(string filename);
  //config(char* filename);
  static vector <string> parser(string line);
  string find(string cond, string type_str);
  void error(int rep, string cond);
  double getFloat(string cond);
  double getDouble(string cond);
  int getInt(string cond);
  bool getBool(string cond);
  string getString(string cond);
  vector <int> getIntArray(string cond);
  vector <float> getFloatArray(string cond);
  vector <double> getDoubleArray(string cond);
  vector <string> getStringArray(string cond);
  vector <vector <int> > getIntMatrix(string cond);

};

#endif

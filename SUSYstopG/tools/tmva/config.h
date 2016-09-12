#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <vector>
#include <ctype.h>


using namespace std;

class config {
public:
  vector <vector <string> > parselist;
  int listsize;
  config(string filename);
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

config::config(string filename) {
  vector <string> rcplist;
  vector <string> parsestring;
  ifstream rcpfile(filename.c_str());
  if(rcpfile.fail()) {
    cout << "Cant open file " << filename << endl;
  }
  string line;

  while(getline(rcpfile,line)) {
    rcplist.push_back(line);
  }
  for(unsigned int i=0; i<rcplist.size(); i++) {
    if((rcplist[i].find("//",0))!=string::npos) {
      rcplist[i].erase(rcplist[i].find("//",0));
    }
    parsestring=parser(rcplist[i]);
    if(parsestring.size()==0) continue;
    parselist.push_back(parsestring);
  }
  listsize=parselist.size();
}

vector <string> config::parser(string line) {
  vector <string> parsestring;
  int start=0;
  int finish=0;
  bool word=false;
  
  for(unsigned int i=0; i<line.size(); i++) {
    if(isspace(line[i])) {
      if(word) {
	finish = i;
	parsestring.push_back(line.substr(start,finish-start));	
	word=false;
      }
    }
    else { 
      if(i==(line.size()-1)) {
	finish=i+1;

	if(!word) start=i; 
	parsestring.push_back(line.substr(start,finish-start));
      }
      if(!word) {
	start=i;
	word=true;
      }      
    }
  }
  return parsestring;
}

string config::find(string cond, string type_str) {
  string var;
  int rep=0;
  for(int i=0; i<listsize; i++) {
    if(parselist[i].size()!=3) continue;
    string type=parselist[i][0];
    string tag=parselist[i][1];
    string svar=parselist[i][2];
    if(tag==cond) {
      if(type!=type_str) 
	cout << "type mismatch for " << tag << endl;
      var=svar;
      rep++;
    }
  }
  error(rep,cond);
  return var;
}

double config::getFloat(string cond) {
  return atof(find(cond,"float").c_str());
}

double config::getDouble(string cond) {
  return atof(find(cond,"double").c_str());
}

int config::getInt(string cond) {
  return atoi(find(cond,"int").c_str());
}

bool config::getBool(string cond) {
  string svar;
  bool var=false;
  svar=find(cond,"bool");
  if(svar=="true") var=true;
  else if(svar=="false") var=false;
  else cout << "boolean not set to true or false for " << cond << endl;
  return var;
}
string config::getString(string cond) {
  string astring = find(cond,"string");
  while(astring.find('"')!=string::npos) {
    astring.replace(astring.find('"'), 1, 1, ' ');
  }
  return astring;
}

vector <int> config::getIntArray(string cond) {
  vector <int> array;
  int rep=0;
  for(int i=0; i<listsize; i++) {
    string type=parselist[i][0];
    string tag=parselist[i][1];
    if(tag==cond) {
      if(type!="int_array") cout << "type mismatch for " << tag << endl;
      rep++;
      if(rep>1) break;
      for(unsigned int j=2; j<parselist[i].size(); j++) {
	array.push_back(atoi(parselist[i][j].c_str()));
      }      
    }
  }
  error(rep,cond);
  return array;  
}

vector <float> config::getFloatArray(string cond) {
  vector <float> array;
  int rep=0;
  for(int i=0; i<listsize; i++) {
    string type=parselist[i][0];
    string tag=parselist[i][1];
    if(tag==cond) {
      if(type!="float_array") cout << "type mismatch for " << tag << endl;
      rep++;
      if(rep>1) break;
      for(unsigned int j=2; j<parselist[i].size(); j++) {
	array.push_back(atof(parselist[i][j].c_str()));
      }      
    }
  }
  error(rep,cond);
  return array;  
}

vector <double> config::getDoubleArray(string cond) {
  vector <double> array;
  int rep=0;
  for(int i=0; i<listsize; i++) {
    string type=parselist[i][0];
    string tag=parselist[i][1];
    if(tag==cond) {
      if(type!="double_array") cout << "type mismatch for " << tag << endl;
      rep++;
      if(rep>1) break;
      for(unsigned int j=2; j<parselist[i].size(); j++) {
	array.push_back(atof(parselist[i][j].c_str()));
      }      
    }
  }
  error(rep,cond);
  return array;  
}

vector <string> config::getStringArray(string cond) {
  vector <string> array;
  int rep=0;
  for(int i=0; i<listsize; i++) {
    string type=parselist[i][0];
    string tag=parselist[i][1];
    if(tag==cond) {
      if(type!="string_array") cout << "type mismatch for " << tag << endl;
      rep++;
      if(rep>1) break;
      for(unsigned int j=2; j<parselist[i].size(); j++) {
        string astring = parselist[i].at(j);
        
        while(astring.find('"')!=string::npos) {
          astring.replace(astring.find('"'), 1, 1, ' ');
        }

        array.push_back(astring);
      }
    }
  }
  error(rep,cond);
  return array;
}

vector <vector <int> > config::getIntMatrix(string cond) {
  vector <int> array;
  vector <vector <int> > matrix;
  int rep=0;
  int pos=0;
  for(int i=0; i<listsize; i++) {
    string type=parselist[i][0];
    string tag=parselist[i][1];
    if(tag==cond) {
      if(type!="int_matrix") cout << "type mismatch for " << tag << endl;
      if(rep<2) {
	for(unsigned int j=2; j<parselist[i].size(); j++) {
	  array.push_back(atoi(parselist[i][j].c_str()));
	}
	matrix.push_back(array);
      }
      pos=i+1;
      rep++;
    }
  }
  array.clear();
  for(int i=pos; i<listsize; i++) {
    string type=parselist[i][0];
    if(type!="\\") break;
    for(unsigned int j=1; j<parselist[i].size(); j++) {
      array.push_back(atoi(parselist[i][j].c_str()));
    }      
    matrix.push_back(array);
    array.clear();
  }
  error(rep,cond);
  return matrix;  
}

void config::error(int rep, string cond) {
  if(rep>1)
    cout << "There is more the one occurrence of \"" << cond << "\" in the config file!" << endl;
  if(rep==0) 
    cout << "I could not find " << cond << " in config file." << endl;
}



#include <assert.h>
#include <stdio.h>
#include "PParameterReader.h"

PParameterReader::PParameterReader(const char *file)
{
  _env = new TEnv(file);

  return;
}


PParameterReader::~PParameterReader()
{
  delete _env;
}


void PParameterReader::_Assert(const char *name)
{
  if( ! _env->Defined(name) ){
    fprintf(stderr,"ERROR: Could not locate parameter with name '%s'\n",name);
    assert( _env->Defined(name) ); 
  }
}

std::string PParameterReader::GetString(const char* name)
{
  this->_Assert(name);
  return this->GetString(name, "(1)");
}


std::string PParameterReader::GetString(const char *name, const std::string& def)
{
  this->_Assert(name);
  return _env->GetValue(name, def.c_str());
}


Float_t PParameterReader::GetFloat(const char *name)
{
  this->_Assert(name);
  return this->GetFloat(name,0.0);
}

Double_t PParameterReader::GetDouble(const char *name)
{
  this->_Assert(name);
  return this->GetDouble(name,0.0);
}

Int_t PParameterReader::GetInt(const char *name)
{
  this->_Assert(name);
  return this->GetInt(name,0);
}

const char* PParameterReader::GetChar(const char *name)
{
  this->_Assert(name);
  return this->GetChar(name,0);
}

Bool_t PParameterReader::GetBool(const char *name)
{
  this->_Assert(name);
  return this->GetBool(name,kTRUE);
}


Float_t PParameterReader::GetFloat(const char *name, const Float_t dflt)
{
  this->_Assert(name);
  Double_t val = _env->GetValue(name,(Double_t)dflt);
  return (Float_t)val;
}

Double_t PParameterReader::GetDouble(const char *name, 
					   const Double_t dflt)
{
  this->_Assert(name);
  return _env->GetValue(name,dflt);
}

Int_t PParameterReader::GetInt(const char *name, const Int_t dflt)
{
  this->_Assert(name);
  return _env->GetValue(name,dflt);
}

const char* PParameterReader::GetChar(const char *name,
					     const char* dflt)
{
  this->_Assert(name);
  return _env->GetValue(name,dflt);
}

Bool_t PParameterReader::GetBool(const char *name, const Bool_t dflt)
{
  this->_Assert(name);
  Bool_t b = _env->GetValue(name,dflt);

  if(b){
    return kTRUE;
  }else{
    return kFALSE;
  }
}

std::vector<float> PParameterReader::GetVFloat(const char *name, const std::string& delim)
{ 
  TString s(this->GetString(name, "").c_str());
  std::vector<float> result;
  TObjArray *tokens = s.Tokenize(delim.c_str());
  TIter iter(tokens);
  while(TObject *p = iter.Next()) {
    TObjString *item = (TObjString*)p;
    result.push_back(atof(item->GetString().Data()));
  }
  delete tokens;
  return result;
}

std::vector<double> PParameterReader::GetVDouble(const char *name, const std::string& delim)
{ 
  TString s(GetString(name, "").c_str());
  std::vector<double> result;
  TObjArray *tokens = s.Tokenize(delim.c_str());
  TIter iter(tokens);
  while(TObject *p = iter.Next()) {
    TObjString *item = (TObjString*)p;
    result.push_back(atof(item->GetString().Data()));
  }
  delete tokens;
  return result;
}

std::vector<std::string> PParameterReader::GetVString(const char *name, const std::string& delim)
{
  TString s(GetString(name, "").c_str());
  std::vector<std::string> result;
  TObjArray *tokens = s.Tokenize(delim.c_str());
  TIter iter(tokens);
  while(TObject *p = iter.Next()) {
    TObjString *item = (TObjString*)p;
    result.push_back(item->GetString().Data());
  }
  delete tokens;
  return result;
}

std::vector<int> PParameterReader::GetVInt(const char *name, const std::string& delim)
{
  TString s(GetString(name, "").c_str());
  std::vector<int> result;
  TObjArray *tokens = s.Tokenize(delim.c_str());
  TIter iter(tokens);
  while(TObject *p = iter.Next()) {
    TObjString *item = (TObjString*)p;
    result.push_back(strtol(item->GetString().Data(),0,0));
  }
  delete tokens;
  return result;
}

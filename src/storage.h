#pragma once

#include <unordered_map>
#include <fstream>
#include <sstream>

#include "define.h"
#include "logical_function.h"

using namespace phil;

template <class T> class storage_t {
 private:
  std::unordered_map<std::string, T> m_storage;

 public:
  storage_t(const std::string &fn) {
    std::ifstream ifsInput(fn.c_str());
    std::string   key, val;
      
    while(ifsInput.good()) {
      getline(ifsInput, key, '\t');
      getline(ifsInput, val);
      
      m_storage[key] = T(val);
    }
  }
  
  const std::unordered_map<std::string, T> &storage() const { return m_storage; };
};

class sparse_vector_storage_t {
 private:
  std::unordered_map<int, float> m_sv;
  
 public:
  sparse_vector_storage_t() {};
  
  sparse_vector_storage_t(const std::string &fv) {    
    std::istringstream iss(fv);
    std::string key, val;
    
    while(iss.good()) {
      getline(iss, key, ':');
      getline(iss, val, ' ');
      m_sv[atoi(key.c_str())] = atof(val.c_str());
    }
  };

  const std::unordered_map<int, float> vector() const { return m_sv; };
  
  std::string toString() const {
    std::string ret;
    for(auto it: m_sv) { ret += phil::format("%d:%.2f ", it.first, it.second); }
    return ret;
  }
};

class logical_function_storage_t {
 private:
  std::vector<lf::logical_function_t> m_lfs;
  
 public:
  logical_function_storage_t() {};
  
  logical_function_storage_t(const std::string &sexp) {
    std::istringstream iss(sexp);
    for(sexp::reader_t sr(iss); !sr.is_end(); sr.read()) {
      if(sr.get_stack()->is_functor("label"))
        for(int i=1; i< sr.get_stack()->children.size(); i++)
          m_lfs.push_back(lf::logical_function_t(*sr.get_stack()->children[i]));
    };
  }

  const std::vector<lf::logical_function_t> &lfs() const { return m_lfs; };
};

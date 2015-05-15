#pragma once

#include "phillip.h"

using namespace phil;

namespace util {
  typedef hash_map<int, float> sparse_vector_t;

  void getEdgesExplaining(const pg::proof_graph_t *graph, const pg::node_t &node, hash_set<pg::edge_idx_t> *pOut);
  bool isExplained(const pg::proof_graph_t *graph, const ilp::ilp_problem_t *pILPP, const ilp::ilp_solution_t &sol, const pg::node_t &node);
  void solutionToLiterals(const pg::proof_graph_t *graph, const ilp::ilp_problem_t *pILPP, const ilp::ilp_solution_t &sol, std::vector<literal_t> *pOut);
  bool doesSolutionContains(const pg::proof_graph_t *graph, const std::vector<literal_t> &lfSol, const lf::logical_function_t &haystack);
  ilp::variable_idx_t createConditionedIndicator(ilp::ilp_problem_t *prob, std::vector<ilp::variable_idx_t> vVars, std::vector<int> vSigns, bool fAnd, const std::string &name = "");
  void getMatchingLiterals(const std::vector<literal_t> &lfTemplate, const std::vector<literal_t> &lfSearchSpace, std::vector<std::vector<literal_t> > *pOut1, std::vector<hash_map<term_t, hash_set<term_t> > > *pOut2);
  void lfToSetOfLiterals(const lf::logical_function_t &haystack, std::vector<literal_t> *pOut);
  void extractRelatedPotentialNodes(const pg::proof_graph_t *graph, const std::vector<literal_t> &lfTemplate, std::vector<pg::node_idx_t> *pOut);

  inline int getIndex(const std::vector<literal_t> &lsLHS, const literal_t &l) {
    for(int i=0; i<lsLHS.size(); i++)
      if(lsLHS[i].predicate == l.predicate && lsLHS[i].terms.size() == l.terms.size())
        return i;
  }
    
  
  inline void forceILPvarval(ilp::ilp_problem_t *prob, ilp::variable_idx_t var, float val) {
    ilp::constraint_t con("force", ilp::OPR_EQUAL, val);
    con.add_term(var, 1.0);
    prob->add_constraint(con);
  }
    
  inline std::string getObsShortName(const std::string &name) { return name.substr(name.find("::")+2); }
  inline std::string literalsToString(const std::vector<literal_t> &literals) {
    std::string ret;
    
    for(auto l: literals) {
      ret += l.to_string() + ",";
    }

    return ret;
  }

  typedef std::pair<term_t, term_t> eq_t;
  
  template <class T, typename F> void enumerateMatchingLiteralsAndEqs(const std::vector<literal_t> &lfTemplate, const std::vector<T> &lfSearchSpace, F FgetLiteral,
                                                                      std::vector<std::vector<int> > *pOut, std::vector<std::vector<std::vector<eq_t> > > *pOutEqs) {
    std::vector<std::vector<literal_t> > out;
    
    const std::string myEqPred = "__EQ";
    
    for(auto lf: lfTemplate) {
      // Obtain set of literals with the specified predicate/arity.
      pOut->push_back(std::vector<int>());
      pOutEqs->push_back(std::vector<std::vector<eq_t> >());

      for(int i=0; i<lfSearchSpace.size(); i++) {
        const literal_t &lfss = FgetLiteral(lfSearchSpace[i]);
        
        if((lfss.predicate == lf.predicate || (lf.predicate[0] == '!' && lfss.predicate == lf.predicate.substr(1)))
           && lfss.terms.size() == lf.terms.size()) {
          std::vector<eq_t> eqs;

          // Add unifiers.
          for(int i=0; i<lfss.terms.size(); i++)
            eqs.push_back(eq_t(lf.terms[i], lfss.terms[i]));
          
          (*pOut)[pOut->size()-1].push_back(i);
          (*pOutEqs)[pOutEqs->size()-1].push_back(eqs);
        }
      }
    }
  }

  template <class T, typename F> bool combination(std::vector<std::vector<T> > &data, F callback) {
    for(std::vector<int> indices(data.size()); indices[indices.size()-1] < data[indices.size()-1].size(); ) {

      // Hey, work on it!
      if(!callback(indices)) return false;
      
      // Go to the next permutation.
      indices[0]++;

      for(int i=0; i<indices.size()-1 && indices[i] >= data[i].size(); i++) {
        indices[i+1]++;
        indices[i] = 0;
      }
    }

    return true;
  }
  
}

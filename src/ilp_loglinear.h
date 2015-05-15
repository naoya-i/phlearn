#pragma once

#include "phillip.h"
#include "storage.h"

#include "util.h"

namespace phil
{
 
  /** A namespace about factories of linear-programming-problems. */
  namespace ilp
  {
    ilp::variable_idx_t _getNodeVar(ilp::ilp_problem_t *prob, pg::node_idx_t ni);
    ilp::variable_idx_t _getHypernodeVar(ilp::ilp_problem_t *prob, pg::hypernode_idx_t ni);
      
    class loglinear_converter_t : public ilp_converter_t
    {
    private:
      const storage_t<sparse_vector_storage_t> &m_stFv;
      const util::sparse_vector_t& m_wv;
    
      mutable hash_map<ilp::variable_idx_t, util::sparse_vector_t> m_fvMap;
      mutable ilp::ilp_problem_t *m_prob;

      mutable hash_map<ilp::variable_idx_t, float> m_scores;
      mutable hash_set<ilp::variable_idx_t>        m_unifyExplanations;
      
    public:
    
    loglinear_converter_t(phillip_main_t *ptr, const util::sparse_vector_t &wv, const storage_t<sparse_vector_storage_t> &stFv)
      : ilp_converter_t(ptr), m_stFv(stFv), m_wv(wv) {}
      virtual ilp_converter_t* duplicate(phillip_main_t *ptr) const;
      virtual ilp::ilp_problem_t* execute() const;
      virtual bool is_available(std::list<std::string>*) const;
      virtual std::string repr() const;

      void reset() const;
      inline ilp::ilp_problem_t *getILPProblem() { return m_prob; };

      const hash_map<ilp::variable_idx_t, float> &scores() const { return m_scores; };
      const hash_map<ilp::variable_idx_t, util::sparse_vector_t> &fvMap() const { return m_fvMap; };
      void getSolutionFeatureVector(ilp::ilp_solution_t &sol, util::sparse_vector_t *pOut) const;
    };

    class loglinear_xml_decorator_t : public ilp::solution_xml_decorator_t
    {
    public:
      loglinear_xml_decorator_t(const hash_map<pg::node_idx_t, ilp::variable_idx_t> &node2asvar,
                                const hash_map<pg::node_idx_t, ilp::variable_idx_t> &node2acvar,
                                const hash_map<ilp::variable_idx_t, util::sparse_vector_t> &fvMap)
        : m_node2asvar(node2asvar), m_node2acvar(node2acvar), m_fvMap(fvMap) {}

      virtual void get_literal_attributes(const ilp_solution_t *sol, pg::node_idx_t idx,
                                          hash_map<std::string, std::string> *out) const;
    private:
      const hash_map<pg::node_idx_t, ilp::variable_idx_t> m_node2asvar, m_node2acvar;
      const hash_map<ilp::variable_idx_t, util::sparse_vector_t>          m_fvMap;
    };
    
  }

}



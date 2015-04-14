
#include "ilp_loglinear.h"

#include "util.h"
#include "eqc.h"

namespace phil
{
  std::string _toString(const util::sparse_vector_t &v) {
    std::string ret;
    
    for(auto it: v) ret += format("%d:%.2f ", it.first, it.second);

    return ret;
  }
  
  namespace ilp {
    inline void _addFeatureVector(util::sparse_vector_t *pOut, const util::sparse_vector_t &sv, float fValue = -9999.0f) {
      for(auto fv: sv)
        (*pOut)[fv.first] += (fValue != -9999.0f ? fValue : fv.second);
    }

    inline void _injectSV(const storage_t<sparse_vector_storage_t> &stFv, const std::string &name, util::sparse_vector_t *pOut, float fValue = -9999.0f) {
      if(stFv.storage().end() != stFv.storage().find(name)) {
        print_console("  FV:     " + stFv.storage().at(name).toString());
        _addFeatureVector(pOut, stFv.storage().at(name).vector(), fValue);
      }
    }

    void _createFeatureVectorOfObservation(const pg::proof_graph_t* graph, const storage_t<sparse_vector_storage_t> &stFv, util::sparse_vector_t *pOut, const pg::node_t &n) {
      print_console("  --");
      print_console("  OBS:    " + n.to_string());
      print_console("  OBS:    " + n.arity());

      _injectSV(stFv, n.arity(), pOut);
    }

    void _createFeatureVectorOfUnification(const pg::proof_graph_t* graph, const storage_t<sparse_vector_storage_t> &stFv, util::sparse_vector_t *pOut, const pg::node_t &nSmaller, const pg::node_t &nLarger) {
      print_console("  --");
      print_console("  UNI:    " + nSmaller.arity());

      _injectSV(stFv, "U:" + nSmaller.arity(), pOut);
    }
    
    void _createFeatureVectorOfAxiom(const pg::proof_graph_t* graph, const storage_t<sparse_vector_storage_t> &stFv, util::sparse_vector_t *pOut, const lf::axiom_t &axiom, const pg::node_t &n) {
      print_console("  --");
      print_console("  SELECT: " + n.to_string());
      print_console("  AXIOM:  " + axiom.name);
      print_console("  LF:     " + axiom.func.to_string());

      // Add per-axiom feature vector.
      _injectSV(stFv, axiom.name, pOut);

      // Add per-predicate feature vector.
      std::vector<literal_t> lsLHS;
      
      assert(axiom.func.is_operator(lf::OPR_IMPLICATION) || axiom.func.is_operator(lf::OPR_PARAPHRASE));
      
      if(axiom.func.branch(0).is_operator(lf::OPR_LITERAL))
        lsLHS.push_back(axiom.func.branch(0).literal());
      
      else {
        for(auto &lf: axiom.func.branch(0).branches()) {
          assert(lf.is_operator(lf::OPR_LITERAL));
          lsLHS.push_back(lf.literal());
        }
      }

      // Identify the index of literal in LHS.
      int idxAxiom = util::getIndex(lsLHS, n.literal());
      _injectSV(stFv, format("%s,%d", axiom.name.c_str(), idxAxiom), pOut);
      print_console(format("  INDEX: %d", idxAxiom));
    }

    float _getImportanceOfLiteralInAxiom(const lf::axiom_t &axiom, const pg::node_t &n) {
      
      // Add per-predicate feature vector.
      std::vector<literal_t> lsLHS;
      
      assert(axiom.func.is_operator(lf::OPR_IMPLICATION) || axiom.func.is_operator(lf::OPR_PARAPHRASE));
      
      if(axiom.func.branch(0).is_operator(lf::OPR_LITERAL))
        lsLHS.push_back(axiom.func.branch(0).literal());
      
      else {
        for(auto &lf: axiom.func.branch(0).branches()) {
          assert(lf.is_operator(lf::OPR_LITERAL));
          lsLHS.push_back(lf.literal());
        }
      }

      return 1.0f / (float)lsLHS.size();
    }
    
    void _createFeatureVector(const pg::proof_graph_t* graph, const storage_t<sparse_vector_storage_t> &stFv, util::sparse_vector_t *pOut, const pg::node_t &n) {
      kb::knowledge_base_t *pKB = kb::knowledge_base_t::instance();

      if(n.is_equality_node() || n.is_non_equality_node() ||
         n.is_transitive_equality_node() ) {
        //_createFeatureVectorOfEq(graph, stFv, pOut, n);

      } else if(pg::NODE_OBSERVABLE == n.type()) {
        _createFeatureVectorOfObservation(graph, stFv, pOut, n);

        _injectSV(stFv, "PER-PRED_IMP", pOut, -1.0f);
        
      } else {
        // Create features of axioms used for n.
        float                       fCounter = -1.0f;
        std::vector<pg::node_idx_t> nstack;
        nstack.push_back(n.index());

        while(nstack.size() > 0) {
          pg::node_idx_t      sn = nstack.back(); nstack.pop_back();
          pg::hypernode_idx_t hn = graph->node(sn).master_hypernode();
          const hash_set<pg::edge_idx_t> *pEdges;

          if(NULL == (pEdges = graph->search_edges_with_hypernode(hn))) continue;

          // I have a favor for nodes containing hn as a head.
          for(auto e: *pEdges) {
            if(-1 == graph->edge(e).head()) continue;
            
            const std::vector<pg::node_idx_t> &nodes = graph->hypernode(graph->edge(e).head());

            // Push the feature vector into the buffer.
            if(nodes.end() != std::find(nodes.begin(), nodes.end(), sn)) {
              float fLocalImp = _getImportanceOfLiteralInAxiom(pKB->get_axiom(graph->edge(e).axiom_id()), graph->node(sn));
              fCounter *= 1.1f*fLocalImp;
              _createFeatureVectorOfAxiom(graph, stFv, pOut, pKB->get_axiom(graph->edge(e).axiom_id()), graph->node(sn));

              print_console(format("  importance: %f, (accum.: %f)", fLocalImp, fCounter));
              
              // Go into deeper...
              for(auto snn: graph->hypernode(graph->edge(e).tail()))
                nstack.insert(nstack.begin(), snn);
              
              break;
            }
          }
        }

        _injectSV(stFv, "PER-PRED_IMP", pOut, fCounter);
      
        for(auto e: n.evidences()) {
          if(pg::NODE_OBSERVABLE == graph->node(e).type()) {
            _createFeatureVectorOfObservation(graph, stFv, pOut, graph->node(e));
          }
        }
      }
    }

    float _getScore(const pg::proof_graph_t *graph, const util::sparse_vector_t &wv, const storage_t<sparse_vector_storage_t> &stFv,
                    const pg::node_t &n, util::sparse_vector_t *pOutFV) {
      float               ret = 0.0;

      assert(NULL != pOutFV);
      
      // Create the feature vetor.
      print_console("-- fv: " + n.to_string());
      _createFeatureVector(graph, stFv, pOutFV, n);
        
      // Calculate the weighted score.
      for(auto it: *pOutFV) {
        auto w = wv.find(it.first);
        ret += (wv.end() != w ? w->second : 0.0) * it.second;
      }

      return ret;
    }

    float _getScoreOfUnification(const pg::proof_graph_t *graph, const util::sparse_vector_t &wv, const storage_t<sparse_vector_storage_t> &stFv,
                                 const pg::node_t &nSmaller, const pg::node_t &nLarger, util::sparse_vector_t *pOutFV) {
      float               ret = 0.0;

      assert(NULL != pOutFV);
      
      // Create the feature vetor.
      print_console("-- fv: " + nSmaller.to_string() + "," + nLarger.to_string());
      _createFeatureVectorOfUnification(graph, stFv, pOutFV, nSmaller, nLarger);
        
      // Calculate the weighted score.
      for(auto it: *pOutFV) {
        auto w = wv.find(it.first);
        ret += (wv.end() != w ? w->second : 0.0) * it.second;
      }

      return ret;
    }
    
    ilp::variable_idx_t _getNodeVar(ilp::ilp_problem_t *prob, pg::node_idx_t ni) {
      ilp::variable_idx_t var = prob->find_variable_with_node(ni);
      return -1 == var ? prob->add_variable_of_node(ni, 0.0) : var;
    }

    ilp::variable_idx_t _getHypernodeVar(ilp::ilp_problem_t *prob, pg::hypernode_idx_t hi) {
      ilp::variable_idx_t var = prob->find_variable_with_hypernode(hi);
      return -1 == var ? prob->add_variable_of_hypernode(hi, 0.0, false) : var;
    }
    
    ilp_converter_t* loglinear_converter_t::duplicate(phillip_main_t *ptr) const {
      return new loglinear_converter_t(ptr, m_wv, m_stFv);
    }

    void loglinear_converter_t::reset() const {
      m_prob           = NULL;
      m_fvMap.clear();
      m_scores.clear();
    }

    void _equalizeEdgeAndHN(ilp::ilp_problem_t *prob, const pg::proof_graph_t* graph, pg::edge_idx_t ei) {
      if(-1 == graph->edge(ei).head()) {
        ilp::constraint_t con("", ilp::OPR_LESS_EQ, 0);
        con.add_term(_getHypernodeVar(prob, graph->edge(ei).tail()), 1.0);
        con.add_term(prob->find_variable_with_edge(ei), -2.0);
        prob->add_constraint(con);
        
      } else {
        ilp::constraint_t con("", ilp::OPR_LESS_EQ, 1);
        con.add_term(_getHypernodeVar(prob, graph->edge(ei).tail()), 1.0);
        con.add_term(_getHypernodeVar(prob, graph->edge(ei).head()), 1.0);
        con.add_term(prob->find_variable_with_edge(ei), -2.0);
        prob->add_constraint(con);
        
      }
    }
    
    void _doubleImplication(ilp::ilp_problem_t *prob, const pg::proof_graph_t* graph, pg::hypernode_idx_t hi) {
      if(-1 == hi) return;
        
      const std::vector<pg::node_idx_t>& hn = graph->hypernode(hi);
      ilp::constraint_t
        con1("hey1", ilp::OPR_LESS_EQ, 0),
        con2("hey2", ilp::OPR_LESS_EQ, hn.size()-1);
      std::string conname;

      conname = prob->variable(_getHypernodeVar(prob, hi)).name() + " <=> ";
      
      con1.add_term(_getHypernodeVar(prob, hi), hn.size());
      con2.add_term(_getHypernodeVar(prob, hi), -1.0);
                    
      for(auto ni: hn) {
        conname += prob->variable(_getNodeVar(prob, ni)).name() + (hn.back() != ni ? " and " : "");
        con1.add_term(_getNodeVar(prob, ni), -1.0);
        con2.add_term(_getNodeVar(prob, ni), 1.0);
      }

      con1.set_name(conname);
      con2.set_name(conname);
      
      prob->add_constraint(con1);
      prob->add_constraint(con2);
    }
    
    ilp::ilp_problem_t* loglinear_converter_t::execute() const {
      reset();

      const pg::proof_graph_t* graph = phillip()->get_latent_hypotheses_set();
      m_prob = new ilp::ilp_problem_t(graph, new ilp::basic_solution_interpreter_t(), true);
      m_prob->disable_economization();

      convert_proof_graph(m_prob);
      if(m_prob->is_timeout()) return m_prob;
      
      hash_map<pg::node_idx_t, ilp::variable_idx_t> node2var;
      hash_map<pg::node_idx_t, ilp::variable_idx_t> uni2var;

      for(auto i=0; i<graph->edges().size(); i++) {
        _doubleImplication(m_prob, graph, graph->edge(i).head());
        _doubleImplication(m_prob, graph, graph->edge(i).tail());
        _equalizeEdgeAndHN(m_prob, graph, i);
      }

      for(auto node: graph->nodes()) {
        // Force this node to be hypothesized if it is observable.
        if(node.type() == pg::NODE_OBSERVABLE) {
          ilp::constraint_t con("",  ilp::OPR_EQUAL, 1);
          con.add_term(_getNodeVar(m_prob, node.index()), 1.0);
          m_prob->add_constraint(con);
        }
        
        // Create an ILP variable F such that:
        //   FV=1 <=> varNode=1 ^ varHN_1=0 ^ varHN_2=0 ^ ...varHN_N=0
        std::vector<ilp::variable_idx_t> vFeatureConditionVars;
        std::vector<int>                 vFeatureConditionSigns;

        vFeatureConditionSigns.push_back(1);
        vFeatureConditionVars.push_back(_getNodeVar(m_prob, node.index()));

        hash_set<pg::edge_idx_t> edges =
          NULL == graph->search_hypernodes_with_node(node.index()) ?
          hash_set<pg::edge_idx_t>() :
          graph->enumerate_edges_with_node(node.index());
        
        for(auto ei: edges) {
          const pg::edge_t                  &e              = graph->edge(ei);
          const std::vector<pg::node_idx_t> &hypernodesTail = graph->hypernode(e.tail());
          
          if(hypernodesTail.end() == std::find(hypernodesTail.begin(), hypernodesTail.end(), node.index())) continue;

          if(e.is_unify_edge()) {
            // e.g., {x=y, u=v} => {p(x,u), p(y,v)}
            assert(hypernodesTail.size() == 2);

            // Make {x=y, u=v, p(x,u)} or {x=y, u=v, p(y,v)} an explanation.
            // Smaller node will be explained by larger nodes.
            util::sparse_vector_t            svT1, svT2;
            pg::node_idx_t                   largerNode =
              _getScore(graph, m_wv, m_stFv, graph->node(hypernodesTail[0]), &svT1) >
              _getScore(graph, m_wv, m_stFv, graph->node(hypernodesTail[1]), &svT2) ?
              hypernodesTail[0] : hypernodesTail[1];
            std::vector<ilp::variable_idx_t> vConditionVars;
            std::vector<int>                 vConditionSigns;
            std::string                      e_str;

            if(largerNode != node.index()) {
              vConditionSigns.push_back(1);
              vConditionVars.push_back(_getNodeVar(m_prob, largerNode));
              e_str += graph->node(largerNode).to_string();
              
              vConditionSigns.push_back(1);
              vConditionVars.push_back(_getNodeVar(m_prob, node.index()));
              e_str += "=>" + graph->node(node.index()).to_string();

              // Add equalities.
              if(-1 != e.head()) { // For the case like: {} => {p(u), p(v)}
                for(auto ni: graph->hypernode(e.head())) {
                  vConditionSigns.push_back(1);
                  vConditionVars.push_back(_getNodeVar(m_prob, ni));
                }
              }

              ilp::variable_idx_t varFI = util::createConditionedIndicator(m_prob, vConditionVars, vConditionSigns, true, "UE_" + e_str);
              vFeatureConditionSigns.push_back(-1);
              vFeatureConditionVars.push_back(varFI);

              float score         = _getScoreOfUnification(graph, m_wv, m_stFv, node, graph->node(largerNode), &m_fvMap[varFI]);
              m_prob->variable(varFI).set_coefficient(score);
              m_scores[varFI]     = score;
              uni2var[largerNode] = varFI;
              
            }
            
          } else {
            vFeatureConditionSigns.push_back(-1);
            vFeatureConditionVars.push_back(_getHypernodeVar(m_prob, e.head()));
            
          }
        }

        // Set the score.
        ilp::variable_idx_t varFI = util::createConditionedIndicator(m_prob, vFeatureConditionVars, vFeatureConditionSigns, true, "F");
        float               score = _getScore(graph, m_wv, m_stFv, node, &m_fvMap[varFI]);
        m_prob->variable(varFI).set_coefficient(score);
        m_scores[varFI] = score;
        node2var[node.index()] = varFI;
      }

      m_prob->add_xml_decorator(new loglinear_xml_decorator_t(node2var, m_fvMap));
      m_prob->add_attributes("converter", "loglinear");
        
      return m_prob;
    }


    bool loglinear_converter_t::is_available(std::list<std::string> *messages) const {
      return true;
    }


    std::string loglinear_converter_t::repr() const {
      return "Log-Linear-Factory";
    }

    void loglinear_converter_t::getSolutionFeatureVector(ilp::ilp_solution_t &sol, util::sparse_vector_t *pOut) const {
      const pg::proof_graph_t* graph = phillip()->get_latent_hypotheses_set();
      const ilp::ilp_problem_t *pILPP = phillip()->get_ilp_problem();

      for(auto it: m_fvMap) {
        if(sol.variable_is_active(it.first)) _addFeatureVector(pOut, it.second);
      }
    }

    void loglinear_xml_decorator_t::get_literal_attributes(
                                                           const ilp_solution_t *sol, pg::node_idx_t idx,
                                                           hash_map<std::string, std::string> *out) const
    {
      auto i = m_node2nodevar.find(idx);
      
      if(i != m_node2nodevar.end()) {
        auto j = m_fvMap.find(i->second);

        if(j != m_fvMap.end()) {
          (*out)["feature_vec"] = _toString(j->second);
        //(*out)["active"] = sol->variable_is_active(nodevar) ? "yes" : "no";
        }
      }
    }
    
  }  
}

    


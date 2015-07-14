
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

    void _createFeatureVectorOfAxiom(const pg::proof_graph_t* graph, const storage_t<sparse_vector_storage_t> &stFv, util::sparse_vector_t *pOut, const lf::axiom_t &axiom, const pg::node_t &n) {
      print_console("  --");
      print_console("  SELECT: " + n.to_string());
      print_console("  AXIOM:  " + axiom.name);
      print_console("  LF:     " + axiom.func.to_string());

      // Add per-axiom feature vector.
      _injectSV(stFv, axiom.name, pOut);

      //
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

    void _createFeatureVector(const pg::proof_graph_t* graph, const storage_t<sparse_vector_storage_t> &stFv, util::sparse_vector_t *pOutAssumed, util::sparse_vector_t *pOutActive, const pg::node_t &n, int sizeOfObservation) {
      kb::knowledge_base_t *pKB               = kb::knowledge_base_t::instance();
      auto                  fIncreasingFactor = 1.2f;

      if(n.is_equality_node() || n.is_non_equality_node() ||
         n.is_transitive_equality_node() ) {
        //_createFeatureVectorOfEq(graph, stFv, pOut, n);

      } else if(pg::NODE_OBSERVABLE == n.type()) {
        _createFeatureVectorOfObservation(graph, stFv, pOutActive, n);

      } else {

        //
        // Create feature vector of axioms hypothesized for n.
        float                       fCounter = -1.0f / sizeOfObservation;
        int                         nVisit   = 0;
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
              if(nVisit++ == 0)
                _createFeatureVectorOfAxiom(graph, stFv, pOutActive, pKB->get_axiom(graph->edge(e).axiom_id()), graph->node(sn));

              float fLocalImp = _getImportanceOfLiteralInAxiom(pKB->get_axiom(graph->edge(e).axiom_id()), graph->node(sn));
              fCounter *= fIncreasingFactor * fLocalImp;

              // print_console(format("  axiom: %s", pKB->get_axiom(graph->edge(e).axiom_id()).func.to_string().c_str()));
              // print_console(format("  importance: %f, (accum.: %f)", fLocalImp, fCounter));

              // Go into deeper...
              for(auto snn: graph->hypernode(graph->edge(e).tail()))
                nstack.insert(nstack.begin(), snn);

              fCounter *= graph->hypernode(graph->edge(e).tail()).size();

              break;
            }
          }
        }

      }
    }

    float _inp(const util::sparse_vector_t &v1, const util::sparse_vector_t &v2) {
      float ret = 0.0f;
      for(auto it: v1) {
        auto w = v2.find(it.first);
        ret += (v2.end() != w ? w->second : 0.0) * it.second;
      }
      return ret;
    }

    void _getScore(const pg::proof_graph_t *graph, const util::sparse_vector_t &wv, const storage_t<sparse_vector_storage_t> &stFv,
                   const pg::node_t &n,
                   util::sparse_vector_t *pOutFVassumed, util::sparse_vector_t *pOutFVactive, float *pOutAssumed, float *pOutActive,
                   int sizeOfObservation) {
      util::sparse_vector_t vas, vac;

      // Create the feature vetor.
      if(!n.is_equality_node() && !n.is_non_equality_node())
        print_console("-- fv: " + n.to_string());

      _createFeatureVector(graph, stFv, &vas, &vac, n, sizeOfObservation);

      // Calculate the weighted score.
      if(NULL != pOutFVassumed) *pOutFVassumed = vas;
      if(NULL != pOutAssumed)   *pOutAssumed   = _inp(vas, wv);
      if(NULL != pOutFVactive)  *pOutFVactive  = vac;
      if(NULL != pOutActive)    *pOutActive    = _inp(vac, wv);
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
      m_prob = NULL;
      m_fvMap.clear();
      m_scores.clear();
    }

    ilp::ilp_problem_t* loglinear_converter_t::execute() const {
      reset();

      const pg::proof_graph_t* graph = phillip()->get_latent_hypotheses_set();
      m_prob = new ilp::ilp_problem_t(graph, new ilp::basic_solution_interpreter_t(), true);
      m_prob->disable_economization();

      convert_proof_graph(m_prob);
      if(m_prob->is_timeout()) return m_prob;

      hash_map<pg::node_idx_t, ilp::variable_idx_t> node2asvar, node2acvar;
      hash_map<pg::node_idx_t, ilp::variable_idx_t> uni2var;
      hash_map<pg::hypernode_idx_t, hash_set<pg::hypernode_idx_t> > hnpairs;

      for(auto i=0; i<graph->edges().size(); i++) {
        pg::hypernode_idx_t
          head = graph->edge(i).head(),
          tail = graph->edge(i).tail();

        _doubleImplication(m_prob, graph, head);
        _doubleImplication(m_prob, graph, tail);
        _equalizeEdgeAndHN(m_prob, graph, i);

        if(!graph->edge(i).is_unify_edge() && -1 != head && -1 != tail)
          hnpairs[tail].insert(head);

      }

      //
      // Make hypernodes having the same parent mutual exclusive. (Explanatory XOR)
      if(phillip()->flag("learn_exp_mex")) {
        for(auto i: hnpairs) {
          if(1 >= i.second.size()) continue;

          ilp::constraint_t con(graph->hypernode2str(i.first) + "<-", ilp::OPR_LESS_EQ, 1);

          print_console("MX: " + graph->hypernode2str(i.first));

          for(auto j: i.second) {
            print_console("MX: <- " + graph->hypernode2str(j));
            con.add_term(_getHypernodeVar(m_prob, j), 1.0);
          }

          m_prob->add_constraint(con);
        }
      }

      int sizeOfObservation = 0;

      struct node_feature_vector_t {
        float fAssumedScore, fActiveScore;
        util::sparse_vector_t vAssumed, vActive;
      };

      hash_map<pg::node_idx_t, node_feature_vector_t> storageNodeFeatureVector;

      // Count the number of observed literals.
      for(auto node: graph->nodes())
        if(node.type() == pg::NODE_OBSERVABLE) sizeOfObservation++;

      // Create the feature vector for each node.
      for(auto node: graph->nodes()) {
        _getScore(graph, m_wv, m_stFv, node,
                  &storageNodeFeatureVector[node.index()].vAssumed,
                  &storageNodeFeatureVector[node.index()].vActive,
                  &storageNodeFeatureVector[node.index()].fAssumedScore,
                  &storageNodeFeatureVector[node.index()].fActiveScore,
                  sizeOfObservation);
      }

      for(auto node: graph->nodes()) {
        // Force this node to be hypothesized if it is observable.
        if(node.type() == pg::NODE_OBSERVABLE) {
          ilp::constraint_t con("obs:" + node.to_string(),  ilp::OPR_EQUAL, 1);
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
            float
              fAssumedScore1 = storageNodeFeatureVector[hypernodesTail[0]].fAssumedScore,
              fAssumedScore2 = storageNodeFeatureVector[hypernodesTail[1]].fAssumedScore;

            pg::node_idx_t                   largerNode = fAssumedScore1 > fAssumedScore2 ? hypernodesTail[0] : hypernodesTail[1],
                                             smallerNode = fAssumedScore1 > fAssumedScore2 ? hypernodesTail[1] : hypernodesTail[0];
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
                vConditionSigns.push_back(1);
                vConditionVars.push_back(_getHypernodeVar(m_prob, e.head()));
              }

              ilp::variable_idx_t varFI = util::createConditionedIndicator(m_prob, vConditionVars, vConditionSigns, true, "UE_" + e_str);
              vFeatureConditionSigns.push_back(-1);
              vFeatureConditionVars.push_back(varFI);

              // float score         = _getScoreOfUnification(graph, m_wv, m_stFv, node, graph->node(largerNode), &m_fvMap[varFI]);
              // m_prob->variable(varFI).set_coefficient(score);
              // m_scores[varFI]     = score;
              // uni2var[largerNode] = varFI;
            }

          } else {
            vFeatureConditionSigns.push_back(-1);
            vFeatureConditionVars.push_back(_getHypernodeVar(m_prob, e.head()));

          }
        }

        // Set the score.
        ilp::variable_idx_t
          varFIas        = util::createConditionedIndicator(m_prob, vFeatureConditionVars, vFeatureConditionSigns, true, "Fas"),
          varFIac        = _getNodeVar(m_prob, node.index());
        float
          fAssumedScore = storageNodeFeatureVector[node.index()].fAssumedScore,
          fActiveScore  = storageNodeFeatureVector[node.index()].fActiveScore;
        m_fvMap[varFIas] = storageNodeFeatureVector[node.index()].vAssumed;
        m_fvMap[varFIac] = storageNodeFeatureVector[node.index()].vActive;

        m_prob->variable(varFIas).set_coefficient(fAssumedScore);
        m_prob->variable(varFIac).set_coefficient(fActiveScore);
        m_scores[varFIas] = fAssumedScore;
        m_scores[varFIac] = fActiveScore;
        node2asvar[node.index()] = varFIas;
        node2acvar[node.index()] = varFIac;
      }

      m_prob->add_xml_decorator(new loglinear_xml_decorator_t(node2asvar, node2acvar, m_fvMap));
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
      auto ias = m_node2asvar.find(idx), iac = m_node2acvar.find(idx);

      if(ias != m_node2asvar.end()) {
        auto j = m_fvMap.find(ias->second);

        if(j != m_fvMap.end()) {
          (*out)["feature_vec"] += _toString(j->second);
        //(*out)["active"] = sol->variable_is_active(nodevar) ? "yes" : "no";
        }
      }

      if(iac != m_node2asvar.end()) {
        auto j = m_fvMap.find(iac->second);

        if(j != m_fvMap.end()) {
          (*out)["feature_vec"] += _toString(j->second);
        //(*out)["active"] = sol->variable_is_active(nodevar) ? "yes" : "no";
        }
      }
    }

  }
}

# include "util.h"

# include "phillip.h"
# include "eqc.h"

using namespace phil;

namespace util {

  void solutionToLiterals(const pg::proof_graph_t *graph, const ilp::ilp_problem_t *pILPP, const ilp::ilp_solution_t &sol, std::vector<literal_t> *pOut) {
    hash_map<term_t, hash_set<term_t> > eqs;
    equivalence_class_t                 eqc;

    // Identify literals and equalities in the solution hypothesis.
    for(auto node: graph->nodes()) {
      if(pILPP->node_is_active(sol, node.index())) {
        if(node.is_equality_node())
          eqc.addEquality(node.literal().terms[0], node.literal().terms[1]);
        else
          pOut->push_back(node.literal());
      }
    }

    eqc.updateClusters();

    // Replace the variables.
    for(auto &l: *pOut) {
      for(int i=0; i<l.terms.size(); i++) {
        l.terms[i] = eqc.classes().end() != eqc.classes().find(l.terms[i]) ? term_t(eqc.classes().at(l.terms[i])) : l.terms[i];
      }
    }
  }

  void lfToSetOfLiterals(const lf::logical_function_t &haystack, std::vector<literal_t> *pOut) {
    assert(haystack.is_operator(lf::OPR_AND));

    // Convert the input to the template.
    for(auto lf: haystack.branches()) {
      literal_t l(lf.literal().predicate);
      
      for(auto t: lf.literal().terms)
        l.terms.push_back(t.is_constant() ? t : term_t("tmpl_" + t.string()));
      
      pOut->push_back(l);
    }
  }
  
  bool doesSolutionContains(const pg::proof_graph_t *graph, const std::vector<literal_t> &lfSol, const lf::logical_function_t &haystack) {
    std::vector<literal_t> lfGold;
    lfToSetOfLiterals(haystack, &lfGold);

    // Assume that lfSol is shrunk.
    
    std::vector<std::vector<int> > ieMatches;
    std::vector<std::vector<std::vector<eq_t> > > ieMatchesEqs;
    enumerateMatchingLiteralsAndEqs(lfGold, lfSol, [](const literal_t &l) { return l; }, &ieMatches, &ieMatchesEqs);

    for(int i=0; i<ieMatches.size(); i++) {
      if(0 == ieMatches[i].size()) return false;
    }
    
    // Check the feasibility of equalities. If we can find at least
    // one feasible hypothesis, it will be good.
    return !combination(ieMatches, [&](const std::vector<int> &indices){
        const std::string                  myEqPred = "__EQ";
        hash_map<term_t,hash_set<term_t> > outEqs;
        
        for(int i=0; i<indices.size(); i++) {
          for(auto eq: ieMatchesEqs[i][indices[i]])
            outEqs[eq.first].insert(eq.second);
        }
        
        for(auto it: outEqs) {
          if(it.second.size() > 1) return true;
          
          if(it.first.is_constant()) {
            if(*it.second.begin() != it.first) return true;
          }
        }

        return false;
      });
  }

  ilp::variable_idx_t createConditionedIndicator(ilp::ilp_problem_t *prob, std::vector<ilp::variable_idx_t> vVars, std::vector<int> vSigns, bool fAnd, const std::string &name) {
    // Create an ILP constraint that enforces:
    //   varSwitch = 1 <=> varLiteral  = 1 & varEq1 = 1 & varEq2 = 1 ... (fAnd = 1)
    //   varSwitch = 1 <=> varLiteral  = 1 v varEq1 = 1 v varEq2 = 1 ... (fAnd = 0)
      
    assert(vVars.size() == vSigns.size());
      
    // Count number of positive signs and negative signs.
    int numPositive = 0, numNegative = 0;

    for(auto i: vSigns) {
      if(i>0) numPositive++; else if(i<0) numNegative++;
    }

    assert(numPositive+numNegative == vSigns.size());

    static int s_sw = 0;
    ilp::variable_idx_t
      varSwitch  = prob->add_variable(ilp::variable_t(format(("sw%d" + name).c_str(), s_sw++), 0.0));

    ilp::constraint_t
      con1("hey1", ilp::OPR_LESS_EQ, numNegative),
      con2("hey2", ilp::OPR_GREATER_EQ, fAnd ? -numPositive + 1 : numNegative);
    std::string conname;
      
    con1.add_term(varSwitch, fAnd ? (numPositive+numNegative) : 1);
    conname += format(("sw%d"+name+"=1 <=> ").c_str(), s_sw-1);
    for(int i=0; i<vVars.size(); i++) {
      con1.add_term(vVars[i], vSigns[i] > 0 ? -1.0 : 1.0);
      conname += prob->variable(vVars[i]).name() + (vSigns[i] > 0 ? "=1" : "=0") + (i < vVars.size()-1 ? (fAnd ? " and " : " or ") : "");
    }

    con2.add_term(varSwitch, fAnd ? 1.0 : (numPositive+numNegative));
    for(int i=0; i<vVars.size(); i++)
      con2.add_term(vVars[i], vSigns[i] > 0 ? -1.0 : 1.0);
    con1.set_name(conname);
    con2.set_name(conname);
      
    prob->add_constraint(con1);
    prob->add_constraint(con2);
        
    return varSwitch;
  }

  void extractRelatedPotentialNodes(const pg::proof_graph_t *graph, const std::vector<literal_t> &lfTemplate, std::vector<pg::node_idx_t> *pOut) {
    for(auto lf: lfTemplate) {
      const hash_set<pg::node_idx_t> *pNodes = graph->search_nodes_with_arity(lf.get_arity());

      if(NULL != pNodes) {
        for(auto ni: *pNodes)
          pOut->push_back(ni);
      }
    }
  }
  
}

// Log-linear abduction

#include <fstream>

#include "./binary.h"
#include "./processor.h"

#include "ilp_loglinear.h"
#include "storage.h"

#include "weight_update.h"

const float fDefaultC = 0.5;
const float fDefaultEta = 0.75;
const float fInitialVar = 100;

using namespace phil;

void _printVector(const util::sparse_vector_t &v, std::ostream *pLog) {
  for(auto it: v)
    (*pLog) << it.first << ":" << it.second << " ";
}

void _printWeightVector(int n, const util::sparse_vector_t &vMean, const util::sparse_vector_t &vVar) {
  for(auto it: vMean) {
    std::cout << (n+1) << "\t" << it.first << "\t" << it.second << std::endl; // << "\t" << vVar.at(it.first) << std::endl;
  }
}

bool _getGoldSetsVars(const pg::proof_graph_t *graph, ilp::ilp_problem_t *prob, const hash_set<pg::hypernode_idx_t> &unihns, const lf::logical_function_t &lfGold,
                      std::vector<ilp::variable_idx_t> *pOut, bool f_exclude_transitiveeq, std::ostream *pLog) {
  assert(lfGold.is_operator(lf::OPR_AND));
  
  // Convert the solution into a set of literals.
  std::vector<literal_t> lsGold;
  util::lfToSetOfLiterals(lfGold, &lsGold);
  
  // Extract related potential literals from the proof graph.
  std::vector<pg::node_idx_t> nsPoSol;
  util::extractRelatedPotentialNodes(graph, lsGold, &nsPoSol);

  // Template matching.
  std::vector<std::vector<int> > matches;
  std::vector<std::vector<std::vector<util::eq_t> > > matchesEqs;

  util::enumerateMatchingLiteralsAndEqs(lsGold, nsPoSol,
                                        [&graph](pg::node_idx_t ni) {
                                          return graph->node(ni).literal();
                                        },
                                        &matches, &matchesEqs);

  for(auto &v: matches) {
    if(0 == v.size()) return false;
  }
  
  std::vector<std::vector<ilp::variable_idx_t> > varGoldNodes;
  int                                            numFound = 0;

  // Work on all possible combinations.
  util::combination(matches, [&](const std::vector<int> &indices){
      hash_map<term_t, hash_set<term_t> > eqMap;
      std::vector<ilp::variable_idx_t>    vars;
      std::vector<pg::node_idx_t>         conditionedNodes;
      
      (*pLog) << "<pattern>" << std::endl;

      //
      // Enumerate set of literals by pattern matching.
      (*pLog) << "<literals>" << std::endl;
      for(int i=0; i<indices.size(); i++) {
        (*pLog) << "<literal>" << graph->node(nsPoSol[matches[i][indices[i]]]).to_string() << "</literal>" << std::endl;
        
        conditionedNodes.push_back(nsPoSol[matches[i][indices[i]]]);
        
        // Create equality mapping (e.g., tmpl_x: {y, z, w, ...}, Child: {x, y, z, ...})
        for(auto eq: matchesEqs[i][indices[i]]) {
          eqMap[eq.first].insert(eq.second);
          
          if(eq.first.is_constant())
            eqMap[eq.first].insert(eq.first);
        }
      }
      (*pLog) << "</literals>" << std::endl;

      //
      // Write variable-matching log.
      bool fPossible = true;
      
      (*pLog) << "<variables>" << std::endl;
      for(auto it: eqMap) {
        (*pLog) << "<match target=\"" << it.first.string() << "\">";
          
        //
        // All the pairwise combination of the variables in it.second
        // must be possible. In the below, check the possibility on
        // proof graph created so far.
        for(auto &it2: it.second) {
          for(auto &it3: it.second) {            
            if(it2.string() >= it3.string()) continue;

            (*pLog) << it2.string() << "=" << it3.string() << ",";
            
            //
            // Check whether the equality conditions are in the proofgraph.
            pg::node_idx_t subnode = graph->find_sub_node(it2, it3);

            if(f_exclude_transitiveeq && -1 != subnode) {
              if(graph->node(subnode).is_transitive_equality_node()) subnode = -1;
            }
            
            conditionedNodes.push_back(subnode);
             
            if(-1 == subnode) {
              (*pLog) << "... impossible., ";
              fPossible = false;
            }
          }
        }
        
        (*pLog) << "</match>" << std::endl;
      }
      (*pLog) << "</variables>" << std::endl;

      //
      // Create the condition.
      if(fPossible) {
        bool fAbort = false;
        
        (*pLog) << "<proofgraph-checking>" << std::endl;
        
        for(auto ni: conditionedNodes) {
          (*pLog) << "<node target=\"" << graph->node(ni).to_string() << "\""
                  << " is_transitive_eq=\"" << graph->node(ni).is_transitive_equality_node() << "\">";
          
          const hash_set<pg::node_idx_t>   *pHypernodes = graph->search_hypernodes_with_node(ni);
          std::vector<ilp::variable_idx_t> varsHN;

          if(NULL == pHypernodes) {
            //
            // If the node is equality and created by a transitivity
            // rule, the hypernode might be NULL.
            if(graph->node(ni).is_equality_node()) {
              if(NULL != pOut)
                varsHN.push_back(ilp::_getNodeVar(prob, ni));
              
            } else {
              (*pLog) << "<no-hypernode-found />";
              fAbort = true;
              
            }
            
          } else {
            for(auto hn: *pHypernodes) {
              if(unihns.end() != unihns.find(hn)) continue;
          
              if(NULL != pOut)
                varsHN.push_back(ilp::_getHypernodeVar(prob, hn));
            }
          }
          
          // Create OR-conditioned indicator that becomes 1 iff at least
          // one hypernode related to the node is hypothesized.
          if(!fAbort && NULL != pOut) {
            vars.push_back(util::createConditionedIndicator(prob, varsHN, std::vector<int>(varsHN.size(), 1), false));
            
            for(auto v: varsHN) {
              (*pLog) << prob->variable(v).name() << " or ";
            }
          }
          
          (*pLog) << "</node>" << std::endl;
        }

        // Push an AND-conditioned indicator that becomes 1 iff all vars
        // are true.
        if(!fAbort && NULL != pOut)
          pOut->push_back(util::createConditionedIndicator(prob, vars, std::vector<int>(vars.size(), 1), true, "cond"));
      
        numFound++;
        (*pLog) << "</proofgraph-checking>" << std::endl;
      }
      
      (*pLog) << "</pattern>" << std::endl;
      
      return true;
    });

  if(0 == numFound)
    (*pLog) << "<no-match />" << std::endl;
  
  return numFound > 0;
}

class app_t : public phillip_main_t {

public:
  void testing(const storage_t<sparse_vector_storage_t> &stFv, const storage_t<logical_function_storage_t> &stLabel,
               ilp::loglinear_converter_t *pLLConv,
               const std::vector<lf::input_t>& parsed_inputs, int idx, util::sparse_vector_t &vecMean) {
    std::ostream *pLog = &std::cout;
    int           K = param_int("kbest", 1);

    // Perform inference.
    reset_for_inference();
    set_input(parsed_inputs[idx]);
    execute_enumerator();
    execute_convertor();

    for(int k=0; k<K; k++) {
      std::vector<ilp::ilp_solution_t> sols;
      ilp_solver()->execute(&sols);
    
      auto sol = sols[0];
      sol.print_graph();

      (*pLog) << "<result k=\""<< (1+k) <<"\">" << std::endl;

      // Check the performance.
      const pg::proof_graph_t *graph = get_latent_hypotheses_set();
      hash_set<pg::hypernode_idx_t>    unihns;

      // Retrieve set of edges representing unification.
      for(auto e: graph->edges()) {
        if(e.is_unify_edge()) unihns.insert(e.tail());
      }
    
      const std::vector<lf::logical_function_t> &lfsGold = stLabel.storage().at(util::getObsShortName(parsed_inputs.at(idx).name)).lfs();
    
      for(auto &lfGold: lfsGold)
        (*pLog) << "<label>" << lfGold.to_string() << "</label>" << std::endl;
    
      // Get the feature vector of our best hypothesis.
      ilp::ilp_problem_t      *prob  = ((ilp::loglinear_converter_t*)ilp_convertor())->getILPProblem();
    
      util::sparse_vector_t  vGold;
      std::vector<literal_t> lfSolGold;
  
      pLLConv->getSolutionFeatureVector(sol, &vGold);
      util::solutionToLiterals(graph, prob, sol, &lfSolGold);

      (*pLog) << "<logical-form>" + util::literalsToString(lfSolGold) << "</logical-form>" << std::endl;
      (*pLog) << "<vector>"; _printVector(vGold, &std::cout); (*pLog) << "</vector>" << std::endl;
  
      // If the best hypothesis = yi, then we want the best hypothesis not
      // including yi. Otherwise (the best hypothesis != yi), then we want
      // the best hypothesis including yi.
      int numCorrects = 0;

      (*pLog) <<  "<label-matching-log>" << std::endl;
      
      for(auto &lfGold: lfsGold) {
        if(util::doesSolutionContains(get_latent_hypotheses_set(),
                                      lfSolGold, lfGold)) numCorrects++;
      }

      (*pLog) << "</label-matching-log>" << std::endl;
      (*pLog) <<  "<correct>" << numCorrects << "</correct>" << std::endl;

      if(k+1 != K) {
        (*pLog) <<  "<constraint-of-next-solution>" << std::endl;
        
        // Impose constraint.
        for(auto node: graph->nodes()) {
          if(node.literal().predicate.substr(0, 5) == "will-" && prob->node_is_active(sol, node.index()) && node.type() != pg::NODE_OBSERVABLE) {
            (*pLog) << "<prohibit>" << node.to_string() << "</prohibit>" << std::endl;
            util::forceILPvarval(prob, ilp::_getNodeVar(prob, node.index()), 0);
          }
        }
        
        (*pLog) <<  "</constraint-of-next-solution>" << std::endl;
      }
      
      (*pLog) << "</result>" << std::endl;
    }
  }

  void learn(const storage_t<sparse_vector_storage_t> &stFv, const storage_t<std::string> &stFt,
             const storage_t<logical_function_storage_t> &stLabel,
            ilp::loglinear_converter_t *pLLConv,
            const std::vector<lf::input_t>& parsed_inputs, int idx, util::sparse_vector_t &vecMean, util::sparse_vector_t &vecVariance, std::ostream *pLog) {

    reset_for_inference();
    set_input(parsed_inputs[idx]);

    execute_enumerator();
  
    // 
    const pg::proof_graph_t          *graph  = get_latent_hypotheses_set();
    hash_set<pg::hypernode_idx_t>    unihns;

    // Retrieve set of edges representing unification.
    for(auto e: graph->edges()) {
      if(e.is_unify_edge()) unihns.insert(e.tail());
    }
  
    // Obtain the ground truth of this observation.
    if(stLabel.storage().end() == stLabel.storage().find(util::getObsShortName(parsed_inputs.at(idx).name))) {
      (*pLog) << "<label-status>no-annotation</label-status>" << std::endl;
      return;
    }

    const std::vector<lf::logical_function_t> &lfsGold = stLabel.storage().at(util::getObsShortName(parsed_inputs.at(idx).name)).lfs();
    int                                        numHits = 0;
    std::ostringstream                         trash;
    
    for(auto &lfGold: lfsGold) {
      (*pLog) << "<label>" << lfGold.to_string() << "</label>" << std::endl;
      if(_getGoldSetsVars(graph, NULL, unihns, lfGold, NULL, flag("learn_exclude_transieq"), &trash)) numHits++;
    }

    if(0 == numHits) {
      (*pLog) << "<label-status>no-potential-gold-literal</label-status>" << std::endl;
      return;
    }

    (*pLog) << "<label-status>potential-gold-literal-found</label-status>" << std::endl;

    if(flag("learn_label_check_only"))
      return;
    
    // Perform inference.
    execute_convertor();
    execute_solver();
    
    (*pLog) << "<current-prediction>" << std::endl;
    ilp::ilp_problem_t               *prob   = ((ilp::loglinear_converter_t*)ilp_convertor())->getILPProblem();
  
    auto sol = get_solutions()[0];
    
    if(flag("learn_print_ilp"))
      prob->print(pLog);
    sol.print_graph(pLog);

    // Get the feature vector of our best hypothesis.
    util::sparse_vector_t  vGold;
    std::vector<literal_t> lfSolGold;
  
    pLLConv->getSolutionFeatureVector(sol, &vGold);
    util::solutionToLiterals(graph, prob, sol, &lfSolGold);

    (*pLog) << "<logical-form>" + util::literalsToString(lfSolGold) << "</logical-form>" << std::endl;
    (*pLog) << "<vector>"; _printVector(vGold, pLog); (*pLog) << "</vector>" << std::endl;
  
    // If the best hypothesis = yi, then we want the best hypothesis not
    // including yi. Otherwise (the best hypothesis != yi), then we want
    // the best hypothesis including yi.
    int numCorrects = 0;

    for(auto &lfGold: lfsGold)
      if(util::doesSolutionContains(get_latent_hypotheses_set(),
                                    lfSolGold, lfGold)) numCorrects++;

    (*pLog) <<  "<result>" << numCorrects << "</result>" << std::endl;
    (*pLog) << "</current-prediction>" << std::endl;

    if("structured_perceptron" == param("learn_algo") && 0 < numCorrects)
      return;
      
    // Create ILP variables expressing inclusion of gold literals.
    (*pLog) << "<latent-variable-completion>" << std::endl;
    
    std::vector<ilp::variable_idx_t> vCondLabelSatisfied;

    for(auto &lfGold: lfsGold) {
      std::vector<ilp::variable_idx_t> vCondGoldSets;
    
      (*pLog) << "<find-label label=\"" << lfGold.to_string() << "\">" << std::endl;
      
      if(_getGoldSetsVars(graph, prob, unihns, lfGold, &vCondGoldSets, flag("learn_exclude_transieq"), pLog))
        vCondLabelSatisfied.push_back(util::createConditionedIndicator(prob, vCondGoldSets, std::vector<int>(vCondGoldSets.size(), 1), false));

      (*pLog) << "</find-label>" << std::endl;
    }

    ilp::variable_idx_t varFI = numCorrects > 0 ?
      util::createConditionedIndicator(prob, vCondLabelSatisfied, std::vector<int>(vCondLabelSatisfied.size(), -1), true) : // Ignore yi.
      util::createConditionedIndicator(prob, vCondLabelSatisfied, std::vector<int>(vCondLabelSatisfied.size(),  1), false);  // Keep watching yi.
    util::forceILPvarval(prob, varFI, 1.0);
    
    //get_ilp_problem()->print(pLog);
  
    // Ok, let us try again with the new constraint.
    std::vector<ilp::ilp_solution_t> sols;
    //pLLConv->adjustScores(-1.0);
    ilp_solver()->execute(&sols);

    if(flag("learn_print_ilp"))
      prob->print(pLog);
    
    sols[0].print_graph(pLog);
    //sols[0].print(pLog);
  
    // Get the feature vector of the competitor.
    util::sparse_vector_t  vCompetitor;
    std::vector<literal_t> lfSolCompetitor;
  
    pLLConv->getSolutionFeatureVector(sols[0], &vCompetitor);
    util::solutionToLiterals(graph, prob, sols[0], &lfSolCompetitor);
  
    (*pLog) << "<logical-form>" + util::literalsToString(lfSolCompetitor) << "</logical-form>" << std::endl;
    (*pLog) << "<vector>"; _printVector(vCompetitor, pLog); (*pLog) << "</vector>" << std::endl;
    (*pLog) << "</latent-variable-completion>" << std::endl;
  
    // Update the weight vector.
    if(0 == numCorrects) std::swap(vCompetitor, vGold);
  
    (*pLog) << "<weight-update>" << std::endl;

    float loss;

    if("exact_scw_" == param("learn_algo")) {
      loss = scw::updateWeightVector(&vecMean, &vecVariance, vGold, vCompetitor,
                                     param_float("learn_C",   fDefaultC),
                                     param_float("learn_eta", fDefaultEta),
                                     param_float("learn_initial_var", fInitialVar),
                                     scw::SCW_I,
                                     stFt,
                                     pLog);
      
    } else if("structured_perceptron" == param("learn_algo")) {
      loss = sp::updateWeightVector(&vecMean, vGold, vCompetitor,
                                    param_float("learn_eta", fDefaultEta),
                                    stFt,
                                    pLog);
      
    }
      
    
    (*pLog) << "<loss>" << loss << "</loss>" << std::endl;
    (*pLog) << "</weight-update>" << std::endl;
  }
};

int main(int argc, char* argv[]) {
  initialize();

  app_t phillip;
  bin::execution_configure_t config;
  std::vector<std::string> inputs;

  //
  print_console("Learnign module for log-linear abduction");
  print_console("  for Phillip ver. " + phillip_main_t::VERSION);
    
  bin::parse_options(argc, argv, (phillip_main_t*)&phillip, &config, &inputs);

  // Force ILP converter to be the log-linear converter.
  if("loglinear" != config.ilp_key) {
    print_error("Not supported: " + config.ilp_key);
    return 0;
  }

  //
  print_console("Loading feature vectors and labels...");

  storage_t<sparse_vector_storage_t>    stFeatureVector(phillip.param("feature_vector"));
  storage_t<std::string>                stFeatureTranslator(phillip.param("feature_translator"));
  storage_t<logical_function_storage_t> stLabel(phillip.param("label"));
  
  print_console("done.");

  //
  util::sparse_vector_t vecMean, vecVariance;
  ilp::loglinear_converter_t *pLLConv = new ilp::loglinear_converter_t((phillip_main_t*)&phillip,
                                                                       vecMean,
                                                                       stFeatureVector);
  phillip.set_ilp_convertor(pLLConv);

  //
  bin::preprocess(config, (phillip_main_t*)&phillip);
  
  std::vector<lf::input_t> parsed_inputs;
  proc::processor_t processor;
  bool do_parallel_inference(phillip.flag("do_parallel_inference"));
  bool do_write_parallel_out(phillip.flag("do_write_parallel_out"));
  bool do_test(phillip.flag("testing"));

  if("" != phillip.param("learn_weight")) {
    print_console("Loading weights from "+ phillip.param("learn_weight") +"...");
    
    std::ifstream ifsWeight(phillip.param("learn_weight").c_str());
    std::string   fk, fm;
    int           n;

    if(do_test) {
      while(ifsWeight >> n >> fk >> fm) {
        if(-1 != phillip.param_int("use_iter", -1) && n != phillip.param_int("use_iter", -1))
          continue;
        
        vecMean[atoi(fk.c_str())] = atof(fm.c_str());
      }
      
    } else {
      while(ifsWeight >> fk >> fm) {
        vecMean[atoi(fk.c_str())] = atof(fm.c_str());
      }
      
    }      
    print_console("done. " + format("Weight vector is now %d dimension.", vecMean.size()));
  }
  
  print_console("Loading observations...");

  processor.add_component(new proc::parse_obs_t(&parsed_inputs));
  processor.process(inputs);

  print_console("Completed to load observations.");
  print_console_fmt("    # of observations: %d", parsed_inputs.size());

  kb::knowledge_base_t::instance()->prepare_query();

  std::ostream *pLog = &std::cout;

  if("" != phillip.param("learn_log")) {
    static std::ofstream ofsLog(phillip.param("learn_log").c_str(), std::ios::out);
    pLog = (std::ostream*)&ofsLog;
  }
  
  (*pLog) << "<phillip-learn>" << std::endl;
  
  if (kb::knowledge_base_t::instance()->is_valid_version()) {
    if(do_test) {      
      phillip.write_header();
      
      for (int i = 0; i < parsed_inputs.size(); ++i) {
        const lf::input_t &ipt = parsed_inputs.at(i);
            
        std::string obs_name = ipt.name;
        if (obs_name.rfind("::") != std::string::npos)
          obs_name = obs_name.substr(obs_name.rfind("::") + 2);
            
        if (phillip.is_target(obs_name) and not phillip.is_excluded(obs_name)) {
          phillip.testing(stFeatureVector, stLabel, pLLConv, parsed_inputs, i, vecMean);
          kb::knowledge_base_t::instance()->clear_distance_cache();
        }
      }

      phillip.write_footer();
      
    } else {
      for(int n=0; n<phillip.param_int("learn_iter", 1); n++) {
        for (int i = 0; i < parsed_inputs.size(); ++i) {
          const lf::input_t &ipt = parsed_inputs.at(i);
            
          std::string obs_name = ipt.name;
          if (obs_name.rfind("::") != std::string::npos)
            obs_name = obs_name.substr(obs_name.rfind("::") + 2);
            
          if (phillip.is_target(obs_name) and not phillip.is_excluded(obs_name)) {
            
            // Step forward!
            (*pLog) << format("<round iteration=\"%d\" observation=\"%s\">", 1+n, obs_name.c_str()) << std::endl;
            phillip.learn(stFeatureVector, stFeatureTranslator, stLabel, pLLConv, parsed_inputs, i, vecMean, vecVariance, pLog);
            (*pLog) << "</round>" << std::endl;
            
            kb::knowledge_base_t::instance()->clear_distance_cache();
          }
        }

        _printWeightVector(n, vecMean, vecVariance);
      }
      
    }
  }

  (*pLog) << "</phillip-learn>" << std::endl;
  
  if(&std::cout != pLog)
    ((std::ofstream*)pLog)->close();  
}

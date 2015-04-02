#include "eqc.h"

using namespace phil;

namespace util {
  
  void equivalence_class_t::updateClusters() {
    hash_set<term_t> visited;

    // Ordinary depth-first search.
    std::function<void(hash_set<term_t>*, const term_t&)> _visitVertex = [&](hash_set<term_t> *pOut, const term_t &v){
      if(visited.end() != visited.find(v)) return;
      else                                 visited.insert(v);

      pOut->insert(v);

      if(m_edges.end() != m_edges.find(v))
        for(auto vDest: m_edges.at(v))
          _visitVertex(pOut, vDest);
    };

    // For each vertex, identify clusters (enumerate connected
    // components on the graph).
    for(auto v: m_vertices) {
      hash_set<term_t> cluster;
      _visitVertex(&cluster, v);

      if(cluster.size() > 0) {
        // Choose the representative.
        term_t reprC, reprV;
        bool   fReprC = false;

        for(auto t: cluster) {
          if(t.is_constant()) { reprC = t; fReprC = true; }
          else                { reprV = t; }
        }

        for(auto t: cluster)
          m_class[t] = fReprC ? reprC : reprV;
      }
    }
  }
  
}

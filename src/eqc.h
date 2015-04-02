#pragma once

#include "phillip.h"

using namespace phil;

namespace util {
  
  class equivalence_class_t {
  private:
    hash_map<term_t, hash_set<term_t> > m_edges;
    hash_set<term_t>                    m_vertices;
    hash_map<term_t, term_t>            m_class;

  public:
    equivalence_class_t() {};

    void updateClusters();

    inline const hash_map<term_t, term_t> &classes() { return m_class; };
    inline const hash_map<term_t, hash_set<term_t> > &edges() { return m_edges; };
    
    inline void addEquality(term_t x, term_t y) {
      m_edges[x].insert(y);
      m_edges[y].insert(x);
      m_vertices.insert(x);
      m_vertices.insert(y);
    }

  };

}

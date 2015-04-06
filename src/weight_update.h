
#include "util.h"
#include "storage.h"

#include <fstream>

namespace scw {
  enum update_type_t {
    SCW_I,
    SCW_II
  };
  
  float updateWeightVector(util::sparse_vector_t       *pOutMean, util::sparse_vector_t *pOutVariance,
                          const util::sparse_vector_t  &vGold, const util::sparse_vector_t &vCompetitor,
                          float                         C, float eta, float initVar, update_type_t typeUpdate,
                          const storage_t<std::string> &stTranslator,
                          std::ostream                 *pLog);
  
  inline float innerProduct(const util::sparse_vector_t &v1, const util::sparse_vector_t &v2) {
    float ret = 0;
    for(auto v: v1)
      ret += v.second * (v2.end() == v2.find(v.first) ? 0.0 : v2.at(v.first));
    return ret;
  }
  
}

namespace sp {
  float updateWeightVector(util::sparse_vector_t       *pOutMean,
                           const util::sparse_vector_t  &vGold, const util::sparse_vector_t &vCompetitor,
                           float eta,
                           const storage_t<std::string> &stTranslator,
                           std::ostream                 *pLog);
  
}

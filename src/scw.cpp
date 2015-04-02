
#include "phillip.h"
#include "util.h"
#include "storage.h"
#include "normsinv.h"

#include "scw.h"

namespace scw {
  float updateWeightVector(util::sparse_vector_t *pOutMean, util::sparse_vector_t *pOutVariance, const util::sparse_vector_t &vGold, const util::sparse_vector_t &vCompetitor, float C, float eta, float initVar, update_type_t typeUpdate, const storage_t<std::string> &stTranslator, std::ostream *pLog) {
    float m, n, v, gamma, phi, u, _u, loss;
    float psi, zeta;
    float alpha, beta;

    util::sparse_vector_t vDiff, vWeightedDiff;
    hash_set<int> featureSet;

    auto get = [](util::sparse_vector_t v, int i, float def = 0.0){ return v.end() != v.find(i) ? v.at(i) : def; };
  
    for(auto it: *pOutMean)     featureSet.insert(it.first);
    for(auto it: *pOutVariance) featureSet.insert(it.first);
    for(auto it: vGold)         featureSet.insert(it.first);
    for(auto it: vCompetitor)   featureSet.insert(it.first);

    for(auto it: featureSet) {
      if(pOutVariance->end() == pOutVariance->find(it)) (*pOutVariance)[it] = initVar;
      
      vDiff[it]         = get(vGold, it) - get(vCompetitor, it);
      vWeightedDiff[it] = get(*pOutVariance, it) * vDiff[it];
    }
    
    v     = innerProduct(vDiff, vWeightedDiff);
    m     = innerProduct(*pOutMean, vDiff);
    loss  = std::max(0.0f, sqrtf(v) - m);

    if(loss == 0.0f) return 0.0f;
    
    phi   = normsinv(eta);
    psi   = 1 + (phi*phi) / 2.0;
    zeta  = 1 + phi*phi;
    n     = v + 1.0 / (2*C);
    gamma = phi * sqrtf(phi*phi*m*m*v*v + 4*n*v*(n+v*phi*phi));
    alpha = SCW_I == typeUpdate ?
      std::min(C, std::max(0.0f, 1.0f / (v*zeta) * (-m * psi + sqrtf(m*m*((phi*phi*phi*phi)/4.0) + v*phi*phi*zeta)))) :
      std::max(0.0f, (-(2*m*v + phi*phi*m*v) + gamma)/(2 * (n*n + n*v*phi*phi)));
    _u    = (-alpha*v*phi + sqrtf(alpha*alpha*v*v*phi*phi + 4*v));
    u     = (1.0 / 4.0) * (_u*_u);
    beta  = (alpha*phi) / (sqrtf(u) + v*alpha*phi);

    for(auto it: featureSet) {
      
      if(0.0f != vDiff[it]) 
        (*pLog) << format("<update id=\"%s\" before_m=\"%.2f\" before_v=\"%.2f\"", format("%d", it).c_str(),
                          // (stTranslator.storage().end() == stTranslator.storage().find(format("%d", it)) ?
                          //  format("%d", it) : stTranslator.storage().find(format("%d", it))->second).c_str(),
                          (*pOutMean)[it], (*pOutVariance)[it]) << " ";
      
      (*pOutMean)[it]     += alpha * (*pOutVariance)[it]*vDiff[it];

      if(vDiff[it] != 0.0)
        (*pOutVariance)[it] -= beta * (*pOutVariance)[it]*v;

      if(0.0f != vDiff[it]) 
        (*pLog) << format("after_m=\"%.2f\" after_v=\"%.2f\" />", (*pOutMean)[it], (*pOutVariance)[it]) << std::endl;
    }

    return loss;
  }
}

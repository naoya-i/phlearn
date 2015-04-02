
namespace scw {
  
  // Powered by: http://www.johndcook.com/blog/cpp_phi/
  inline float ndPhi(float x) {
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
      sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return (float)(0.5*(1.0 + sign*y));
  }
  
  inline float innerProduct(const util::sparse_vector_t &v1, const util::sparse_vector_t &v2) {
    float ret = 0;
    for(auto v: v1)
      ret += v.second * (v2.end() == v2.find(v.first) ? 0.0 : v2.at(v.first));
    return ret;
  }

  enum update_type_t {
    SCW_I,
    SCW_II
  };
  
  inline void updateWeightVector(util::sparse_vector_t *pOutMean, util::sparse_vector_t *pOutVariance, const util::sparse_vector_t &vGold, const util::sparse_vector_t &vCompetitor, float C, float eta, update_type_t typeUpdate, std::ofstream &ofsLog) {
    float m, n, v, gamma, phi, u, _u;
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
      vDiff[it]         = get(vGold, it) - get(vCompetitor, it);
      vWeightedDiff[it] = get(*pOutVariance, it, 1) * (get(vGold, it) - get(vCompetitor, it));
    }
  
    v     = innerProduct(vDiff, vWeightedDiff);
    m     = innerProduct(*pOutMean, vDiff);

    if(0.0f == m) return;
    
    phi   = 1.0 / ndPhi(eta);
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
      if(pOutVariance->end() == pOutVariance->find(it))
        (*pOutVariance)[it] = 1.0;

      ofsLog << format("%d: %.2f, %.2f", it, (*pOutMean)[it], (*pOutVariance)[it]) << " ";
      
      (*pOutMean)[it]     += alpha * (*pOutVariance)[it]*(get(vGold,it) - get(vCompetitor,it));
      (*pOutVariance)[it] -= beta * (get(vGold,it) - get(vCompetitor,it)) * (get(vGold,it) - get(vCompetitor,it)); // NOT SURE THIS IS GOOD...
      
      ofsLog << format(" => %.2f, %.2f", (*pOutMean)[it], (*pOutVariance)[it]) << std::endl;
    }
  }
}

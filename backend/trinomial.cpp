#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <algorithm>

/**
 * ============================================================================
 * TRINOMIAL LOCAL VOLATILITY MODEL
 * ============================================================================
 * This engine calculates option prices (European/American) using a
 * trinomial tree where volatility can change over time.
 */

enum OptionType { CALL = 0, PUT = 1 };

// Main configuration structure sent from Python
typedef struct {
  double S0;      // Initial Spot
  double r;       // Risk-free rate
  double T;       // Time to maturity (years)
  double K;       // Strike price
  int isAmerican; // 1: American, 0: European
  enum OptionType type;
  int N;         // Number of time steps
  int M;         // Number of volatility buckets
  double *Theta; // Sigma values for each bucket
  double *tau;   // Time boundaries for buckets
} Params;

// Results to export to Python
typedef struct {
  double price;      // Price at t=0
  int N;             // Steps used
  double *priceTree; // Asset price tree (flattened)
  double *valueTree; // Option value tree (flattened)
} TrinomialResults;

// Internal structure to centralize Grid configuration
// Avoids recalculating u, d, dt, and sigma_ref in multiple functions.
typedef struct {
  double dt;        // Delta t (T/N)
  double u;         // Up factor
  double d;         // Down factor
  double disc;      // Discount factor e^(-r*dt)
  double sigma_ref; // Reference volatility for the fixed grid
} GridConfig;

// =========================================================
// HELPER FUNCTIONS
// =========================================================

/**
 * idx: Calculates linear index in a 1D array to simulate a 2D matrix.
 * i: Time step (0 to N). Horizontal axis.
 * j: Price vertical level (-i to i). Vertical axis.
 * N: Total resolution.
 *
 * Geometry: At step i, the tree has 2*i + 1 nodes.
 * We use a fixed width of 2*N + 1 to simplify access.
 */
static inline int idx(int i, int j, int N) { return i * (2 * N + 1) + (j + N); }

/**
 * sigma_at_step: Selects local volatility Theta[m] for current time t.
 * Implements the local volatility term structure.
 */
static double sigma_at_step(int i, Params p) {
  double t = i * (p.T / p.N);
  for (int m = 0; m < p.M; m++) {
    if (t >= p.tau[m] && t < p.tau[m + 1])
      return p.Theta[m];
  }
  return p.Theta[p.M - 1]; // Fallback to last bucket
}

/**
 * calculate_sigma_ref: Calculates weighted average of volatility.
 * Used to define a fixed grid (Fixed Grid) that is stable.
 */
static double calculate_sigma_ref(Params p) {
  double sum_sigma = 0.0, total_time = 0.0;
  if (p.M <= 0)
    return 0.2;
  for (int m = 0; m < p.M; m++) {
    double dt_interval = p.tau[m + 1] - p.tau[m];
    sum_sigma += p.Theta[m] * (dt_interval < 0 ? 0 : dt_interval);
    total_time += dt_interval;
  }
  return (total_time <= 0.0) ? p.Theta[0] : sum_sigma / total_time;
}

/**
 * get_grid_config: Precalculates constant grid parameters.
 * u = e^(sigma_ref * sqrt(2*dt)) ensures the tree is recombinant and trinomial.
 */
static GridConfig get_grid_config(Params p) {
  GridConfig cfg;
  cfg.dt = p.T / p.N;
  cfg.disc = exp(-p.r * cfg.dt);
  cfg.sigma_ref = calculate_sigma_ref(p);
  cfg.u = exp(cfg.sigma_ref * sqrt(2.0 * cfg.dt));
  cfg.d = 1.0 / cfg.u;
  return cfg;
}

/**
 * trinomial_probs: Solves probability math at each node.
 * pu, pm, pd are calculated to match mean (r) and variance (sigma_local).
 */
static void trinomial_probs(double sigma_i, double r, double dt, double u,
                            double d, double *pu, double *pm, double *pd) {
  double M_factor = exp(r * dt);
  double V = sigma_i * sigma_i * dt + M_factor * M_factor; // Expected E[S^2]

  // Risk-neutral consistency formulas
  *pu = (V - M_factor * (d + 1.0) + d) / ((u - 1.0) * (u - d));
  *pd = (V - M_factor * (u + 1.0) + u) / ((d - 1.0) * (d - u));
  

  // Numerical corrections to keep probabilities in [0, 1]
  if (*pu < 0)
    *pu = 0;
  if (*pu > 1)
    *pu = 1;
  if (*pd < 0)
    *pd = 0;
  if (*pd > 1)
    *pd = 1;

  *pm = 1.0 - *pu - *pd;

  if (*pm < 0.0) {
    *pm = 0.0;
    double sum = *pu + *pd;
    if (sum > 0) {
      *pu /= sum;
      *pd /= sum;
    }
  }
}

static inline double calculate_payoff(double S, double K,
                                      enum OptionType type) {
  double val = (type == CALL ? (S - K) : (K - S));
  return (val > 0) ? val : 0.0;
}

// =========================================================
// INDUCTION ENGINE
// =========================================================

/**
 * computeValueInduction: The heart of the model.
 * Used for both generating the full tree and fast calibration.
 *
 * PT (PriceTree): If null, prices calculated on the fly (saves memory).
 * VT (ValueTree): Array where option values will be stored.
 */
static void computeValueInduction(Params p, GridConfig cfg, double *VT,
                                  double *PT) {
  // 1. Final Condition at T (Maturity)
  for (int j = -p.N; j <= p.N; j++) {
    double S = PT ? PT[idx(p.N, j, p.N)] : (p.S0 * pow(cfg.u, j));
    VT[idx(p.N, j, p.N)] = calculate_payoff(S, p.K, p.type);
  }

  // 2. Backward Induction (Future to Present)
  for (int i = p.N - 1; i >= 0; i--) {
    double sigma_local = sigma_at_step(i, p);
    double pu, pm, pd;
    trinomial_probs(sigma_local, p.r, cfg.dt, cfg.u, cfg.d, &pu, &pm, &pd);

    for (int j = -i; j <= i; j++) {
      // Continuation value (Discounted Risk-Neutral Expectation)
      double cont = cfg.disc * (pu * VT[idx(i + 1, j + 1, p.N)] +
                                pm * VT[idx(i + 1, j, p.N)] +
                                pd * VT[idx(i + 1, j - 1, p.N)]);

      if (p.isAmerican) {
        double S = PT ? PT[idx(i, j, p.N)] : (p.S0 * pow(cfg.u, j));
        double exercise = calculate_payoff(S, p.K, p.type);
        VT[idx(i, j, p.N)] = (cont > exercise) ? cont : exercise;
      } else {
        VT[idx(i, j, p.N)] = cont;
      }
    }
  }
}

// Lightweight pricing version for calibration (returns only price at t=0)
static double trinomialPricerOnly(Params p) {
  GridConfig cfg = get_grid_config(p);
  // Use std::vector for automatic memory management (RAII)
  std::vector<double> V((p.N + 1) * (2 * p.N + 1));
  
  computeValueInduction(
      p, cfg, V.data(),
      nullptr); // Pass nullptr to avoid dependency on physical PriceTree

  return V[idx(0, 0, p.N)];
}

extern "C" {

// Main entry point for /tree endpoint
TrinomialResults runTrinomialModel(Params p) {
  TrinomialResults res = {0, p.N, nullptr, nullptr};
  GridConfig cfg = get_grid_config(p);
  int size = (p.N + 1) * (2 * p.N + 1);

  try {
    res.priceTree = new double[size];
    res.valueTree = new double[size];
  } catch(const std::bad_alloc&) {
    if (res.priceTree) delete[] res.priceTree;
    if (res.valueTree) delete[] res.valueTree;
    res.priceTree = nullptr;
    res.valueTree = nullptr;
    return res;
  }

  if (res.priceTree && res.valueTree) {
    // A. Price Tree Construction (Forward)
    for (int i = 0; i <= p.N; i++) {
      for (int j = -i; j <= i; j++) {
        res.priceTree[idx(i, j, p.N)] = p.S0 * pow(cfg.u, j);
      }
    }
    // B. Value Tree Construction (Backward)
    computeValueInduction(p, cfg, res.valueTree, res.priceTree);
    res.price = res.valueTree[idx(0, 0, p.N)];
  }
  return res;
}

void freeTrinomialResults(TrinomialResults *res) {
  if (res->priceTree) {
    delete[] res->priceTree;
    res->priceTree = nullptr;
  }
  if (res->valueTree) {
    delete[] res->valueTree;
    res->valueTree = nullptr;
  }
}

// =========================================================
// CALIBRATION (NELDER-MEAD)
// =========================================================


static double calibrationError(double *Theta, int M, double lambda,
                               Params *tmpl, double *Klist, double *Tlist,
                               double *Vmarket, int nOptions) {
  double sse = 0.0, penalty = 0.0;
  Params p = *tmpl;
  p.Theta = Theta;

  for (int m = 0; m < M; m++)
    if (Theta[m] < 0.0)
      return 1e9; // Penalty for negative sigma

  for (int k = 0; k < nOptions; k++) {
    p.K = Klist[k];
    p.T = Tlist[k];
    double model_price = trinomialPricerOnly(p);
    double diff = model_price - Vmarket[k];
    sse += diff * diff;
  }

  // Penalize sharp jumps in temporal volatility
  for (int m = 0; m < M - 1; m++) {
    double d = Theta[m + 1] - Theta[m];
    penalty += d * d;
  }
  return sse/nOptions + lambda * penalty;
}

void nelderMead(double *ThetaStart, int M, double lambda, Params *tmpl,
                double *Klist, double *Tlist, double *Vmarket,
                int nOptions, int maxIter, double tol) {
  int n = M, n_pts = n + 1;
  double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;

  std::vector<std::vector<double>> simplex(n_pts, std::vector<double>(n));
  std::vector<double> scores(n_pts);

  // 1. Simplex Initialization
  for (int i = 0; i < n_pts; i++) {
    // Copy from initial array or perturb
    if (i == 0) {
        for(int k=0; k<n; k++) simplex[i][k] = ThetaStart[k];
    } else {
        for(int k=0; k<n; k++) {
            simplex[i][k] = ThetaStart[k];
            if (k == i - 1) {
                simplex[i][k] += (simplex[i][k] != 0 ? simplex[i][k] * 0.05 : 0.01);
            }
        }
    }
    scores[i] = calibrationError(simplex[i].data(), n, lambda, tmpl, Klist, Tlist, Vmarket, nOptions);
  }

  std::vector<double> centroid(n);
  std::vector<double> reflected(n);
  std::vector<double> expanded(n);
  std::vector<double> contracted(n);

  // 2. Optimization Loop
  std::vector<int> idxs(n_pts);

  for (int iter = 0; iter < maxIter; iter++) {
    // Sort candidates
    for (int i = 0; i < n_pts; i++)
      idxs[i] = i;
    
    std::sort(idxs.begin(), idxs.end(), [&](int a, int b) {
        return scores[a] < scores[b];
    });

    if (fabs(scores[idxs[n]] - scores[idxs[0]]) < tol)
      break;

    // Calculate Centroid (average of best n points)
    for (int j = 0; j < n; j++) {
      centroid[j] = 0;
      for (int i = 0; i < n; i++)
        centroid[j] += simplex[idxs[i]][j];
      centroid[j] /= n;
    }

    int worst = idxs[n];
    // Reflect worst point
    for (int j = 0; j < n; j++)
      reflected[j] = centroid[j] + alpha * (centroid[j] - simplex[worst][j]);
    double r_score = calibrationError(reflected.data(), n, lambda, tmpl, Klist, Tlist,
                                      Vmarket, nOptions);

    if (r_score < scores[idxs[n - 1]] && r_score >= scores[idxs[0]]) {
      simplex[worst] = reflected; // Direct copy
      scores[worst] = r_score;
    } else if (r_score < scores[idxs[0]]) {
      // Try Expansion
      for (int j = 0; j < n; j++)
        expanded[j] = centroid[j] + gamma * (reflected[j] - centroid[j]);
      double e_score = calibrationError(expanded.data(), n, lambda, tmpl, Klist, Tlist,
                                        Vmarket, nOptions);
      if (e_score < r_score) {
        simplex[worst] = expanded;
        scores[worst] = e_score;
      } else {
        simplex[worst] = reflected;
        scores[worst] = r_score;
      }
    } else {
      // Try Contraction
      for (int j = 0; j < n; j++) {
        double val = (r_score < scores[worst] ? reflected[j] : simplex[worst][j]);
        contracted[j] = centroid[j] + rho * (val - centroid[j]);
      }
      double c_score = calibrationError(contracted.data(), n, lambda, tmpl, Klist,
                                        Tlist, Vmarket, nOptions);
      if (c_score < scores[idxs[n]]) {
        simplex[worst] = contracted;
        scores[worst] = c_score;
      } else {
        // Shrink
        for (int i = 1; i < n_pts; i++) {
            // idxs[0] is best
            int current_idx = idxs[i];
            int best_idx = idxs[0];
            for (int j = 0; j < n; j++) {
                 simplex[current_idx][j] = simplex[best_idx][j] +
                                          sigma * (simplex[current_idx][j] - simplex[best_idx][j]);
            }
          scores[current_idx] =
              calibrationError(simplex[current_idx].data(), n, lambda, tmpl, Klist, Tlist,
                               Vmarket, nOptions);
        }
      }
    }
  }

  // Return best result found
  int best = 0;
  for (int i = 1; i < n_pts; i++)
    if (scores[i] < scores[best])
      best = i;
  
  for(int k=0; k<n; k++) ThetaStart[k] = simplex[best][k];

}

} 

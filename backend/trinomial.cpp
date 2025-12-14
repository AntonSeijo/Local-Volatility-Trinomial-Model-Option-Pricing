#include "trinomial.hpp"
#include <cfloat>  // For NAN
#include <cmath>   // For exp, sqrt, pow, fabs
#include <cstdio>  // For printf (if needed)
#include <cstdlib> // For malloc, free
#include <cstring> // For memcpy

// =========================================================
// Helper Functions
// =========================================================

// Calculates linear index for a node (i, j) in a flattened array.
// i: time step (0 to N)
// j: vertical position (-i to i)
// The trinomial tree has width 2*i + 1 at step i.
// For memory simplicity, we use a fixed maximum width of 2*N + 1.
static inline int idx(int i, int j, int N) { return i * (2 * N + 1) + (j + N); }

// Gets the local volatility sigma(t) corresponding to step i.
// Finds in which time interval [tau_m, tau_{m+1}) the time t = i*dt falls.
static double sigma_at_step(int i, Params p) {
  double dt = p.T / p.N;
  double t = i * dt;

  for (int m = 0; m < p.M; m++) {
    if (t >= p.tau[m] && t < p.tau[m + 1]) {
      return p.Theta[m];
    }
  }
  return p.Theta[p.M - 1];
}

// Calculates a reference volatility (time-weighted average of Theta) to define
// the fixed grid. This ensures u*d = 1 and provides a stable grid.
static double calculate_sigma_ref(Params p) {
  double sum_sigma = 0.0;
  double total_time = 0.0;

  if (p.M <= 0 || !p.Theta || !p.tau)
    return 0.2; // Default fallback

  for (int m = 0; m < p.M; m++) {
    // Current interval duration: tau[m+1] - tau[m]
    // Note: tau array has size M+1
    double dt_interval = p.tau[m + 1] - p.tau[m];

    if (dt_interval < 0)
      dt_interval = 0; // Should not happen if tau is sorted

    sum_sigma += p.Theta[m] * dt_interval;
    total_time += dt_interval;
  }

  // If total time is 0 check p.T or return first theta
  if (total_time <= 0.0) {
    return (p.Theta[0] > 0) ? p.Theta[0] : 0.2;
  }

  return sum_sigma / total_time;
}

// Calculates trinomial transition probabilities (pu, pm, pd).
// Ensures probabilities sum to 1 and are within [0, 1].
// Based on standard discretization where V = sigma^2 * dt.
static void trinomial_probs(double sigma_i, double r, double dt, double u,
                            double d, double *pu, double *pm, double *pd) {
  double M_factor = exp(r * dt);

  // V is the variance over dt.
  double V = sigma_i * sigma_i * dt + M_factor * M_factor;

  // prob_u = (V - M(d+1) + d) / ((u-1)(u-d))
  // prob_d = (V - M(u+1) + u) / ((d-1)(d-u))
  // prob_m = 1 - prob_u - prob_d

  double den_ud = (u - 1.0) * (u - d);
  double den_du = (d - 1.0) * (d - u);

  *pu = (V - M_factor * (d + 1.0) + d) / den_ud;
  *pd = (V - M_factor * (u + 1.0) + u) / den_du;
  *pm = 1.0 - *pu - *pd;

  // Numerical correction to ensure valid probabilities
  // With sigma_ref >= sigma_i, probabilities should generally be positive,
  // but floating point issues or extreme r can cause slight violations.
  if (*pu < 0.0)
    *pu = 0.0;
  if (*pu > 1.0)
    *pu = 1.0;
  if (*pd < 0.0)
    *pd = 0.0;
  if (*pd > 1.0)
    *pd = 1.0;

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

// =========================================================
// Trinomial Tree Logic
// =========================================================

// Builds the Price Tree (Forward Induction).
// S_{i, j} = S0 * u^j using fixed u derived from sigma_ref
static void buildPriceTree(Params p, double *PT) {
  if (!PT)
    return;

  double dt = p.T / p.N;

  double sigma_ref = calculate_sigma_ref(p);
  double u = exp(sigma_ref * sqrt(2.0 * dt));

  for (int i = 0; i <= p.N; i++) {
    for (int j = -i; j <= i; j++) {
      PT[idx(i, j, p.N)] = p.S0 * pow(u, j);
    }
  }
}

// Calculates the Option Value Tree (Backward Induction).
// V_{i, j} = DF * (pu*V_{i+1, j+1} + pm*V_{i+1, j} + pd*V_{i+1, j-1})
static void buildValueTree(Params p, double *VT, double *PT) {

  if (!VT)
    return;

  double dt = p.T / p.N;
  double disc = exp(-p.r * dt); // Discount factor

  double sigma_ref = calculate_sigma_ref(p);
  double u = exp(sigma_ref * sqrt(2.0 * dt));
  double d = 1.0 / u;

  // We need the Price Tree for:
  // 1. Payoff at maturity (t=T)
  // 2. Early exercise payoff (American)
  double *localPT = nullptr;
  if (!PT) {
    localPT = (double *)malloc((p.N + 1) * (2 * p.N + 1) * sizeof(double));
    buildPriceTree(p, localPT);
    PT = localPT;
  }

  // 1. Initialize final conditions at t=T (Maturity)
  {
    for (int j = -p.N; j <= p.N; j++) {

      double S = PT[idx(p.N, j, p.N)];

      double val = (p.type == CALL ? (S - p.K) : (p.K - S));
      VT[idx(p.N, j, p.N)] = (val > 0) ? val : 0.0;
    }
  }

  // 2. Backward Induction from T-dt to 0
  for (int i = p.N - 1; i >= 0; i--) {

    double sigma_local = sigma_at_step(i, p);

    double pu, pm, pd;
    // CRITICAL: Compute probabilities using LOCAL sigma but FIXED grid (u, d)
    trinomial_probs(sigma_local, p.r, dt, u, d, &pu, &pm, &pd);

    for (int j = -i; j <= i; j++) {
      // Continuation value (European)
      double cont = disc * (pu * VT[idx(i + 1, j + 1, p.N)] +
                            pm * VT[idx(i + 1, j, p.N)] +
                            pd * VT[idx(i + 1, j - 1, p.N)]);

      // Check early exercise for American Option
      if (p.isAmerican) {
        double S_node = PT[idx(i, j, p.N)];
        double intr = (p.type == CALL ? (S_node - p.K) : (p.K - S_node));
        if (intr < 0)
          intr = 0;
        VT[idx(i, j, p.N)] = (cont > intr) ? cont : intr;
      } else {
        VT[idx(i, j, p.N)] = cont;
      }
    }
  }

  if (localPT)
    free(localPT);
}

// Simplified pricer version for calibration.

static double trinomialPricerOnly(Params p) {
  double dt = p.T / p.N;
  double disc = exp(-p.r * dt);

  int W = 2 * p.N + 1;
  double *V = (double *)malloc((p.N + 1) * W * sizeof(double));

  if (!V)
    return NAN;

  double sigma_ref = calculate_sigma_ref(p);
  double u = exp(sigma_ref * sqrt(2.0 * dt));
  double d = 1.0 / u;

  // Final conditions
  {
    // At T, calculating S is easy: S0 * u^j
    for (int j = -p.N; j <= p.N; j++) {
      double S = p.S0 * pow(u, j); // Calculate S on the fly using FIXED u
      double val = (p.type == CALL ? (S - p.K) : (p.K - S));
      V[idx(p.N, j, p.N)] = (val > 0) ? val : 0.0;
    }
  }

  // Backward steps
  for (int i = p.N - 1; i >= 0; i--) {

    double sigma_local = sigma_at_step(i, p);

    double pu, pm, pd;
    // Probabilities adapt to local sigma on fixed grid
    trinomial_probs(sigma_local, p.r, dt, u, d, &pu, &pm, &pd);

    for (int j = -i; j <= i; j++) {

      double cont =
          disc * (pu * V[idx(i + 1, j + 1, p.N)] + pm * V[idx(i + 1, j, p.N)] +
                  pd * V[idx(i + 1, j - 1, p.N)]);

      if (p.isAmerican) {
        double S = p.S0 * pow(u, j);
        double intr = (p.type == CALL ? (S - p.K) : (p.K - S));
        if (intr < 0)
          intr = 0;
        V[idx(i, j, p.N)] = (cont > intr) ? cont : intr;
      } else {
        V[idx(i, j, p.N)] = cont;
      }
    }
  }

  double res = V[idx(0, 0, p.N)];
  free(V);
  return res;
}

// =========================================================
// Exported Functions (API)
// =========================================================

// Main Function: Runs the full model and returns structures.
// Used by /price and /tree endpoints.
TrinomialResults runTrinomialModel(Params p) {
  TrinomialResults res;
  res.N = p.N;
  int size = (p.N + 1) * (2 * p.N + 1);

  // Allocate memory for trees
  res.priceTree = (double *)malloc(size * sizeof(double));
  res.valueTree = (double *)malloc(size * sizeof(double));

  if (!res.priceTree || !res.valueTree) {
    if (res.priceTree)
      free(res.priceTree);
    if (res.valueTree)
      free(res.valueTree);
    res.price = NAN;
    res.priceTree = nullptr;
    res.valueTree = nullptr;
    return res;
  }

  // Build trees
  buildPriceTree(p, res.priceTree);
  buildValueTree(p, res.valueTree, res.priceTree);

  // Option price is the value at the root (t=0, j=0)
  res.price = res.valueTree[idx(0, 0, p.N)];

  return res;
}

// Releases memory of results returned to Python/API.
void freeTrinomialResults(TrinomialResults *res) {
  if (res && res->priceTree) {
    free(res->priceTree);
    res->priceTree = nullptr;
  }
  if (res && res->valueTree) {
    free(res->valueTree);
    res->valueTree = nullptr;
  }
}

// =========================================================
// Calibration (Local Volatility)
// =========================================================

// Objective Function: Squared Error + Penalty (Smoothing).
// E(Theta) = sum(w * (Model - Market)^2) + lambda * penalty
static double calibrationError(double *Theta, int M, double lambda,
                               Params *tmpl, double *Klist, double *Tlist,
                               double *Vmarket, double *w, int nOptions) {
  double sse = 0.0;
  Params p = *tmpl;
  p.Theta = Theta; // Use candidate Thetas

  for (int k = 0; k < nOptions; k++) {
    p.K = Klist[k];
    p.T = Tlist[k];
    p.isAmerican = 0; // Assume calibration with liquid European options

    // Use optimized pricer (no tree storage)
    double modelPrice = trinomialPricerOnly(p);
    double diff = modelPrice - Vmarket[k];
    sse += w[k] * (diff * diff);
  }

  // Penalty for excessive volatility variation (Tikhonov variational
  // regularization)
  double penalty = 0.0;
  for (int m = 0; m < M - 1; m++) {
    double d = Theta[m + 1] - Theta[m];
    penalty += d * d;
  }

  return sse + lambda * penalty;
}

// Nelder-Mead Algorithm (Simplex) to minimize calibration error.
// Finds optimal Theta vector.
void nelderMead(double *ThetaStart, int M, double lambda, Params *tmpl,
                double *Klist, double *Tlist, double *Vmarket, double *w,
                int nOptions, int maxIter, double tol) {

  // Normalize weights so they sum to 1.0
  double sum_w = 0.0;
  for (int i = 0; i < nOptions; i++) {
    sum_w += w[i];
  }

  if (sum_w > 0.0) {
    for (int i = 0; i < nOptions; i++) {
      w[i] /= sum_w;
    }
  } else {
    // Handle zero sum: Set uniform weights
    double uniform = 1.0 / nOptions;
    for (int i = 0; i < nOptions; i++)
      w[i] = uniform;
  }

  // Standard Nelder-Mead parameters
  double alpha = 1.0, gamma = 2.0, rho = 0.5, sigma = 0.5;

  int n = M;            // Problem dimension (number of volatility buckets)
  int n_points = n + 1; // Points in simplex

  // Memory allocation for simplex
  double **simplex = (double **)malloc(n_points * sizeof(double *));
  double *scores = (double *)malloc(n_points * sizeof(double));

  // Initial point (ThetaStart provided by user)
  simplex[0] = (double *)malloc(n * sizeof(double));
  memcpy(simplex[0], ThetaStart, n * sizeof(double));
  scores[0] = calibrationError(simplex[0], n, lambda, tmpl, Klist, Tlist,
                               Vmarket, w, nOptions);

  // Generate remaining initial points by perturbing each dimension
  for (int i = 1; i < n_points; i++) {
    simplex[i] = (double *)malloc(n * sizeof(double));
    memcpy(simplex[i], ThetaStart, n * sizeof(double));
    // Perturbation of 5% or 0.01 absolute if zero
    double step = (simplex[i][i - 1] != 0) ? simplex[i][i - 1] * 0.05 : 0.01;
    simplex[i][i - 1] += step;
    scores[i] = calibrationError(simplex[i], n, lambda, tmpl, Klist, Tlist,
                                 Vmarket, w, nOptions);
  }

  // Temp buffers for transformed points
  double *centroid = (double *)malloc(n * sizeof(double));
  double *reflected = (double *)malloc(n * sizeof(double));
  double *expanded = (double *)malloc(n * sizeof(double));
  double *contracted = (double *)malloc(n * sizeof(double));

  // Main Optimization Loop
  for (int iter = 0; iter < maxIter; iter++) {
    // 1. Sort by score (Bubble sort for simplicity with small N)
    int indices[n_points];
    for (int i = 0; i < n_points; i++)
      indices[i] = i;

    for (int i = 0; i < n_points - 1; i++) {
      for (int j = 0; j < n_points - i - 1; j++) {
        if (scores[indices[j]] > scores[indices[j + 1]]) {
          int temp = indices[j];
          indices[j] = indices[j + 1];
          indices[j + 1] = temp;
        }
      }
    }

    // Convergence check
    double best_score = scores[indices[0]];
    double worst_score = scores[indices[n]];
    if (fabs(worst_score - best_score) < tol)
      break;

    // 2. Calculate Centroid (excluding worst point)
    for (int j = 0; j < n; j++)
      centroid[j] = 0.0;
    for (int i = 0; i < n; i++) {
      int idx_p = indices[i];
      for (int j = 0; j < n; j++)
        centroid[j] += simplex[idx_p][j];
    }
    for (int j = 0; j < n; j++)
      centroid[j] /= n;

    int worst_idx = indices[n];

    // 3. Reflect
    for (int j = 0; j < n; j++)
      reflected[j] =
          centroid[j] + alpha * (centroid[j] - simplex[worst_idx][j]);
    double reflected_score = calibrationError(reflected, n, lambda, tmpl, Klist,
                                              Tlist, Vmarket, w, nOptions);

    if (reflected_score >= scores[indices[0]] &&
        reflected_score < scores[indices[n - 1]]) {
      memcpy(simplex[worst_idx], reflected, n * sizeof(double));
      scores[worst_idx] = reflected_score;
      continue;
    }

    // 4. Expand
    if (reflected_score < scores[indices[0]]) {
      for (int j = 0; j < n; j++)
        expanded[j] = centroid[j] + gamma * (reflected[j] - centroid[j]);
      double expanded_score = calibrationError(expanded, n, lambda, tmpl, Klist,
                                               Tlist, Vmarket, w, nOptions);

      if (expanded_score < reflected_score) {
        memcpy(simplex[worst_idx], expanded, n * sizeof(double));
        scores[worst_idx] = expanded_score;
      } else {
        memcpy(simplex[worst_idx], reflected, n * sizeof(double));
        scores[worst_idx] = reflected_score;
      }
      continue;
    }

    // 5. Contract
    if (reflected_score < scores[indices[n]]) {
      // Outside Contraction
      for (int j = 0; j < n; j++)
        contracted[j] = centroid[j] + rho * (reflected[j] - centroid[j]);
      double contracted_score = calibrationError(
          contracted, n, lambda, tmpl, Klist, Tlist, Vmarket, w, nOptions);
      if (contracted_score < reflected_score) {
        memcpy(simplex[worst_idx], contracted, n * sizeof(double));
        scores[worst_idx] = contracted_score;
        continue;
      }
    } else {
      // Inside Contraction
      for (int j = 0; j < n; j++)
        contracted[j] =
            centroid[j] + rho * (simplex[worst_idx][j] - centroid[j]);
      double contracted_score = calibrationError(
          contracted, n, lambda, tmpl, Klist, Tlist, Vmarket, w, nOptions);
      if (contracted_score < scores[indices[n]]) {
        memcpy(simplex[worst_idx], contracted, n * sizeof(double));
        scores[worst_idx] = contracted_score;
        continue;
      }
    }

    // 6. Shrink entire simplex towards best point
    int best_idx = indices[0];
    for (int i = 1; i < n_points; i++) {
      int idx_p = indices[i];
      for (int j = 0; j < n; j++) {
        simplex[idx_p][j] = simplex[best_idx][j] +
                            sigma * (simplex[idx_p][j] - simplex[best_idx][j]);
      }
      scores[idx_p] = calibrationError(simplex[idx_p], n, lambda, tmpl, Klist,
                                       Tlist, Vmarket, w, nOptions);
    }
  }

  // Copy best result to ThetaStart (output)
  int final_best_idx = 0;
  double min_s = scores[0];
  for (int i = 1; i < n_points; i++) {
    if (scores[i] < min_s) {
      min_s = scores[i];
      final_best_idx = i;
    }
  }
  memcpy(ThetaStart, simplex[final_best_idx], n * sizeof(double));

  // Memory Cleanup
  for (int i = 0; i < n_points; i++)
    free(simplex[i]);
  free(simplex);
  free(scores);
  free(centroid);
  free(reflected);
  free(expanded);
  free(contracted);
}

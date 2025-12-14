// =========================================================
// TRINOMIAL.H
// Definitions for the Trinomial Local Volatility Model.
// Compatible with C and C++.
// =========================================================

#ifndef TRINOMIAL_H
#define TRINOMIAL_H

// Option Type: Call (Buy Right) or Put (Sell Right).
enum OptionType { CALL = 0, PUT = 1 };

// Main structure containing ALL model configuration.
// Passed directly from Python (ctypes).
typedef struct {
  double S0; // Initial Spot price of the underlying asset.
  double r;  // Risk-free interest rate (annualized).

  // Option Contract Parameters
  double T;             // Time to Maturity (years).
  double K;             // Strike Price.
  int isAmerican;       // 1: American Option (early exercise), 0: European.
  enum OptionType type; // CALL or PUT.

  // Numerical Parameters
  int N; // Number of time steps in the trinomial tree.

  // Local Volatility (term structure)
  // Allows defining different volatilities for time buckets.
  int M;         // Number of volatility buckets.
  double *Theta; // Vector of volatilities (sigma) for each bucket.
  double *tau;   // Vector of times where volatility changes (Boundaries).

} Params;

// Result structure returned to the client.
// Includes price and full trees for visualization.
typedef struct {
  double price;      // Calculated Option Price (at t=0).
  int N;             // Steps used (to dimension arrays in Python).
  double *priceTree; // Underlying Price Tree (Flattened Vector).
  double *valueTree; // Option Value Tree (Flattened Vector).
} TrinomialResults;

#ifdef __cplusplus
extern "C" {
#endif

// =========================================================
// Public Functions (API Exports)
// =========================================================

// Runs the full Trinomial Model.
// - Allocates memory for price and value trees.
// - Calculates forward nodes (prices) and backward nodes (values).
// - Returns TrinomialResults structure (Must be freed with
// freeTrinomialResults!).
TrinomialResults runTrinomialModel(Params p);

// Frees the memory of the trees within TrinomialResults.
// Must be called from Python after consuming data to avoid leaks.
void freeTrinomialResults(TrinomialResults *res);

// Executes Local Volatility Calibration using Nelder-Mead.
// - Adjusts Theta vector to minimize difference between model and market
// prices.
// - Theta0: Initial point. Updated in-place with optimal result.
// - maxIter/tol: Optimization stopping criteria.
void nelderMead(double *Theta0, int M, double lambda, Params *tmpl,
                double *Klist, double *Tlist, double *Vmarket, double *w,
                int nOptions, int maxIter, double tol);

#ifdef __cplusplus
}
#endif

#endif

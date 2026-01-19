// =========================================================
// TRINOMIAL.H
// Definitions for the Trinomial Local Volatility Model.
// Compatible with C and C++.
// =========================================================

#ifndef TRINOMIAL_H
#define TRINOMIAL_H

// Option Type: CALL (0) or PUT (1).
enum OptionType { CALL = 0, PUT = 1 };

/**
 * Params: Main configuration structure.
 * This is the bridge between Python (ctypes) and C++.
 */
typedef struct {
  double S0;            // Initial Spot price of the asset.
  double r;             // Annualized risk-free interest rate.
  double T;             // Time to maturity in years.
  double K;             // Strike price of the option.
  int isAmerican;       // 1 for American (early exercise), 0 for European.
  enum OptionType type; // CALL or PUT.
  int N;                // Number of time steps (Resolution).

  // Local Volatility Term Structure
  int M;         // Number of volatility time-buckets.
  double *Theta; // Sigma values for each bucket.
  double *tau;   // Time boundaries for volatility changes.
} Params;

/**
 * TrinomialResults: Output structure.
 * Contains the final price and the flattened trees for visualization.
 */
typedef struct {
  double price;      // Option price at t=0.
  int N;             // Number of steps (used for array dimensioning).
  double *priceTree; // Flattened vector of the underlying price tree.
  double *valueTree; // Flattened vector of the option value tree.
} TrinomialResults;

#ifdef __cplusplus
extern "C" {
#endif

// =========================================================
// PUBLIC API EXPORTS
// =========================================================

/**
 * runTrinomialModel
 * Executes the full engine:
 * 1. Pre-calculates the fixed grid geometry (GridConfig).
 * 2. Builds the price tree (Forward induction).
 * 3. Builds the value tree (Backward induction).
 *
 * Returns a structure with the calculated price and pointers to the trees.
 * IMPORTANT: Memory must be released using freeTrinomialResults.
 */
TrinomialResults runTrinomialModel(Params p);

/**
 * freeTrinomialResults
 * Releases the heap memory allocated for priceTree and valueTree.
 * Must be called from the client (e.g., Python) to prevent memory leaks.
 */
void freeTrinomialResults(TrinomialResults *res);

/**
 * nelderMead
 * Calibration algorithm to find optimal Local Volatility (Theta).
 *
 * @param Theta0: Initial sigma vector. Will be updated with optimal results.
 * @param M: Dimension of the optimization problem.
 * @param lambda: Tikhonov regularization parameter (smoothing).
 * @param tmpl: Template parameters (Spot, rate, steps).
 * @param Klist/Tlist/Vmarket: Market data for calibration.
 * @param w: Weights for each market option.
 * @param nOptions: Size of market data arrays.
 */
void nelderMead(double *Theta0, int M, double lambda, Params *tmpl,
                double *Klist, double *Tlist, double *Vmarket,
                int nOptions, int maxIter, double tol);

#ifdef __cplusplus
}
#endif

#endif
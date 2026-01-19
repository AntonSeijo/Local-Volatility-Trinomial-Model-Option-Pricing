from fastapi import FastAPI
from fastapi.responses import JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Optional
import ctypes
import numpy as np
import os
from dotenv import load_dotenv

# ==========================================================
# Configuration and Environment
# ==========================================================
load_dotenv()

# Get frontend URL from environment for CORS policy
FRONTEND_URL = os.getenv("FRONTEND_URL", "http://localhost:5173")

app = FastAPI(title="Trinomial Local Volatility Engine")

# Configure CORS
allowed_origins = [
    FRONTEND_URL,
    "http://localhost:5173",
    "http://127.0.0.1:5173",
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=allowed_origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ==========================================================
# Load C Library (Shared Object)
# ==========================================================
try:
    # Standard path for Docker/Production
    lib = ctypes.CDLL("./trinomial_model.so")
except OSError:
    # Fallback path for local development
    lib = ctypes.CDLL("./backend/trinomial_model.so")


# ==========================================================
# C-Structure Mappings (ctypes)
# ==========================================================

class Params(ctypes.Structure):
    """Matches the 'Params' struct in C++."""
    _fields_ = [
        ("S0", ctypes.c_double),
        ("r", ctypes.c_double),
        ("T", ctypes.c_double),
        ("K", ctypes.c_double),
        ("isAmerican", ctypes.c_int),
        ("type", ctypes.c_int),  # 0=CALL, 1=PUT
        ("N", ctypes.c_int),
        ("M", ctypes.c_int),
        ("Theta", ctypes.POINTER(ctypes.c_double)),
        ("tau", ctypes.POINTER(ctypes.c_double))
    ]

class TrinomialResults(ctypes.Structure):
    """Matches the 'TrinomialResults' struct in C++."""
    _fields_ = [
        ("price", ctypes.c_double),
        ("N", ctypes.c_int),
        ("priceTree", ctypes.POINTER(ctypes.c_double)),
        ("valueTree", ctypes.POINTER(ctypes.c_double))
    ]

# Define Function Signatures for C Library
lib.runTrinomialModel.argtypes = [Params]
lib.runTrinomialModel.restype = TrinomialResults

lib.freeTrinomialResults.argtypes = [ctypes.POINTER(TrinomialResults)]
lib.freeTrinomialResults.restype = None

lib.nelderMead.argtypes = [
    ctypes.POINTER(ctypes.c_double),     # ThetaStart
    ctypes.c_int,                       # M
    ctypes.c_double,                    # lambda
    ctypes.POINTER(Params),             # template
    ctypes.POINTER(ctypes.c_double),     # Klist
    ctypes.POINTER(ctypes.c_double),     # Tlist
    ctypes.POINTER(ctypes.c_double),     # Vmarket
    ctypes.c_int,                       # nOptions
    ctypes.c_int,                       # maxIter
    ctypes.c_double                     # tol
]
lib.nelderMead.restype = None


# ==========================================================
# Pydantic Schemas for Request Validation
# ==========================================================

class ModelParams(BaseModel):
    S0: float
    K: float
    T: float
    r: float
    N: int
    type: str  # "CALL" or "PUT"
    isAmerican: bool
    Theta: List[float]
    tau: List[float]

class CalibrationParams(BaseModel):
    S0: float
    r: float
    N: int
    M: int
    Theta_initial: List[float]
    tau: List[float]
    Klist: List[float]
    Tlist: List[float]
    Vmarket: List[float]
    lambda_penalty: float = 0.01
    max_iter: int = 500
    tolerance: float = 1e-6


# ==========================================================
# Helper Functions
# ==========================================================

def get_c_params(input_data: ModelParams) -> Params:
    """Converts Pydantic Model to C-Compatible Structure."""
    M = len(input_data.Theta)
    
    # Create C-compatible double arrays
    ThetaArray = (ctypes.c_double * M)(*input_data.Theta)
    TauArray = (ctypes.c_double * len(input_data.tau))(*input_data.tau)
    
    option_type = 0 if input_data.type.upper() == "CALL" else 1
    
    p = Params()
    p.S0 = input_data.S0
    p.r = input_data.r
    p.T = input_data.T
    p.K = input_data.K
    p.isAmerican = 1 if input_data.isAmerican else 0
    p.type = option_type
    p.N = input_data.N
    p.M = M
    p.Theta = ctypes.cast(ThetaArray, ctypes.POINTER(ctypes.c_double))
    p.tau = ctypes.cast(TauArray, ctypes.POINTER(ctypes.c_double))
    
    # CRITICAL: Keep reference to arrays to prevent Python Garbage Collector 
    # from freeing memory while C is still reading it.
    p._refs = (ThetaArray, TauArray) 
    
    return p


# ==========================================================
# API Endpoints
# ==========================================================

@app.post("/tree")
def calculate_tree(params: ModelParams):
    """Runs the pricing model and returns the full nested tree structure."""
    c_params = get_c_params(params)
    
    # Execute C++ Engine
    results = lib.runTrinomialModel(c_params)
    
    if np.isnan(results.price):
        lib.freeTrinomialResults(ctypes.byref(results))
        return JSONResponse(status_code=400, content={"error": "Calculation failed (NaN)"})
    
    N = results.N
    width = 2 * N + 1
    total_size = (N + 1) * width
    
    # Safely convert C pointers to Numpy arrays and copy them to Python managed memory
    price_flat = np.ctypeslib.as_array(results.priceTree, shape=(total_size,)).copy()
    value_flat = np.ctypeslib.as_array(results.valueTree, shape=(total_size,)).copy()
    
    # Immediately free C++ memory after copying
    lib.freeTrinomialResults(ctypes.byref(results))
    
    # Reconstruct the 1D flat arrays into a nested list (Tree structure)
    # Level i contains 2*i + 1 nodes.
    price_tree = []
    value_tree = []
    
    for i in range(N + 1):
        level_prices = []
        level_values = []
        for j in range(-i, i + 1):
            # idx(i, j) mapping logic
            idx_mapped = i * width + (j + N)
            level_prices.append(float(price_flat[idx_mapped]))
            level_values.append(float(value_flat[idx_mapped]))
        price_tree.append(level_prices)
        value_tree.append(level_values)

    return {
        "price": float(results.price),
        "priceTree": price_tree,
        "valueTree": value_tree
    }


@app.post("/calibrate")
def calibrate(params: CalibrationParams):
    """Executes Nelder-Mead optimization to find the best Volatility Theta values."""
    
    # Validation
    if len(params.Theta_initial) != params.M:
        return JSONResponse(status_code=400, content={"error": "Theta_initial length must equal M"})
    if len(params.tau) != params.M + 1:
        return JSONResponse(status_code=400, content={"error": "tau length must equal M+1"})
    
    nOptions = len(params.Klist)
    
    # Prepare C Arrays
    ThetaArray = (ctypes.c_double * params.M)(*params.Theta_initial)
    TauArray = (ctypes.c_double * len(params.tau))(*params.tau)
    KlistArray = (ctypes.c_double * nOptions)(*params.Klist)
    TlistArray = (ctypes.c_double * nOptions)(*params.Tlist)
    VmarketArray = (ctypes.c_double * nOptions)(*params.Vmarket)
    
    # Create template Params for the pricer within the calibration loop
    template = Params()
    template.S0, template.r, template.N, template.M = params.S0, params.r, params.N, params.M
    template.tau = ctypes.cast(TauArray, ctypes.POINTER(ctypes.c_double))
    template.type, template.isAmerican = 0, 0 # Calibrate using Euro Calls typically
    
    # Run Optimization (Nelder-Mead)
    # ThetaArray will be modified in-place by the C++ function
    lib.nelderMead(
        ThetaArray,
        params.M,
        params.lambda_penalty,
        ctypes.byref(template),
        KlistArray,
        TlistArray,
        VmarketArray,
        nOptions,
        params.max_iter,
        params.tolerance
    )
    
    return {
        "calibrated_theta": [ThetaArray[i] for i in range(params.M)],
        "tau": params.tau,
        "message": "Calibration completed successfully"
    }
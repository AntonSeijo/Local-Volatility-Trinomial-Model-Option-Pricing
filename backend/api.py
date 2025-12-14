from fastapi import FastAPI
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import List, Optional
import ctypes
import numpy as np
import os
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Get frontend URL from environment (for CORS)
FRONTEND_URL = os.getenv("FRONTEND_URL", "http://localhost:5173")

app = FastAPI()

from fastapi.middleware.cors import CORSMiddleware

# Configure CORS with environment variable
allowed_origins = [
    FRONTEND_URL,
    "http://localhost:5173",  # Fallback for local dev
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
# Load C Library
# ==========================================================
try:
    lib = ctypes.CDLL("./trinomial_model.so")
except OSError:
    # Fallback for local testing if not in docker or different path
    lib = ctypes.CDLL("./backend/trinomial_model.so")


# ==========================================================
# C Structures Definition
# ==========================================================

class Params(ctypes.Structure):
    _fields_ = [
        ("S0", ctypes.c_double),
        ("r", ctypes.c_double),
        ("T", ctypes.c_double),
        ("K", ctypes.c_double),
        ("isAmerican", ctypes.c_int),
        ("type", ctypes.c_int), # 0=CALL, 1=PUT
        ("N", ctypes.c_int),
        ("M", ctypes.c_int),
        ("Theta", ctypes.POINTER(ctypes.c_double)),
        ("tau", ctypes.POINTER(ctypes.c_double))
    ]

class TrinomialResults(ctypes.Structure):
    _fields_ = [
        ("price", ctypes.c_double),
        ("N", ctypes.c_int),
        ("priceTree", ctypes.POINTER(ctypes.c_double)),
        ("valueTree", ctypes.POINTER(ctypes.c_double))
    ]

# Configure function signatures
lib.runTrinomialModel.argtypes = [Params]
lib.runTrinomialModel.restype = TrinomialResults

lib.freeTrinomialResults.argtypes = [ctypes.POINTER(TrinomialResults)]
lib.freeTrinomialResults.restype = None

# ==========================================================
# Pydantic Models for API
# ==========================================================

class ModelParams(BaseModel):
    S0: float
    K: float
    T: float
    r: float
    N: int
    type: str # "CALL" or "PUT"
    isAmerican: bool
    tm: Optional[int] = None # M, renaming to avoid confusion but mapping to M
    Theta: List[float]
    tau: List[float]

    class Config:
        json_schema_extra = {
            "example": {
                "S0": 100.0, "K": 100.0, "T": 1.0, "r": 0.05, "N": 50,
                "type": "CALL", "isAmerican": True,
                "u": 1.02, "d": 0.98039,
                "Theta": [0.20, 0.25, 0.22], "tau": [0.0, 0.4, 0.7, 1.0]
            }
        }

# Helpers

def get_c_params(input_data: ModelParams) -> Params:
    # Prepare arrays
    M = len(input_data.Theta)
    
    # Validate tau length (should be M+1 according to C logic, or handled there)
    
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
    p.type = option_type
    p.N = input_data.N
    p.M = M
    p.Theta = ctypes.cast(ThetaArray, ctypes.POINTER(ctypes.c_double))
    p.tau = ctypes.cast(TauArray, ctypes.POINTER(ctypes.c_double))
    
    # Important: Keep references to arrays to prevent GC
    p._refs = (ThetaArray, TauArray) 
    
    return p

# Endpoints

@app.post("/tree")
def calculate_tree(params: ModelParams):
    c_params = get_c_params(params)
    
    results = lib.runTrinomialModel(c_params)
    
    if np.isnan(results.price):
        lib.freeTrinomialResults(ctypes.byref(results))
        return JSONResponse(status_code=400, content={"error": "Calculation failed (NaN)"})
    
    # Extract trees
    # Tree size is (N+1)*(2N+1) according to C, but indexing is idx(i, j)
    # idx(i, j) = i * (2N+1) + (j+N)
    # We want to convert this to a friendly JSON structure.
    # List of lists where each level i has its j nodes.
    
    N = results.N
    width = 2 * N + 1
    total_size = (N + 1) * width
    
    # Convert pointers to numpy arrays (copy)
    # Note: as_array does not own memory, so we copy before free
    price_flat = np.ctypeslib.as_array(results.priceTree, shape=(total_size,)).copy()
    value_flat = np.ctypeslib.as_array(results.valueTree, shape=(total_size,)).copy()
    
    lib.freeTrinomialResults(ctypes.byref(results))
    
    # Structure into levels
    # For level i, j goes from -i to i.
    price_tree = []
    value_tree = []
    
    for i in range(N + 1):
        level_prices = []
        level_values = []
        for j in range(-i, i + 1):
            # Indice lineal
            idx = i * width + (j + N)
            level_prices.append(float(price_flat[idx]))
            level_values.append(float(value_flat[idx]))
        price_tree.append(level_prices)
        value_tree.append(level_values)

    return {
        "price": float(results.price),
        "priceTree": price_tree,
        "valueTree": value_tree
    }

# Endpoint 3: CALIBRATION

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
    weights: Optional[List[float]] = None
    lambda_penalty: float = 0.01
    max_iter: int = 500
    tolerance: float = 1e-6

# Configure nelderMead signature
lib.nelderMead.argtypes = [
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_double,
    ctypes.POINTER(Params),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_double),
    ctypes.c_int,
    ctypes.c_int,
    ctypes.c_double
]
lib.nelderMead.restype = None

@app.post("/calibrate")
def calibrate(params: CalibrationParams):
    if len(params.Theta_initial) != params.M:
        return JSONResponse(status_code=400, content={"error": "Theta_initial length must equal M"})
    
    if len(params.tau) != params.M + 1:
        return JSONResponse(status_code=400, content={"error": "tau length must equal M+1"})
    
    nOptions = len(params.Klist)
    if len(params.Tlist) != nOptions or len(params.Vmarket) != nOptions:
        return JSONResponse(status_code=400, content={"error": "Klist, Tlist, and Vmarket must have same length"})
    
    weights = params.weights if params.weights else [1.0] * nOptions
    if len(weights) != nOptions:
        return JSONResponse(status_code=400, content={"error": "weights length must match number of options"})
    
    M = params.M
    ThetaArray = (ctypes.c_double * M)(*params.Theta_initial)
    TauArray = (ctypes.c_double * len(params.tau))(*params.tau)
    KlistArray = (ctypes.c_double * nOptions)(*params.Klist)
    TlistArray = (ctypes.c_double * nOptions)(*params.Tlist)
    VmarketArray = (ctypes.c_double * nOptions)(*params.Vmarket)
    WeightsArray = (ctypes.c_double * nOptions)(*weights)
    
    template = Params()
    template.S0 = params.S0
    template.r = params.r
    template.N = params.N
    template.M = M
    template.Theta = ctypes.cast(ThetaArray, ctypes.POINTER(ctypes.c_double))
    template.tau = ctypes.cast(TauArray, ctypes.POINTER(ctypes.c_double))
    template.type = 0
    template.isAmerican = 0
    template.K = 0.0
    template.T = 0.0
    
    lib.nelderMead(
        ThetaArray,
        M,
        params.lambda_penalty,
        ctypes.byref(template),
        KlistArray,
        TlistArray,
        VmarketArray,
        WeightsArray,
        nOptions,
        params.max_iter,
        params.tolerance
    )
    
    calibrated_theta = [ThetaArray[i] for i in range(M)]
    
    return {
        "calibrated_theta": calibrated_theta,
        "tau": params.tau,
        "message": "Calibration completed successfully"
    }


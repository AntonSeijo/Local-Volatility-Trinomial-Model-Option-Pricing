import React, { useState } from 'react';
import { motion } from 'framer-motion';
import axios from 'axios';

export default function CalibrationPanel({ params, setParams, style }) {
    const [calibrating, setCalibrating] = useState(false);
    const [calibrationResult, setCalibrationResult] = useState(null);
    const [marketData, setMarketData] = useState({
        Klist: "95, 100, 105",
        Tlist: "1.0, 1.0, 1.0",
        Vmarket: "8.5, 5.2, 3.1"
    });
    const [calibParams, setCalibParams] = useState({
        lambda_penalty: 0.01,
        max_iter: 500,
        tolerance: 0.000001
    });

    const handleMarketDataChange = (e) => {
        setMarketData({ ...marketData, [e.target.name]: e.target.value });
    };

    const handleCalibParamsChange = (e) => {
        if (e.target.type === 'number') {
            setCalibParams({ ...calibParams, [e.target.name]: parseFloat(e.target.value) });
        } else {
            setCalibParams({ ...calibParams, [e.target.name]: e.target.value });
        }
    };

    const runCalibration = async () => {
        setCalibrating(true);
        setCalibrationResult(null);

        try {
            // Parse inputs
            const Theta = params.ThetaStr.split(',').map(s => Number(s.trim()));
            const tau = params.tauStr.split(',').map(s => Number(s.trim()));
            const Klist = marketData.Klist.split(',').map(s => Number(s.trim()));
            const Tlist = marketData.Tlist.split(',').map(s => Number(s.trim()));
            const Vmarket = marketData.Vmarket.split(',').map(s => Number(s.trim()));

            // Basic validation before sending
            if (Theta.length + 1 !== tau.length) {
                // Mismatch in Theta and tau lengths
            }

            const payload = {
                S0: Number(params.S0),
                r: Number(params.r),
                N: Number(params.N),
                M: Theta.length,
                Theta_initial: Theta,
                tau: tau,
                Klist: Klist,
                Tlist: Tlist,
                Vmarket: Vmarket,
                lambda_penalty: calibParams.lambda_penalty,
                max_iter: calibParams.max_iter,
                tolerance: calibParams.tolerance
            };

            const apiUrl = import.meta.env.VITE_API_URL || 'http://localhost:8000';
            const response = await axios.post(`${apiUrl}/calibrate`, payload);

            setCalibrationResult(response.data);

            // Update main params with calibrated values
            const newThetaStr = response.data.calibrated_theta.map(v => v.toFixed(6)).join(', ');
            setParams({ ...params, ThetaStr: newThetaStr });

        } catch (err) {
            console.error(err);
            const msg = err.response?.data?.error || err.message;
            alert("Calibration error: " + msg);
        }
        setCalibrating(false);
    };

    return (
        <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            className="glass-panel"
            style={{ padding: '1.5rem', ...style }} // Merge style prop with default padding
        >
            <h3 style={{ marginTop: 0, marginBottom: '1rem', fontSize: '1.2rem', color: '#fff' }}>
                Market Calibration
            </h3>
            <p style={{ fontSize: '0.8rem', color: '#94a3b8', marginBottom: '1rem' }}>
                Calibrate Local Volatility (Θ) to market prices using Nelder-Mead optimization.
            </p>

            <div style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>

                <div>
                    <label>Market Strikes (K) [comma separated]</label>
                    <input
                        type="text"
                        name="Klist"
                        value={marketData.Klist}
                        onChange={handleMarketDataChange}
                        placeholder="95, 100, 105"
                    />
                </div>

                <div>
                    <label>Market Maturities (T) [comma separated]</label>
                    <input
                        type="text"
                        name="Tlist"
                        value={marketData.Tlist}
                        onChange={handleMarketDataChange}
                        placeholder="1.0, 1.0, 1.0"
                    />
                </div>

                <div>
                    <label>Market Prices (V) [comma separated]</label>
                    <input
                        type="text"
                        name="Vmarket"
                        value={marketData.Vmarket}
                        onChange={handleMarketDataChange}
                        placeholder="8.5, 5.2, 3.1"
                    />
                </div>

                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '0.5rem' }}>
                    <div>
                        <label style={{ fontSize: '0.75rem' }}>Lambda (λ)</label>
                        <input
                            type="number"
                            name="lambda_penalty"
                            value={calibParams.lambda_penalty}
                            onChange={handleCalibParamsChange}
                            step="0.001"
                        />
                    </div>
                    <div>
                        <label style={{ fontSize: '0.75rem' }}>Max Iter</label>
                        <input
                            type="number"
                            name="max_iter"
                            value={calibParams.max_iter}
                            onChange={handleCalibParamsChange}
                        />
                    </div>
                    <div>
                        <label style={{ fontSize: '0.75rem' }}>Tolerance</label>
                        <input
                            type="number"
                            name="tolerance"
                            value={calibParams.tolerance}
                            onChange={handleCalibParamsChange}
                            step="0.000001"
                        />
                    </div>
                </div>

                <button
                    onClick={runCalibration}
                    className="btn-primary"
                    disabled={calibrating}
                    style={{ marginTop: '0.5rem', background: calibrating ? '#475569' : '' }}
                >
                    {calibrating ? 'Running Optimization...' : 'Run Calibration'}
                </button>

                {calibrationResult && (
                    <motion.div
                        initial={{ opacity: 0, y: -10 }}
                        animate={{ opacity: 1, y: 0 }}
                        style={{
                            marginTop: '0.5rem',
                            padding: '0.5rem',
                            background: 'rgba(16, 185, 129, 0.1)',
                            borderRadius: '8px',
                            border: '1px solid rgba(16, 185, 129, 0.3)'
                        }}
                    >
                        <h4 style={{ margin: '0 0 0.5rem 0', color: '#34d399', fontSize: '0.95rem' }}>
                            ✓ Calibration Successful
                        </h4>
                        <div style={{ fontSize: '0.85rem' }}>
                            <p style={{ margin: '0.25rem 0' }}>
                                <strong>Theta Updated:</strong>
                            </p>
                            <code style={{ display: 'block', background: 'rgba(0,0,0,0.3)', padding: '0.5rem', borderRadius: '4px', wordBreak: 'break-all' }}>
                                {calibrationResult.calibrated_theta.map(v => v.toFixed(6)).join(', ')}
                            </code>
                        </div>
                    </motion.div>
                )}
            </div>
        </motion.div >
    );
}

import React, { useEffect } from 'react';
import { motion } from 'framer-motion';

export default function ControlPanel({ params, handleChange, calculate, loading, setParams }) {

    return (
        <motion.div
            initial={{ x: -20, opacity: 0 }}
            animate={{ x: 0, opacity: 1 }}
            className="glass-panel"
            style={{ padding: '1.5rem', height: 'fit-content' }}
        >
            <form onSubmit={calculate} style={{ display: 'flex', flexDirection: 'column', gap: '1rem' }}>
                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem' }}>
                    <div>
                        <label>Spot Price (S0)</label>
                        <input type="number" name="S0" value={params.S0} onChange={handleChange} step="0.1" required />
                    </div>
                    <div>
                        <label>Strike Price (K)</label>
                        <input type="number" name="K" value={params.K} onChange={handleChange} step="0.1" required />
                    </div>
                </div>

                <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '1rem' }}>
                    <div>
                        <label>Time (T)</label>
                        <input type="number" name="T" value={params.T} onChange={handleChange} step="0.01" required />
                    </div>
                    <div>
                        <label>Risk-free Rate (r)</label>
                        <input type="number" name="r" value={params.r} onChange={handleChange} step="0.01" required />
                    </div>
                </div>

                <div>
                    <label>Steps (N)</label>
                    <input type="number" name="N" value={params.N} onChange={handleChange} max="100" required />
                    <small style={{ display: 'block', color: '#64748b', fontSize: '0.7rem', marginTop: 4 }}>
                        {Number(params.N) > 100
                            ? "Max N is 100"
                            : "P&L Heatmap Visualization"}
                    </small>
                </div>

                <div>
                    <label>Local Volatility (Theta) [comma separated]</label>
                    <input type="text" name="ThetaStr" value={params.ThetaStr} onChange={handleChange} placeholder="0.2, 0.25, 0.22" required />
                </div>

                <div>
                    <label>Time Intervals (tau) [comma separated]</label>
                    <input type="text" name="tauStr" value={params.tauStr} onChange={handleChange} placeholder="0.0, 0.5, 1.0" required />
                </div>

                <div style={{ display: 'flex', gap: '1rem', alignItems: 'center' }}>
                    <div style={{ flex: 1 }}>
                        <label>Option Type</label>
                        <select name="type" value={params.type} onChange={handleChange} style={{ width: '100%' }}>
                            <option value="CALL">Call</option>
                            <option value="PUT">Put</option>
                        </select>
                    </div>
                    <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem', marginTop: '1.2rem' }}>
                        <input
                            type="checkbox"
                            name="isAmerican"
                            checked={params.isAmerican}
                            onChange={(e) => setParams({ ...params, isAmerican: e.target.checked })}
                            style={{ width: 'auto' }}
                        />
                        <span style={{ fontSize: '0.9rem' }}>American</span>
                    </div>
                </div>

                <button type="submit" className="btn-primary" disabled={loading}>
                    {loading ? 'Calculating...' : 'Calculate'}
                </button>
            </form>
        </motion.div>
    );
}

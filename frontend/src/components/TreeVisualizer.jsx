import React, { useMemo } from 'react';
import { motion } from 'framer-motion';
import TreeNodes from './TreeNodes';

export default function TreeVisualizer({ priceTree, valueTree, nSteps }) {

    // Helper: Interpolate color from Red (low) to Green (high)
    // Using cubic root (1/3) to bias gradient towards green VERY early
    const getColor = (val, maxVal) => {
        if (maxVal === 0) return '#ef4444';
        let ratio = Math.max(0, Math.min(1, val / maxVal));

        // Strong non-linear scaling
        ratio = Math.pow(ratio, 0.4);

        const r = Math.round(239 + (34 - 239) * ratio);
        const g = Math.round(68 + (197 - 68) * ratio);
        const b = Math.round(68 + (94 - 68) * ratio);
        return `rgb(${r}, ${g}, ${b})`;
    };

    const visualizationData = useMemo(() => {
        if (!priceTree || !priceTree.length || !valueTree || !valueTree.length) return null;

        const N = Number(nSteps);

        // Find max value for heatmap normalization
        let maxOptionValue = 0;
        valueTree.forEach(level => {
            level.forEach(val => {
                if (val > maxOptionValue) maxOptionValue = val;
            });
        });

        return {
            N,
            maxOptionValue,
            priceTree,
            valueTree
        };
    }, [priceTree, valueTree, nSteps]);

    if (!visualizationData) {
        return (
            <motion.div
                initial={{ opacity: 0 }}
                animate={{ opacity: 1 }}
                className="glass-panel"
                style={{ padding: '1.5rem', minHeight: '300px', display: 'flex', alignItems: 'center', justifyContent: 'center' }}
            >
                <p style={{ textAlign: 'center', color: '#64748b', fontSize: '0.95rem' }}>
                    Click "Calculate" to generate the tree visualization
                </p>
            </motion.div>
        );
    }

    // Render Tree Nodes (Filtered if N > 10)
    return (
        <TreeNodes
            priceTree={priceTree}
            valueTree={valueTree}
            nSteps={visualizationData.N}
            maxOptionValue={visualizationData.maxOptionValue}
            getColor={getColor}
            filterSteps={visualizationData.N > 10}
        />
    );
}

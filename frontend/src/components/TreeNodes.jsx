import React, { useMemo } from 'react';
import { motion } from 'framer-motion';

export default function TreeNodes({ priceTree, valueTree, nSteps, maxOptionValue, filterSteps, getColor }) {

    const treeData = useMemo(() => {
        let displaySteps = [];
        const totalSteps = priceTree.length;

        if (!filterSteps) {
            for (let i = 0; i < totalSteps; i++) displaySteps.push(i);
        } else {
            for (let i = 0; i < totalSteps; i += 10) {
                displaySteps.push(i);
            }
            if (displaySteps[displaySteps.length - 1] !== totalSteps - 1) {
                displaySteps.push(totalSteps - 1);
            }
        }

        const nodes = [];
        const links = [];
        const stepLabels = [];

        const isLargeN = nSteps > 10;
        const nodeRadius = isLargeN ? 12 : 16;
        const fontSize = isLargeN ? 7 : 9;
        const priceSize = isLargeN ? 6 : 7;

        const levelWidth = 140;
        const width = displaySteps.length * levelWidth + 100;

        const lastRealStep = displaySteps[displaySteps.length - 1];
        const maxVerticalNodes = 2 * lastRealStep + 1;
        const spacingY = isLargeN ? 35 : 50;
        const height = Math.max(600, maxVerticalNodes * spacingY + 100);
        const centerY = height / 2;

        displaySteps.forEach((realStepIndex, visualIndex) => {
            const levelNodes = priceTree[realStepIndex];
            const levelValues = valueTree[realStepIndex];

            const x = visualIndex * levelWidth + 50;

            stepLabels.push({ x: x, y: 30, step: realStepIndex });

            if (visualIndex > 0) {
                const prevRealStep = displaySteps[visualIndex - 1];
                if (realStepIndex - prevRealStep > 1) {
                    const midX = x - levelWidth / 2;
                    const dotSpacing = 15;
                    for (let d = 0; d < 3; d++) {
                        nodes.push({
                            id: `separator-${visualIndex}-${d}`,
                            x: midX,
                            y: centerY - dotSpacing + d * dotSpacing,
                            isSeparator: true
                        });
                    }
                }
            }

            for (let k = 0; k < levelNodes.length; k++) {
                const j = k - realStepIndex;
                const y = centerY - j * spacingY;

                const val = levelValues[k];
                const price = levelNodes[k];
                const color = getColor(val, maxOptionValue);

                nodes.push({
                    id: `${realStepIndex}-${k}`,
                    x, y, price, val, color,
                    isSeparator: false
                });

                if (visualIndex < displaySteps.length - 1) {
                    const nextRealStep = displaySteps[visualIndex + 1];
                    if (nextRealStep === realStepIndex + 1) {
                        const nextX = (visualIndex + 1) * levelWidth + 50;
                        const targets = [k, k + 1, k + 2];
                        targets.forEach(targetK => {
                            const targetJ = targetK - nextRealStep;
                            const targetY = centerY - targetJ * spacingY;
                            links.push({ x1: x, y1: y, x2: nextX, y2: targetY });
                        });
                    }
                }
            }
        });

        return { nodes, links, width, height, stepLabels, nodeRadius, fontSize, priceSize };
    }, [filterSteps, nSteps, priceTree, valueTree, maxOptionValue, getColor]);

    const { nodes, links, width, height, stepLabels, nodeRadius, fontSize, priceSize } = treeData;

    return (
        <motion.div
            initial={{ opacity: 0 }}
            animate={{ opacity: 1 }}
            className="glass-panel"
            style={{ padding: '1.5rem', minHeight: '400px', overflowX: 'auto', overflowY: 'hidden' }}
        >
            <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem' }}>
                <h3 style={{ margin: 0 }}>Trinomial Heatmap Tree</h3>
                <div style={{ fontSize: '0.8rem', color: '#94a3b8' }}>
                    Showing steps: {nSteps > 10 ? '0, 10, 20, 30...' : 'All'}
                </div>
            </div>

            <div style={{ overflow: 'auto', textAlign: 'center', maxHeight: '800px' }}>
                <svg width={width} height={height}>
                    {stepLabels.map((label, idx) => (
                        <text
                            key={`label-${idx}`}
                            x={label.x} y={label.y}
                            textAnchor="middle" fontSize="11" fill="#94a3b8" fontWeight="600"
                        >
                            Step {label.step}
                        </text>
                    ))}

                    {links.map((link, idx) => (
                        <line
                            key={idx}
                            x1={link.x1} y1={link.y1} x2={link.x2} y2={link.y2}
                            stroke="#475569" strokeWidth="1" opacity="0.3"
                        />
                    ))}

                    {nodes.map((node) => {
                        if (node.isSeparator) {
                            return <circle key={node.id} cx={node.x} cy={node.y} r={3} fill="#64748b" />;
                        }

                        return (
                            <g key={node.id}>
                                <circle
                                    cx={node.x} cy={node.y} r={nodeRadius}
                                    fill={node.color}
                                    stroke={node.val === maxOptionValue ? "#fff" : "none"}
                                    strokeWidth={node.val === maxOptionValue ? 2 : 0}
                                />
                                <text x={node.x} y={node.y - 3} textAnchor="middle" fontSize={fontSize} fill="white" fontWeight="bold" style={{ textShadow: '0 1px 2px rgba(0,0,0,0.8)' }}>
                                    {node.val.toFixed(2)}
                                </text>
                                <text x={node.x} y={node.y + 7} textAnchor="middle" fontSize={priceSize} fill="rgba(255,255,255,0.8)" style={{ textShadow: '0 1px 2px rgba(0,0,0,0.8)' }}>
                                    {node.price.toFixed(1)}
                                </text>
                            </g>
                        );
                    })}
                </svg>
            </div>

            <div style={{ display: 'flex', gap: '1rem', justifyContent: 'center', marginTop: '1rem', fontSize: '0.8rem', color: '#94a3b8' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                    <span style={{ width: 12, height: 12, background: '#ef4444', borderRadius: '50%', display: 'inline-block' }}></span>
                    <span>Low Value (0)</span>
                </div>
                <div style={{ display: 'flex', alignItems: 'center', gap: '0.5rem' }}>
                    <span style={{ width: 12, height: 12, background: '#22c55e', borderRadius: '50%', display: 'inline-block' }}></span>
                    <span>High Value</span>
                </div>
            </div>
        </motion.div>
    );
}

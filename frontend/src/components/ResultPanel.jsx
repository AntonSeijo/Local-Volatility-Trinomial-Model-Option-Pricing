import React from 'react';
import { motion } from 'framer-motion';

export default function ResultPanel({ result }) {
    return (
        <motion.div
            initial={{ y: 20, opacity: 0 }}
            animate={{ y: 0, opacity: 1 }}
            className="glass-panel"
            style={{ padding: '2rem', textAlign: 'center' }}
        >
            <h2 style={{ margin: 0, color: '#94a3b8', fontSize: '1rem', textTransform: 'uppercase', letterSpacing: '1px' }}>Option Price</h2>
            <div style={{ fontSize: '3rem', fontWeight: 'bold', margin: '0.5rem 0' }}>
                {result ? `$${result.price.toFixed(4)}` : '---'}
            </div>
        </motion.div>
    );
}

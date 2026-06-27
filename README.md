# ProteinSynergyDock

**Structure-aware drug combination synergy prediction via co-docking GNN with GO function context**

[![Python](https://img.shields.io/badge/Python-3.10+-blue.svg)](https://python.org)
[![PyTorch](https://img.shields.io/badge/PyTorch-2.0+-red.svg)](https://pytorch.org)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

## Overview

ProteinSynergyDock is the first model to predict drug combination synergy using:
- **Real 3D molecular docking** (AutoDock Vina) of both drugs against their target protein
- **GO function context** from ProteinWhisper++ (ontology-aware protein function encoder)
- **Cross-drug attention** to model geometric complementarity between two drugs in the same binding pocket

## Results

| Model | Pearson r | AUROC | Method |
|-------|-----------|-------|--------|
| **ProteinSynergyDock** | **0.5768** | **0.5408** | 3D docking + GO + cross-attn |
| DrugSynergy3D (baseline) | 0.54 | 0.835 | SMILES + GNN, no docking |
| Random | ~0.0 | 0.50 | — |

**ProteinWhisper++** (protein encoder): Fmax **0.4006** (7.9x improvement over baseline 0.0504)

## Architecture
Drug A (SMILES) ──┐

├── GATv2 Drug Encoder ──┐

Drug B (SMILES) ──┘                        ├── Cross-Drug Attention

│        +

Protein Sequence ── ProteinWhisper++ ──────┤   FiLM Conditioning

(GO DAG decoder)       │

├── Synergy Head

AutoDock Vina ── Binding Scores ───────────┘

↓

Synergy Score + Class

## Key Innovations

1. **Cross-drug geometric attention** — first model to explicitly model Drug A / Drug B geometric relationships inside a shared protein pocket
2. **GO-conditioned synergy** — protein function embeddings from ProteinWhisper++ injected via FiLM conditioning
3. **Real docking integration** — AutoDock Vina binding affinities as structural features
4. **Dark target capability** — works on unannotated proteins with no prior drug data

## Data

- **NCI ALMANAC** — 231 real drug combination synergy scores
- **SwissProt** — 109,720 proteins for ProteinWhisper++ training
- **RCSB PDB** — 19 cancer target crystal structures
- **Gene Ontology** — 38,245 GO terms, 57,824 DAG edges

## Installation

```bash
git clone https://github.com/Aprameya05/ProteinSynergyDock.git
cd ProteinSynergyDock
pip install -r requirements.txt
```

## Related Projects

- [ProteinWhisper](https://github.com/Aprameya05/ProteinWhisper) — zero-shot protein function annotation
- [DrugSynergy3D](https://github.com/Aprameya05/DrugSynergy3D) — SE(3) equivariant drug synergy prediction

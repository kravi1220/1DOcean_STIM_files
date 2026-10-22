#!/bin/bash
# Regenerate every table and figure from the model output in experiments/ (or below $KG6_ROOT).
# Tables go to results/tables/, figures to results/figures/. Run after all model runs have finished.
# Usage: scripts/run_analysis.sh
set -e
cd "$(dirname "$0")/../analysis"
TAB=$(python3 -c 'import paths; print(paths.TAB)'); G=$(python3 -c 'import paths; print(paths.GOTM)')
mkdir -p "$TAB"
# tables (text); the figure scripts read blockA_score, blockW_compare and blockD_table
python3 blockA_score.py        > "$TAB/blockA_score.txt"         # Section 4.2, Fig. 3
python3 compare_fullyear.py    > "$TAB/compare_fullyear.txt"     # Sections 4.3.3, 4.6
python3 compare_oct15.py       > "$TAB/compare_oct15.txt"        # Section 4.3.2
python3 blockW_compare.py      > "$TAB/blockW_compare.txt"       # Section 4.3.1, Fig. 7
python3 blockD_table.py        > "$TAB/blockD_table.txt"         # Section 4.5, Fig. 9
python3 adc_dt_ladder.py       > "$TAB/adc_dt_ladder.txt"        # Section 4.2.3
python3 les_convergence.py     > "$TAB/les_convergence.txt"      # Section 3.5, Appendix A
python3 blockR_metrics.py $G/blockR/dt*/R* $G/blockR/R5_dt112.5/m? > "$TAB/blockR_metrics.txt"   # Section 4.1
python3 icepack_salt_budget.py $G/blockC/C1_bl99 $G/blockC/C2_mushy > "$TAB/icepack_salt_budget.txt"   # Appendix A
# figures
python3 fig_map.py             # Fig. 1  (K0)
python3 fig_main.py            # Figs. 2-6 (K1-K5)
python3 fig_summary.py         # Figs. 7-11 (K7-K11)

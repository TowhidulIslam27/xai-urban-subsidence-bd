"""Section 5.1, Table 4  LaTeX body of the screening table, from analysis/screening/screen_<c>.csv.

Kept candidates show Spearman's rho with vU at PS cells; dropped ones show d (reason):
rho = |rho| > 0.8 with a factor of higher priority, VIF, const. = one value in > 95 % of cells,
sparse = a value in < 5 % of cells, coord. = acts as a coordinate, res. = native resolution > 1 km. A last block lists the check-only factors.
Two blocks of factors side by side (6 columns: factor, Dhaka, Chattogram, twice).
Output    analysis/screening/table4_screening.tex (paste between \\midrule and \\bottomrule)
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import pandas as pd
from factor_groups import PRIORITY, GROUP, LC_CHECK, PS_ONLY
from figlib import LABEL
from paths import SCREEN

S = {c: pd.read_csv(SCREEN / f'screen_{c}.csv', index_col=0) for c in ['dha', 'ctg']}
TAG = [('near', 'const.'), ('too', 'sparse'), ('VIF', 'VIF'), ('acts', 'coord.'), ('native', 'res.')]


def num(v):
    return ('$+$' if v >= 0 else '$-$') + f'{abs(v):.2f}'


def cell(c, k):
    if k not in S[c].index:
        return '--'
    r = S[c].loc[k]
    if r.status == 'dropped':
        return 'd (' + next((t for p, t in TAG if r.reason.startswith(p)), '$\\rho$') + ')'
    return num(r.rho_vU)


def lab(k):
    return (LABEL.get(k, k).replace('%', '\\%').replace('km²', 'km$^2$').replace('m³', 'm$^3$')
            .replace('°', '$^\\circ$').replace('²', '$^2$'))


rows, last = [], None                      # (kind, label, dha, ctg); kind 'g' = group header
from figlib import display_order
for k in display_order([k for k in PRIORITY if k in S['dha'].index or k in S['ctg'].index]):
    if GROUP[k] != last:
        rows.append(('g', GROUP[k], '', '')); last = GROUP[k]
    rows.append(('f', lab(k), cell('dha', k), cell('ctg', k)))
rows.append(('g', 'Checks only (partial coverage or PS cells only)', '', ''))
for k in LC_CHECK + PS_ONLY:
    if k in S['dha'].index and k in S['ctg'].index:
        rows.append(('f', lab(k) + (' (PS)' if k in PS_ONLY else ''), num(S['dha'].loc[k].rho_vU), num(S['ctg'].loc[k].rho_vU)))

# two blocks side by side, split at the group header closest to the middle
heads = [i for i, r in enumerate(rows) if r[0] == 'g']
cut = min(heads[1:], key=lambda i: abs(i - len(rows) / 2))
L, R = rows[:cut], rows[cut:]


def fmt(r):
    if r is None:
        return ' & & '
    if r[0] == 'g':
        return f'\\multicolumn{{3}}{{@{{}}l}}{{\\textit{{{r[1]}}}}}'
    return f'\\quad {r[1]} & {r[2]} & {r[3]}'


out = [f'{fmt(l)} & {fmt(r)} \\\\' for l, r in
       ((L[i] if i < len(L) else None, R[i] if i < len(R) else None) for i in range(max(len(L), len(R))))]
open(SCREEN / 'table4_screening.tex', 'w').write('\n'.join(out) + '\n')
print('written', SCREEN / 'table4_screening.tex', len(out), 'lines (two blocks of', len(L), 'and', len(R), 'rows)')

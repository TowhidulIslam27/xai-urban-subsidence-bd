"""Section 3.2 / 5.1, Table 4  Screening of the candidate factors (PS cells, per city).

Reads analysis/factors/<c>_factors_full.csv (after the GEE merge) or <c>_factors.csv.
Candidates are the factors of common/factor_groups.PRIORITY present in the table; the coverage-
limited land-cover factors (LC_CHECK) and the PS-only building factors (PS_ONLY) are reported
with status 'check only' and never enter the models.
1. native resolution coarser than 1 km             -> dropped (a few values per city = location code)
2. near-constant: one value in more than 95 % of PS cells, or < 3 distinct values -> dropped;
   too sparse: a value in fewer than 5 % of PS cells (e.g. fill age, defined only where water
   was lost) -> dropped
   (not IQR = 0: zero-inflated factors such as the share of a cell that lost its water are 0 in
   most cells but carry the signal in the rest)
3. acts as a coordinate: R2 > 0.7 against a quadratic surface in E, N -> dropped
   (the rule that already excluded distance to faults, now applied to every factor)
4. Spearman |rho| > 0.8: keep the factor earlier in PRIORITY, drop the other
5. iterative VIF > 10 on the standardised survivors (complete cases), following PRIORITY as in
   step 4: take the factor with the largest VIF, regress it on the others, and drop its main
   contributor that comes later in PRIORITY (|standardised coefficient| >= 0.25; the latest such
   one); only if it has none is the factor itself dropped. So elevation is not lost to TPI 1 km
   (elevation minus its 1 km mean) or to HAND, water loss is kept over the 1990 water share it is
   built from, and the closure of the Dynamic World probabilities (they sum to 1) costs the least
   useful classes rather than built-up. Factors that are empty by
   design in part of a city (fill age where nothing was filled, height where nothing is built;
   < 90 % of PS cells) are left out of the VIF so they do not shrink the sample; their overlap
   with the others is still covered by the Spearman step.
vU is never used to choose factors; rho_vU is reported only.
Writes analysis/screening/screen_<c>.csv and kept_factors.csv (with the transfer set: kept in
both cities and meaning the same thing in both).
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import numpy as np, pandas as pd
from scipy.stats import spearmanr

from paths import FACTORS, SCREEN
from factor_groups import (PRIORITY, GROUP, SET, RES, TRANSFER, LC_CHECK, PS_ONLY, WELL_CHECK, MAX_RES, TREND_R2_MAX,
                           present)
D, O = str(FACTORS) + os.sep, str(SCREEN) + os.sep


def table(c):
    f = D + f'{c}_factors_full.csv'
    return pd.read_csv(f if os.path.exists(f) else D + f'{c}_factors.csv', low_memory=False)


def vif(X):
    Z = (X - X.mean(0)) / X.std(0)
    C = np.corrcoef(Z, rowvar=False)
    return np.diag(np.linalg.pinv(C))


def trend_r2(p, k):
    """Share of a factor's variance explained by a quadratic surface in E, N (PS cells)."""
    e = (p.E - p.E.mean()) / p.E.std(); n = (p.N - p.N.mean()) / p.N.std()
    y = p[k].to_numpy(float); ok = np.isfinite(y)
    if ok.sum() < 10 or np.nanstd(y) == 0:
        return np.nan
    X = np.c_[np.ones(len(p)), e, n, e * e, n * n, e * n][ok]
    b = np.linalg.lstsq(X, y[ok], rcond=None)[0]
    return 1 - np.var(y[ok] - X @ b) / np.var(y[ok])


def screen(c):
    df = table(c)
    p = df[df.has_ps.astype(bool)].copy()
    cand = present(df.columns)
    checks = present(df.columns, LC_CHECK + PS_ONLY + WELL_CHECK)
    allk = cand + checks
    rec = pd.DataFrame({'factor': allk, 'group': [GROUP[k] for k in allk],
                        'set': [SET.get(k, '-') for k in allk], 'res_m': [RES.get(k, np.nan) for k in allk],
                        'transfer': [TRANSFER.get(k, False) for k in allk]}).set_index('factor')
    rec['n'] = p[allk].notna().sum()
    rec['rho_vU'] = [spearmanr(p[k], p.vU, nan_policy='omit')[0] for k in allk]
    rec['trend_r2'] = [trend_r2(p, k) for k in allk]
    rec['status'] = 'kept'
    rec['reason'] = ''
    rec.loc[checks, ['status', 'reason']] = 'check only', 'partial coverage or PS cells only'
    rec.loc[[k for k in checks if k in WELL_CHECK], 'reason'] = 'borehole survey network, not pumping (Dhaka only)'

    def drop(k, why):
        rec.loc[k, ['status', 'reason']] = 'dropped', why

    alive = lambda: [k for k in cand if rec.loc[k, 'status'] == 'kept']
    for k in cand:
        if RES[k] > MAX_RES:
            drop(k, f'native resolution {RES[k] / 1000:g} km')
    for k in alive():
        q = p[k].dropna()
        if len(q) < 0.05 * len(p):
            drop(k, f'too sparse (value in {len(q) / len(p):.1%} of cells)'); continue
        top = q.value_counts(normalize=True).iloc[0] if len(q) else 1.0
        rec.loc[k, 'mode_share'] = top
        if top > 0.95 or q.nunique() < 3:
            drop(k, f'near-constant ({top:.0%} of cells at one value)')
    for k in alive():
        if rec.loc[k, 'trend_r2'] > TREND_R2_MAX:
            drop(k, f'acts as a coordinate (trend R2 = {rec.loc[k, "trend_r2"]:.2f})')
    a = alive()
    R = p[a].rank().corr()                     # Spearman via ranks, pairwise complete
    for i, k1 in enumerate(a):
        if rec.loc[k1, 'status'] != 'kept':
            continue
        for k2 in a[i + 1:]:
            if rec.loc[k2, 'status'] == 'kept' and abs(R.loc[k1, k2]) > 0.8:
                drop(k2, f'|rho| = {abs(R.loc[k1, k2]):.2f} with {k1}')
    a = [k for k in alive() if rec.loc[k, 'n'] >= 0.9 * len(p)]
    for k in alive():
        if k not in a:
            rec.loc[k, 'reason'] = 'VIF not computed (empty by design in part of the city)'
    X = p[a].dropna()
    rec['vif'] = np.nan
    while True:
        v = pd.Series(vif(X[a].to_numpy(float)), index=a)
        rec.loc[a, 'vif'] = v
        if v.max() <= 10:
            break
        h = v.idxmax()
        Z = (X[a] - X[a].mean()) / X[a].std()
        others = [f for f in a if f != h]
        beta = pd.Series(np.linalg.lstsq(np.c_[np.ones(len(Z)), Z[others]], Z[h], rcond=None)[0][1:], index=others)
        later = [f for f in others if a.index(f) > a.index(h) and abs(beta[f]) >= 0.25]
        k = later[-1] if later else h
        drop(k, f'VIF = {v[h]:.1f} for {h}' if k != h else f'VIF = {v[k]:.1f}')
        a.remove(k)
    rec['n_complete'] = len(X)
    rec.to_csv(O + f'screen_{c}.csv')
    return rec


if __name__ == '__main__':
    out = {}
    for c in ['dha', 'ctg']:
        rec = screen(c)
        out[c] = rec
        cand = rec[rec.status != 'check only']
        print(c, 'kept', (cand.status == 'kept').sum(), 'of', len(cand), 'candidates | complete cases',
              rec.n_complete.iloc[0], 'of', int(rec.n.max()))
        print(cand[cand.status != 'kept'][['reason']].to_string())
    idx = list(dict.fromkeys(list(out['dha'].index) + list(out['ctg'].index)))
    k = pd.DataFrame({c: out[c].status.reindex(idx).fillna('not available') for c in out})
    k['group'] = [GROUP[i] for i in idx]
    k['set'] = [SET.get(i, '-') for i in idx]
    k['transfer'] = [TRANSFER.get(i, False) for i in idx]
    k['kept_both'] = (k.dha == 'kept') & (k.ctg == 'kept')
    k['transfer_set'] = k.kept_both & k.transfer
    k.to_csv(O + 'kept_factors.csv')
    print(k.to_string())
    print('transfer set:', list(k.index[k.transfer_set]))

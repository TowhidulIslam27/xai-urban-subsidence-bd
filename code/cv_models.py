"""Section 4.2 / 5.2, Tables 5-6, Fig. 4  RQ1: how much of vU do the nested factor sets explain?

Data      analysis/factors/<c>_factors_full.csv, PS cells with vU; factors kept by the screening
          (analysis/screening/screen_<c>.csv, status 'kept').
Sets      nested: M0 terrain + hydrology + soil; M1 = M0 + geology (contact distance and the shared
          geology class, one-hot); M2 = M1 + urban load + land cover (incl. the city-specific
          objects); M3 = M2 + fill history + groundwater proxies. Measured groundwater heads are not
          in yet: this M3 is 'M3a'; Dhaka gets an M3b once the kriged heads exist.
Models    gbm   gradient-boosted trees: LightGBM when installed, otherwise scikit-learn's
                HistGradientBoostingRegressor (same histogram algorithm, leaf-wise growth, NaN
                handled natively); fixed settings, no early stopping (no inner random split)
          rf    random forest (200 trees, max_features 1/3, min leaf 5, 30 % subsample per tree)
          ols   linear regression on standardised factors (NaN -> median + missing flag)
          xy    gbm on the cell coordinates only: what spatial position alone predicts under the
                same CV (reference line, not a candidate model)
CV        spatial block CV, 5 folds x 10 repeats. Block size 2 km: the residual semivariogram of a
          1 km-block run has an effective range of ~2.0 km in both cities (range_check.csv), so
          1 km blocks would let neighbouring cells leak. Blocks are shuffled and dealt to the fold
          with the fewest cells. Same folds for every set and model in a repeat (paired).
          Sensitivity: gbm-M3 with 1 km and 4 km blocks.
Metrics   per repeat, on the pooled out-of-fold predictions: R2, RMSE, MAE, Spearman rho.
          Gain of a step = paired difference of R2 between consecutive sets in the same repeat.
Output    analysis/models/cv_<c>.csv (one row per city, set, model, repeat), cv_summary.csv,
          step_gain.csv, block_sensitivity.csv, range_check.csv, oof_<c>.csv (gbm, every set,
          repeat 0, with fold ids), table5_cv.tex, table6_gain.tex; Figures/Fig04_cv_skill.png
usage     python cv_models.py [--quick]   (--quick: 2 repeats, for testing)
"""
import os, sys, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "common"))
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from scipy.spatial.distance import pdist
from scipy.optimize import curve_fit
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression

from paths import ROOT, FACTORS, SCREEN, FIGS, CITIES
from factor_groups import SET

OUT = ROOT / 'analysis' / 'models'
OUT.mkdir(parents=True, exist_ok=True)
QUICK = '--quick' in sys.argv
N_REP, N_FOLD, BLOCK_CELLS = (2 if QUICK else 10), 5, 20          # 20 cells = 2 km
SETS = ['M0', 'M1', 'M2', 'M3']
SET_LABEL = {'M0': 'M0 terrain, hydrology, soil', 'M1': 'M1 + geology', 'M2': 'M2 + urban load, land cover',
             'M3': 'M3a + fill history, GW proxies'}
MODELS = ['gbm', 'rf', 'ols']
SEED = 20261020

try:
    import lightgbm as lgb
    ENGINE = 'LightGBM ' + lgb.__version__

    def gbm():
        return lgb.LGBMRegressor(n_estimators=500, learning_rate=0.05, num_leaves=31, min_child_samples=40,
                                 subsample=0.8, subsample_freq=1, colsample_bytree=0.8, reg_lambda=1.0,
                                 n_jobs=2, verbose=-1, random_state=0)
except ImportError:
    from sklearn.ensemble import HistGradientBoostingRegressor
    import sklearn
    ENGINE = 'scikit-learn HistGradientBoosting ' + sklearn.__version__

    def gbm():
        return HistGradientBoostingRegressor(max_iter=500, learning_rate=0.05, max_leaf_nodes=31,
                                             min_samples_leaf=40, l2_regularization=1.0, max_features=0.8,
                                             early_stopping=False, random_state=0)


class OLS:
    """Linear baseline: median-impute (+ missing flag where a factor has gaps), standardise."""
    def fit(self, X, y):
        self.med = np.nanmedian(X, 0)
        self.flag = np.isnan(X).any(0)
        Z = self._z(X, fit=True)
        self.m = LinearRegression().fit(Z, y)
        return self

    def _z(self, X, fit=False):
        Z = np.where(np.isnan(X), self.med, X)
        Z = np.c_[Z, np.isnan(X[:, self.flag]).astype(float)]
        if fit:
            self.mu, self.sd = Z.mean(0), Z.std(0)
            self.sd[self.sd == 0] = 1
        return (Z - self.mu) / self.sd

    def predict(self, X):
        return self.m.predict(self._z(X))


def make(name):
    if name in ('gbm', 'xy'):
        return gbm()
    if name == 'rf':
        return RandomForestRegressor(n_estimators=200, max_features=1 / 3, min_samples_leaf=5, max_samples=0.3,
                                     n_jobs=2, random_state=0)
    return OLS()


def fit_predict(name, X, y, tr, te):
    m = make(name)
    Xt = X[tr]
    if name == 'rf':                       # RF in scikit-learn takes NaN, but impute for a clean baseline
        med = np.nanmedian(Xt, 0)
        Xt, Xe = np.where(np.isnan(Xt), med, Xt), np.where(np.isnan(X[te]), med, X[te])
        return m.fit(Xt, y[tr]).predict(Xe)
    return m.fit(Xt, y[tr]).predict(X[te])


def load(c):
    df = pd.read_csv(FACTORS / f'{c}_factors_full.csv', low_memory=False)
    p = df[df.has_ps.astype(bool) & df.vU.notna()].reset_index(drop=True)
    s = pd.read_csv(SCREEN / f'screen_{c}.csv', index_col=0)
    kept = [k for k in s.index[s.status == 'kept']]
    geo = pd.get_dummies(p.geo_shared, prefix='geo', dtype=float)
    geo = geo.loc[:, geo.sum() >= 30]                    # classes with >= 30 PS cells
    p = pd.concat([p, geo], axis=1)
    cols = {}
    for st in SETS:
        add = [k for k in kept if SET.get(k) == st] + (list(geo.columns) if st == 'M1' else [])
        cols[st] = (cols[SETS[SETS.index(st) - 1]] if st != 'M0' else []) + add
    return p, cols


def folds(p, block_cells, rep):
    """Shuffle the blocks and deal them to the fold with the fewest cells (balanced, random)."""
    b = (p.row // block_cells).astype(str) + '_' + (p.col // block_cells).astype(str)
    cnt = b.value_counts()
    order = cnt.index.to_numpy()
    rng = np.random.default_rng(SEED + rep)
    order = order[rng.permutation(len(order))]
    load_, f_of = np.zeros(N_FOLD), {}
    for blk in order:
        k = int(np.argmin(load_)); f_of[blk] = k; load_[k] += cnt[blk]
    return b.map(f_of).to_numpy()


def scores(y, yh):
    r = y - yh
    return dict(r2=1 - (r ** 2).sum() / ((y - y.mean()) ** 2).sum(), rmse=np.sqrt((r ** 2).mean()),
                mae=np.abs(r).mean(), rho=spearmanr(y, yh)[0])


def cv_run(p, X, y, name, block_cells, rep):
    f = folds(p, block_cells, rep)
    oof = np.full(len(y), np.nan)
    for k in range(N_FOLD):
        tr, te = np.where(f != k)[0], np.where(f == k)[0]
        oof[te] = fit_predict(name, X, y, tr, te)
    return oof, f


def vario_range(xy, z, rng):
    i = rng.choice(len(z), min(5000, len(z)), replace=False)
    d, g = pdist(xy[i]), 0.5 * pdist(z[i, None]) ** 2
    b = np.linspace(0, 6000, 25); k = np.digitize(d, b) - 1; m = (k >= 0) & (k < 24)
    h = (b[:-1] + b[1:]) / 2; gm = np.array([g[m & (k == j)].mean() for j in range(24)])
    f = lambda h, n, s, r: n + s * (1 - np.exp(-3 * h / r))
    (n, s, r), _ = curve_fit(f, h, gm, p0=[gm[0], gm[-1], 1500], bounds=([0, 0, 50], [np.inf, np.inf, 50000]))
    return n, s, r


def run_city(c):
    t0 = time.time()
    p, cols = load(c)
    y = p.vU.to_numpy(float)
    print(c, len(p), 'PS cells |', ', '.join(f'{s} {len(cols[s])}' for s in SETS), '|', ENGINE, flush=True)
    rows, oof_tab = [], p[['cell_id', 'row', 'col', 'E', 'N', 'vU']].copy()
    for rep in range(N_REP):
        for st in SETS:
            X = p[cols[st]].to_numpy(float)
            for mdl in MODELS:
                yh, f = cv_run(p, X, y, mdl, BLOCK_CELLS, rep)
                rows.append(dict(city=c, set=st, model=mdl, rep=rep, n_factors=len(cols[st]), **scores(y, yh)))
                if rep == 0 and mdl == 'gbm':
                    oof_tab[f'oof_{st}'] = yh
                    oof_tab['fold'] = f
        yh, _ = cv_run(p, p[['E', 'N']].to_numpy(float), y, 'xy', BLOCK_CELLS, rep)
        rows.append(dict(city=c, set='XY', model='xy', rep=rep, n_factors=2, **scores(y, yh)))
        print(f'  {c} repeat {rep + 1}/{N_REP}  {time.time() - t0:5.0f}s', flush=True)
    cv = pd.DataFrame(rows)
    cv.to_csv(OUT / f'cv_{c}.csv', index=False)
    oof_tab.to_csv(OUT / f'oof_{c}.csv', index=False)

    # block-size sensitivity (gbm, M3)
    X = p[cols['M3']].to_numpy(float)
    sens = []
    for bc in (10, 20, 40):
        for rep in range(N_REP):
            yh, _ = cv_run(p, X, y, 'gbm', bc, rep)
            sens.append(dict(city=c, block_km=bc / 10, rep=rep, **scores(y, yh)))
    # residual range of the 1 km and 2 km runs (repeat 0)
    rng = np.random.default_rng(SEED)
    xy = p[['E', 'N']].to_numpy(float)
    rc = []
    for bc in (10, 20):
        yh, _ = cv_run(p, X, y, 'gbm', bc, 0)
        n, s, r = vario_range(xy, y - yh, rng)
        rc.append(dict(city=c, block_km=bc / 10, nugget=n, partial_sill=s, eff_range_m=r))
    n, s, r = vario_range(xy, y, rng)
    rc.append(dict(city=c, block_km=np.nan, nugget=n, partial_sill=s, eff_range_m=r, note='vU itself'))
    print(f'  {c} done {time.time() - t0:5.0f}s', flush=True)
    return cv, pd.DataFrame(sens), pd.DataFrame(rc), cols


def summarise(cv, sens):
    g = cv.groupby(['city', 'set', 'model'])
    summ = g[['r2', 'rmse', 'mae', 'rho']].agg(['mean', 'std'])
    summ.columns = [f'{a}_{b}' for a, b in summ.columns]
    summ['n_factors'] = g.n_factors.first()
    summ = summ.reset_index()
    summ.to_csv(OUT / 'cv_summary.csv', index=False)
    gain = []
    for (c, mdl), d in cv[cv.model != 'xy'].groupby(['city', 'model']):
        w = d.pivot(index='rep', columns='set', values='r2')
        for a, b in zip(SETS[:-1], SETS[1:]):
            dd = w[b] - w[a]
            gain.append(dict(city=c, model=mdl, step=f'{a}->{b}', added=SET_LABEL[b].split(' + ', 1)[-1],
                             d_r2_mean=dd.mean(), d_r2_lo=dd.quantile(0.025), d_r2_hi=dd.quantile(0.975),
                             share_reps_positive=(dd > 0).mean()))
    gain = pd.DataFrame(gain)
    gain.to_csv(OUT / 'step_gain.csv', index=False)
    s = sens.groupby(['city', 'block_km'])[['r2', 'rmse']].agg(['mean', 'std'])
    s.columns = [f'{a}_{b}' for a, b in s.columns]
    s.reset_index().to_csv(OUT / 'block_sensitivity.csv', index=False)
    return summ, gain


CNAME = {'dha': 'Dhaka', 'ctg': 'Chattogram'}
MNAME = {'gbm': 'GBM', 'rf': 'RF', 'ols': 'OLS'}


def tex(summ, gain):
    L = [r'% Table 5 body: spatial block CV (2 km blocks, 5 folds x 10 repeats), mean (s.d.) over repeats',
         r'Set & Model & \multicolumn{3}{c}{Dhaka} & \multicolumn{3}{c}{Chattogram} \\',
         r' & & $R^2$ & RMSE & $\rho$ & $R^2$ & RMSE & $\rho$ \\ \midrule']
    for st in SETS + ['XY']:
        for mdl in (MODELS if st != 'XY' else ['xy']):
            cells = []
            for c in CITIES:
                r = summ[(summ.city == c) & (summ.set == st) & (summ.model == mdl)]
                if len(r):
                    r = r.iloc[0]
                    cells += [f'{r.r2_mean:.2f} ({r.r2_std:.2f})', f'{r.rmse_mean:.2f}', f'{r.rho_mean:.2f}']
                else:
                    cells += ['--'] * 3
            lab = st if st != 'XY' else 'Coordinates only'
            L.append(f'{lab} & {MNAME.get(mdl, "GBM")} & ' + ' & '.join(cells) + r' \\')
        L.append(r'\addlinespace')
    open(OUT / 'table5_cv.tex', 'w').write('\n'.join(L) + '\n')
    L = [r'% Table 6 body: gain in CV R2 per step, paired over repeats: mean [2.5, 97.5 %]',
         r'Step & Added & Model & Dhaka & Chattogram \\ \midrule']
    for a, b in zip(SETS[:-1], SETS[1:]):
        for mdl in MODELS:
            cells = []
            for c in CITIES:
                r = gain[(gain.city == c) & (gain.model == mdl) & (gain.step == f'{a}->{b}')].iloc[0]
                cells.append(f'{r.d_r2_mean:+.3f} [{r.d_r2_lo:+.3f}, {r.d_r2_hi:+.3f}]')
            L.append(f'{a}$\\to${b} & {SET_LABEL[b].split(" + ", 1)[-1]} & {MNAME[mdl]} & ' + ' & '.join(cells) + r' \\')
        L.append(r'\addlinespace')
    open(OUT / 'table6_gain.tex', 'w').write('\n'.join(L) + '\n')


def figure(summ):
    import matplotlib.pyplot as plt
    from figstyle import set_style, save, panel_label, W2, CITY_NAME
    set_style()
    fig, axs = plt.subplots(1, 2, figsize=(W2, 62 / 25.4), sharey=True)
    col = {'gbm': '#0072B2', 'rf': '#009E73', 'ols': '#7f7f7f'}
    mk = {'gbm': 'o', 'rf': 's', 'ols': '^'}
    x = np.arange(len(SETS))
    for i, (ax, c) in enumerate(zip(axs, CITIES)):
        s = summ[summ.city == c]
        for j, mdl in enumerate(MODELS):
            d = s[s.model == mdl].set_index('set').loc[SETS]
            ax.errorbar(x + (j - 1) * 0.06, d.r2_mean, yerr=d.r2_std, color=col[mdl], marker=mk[mdl], ms=3.2,
                        lw=0.9, capsize=1.5, label=MNAME[mdl])
        xy = s[s.model == 'xy'].iloc[0]
        ax.axhline(xy.r2_mean, color='k', ls='--', lw=0.6)
        ax.text(len(SETS) - 0.65, xy.r2_mean, 'coordinates only', ha='right', va='bottom', fontsize=6)
        ax.set_xticks(x, ['M0\nterrain,\nhydrology, soil', 'M1\n+ geology', 'M2\n+ urban load,\nland cover',
                          'M3a\n+ fill, GW\nproxies'])
        ax.set_xlim(-0.35, len(SETS) - 0.65)
        ax.set_title(CITY_NAME[c])
        ax.grid(axis='y', lw=0.3, color='0.85')
        panel_label(ax, '(' + 'ab'[i] + ')')
        if i == 0:
            ax.set_ylabel('Cross-validated $R^2$ of $v_U$ (2 km spatial blocks)')
        else:
            ax.legend(loc='lower right')
    save(fig, str(FIGS), 'Fig04_cv_skill', formats=('png',))


if __name__ == '__main__':
    cvs, senss, rcs = [], [], []
    for c in CITIES:
        cv, sens, rc, cols = run_city(c)
        cvs.append(cv); senss.append(sens); rcs.append(rc)
        pd.Series({s: ';'.join(v) for s, v in cols.items()}).to_csv(OUT / f'sets_{c}.csv', header=['factors'])
    cv, sens = pd.concat(cvs), pd.concat(senss)
    pd.concat(rcs).to_csv(OUT / 'range_check.csv', index=False)
    summ, gain = summarise(cv, sens)
    tex(summ, gain)
    figure(summ)
    open(OUT / 'engine.txt', 'w').write(ENGINE + '\n')
    pd.set_option('display.width', 200)
    print(summ[['city', 'set', 'model', 'n_factors', 'r2_mean', 'r2_std', 'rmse_mean', 'rho_mean']].round(3).to_string(index=False))
    print(gain.round(3).to_string(index=False))
    print(pd.read_csv(OUT / 'block_sensitivity.csv').round(3).to_string(index=False))
    print(pd.read_csv(OUT / 'range_check.csv').round(1).to_string(index=False))

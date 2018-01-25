from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from sklearn import mixture


def _fit_mixture_models(s, max_gaussians=5):
    """
    Fit mixture models to a numeric Series column.

    :param s:
    :param max_gaussians:
    :return:
    """
    x = s.values.reshape(-1, 1)

    # Fit GMM with a range of mixture components
    models = [mixture.GaussianMixture(n, n_init=3) for n in range(1, max_gaussians+1)]
    models = [m.fit(x) for m in models]
    return {'models': models, 'data': x}


def _score_models(models, x):
    """Return AIC and BIC for list of GMMs.

    :param models:
    :return:
    """
    # Compute AIC and BIC
    AIC = [m.aic(x) for m in models]
    BIC = [m.bic(x) for m in models]

    return {'AIC': AIC, 'BIC': BIC}


def _pick_best(models, x):
    """Pick best model from list of GMMs.

    :param models:
    :return:
    """
    # Pick best
    aic = _score_models(models, x)['AIC']
    best = np.argmin(aic)
    print('Best model has {} components'.format(best + 1))
    return models[best]


def _gmm_plot(models, x):
    """Plot some diagnostics for the GMM fitting.

    :param models:
    :return:
    """
    # Adapted from: http://www.astroml.org/book_figures/chapter4/fig_GMM_1D.html
    # ------------------------------------------------------------
    # Plot the results
    #  We'll use three panels:
    #   1) data + best-fit mixture
    #   2) AIC and BIC vs number of components
    #   3) probability that a point came from each component

    fig = plt.figure(figsize=(10, 5))
    fig.subplots_adjust(left=0.12, right=0.97,
                        bottom=0.21, top=0.9, wspace=0.5)

    # plot 1: data + best-fit mixture
    best_model = _pick_best(models, x)

    ax = fig.add_subplot(121)

    grid = np.linspace(1, max(x), 1000)  # Can't go from 0 because get p > 1 !!!
    grid = grid.reshape(-1, 1)

    logprob = best_model.score_samples(grid)
    responsibilities = best_model.predict_proba(grid)
    pdf = np.exp(logprob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]

    ax.hist(x, 30, normed=True, histtype='stepfilled', alpha=0.4)
    # ax.plot(x, pdf, '-k')
    # ax.plot(x, pdf_individual, '--k')
    ax.text(0.04, 0.96, "Best-fit Mixture",
            ha='left', va='top', transform=ax.transAxes)
    ax.set_xlabel('$x$')
    ax.set_ylabel('$p(x)$')
    ax.plot(grid, pdf_individual, '--k')

    # plot 2: AIC and BIC
    scores = _score_models(models, x)

    ax = fig.add_subplot(122)
    ax.plot(list(range(1, len(scores['AIC']) + 1)), scores['AIC'], '-k', label='AIC')
    ax.plot(list(range(1, len(scores['BIC']) + 1)), scores['BIC'], '--k', label='BIC')
    ax.set_xlabel('n. components')
    ax.set_ylabel('information criterion')
    ax.legend(loc=2)

    fig.suptitle('Residue Occupancy GMM Diagnostics')

    return None


def _core_column_mask(model, x, n_groups=1):
    """Return a row mask indicating membership of most occupied column set.

    :param model:
    :param x:
    :return:
    """
    return pd.Series(model.predict(x)).isin(np.argsort(model.means_.reshape(-1))[-n_groups:]).values


import pandas as pd
from scipy import stats


def _comparative_regression(column_variant_counts, filter_mask=None, regressor ='Human_res_occupancy'):
    """Return table of regression parameters for missense / synonymous counts vs. a regressor.

    :param column_variant_counts:
    :param filter_mask:
    :param regressor:
    :return:
    """
    regressions = []
    # Regressions for all columns
    regressions.append(stats.linregress(x=column_variant_counts[regressor],
                                        y=column_variant_counts['missense_variant'])._asdict())
    regressions.append(stats.linregress(x=column_variant_counts[regressor],
                                        y=column_variant_counts['synonymous_variant'])._asdict())
    row_names = ['missense', 'synonymous']
    if filter_mask is not None:
        # Regressions for filtered columns
        regressions.append(stats.linregress(x=column_variant_counts[filter_mask][regressor],
                                            y=column_variant_counts[filter_mask]['missense_variant'])._asdict())
        regressions.append(stats.linregress(x=column_variant_counts[filter_mask][regressor],
                                            y=column_variant_counts[filter_mask]['synonymous_variant'])._asdict())
        row_names += ['filtered_missense', 'filtered_synonymous']

    results = pd.DataFrame(regressions, index=row_names)
    results.columns.name = regressor
    return results


def _column_variant_scores(column_variant_counts, variant_class='missense_variant'):
    """Calculate variation scores from a DataFrame of column variant totals.

    :param column_variant_counts:
    :param variant_class:
    :return:
    """
    # column_data.query('Human_res_occupancy >= 5', inplace=True)
    alignment_totals = column_variant_counts.loc[:, [variant_class, 'Human_res_occupancy']].sum()

    # Calculate missense scores
    missense_scores = []
    for mis_col, occ_col in column_variant_counts[[variant_class, 'Human_res_occupancy']].itertuples(index=False):
        mis_other = alignment_totals[variant_class] - mis_col
        occ_other = alignment_totals['Human_res_occupancy'] - occ_col
        oddsratio, pvalue = stats.fisher_exact([[mis_col, mis_other], [occ_col, occ_other]])
        missense_scores.append((oddsratio, pvalue))

    # Parse to dataframe
    return pd.DataFrame(missense_scores, columns=['oddsratio', 'pvalue'], index=column_variant_counts.index)


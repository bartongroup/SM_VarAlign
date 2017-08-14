import pandas as pd
import seaborn as sns
from scipy import stats


def _aggregate_annotation(aligned_variant_table, annotation_column, aggregate_by=[('Alignment', 'Column')],
                          fill_value=0):
    """Count category frequencies over a set of aggregation columns.

    Create frequency table for a categorical column over other columns and/or index levels. Index names can only be used
    in Pandas >= 0.20.1.

    :param aligned_variant_table: DataFrame
    :param annotation_column: Column name containing categories to count.
    :param aggregate_by: List of column and/or index names to group by.
    :param fill_value: Value to fill unobserved categories.
    :return: DataFrame of counts indexed by `aggregate_by` with the unique values of `annotation_column` as columns.
    """
    groupby = aggregate_by + [annotation_column]
    grouped = aligned_variant_table.groupby(groupby)
    return grouped.size().unstack(fill_value=fill_value).rename_axis('', 1)


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


def _column_variant_scores(column_variant_counts, variant_class='missense_variant', occupancy='Human_res_occupancy'):
    """Calculate variation scores from a DataFrame of column variant totals.

    :param column_variant_counts:
    :param variant_class:
    :return:
    """
    # column_data.query('Human_res_occupancy >= 5', inplace=True)
    alignment_totals = column_variant_counts.loc[:, [variant_class, occupancy]].sum()

    # Calculate missense scores
    missense_scores = []
    for mis_col, occ_col in column_variant_counts[[variant_class, occupancy]].itertuples(index=False):
        mis_other = alignment_totals[variant_class] - mis_col
        occ_other = alignment_totals[occupancy] - occ_col
        oddsratio, pvalue = stats.fisher_exact([[mis_col, mis_other], [occ_col, occ_other]])
        missense_scores.append((oddsratio, pvalue))

    # Parse to dataframe
    return pd.DataFrame(missense_scores, columns=['oddsratio', 'pvalue'], index=column_variant_counts.index)


def _missense_per_residue_plot(aligned_variants_table, mapping_table):
    """

    :param aligned_variants_table:
    :return:
    """
    residue_counts = aligned_variants_table.pipe(_aggregate_annotation,
                                                 ('VEP', 'Consequence'),
                                                 aggregate_by=['SOURCE_ID', 'Protein_position'])
    residue_counts = residue_counts.reindex(mapping_table.index).fillna(0)  # Fill in residues with no variants
    ax = residue_counts['missense_variant'].astype(int).value_counts().plot.bar(width=0.9, facecolor='black',
                                                                                edgecolor='black')

    return residue_counts.rename_axis('', 1)


def _variant_per_protein_plot(aligned_variants_table):
    """

    :param aligned_variants_table:
    :return:
    """
    protein_consequences = _aggregate_annotation(aligned_variants_table, ('VEP', 'Consequence'),
                                                 aggregate_by=['SOURCE_ID'])
    ax = protein_consequences.loc[:, protein_consequences.sum() > 100].hist(facecolor='black', edgecolor='black',
                                                                            figsize=(10, 10))

    return protein_consequences.rename_axis('', 1)


def _variants_vs_length_plot(protein_variant_counts, alignment_info):
    """

    :param protein_variant_counts:
    :param alignment_info:
    :return:
    """
    # Calculate sequence lengths
    protein_variant_counts = protein_variant_counts.join(alignment_info['lengths'])

    # Plot
    plot_data = pd.melt(protein_variant_counts.loc[:, protein_variant_counts.sum() > 100],
                        id_vars=['length'],
                        var_name='Variant_Effect', value_name='Count')
    sns.lmplot(x='length', y='Count', col='Variant_Effect', hue='Variant_Effect',
               data=plot_data, fit_reg=True, sharey=False, col_wrap=3)

    return None



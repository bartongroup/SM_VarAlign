

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
    return grouped.size().unstack(fill_value=fill_value)

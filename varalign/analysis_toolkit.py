

def _aggregate_annotation(aligned_variant_table, annotation_column, fill_value=0):
    """

    :param aligned_variant_table:
    :param annotation_column:
    :return:
    """
    grouped = aligned_variant_table.groupby([('Alignment', 'Column'), annotation_column])
    return grouped.size().unstack(fill_value=fill_value)

"""
This module contains functions that write Jalview annotation files.
"""
import logging

log = logging.getLogger(__name__)


def write_jalview_annotation(ordered_values, file_name, title, description, append=False, tooltips=None):
    """
    Create and/or append tracks to a Jalview alignment annotation file.

    :param ordered_values: Tuple or list of tuples containing annotation scores for each column in alignment.
    :param file_name: Filename for Jalview output.
    :param title: Tuple or list of tuples of titles for each annotation track.
    :param description: As above containing descriptions.
    :param append: Whether to append to an existing Jalview annotation file.
    :return:
    """
    file_mode = 'w'
    if append:
        file_mode = 'a'

    if not tooltips:
        tooltips = ordered_values

    with open(file_name, file_mode) as results_file:
        if not append:
            results_file.write('JALVIEW_ANNOTATION\n')
        if isinstance(ordered_values, tuple) and all(map(lambda x: isinstance(x, str), [title, description])):
            results_file.write('BAR_GRAPH\t{}\t{}\t'.format(title, description) +
                               '|'.join('{},,{}'.format(str(x), str(y)) for x, y in zip(ordered_values, tooltips)))
            results_file.write('\n')
        elif all(map(lambda x: isinstance(x, list), [ordered_values, title, description])):
            arg_lengths = map(len, [ordered_values, title, description])
            if len(set(arg_lengths)) == 1:
                for v, vv, t, d in zip(ordered_values, tooltips, title, description):
                    results_file.write('BAR_GRAPH\t{}\t{}\t'.format(t, d) +
                                       '|'.join('{},,{}'.format(str(x), str(y)) for x, y in zip(v, vv)))
                    results_file.write('\n')
            else:
                log.error('List arguments must be of same length')
                raise TypeError
        else:
            log.error('Must provide same number of titles/descriptions as lists of values.')
            raise TypeError

    return 0


def append_jalview_variant_features(seq_id, positions, descriptions, feature_type, file_name):
    """
    Append features to a Jalview features file.

    :param seq_id: The sequence ID that the features are associated with.
    :param positions: The sequence positions of the features.
    :param descriptions: The description for each feature.
    :param feature_type: The Jalview featureType.
    :param file_name: Jalview feature file to append.
    :return:

    feature_list format: [description, sequenceId, sequenceIndex, start, end, featureType]
    """
    with open(file_name, 'a') as output:
        for pos, desc in zip(positions, descriptions):
            feature_list = [desc, seq_id, '-1', str(pos), str(pos), feature_type, '0.0']
            output.write('\t'.join(feature_list) + '\n')


def create_jalview_feature_file(feature_dict, file_name):
    """
    Write feature definitions to a Jalview features file.

    :param feature_dict: {'featureType': colour}
    :param file_name: Filename for Jalview output.
    :return:
    """
    with open(file_name, 'w') as output:
        for feature, colour in feature_dict.items():
            output.write('{}\t{}\n'.format(feature, colour))


def marked_columns_track(mask, title, description, filename):
    """Write NO_GRAPH jalview track for a mask with '*' for True and '' for False.

    :param mask:
    :param title:
    :param description:
    :param filename:
    :return:
    """
    values = mask.copy()
    values.loc[mask] = '*'
    values.loc[~mask] = ''

    with open(filename, 'wb') as f:
        f.write('JALVIEW_ANNOTATION\n')
        f.write('NO_GRAPH\t{}\t{}\t'.format(title, description) + '|'.join(values.tolist()))
        f.write('\n')


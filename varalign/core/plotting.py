# This module contains all plotting-related functions.

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pandas as pd


def _missense_per_residue_plot(aligned_variants_table, mapping_table):
    """

    :param aligned_variants_table:
    :return:
    """
    residue_counts = aligned_variants_table.pipe(
        _aggregate_annotation,
        ("VEP", "Consequence"),
        aggregate_by=["SOURCE_ID", "Protein_position"],
    )
    residue_counts = residue_counts.reindex(mapping_table.index).fillna(
        0
    )  # Fill in residues with no variants
    ax = (
        residue_counts["missense_variant"]
        .astype(int)
        .value_counts()
        .plot.bar(width=0.9, facecolor="black", edgecolor="black")
    )

    return residue_counts.rename_axis("", 1)


def _variant_per_protein_plot(aligned_variants_table):
    """

    :param aligned_variants_table:
    :return:
    """
    protein_consequences = _aggregate_annotation(
        aligned_variants_table, ("VEP", "Consequence"), aggregate_by=["SOURCE_ID"]
    )
    ax = protein_consequences.loc[:, protein_consequences.sum() > 100].hist(
        facecolor="black", edgecolor="black", figsize=(10, 10)
    )

    return protein_consequences.rename_axis("", 1)


def _variants_vs_length_plot(protein_variant_counts, alignment_info):
    """

    :param protein_variant_counts:
    :param alignment_info:
    :return:
    """
    # Calculate sequence lengths
    protein_variant_counts = protein_variant_counts.join(alignment_info["length"])

    # Plot
    plot_data = pd.melt(
        protein_variant_counts.loc[:, protein_variant_counts.sum() > 100],
        id_vars=["length"],
        var_name="Variant_Effect",
        value_name="Count",
    )
    sns.lmplot(
        x="length",
        y="Count",
        col="Variant_Effect",
        hue="Variant_Effect",
        data=plot_data,
        fit_reg=True,
        sharey=False,
        col_wrap=3,
    )

    return None


def align_variants_plot_function_1():
    pdf = PdfPages(results_prefix + ".figures.pdf", metadata={"creationDate": None})
    # PDF metadata
    d = pdf.infodict()
    d["Title"] = "Aligned Variant Diagnostics Plots for {}".format(path_to_alignment)
    d["Author"] = "align_variants.py"

    # Plot GMM diagnostics
    occ_gmm._gmm_plot(M_best, gmms["models"], gmms["data"])
    pdf.attach_note("Residue Occupancy GMM Diagnostics")
    pdf.savefig(metadata={"creationDate": None})
    plt.close()


def align_variants_plot_function_2():
    # Plot 1
    fig, axs = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)
    column_summary.plot.scatter("occupancy", "missense_variant", ax=axs[0])
    column_summary.plot.scatter("occupancy", "synonymous_variant", ax=axs[1])
    axs[0].axvline(column_summary[subset_mask_gmm]["occupancy"].min())
    axs[1].axvline(column_summary[subset_mask_gmm]["occupancy"].min())
    fig.suptitle("N Variants vs. Occupancy")
    pdf.attach_note("N Variants vs. Occupancy")
    pdf.savefig(metadata={"creationDate": None})
    plt.close()


def align_variants_plot_function_3():
    # Conservation plane plot: Variant counts vs. Shenkin
    fig, axs = plt.subplots(1, 2, figsize=(15, 5), sharex=True, sharey=True)
    sns.regplot(
        x="shenkin",
        y="missense_variant",
        data=column_summary[subset_mask_gmm],
        ax=axs[0],
    )
    pd.plotting.table(
        axs[0],
        shenkin_regressions.loc[["missense", "filtered_missense"]].round(2),
        loc="upper right",
        colWidths=[0.12] * 5,
        zorder=100,
    )
    sns.regplot(
        x="shenkin",
        y="synonymous_variant",
        data=column_summary[subset_mask_gmm],
        ax=axs[1],
    )
    pd.plotting.table(
        axs[1],
        shenkin_regressions.loc[["synonymous", "filtered_synonymous"]].round(2),
        loc="upper right",
        colWidths=[0.12] * 5,
        zorder=100,
    )
    plt.title("N Variants vs. Shenkin")
    pdf.attach_note("N Variants vs. Shenkin")
    pdf.savefig(metadata={"creationDate": None})
    plt.close()


def align_variants_plot_function_4():
    # Conservation plane plot: Missense Scores vs. Shenkin
    plot_data = column_summary[subset_mask_gmm]
    plot_data = plot_data.assign(pass_alpha=plot_data["pvalue"] < 0.1)
    log.info("plot_data:\n%s", plot_data.head().to_string())
    # TODO: plot.scatter throws AttributeError with pandas 0.22.0 or matplotlib 2.1.2
    # ax = plot_data.plot.scatter('shenkin', 'oddsratio', c='pass_alpha',  # Valdar is well correlated...
    #                            colorbar=False,
    #                            logy=True, figsize=(10, 10))
    # _ = plt.setp(ax.get_xticklabels(), visible=True)
    # plt.title('Missense Score vs. Shenkin')
    # pdf.attach_note('Missense Score vs. Shenkin')
    # pdf.savefig(metadata={'creationDate': None})


def align_variants_plot_function_5():
    # Variants per sequence histogram
    protein_consequences = analysis_toolkit._aggregate_annotation(
        alignment_variant_table, ("VEP", "Consequence"), aggregate_by=["SOURCE_ID"]
    )
    protein_consequences.hist(facecolor="black", edgecolor="black")
    plt.title("Variants per Sequence")
    pdf.attach_note("Distribution of variants over alignment sequences")
    pdf.savefig(metadata={"creationDate": None})
    # Plotting logic has been moved to 'plotting.py'
    # Variants per residue and column histograms
    fig, axes = plt.subplots(1, 2)
    residue_counts = alignment_variant_table.pipe(
        analysis_toolkit._aggregate_annotation,
        ("VEP", "Consequence"),
        aggregate_by=["SOURCE_ID", "Protein_position"],
    )
    residue_counts = residue_counts.reindex(indexed_mapping_table.index).fillna(
        0
    )  # Fill in residues with no variants
    residue_counts["missense_variant"].astype(int).value_counts().plot.bar(
        ax=axes[0], width=1, facecolor="black", edgecolor="black"
    )
    axes[0].set_title("Missense Variants per Residue")
    # column_variant_counts['missense_variant'].hist(ax=axes[0])
    column_summary.loc[subset_mask_gmm, "missense_variant"].hist(
        ax=axes[1], facecolor="black", edgecolor="black"
    )
    axes[1].set_title("Missense Variants per Column")
    pdf.attach_note("Distribution of variants over residues and alignment columns")
    pdf.savefig(metadata={"creationDate": None})
    plt.close()
    pdf.close()


def prointvar_analysis_plot_function_1():
    plt.setp(markers, marker="o", markersize=3)
    _ = plt.setp(stemlines, linewidth=0.5)
    _ = plt.setp(baseline, linewidth=1, visible=True)


def prointvar_analysis_plot_function_2():
    plt.subplots()
    # _ = ax.plot(plot_data['sifts_index'], '.')
    # ax.xaxis.set_major_locator(IndexLocator(10, 0))
    # # Label top three
    # for i, row in plot_data.iloc[-3:].iterrows():
    #     _ = ax.annotate(row.query, (i, row.sifts_index),
    #                     xytext=(0.85 * i, row.sifts_index))
    # pdf.attach_note('Available PDBs from SIFTS')
    # pdf.savefig()
    # plt.close()
    # Structure analyses
    structure_stats = prointvar_stats.collect_column_structure_stats(structure_table)
    # Add variant column stats
    column_stats.rename(
        columns={"('Alignment', 'Column')": "Alignment_column"}, inplace=True
    )
    # Add Shenkin percentile score
    if "shenkin_percentile" not in column_stats.columns:
        # subset_mask_gmm??
        column_stats = column_stats.join(
            column_stats["shenkin"].rank(pct=True), rsuffix="_percentile"
        )
    column_annotations = column_stats.set_index("Alignment_column").join(
        structure_stats
    )
    column_annotations.to_csv(results_prefix + "_column_data.csv")
    # Plot comparing structural features on the alignment
    fig, axs = plt.subplots(2, 1, sharex=True, figsize=(9, 8))
    alignment_ligand_plot(structure_stats, axs[0])
    alignment_ppi_plot(structure_stats, axs[1])
    _ = axs[0].set_ylabel("Protein-Ligand Interactions")
    _ = axs[1].set_ylabel("Protein-Protein Interactions")
    _ = axs[1].set_xlabel("Alignment column")
    _ = plt.suptitle("Distirbution of Structural Features on Alignment")
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    pdf.attach_note("Ligand and PPI residue distribution.")
    pdf.savefig()
    plt.close()
    pdf.close()
    log.info("DONE.")


def occ_gmm_plot_function_1():
    plt.figure(figsize=(10, 5))
    fig.subplots_adjust(left=0.12, right=0.97, bottom=0.21, top=0.9, wspace=0.5)

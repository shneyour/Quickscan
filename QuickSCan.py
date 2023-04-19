import gseapy as gp
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np


def enrichr_plot(gene_list, color, output_path, cell1=None, cell2=None, i=None, cluster=None) -> None:
    pointer = 'Upregulated' if i == 0 else 'Downregulated'
    try:
        enr = gp.enrichr(gene_list=gene_list,  # or "./tests/data/gene_list.txt",
                         gene_sets='GO_Biological_Process_2021',
                         organism='human',
                         outdir=None)
        enr.results.Term = enr.results.Term.str.split(" \(GO").str[0]
        terms = enr.results[enr.results['Adjusted P-value'] < 0.05]
        terms['-log10(P value)'] = -np.log10(terms['Adjusted P-value'])
        terms = terms.sort_values(by='-log10(P value)', ascending=True)
        terms = terms.tail(5)
        y = terms['Term']
        x = terms['-log10(P value)']
        fig, ax = plt.subplots(figsize=(5.3, 3))
        ax.spines['top'].set_visible(True)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(True)
        ax.xaxis.tick_top()
        width = 0.6
        ind = np.arange(len(x))
        ax.barh(ind, x, width, color=color,align='edge')
        plt.tick_params(left=False, labelleft=False, bottom=False)
        if (x.max().round() / 10).round() != 0:
            plt.xticks(np.arange(0, x.max().round() + 1, (x.max().round() / 10).round()))
        else:
            plt.xticks(np.arange(0, x.max().round() + 1, 1))
        plt.margins(0, 0.05)
        plt.xlabel(r'$- \log_{10}$(Adj. $P$-Value)')
        ax.xaxis.set_label_coords(0.5, 1.2)
        for bar, term in zip(ax.patches, list(y)):
            ax.text(0.1, bar.get_y() + bar.get_height() / 2, term, color='black', ha='left', va='center')
        if cell1:
            plt.title(f'{cell1} vs {cell2} - {pointer}', loc='center')
            plt.tight_layout()
            plt.savefig(output_path + cell1 + 'vs' + cell2 + pointer + '.png', dpi=300)
        else:
            plt.title(f'{cluster}', loc='center')
            plt.tight_layout()
            plt.savefig(output_path + cluster + '.png', dpi=300)

    except ValueError:
        if cell1:
            print(f'No processes for {cell1} vs {cell2} - {pointer}')
        else:
            print(f'No processes for {cluster}')


def preprocess_enrichr(data, output_path, groups, cell1=None, cell2=None) -> None:
    if groups:
        color_list = ['bisque', 'lightsalmon', 'lavender', 'thistle', 'khaki', 'lightsteelblue', 'darkseagreen',
                      'lightcyan', 'mediumaquamarine', 'honeydew',
                      'cornsilk', 'mistyrose', 'whitesmoke', 'lightgrey'] # nice colors I chose
        clusters = data['group'].unique().tolist()
        for i, cluster in enumerate(clusters):
            gene_list_enrichr = data.index[(data['group'] == cluster) & (data['logfoldchanges'] > 0.4)].tolist()
            color = color_list[i]
            enrichr_plot(gene_list_enrichr, color, output_path, cluster=cluster)
    else:
        colors = ['lightcoral', 'powderblue']
        for i in range(0, 2):
            gene_list_enrichr = data.index[data['logfoldchanges'] > 0.4].tolist() if i == 0 else data.index[data['logfoldchanges'] < -0.4].tolist()
            color = colors[i]
            enrichr_plot(gene_list_enrichr, color, output_path, cell1=cell1, cell2=cell2, i=i)


def analyze_sc(data_path, output_path, ann_column, cell1=None, cell2=None, general=True, pathways=True) -> None:
    adata = sc.read_h5ad(data_path)
    data = adata.obs
    data_dict = {}
    for cluster in data[ann_column].unique():
        data_dict[cluster] = len(data[data[ann_column] == cluster])
    fig, ax = plt.subplots(figsize=(6, 3), subplot_kw=dict(aspect="equal"))
    sizes = np.array(list(data_dict.values()))
    percent = 100. * sizes / sizes.sum()
    labels = list(data_dict.keys())
    colors = ['yellowgreen', 'red', 'gold', 'lightskyblue', 'purple', 'lightcoral', 'blue', 'pink', 'darkgreen',
              'yellow', 'grey', 'violet', 'magenta', 'cyan']
    labels_fixed = ['{0} - {1:1.2f}%'.format(i, j) for i, j in zip(labels, percent)]
    patches, texts = plt.pie(sizes, colors=colors, startangle=90, radius=1)
    sort_legend = True
    if sort_legend:
        patches, labels, dummy = zip(*sorted(zip(patches, labels_fixed, sizes),
                                             key=lambda x: x[2],
                                             reverse=True))

    ax.legend(patches, labels, loc='upper left', bbox_to_anchor=(1, 0, 0.5, 1),
              fontsize=8)
    ax.set_title(ann_column)
    fig.savefig(output_path + ann_column + '_pie.png',
                bbox_inches='tight', dpi=300)

    if general:
        sc.tl.rank_genes_groups(adata, ann_column, method='wilcoxon')
        sc.pl.rank_genes_groups_matrixplot(
            adata,
            n_genes=4,
            values_to_plot="logfoldchanges",
            cmap='bwr',
            vmin=-4,
            vmax=4,
            min_logfoldchange=3,
            colorbar_title='log fold change',
            dendrogram=False,
            show=False
        )
        plt.savefig(output_path + ann_column + '_DEG_heatmap.png', dpi=300)
        plt.clf()
        sc.pl.rank_genes_groups(adata, show=False)
        plt.savefig(output_path + f'{ann_column}.png', dpi=300)

    else:
        sc.tl.rank_genes_groups(adata, ann_column, groups=[cell1], reference=cell2, method='wilcoxon')
        sc.pl.rank_genes_groups(adata, show=False)
        plt.savefig(output_path + f'{cell1}_vs_{cell2}_rank.png', dpi=300)

    if pathways:
        test_data = sc.get.rank_genes_groups_df(adata, group=None)
        test_data = test_data[test_data['pvals_adj'] < 0.05].sort_values(by='pvals', ascending=True).set_index('names')
        if 'group' in test_data.columns:
            preprocess_enrichr(data=test_data, groups=True, output_path=output_path)

        else:
            preprocess_enrichr(data=test_data, groups=False, output_path=output_path, cell1=cell1, cell2=cell2)


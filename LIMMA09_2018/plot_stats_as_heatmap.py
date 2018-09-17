import pandas, os, glob, seaborn
import matplotlib.pyplot as plt

"""
We are primarily interest in three questions:
   1) Does gene g respond to TGFb? - Compare TGFb vs Control groups
   2) Is the gene g response different in adult/senescent compared to neonatal cell lines? - Compare Control groups for cell lines 
   3) Is the gene g response to TGFb different in adult/senescent compared to neonatal cell lines? - Compare treated groups for cell lines

Analysis has been conducted in R, with LIMMA. This script is only for plotting. 

"""
PLOT_WITHIN = False
PLOT_BETWEEN = False
PLOT_PVAL_GRAPH = True

seaborn.set_context('talk', font_scale=2)
seaborn.set_style('white')

directory = r'/home/b3053674/Documents/LargeStudy/LIMMA09-2018'
saved_objects_path = os.path.join(directory, 'SavedObjects')

pval_settings_file = os.path.join(directory, 'pval')
with open(pval_settings_file) as f:
    PVAL = f.read().strip()

print('pval is "{}"'.format(PVAL))
pval_path = os.path.join(saved_objects_path, 'pval_less_than_{}'.format(str(PVAL).replace('.', '_')))
if not os.path.isdir(pval_path):
    print('makeing file', pval_path)
    os.mkdir(pval_path)  # if os.path.isdir(pval_path) is False else None

between_control_path = os.path.join(pval_path, 'between_control_statistics.csv')
between_tgfb_path = os.path.join(pval_path, 'between_tgfb_statistics.csv')
within_neonatal_path = os.path.join(pval_path, 'within_neonatal.csv')
within_adult_path = os.path.join(pval_path, 'within_adult.csv')
within_senescent_path = os.path.join(pval_path, 'within_senescent.csv')

between_control_data = pandas.read_csv(between_control_path, index_col=0)
between_tgfb_data = pandas.read_csv(between_tgfb_path, index_col=0)
within_neonatal_data = pandas.read_csv(within_neonatal_path, index_col=0)
within_adult_data = pandas.read_csv(within_adult_path, index_col=0)
within_senescent_data = pandas.read_csv(within_senescent_path, index_col=0)

between_control_data.index.name = 'Gene'
between_tgfb_data.index.name = 'Gene'
within_neonatal_data.index.name = 'Gene'
within_adult_data.index.name = 'Gene'
within_senescent_data.index.name = 'Gene'

# print(between_tgfb_data)
between_control_data = between_control_data[['sen_perc', 'ad_perc']]
between_tgfb_data = between_tgfb_data[['sen_perc', 'ad_perc']]

within_neonatal_data.columns = pandas.MultiIndex.from_product([['within_neonatal'], list(within_neonatal_data.columns)])
within_adult_data.columns = pandas.MultiIndex.from_product([['within_adult'], list(within_adult_data.columns)])
within_senescent_data.columns = pandas.MultiIndex.from_product(
    [['within_senescent'], list(within_senescent_data.columns)])
between_control_data.columns = pandas.MultiIndex.from_product([['between_control'], list(between_control_data.columns)])
between_tgfb_data.columns = pandas.MultiIndex.from_product([['between_tgfb'], list(between_tgfb_data.columns)])

df = pandas.concat(
    [within_neonatal_data, within_senescent_data, within_adult_data, between_control_data, between_tgfb_data], axis=1
)
df = df.stack().stack()
df = pandas.DataFrame(df)
df.index.names = ['gene', 'group1', 'group2']

df.columns = ['percentage']

df = df.reset_index()
tgf = r'TGF$\beta$'

label_down = -0.095

if PLOT_WITHIN:
    within = df.query("group2 in ['within_adult', 'within_senescent', 'within_neonatal']")
    # within = df[df['group2'].any() in ]

    new_col = ["{}_{}".format(
        within['group2'].iloc[i], within['group1'].iloc[i]) for i in range(within.shape[0])
    ]
    within['group'] = new_col
    within = within.pivot(index='gene', columns='group', values='percentage')
    group_labels = ['adult_control', 'adult_tgf',
                    'neonatal_control', 'neonatal_tgf',
                    'senescent_control', 'senescent_tgf']
    within.columns = group_labels
    new_order = [group_labels[2], group_labels[0], group_labels[4],
                 group_labels[3], group_labels[1], group_labels[5]]
    within = within[new_order]
    fig, ax = plt.subplots(figsize=(10, 20), dpi=350)

    # print(within['neonatal_tgf'])
    print('within, pval ({}) count below 60%'.format(PVAL))
    count = pandas.DataFrame(within[within > 60].count())
    count.columns = ['Count']
    cell, treat = zip(*[i.split('_') for i in list(count.index)])
    count['comparison'] = cell
    count['treatment'] = treat
    count = count.set_index(['comparison', 'treatment'])
    count = count.unstack(level=1)
    count_fname = os.path.join(pval_path, 'within_count_greater_than_60_percent.csv')
    count.to_csv(count_fname)
    print(count)

    ax = seaborn.heatmap(within, linewidths=1, cmap='YlGnBu',
                         ax=ax, linecolor='black',
                         cbar_kws={'label': '% p-value < {}'.format(PVAL)})

    ax.set_xticklabels(['Neo', 'Adult', 'Sen'] * 2, rotation=90)
    plt.yticks(rotation=0, fontsize=16)
    plt.ylabel('')
    plt.xlabel('')

    trans = ax.get_xaxis_transform()
    down_loc = -380

    ax.annotate('Control', xy=(1.0, label_down), xycoords=trans)
    ax.annotate('', xycoords='axes pixels', xytext=(100, down_loc), xy=(1000, down_loc),
                arrowprops=dict(arrowstyle='-',
                                color='black',
                                linewidth=5)
                )

    ax.annotate(tgf, xy=(4.0, label_down), xycoords=trans)
    ax.annotate('', xycoords='axes pixels', xytext=(1200, down_loc), xy=(2100, down_loc),
                arrowprops=dict(arrowstyle='-',
                                color='black',
                                linewidth=5)
                )

    fname = os.path.join(pval_path, 'within_heatmap_{}'.format(str(PVAL).replace('.', '_')))
    fig.savefig(fname, dpi=350, bbox_inches='tight')

if PLOT_BETWEEN:
    between = df.query("group2 in ['between_control', 'between_tgfb']")

    new_col = [
        "{}_{}".format(
            between['group1'].iloc[i], between['group2'].iloc[i]
        ) for i in range(between.shape[0])
    ]

    between['group'] = new_col

    between = between.pivot_table(index='gene', columns='group', values='percentage')
    new_labels = ['adult_control', 'adult_tgfb', 'senescent_control', 'senescent_tgfb']
    between.columns = new_labels
    new_order = [new_labels[0], new_labels[2], new_labels[1], new_labels[3]]
    between = between[new_order]

    print('between, pval ({}) count below 60%'.format(PVAL))
    count = pandas.DataFrame(between[between > 60].count())
    count.columns = ['Count']
    cell, treat = zip(*[i.split('_') for i in list(count.index)])
    count['comparison'] = cell
    count['treatment'] = treat
    count = count.set_index(['comparison', 'treatment'])
    count = count.unstack(level=1)
    print(count)
    count_fname = os.path.join(pval_path, 'between_count_greater_than_60_percent.csv')
    count.to_csv(count_fname)

    fig, ax = plt.subplots(figsize=(10, 20), dpi=350)
    # [YlOrRd,YlOrBr, YlGnBu]
    ax = seaborn.heatmap(between, linewidths=1, cmap='YlGnBu',
                         ax=ax, linecolor='black',
                         cbar_kws={'label': '% p-value < {}'.format(PVAL)})

    ax.set_xticklabels(['Adult', 'Sen'] * 2, rotation=90)
    plt.yticks(rotation=0, fontsize=16)
    plt.ylabel('')
    plt.xlabel('')
    trans = ax.get_xaxis_transform()

    down_loc = -370
    ax.annotate('Control', xy=(0.5, label_down), xycoords=trans)
    ax.annotate('', xycoords='axes pixels', xytext=(100, down_loc), xy=(1000, down_loc),
                arrowprops=dict(arrowstyle='-',
                                color='black',
                                linewidth=5)
                )

    ax.annotate(tgf, xy=(2.75, label_down), xycoords=trans)
    ax.annotate('', xycoords='axes pixels', xytext=(1200, down_loc), xy=(2100, down_loc),
                arrowprops=dict(arrowstyle='-',
                                color='black',
                                linewidth=5)
                )

    fname = os.path.join(pval_path, 'between_heatmap{}'.format(str(PVAL).replace('.', '_')))
    fig.savefig(fname, dpi=350, bbox_inches='tight')

if PLOT_PVAL_GRAPH:
    import glob

    dirs = glob.glob(os.path.join(saved_objects_path, 'pval_less_than_*'))
    all_pvals = ['{}.{}'.format(i.split('_')[-2:][0], i.split('_')[-2:][1]) for i in dirs]
    all_pvals = [float(i) for i in all_pvals]
    dirs = zip(dirs, all_pvals)
    between_df_dct = {}
    within_df_dct = {}
    for (d, p) in dirs:
        assert os.path.isdir(d)
        print(d)
        between_file = glob.glob(os.path.join(d, 'between_count*'))
        within_file = glob.glob(os.path.join(d, 'within_count*'))
        assert len(between_file) == 1, between_file
        assert len(within_file) == 1
        between_count = pandas.read_csv(between_file[0], index_col=[0], skiprows=[0]).dropna()
        within_count = pandas.read_csv(within_file[0], index_col=[0], skiprows=[0]).dropna()
        between_df_dct[p] = between_count
        within_df_dct[p] = within_count
        # print(within_count)
    between = pandas.concat(between_df_dct)
    within = pandas.concat(within_df_dct)
    within.index.names = ['p_val', 'comparison']
    between.index.names = ['p_val', 'comparison']
    between = between.unstack(level=1).sort_index(ascending=False)
    within = within.unstack(level=1).sort_index(ascending=False)

    within = within.stack().stack().reset_index()
    new_col = []
    for i in range(within.shape[0]):
        new_col.append("{}_{}".format(
            within['comparison'].iloc[i],
            within['level_2'].iloc[i])
        )
    within['label'] = new_col
    # print(within)
    fig, ax = plt.subplots()
    seaborn.barplot(data=within, x='p_val', y=0, hue='label', ax=ax,
                    hue_order=[
                        'neonatal_control', 'adult_control', 'senescent_control',
                        'neonatal_tgf', 'adult_tgf', 'senescent_tgf'
                    ],
                    order=[
                        0.01, 0.001, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10
                    ])
    seaborn.despine(fig=fig, top=True, right=True)
    L = plt.legend(loc=(1, 0.1))
    L.get_texts()[0].set_text('neonatal, control')
    L.get_texts()[1].set_text('adult, control')
    L.get_texts()[1].set_text('senescent, control')
    L.get_texts()[2].set_text(r'neonatal, TGF$\beta$')
    L.get_texts()[3].set_text(r'adult, TGF$\beta$')
    L.get_texts()[3].set_text(r'senescent, TGF$\beta$')
    plt.xticks(rotation=90)
    plt.ylabel('Count >60%')
    plt.xlabel('FDR corrected p-value cut-off')
    fname = os.path.join(saved_objects_path, 'within_pvalue_counts.png')
    fig.savefig(fname, bbox_inches='tight', dpi=100)

    between = between.stack().stack().reset_index()
    new_col = []
    for i in range(between.shape[0]):
        new_col.append("{}_{}".format(
            between['comparison'].iloc[i],
            between['level_2'].iloc[i])
        )
    between['label'] = new_col
    print('between')
    print(between)
    fig, ax = plt.subplots()
    seaborn.barplot(data=between, x='p_val', y=0, hue='label', ax=ax,
                    order=[
                        0.01, 0.001, 0.0001, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10
                    ],
                    hue_order=[
                        'adult_control', 'senescent_control',
                        'adult_tgfb', 'senescent_tgfb'
                    ]
                    )
    seaborn.despine(fig=fig, top=True, right=True)
    L = plt.legend(loc=(1, 0.1))
    L.get_texts()[0].set_text('adult, control')
    L.get_texts()[1].set_text('senescent, control')
    L.get_texts()[2].set_text(r'adult, TGF$\beta$')
    L.get_texts()[3].set_text(r'adult, TGF$\beta$')
    # plt.legend(loc=(1, 0.1))
    plt.xticks(rotation=90)
    plt.ylabel('Count >60%')
    plt.xlabel('FDR corrected p-value cut-off')
    fname = os.path.join(saved_objects_path, 'between_pvalue_counts.png')
    fig.savefig(fname, bbox_inches='tight', dpi=100)
    # print(between)

# if MODE == 1:
#
#     if PLOT_WITHIN:
#         within = df.query("group2 in ['within_adult', 'within_senescent', 'within_neonatal']")
#         # within = df[df['group2'].any() in ]
#
#         new_col = ["{}_{}".format(
#             within['group2'].iloc[i], within['group1'].iloc[i]) for i in range(within.shape[0])
#         ]
#         within['group'] = new_col
#         within = within.pivot(index='gene', columns='group', values='percentage')
#         group_labels = ['adult, control', 'adult, {}'.format(tgf),
#                         'neonatal, control', 'neonatal, {}'.format(tgf),
#                         'senescent, control', 'senescent, {}'.format(tgf)
#                         ]
#         within.columns = group_labels
#         new_order = [group_labels[2], group_labels[3],
#                      group_labels[0], group_labels[1],
#                      group_labels[4], group_labels[5]]
#         within = within[new_order]
#
#         fig, ax = plt.subplots(figsize=(10, 20), dpi=350)
#         ax = seaborn.heatmap(within, linewidths=1, cmap='YlGnBu',
#                              ax=ax, linecolor='black',
#                              cbar_kws={'label': '% p-value < 0.001'})
#
#         ax.set_xticklabels(['Control', tgf] * 3, rotation=90)
#         plt.yticks(rotation=0, fontsize=16)
#         plt.ylabel('')
#         plt.xlabel('')
#
#         trans = ax.get_xaxis_transform()
#         ax.annotate('Neonatal', xy=(0.25, label_down), xycoords=trans)
#         down_loc = -480
#         gap = 40
#         length = 160
#         start_loc = 20
#         ax.annotate('', xycoords='axes pixels', xytext=(50, down_loc), xy=(700, down_loc),
#                     arrowprops=dict(arrowstyle='-',
#                                     color='black',
#                                     linewidth=5)
#                     )
#
#         ax.annotate('Adult', xy=(2.5, label_down), xycoords=trans)
#         ax.annotate('', xycoords='axes pixels', xytext=(800, down_loc), xy=(1400, down_loc),
#                     arrowprops=dict(arrowstyle='-',
#                                     color='black',
#                                     linewidth=5)
#                     )
#
#         ax.annotate('Senescent', xy=(4.1, label_down), xycoords=trans)
#         ax.annotate('', xycoords='axes pixels', xytext=(1500, down_loc), xy=(2100, down_loc),
#                     arrowprops=dict(arrowstyle='-',
#                                     color='black',
#                                     linewidth=5)
#                     )
#         fname = os.path.join(saved_objects_path, 'within_heatmap0_001')
#         fig.savefig(fname, dpi=350, bbox_inches='tight')
#
#     if PLOT_BETWEEN:
#         between = df.query("group2 in ['between_control', 'between_tgfb']")
#
#         new_col = [
#             "{}_{}".format(
#                 between['group1'].iloc[i], between['group2'].iloc[i]
#             ) for i in range(between.shape[0])
#         ]
#
#         between['group'] = new_col
#
#         between = between.pivot_table(index='gene', columns='group', values='percentage')
#         between.columns = ['adult_control', 'adult_tgfb', 'senescent_control', 'senescent_tgfb']
#
#         fig, ax = plt.subplots(figsize=(10, 20), dpi=350)
#         # [YlOrRd,YlOrBr, YlGnBu]
#         ax = seaborn.heatmap(between, linewidths=1, cmap='YlGnBu',
#                              ax=ax, linecolor='black',
#                              cbar_kws={'label': '% p-value < 0.001'})
#
#         ax.set_xticklabels(['Control', tgf] * 3, rotation=90)
#         plt.yticks(rotation=0, fontsize=16)
#         plt.ylabel('')
#         plt.xlabel('')
#         trans = ax.get_xaxis_transform()
#
#         down_loc = -480
#         gap = 40
#         length = 160
#         start_loc = 20
#
#         ax.annotate('Adult', xy=(0.5, label_down), xycoords=trans)
#         ax.annotate('', xycoords='axes pixels', xytext=(100, down_loc), xy=(1000, down_loc),
#                     arrowprops=dict(arrowstyle='-',
#                                     color='black',
#                                     linewidth=5)
#                     )
#
#         ax.annotate('Senescent', xy=(2.3, label_down), xycoords=trans)
#         ax.annotate('', xycoords='axes pixels', xytext=(1200, down_loc), xy=(2100, down_loc),
#                     arrowprops=dict(arrowstyle='-',
#                                     color='black',
#                                     linewidth=5)
#                     )
#
#         fname = os.path.join(saved_objects_path, 'between_heatmap')
#         fig.savefig(fname, dpi=350, bbox_inches='tight')

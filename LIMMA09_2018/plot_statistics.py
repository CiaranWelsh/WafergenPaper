import pandas, os, glob, seaborn
import matplotlib.pyplot as plt

"""
We are primarily interest in three questions:
   1) Does gene g respond to TGFb? - Compare TGFb vs Control groups
   2) Is the gene g response different in adult/senescent compared to neonatal cell lines? - Compare Control groups for cell lines 
   3) Is the gene g response to TGFb different in adult/senescent compared to neonatal cell lines? - Compare treated groups for cell lines

Analysis has been conducted in R, with LIMMA. This script is only for plotting. 

"""
PVAL = 0.001


seaborn.set_context('talk', font_scale=2)
seaborn.set_style('white')

directory = r'/home/b3053674/Documents/LargeStudy/LIMMA09-2018'
saved_objects_path = os.path.join(directory, 'SavedObjects')

## first do the 'between' data plots

between_control_path = os.path.join(saved_objects_path, 'between_control_statistics.csv')
between_tgfb_path = os.path.join(saved_objects_path, 'between_tgfb_statistics.csv')

between_control_data = pandas.read_csv(between_control_path, index_col=0)
between_tgfb_data = pandas.read_csv(between_tgfb_path, index_col=0)

between_control_data.index.name = 'Gene'
between_tgfb_data.index.name = 'Gene'

between_control_data.reset_index(inplace=True)
between_tgfb_data.reset_index(inplace=True)


def plot(data, y=None, fname=None):
    if y is None:
        print(data.head())
        raise ValueError('y cannot be None')

    import matplotlib
    cmap = matplotlib.cm.get_cmap('gist_rainbow')
    colours = []
    count = 0
    num_colours_needed = data.shape[0]
    for i in range(num_colours_needed):
        colours.append(cmap(count))
        count += 1.0 / num_colours_needed

    fig, ax = plt.subplots(figsize=(30, 10))
    seaborn.barplot(x='Gene',
                    y=y,
                    data=data,
                    linewidth=2,
                    edgecolor='black',
                    palette=reversed(colours)
                    )

    plt.xticks(rotation=90)
    plt.xlabel('')
    plt.ylabel('% < {}'.format(PVAL))

    seaborn.despine(fig=fig, top=True, right=True)
    if fname is not None:
        fig.savefig(fname, dpi=350, bbox_inches='tight')
    plt.show()
    # print(q2_data)


between_control_adult_fname = os.path.join(saved_objects_path, 'between_control_adult.png')
between_control_senescent_fname = os.path.join(saved_objects_path, 'between_control_sen.png')
between_tgfb_adult_fname = os.path.join(saved_objects_path, 'between_tgf_adult.png')
between_tgfb_senescent_fname = os.path.join(saved_objects_path, 'between_tgf_sen.png')

plot(between_control_data, 'ad_perc', fname=between_control_adult_fname)
plot(between_control_data, 'sen_perc', fname=between_control_senescent_fname)
plot(between_tgfb_data, 'ad_perc', fname=between_tgfb_adult_fname)
plot(between_tgfb_data, 'sen_perc', fname=between_tgfb_senescent_fname)



## Similarly, do the 'within' data plots

within_neonatal_path = os.path.join(saved_objects_path, 'within_neonatal.csv')
within_adult_path = os.path.join(saved_objects_path, 'within_adult.csv')
within_senescent_path = os.path.join(saved_objects_path, 'within_senescent.csv')

within_neonatal_data = pandas.read_csv(within_neonatal_path, index_col=0)
within_adult_data = pandas.read_csv(within_adult_path, index_col=0)
within_senescent_data = pandas.read_csv(within_senescent_path, index_col=0)

within_neonatal_data.index.name = 'Gene'
within_adult_data.index.name = 'Gene'
within_senescent_data.index.name = 'Gene'

within_neonatal_data.reset_index(inplace=True)
within_adult_data.reset_index(inplace=True)
within_senescent_data.reset_index(inplace=True)

# print(within_adult_data)
def plot(data, y=None, fname=None):
    if y is None:
        print(data.head())
        raise ValueError('y cannot be None')

    import matplotlib
    cmap = matplotlib.cm.get_cmap('gist_rainbow')
    colours = []
    count = 0
    num_colours_needed = data.shape[0]
    for i in range(num_colours_needed):
        colours.append(cmap(count))
        count += 1.0 / num_colours_needed

    fig, ax = plt.subplots(figsize=(30, 10))

    seaborn.barplot(x='Gene',
                    y=y,
                    data=data,
                    linewidth=2,
                    edgecolor='black',
                    palette=reversed(colours)
                    )

    plt.xticks(rotation=90)
    plt.xlabel('')
    plt.ylabel('Percentage')

    seaborn.despine(fig=fig, top=True, right=True)
    if fname is not None:
        fig.savefig(fname, dpi=350, bbox_inches='tight')
    plt.show()
    # print(q2_data)


within_neonatal_control_fname = os.path.join(saved_objects_path, 'within_neonatal_control.png')
within_neonatal_tgfb_fname = os.path.join(saved_objects_path, 'within_neonatal_tgfb.png')


within_adult_control_fname = os.path.join(saved_objects_path, 'within_adult_control.png')
within_adult_tgfb_fname = os.path.join(saved_objects_path, 'within_adult_tgfb.png')


within_senescent_control_fname = os.path.join(saved_objects_path, 'within_senescent_control.png')
within_senescent_tgfb_fname = os.path.join(saved_objects_path, 'within_senescent_tgfb.png')



plot(within_neonatal_data, 'tgfb_perc', fname=within_neonatal_control_fname)
plot(within_neonatal_data, 'ctrl_perc', fname=within_neonatal_tgfb_fname)

plot(within_adult_data, 'tgfb_perc', fname=within_adult_control_fname)
plot(within_adult_data, 'ctrl_perc', fname=within_adult_tgfb_fname)

plot(within_senescent_data, 'tgfb_perc', fname=within_senescent_control_fname)
plot(within_senescent_data, 'ctrl_perc', fname=within_senescent_tgfb_fname)

import click
import pandas as pd
import numpy as np
import os
import pickle
import glob

import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("white")
sns.set_context("talk")


from dreem import logger

# universal logger
log = logger.log

def get_trace(df):
    xaxis_coordinates = [i for i in range(1, len(df) + 1)]
    colors = []
    ref_bases = []
    data = []
    cmap = {"A": "red", "T": "green", "G": "orange", "C": "blue"}  # Color map
    for i, row in df.iterrows():
        colors.append(cmap[row['nuc']])
        ref_bases.append(row['nuc'])
        data.append(row['mismatches'])
    trace = go.Bar(
            x=xaxis_coordinates,
            y=data,
            text=ref_bases,
            marker=dict(color=colors),
            showlegend=False,
    )
    return trace


def normalize_reactivity(df):
    df_sub = df[(df['nuc'] == 'A') | (df['nuc'] == 'C') ]
    avg = df['mismatches'].mean()
    df['mismatches'] = df['mismatches'] / avg

def calc_avg(df):
    df = df[20:-20]
    df_sub = df[(df['nuc'] == 'A') | (df['nuc'] == 'C')]
    return df['mismatches'].mean()


@click.group()
def cli():
    pass


@cli.command()
@click.option('--normalize', is_flag=True)
@click.argument("dirs")
def concat(**args):
    spl = args['dirs'].split(",")
    traces = []
    names = []
    for s in spl:
        if not os.path.isdir(s):
            print(s + " is not a directory!")
            exit()
        p = s + "/output/BitVector_Files"
        if not os.path.isdir(p):
            print(p + " does not exist! Must be a valid dreem directory")
            exit()
        mut_histos = {}
        with open(p + "/mutation_histos.p", "rb") as handle:
            mut_histos = pickle.load(handle)
        keys = list(mut_histos.keys())
        keys.sort()
        for k in keys:
            mh = mut_histos[k]
            df = mh.to_pop_avg_data_frame()
            traces.append(get_trace(df))
            names.append(mh.name)
    mut_fig = make_subplots(
            rows=len(traces), cols=1, subplot_titles=names)
    for i, t in enumerate(traces):
        mut_fig.append_trace(t, i + 1, 1)
    mut_fig.update_layout(height=200 * len(traces))
    plotly.offline.plot(
            mut_fig, filename="test.html", auto_open=False
    )

@cli.command()
@click.option('--normalize', is_flag=True)
@click.argument("base_dir")
@click.argument("name")
def scan(**args):
    d = args['base_dir']
    name = args['name']
    traces = []
    names = []
    files = glob.glob(f"{d}/*/output/BitVector_Files/{name}_*.csv")
    files.sort()
    for file in files:
        df = pd.read_csv(file)
        normalize_reactivity(df)
        names.append(file)
        traces.append(get_trace(df))
    mut_fig = make_subplots(
            rows=len(traces), cols=1, subplot_titles=names)
    for i, t in enumerate(traces):
        mut_fig.append_trace(t, i + 1, 1)
    mut_fig.update_layout(height=200 * len(traces))
    plotly.offline.plot(
            mut_fig, filename="test.html", auto_open=False
    )

@cli.command()
@click.option('--normalize', is_flag=True)
@click.argument("base_dir")
@click.argument("name")
@click.argument("res")
def titration(**args):
    d = args['base_dir']
    name = args['name']
    res = [int(x) for x in args['res'].split(",")]
    files = glob.glob(f"{d}/*/output/BitVector_Files/{name}_*.csv")
    files.sort()
    data = []
    for file in files:
        df = pd.read_csv(file)
        #max_v = df["mismatches"].mean()
        normalize_reactivity(df)
        max_v = df[17:20]["mismatches"].mean()
        cur_data = []
        for r in res:
            cur_data.append((df.loc[r-1]["mismatches"]/max_v))
        data.append(sum(cur_data)/len(res))
    data = np.array(data)
    min_d = np.min(data)
    #data -= min_d
    #max_d = np.max(data)
    max_d = 1
    data /= max_d
    for d in data:
        print(d)

@cli.command()
@click.argument("base_dir")
@click.argument('name')
def temptitration(**args):
    d = args['base_dir']
    name = args['name']
    files = glob.glob(f"{d}/*/output/BitVector_Files/{name}_*.csv")
    files.sort()
    temps = []
    for f in files:
        spl = f.split("/")
        spl2 = spl[-4].split("_")
        temp = spl2[-1][:-1]
        temps.append(float(temp))
    all_data = []
    seq = None
    avgs = []
    for file in files:
        df = pd.read_csv(file)
        if seq is None:
            seq = "".join(list(df['nuc'])[20:-20])
        # max_v = df["mismatches"].mean()
        #normalize_reactivity(df)
        max_v = df[20:-20]["mismatches"].mean()
        data = list(df[20:-20]['mismatches'] / 1)
        #print(file, np.mean(data))
        avgs.append(calc_avg(df))
        all_data.append(data)

    min_a, max_a = min(avgs), max(avgs)
    scale = [x / min_a for x in avgs ]

    if not os.path.isdir('plots'):
        os.mkdir('plots')
    for i in range(0, len(seq)):
        titration = []
        for j, data in enumerate(all_data):
            titration.append(data[i]/scale[j])
        plt.plot(temps, titration)
        plt.savefig(f'plots/{i+21}{seq[i]}.png')
        plt.clf()


if __name__ == "__main__":
    cli()


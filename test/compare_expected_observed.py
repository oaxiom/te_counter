
import numpy
import matplotlib.pyplot as plot

def open_mat(filename):
    tab = []

    with open(filename, 'rt') as oh:
        for line in oh:
            if 'name' in line:
                continue
            t = [int(i) for i in line.strip().split('\t')[1:]]
            tab.append(t)
    np = numpy.array(tab)
    np = np.T

    return np

obs = open_mat('run_results/single_cell_out.tsv')
exp = open_mat('run_results/single_cell_out-expected.tsv')

fig = plot.figure()
ax = fig.add_subplot(111)

for p1, p2 in zip(obs, exp):
    ax.plot(p1, p2)

fig.savefig('single_cell_out-plot.pdf')

obs = open_mat('run_results/single_cell_strand_out.tsv')
exp = open_mat('run_results/single_cell_strand_out-expected.tsv')

fig = plot.figure()
ax = fig.add_subplot(111)

for p1, p2 in zip(obs, exp):
    ax.plot(p1, p2)

fig.savefig('single_cell_strand_out-plot.pdf')

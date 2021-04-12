
import numpy
import matplotlib.pyplot as plot

def open_mat(filename):
    tab = []

    oh = open(filename, 'rt')
    for line in oh:
        if 'name' in line:
            continue
        t = [int(i) for i in line.strip().split('\t')[1:]]
        tab.append(t)
    oh.close()

    np = numpy.array(tab)
    np = np.T

    return np

obs = open_mat('single_cell_out.tsv')
exp = open_mat('single_cell_out-expected.tsv')

fig = plot.figure()
ax = fig.add_subplot(111)

for p1, p2 in zip(obs, exp):
    ax.plot(p1, p2)

fig.savefig('single_cell_out-plot.pdf')

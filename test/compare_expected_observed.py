
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

    np = numpy.mean(np, axis=1)

    return np

print('Loading')
obs = open_mat('run_results/single_cell_out.tsv')
exp = open_mat('run_results/single_cell_out-expected.tsv')

print('Plotting')
fig = plot.figure()
ax = fig.add_subplot(111)

ax.scatter(obs, exp, alpha=0.1)
ax.set_xlabel('Observed')
ax.set_ylabel('Expected')

fig.savefig('single_cell_out-plot.pdf')

print('Loading')
obs = open_mat('run_results/single_cell_strand_out.tsv')
exp = open_mat('run_results/single_cell_strand_out-expected.tsv')

fig = plot.figure()
ax = fig.add_subplot(111)

print('Plotting')
ax.scatter(obs, exp, alpha=0.1)
ax.set_xlabel('Observed')
ax.set_ylabel('Expected')
fig.savefig('single_cell_strand_out-plot.pdf')

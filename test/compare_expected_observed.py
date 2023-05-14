
import numpy, gzip
import matplotlib.pyplot as plot

def open_mat(filename, gzip_file=False):
    tab = []

    if gzip_file:
        fh = gzip.open(filename, 'rt')
    else:
        fh = open(filename, 'rt')

    for line in fh:
        if 'name' in line:
            continue
        t = [int(i) for i in line.strip().split('\t')[1:]]
        tab.append(t)

    fh.close()

    np = numpy.array(tab)
    np = np.T

    np = numpy.mean(np, axis=1)

    return np

print('Loading')
obs = open_mat('run_results/single_cell_out.tsv')
exp = open_mat('run_results/single_cell_out-expected.tsv.gz', True)

print('Plotting')
fig = plot.figure()
ax = fig.add_subplot(111)

ax.scatter(obs, exp, alpha=0.1)
ax.set_xlabel('Observed')
ax.set_ylabel('Expected')

fig.savefig('single_cell_out-plot.pdf')

print('Loading')
obs = open_mat('run_results/single_cell_strand_out.tsv')
exp = open_mat('run_results/single_cell_strand_out-expected.tsv.gz', True)

fig = plot.figure()
ax = fig.add_subplot(111)

print('Plotting')
ax.scatter(obs, exp, alpha=0.1)
ax.set_xlabel('Observed')
ax.set_ylabel('Expected')
fig.savefig('single_cell_strand_out-plot.pdf')

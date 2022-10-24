"""
Think of me as a delayed version of geneList

* this is the only special case of geneList (delayed)
"""
import sys, os, time, copy, csv, gzip

from . import config
from . import utils
from .genelist import genelist
from .location import location

ignorekeys = frozenset( # these are functional tags - so I should ignore them.
    ["dialect",
    "duplicates_key",
    "skiplines",
    "debug",
    "special",
    "skiptill",
    "force_tsv",
    "gtf_decorators",
    "endwith",
    "__description__",
    "commentlines",
    "keepifxin",
    '__column_must_be_used',
    '__ignore_empty_columns'
    ])

typical_headers = frozenset(["chipseq_loc", "loc", "chr", "#",
    "Gene Name", "", "GenBank", "RefSeq",
    "Systematic", "mm8.refGene.chrom", "mm8", "loc", 'chromosome',
    "mm9.refGene.chrom", "mm9",
    "======================================================================", # stupid sissrs format garbage.
    "=====================================================================", # stupid sissrs format garbage.
    "======================================================================="] # stupid sissrs format garbage.
    ) # typical header labels;

class delayedlist(genelist):
    """
    **Purpose**

        delayedlist is very similar to a genelist - except the data never
        makes it into memory. This allows you to iterate over huge files from disk.

    **Arguments**
        name (Optional)
            Will default to the name of the file if a file is loaded,
            otherwise name will be set to "Generic List" by default.
            us the name argument to give it a custom nam,e

        filename (Optional)
            the name of a file to attempt to load.

        force_tsv (Optional)
            if you specify a filename to load then
            setting this argument to True will force the file to be treated
            as a tsv (tab-separated) rather than the default csv (comma
            separated).

        gzip (Optional, default=False)
            The input file is a gzip file.

        format
            format specifier. (see docs... complex)

    """

    def __init__(self, filename=None, format=None, force_tsv=False, gzip=False, **kargs):
        genelist.__init__(self) # no kargs though. I want the mpty version.

        self.__len_estimate = None

        assert filename, "No Filename"
        assert os.path.exists(filename), "%s not found" % (filename)
        assert format, "You must provide a format for delayedlist. I cannot guess its format."

        self.path = os.path.split(os.path.realpath(filename))[0]
        self.filename = os.path.split(os.path.realpath(filename))[1]
        self.fullpath = filename
        self.filehandle = None
        self.format = format # override default
        self.gzip = gzip

        if force_tsv:
            self.format["force_tsv"] = True

        config.log.info("Bound '%s' as a delayedlist" % filename)
        self._optimiseData()

    def __repr__(self):
        return("glbase.delayedlist")

    def __len__(self):
        # I need to collect an estimate
        if not self.__len_estimate:
            lines = 0
            if not self.gzip:
                f = open(self.fullpath, 'rb')
                for _ in f: lines += 1

            else: # gzipped file variant
                f = gzip.open(self.fullpath, 'rb') # must be rb :(
                for _ in f.readlines(): lines += 1
            self.__len_estimate = lines-1 # start from 0

        return(self.__len_estimate)

    def __getitem__(self, index):
        """
        (Override)
        confers a = geneList[0] behaviour
        This is broken. It returns only the very first entry,
        whatever 'index' is sent to it, it disregards it.
        This only continues to exist for compatability with some internal
        routines.
        """
        self._optimiseData()
        return(next(self.__iter__()))

    def __iter__(self):
        """
        (Override)
        make the geneList behave like a normal iterator (list)
        """
        try:
            column = next(self.__reader) # get started
            self.cindex += 1

            while column:
                d = None
                while not d:
                    if column: # list is empty, so omit.
                        if "commentlines" in self.format:
                            if column[0][0] == self.format["commentlines"]: # csv will split the table and returns a list
                                column = None # force a skip of this row, don't use continue it will just hang
                                d = None
                            elif column[0].startswith(self.format["commentlines"]):
                                column = None
                                d = None

                    if column:
                        if "keepifxin" in self.format:
                            if True in [self.format["keepifxin"] in i for i in column]:
                                if (not (column[0] in typical_headers)):
                                    d = self._processKey(self.format, column)
                            else:
                                d = None # not present, skip this line

                        else: # just do normally
                            if (not (column[0] in typical_headers)):
                                d = self._processKey(self.format, column)

                    if not d: # d is bad, grab another
                        column = next(self.__reader)
                        self.cindex += 1

                # I do quoting = NONE now, so I need to manually deal with containing quotes.
                for k in d: # d must be valid to get here
                    if isinstance(d[k], str): # only for strings
                        if d[k][0] == "'" and d[k][-1] == "'":
                            d[k] = d[k].strip("'")
                        if d[k][0] == '"' and d[k][-1] == '"':
                            d[k] = d[k].strip('"')

                yield d # d should be valid
                column = next(self.__reader) # get the next item
                self.cindex += 1

        except StopIteration:
            self._optimiseData()
            return # py3.7 new

    def _optimiseData(self):
        """
        (Override)
        Impossible to optimise the data.
        so we just reset the entry point.
        This makes the iterator work like you would expect:
        a new iterator will go back to the beginning of the list.
        """
        if self.filehandle:
            self.filehandle.close()

        if not self.gzip:
            self.filehandle = open(self.fullpath, "rt")
        else:
            self.filehandle = gzip.open(self.fullpath, 'rt') # must be rb :(

        if "force_tsv" in self.format and self.format["force_tsv"]:
            self.__reader = csv.reader(self.filehandle, dialect=csv.excel_tab, quoting=csv.QUOTE_NONE)
        elif "dialect" in self.format:
            self.__reader = csv.reader(self.filehandle, dialect=self.format["dialect"], quoting=csv.QUOTE_NONE)
        else:
            self.__reader = csv.reader(self.filehandle, quoting=csv.QUOTE_NONE)

        if "skiplines" in self.format:
            if self.format["skiplines"] != -1: # no skipped lines, good to go.
                for i, x in enumerate(self.__reader):
                    if i == self.format["skiplines"]:
                        break
        else: # default behaviour of genelist is to always skip the first line
            next(self.__reader)

        if "skiptill" in self.format:
            done = False
            while not done:
                for line in self.__reader:
                    if self.format["skiptill"] in line:
                        done = True
                        break

        self.linearData = self.__iter__()
        self.cindex = 0
        return True

    def __str__(self):
        self._optimiseData()

        # load in a bunch of data and dump it into self.linearData
        temp_data = []
        for index, item in enumerate(self):
            temp_data.append(item)
            if index > config.NUM_ITEMS_TO_PRINT-2:
                break # only get the first n data.s
        self.linearData = temp_data
        ret = genelist.__str__(self)
        self._optimiseData()
        return("%s\nThis is a delayedlist - only the first %s entries are shown" %(ret, config.NUM_ITEMS_TO_PRINT))

    def reset(self):
        """
        **Purpose**
            reset the delayedlist to the 0th element.
            This is a bit of a hack. for most purposes
            delayedlist will be correctly reset. The exception is this case:

            for item in delayedlist:
                ... some code
                break

            for item in delayedlist:
                !!! Error! The list
                continues from where it left off, not from the zeroth
                element as expected.

            Actually, (not tested) I think iterating over a delayedlist
            twice is generally broken, and you should reset.
            However, there is no way for delayedlist to know if the next
            iteration is actually the first iteration of the list
            and not a continuing iteration.

        **Arguments**
            None

        **Results**
            resets the list to the zeroth entry.
        """
        self._optimiseData()

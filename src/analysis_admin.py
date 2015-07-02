#!/usr/bin/env python
# Authors: Gao Wang <gaow@uchicago.edu>
# Created June 03, 2015
"""
Administrative program for data (GTEx) analysis
"""
import os, sys, argparse, random, gzip, re, warnings, copy, math
from utils import env, get_tb_grps, get_gs_pairs
from utils import TBData as SSData
import tables as tb
import numpy as np
from subprocess import check_output
from pandas import DataFrame, concat

def get_batch(nlines, nbatches, batch_id):
    '''Note: 1. The last batch may have smaller number of entries
    2. It is a little wasteful to compute the entire batch index each time but only take
    one batch. But it is trivial compare to the eqtlbma step'''
    per_batch = int(nlines / nbatches + 1)
    batches = {}
    for idx, value in enumerate(range(1, nlines, per_batch)):
        batches[idx + 1] = (value, min(value + per_batch - 1, nlines))
    return batches[batch_id]

class CoordChopper:
    def __init__(self, gene_coords, snp_coords, start, end, batch_idstr):
        self.gene_coords = gene_coords
        self.snp_coords = snp_coords
        self.start = start
        self.end = end
        self.gf = '%s/genes_%s.bed.gz' % (env.batch_file_dir, batch_idstr)
        self.sf = '%s/snps_%s.bed.gz' % (env.batch_file_dir, batch_idstr)

    def ChopGeneCoords(self):
        env.log('Creating gene coords file ...')
        os.system('zcat {} | sed -n {},{}p | bgzip > {}'.\
                     format(self.gene_coords, self.start, self.end, self.gf))
        os.system('tabix -p bed {}'.format(self.gf))

    def ChopSNPCoords(self, window):
        env.log('Creating SNP coords file ...')
        os.system('bedtools window -w {} -a {} -b {} | cut -f7,8,9,10 | sort -k1,1g -k2,2g | bgzip > {}'.\
                  format(window, self.gf, self.snp_coords, self.sf))
        os.system('tabix -p bed {}'.format(self.sf))

    def Clean(self):
        os.system('rm -f %s*' % self.gf)
        os.system('rm -f %s*' % self.sf)
        
class EQTLBMARunner:
    def __init__(self, args_file, exe_path, seed, batch_id, batch_idstr):
        self.args = ' '.join([line.strip('\n').strip('\\') for line in open(args_file).readlines()]).split()
        self.exe = exe_path
        if seed:
            random.seed(seed + batch_id)
            self.seed = random.randint(1,100000)
        else:
            self.seed = None
        self.name = batch_idstr

    def GenerateCommand(self, gf, sf):
        try:
            self.args[self.args.index('--gcoord') + 1] = gf
        except ValueError:
            self.args.extend(['--gcoord', gf])
        try:
            self.args[self.args.index('--scoord') + 1] = sf
        except ValueError:
            self.args.extend(['--scoord', sf])
        #
        try:
            try:
                base, ext = self.args[self.args.index('--out') + 1].split(os.extsep, 1)
            except:
                base = self.args[self.args.index('--out') + 1]
                ext = ''
            self.name = '{}_{}'.format(base, self.name)
            self.args[self.args.index('--out') + 1] = '{0}/{0}{1}'.\
              format(self.name, ('.' + ext) if ext else '')
        except ValueError:
            raise ValueError('Please specify "--out" in argument file!')
        # 
        if self.seed:
            try:
                self.args[self.args.index('--seed') + 1] = str(self.seed)
            except ValueError:
                self.args.extend(['--seed', str(self.seed)])
        return ' '.join(self.args)

    def Run(self, cmd):
        # simply use os.system for now. will see if there is a need to use subprocess call
        env.log('Running eqtlbma ...')
        os.system('mkdir -p {}'.format(self.name))
        os.system('echo {} > {}.log'.format(' '.join(sys.argv), self.name))
        os.system('{0} {1} >> {2}.log 2> {2}.err'.format(self.exe, cmd, self.name))

    def Print(self, cmd):
        print('{} {}'.format(self.exe, cmd))

class SSDataParser:
    def __init__(self, header = False):
        self.data = {'buffer':{'data':[], 'rownames':[]}, 'output':{}}
        self.header = header
        self.previous_name = self.current_name = None
        self.count = -1

    def parse(self, line):
        # input line is snp, gene, beta, t, pval
        if not line:
            self.__reset()
            self.current_name = None
            return 1
        line = line.strip().split()
        self.count += 1
        if self.header and self.count == 0:
            return 0
        #
        if self.previous_name is None:
            self.previous_name = line[1]
        self.current_name = line[1]
        if self.current_name != self.previous_name:
            self.__reset()
        self.data['buffer']['data'].append([line[2], line[3], line[4]])
        self.data['buffer']['rownames'].append(line[0])
        return 0

    def __reset(self):
        self.data['buffer']['data'] = np.array(self.data['buffer']['data'], dtype = env.float)
        self.data['buffer']['rownames'] = np.array(self.data['buffer']['rownames'])
        self.data['buffer']['colnames'] = np.array(['beta','t-stat','p-value'])
        self.data['output'] = copy.deepcopy(self.data['buffer'])
        self.data['buffer'] = {'data':[], 'rownames':[]}
        
    def dump(self):
        return self.data['output']
    
def eqtlbma_batch(args):
    os.system('mkdir -p %s' % env.batch_file_dir)
    n_tss = int(check_output('zcat {} | wc -l'.format(args.gene_coords), shell = True))
    start, end = get_batch(n_tss, args.n_batches, args.batch_id)
    batch_idstr = '{}_{}_{}'.format(args.n_batches, start, end)
    cc = CoordChopper(args.gene_coords, args.snp_coords, start, end, batch_idstr)
    cc.ChopGeneCoords()
    cc.ChopSNPCoords(args.window)
    er = EQTLBMARunner(args.args_file, args.eqtlbma_path, args.seed, args.batch_id, batch_idstr)
    cmd = er.GenerateCommand(cc.gf, cc.sf)
    if not args.dry_run:
        er.Run(cmd)
        if args.clean:
            cc.Clean()
    else:
        er.Print(cmd)
    return

class SSDataMerger(SSData):
    def __init__(self, files, name, msg = None):
        SSData.__init__(self, {}, name, msg)
        self.files = sorted(files)
        self.__group = name

    def merge(self):
        data = {}
        one_snp = None
        failure_ct = 0
        # Collect data
        for item in self.files:
            tissue = re.sub(r'{}$'.format(env.common_suffix), '', os.path.basename(item))
            try:
                data[tissue] = SSData(item, self.__group)
                if one_snp is None: one_snp = data[tissue]['rownames'][0]
            except ValueError:
                data[tissue] = {'data' : np.array([[np.nan, np.nan, np.nan]]), 'rownames': None}
                failure_ct += 1
            # Fix row name
            # Because in GTEx data file there are duplicated gene-snp pairs having different sumstats!!
            if data[tissue]['rownames'] is not None:
                data[tissue]['rownames'] = self.__dedup(data[tissue]['rownames'], item)
        if failure_ct == len(self.files):
            return 1
        # Merge data
        for idx, item in enumerate(['beta','t-stat','p-value']):
            self[item] = concat([DataFrame(
                {tissue : data[tissue]['data'][:,idx]},
                index = data[tissue]['rownames'] if data[tissue]['rownames'] is not None else [one_snp]
                ) for tissue in sorted(data.keys())], axis = 1)
            if 'rownames' not in self:
                self['rownames'] = np.array(self[item].index, dtype = str)
            if 'colnames' not in self:
                self['colnames'] = np.array(self[item].columns.values.tolist(), dtype = str)
            self[item] = np.array(self[item].as_matrix(), dtype = env.float)
        # np.savetxt(sys.stdout, self['p-value'], fmt='%10.5f')
        # print(self['rownames'])
        # print(self['colnames'])
        return 0

    def __dedup(self, seq, filename):
        seen = {}
        dups = set()
        def __is_seen(x, seen):
            if x not in seen:
                seen[x] = 0
                return 0
            else:
                seen[x] += 1
                dups.add(x)
                return 1
        # Tag them
        obs = [x if not __is_seen(x, seen) else '%s%s%s' % (x, env.duplicate_tag, seen[x]) for x in seq]
        # Log them
        if len(dups):
            filename = os.path.splitext(filename)[0]
            for item in dups: 
                env.error('{}:{} appeared {} times in {}'.\
                          format(self.__group, item, seen[item] + 1, filename),
                          to_file = filename + '.error')
        return obs
    
def ss_to_h5(args):
    if args.action == 'convert':
        # Benchmark:
        # 118min to convert a 5.0G gz sumstats file with 179083485 lines to 3.8G h5 file
        # when TBData format bzip2 compression is used and env.float = numpy.float32
        # file size is 5.1G when env.float=numpy.float64
        # 138min to convert to 4.9G h5 when zlib is used instead of bzip2
        # 39min to convert the same file to 6.9G h5 file
        # when HPData format is used with zlib and env.float = numpy.float64
        for item in args.input:
            ssp = SSDataParser(header = True)
            group_counts = 0
            bname = os.path.basename(item)
            with gzip.open(item) as f:
                while True:
                    line = f.readline()
                    quit = ssp.parse(line)
                    if ssp.current_name != ssp.previous_name:
                        group_counts += 1.0
                        with warnings.catch_warnings():
                            warnings.filterwarnings("ignore", category = tb.NaturalNameWarning)
                            # warnings.filterwarnings("ignore", category = tb.PerformanceWarning)
                            data = SSData(ssp.dump(), ssp.previous_name, args.message)
                            data.sink(os.path.join(args.output, (bname.split(os.extsep, 1)[0] + ('_%i.h5' % (math.ceil(group_counts / args.maxsize)) if args.maxsize else '.h5'))))
                        ssp.previous_name = ssp.current_name
                    if quit:
                        env.log("[%s] Processed %s lines\n" % (bname, ssp.count), flush = True)
                        break
                    if ssp.count % env.batch['lines'] == 0:
                        env.log("[%s] Processed %s lines" % (bname, ssp.count), flush = True)
    if args.action == 'merge':
        if args.gene_list:
            gene_names = sorted(set([line.strip('\n') for line in open(args.gene_list).readlines()]))
        else:
            gene_names = get_tb_grps(args.input)
        failure_ct = 0
        for idx, item in enumerate(gene_names):
            env.log("Merging group #%s ..." % (idx + 1), flush = True)
            ssm = SSDataMerger(args.input, item, args.message)
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category = tb.NaturalNameWarning)
                if ssm.merge() == 0:
                    ssm.sink(args.output + ('_%i.h5' % (math.ceil((idx + 1.0) / args.maxsize)) if args.maxsize else '.h5'))
                else:
                    failure_ct += 1
        env.log("%s groups merged!\n" % (idx + 1 - failure_ct), flush = True)
    if args.action == 'cat':
        for item in sorted(args.input):
            for name in get_tb_grps([item]):
                os.system('h5copy -i {0} -o {1} -s "/{2}" -d "/{2}"'.format(item, args.output, name))
    if args.action == 'summary':
        names = get_tb_grps(args.input)
        for item in names:
            print(item)
    if args.action in ['max', 'null']:
        assert len(args.input) == 1
        if args.gene_list:
            names = sorted(set([line.strip('\n') for line in open(args.gene_list).readlines()]))
        else:
            names = get_tb_grps(args.input)
        data = {'colnames': None, 'rownames': [], 'beta': None, 't-stat': None, 'p-value': None}
        failure_ct = 0
        for idx, name in enumerate(names):
            env.log("Calculating #%s ..." % (idx + 1), flush = True)
            # extract the best gene-snp pair or some null gene-snp pairs
            best = get_gs_pairs(SSData(args.input[0], name), name,
                                0 if args.action == 'max' else env.nb_null_pairs)
            #
            if best is None:
                failure_ct += 1
                continue
            for k in data:
                if data[k] is None:
                    data[k] = best[k]
                    continue
                if k == 'colnames':
                    continue
                elif k == 'rownames':
                    data[k].extend(best[k])
                else:
                    data[k] = np.vstack((data[k], best[k]))
        #
        if failure_ct < len(names):
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category = tb.NaturalNameWarning)
                SSData(dict(data), args.action, msg = "%s gene-snp pair, GTEx v6" % args.action).\
                  sink(args.output)
        env.log("%s groups processed!\n" % (idx + 1 - failure_ct), flush = True)

def sumstat_query(args):
    env.log("Loading data ...")
    try:
        data = SSData(args.data, args.gene)
    except Exception as e:
        env.error(e, exit=1)
    if not args.merged:
        df = data.dump('data', True)
        try:
            print(df.loc[args.snp])
        except KeyError:
            env.error("{} not found in {}!".format(args.snp, args.gene))
            env.log("Available SNP names are:")
            for item in data['rownames']:
                print(item)

if __name__ == '__main__':
    def uint(value):
        ivalue = int(value)
        if ivalue < 0:
            raise argparse.ArgumentTypeError("%s is an invalid positive int value" % value)
        return ivalue
        
    def pathstr(value):
        return os.path.abspath(os.path.expanduser(value))
    #
    parser = argparse.ArgumentParser(description = __doc__)
    subparsers = parser.add_subparsers()
    p = subparsers.add_parser('eqtlbma_batch', help = 'Divide data into gene batches and run eqtlbma',
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-g', dest = 'gene_coords', type = pathstr, required = True,
                   help = 'Gene (TSS) coordinate file, in bed.gz format.')
    p.add_argument('-s', dest = 'snp_coords', type = pathstr, required = True,
                   help = 'SNP coordinate file, in bed.gz format.')
    p.add_argument('-w', metavar = 'N', dest = 'window', type = uint, default = 100000, help = 'Window size.')
    p.add_argument('-n', dest = 'n_batches', type = uint,
                   default = 1000, help = 'Total number of batches.')
    p.add_argument('-b', dest = 'batch_id', type = uint, default = 1,
                   help = 'Execute the i-th batch. Program will quit if invalid ID is provided.')
    p.add_argument('-a', dest = 'args_file', type = pathstr, required = True,
                   help = 'Path to file containing additional eqtlbma command arguments.')
    p.add_argument('--seed', metavar = 'N', type = uint,
                   help = '''If specified, a random number will be generated using (N + batch_id) as seed,
                   and the eqtlbma command will be appended a "--seed" argument with the number generated
                   here.''')
    p.add_argument('-e', dest = 'eqtlbma_path', type = pathstr, required = True,
                   help = 'Path to an eqtlbma_* executable.')
    p.add_argument('--dry-run', dest = 'dry_run', action = 'store_true',
                   help = 'Only generate and save the batch data & commands without performing analysis.')
    p.add_argument('--clean', dest = 'clean', action = 'store_true',
                   help = 'Remove batch gene / snp coords file upon finishing the analysis.')
    p.set_defaults(func=eqtlbma_batch)
    #
    p = subparsers.add_parser('ss_to_h5', help = 'Convert summary statistics output to HDF5 format',
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('input', nargs = '+', type = pathstr, help = 'Input data files.')
    p.add_argument('--action', required = True,
                   choices = ['convert', 'merge', 'cat', 'summary', 'max', 'null'],
                   help = '''Convert: deals with converting summary statistic text files to HDF5 file;
                   Merge: merging multiple HDF5 files to one
                   (both for all data as well as for most significant SNP per gene data);
                   Cat: concatenate multiple HDF5 files;
                   Summary: output unique groups in an HDF5 file;
                   max: output the "best" gene-snp pair;
                   null: output %s null gene-snp pairs''' % env.nb_null_pairs)
    p.add_argument('--maxsize', type = uint, help = 'Maximum number of groups per HDF5 file.')
    p.add_argument('--output', required = True, type = pathstr, help = 'Output data dir / base name.')
    p.add_argument('--gene-list', dest = 'gene_list', type = pathstr, help = 'Path to gene list file. If specified, genes for the merger will be read from this list rather than from input data (which will be much slower).')
    p.add_argument('--message', help = 'A message string of data description.')
    p.set_defaults(func=ss_to_h5)
    p = subparsers.add_parser('sumstat_query',
                              help = 'Query information from sumstats file in HDF5 format',
                              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('data', type = pathstr, help = 'Input data file name')
    p.add_argument('-g', dest = 'gene', help = 'Gene name')
    p.add_argument('-s', dest = 'snp', help = 'SNP id')
    p.add_argument('--merged', action = 'store_true', help = 'Input data has been merged across tissues.')
    p.set_defaults(func=sumstat_query)
    args = parser.parse_args()
    try:
        args.func(args)
    except Exception as e:
        raise
        sys.exit(e)

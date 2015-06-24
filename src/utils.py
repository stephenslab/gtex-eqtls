from __future__ import print_function
import sys
import numpy as np
from pandas import DataFrame
import tables as tb
import h5py as hp
tb.parameters.MAX_GROUP_WIDTH = 51200
# tb.parameters.NODE_CACHE_SLOTS = -51200
# tb.parameters.METADATA_CACHE_SIZE = 1048576 * 100000
# tb.parameters.CHUNK_CACHE_SIZE = 2097152 * 100000
# tb.parameters.CHUNK_CACHE_NELMTS = 521

class TBData(dict):
    def __init__(self, data, name, msg = None):
        self.__group = name
        self.__msg = msg
        try:
            if type(data) is dict:
                self.update(data)
            elif type(data) is str:
                # is file name
                self.__load(tb.open_file(data))
            else:
                # is file stream
                self.__load(data)
        except tb.exceptions.NoSuchNodeError:
            raise ValueError('Cannot find dataset {}!'.format(name))
        # self.tb_filters = tb.Filters(complevel=9, complib='bzip2') # bzip2 is not compatible with other hdf5 applications
        self.tb_filters = tb.Filters(complevel=9, complib='zlib')

    def sink(self, filename):
        with tb.open_file(filename, 'a') as f:
            try:
                # there is existing data -- have to merge with current data
                # have to do this because the input file lines are not grouped by gene names!!
                # use try ... except to hopefully faster than if ... else
                # e.g., if not f.__contains__('/{}'.format(self.__group)) ... else ...
                for element in f.list_nodes('/{}'.format(self.__group)):
                    if element.name != 'colnames':
                        self[element.name] = np.concatenate((element[:], self[element.name]))
            except tb.exceptions.NoSuchNodeError:
                f.create_group("/", self.__group, self.__msg if self.__msg else self.__group)
            for key in self:
                self.__store_array(key, f)
            f.flush()

    def dump(self, table, output = False):
        if output:
            DataFrame({self['colnames'][i] : self[table][:,i] for i in range(len(self['colnames']))}, index = self['rownames']).to_csv(sys.stdout, na_rep = 'NA')
            return None
        else:
            return DataFrame({self['colnames'][i] : self[table][:,i] for i in range(len(self['colnames']))}, index = self['rownames'])

    def __load(self, fstream):
        try:
            for element in fstream.list_nodes('/{}'.format(self.__group)):
                self[element.name] = element[:]
            fstream.close()
        except:
            fstream.close()
            raise

    def __roll_back(self, group, name):
        try:
            n = getattr(group, name)
            n._f_remove()
        except AttributeError:
            pass

    def __store_array(self, name, fstream):
        element = getattr(fstream.root, self.__group)
        arr = self[name]
        if type(arr) is list:
            arr = np.array(arr)
        self.__roll_back(element, name)
        #
        if arr.shape != (0,):
            ds = fstream.create_carray(element, name, tb.Atom.from_dtype(arr.dtype), arr.shape,
                                       filters = self.tb_filters)
            ds[:] = arr

class HPData(dict):
    def __init__(self, data, name, msg = None):
        self.__group = name
        self.__msg = msg
        try:
            if type(data) is dict:
                self.update(data)
            elif type(data) is str:
                # is file name
                self.__load(hp.File(data))
            else:
                # is file stream
                self.__load(data)
        except Exception as e:
            raise ValueError(e)
        self.hp_filters = dict(compression = 'gzip', compression_opts = 9, shuffle = True, fletcher32 = True)    

    def sink(self, filename):
        with hp.File(filename, 'a') as f:
            if f.__contains__('{}'.format(self.__group)):
                g = f.get(self.__group)
                # there is existing data -- have to merge with current data
                # have to do this because the input file lines are not grouped by gene names!!
                for key, value in g.iteritems():
                    if key != 'colnames':
                        self[key] = np.concatenate((value[:], self[key]))
                # raise ValueError('Group {} already exists in {}!'.format(self.__group, filename))
            else:
                g = f.create_group(self.__group)
            for key in self:
                self.__store_array(key, g)
            f.flush()

    def dump(self, table, output = False):
        if output:
            DataFrame({self['colnames'][i] : self[table][:,i] for i in range(len(self['colnames']))}, index = self['rownames']).to_csv(sys.stdout, na_rep = 'NA')
            return None
        else:
            return DataFrame({self['colnames'][i] : self[table][:,i] for i in range(len(self['colnames']))}, index = self['rownames'])

    def __load(self, fstream):
        if not fstream.__contains__(self.__group):
            fstream.close()
            raise ValueError('Group {} does not exist!'.format(self.__group))
        #
        gob = fstream.get(self.__group)
        for key, value in gob.iteritems():
            self[key] = value[:]
        fstream.close()

    def __store_array(self, name, gob):
        arr = self[name]
        if type(arr) is list:
            arr = np.array(arr)
        if name in gob:
            del gob[name]
        if arr.shape != (0,):
            gob.create_dataset(name, shape = arr.shape, dtype = arr.dtype, data = arr, **self.hp_filters)


class Environment:
    '''Class of "global" variables'''
    def __init__(self):
        self.float = np.float64
        self.debug = False
        self.quiet = False
        self.__width_cache = 1
        self.prog = sys.argv[0]
        self.batch = {'lines': 10000}
        self.common_suffix = '_Analysis.h5'
        self.duplicate_tag = '_duplicated_'
        self.batch_file_dir = 'batch_files'
        
    def error(self, msg = None, show_help = False, exit = False, to_file = None):
        if to_file:
            lf = open(to_file, 'a')
        else:
            lf = sys.stderr
        if msg is None:
            print('\n', file = lf)
            if to_file: lf.close()
            return
        if type(msg) is list:
            msg = ' '.join(map(str, msg))
        else:
            msg = str(msg)
        start = '\n' if msg.startswith('\n') else ''
        end = '\n' if msg.endswith('\n') else ''
        msg = msg.strip()
        print(start + "\033[1;40;33mERROR: {}\033[0m".format(msg) + end, file = lf)
        if to_file: lf.close()
        if show_help:
            self.log("Type '{} -h' for help message".format(env.prog))
            sys.exit()
        if exit:
            sys.exit()
        
    def log(self, msg = None, flush=False):
        if self.debug or self.quiet:
            return
        if msg is None:
            sys.stderr.write('\n')
            return
        if type(msg) is list:
            msg = ' '.join(map(str, msg))
        else:
            msg = str(msg)
        start = "{0:{width}}".format('\r', width = self.__width_cache + 10) + "\r" if flush else ''
        end = '' if flush else '\n'
        start = '\n' + start if msg.startswith('\n') else start
        end = end + '\n' if msg.endswith('\n') else end
        msg = msg.strip()
        sys.stderr.write(start + "\033[1;40;32mMESSAGE: {}\033[0m".format(msg) + end)
        self.__width_cache = len(msg)

env = Environment()

def get_tb_grps(filenames):
    env.log('Collecting group names from input files ...')
    names = set()
    for filename in filenames:
        env.log(filename, flush=True)
        with tb.open_file(filename) as f:
            names.update([node._v_name for node in f.root])
    env.log('%s unique groups identified from %s files\n' % (len(names), len(filenames)), flush = True)
    return sorted(names)

def get_best_pair(data, name):
    output = {'colnames' : data['colnames']}
    # Find max SNP-gene pair
    t = data.dump('t-stat')
    t = t[np.all(np.isfinite(t), axis=1)]
    #
    if t.empty:
        return None
    rowidx = np.where(data['rownames'] == t.abs().max(axis=1).idxmax())[0]
    output['rownames'] = '%s_%s' % (name, data['rownames'][rowidx][0])
    for k in ['beta', 'p-value', 't-stat']:
        output[k] = data[k][rowidx, :][0]
    return output

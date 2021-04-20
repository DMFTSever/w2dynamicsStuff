import h5py as hdf5
import numpy as np


class MissmatchError(Exception):

    def __init__(self, message, filename):
        self.message  = message
        self.filename = filename

    def __str__(self):
        string = "Concerning processing of file: {}\n".format(self.filename)
        string += self.message
        return string


class W2dynfile():

    def __init__(self, filename):

        self.filename = filename
        self.file     = hdf5.File(filename, 'r')
        self.dmftiter = int(self.file['/.config']
                            .attrs.get('general.dmftsteps'))
        self.statiter = int(self.file['/.config']
                                 .attrs.get('general.statisticsteps'))

        if '/finish' in self.file:
            self.lastdmftiter = self.dmftiter
            self.laststatiter = self.statiter
        else:
            for i in range(self.dmftiter, 0, -1):
                if i == 1:
                    print("Found no valid dmft iteration in file: {}"
                          .format(self.filename))
                    print("Here is what is in the file:")
                    self.walk_file()
                    raise RuntimeError('Aborting...')
                elif '/dmft-%03i' % i in self.file:
                    self.lastdmftiter = i-1  # str('%03i'%(i-1))
                    print('Last valid dmft iteration in file {}: {}'
                          .format(self.filename, i-1))
                    break
                else:
                    pass
            for i in range(self.statiter, 0, -1):
                if i == 1:
                    print("Found no valid stat iteration in file: {}"
                          .format(self.filename))
                    print("Here is what is in the file:")
                    self.walk_file()
                    raise RuntimeError('Aborting...')
                elif '/stat-%03i' % i in self.file:
                    self.lastdmftiter = i-1  # str('%03i'%(i-1))
                    print('Last valid stat iteration in file {}: {}'
                          .format(self.filename, i-1))
                    break
                else:
                    pass

        if self.laststatiter != 0:
            self.ref = 'stat-{:03d}'.format(self.laststatiter)
        else:
            self.ref = 'dmft-{:03d}'.format(self.lastdmftiter)

        self.natoms = int(self.file['/.config'].attrs.get('general.nat'))
        self.nineq  = 0

        try:
            for element in self.file['/dmft-001'].keys():
                if "ineq" in element:
                    self.nineq += 1
        except KeyError:
            self.ninew = 0
            for element in self.file['/stat-001'].keys():
                if "ineq" in element:
                    self.nineq += 1

        self.equiv = self.file['/.config'].attrs.get('general.equiv')

        if self.equiv is None and self.natoms != self.nineq:
            raise MissmatchError("Unknown quivalence of atoms! Aborting...",
                                 self.filename)
        elif self.equiv is None:
            self.equiv = np.arange(1, self.natoms+1)
        else:
            self.equiv = self.equiv.astype(int)+1

        self.niw = self.file['/.config'].attrs.get('qmc.niw')
        self.ntau  = self.file['/.config'].attrs.get('qmc.ntau')
        self.n4iwf = self.file['/.config'].attrs.get('qmc.n4iwf')
        self.n4iwb = self.file['/.config'].attrs.get('qmc.n4iwb')
        if self.n4iwb == 0: self.n4iwb = 1
        self.n4tau  = self.file['/.config'].attrs.get('qmc.n4tau')

    def walk_file(self, iterations='all', stat_iter=False):

        def to_int(string):
            try:
                iteration = int(string.strip("dmft-").strip("stat-"))
            except:
                return -1
            return iteration

        if iterations != 'all' and type(iterations) != list:
            raise MissmatchError("iterations need of to be of type list!",
                                 self.filename)

        if stat_iter is False and self.lastdmftiter == 0:
            stat_iter = True

        if stat_iter:
            keys = ['stat', 'dmft']
        else:
            keys = ['dmft', 'stat']

        for key_level1 in self.file.keys():
            if iterations != 'all' and keys[1] in key_level1:
                pass
            elif iterations != 'all' and keys[0] in key_level1 \
                                     and to_int(key_level1) not in iterations:
                pass
            else:
                print("{:20}".format(key_level1), self.file[key_level1])
                try:
                    for key_level2 in self.file[key_level1].keys():
                        print("  {:18}".format(key_level2),
                              self.file[key_level1][key_level2])
                        try:
                            for key_level3 in self.file[key_level1][key_level2]\
                                                  .keys():
                                print("    {:16}".format(key_level3),
                                      self.file[key_level1][key_level2]
                                               [key_level3])
                                try:
                                    for key_level4 in self.file[key_level1][
                                                                key_level2][
                                                                key_level3]\
                                                                .keys():
                                        print("      {:14}".format(key_level4),
                                              self.file[key_level1][key_level2]
                                                       [key_level3][key_level4])
                                except:
                                    pass
                        except:
                            pass
                except:
                    pass
                print("")

        for atr in self.file.attrs.keys():
            print(atr, self.file.attrs.get(atr))

        if iterations != 'all':
            for x in iterations:
                if not stat_iter and int(x) > self.lastdmftiter:
                    print("ITERATION {} LARGER THAN MAXIMAL ITERATION IN FILE"
                          .format(x))
                elif int(x) > self.laststatiter:
                    print("ITERATION {} LARGER THAN MAXIMAL ITERATION IN FILE"
                          .format(x))

    def get_tauax(self, ntau='all'):

        if ntau == 'all':
            ntau = self.ntau
        elif ntau > self.ntau or type(ntau) != int:
            raise MissmatchError("ntau hast to be of type integer < self.ntau",
                                 self.filename)

        return self.file['/.axes/taubin'][0:ntau]

    def get_iwax(self, niw='all'):

        if niw == 'all':
            niw = self.niw
        elif niw > self.niw or type(niw) != int:
            raise MissmatchError("niw hast to be of type integer < self.niw",
                                 self.filename)

        return self.file['/.axes/iw'][self.niw-niw:self.niw+niw]

    def get_iw(self, quant, atom, iterations='all', niw='all',
               get_ax=False, stat_iter=False):

        if stat_iter is False and self.lastdmftiter == 0:
            stat_iter = True

        if stat_iter:
            key = 'stat'
        else:
            key = 'dmft'

        if iterations == 'all' and stat_iter:
            iterations = np.arange(1, self.laststatiter+1)
        elif iterations == 'all':
            iterations = np.arange(1, self.lastdmftiter+1)
        elif type(iterations) != list:
            raise MissmatchError("iterations have to be of type list",
                                 self.filename)

        if niw == 'all':
            niw = self.niw
        elif niw > self.niw or type(niw) != int:
            raise MissmatchError("niw hast to be of type integer < self.niw",
                                 self.filename)

        atstr      = str('%03i' % self.equiv[atom-1])
        innershape = self.file['/{}/ineq-{}/{}/value'
                               .format(self.ref, atstr, quant)].shape

        if innershape[-3:] == (2*self.n4iwf, 2*self.n4iwf, self.n4iwb):
            val = np.zeros((len(iterations),) + innershape[:-3]
                           + (2*niw, 2*niw, self.n4iwb),
                           dtype=complex)
            for i in iterations:
                val[i-1, ...] = self.file['/{}-{:03d}/ineq-{}/{}/value'
                                          .format(key, i, atstr, quant)
                                          ][..., self.n4iwf-niw:self.n4iwf+niw,
                                                 self.n4iwf-niw:self.n4iwf+niw,
                                                 :]
        else:
            val = np.zeros((len(iterations),) + innershape[:-1] + (2*niw,),
                           dtype=complex)
            for i in iterations:
                val[i-1, ...] = self.file['/{}-{:03d}/ineq-{}/{}/value'
                                          .format(key, i, atstr, quant)
                                          ][..., self.niw-niw:self.niw+niw]

        if get_ax:
            return self.get_iwax(niw=niw), val
        else:
            return val

    def get_tau(self, quant, atom, iterations='all', ntau='all',
                get_ax=False, stat_iter=False):

        if stat_iter is False and self.lastdmftiter == 0:
            stat_iter = True

        if stat_iter:
            key = 'stat'
        else:
            key = 'dmft'

        if iterations == 'all' and stat_iter:
            iterations = np.arange(1, self.laststatiter+1)
        elif iterations == 'all':
            iterations = np.arange(1, self.lastdmftiter+1)
        elif type(iterations) != list:
            raise MissmatchError("iterations have to be of type list",
                                 self.filename)

        if ntau == 'all':
            ntau = self.ntau
        elif ntau > self.ntau or type(ntau) != int:
            raise MissmatchError("ntau hast to be of type integer < self.ntau",
                                 self.filename)

        atstr      = str('%03i' % self.equiv[atom-1])
        innershape = self.file['/{}/ineq-{}/{}/value'
                               .format(self.ref, atstr, quant)].shape

        if innershape[-2:] == (self.ntau, self.ntau):
            val = np.zeros((len(iterations),) + innershape[:-2] + (ntau, ntau),
                           dtype=complex)

            for i in iterations:
                val[i-1, ...] = self.file['/{}-{:03d}/ineq-{}/{}/value'
                                          .format(key, i, atstr, quant)
                                          ][..., 0:ntau, 0:ntau]
        else:
            val = np.zeros((len(iterations),) + innershape[:-1] + (ntau,),
                           dtype=complex)

            for i in iterations:
                val[i-1, ...] = self.file['/{}-{:03d}/ineq-{}/{}/value'
                                          .format(key, i, atstr, quant)
                                          ][..., 0:ntau]

        if get_ax:
            return self.get_tauax(ntau=ntau), val
        else:
            return val

    def get_noax(self, quant, atom, iterations='all', stat_iter=False):

        if stat_iter is False and self.lastdmftiter == 0:
            stat_iter = True

        if stat_iter:
            key = 'stat'
        else:
            key = 'dmft'

        if iterations == 'all' and stat_iter:
            iterations = np.arange(1, self.laststatiter+1)
        elif iterations == 'all':
            iterations = np.arange(1, self.lastdmftiter+1)
        elif type(iterations) != list:
            raise MissmatchError("iterations have to be of type list",
                                 self.filename)

        atstr      = str('%03i' % self.equiv[atom-1])
        innershape = self.file['/{}/ineq-{}/{}/value'
                               .format(self.ref, atstr, quant)].shape
        val = np.zeros((len(iterations),) + innershape, dtype=complex)
        for i in iterations:
            val[i-1, ...] = self.file['/{}-{:03d}/ineq-{}/{}/value'
                                      .format(key, i, atstr, quant)
                                      ][...]

        return val

    def get_giw(self, atom, iterations='all', niw='all', get_ax=False,
                 stat_iter=False):

        return self.get_iw(quant='giw', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_giwfull(self, atom, iterations='all', niw='all', get_ax=False,
                    stat_iter=False):

        return self.get_iw(quant='giw-full', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_g0iw(self, atom, iterations='all', niw='all', get_ax=False,
                 stat_iter=False):

        return self.get_iw(quant='g0iw', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_g0iwfull(self, atom, iterations='all', niw='all', get_ax=False,
                     stat_iter=False):

        return self.get_iw(quant='g0iwfull', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_fiw(self, atom, iterations='all', niw='all', get_ax=False,
                 stat_iter=False):

        return self.get_iw(quant='fiw', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_fiwfull(self, atom, iterations='all', niw='all', get_ax=False,
                    stat_iter=False):

        return self.get_iw(quant='fiw-full', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_siw(self, atom, iterations='all', niw='all', get_ax=False,
                 stat_iter=False):

        return self.get_iw(quant='siw', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_siwfull(self, atom, iterations='all', niw='all', get_ax=False,
                    stat_iter=False):

        return self.get_iw(quant='siw-full', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_g4iw(self, atom, iterations='all', niw='all', get_ax=False,
                 stat_iter=False):

        return self.get_iw(quant='g4iw', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_g4iwpp(self, atom, iterations='all', niw='all', get_ax=False,
                   stat_iter=False):

        return self.get_iw(quant='g4iw-pp', atom=atom, iterations=iterations,
                           niw=niw, get_ax=get_ax, stat_iter=stat_iter)

    def get_gtau(self, atom, iterations='all', ntau='all', get_ax=False,
                 stat_iter=False):

        return self.get_tau(quant='gtau', atom=atom, iterations=iterations,
                            ntau=ntau, get_ax=get_ax, stat_iter=stat_iter)

    def get_gtaufull(self, atom, iterations='all', ntau='all', get_ax=False,
                     stat_iter=False):

        return self.get_tau(quant='gtau-full', atom=atom, iterations=iterations,
                            ntau=ntau, get_ax=get_ax, stat_iter=stat_iter)

    def get_ftau(self, atom, iterations='all', ntau='all', get_ax=False,
                 stat_iter=False):

        return self.get_tau(quant='ftau', atom=atom, iterations=iterations,
                            ntau=ntau, get_ax=get_ax, stat_iter=stat_iter)

    def get_ftaufull(self, atom, iterations='all', ntau='all', get_ax=False,
                     stat_iter=False):

        return self.get_tau(quant='ftau-full', atom=atom, iterations=iterations,
                            ntau=ntau, get_ax=get_ax, stat_iter=stat_iter)

    def get_occ(self, atom, iterations='all', stat_iter=False):

        return self.get_noax(quant='occ', atom=atom, iterations=iterations,
                             stat_iter=stat_iter)

    def get_rho1(self, atom, iterations='all', stat_iter=False):

        return self.get_noax(quant='rho1', atom=atom, iterations=iterations,
                             stat_iter=stat_iter)

    def get_rho2(self, atom, iterations='all', stat_iter=False):

        return self.get_noax(quant='rho2', atom=atom, iterations=iterations,
                             stat_iter=stat_iter)

    def get_dc(self, atom, iterations='all', stat_iter=False):

        return self.get_noax(quant='dc', atom=atom, iterations=iterations,
                             stat_iter=stat_iter)

    def get_muimp(self, atom, iterations='all', stat_iter=False):

        return self.get_noax(quant='muimp', atom=atom,
                             iterations=iterations, stat_iter=stat_iter)

    def get_muimpfull(self, atom, iterations='all', stat_iter=False):

        return self.get_noax(quant='muimp-full', atom=atom,
                             iterations=iterations, stat_iter=stat_iter)

    def get_all_atoms(self, function, **args):

        result = []
        ax = None

        if 'get_ax' in args.keys():
            if args['get_ax'] is True:
                ax, first = function(atom=1, **args)
                args['get_ax'] = False
            else:
                first = function(atom=1, **args)
        else:
            first = function(atom=1, **args)
        result.append(first)

        for i in range(2, self.natoms+1):
            result.append(function(atom=i, **args))

        if ax is not None:
            return ax, result
        else:
            return result

    def get_beta(self):
        return self.file['/.config'].attrs.get('general.beta')

    def check_for_statiter(self):
        return self.laststatiter > 0


def atomloop(function, atoms, **args):

    result = []
    ax = None

    if 'get_ax' in args.keys():
        if args['get_ax'] is True:
            ax, first = function(atom=atoms[0], **args)
            args['get_ax'] = False
        else:
            first = function(atom=atoms[0], **args)
    else:
        first = function(atom=atoms[0], **args)
    result.append(first)

    for i in atoms[1:]:
        result.append(function(atom=i, **args))

    if ax is not None:
        return ax, result
    else:
        return result

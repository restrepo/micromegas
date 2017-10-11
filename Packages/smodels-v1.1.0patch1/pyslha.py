#! /usr/bin/env python

"""\
A simple but flexible handler of the SUSY Les Houches Accord (SLHA) data format.

pyslha is a parser/writer module for particle physics SUSY Les Houches Accord
(SLHA) supersymmetric spectrum/decay files, and a collection of scripts which
use the interface, e.g. for conversion to and from the legacy ISAWIG format, or
to plot the mass spectrum and decay chains.

The current release supports SLHA version 1, and as far as I'm aware is also
fully compatible with SLHA2: the block structures are read and accessed
generically. If you have any problems, please provide an example input file and
I'll happily investigate. SLHA3 is not yet supported (or standardised) but in
recent releases the new structures will not crash the parser. Support will be
added once the format is standardised (and in response to demand!)

The plotting script provides output in PDF, EPS and PNG via LaTeX and the TikZ
graphics package, and as LaTeX/TikZ source for direct embedding into documents or
user-tweaking of the generated output.

Users of version 2.x should note that the interface has changed a little in
version 3.x: there are now generic read()/write() functions which
can operate on filenames or file objects, and all I/O functions now
return/accept a single Doc object rather than a tuple of blocks, decays,
etc. dicts. This single-object interface to the SLHA document allows for more
coherent handling of the data, as well as better robustness against future
changes in the format and planned support for documentation comments.

API example
-----------

>>> import pyslha
>>> # generic read from file:
>>> d = pyslha.read('spcfiles/sps1a.spc')
>>> d
<PySLHA Doc: 22 blocks, 35 decays, 0 xsections>
>>> # or, to ignore blocks known to be e.g. badly formatted:
>>> d = pyslha.read('spcfiles/sps1a.spc', ignoreblocks=['DCINFO'])
>>> d
<PySLHA Doc: 21 blocks, 35 decays, 0 xsections>
>>> d.blocks
Blocks
  SPINFO { 1 : SOFTSUSY; 2 : 2.0.5 }
  MODSEL { 1,1 : sugra }
  SMINPUTS { 1 :  1.27934000e+02; 2 :  1.16637000e-05; ...
  ...
>>> d.blocks.has_key('MODSEL')
True
>>> d.blocks['MODSEL']
MODSEL { 1,1 : sugra }
>>> d.blocks['MODSEL'][1,1]
'sugra'
>>> d.blocks['MODSEL'][1,2] = 'foo'
>>> d.blocks['MODSEL']
MODSEL { 1,1 : sugra; 1,2 : foo }


Blocks
------

The Block interface also supplies dict-like has_key(), keys(), and items()
methods, as well as more specialist value(), set_value() and is_single_valued()
methods for improved access to ALPHA and any other unindexed blocks.


Decays
------

The decay system has a similar interface to that of the generic Block for
accessing particle decay channel information:

>>> d.decays.keys()
[23, 24, 6, 25, 35, 36, 37, 1000021, 1000006, 2000006, 1000005, 2000005, ...
>>> d.decays[1000021]
1000021 : mass = 6.07713704e+02 GeV : total width = 5.50675438e+00 GeV
   1.05840237e-01 [1000005, -5]
   1.05840237e-01 [-1000005, 5]
   ...

The DECAY blocks in SLHA files are in fact mapped into Particle objects, each of
which contains multiple Decay objects. Check the Python documentation for Doc,
Block, Particle and Decay, plus the pyslha module free functions for I/O.


Cross-sections
--------------

Cross-section information, cf. the new XSECTION block type, is available via
Doc.xsections. As for decays, the API structure is a bit different from the text
format: xsections is a dict of Process objects, each of which contains all
XSECTION lines for a given list of (sorted) initial + final state particle
IDs. This is different from the text format in that a single Process contains
cross-sections, as XSec objects, for multiple centre-of-mass energies
(sqrts). For example:

>>> print d.xsections
[(2212, 2212, 1000001, 1000003), (2212, 2212, -1000002, 2000002), ...]
>>> myproc = d.xsections[2212,2212,1000001,1000003]

The Process interface supplies a convenient method for filtering the
contained XSecs on any of their defining attributes, including scale scheme, QCD
and EW orders, multiplicative scale factors, PDF ID code, and
generator/integrator.

>>> myproc.get_xsecs(sqrts=13000., kappa_r=2., code='Prospino')


Citation
--------

If you use PySLHA, for either model data handling or spectrum visualisation,
please cite the paper: http://arxiv.org/abs/1305.4194

TODOs:

  For 3.2.x:
   * In set_value, if first item is non-int, treat as None-indexed.
   * Refine value string heuristic for strings with ints in them?

  For 3.3.0:
   * Use Doc to handle document-level header comments.
   * Use _dict to handle block and decay summary comments.
   * Preserve _inline_ comments from read -> write (needs full-line/inline
     comment separation). Can use separate comment dicts in Block, Decay,
     etc. and attach a multiline .comment attr to the returned/written dicts.

  Later, maybe:
   * Identify HERWIG decay matrix element to use in ISAWIG.
   * Handle RPV SUSY in ISAWIG.
"""

__author__ = "Andy Buckley <andy.buckley@cern.ch>"
__version__ = "3.2.0"


## Python version: require >= 2.6 (including Py3) for "as" exception handling syntax and ternary "x if a else y" syntax
import sys
if sys.version_info < (2,6):
    raise Exception("pyslha requires Python >= 2.6")


###############################################################################
## Private utility functions

def _mkdict():
    """Return an OrderedDict if possible, or a normal dict if not."""
    try:
        from collections import OrderedDict
        return OrderedDict()
    except:
        try:
            from ordereddict import OrderedDict
            return OrderedDict()
        except:
            return dict()

_d = _mkdict()
class _dict(type(_d)):
    """
    A cosmetic wrapper on an OrderedDict if possible, or a normal dict if not.

    TODO: Add (optional) comment support
    """

    def __init__(self, name=None):
        super(_dict, self).__init__()
        self.name = name

    def __getitem__(self, key):
        return super(_dict, self).__getitem__(_autotuple(key))

    def __setitem__(self, key, val):
        return super(_dict, self).__setitem__(_autotuple(key), val)

    def __repr__(self):
        s = ""
        if self.name:
            s += self.name + "\n  "
        if self.values():
            s += "\n  ".join(str(val) for val in self.values())
        else:
            s += "No entries"
        return s

def _autotype(var):
    """Automatically convert strings to numerical types if possible."""
    if type(var) is not str:
        return var
    if var.isdigit() or (var.startswith("-") and var[1:].isdigit()):
        return int(var)
    try:
        f = float(var)
        return f
    except ValueError:
        return var

def _autostr(var, precision=8):
    """Automatically format numerical types as the right sort of string."""
    # print("@", type(var), "@")
    if type(var) is str:
        return var
    elif type(var) is float:
        return ("% ." + str(precision) + "g") % var
    elif not hasattr(var, "__iter__"):
        return str(var)
    else:
        return ",".join(_autostr(subval) for subval in var)

def _autotuple(a):
    """Automatically convert the supplied iterable to a scalar or tuple as appropriate."""
    if not hasattr(a, "__iter__"):
        return a
    if len(a) == 1:
        return a[0]
    return tuple(a)

def _read(f):
    """Read a file's contents, autodetecting whether the arg is a file or filename,
    and treating '-' as as indication to read from stdin."""
    if type(f) is str:
        if f == "-":
            return sys.stdin.read()
        else:
            with open(f, "r") as ff:
                return ff.read()
    else:
        return f.read()

def _write(f, txt):
    """Write to a file, autodetecting whether the arg is a file or filename,
    and treating '-' as as indication to write to stdout."""
    if type(f) is str:
        if f == "-":
            sys.stdout.write(txt)
        else:
            with open(f, "w") as ff:
                ff.write(txt)
    else:
        f.write(txt)


###############################################################################
## Exceptions

class AccessError(Exception):
    "Exception object to be raised when a SLHA block is accessed in an invalid way"
    def __init__(self, errmsg):
        self.msg = errmsg
    def __str__(self):
        return self.msg

class ParseError(Exception):
    "Exception object to be raised when a spectrum file/string is malformed"
    def __init__(self, errmsg):
        self.msg = errmsg
    def __str__(self):
        return self.msg


###############################################################################
## The Doc top-level container object


class Doc(object):
    "Top level container for everything in an SLHA record"
    def __init__(self, blocks, decays=None, xsections=None):
        self.blocks = blocks
        self.decays = decays
        self.xsections = xsections
        self.comment = ""

    def write(self, filename=None, ignorenobr=True, precision=8):
        """
        Convenient method for converting an SLHA Doc object to SLHA format,
        either returned as a string or written to a file depending on whether
        the filename variable is None.
        """
        if filename is None:
            return writeSLHA(self, ignorenobr, precision)
        else:
            write(filename, self, ignorenobr, precision)

    def __repr__(self):
        s = "<PySLHA Doc: %d blocks, %d decays, %d xsections" % \
            (len(self.blocks or []), len(self.decays or []), len(self.xsections or []))
        if self.comment:
            s += " %s" % self.comment
        s += ">"
        return s


###############################################################################
## The data block, decay and particle classes

class Block(object):
    """
    Object representation of any BLOCK elements read from an SLHA file.

    Blocks have a name, may have an associated Q value, and contain a collection
    of data entries, each indexed by one or more keys. Entries in the dictionary
    are stored as numeric types (int or float) when a cast from the string in
    the file has been possible.

    Block is closely related to a Python dict (and, in fact, is implemented via
    an OrderedDict when possible). The preferred methods of entry access use the
    dict-like [] operator for getting and setting, and the keys() and items()
    methods for iteration. Purely iterating over the object behaves like keys(),
    as for an ordinary dict.

    Multiple (integer) indices are possible, especially for entries in mixing
    matrix blocks. These are now implemented in the natural way, e.g. for access
    to the (1,2) element of a mixing matrix block, use 'bmix[1,2] = 0.123' and
    'print bmix[1,2]'. The value() and set_value() functions behave
    similarly. Multi-element values are also permitted.

    It is possible, although not usual, to store unindexed values in a
    block. This is only supported when that entry is the only one in the block,
    and it is stored in the normal fashion but with None as the lookup key. The
    value() method may be used without a key argument to retrieve this value, if
    the is_single_valued() method returns True, and similarly the set_value()
    method may be used to set it if only one argument is supplied and the object
    is compatible.
    """
    def __init__(self, name, q=None, entries=None):
        self.name = name
        self.entries = _dict()
        if entries is not None:
            self.entries.update(entries)
        self.q = _autotype(q)

    # TODO: Rename? To what?
    def add_entry(self, args):
        """Add an entry to the block from an iterable (i.e. list or tuple) of
        strings, or from a whitespace-separated string.

        This method is just for convenience: it splits the single string
        argument if necessary and converts the list of strings into numeric
        types when possible. For the treatment of the resulting iterable see the
        set_value method.
        """
        ## If the argument is a single string, split it and proceed
        if type(args) is str:
            args = args.split()
        ## Check that the arg is an iterable
        if not hasattr(args, "__iter__"):
            raise AccessError("Block entries must be iterable")
        ## Auto-convert the types in the list
        args = [_autotype(a) for a in args]
        ## Re-join consecutive strings into single entries
        i = 0
        while i < len(args)-1:
            if type(args[i]) is str and type(args[i+1]) is str:
                args[i] += " " + args[i+1]
                del args[i+1]
                continue
            i += 1
        ## Add the entry to the map, with appropriate indices
        self.set_value(*args)

    def is_single_valued(self):
        """Return true if there is only one entry, and it has no index: the
        'value()' attribute may be used in that case without an argument."""
        return len(self.entries) == 1 and list(self.entries.keys())[0] is None

    def value(self, key=None, default=1):
        """Get a value from the block with the supplied key.

        If no key is given, then the block must contain only one non-indexed
        value otherwise an AccessError exception will be raised.\
        """
        if key is None and not self.is_single_valued():
            raise AccessError("Tried to access unique value of multi-value block " + self.name)
        if not self.has_key(key):
            return default
        return self.entries[key]

    "Alias"
    get = value

    def set_value(self, *args):
        """Set a value in the block via supplied key/value arguments.

        Indexing is determined automatically: any leading integers will be
        treated as a multi-dimensional index, with the remaining entries being a
        (potentially multi-dimensional) value. If all N args are ints, then the
        first N-1 are treated as the index and the Nth as the value.

        If there is only one arg it will be treated as the value of a
        single-valued block. In this case the block must already contain at most
        one non-indexed value otherwise an AccessError exception will be
        raised.\
        """
        if len(args) == 0:
            raise AccessError("set_value() called without arguments")
        elif len(args) == 1:
            if len(self.entries) > 0 and not self.is_single_valued():
                raise AccessError("Tried to set a unique value on a multi-value block " + self.name)
            self.entries[None] = args[0]
        else:
            ## Find the first non-integer -- all previous items are indices
            i_first_nonint = -1
            for i, x in enumerate(args):
                if type(x) is not int:
                    i_first_nonint = i
                    break
            # if i_first_nonint == 0:
            #     raise AccessError("Attempted to set a block entry with a non-integer(s) index: %s" % str(args))
            # else:
            self.entries[_autotuple(args[:i_first_nonint])] = _autotuple(args[i_first_nonint:])

    "Alias"
    set = set_value

    def has_key(self, key):
        """Does the block have the given key?"""
        return key in self.entries

    def keys(self):
        """Access the block item keys."""
        return self.entries.keys()

    def values(self):
        """Access the block item values."""
        return self.entries.values()

    def items(self, key=None):
        """Access the block items as (key,value) tuples."""
        return self.entries.items()

    def __len__(self):
        return len(self.entries)

    def __contains__(self, key):
        return key in self.entries

    def __iter(self):
        return self.entries.__iter__()

    def __getitem__(self, key):
        return self.entries[key]

    def __setitem__(self, key, value):
        if key is not None and type(key) is not int and not all(type(x) is int for x in key):
            raise AccessError("Attempted to set a block entry with a non-integer(s) index")
        self.entries[key] = value

    def __lt__(self, other):
        return (self.name < other.name) and (self.entries < other.entries)

    def __repr__(self):
        s = self.name
        if self.q is not None:
            s += " (Q=%s)" % self.q
            s += " { " + "; ".join(_autostr(k,3) + " : " + _autostr(v,3) for k, v in self.items()) + " }"
        return s


class Decay(object):
    """
    Object representing a decay entry on a particle decribed by the SLHA file.
    'Decay' objects are not a direct representation of a DECAY block in an SLHA
    file... that role, somewhat confusingly, is taken by the Particle class.

    Decay objects have three properties: a branching ratio, br, an nda number
    (number of daughters == len(ids)), and a tuple of PDG PIDs to which the
    decay occurs. The PDG ID of the particle whose decay this represents may
    also be stored, but this is normally known via the Particle in which the
    decay is stored.
    """
    def __init__(self, br, nda, ids, parentid=None):
        self.parentid = parentid
        self.br = br
        self.nda = nda
        self.ids = list(ids)
        assert(self.nda == len(self.ids))

    def __lt__(self, other):
        return (other.br < self.br)

    def __repr__(self):
        return "% .2g %s" % (self.br, self.ids)


class Particle(object):
    """
    Representation of a single, specific particle, decay block from an SLHA
    file.  These objects are not themselves called 'Decay', since that concept
    applies more naturally to the various decays found inside this
    object. Particle classes store the PDG ID (pid) of the particle being
    represented, and optionally the mass (mass) and total decay width
    (totalwidth) of that particle in the SLHA scenario. Masses may also be found
    via the MASS block, from which the Particle.mass property is filled, if at
    all. They also store a list of Decay objects (decays) which are probably the
    item of most interest.
    """
    def __init__(self, pid, totalwidth=None, mass=None):
        self.pid = pid
        self.totalwidth = totalwidth
        self.mass = mass
        self.decays = []

    def add_decay(self, br, nda, ids):
        self.decays.append(Decay(br, nda, ids))
        self.decays.sort()

    def __lt__(self, other):
        if abs(self.pid) == abs(other.pid):
            return (self.pid < other.pid)
        return (abs(self.pid) < abs(other.pid))

    def __repr__(self):
        s = str(self.pid)
        if self.mass is not None:
            s += " : mass = %.3g GeV" % self.mass
        if self.totalwidth is not None:
            s += " : total width = %.3g GeV" % self.totalwidth
        for d in self.decays:
            if d.br > 0.0:
                s += "\n  %s" % d
        return s


class XSec(object):
    """
    A cross-section value for a specific combination of energy, scheme, scale
    choice, PDF, and the computational code (+ version) that made the calculation.
    """
    def __init__(self, sqrts, scale_scheme, qcd_order, ew_order, kappa_f, kappa_r, pdf_id, value, code=None):
        self.sqrts = sqrts
        self.scale_scheme = scale_scheme
        self.qcd_order = qcd_order
        self.ew_order = ew_order
        self.kappa_f = kappa_f
        self.kappa_r = kappa_r
        self.pdf_id = pdf_id
        self.value = value
        ## Make sure that code is always a 2-list
        self.code = code
        if self.code is None:
            self.code = [None, None]
        elif type(self.code) is str:
            self.code = [self.code, None]
        assert len(self.code) == 2

    @property
    def scale_scheme_str(self):
        return ("avg mass", "fixed", "s_hat", "mT")[self.scale_scheme]

    @property
    def qcd_order_str(self):
        return ("Born", "NLO", "NLO+LL")[self.qcd_order] if self.qcd_order is not None else ""

    @property
    def ew_order_str(self):
        return ("Born", "NLO", "NLO+LL")[self.ew_order] if self.ew_order is not None else ""

    def __cmp__(self, other):
        return \
            (self.sqrts < other.sqrts) and \
            (self.scale_scheme < other.scale_scheme) and \
            (self.qcd_order < other.qcd_order) and \
            (self.ew_order < other.ew_order) and \
            (self.kappa_f < other.kappa_f) and \
            (self.kappa_r < other.kappa_r) and \
            (self.pdf_id < other.pdf_id) and \
            (self.value < other.value) and \
            (self.code < other.code)

    def __repr__(self):
        s = "sqrt(s) = %g GeV" % self.sqrts
        if self.scale_scheme_str:
            s += ", " + self.scale_scheme_str + " scale scheme"
        if self.qcd_order_str:
            s += ", " + self.qcd_order_str + " QCD"
        if self.ew_order_str:
            s += ", " + self.ew_order_str + " EW"
        s += "; "
        s += "K_F = %g, " % self.kappa_f
        s += "K_R = %g, " % self.kappa_r
        s += "PDF = %d " % self.pdf_id
        s += ": xsec = %g pb" % self.value
        if self.code[0]:
            vstr = ""
            if self.code[1]:
                vstr = " " + self.code[1]
            s += " (%s%s)" % (self.code[0], vstr)
        return s


class Process(object):
    """
    Representation of a single a b -> x y z ... scattering process, containing
    cross-sections for that topology with different energies, schemes, scale
    choices, PDFs, and origins. A method is provided to filter the contained
    cross-section values by these criteria.
    """
    def __init__(self, pidsinitial, pidsfinal):
        self.pidsinitial = sorted(pidsinitial)
        self.pidsfinal = sorted(pidsfinal)
        self.xsecs = []

    def add_xsec(self, sqrts, scale_scheme, qcd_order, ew_order, kappa_f, kappa_r, pdf_id, value, code=None):
        "Add an XSec object by passing the XSec constructor arguments directly"
        self.xsecs.append(XSec(sqrts, scale_scheme, qcd_order, ew_order, kappa_f, kappa_r, pdf_id, value, code))
        #self.xsecs.sort()

    def get_xsecs(self, sqrts=None, scale_scheme=None, qcd_order=None, ew_order=None, kappa_f=None, kappa_r=None, pdf_id=None, code=None):
        "Get all contained XSec objects matching the specified (all optional) filter values"
        rtn = self.xsecs
        if sqrts:
            rtn = [x for x in rtn if x.sqrts == sqrts]
        if scale_scheme:
            rtn = [x for x in rtn if x.scale_scheme == scale_scheme]
        if qcd_order:
            rtn = [x for x in rtn if x.qcd_order == qcd_order]
        if ew_order:
            rtn = [x for x in rtn if x.ew_order == ew_order]
        if kappa_f:
            rtn = [x for x in rtn if x.kappa_f == kappa_f]
        if kappa_r:
            rtn = [x for x in rtn if x.kappa_r == kappa_r]
        if pdf_id:
            rtn = [x for x in rtn if x.pdf_id == pdf_id]
        if code:
            if type(code) is str or code[1] is None:
                rtn = [x for x in rtn if x.code[0] == code[0]]
            else:
                rtn = [x for x in rtn if x.code == code]
        return rtn

    @property
    def sqrtses(self):
        # TODO: could just use an inline sorted set literal, but that would require Python >= 2.7... can't do that just yet
        rtn = []
        for x in self.xsecs:
            if x.sqrts in rtn:
                continue
            rtn.append(x.sqrts)
        rtn.sort()
        return rtn

    def __lt__(self, other):
        if self.pidsinitial != other.pidsinitial:
            return (self.pidsinitial < other.pidsinitial)
        if self.pidsfinal != other.pidsfinal:
            return (self.pidsfinal < other.pidsfinal)
        return (self.xsecs < other.xsecs)

    def __repr__(self):
        s = " ".join(str(pid) for pid in self.pidsinitial)
        s += "  ->  "
        s += " ".join(str(pid) for pid in self.pidsfinal)
        # s += " %d " % len(self.xsecs)
        for x in self.xsecs:
            s += "\n    %s" % x
        return s




###############################################################################
## SLHA parsing and writing functions

def readSLHA(spcstr, ignorenobr=False, ignorenomass=False, ignoreblocks=[]):
    """
    Read an SLHA definition from a string, returning dictionaries of blocks,
    decays ('Particle's), and cross-sections ('Process'es).

    If the ignorenobr parameter is True, do not store decay entries with a
    branching ratio of zero.

    If the ignorenomass parameter is True, parse file even if mass block is
    absent in the file (default is to raise a ParseError).
    """
    blocks = _dict("Blocks")
    decays = _dict("Decays")
    xsections = _dict("Cross-sections")
    #
    currentblock = None
    currentdecay = None
    currentsqrts = None
    currentproc = None
    #
    import re
    re_indented = re.compile(r"^\s+")
    for line in spcstr.splitlines():
        ## Handle (ignore) comment lines
        # TODO: Store block/entry comments
        if line.startswith("#"):
            continue
        if "#" in line:
            line = line[:line.index("#")]
        ## Ignore empty lines (after comment removal and whitespace trimming)
        if not line.strip():
            continue

        ## Section header lines start with a non-whitespace character, data lines have a whitespace indent
        ## NOTE: we'll accept any whitespace rather than just spaces cf. the SLHA standard.
        #if not line.startswith(" "):
        if not re_indented.match(line):
            line = line.strip()

            ## Handle BLOCK start lines
            if line.upper().startswith("BLOCK"):
                #print line
                match = re.match(r"BLOCK\s+(\w+)(\s+Q\s*=\s*.+)?", line.upper())
                if not match:
                    continue
                blockname = match.group(1)
                if blockname in ignoreblocks:
                    currentblock = None
                    currentdecay = None
                else:
                    qstr = match.group(2)
                    if qstr is not None:
                        qstr = qstr[qstr.find("=")+1:].strip()
                    currentblock = blockname
                    currentdecay = None
                    currentsqrts = None
                    currentproc = None
                    blocks[blockname] = Block(blockname, q=qstr)
                    # TODO: Warn or combine if there are multiple blocks with the same name etc.?
                    #blocks.setdefault(blockname, Block(blockname, q=qstr))

            ## Handle DECAY start lines
            elif line.upper().startswith("DECAY"):
                match = re.match(r"DECAY\s+(-?\d+)\s+([\d\.E+-]+|NAN).*", line.upper())
                if not match:
                    continue
                pdgid = int(match.group(1))
                width = float(match.group(2)) if match.group(2) != "NAN" else None
                currentblock = "DECAY"
                currentdecay = pdgid
                currentsqrts = None
                currentproc = None
                decays[pdgid] = Particle(pdgid, width)
                # TODO: Warn or combine if there are multiple decay blocks with the same PID?
                #decays.setdefault(pdgid, Particle(pdgid, width))

            ## Handle XSECTION start lines
            elif line.upper().startswith("XSECTION"):
                match = re.match(r"XSECTION\s+([\d+\.E+-]+)\s+(-?\d+)\s+(-?\d+)\s+(\d+)\s+(.*)", line.upper())
                if not match:
                    continue
                sqrts = float(match.group(1))
                pdgid1 = int(match.group(2))
                pdgid2 = int(match.group(3))
                nf = int(match.group(4))
                fspdgids = []
                if match.group(4):
                    fspdgids = [int(i) for i in match.group(5).strip().split()]
                fspdgids.sort()
                assert len(fspdgids) == nf
                currentblock = "XSECTION"
                currentdecay = None
                currentsqrts = sqrts
                currentproc = [pdgid1,pdgid2] + fspdgids
                ## There might be more than one energy so we need to combine:
                # TODO: need to do something about the tuplification for UX reasons?
                xsections.setdefault(_autotuple(currentproc), Process([pdgid1,pdgid2], fspdgids))

            ## Handle unknown section type start lines (and continue ignoring until a non-header line is found)
            elif type(_autotype(line.split()[0])) is str:
                sys.stderr.write("Ignoring unknown section type: %s\n" % line.split()[0])
                currentblock = None
                currentdecay = None
                currentsqrts = None
                currentproc = None

        ## This non-empty line starts with an indent, hence must be an in-block data line (provided we _are_ in a block, otherwise ignore)
        elif currentblock is not None:
            # print "@", currentblock
            items = line.split()
            if not items: #< Shouldn't be possible, but whatever...
                continue
            # print items
            if currentblock == "DECAY":
                br = float(items[0]) if items[0].upper() != "NAN" else None
                nda = int(items[1])
                ids = map(int, items[2:])
                if br or not ignorenobr:
                    if not br or br <= 0:
                        br = 0.0
                    decays[currentdecay].add_decay(br, nda, ids)
            elif currentblock == "XSECTION":
                scale_scheme = int(items[0])
                qcd_order = int(items[1])
                ew_order = int(items[2])
                kappa_f = float(items[3])
                kappa_r = float(items[4])
                pdf_id = int(items[5])
                value = float(items[6])
                code = items[7] if len(items) >= 8 else None
                version = items[8] if len(items) >= 9 else None
                xsections[currentproc].add_xsec(currentsqrts, scale_scheme, qcd_order, ew_order, kappa_f, kappa_r, pdf_id, value, [code, version])
            else: # we're in a generic BLOCK
                blocks[currentblock].add_entry(items)

    ## Try to populate Particle masses from the MASS block
    # print blocks.keys()
    try:
        for pid in blocks["MASS"].keys():
            ## Set zero width for particles that have mass but no decays, assuming them to be stable
            decays.setdefault(pid, Particle(pid, 0)).mass = blocks["MASS"][pid]
    except:
        if not ignorenomass:
            raise ParseError("No MASS block found: cannot populate particle masses")

    rtn = Doc(blocks, decays, xsections)
    return rtn


def writeSLHABlocks(blocks, precision=8):
    """Return an SLHA definition as a string, from the supplied blocks dict."""
    if not blocks:
        return ""
    sep = 3 * " "
    blockstrs = []
    for bname, b in blocks.items():
        namestr = b.name
        if b.q is not None:
            namestr += " Q= " + _autostr(float(b.q), precision)
        blockstr = "BLOCK %s\n" % namestr
        entrystrs = []
        for k, v in b.items():
            entrystr = sep
            if type(k) == tuple:
                entrystr += sep.join(_autostr(i, precision) for i in k)
            elif k is not None:
                entrystr += _autostr(k, precision)
            entrystr += sep + _autostr(v, precision)
            entrystrs.append(entrystr)
        blockstr += "\n".join(entrystrs)
        blockstrs.append(blockstr)
    return "\n\n".join(blockstrs)



def writeSLHADecays(decays, ignorenobr=False, precision=8):
    """Return an SLHA decay definition as a string, from the supplied decays dict."""
    if not decays:
        return ""
    sep = 3 * " "
    blockstrs = []
    for pid, particle in decays.items():
        blockstr = ("DECAY %d " % particle.pid) + _autostr(particle.totalwidth or -1, precision) + "\n"
        decaystrs = []
        for d in particle.decays:
            if d.br > 0.0 or not ignorenobr:
                products_str = sep.join("% d" % i for i in d.ids)
                decaystr = sep + _autostr(d.br, precision) + sep + str(len(d.ids)) + sep + products_str
                decaystrs.append(decaystr)
        blockstr += "\n".join(decaystrs)
        blockstrs.append(blockstr)
    return "\n\n".join(blockstrs)



def writeSLHAXSections(xsections, precision=8):
    """Return an SLHA cross-section definition as a string, from the supplied xsections dict."""
    if not xsections:
        return ""
    sep = 3 * " "
    blockstrs = []
    for _, proc in xsections.items():
        for sqrts in proc.sqrtses:
            blockstr = "XSECTION " + _autostr(sqrts) + " "
            blockstr += " ".join(str(pid) for pid in proc.pidsinitial)
            blockstr += " %d " % len(proc.pidsfinal)
            blockstr += " ".join(str(pid) for pid in proc.pidsfinal)
            for x in proc.get_xsecs(sqrts=sqrts):
                blockstr += "\n" + sep
                blockstr += "%d %d %d " % (x.scale_scheme, x.qcd_order, x.ew_order)
                blockstr += "%s %s " % (_autostr(x.kappa_f), _autostr(x.kappa_r))
                blockstr += "%d %s" % (x.pdf_id, _autostr(x.value))
                if x.code[0]:
                    blockstr += " " + x.code[0]
                if x.code[1]:
                    blockstr += " " + x.code[1]
            blockstrs.append(blockstr)
    return "\n\n".join(blockstrs)



def writeSLHA(doc, ignorenobr=False, precision=8):
    """Return an SLHA definition as a string, from the supplied Doc."""
    ss = [x for x in (writeSLHABlocks(doc.blocks, precision), writeSLHADecays(doc.decays, ignorenobr, precision), writeSLHAXSections(doc.xsections, precision)) if x]
    return "\n\n".join(ss)



###############################################################################
## PDG <-> HERWIG particle ID code translations for ISAWIG handling

## Static array of HERWIG IDHW codes mapped to PDG MC ID codes, based on
## http://www.hep.phy.cam.ac.uk/~richardn/HERWIG/ISAWIG/susycodes.html
## + the IDPDG array and section 4.13 of the HERWIG manual.
_HERWIGID2PDGID = {}
_HERWIGID2PDGID[7]   = -1
_HERWIGID2PDGID[8]   = -2
_HERWIGID2PDGID[9]   = -3
_HERWIGID2PDGID[10]  = -4
_HERWIGID2PDGID[11]  = -5
_HERWIGID2PDGID[12]  = -6
_HERWIGID2PDGID[13]  =  21
_HERWIGID2PDGID[59]  =  22
_HERWIGID2PDGID[121] =  11
_HERWIGID2PDGID[122] =  12
_HERWIGID2PDGID[123] =  13
_HERWIGID2PDGID[124] =  14
_HERWIGID2PDGID[125] =  15
_HERWIGID2PDGID[126] =  16
_HERWIGID2PDGID[127] = -11
_HERWIGID2PDGID[128] = -12
_HERWIGID2PDGID[129] = -13
_HERWIGID2PDGID[130] = -14
_HERWIGID2PDGID[131] = -15
_HERWIGID2PDGID[132] = -16
_HERWIGID2PDGID[198] =  24 # W+
_HERWIGID2PDGID[199] = -24 # W-
_HERWIGID2PDGID[200] =  23 # Z0
_HERWIGID2PDGID[201] =  25 ## SM HIGGS
_HERWIGID2PDGID[203] =  25 ## HIGGSL0 (== PDG standard in this direction)
_HERWIGID2PDGID[204] =  35 ## HIGGSH0
_HERWIGID2PDGID[205] =  36 ## HIGGSA0
_HERWIGID2PDGID[206] =  37 ## HIGGS+
_HERWIGID2PDGID[207] = -37 ## HIGGS-
_HERWIGID2PDGID[401] =  1000001 ## SSDLBR
_HERWIGID2PDGID[407] = -1000001 ## SSDLBR
_HERWIGID2PDGID[402] =  1000002 ## SSULBR
_HERWIGID2PDGID[408] = -1000002 ## SSUL
_HERWIGID2PDGID[403] =  1000003 ## SSSLBR
_HERWIGID2PDGID[409] = -1000003 ## SSSL
_HERWIGID2PDGID[404] =  1000004 ## SSCLBR
_HERWIGID2PDGID[410] = -1000004 ## SSCL
_HERWIGID2PDGID[405] =  1000005 ## SSB1BR
_HERWIGID2PDGID[411] = -1000005 ## SSB1
_HERWIGID2PDGID[406] =  1000006 ## SST1BR
_HERWIGID2PDGID[412] = -1000006 ## SST1
_HERWIGID2PDGID[413] =  2000001 ## SSDR
_HERWIGID2PDGID[419] = -2000001 ## SSDRBR
_HERWIGID2PDGID[414] =  2000002 ## SSUR
_HERWIGID2PDGID[420] = -2000002 ## SSURBR
_HERWIGID2PDGID[415] =  2000003 ## SSSR
_HERWIGID2PDGID[421] = -2000003 ## SSSRBR
_HERWIGID2PDGID[416] =  2000004 ## SSCR
_HERWIGID2PDGID[422] = -2000004 ## SSCRBR
_HERWIGID2PDGID[417] =  2000005 ## SSB2
_HERWIGID2PDGID[423] = -2000005 ## SSB2BR
_HERWIGID2PDGID[418] =  2000006 ## SST2
_HERWIGID2PDGID[424] = -2000006 ## SST2BR
_HERWIGID2PDGID[425] =  1000011 ## SSEL-
_HERWIGID2PDGID[431] = -1000011 ## SSEL+
_HERWIGID2PDGID[426] =  1000012 ## SSNUEL
_HERWIGID2PDGID[432] = -1000012 ## SSNUELBR
_HERWIGID2PDGID[427] =  1000013 ## SSMUL-
_HERWIGID2PDGID[433] = -1000013 ## SSMUL+
_HERWIGID2PDGID[428] =  1000014 ## SSNUMUL
_HERWIGID2PDGID[434] = -1000014 ## SSNUMLBR
_HERWIGID2PDGID[429] =  1000015 ## SSTAU1-
_HERWIGID2PDGID[435] = -1000015 ## SSTAU1+
_HERWIGID2PDGID[430] =  1000016 ## SSNUTL
_HERWIGID2PDGID[436] = -1000016 ## SSNUTLBR
_HERWIGID2PDGID[437] =  2000011 ## SSEL-
_HERWIGID2PDGID[443] = -2000011 ## SSEL+
_HERWIGID2PDGID[438] =  2000012 ## SSNUEL
_HERWIGID2PDGID[444] = -2000012 ## SSNUELBR
_HERWIGID2PDGID[439] =  2000013 ## SSMUL-
_HERWIGID2PDGID[445] = -2000013 ## SSMUL+
_HERWIGID2PDGID[440] =  2000014 ## SSNUMUL
_HERWIGID2PDGID[446] = -2000014 ## SSNUMLBR
_HERWIGID2PDGID[441] =  2000015 ## SSTAU1-
_HERWIGID2PDGID[447] = -2000015 ## SSTAU1+
_HERWIGID2PDGID[442] =  2000016 ## SSNUTL
_HERWIGID2PDGID[448] = -2000016 ## SSNUTLBR
_HERWIGID2PDGID[449] =  1000021 ## GLUINO
_HERWIGID2PDGID[450] =  1000022 ## NTLINO1
_HERWIGID2PDGID[451] =  1000023 ## NTLINO2
_HERWIGID2PDGID[452] =  1000025 ## NTLINO3
_HERWIGID2PDGID[453] =  1000035 ## NTLINO4
_HERWIGID2PDGID[454] =  1000024 ## CHGINO1+
_HERWIGID2PDGID[456] = -1000024 ## CHGINO1-
_HERWIGID2PDGID[455] =  1000037 ## CHGINO2+
_HERWIGID2PDGID[457] = -1000037 ## CHGINO2-
_HERWIGID2PDGID[458] =  1000039 ## GRAVTINO

def herwigid2pdgid(hwid):
    """
    Convert a particle ID code in the HERWIG internal IDHW format (as used by
    ISAWIG) into its equivalent in the standard PDG ID code definition.
    """
    return _HERWIGID2PDGID.get(hwid, hwid)


## PDG MC ID codes mapped to HERWIG IDHW codes, based on
## http://www.hep.phy.cam.ac.uk/~richardn/HERWIG/ISAWIG/susycodes.html
## + the IDPDG array and section 4.13 of the HERWIG manual.
_PDGID2HERWIGID = {}
_PDGID2HERWIGID[      -1] = 7
_PDGID2HERWIGID[      -2] = 8
_PDGID2HERWIGID[      -3] = 9
_PDGID2HERWIGID[      -4] = 10
_PDGID2HERWIGID[      -5] = 11
_PDGID2HERWIGID[      -6] = 12
_PDGID2HERWIGID[      21] = 13
_PDGID2HERWIGID[      22] = 59
_PDGID2HERWIGID[      11] = 121
_PDGID2HERWIGID[      12] = 122
_PDGID2HERWIGID[      13] = 123
_PDGID2HERWIGID[      14] = 124
_PDGID2HERWIGID[      15] = 125
_PDGID2HERWIGID[      16] = 126
_PDGID2HERWIGID[     -11] = 127
_PDGID2HERWIGID[     -12] = 128
_PDGID2HERWIGID[     -13] = 129
_PDGID2HERWIGID[     -14] = 130
_PDGID2HERWIGID[     -15] = 131
_PDGID2HERWIGID[     -16] = 132
_PDGID2HERWIGID[      24] = 198 ## W+
_PDGID2HERWIGID[     -24] = 199 ## W-
_PDGID2HERWIGID[      23] = 200 ## Z
_PDGID2HERWIGID[      25] = 203 ## HIGGSL0 (added for PDG standard -> HERWIG IDHW) # TODO: should be 201?
_PDGID2HERWIGID[      26] = 203 ## HIGGSL0
_PDGID2HERWIGID[      35] = 204 ## HIGGSH0
_PDGID2HERWIGID[      36] = 205 ## HIGGSA0
_PDGID2HERWIGID[      37] = 206 ## HIGGS+
_PDGID2HERWIGID[     -37] = 207 ## HIGGS-
_PDGID2HERWIGID[ 1000001] = 401 ## SSDLBR
_PDGID2HERWIGID[-1000001] = 407 ## SSDLBR
_PDGID2HERWIGID[ 1000002] = 402 ## SSULBR
_PDGID2HERWIGID[-1000002] = 408 ## SSUL
_PDGID2HERWIGID[ 1000003] = 403 ## SSSLBR
_PDGID2HERWIGID[-1000003] = 409 ## SSSL
_PDGID2HERWIGID[ 1000004] = 404 ## SSCLBR
_PDGID2HERWIGID[-1000004] = 410 ## SSCL
_PDGID2HERWIGID[ 1000005] = 405 ## SSB1BR
_PDGID2HERWIGID[-1000005] = 411 ## SSB1
_PDGID2HERWIGID[ 1000006] = 406 ## SST1BR
_PDGID2HERWIGID[-1000006] = 412 ## SST1
_PDGID2HERWIGID[ 2000001] = 413 ## SSDR
_PDGID2HERWIGID[-2000001] = 419 ## SSDRBR
_PDGID2HERWIGID[ 2000002] = 414 ## SSUR
_PDGID2HERWIGID[-2000002] = 420 ## SSURBR
_PDGID2HERWIGID[ 2000003] = 415 ## SSSR
_PDGID2HERWIGID[-2000003] = 421 ## SSSRBR
_PDGID2HERWIGID[ 2000004] = 416 ## SSCR
_PDGID2HERWIGID[-2000004] = 422 ## SSCRBR
_PDGID2HERWIGID[ 2000005] = 417 ## SSB2
_PDGID2HERWIGID[-2000005] = 423 ## SSB2BR
_PDGID2HERWIGID[ 2000006] = 418 ## SST2
_PDGID2HERWIGID[-2000006] = 424 ## SST2BR
_PDGID2HERWIGID[ 1000011] = 425 ## SSEL-
_PDGID2HERWIGID[-1000011] = 431 ## SSEL+
_PDGID2HERWIGID[ 1000012] = 426 ## SSNUEL
_PDGID2HERWIGID[-1000012] = 432 ## SSNUELBR
_PDGID2HERWIGID[ 1000013] = 427 ## SSMUL-
_PDGID2HERWIGID[-1000013] = 433 ## SSMUL+
_PDGID2HERWIGID[ 1000014] = 428 ## SSNUMUL
_PDGID2HERWIGID[-1000014] = 434 ## SSNUMLBR
_PDGID2HERWIGID[ 1000015] = 429 ## SSTAU1-
_PDGID2HERWIGID[-1000015] = 435 ## SSTAU1+
_PDGID2HERWIGID[ 1000016] = 430 ## SSNUTL
_PDGID2HERWIGID[-1000016] = 436 ## SSNUTLBR
_PDGID2HERWIGID[ 2000011] = 437 ## SSEL-
_PDGID2HERWIGID[-2000011] = 443 ## SSEL+
_PDGID2HERWIGID[ 2000012] = 438 ## SSNUEL
_PDGID2HERWIGID[-2000012] = 444 ## SSNUELBR
_PDGID2HERWIGID[ 2000013] = 439 ## SSMUL-
_PDGID2HERWIGID[-2000013] = 445 ## SSMUL+
_PDGID2HERWIGID[ 2000014] = 440 ## SSNUMUL
_PDGID2HERWIGID[-2000014] = 446 ## SSNUMLBR
_PDGID2HERWIGID[ 2000015] = 441 ## SSTAU1-
_PDGID2HERWIGID[-2000015] = 447 ## SSTAU1+
_PDGID2HERWIGID[ 2000016] = 442 ## SSNUTL
_PDGID2HERWIGID[-2000016] = 448 ## SSNUTLBR
_PDGID2HERWIGID[ 1000021] = 449 ## GLUINO
_PDGID2HERWIGID[ 1000022] = 450 ## NTLINO1
_PDGID2HERWIGID[ 1000023] = 451 ## NTLINO2
_PDGID2HERWIGID[ 1000025] = 452 ## NTLINO3
_PDGID2HERWIGID[ 1000035] = 453 ## NTLINO4
_PDGID2HERWIGID[ 1000024] = 454 ## CHGINO1+
_PDGID2HERWIGID[-1000024] = 456 ## CHGINO1-
_PDGID2HERWIGID[ 1000037] = 455 ## CHGINO2+
_PDGID2HERWIGID[-1000037] = 457 ## CHGINO2-
_PDGID2HERWIGID[ 1000039] = 458 ## GRAVITINO


###############################################################################
## ISAWIG format reading/writing

def readISAWIG(isastr, ignorenobr=False):
    """
    Read a spectrum definition from a string in the ISAWIG format, returning
    dictionaries of blocks and decays. While this is not an SLHA format, it is
    informally supported as a useful mechanism for converting ISAWIG spectra to
    SLHA.

    ISAWIG parsing based on the HERWIG SUSY specification format, from
    http://www.hep.phy.cam.ac.uk/~richardn/HERWIG/ISAWIG/file.html

    If the ignorenobr parameter is True, do not store decay entries with a
    branching ratio of zero.
    """

    blocks = _dict("Blocks")
    decays = _dict("Decays")
    LINES = isastr.splitlines()

    def getnextvalidline():
        while LINES:
            s = LINES.pop(0).strip()
            # print "*", s, "*"
            ## Return None if EOF reached
            if len(s) == 0:
                continue
            ## Strip comments
            if "#" in s:
                s = s[:s.index("#")].strip()
            ## Return if non-empty
            if len(s) > 0:
                return s

    def getnextvalidlineitems():
        return [_autotype(x) for x in getnextvalidline().split()]

    ## Populate MASS block and create decaying particle objects
    masses = Block("MASS")
    numentries = int(getnextvalidline())
    for i in range(numentries):
        hwid, mass, lifetime = getnextvalidlineitems()
        width = 1.0/(lifetime * 1.51926778e24) ## width in GeV == hbar/lifetime in seconds
        pdgid = herwigid2pdgid(hwid)
        masses[pdgid] = mass
        decays[pdgid] = Particle(pdgid, width, mass)
        #print pdgid, mass, width
    blocks["MASS"] = masses

    ## Populate decays
    for n in range(numentries):
        numdecays = int(getnextvalidline())
        #print n+1, "/", numentries, numdecays
        for d in range(numdecays):
            #print n, numentries-1, d, numdecays-1
            decayitems = getnextvalidlineitems()
            hwid = decayitems[0]
            pdgid = herwigid2pdgid(hwid)
            br = decayitems[1]
            nme = decayitems[2]
            daughter_hwids = decayitems[3:]
            daughter_pdgids = []
            for hw in daughter_hwids:
                if hw != 0:
                    daughter_pdgids.append(herwigid2pdgid(hw))
            if pdgid not in decays:
                #print "Decay for unlisted particle %d, %d" % (hwid, pdgid)
                decays[pdgid] = Particle(pdgid)
            decays[pdgid].add_decay(br, len(daughter_pdgids), daughter_pdgids)


    ## Now the SUSY parameters
    TANB, ALPHAH = getnextvalidlineitems()
    blocks["MINPAR"] = Block("MINPAR")
    blocks["MINPAR"][3] = TANB
    blocks["ALPHA"] = Block("ALPHA")
    blocks["ALPHA"].set_value(ALPHAH)
    #
    ## Neutralino mixing matrix
    blocks["NMIX"] = Block("NMIX")
    for i in range(1, 5):
        nmix_i = getnextvalidlineitems()
        for j, v in enumerate(nmix_i):
            blocks["NMIX"][i, j+1] = v
    #
    ## Chargino mixing matrices V and U
    blocks["VMIX"] = Block("VMIX")
    vmix = getnextvalidlineitems()
    blocks["VMIX"][1, 1] = vmix[0]
    blocks["VMIX"][1, 2] = vmix[1]
    blocks["VMIX"][2, 1] = vmix[2]
    blocks["VMIX"][2, 2] = vmix[3]
    blocks["UMIX"] = Block("UMIX")
    umix = getnextvalidlineitems()
    blocks["UMIX"][1, 1] = umix[0]
    blocks["UMIX"][1, 2] = umix[1]
    blocks["UMIX"][2, 1] = umix[2]
    blocks["UMIX"][2, 2] = umix[3]
    #
    THETAT, THETAB, THETAL = getnextvalidlineitems()
    import math
    blocks["STOPMIX"] = Block("STOPMIX")
    blocks["STOPMIX"][1, 1] =  math.cos(THETAT)
    blocks["STOPMIX"][1, 2] = -math.sin(THETAT)
    blocks["STOPMIX"][2, 1] =  math.sin(THETAT)
    blocks["STOPMIX"][2, 2] =  math.cos(THETAT)
    blocks["SBOTMIX"] = Block("SBOTMIX")
    blocks["SBOTMIX"][1, 1] =  math.cos(THETAB)
    blocks["SBOTMIX"][1, 2] = -math.sin(THETAB)
    blocks["SBOTMIX"][2, 1] =  math.sin(THETAB)
    blocks["SBOTMIX"][2, 2] =  math.cos(THETAB)
    blocks["STAUMIX"] = Block("STAUMIX")
    blocks["STAUMIX"][1, 1] =  math.cos(THETAL)
    blocks["STAUMIX"][1, 2] = -math.sin(THETAL)
    blocks["STAUMIX"][2, 1] =  math.sin(THETAL)
    blocks["STAUMIX"][2, 2] =  math.cos(THETAL)
    #
    ATSS, ABSS, ALSS = getnextvalidlineitems()
    blocks["AU"] = Block("AU")
    blocks["AU"][3, 3] = ATSS
    blocks["AD"] = Block("AD")
    blocks["AD"][3, 3] = ABSS
    blocks["AE"] = Block("AE")
    blocks["AE"][3, 3] = ALSS
    #
    MUSS = getnextvalidlineitems()[0]
    blocks["MINPAR"][4] = MUSS
    #

    # TODO: Parse RPV boolean and couplings into SLHA2 blocks

    return Doc(blocks, decays)


def writeISAWIG(doc, ignorenobr=True, precision=8):
    """
    Return a SUSY spectrum definition in the format produced by ISAWIG for input to HERWIG
    as a string, from the supplied SLHA blocks and decays dicts.

    ISAWIG parsing based on the HERWIG SUSY specification format, from
    http://www.hep.phy.cam.ac.uk/~richardn/HERWIG/ISAWIG/file.html

    If the ignorenobr parameter is True, do not write decay entries with a
    branching ratio of zero.
    """

    blocks = doc.blocks
    decays = doc.decays

    masses = blocks["MASS"]

    ## Init output string
    out = ""

    ## First write out masses section:
    ##   Number of SUSY + top particles
    ##   IDHW, RMASS(IDHW), RLTIM(IDHW)
    ##   repeated for each particle
    ## IDHW is the HERWIG identity code.
    ## RMASS and RTLIM are the mass in GeV, and lifetime in seconds respectively.
    massout = ""
    for pid in masses.keys():
        lifetime = -1
        try:
            width = decays[pid].totalwidth
            if width and width > 0:
                lifetime = 1.0/(width * 1.51926778e24) ## lifetime in seconds == hbar/width in GeV
            if lifetime < 0: ## hack because HERWIG doesn't take negative lifetimes, but sets large positive ones stable
                lifetime = -lifetime
        except:
            pass
        if -pid in _PDGID2HERWIGID.keys():
            massout += ("%4d %11.4f %14.5e \n" % (_PDGID2HERWIGID.get(pid,pid), masses[pid], lifetime))
            massout += ("%4d %11.4f %14.5e \n" % (_PDGID2HERWIGID.get(-pid,-pid), masses[pid], lifetime))
        else:
            # print pid
            massout += ("%4d %11.4f %14.5e \n" % (_PDGID2HERWIGID.get(pid,pid), masses[pid], lifetime))
        # massout += ("%d " % _PDGID2HERWIGID.get(pid,pid)) + _autostr(masses[pid], precision) + " " + _autostr(lifetime, precision) + "\n"
    out += "%d\n" % massout.count("\n")
    out += massout

    # TODO: it's possible for there to be entries in DECAYS which aren't in MASS (sigh)
    # print len(masses), len(decays)
    # print zip(sorted(masses.keys()), sorted(decays.keys()))
    # assert(len(masses) == len(decays))

    ## Next each particles decay modes together with their branching ratios and matrix element codes
    ##   Number of decay modes for a given particle (IDK)
    ##     IDK(*), BRFRAC(*), NME(*) & IDKPRD(1-5,*)
    ##     repeated for each mode.
    ##   Repeated for each particle.
    ## IDK is the HERWIG code for the decaying particle, BRFRAC is the branching ratio of
    ## the decay mode. NME is a code for the matrix element to be used, either from the
    ## SUSY elements or the main HERWIG MEs. IDKPRD are the HERWIG identity codes of the decay products.
    for i, pid in enumerate(decays.keys()):
        if not pid in masses.keys():
            continue
        # if not decays.has_key(pid):
        #     continue
        hwid = _PDGID2HERWIGID.get(pid,pid)
        decayout = ""
        hasantiparticle = False
        if -pid in _PDGID2HERWIGID.keys():
            hasantiparticle = True
            hwidanti = _PDGID2HERWIGID.get(-pid,-pid)
            decayoutanti = ""
        #decayout += "@@@@ %d %d %d\n" % (i, pid, hwid)
        for i_d, d in enumerate(decays[pid].decays):
            ## Skip decay if it has no branching ratio
            if ignorenobr and d.br == 0:
                continue

            ## Identify decay matrix element to use
            ## From std HW docs, or from this pair:
            ## Two new matrix element codes have been added for these new decays:
            ##    NME =	200 	3 body top quark via charged Higgs
            ##    	300 	3 body R-parity violating gaugino and gluino decays
            nme = 0
            # TODO: Get correct condition for using ME 100... this guessed from some ISAWIG output
            if abs(pid) in (6, 12):
                nme = 100
            ## Extra SUSY MEs
            if len(d.ids) == 3:
                # TODO: How to determine the conditions for using 200 and 300 MEs? Enumeration of affected decays?
                pass
            decayout += ("%5d %15.8f %5d " % (hwid, d.br, nme))
            if hasantiparticle:
                decayoutanti += ("%5d %15.8f %5d " % (hwidanti, d.br, nme))
            # decayout += ("%d " % hwid) + _autostr(d.br, precision) + (" %d " % nme)

            def is_quark(pid):
                return (abs(pid) in range(1, 7))

            def is_lepton(pid):
                return (abs(pid) in range(11, 17))

            def is_squark(pid):
                if abs(pid) in range(1000001, 1000007):
                    return True
                if abs(pid) in range(2000001, 2000007):
                    return True
                return False

            def is_slepton(pid):
                if abs(pid) in range(1000011, 1000017):
                    return True
                if abs(pid) in range(2000011, 2000016, 2):
                    return True
                return False

            def is_gaugino(pid):
                if abs(pid) in range(1000022, 1000026):
                    return True
                if abs(pid) in (1000035, 1000037):
                    return True
                return False

            def is_susy(pid):
                return (is_squark(pid) or is_slepton(pid) or is_gaugino(pid) or pid == 1000021)

            absids = [abs(i) for i in d.ids]

            ## Order decay products as required by HERWIG
            ## Top
            if abs(pid) == 6:
                def key_bottomlast(item):
                    """Comparison function which always puts b/bar last"""
                    answer = (item == 5)*1e9 #< I never liked this ad-hocness in the 'key' approach
                    return answer if answer != 0 else item
                if len(absids) == 2:
                    ## 2 body mode, to Higgs: Higgs; Bottom
                    if (25 in absids or 26 in absids) and 5 in absids:
                        d.ids = sorted(d.ids, key=key_bottomlast)
                elif len(absids) == 3:
                    ## 3 body mode, via charged Higgs/W: quarks or leptons from W/Higgs; Bottom
                    if 37 in absids or 23 in absids:
                        d.ids = sorted(d.ids, key=key_bottomlast)
            ## Gluino
            elif abs(pid) == 1000021:
                if len(absids) == 2:
                    ## 2 body mode
                    ## without gluon: any order
                    ## with gluon: gluon; colour neutral
                    if 21 in absids:
                        def key_gluonfirst(item):
                            """Comparison function which always puts gluon first"""
                            answer = (item == 21)*-1e9 #< I never liked this ad-hocness in the 'key' approach
                            return answer if answer != 0 else item
                        d.ids = sorted(d.ids, key=key_gluonfirst)
                elif len(absids) == 3:
                    ## 3-body modes, R-parity conserved: colour neutral; q or qbar
                    def key_quarkslast(item):
                        """Comparison function which always puts quarks last"""
                        answer = is_quark(item)*1e9 #< I never liked this ad-hocness in the 'key' approach
                        return answer if answer != 0 else item
                    d.ids = sorted(d.ids, key=key_quarkslast)
            ## Squark/Slepton
            elif is_squark(pid) or is_slepton(pid):
                def key_susy_quark_lepton(item):
                    answer = (is_susy(item) + is_quark(item))*-1e9 #< I never liked this ad-hocness in the 'key' approach
                    return answer if answer != 0 else item

                ##   2 body modes: Gaugino/Gluino with Quark/Lepton     Gaugino      quark
                ##                                                      Gluino       lepton
                ##   3 body modes: Weak                                 sparticle    particles from W decay
                ## Squark
                ##   2 body modes: Lepton Number Violated               quark     lepton
                ##                 Baryon number violated               quark     quark
                ## Slepton
                ##   2 body modes: Lepton Number Violated               q or qbar
                # d.ids = sorted(d.ids, key=key_bottomlast)
                d.ids = sorted(d.ids, key=key_susy_quark_lepton)
            ## Higgs
            elif pid in (25, 26):
                # TODO: Includes SUSY Higgses?
                ## Higgs
                ##   2 body modes: (s)quark-(s)qbar                     (s)q or (s)qbar
                ##                 (s)lepton-(s)lepton                  (s)l or (s)lbar
                ##   3 body modes:                                      colour neutral       q or qbar
                if len(absids) == 3:
                    def key_quarkslast(item):
                        """Comparison function which always puts quarks last"""
                        answer = is_quark(item)*1e9 #< I never liked this ad-hocness in the 'key' approach
                        return answer if answer != 0 else item
                    d.ids = sorted(d.ids, key=key_quarkslast)
            elif is_gaugino(pid):
                # TODO: Is there actually anything to do here?
                ## Gaugino
                ##   2 body modes: Squark-quark                         q or sq
                ##                 Slepton-lepton                       l or sl
                ##
                ##   3 body modes: R-parity conserved                   colour neutral       q or qbar
                ##                                                                           l or lbar
                if len(absids) == 3:
                    def key_quarkslast(item):
                        """Comparison function which always puts quarks last"""
                        answer = is_quark(item)*1e9 #< I never liked this ad-hocness in the 'key' approach
                        return answer if answer != 0 else item
                    d.ids = sorted(d.ids, key=key_quarkslast)

            # TODO: Gaugino/Gluino
            ##   3 body modes:  R-parity violating:   Particles in the same order as the R-parity violating superpotential

            ## Pad out IDs list with zeros
            ids = [0,0,0,0,0]
            for i, pid in enumerate(d.ids):
                ids[i] = pid
            for id_i in ids:
                decayout += "%5d " % _PDGID2HERWIGID.get(id_i, id_i)
                if hasantiparticle:
                    id_anti_i = -id_i if -id_i in _PDGID2HERWIGID.keys() else id_i
                    decayoutanti += "%5d " % _PDGID2HERWIGID.get(id_anti_i, id_anti_i)
            decayout += "\n"
            if hasantiparticle:
                decayoutanti += "\n"
            # decayout += " ".join(ids) + "\n"
        ndecayout = decayout.count("\n")
        # if ndecayout:
        decayout = ("%d\n" % ndecayout) + decayout
        out += decayout
        if hasantiparticle:
            decayoutanti = ("%d\n" % ndecayout) + decayoutanti
            out += decayoutanti

    ## Now the SUSY parameters
    ## TANB, ALPHAH:
    out += ("%15.8f %15.8f \n" %(blocks["MINPAR"][3], blocks["ALPHA"].value()))
    # out += _autostr(blocks["MINPAR"][3], precision) + " " + _autostr(blocks["ALPHA"].value(), precision) + "\n"
    ## Neutralino mixing matrix
    nmix = blocks["NMIX"]
    for i in range(1, 5):
        out += ("%15.8f %15.8f %15.8f %15.8f \n" %(nmix[i,1], nmix[i,2], nmix[i,3], nmix[i,4]))
        # out += _autostr(nmix[i,1], precision) + " " + \
        #        _autostr(nmix[i,2], precision) + " " + \
        #        _autostr(nmix[i,3], precision) + " " + \
        #        _autostr(nmix[i,4], precision) + "\n"
    ## Chargino mixing matrices V and U
    vmix = blocks["VMIX"]
    out += ("%15.8f %15.8f %15.8f %15.8f \n" %(vmix[1,1], vmix[1,2], vmix[2,1], vmix[2,2]))
    # out += _autostr(vmix[1,1], precision) + " " + \
    #        _autostr(vmix[1,2], precision) + " " + \
    #        _autostr(vmix[2,1], precision) + " " + \
    #        _autostr(vmix[2,2], precision) + "\n"
    umix = blocks["UMIX"]
    out += ("%15.8f %15.8f %15.8f %15.8f \n" %(umix[1,1], umix[1,2], umix[2,1], umix[2,2]))
    # out += _autostr(umix[1,1], precision) + " " + \
    #        _autostr(umix[1,2], precision) + " " + \
    #        _autostr(umix[2,1], precision) + " " + \
    #        _autostr(umix[2,2], precision) + "\n"
    ## THETAT,THETAB,THETAL
    import math
    out += ("%15.8f %15.8f %15.8f \n" %(blocks["STOPMIX"][1,1], blocks["SBOTMIX"][1,1], blocks["STAUMIX"][1,1]))
    # out += _autostr(math.acos(blocks["STOPMIX"][1,1]), precision) + " " + \
    #        _autostr(math.acos(blocks["SBOTMIX"][1,1]), precision) + " " + \
    #        _autostr(math.acos(blocks["STAUMIX"][1,1]), precision) + "\n"
    ## ATSS,ABSS,ALSS
    out += ("%15.8f %15.8f %15.8f \n" %(blocks["AU"][3,3], blocks["AD"][3,3], blocks["AE"][3,3]))
    # out += _autostr(blocks["AU"][3,3], precision) + " " + \
    #        _autostr(blocks["AD"][3,3], precision) + " " + \
    #        _autostr(blocks["AE"][3,3], precision) + "\n"
    ## MUSS == sign(mu)
    out += "%15.8f\n" % blocks["MINPAR"][4]

    ## RPV SUSY
    # JEM EDIT: HERWIG reads the logical variable isRPC, not RPV. This is the opposite.
    # isRPV = False
    isRPC = True
    if isRPC:
        out += "  T  \n"
    else:
        out += "  F  \n"
    # TODO: Write RPV couplings if RPV is True (lambda1,2,3; 27 params in each, sci format.
    # TODO: Get the index orderings right
    # if isRPV: ...

    return out



###############################################################################
## File-level read/write functions


def read(spcfile, **kwargs):
    """
    Read an SLHA or ISAWIG file (or stdin).

    spcfile may either be a string filename or a file object.
    If a s string, the assumed file format is based from the
    filename; if a file it is assumed to be SLHA format.

    Other keyword parameters are passed to readSLHA/readISAWIG.
    """
    txt = _read(spcfile)
    if type(spcfile) is str and spcfile.endswith(".isa"):
        return readISAWIG(txt, **kwargs)
    else:
        return readSLHA(txt, **kwargs)

def readSLHAFile(spcfile, **kwargs):
    """
    Read an SLHA file, returning dictionaries of blocks and decays.

    Other keyword parameters are passed to readSLHA.
    """
    return readSLHA(_read(spcfile), **kwargs)

def readISAWIGFile(isafile, **kwargs):
    """
    Read a spectrum definition from a file in the ISAWIG format, returning
    dictionaries of blocks and decays. While this is not an SLHA format, it is
    informally supported as a useful mechanism for converting ISAWIG spectra to
    SLHA.

    isafile may either be a string filename or a file object.

    Other keyword parameters are passed to readISAWIG.
    """
    return readISAWIG(_read(isafile), **kwargs)


def write(spcfile, doc, **kwargs):
    """
    Write to an SLHA or ISAWIG file (or stdout).

    spcfile may either be a string filename or a file object.
    If a s string, the assumed file format is based from the
    filename; if a file it is assumed to be SLHA format.

    Other keyword parameters are passed to writeSLHA/writeISAWIG.
    """
    if type(spcfile) is str and spcfile.endswith(".isa"):
        writeISAWIGFile(spcfile, doc, **kwargs)
    else:
        writeSLHAFile(spcfile, doc, **kwargs)

def writeSLHAFile(spcfile, doc, **kwargs):
    """
    Write an SLHA file from the supplied blocks and decays dicts.

    Other keyword parameters are passed to writeSLHA.
    """
    _write(spcfile, writeSLHA(doc, **kwargs))

def writeISAWIGFile(isafile, doc, **kwargs):
    """
    Write an ISAWIG file from the supplied blocks and decays dicts.

    isafile may either be a string filename or a file object.

    Other keyword parameters are passed to writeISAWIG.
    """
    _write(isafile, writeISAWIG(doc, **kwargs))


###############################################################################
## Main function for module testing

if __name__ == "__main__":

    for a in sys.argv[1:]:
        doc = read(a)
        print(doc)
        print("")

        for bname, b in sorted(doc.blocks.items()):
            print(b)
            print("")

        print(doc.blocks.keys())

        print(doc.blocks["MASS"].get(25))
        print("")

        for p in sorted(doc.decays.values()):
            print(p)
            print("")

        print(writeSLHA(doc, ignorenobr=True))

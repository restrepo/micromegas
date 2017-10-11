##########################################################################
#
#  This file is part of Lilith
#  made by J. Bernon and B. Dumont
#
#  Web page: http://lpsc.in2p3.fr/projects-th/lilith/
#
#  In case of questions email bernon@lpsc.in2p3.fr 
#
#
#    Lilith is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Lilith is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Lilith.  If not, see <http://www.gnu.org/licenses/>.
#
##########################################################################

import sys
try:
    from lxml import etree
except:
    import xml.etree.ElementTree as etree
from ..errors import ExpInputError, ExpInputIOError
from scipy import interpolate
import numpy as np
import math
from . import brsm as BR_SM
from warnings import warn

class ReadExpInput:
    """Read the experimental input in XML and extracts all information."""

    def __init__(self):
        """Initialize the VBF, WH and ZH cross section ratios."""

        self.eff_VVH = BR_SM.geteffVVHfunctions()
        self.mu = []
        self.filepath = ""

    def warning(self, message):
        """Customized warnings."""

        warn("in the file " + self.filepath + ": " + message, Warning,
                      stacklevel=3)

    def get_filelist(self, filepath):
        """get list of files from .list experimental input"""

        expfiles = [] # list of all paths to XML input files

        filepath_split = filepath.split("/")
        expdata_dir = "/".join(filepath_split[:-1])

        try:
            with open(filepath) as linput:
                expfiles = [] # list of all paths to XML input files
                for line in linput:
                    # remove comment, space and new lines from the line
                    line = line.split("#")[0].rstrip("\n").strip()
                    if line == "": # empty line or comment
                        continue
                    if line[0] == "/": # absolute filepath
                        expfiles.append(line)
                    else: # relative filepath
                        expfiles.append(expdata_dir+"/"+line)
        except IOError as e:
            raise ExpInputIOError(
                'I/O error({0}): {1}'.format(e.errno, e.strerror) + '; cannot' +
                ' open the experimental list file "' + filepath + '".')
        
        return expfiles
    

    def read_file(self, filepath):
        """read individual xml files"""

        self.filepath = filepath

        root = self.produce_tree()
        if root.tag != "expmu":
            raise ExpInputError(self.filepath, "root tag is not <expmu>")

        (dim, decay, type) = self.get_mode(root)
        
        self.get_mass(root)
        (experiment, source, sqrts) = self.get_metadata(root)
        
        eff = self.read_eff(root, dim, decay)
        (bestfit, param, grid, Lxy, LChi2min) = self.read_mus(root, dim, type)

        self.mu.append({"filepath": self.filepath,
                        "dim": dim, "type": type,
                        "eff": eff,
                        "bestfit": bestfit, "param": param, "grid": grid,
                        "Lxy": Lxy, "LChi2min": LChi2min,
                        "experiment": experiment, "source": source,
                        "sqrts": sqrts, "eff": eff})

    def produce_tree(self):
        """Produce the XML tree with ElementTree."""

        try:
            with open(self.filepath) as f:
                tree = etree.parse(f)
        except IOError as e:
            raise ExpInputIOError(
                'I/O error({0}): {1}'.format(e.errno, e.strerror) + '; cannot' +
                'open the experimental input file "' + self.filepath + '".')

        return tree.getroot()

    def get_mode(self, root):
        """Get the dimension, decay and type of the experimental mu."""
        
        allowed_decays = ["gammagamma", "ZZ", "WW", "Zgamma",
                          "tautau", "bb", "cc", "mumu", "invisible"]

        mandatory_attribs = {"dim":["1", "2"],
                             "type":["n", "f"]}
        optional_attribs = {"decay": allowed_decays}

        for mandatory_attrib, allowed_values in mandatory_attribs.items():
            if mandatory_attrib not in root.attrib:
                # if "dim" or "type" not in attribute
                raise ExpInputError(self.filepath,
                                    'mandatory attribute of root tag "' +
                                    mandatory_attrib + '" is not present.')
            
            if root.attrib[mandatory_attrib] not in allowed_values:
                # if dim="3" or type="z" for instance
                raise ExpInputError(self.filepath,
                                    'mandatory attribute of root tag "' +
                                    mandatory_attrib + '" has value "' +
                                   root.attrib[mandatory_attrib] +
                                   '" which is unknown. Allowed values are : ' +
                                   str(allowed_values))

        dim = int(root.attrib["dim"])
        type = root.attrib["type"]

        decay = "mixture"

        for optional_attrib, allowed_values in optional_attribs.items():
            if optional_attrib in root.attrib:
                # if "decay" in attribute
                if root.attrib[optional_attrib] not in allowed_values:
                    # if decay="yy" for instance
                    raise ExpInputError(self.filepath,
                                        'optional attribute of root tag "' +
                                        optional_attrib + '" has value "' +
                                        root.attrib[optional_attrib] +
                                        '" which is unknown. Allowed values ' +
                                        + 'are: ' + str(allowed_values))
                else:
                    decay = root.attrib["decay"]

        return (dim, decay, type)
    
    def get_mass(self, root):
        def_mass = 125. # default value
        mass = def_mass
    
        for child in root:
            if child.tag == "mass":
                try:
                    mass = float(child.text)
                except TypeError: # empty tag is of type NULL
                    self.warning('<mass> tag is empty; ' +
                                 'setting the mass to ' + def_mass + ' GeV')
                    mass = def_mass
                except ValueError:
                    raise ExpInputError(self.filepath,
                                        "value of <mass> tag is not a number")
        self.mass = mass

    def get_metadata(self, root):
        experiment = ""
        source = ""
        sqrts = ""

        for child in root:
            if child.tag == "experiment":
                experiment = child.text
            if child.tag == "source":
                source = child.text
            if child.tag == "sqrts":
                sqrts = child.text
    
        return (experiment, source, sqrts)

    def read_eff(self, root, dim, decay):
        allowed_decays = ["gammagamma", "ZZ", "WW", "Zgamma",
                          "tautau", "bb", "cc", "mumu", "invisible"]

        # read the efficiencies
        if dim == 1: # 1D signal strength
            eff = {"x": {}}
            axis_label = "x"

            mandatory_attribs = {"prod": ["ggH", "VVH", "VBF", "VH", "WH", "ZH", "ttH"]}
            if decay == "mixture":
                mandatory_attribs["decay"] = allowed_decays

            for child in root:
                if child.tag == "eff":
                    for mandatory_attrib, allowed_values in mandatory_attribs.items():
                        if mandatory_attrib not in child.attrib:
                            # if "axis" or "prod" not in attribute
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" is not present.')
                        if child.attrib[mandatory_attrib] not in allowed_values:
                            # if axis="z" or prod="yy"
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" has value "' +
                                                child.attrib[mandatory_attrib] + '" which is unknown. Allowed values are : ' + str(allowed_values))

                    prod_label = child.attrib["prod"]
                    if decay == "mixture":
                        decay_label = child.attrib["decay"]
                    else:
                        decay_label = decay

                    if (prod_label,decay_label) in eff["x"]:
                        self.warning('<eff> tag with prod="' + prod_label +
                                     '" and decay="' + decay_label +
                                     '" is being redefined.')
                        
                    try:
                        eff[axis_label][prod_label,decay_label] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        self.warning('<eff> tag for axis="' + axis_label +
                                     '", prod="' + prod_label + '" and decay="' +
                                     decay_label + '" is empty; setting to ' +
                                     'default value of 0')
                        eff[axis_label][prod_label,decay_label] = 0.
                    except ValueError:
                        raise ExpInputError(self.filepath,
                                            'value of <eff> tag with axis="' + axis_label +
                                            '" and prod="' + prod_label + '" and decay="' + decay_label + '" is not a number')

        else: # 2D signal strength
            eff = {"x": {}, "y": {}}

            mandatory_attribs = {"axis": ["x", "y"],
                             "prod": ["ggH", "VVH", "VBF", "VH", "WH", "ZH", "ttH"]}
            if decay == "mixture":
                mandatory_attribs["decay"] = allowed_decays
            
            for child in root:
                if child.tag == "eff":
                    for mandatory_attrib, allowed_values in mandatory_attribs.items():
                        if mandatory_attrib not in child.attrib:
                            # if "axis" or "prod" not in attribute
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" is not present.')
                        if child.attrib[mandatory_attrib] not in allowed_values:
                            # if axis="z" or prod="yy"
                            raise ExpInputError(self.filepath,
                                                'mandatory attribute of <eff> tag "' +
                                                mandatory_attrib + '" has value "' +
                                                child.attrib[mandatory_attrib] + '" which is unknown. Allowed values are : ' + str(allowed_values))

                        axis_label = child.attrib["axis"]
                        prod_label = child.attrib["prod"]
                        if decay == "mixture":
                            decay_label = child.attrib["decay"]
                        else:
                            decay_label = decay


                    if (prod_label,decay_label) in eff[axis_label]:
                        self.warning('<eff> tag with axis="' + axis_label +
                                     '", prod="' + prod_label +
                                     '" and decay="' +decay_label +
                                     '" is being redefined.')
                        
                    try:
                        eff[axis_label][prod_label,decay_label] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        self.warning('<eff> tag for axis="' + axis_label +
                                     '", prod="' + prod_label + '" and decay="' +
                                     decay_label + '" is empty; setting to ' +
                                     'default value of 0')
                        eff[axis_label][prod_label,decay_label] = 0.
                    except ValueError:
                        raise ExpInputError(self.filepath,
                                            'value of <eff> tag with axis="' + axis_label +
                                            '" and prod="' + prod_label + '" and decay="' + decay_label + '" is not a number')

        effWH = self.eff_VVH["eff_WH"](self.mass)
        effZH = self.eff_VVH["eff_ZH"](self.mass)
        
        effVBF = self.eff_VVH["eff_VBF"](self.mass)
        effVH = self.eff_VVH["eff_VH"](self.mass)
        effVWH = effVH * effWH
        effVZH = effVH * effZH

        multiprod = {"VH": {"WH": effWH, "ZH": effZH}, "VVH": {"VBF": effVBF, "WH": effVWH, "ZH": effVZH}}

        self.check_multiprod(eff["x"], multiprod)
        if dim == 2:
            self.check_multiprod(eff["y"], multiprod)

        # now all reduced couplings have been properly defined, one can
        # delete all multiparticle labels
        effCleanX = eff["x"].copy()
        for (p,decay) in eff["x"]:
            if p in multiprod:
                del effCleanX[p,decay]

        if dim == 2:
            effCleanY = eff["y"].copy()

            for (p,decay) in eff["y"]:
                if p in multiprod:
                    del effCleanY[p,decay]

        eff["x"] = effCleanX
        if dim == 2:
            eff["y"] = effCleanY

        # check that efficiencies add up to 1, otherwise issue a warning
        # or an error
        for axis in eff:
            sumeff = 0
            for prod in eff[axis]:
                sumeff += eff[axis][prod]
            
            if sumeff == 0:
                raise ExpInputError(self.filepath,
                                    "no <eff> tag found for " + axis + " axis")

            if sumeff < 0.99:
                self.warning('the sum of efficiencies for axis="' +
                             axis + '" is less than 1 (value: ' +
                             str(sumeff) + ')')
            elif sumeff > 1.01:
                raise ExpInputError(self.filepath,
                                    'the sum of efficiencies for axis="' +
                                    axis + '" is greater than 1 (value: ' +
                                    str(sumeff) + ')')

        return eff


    def read_mus(self, root, dim, type):
        # first, read the bestfit
        bestfit = {}
        LChi2min = 0
                
        for child in root:
            if child.tag == "bestfit":
                if type == "f":
                    self.warning('block <bestfit> in experimental mu of ' +
                                 'type "full"... skipping.')
                    continue

                if dim == 1:
                    # read directly the number
                    
                    if "x" in bestfit:
                        self.warning("redefinition of the bestfit...")
                    
                    try:
                        bestfit["x"] = float(child.text)
                    except TypeError: # empty tag is of type NULL
                        self.warning('<x> tag in <bestfit> block is empty; ' +
                                     'setting to 0')
                        bestfit["x"] = 0.
                    except ValueError:
                        raise ExpInputError(self.filepath,
                                            "value of <besfit> tag is not a number")
            
                elif dim == 2:
                    bestfit_allowedsubtags = ["x", "y"]
                
                    for bfit in child:
                        if bfit.tag in bestfit_allowedsubtags:
                            
                            if bfit.tag in bestfit:
                                self.warning("redefinition of the bestfit...")
                        
                            try:
                                bestfit[bfit.tag] = float(bfit.text)
                            except TypeError: # empty tag is of type NULL
                                self.warning('<' + bfit.tag + '> tag in ' +
                                             '<bestfit> block is empty; ' +
                                             'setting to 0')
                                bestfit[bfit.tag] = 0.
                            except ValueError:
                                raise ExpInputError(self.filepath,
                                                    "value of <besfit> tag is not a number")
                        else:
                            raise ExpInputError(self.filepath,
                                                "subtag in bestfit not known")
                    
                if dim == 1 and "x" not in bestfit:
                    raise ExpInputError(self.filepath,
                                        "best fit point should be specified.")
                if dim == 2 and ("x" not in bestfit or "y" not in bestfit):
                    raise ExpInputError(self.filepath,
                                        "best fit point should be specified for x and y.")
        
        # then, read the param...
        param = {}
        
        for child in root:
            if child == "param":
                break
        param_tag = child

        param["uncertainty"] = {}
        
        for child in param_tag:
            if child.tag is etree.Comment:
                # ignore all comments
                continue

            if dim == 1:
                if child.tag == "uncertainty":
                    if "side" not in child.attrib:
                        try:
                            unc_value = float(child.text)
                        except TypeError: # empty tag is of type NULL
                            self.warning('<uncertainty> tag is empty; ' +
                                         'setting to 0')
                            unc_value = 0.
                        except ValueError:
                            raise ExpInputError(self.filepath,
                                                "value of <uncertainty> tag is not a number")

                        param["uncertainty"]["left"] = unc_value
                        param["uncertainty"]["right"] = unc_value
                    else:
                        if child.attrib["side"] not in ["left", "right"]:
                            raise ExpInputError(self.filepath,
                                                "attribute of uncertainty is not left nor right")
                        else:
                            try:
                                unc_value = float(child.text)
                            except TypeError: # empty tag is of type NULL
                                self.warning('<uncertainty> tag is empty; ' +
                                             'setting to 0')
                                unc_value = 0.
                            except ValueError:
                                raise ExpInputError(self.filepath,
                                                    "value of <uncertainty> tag is " +
                                                    "not a number")
                                
                        param["uncertainty"][child.attrib["side"]] = unc_value
                else:
                    raise ExpInputError(self.filepath,
                                        "subtag or param should be uncertainty")
                    
            elif dim == 2:
                allowed_tags = ["a", "b", "c"]
                if child.tag not in allowed_tags:
                    raise ExpInputError(self.filepath,
                                        "only allowed tags are <a>, <b> and <c> in " +
                                        "block param in 2D normal mode")

                if child.tag in param:
                    self.warning("redefinition of tag <" + child.tag + ">")
                    
                try:
                    param_value = float(child.text)
                except TypeError: # empty tag is of type NULL
                    self.warning('<' + child.tag + '> tag is empty; ' +
                                 'setting to 0')
                    param_value = 0.
                except ValueError:
                    raise ExpInputError(self.filepath,
                                        "value of <" + child.tag + "> tag is not a number")
                
                param[child.tag] = param_value
        
        # check that everything is there
        if type == "n" and dim == 1:
            if ("uncertainty" not in param or
                "left" not in param["uncertainty"] or
                "right" not in param["uncertainty"]):
                raise ExpInputError(self.filepath,
                                    "uncertainties are not given consistently in block param")
        elif type == "n" and dim == 2:
            if "a" not in param or "b" not in param or "c" not in param:
                raise ExpInputError(self.filepath,
                                    "a, b, c tags are not given in block param")
        
        # or the grid
        grid = {}
        Lxy = None
        
        for child in root:
            if child == "grid":
                break
        grid_raw = child.text

        if type == "f" and dim == 1:
            x = []
            L = []
        
            grid_raw = grid_raw.strip("\n").strip().split("\n")
        
            i = -1
        
            for line in grid_raw:
                tab = line.split()
                if len(tab) != 2:
                    raise ExpInputError(self.filepath,
                                        'incorrect <grid> entry on line "' + line + '"')

                cur_x = float(tab[0])
                cur_L = float(tab[1])

                if cur_x not in x:
                    x.append(cur_x)
                    L.append(cur_L)
                    i += 1
                else:
                    i = x.index(cur_x)

            grid["x"] = x
            grid["L"] = L
            LChi2min = min(grid["L"])

            Lxy = interpolate.UnivariateSpline(grid["x"], grid["L"], k = 3, s = 0)
        
        elif type == "f" and dim == 2:
            x = []
            y = []
            L = []
            
            grid_raw = grid_raw.strip("\n").strip().split("\n")
            
            i = -1
            
            for line in grid_raw:
                tab = line.split()
                if len(tab) != 3:
                    raise ExpInputError(self.filepath,
                                        'incorrect <grid> entry on line "' + line + '"')

                cur_x = float(tab[0])
                cur_y = float(tab[1])
                cur_L = float(tab[2])

                if cur_x not in x:
                    x.append(cur_x)
                    L.append([])
                    i += 1
                else:
                    i = x.index(cur_x)

                if cur_y not in y:
                    y.append(cur_y)

                L[i].append(cur_L)

            grid["x"] = np.array(x)
            grid["y"] = np.array(y)
            grid["L"] = np.array(L)
            
            LChi2min = min(min(p[1:]) for p in grid["L"])

            Lxy = interpolate.RectBivariateSpline(grid["x"],
                    grid["y"], grid["L"])
            
        return (bestfit, param, grid, Lxy, LChi2min)

    def check_multiprod(self, eff_dict, multiprod):
        """..."""

        # check consistency in the definition of multi-particle labels
        for (prod,decay) in eff_dict:
            if prod in multiprod:
                # there is a multi-particle label
                # in that case, check if individual particle are also defined
                for label in multiprod[prod]:
                    if (label,decay) in eff_dict:
                        raise ExpInputError(self.filepath,
                                            '<eff> tags for "' + label + '" and "' +
                                            prod + '" cannot both be defined')
                # also, only one multi-particle label can be used (VH or VVH),
                # not both
                for label in multiprod:
                    if label != prod and (label,decay) in eff_dict:
                        raise ExpInputError(self.filepath,
                                            '<eff> tags for "' + label + '" and "' +
                                             prod + '" cannot both be defined')

        # it is consistent, resolve multi-particle labels
        new_eff = {}
        for (prod,decay) in eff_dict:
            if prod in multiprod:
                for label in multiprod[prod]:
                    new_eff[label,decay] = eff_dict[prod,decay]*multiprod[prod][label]

        for elem in new_eff:
            eff_dict[elem] = new_eff[elem]

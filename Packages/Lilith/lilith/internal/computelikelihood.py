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

from ..errors import LikelihoodComputationError

def compute_likelihood(exp_mu, user_mu):
    """Computes the likelihood from experimental mu and user mu."""

    likelihood_results = []
    l = 0. # actually -2log(likelihood)

    for mu in exp_mu:
        # compute user mu value scaled to efficiencies
        user_mu_effscaled = {}

        try:
            user_mu_effscaled["x"] = 0.
            for (prod,decay),eff_prod in mu["eff"]["x"].items():
                user_mu_effscaled["x"] += eff_prod*user_mu[prod,decay]

            if mu["dim"] == 2:
                user_mu_effscaled["y"] = 0.
                for (prod,decay),eff_prod in mu["eff"]["y"].items():
                    user_mu_effscaled["y"] += eff_prod*user_mu[prod,decay]
        except KeyError as s:
            if "s" in ["eff", "x", "y"]:
                # the experimental mu dictionnary is not filled correctly
                raise LikelihoodComputationError(
                    'there are missing elements in exp_mu: key "' + str(s) +
                    '" is not found')
            else:
                # the total user mu dictionnary is not filled correctly
                raise LikelihoodComputationError(
                    'there are missing elements in user_mu_tot: key "' +
                    str(s) + '" is not found')

        try:
            # likelihood computation in case of a type="normal" (Gaussian)
            if mu["type"] == "n":
                if mu["dim"] == 1:
                    if user_mu_effscaled["x"] < mu["bestfit"]["x"]:
                        unc = mu["param"]["uncertainty"]["left"]
                    else:
                        unc = mu["param"]["uncertainty"]["right"]
                    cur_l = ((mu["bestfit"]["x"] - user_mu_effscaled["x"])**2
                             / unc**2)

                elif mu["dim"] == 2:
                    a = mu["param"]["a"]
                    b = mu["param"]["b"]
                    c = mu["param"]["c"]

                    cur_l = a*(mu["bestfit"]["x"] - user_mu_effscaled["x"])**2
                    cur_l += c*(mu["bestfit"]["y"] - user_mu_effscaled["y"])**2
                    cur_l += (2*b*(mu["bestfit"]["x"] - user_mu_effscaled["x"])
                             * (mu["bestfit"]["y"] - user_mu_effscaled["y"]))

            if mu["type"] == "f":
                if mu["dim"] == 1:
                    cur_l = mu["Lxy"](user_mu_effscaled["x"]) - mu["LChi2min"]
                elif mu["dim"] == 2:
                    cur_l = (mu["Lxy"](user_mu_effscaled["x"],
                                       user_mu_effscaled["y"])[0][0]
                            - mu["LChi2min"])
        except KeyError as s:
            raise LikelihoodComputationError(
                'there are missing elements in exp_mu: key "' + s +
                '" is not found')

        l += cur_l
        likelihood_results.append(
            {"experiment": mu["experiment"], "source": mu["source"],
            "sqrts": mu["sqrts"], "dim": mu["dim"],
            "type": mu["type"], "eff": mu["eff"], "l": cur_l})

    return likelihood_results, l


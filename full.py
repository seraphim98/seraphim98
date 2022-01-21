import ytree
import math as m
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.stats import beta as beta_fitting
from scipy import stats
sys.path.insert(1, '/disk01/jhorlock')
from Halo import Halo
import random

print("version 27012021full")

def graph(data):

    """" Plots the Histogram of the data retrieved from merger tree data"""

    counts, bins, plot = plt.hist(data, bins='auto', orientation='horizontal', histtype='step')
    # arguments are passed to np.histogram
    bins = bins[:-1]
    plt.clf()
    return counts, bins


def mass_calc(halo_ancestors):

    masses = []
    x_list = []

    for i in range(len(halo_ancestors)):
        masses.append(halo_ancestors[i].node['mass'].to_value())
    for i in range(len(masses)):
        x_list.append(masses[i] / np.sum(masses))

    return x_list


def h_calc(x_list, alpha, f):

    calc_list = []

    for i in range(len(x_list)):
        calc_list.append((x_list[i] ** alpha) * np.log(x_list[i]))

    return -f * np.sum(calc_list)


def entropy_calc(halo, x_list, H, ancestors, a, b, gamma):

    sum_term = []
    masses = []
    halo_mass = halo.node['mass'].to_value()
    for i in range(len(ancestors)):
        calc = (x_list[i]**2) * (ancestors[i].entropy - H)
        sum_term.append(calc)
        masses.append(ancestors[i].node['mass'].to_value())

    entropy = H + (1 + a*H + b*H**2)*np.sum(sum_term)
    delta_m = halo_mass - np.sum(masses)
    snew = entropy * (np.sum(masses) / halo_mass)**gamma
    if snew > 1:
        pass
    else:
        entropy = snew

    halo.entropy_update(entropy)
    halo.accrete_mass_update(delta_m)


def one_ancestor(halo, ancestor, gamma):

    halo_mass = halo.node['mass'].to_value()
    mass = ancestor.node['mass'].to_value()
    entropy = ancestor.entropy
    delta_m = halo_mass - np.sum(mass)
    snew = entropy * (mass/halo_mass) ** gamma
    if snew > 1:
        pass
    else:
        entropy = snew

    halo.entropy_update(entropy)
    halo.accrete_mass_update(delta_m)

def main():
    treefile = "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/rockstar/trees/tree_0_0_0.dat"
    tree = ytree.load(treefile)
    fn = tree.save_arbor()
    tree = ytree.load(fn)
    halos_0 = tree.select_halos('tree["tree", "redshift"] == 0 ')
    halos = random.sample(list(halos_0), 10000)

    nc = 2
    beta = 0.75
    gamma = 1/3
    alpha = 1 / np.log(nc) + 1
    f = (alpha - 1) * m.exp(1)
    a = (2 - gamma) / f
    b = nc * (1 - beta) - 1 - a
    less = 0
    other = 0

    initial_entropy = 0
    generations = []
    generations_id = []
    no_sat = []
    calc_list = []
    num_gen = len(list(halos_0[0]['prog']))
    no_res = 0
    res = 0
    resolution = 10 ** 10
    indices = []

    for i in range(num_gen):
        generations.append([])
        generations_id.append([])

    for i in halos:
        prog_list = list(i['prog'])
        index = num_gen - len(prog_list)
        index2 = len(prog_list) - 1
        generations[index].append(Halo(prog_list[index2], initial_entropy, prog_list[index2]['mass'].to_value()))
        generations_id[index].append(prog_list[index2]['uid'])

    for i in range(len(generations)):

        test_list = []
        calc_gen = generations[i]
        if i == 0:
            initial_redshift = calc_gen[0].node['redshift']
        if i != (len(generations) - 1):
            for j in range(len(calc_gen)):
                halo_calc = calc_gen[j].node.descendent
                if j == 0:
                    if 4.95 < halo_calc['redshift'] < 5.05:
                        indices.append(i + 1)
                    if 3.95 < halo_calc['redshift'] < 4.05:
                        indices.append(i + 1)
                    if 2.95 < halo_calc['redshift'] < 3.05:
                        indices.append(i + 1)
                    if 1.95 < halo_calc['redshift'] < 2.05:
                        indices.append(i + 1)
                    if 0.95 < halo_calc['redshift'] < 1.05:
                        indices.append(i + 1)

                if halo_calc is None:
                    no_res = no_res + 1
                else:
                    res = res + 1
                    if halo_calc['uid'] not in test_list:

                        halo_calc_ancestors = list(halo_calc.ancestors)

                        if len(halo_calc_ancestors) == 1:
                            k = halo_calc_ancestors[0]
                            x = generations_id[i].index(k['uid'])
                            initial_mass = calc_gen[x].initial_mass
                            halo_object = Halo(halo_calc, initial_entropy, initial_mass)
                            one_ancestor(halo_object, calc_gen[x], gamma)
                            generations[i + 1].append(halo_object)
                            generations_id[i + 1].append(halo_object.node['uid'])
                            less = less + 1
                        else:
                            other = other + 1

                            test_list.append(halo_calc['uid'])

                            for k in halo_calc_ancestors:
                                if k['uid'] not in generations_id[i]:
                                    calc_list.append(Halo(k, initial_entropy, k['mass'].to_value()))
                                else:
                                    x = generations_id[i].index(k['uid'])
                                    calc_list.append(calc_gen[x])

                            initial_mass = calc_list[0].initial_mass
                            halo_object = Halo(halo_calc, initial_entropy, initial_mass)
                            x_list = mass_calc(calc_list)
                            H = h_calc(x_list, alpha, f)
                            entropy_calc(halo_object, x_list, H, calc_list, a, b, gamma)
                            generations[i + 1].append(halo_object)
                            generations_id[i + 1].append(halo_object.node['uid'])

        calc_list = []
        print("Calculations...")
        print(str(i * 100 / len(generations)) + '%')
    print("Calculations...")
    print(str(100) + '%')

    no_sat = []
    no_sat_small = []
    no_sat_large = []
    no_sat_10 = []
    no_sat_11 = []
    no_sat_12 = []
    no_sat_13 = []
    tot_ang_small = []
    tot_ang_large = []
    mass_10 = []
    tot_ang_10 = []
    mass_11 = []
    tot_ang_11 = []
    mass_12 = []
    tot_ang_12 = []
    mass_13 = []
    tot_ang_13 = []


    for i in generations[num_gen - 1]:

        ancestors = list(i.node.ancestors)
        x = generations_id[num_gen - 1].index(i.node['uid'])
        halo = generations[num_gen - 1][x]
        halo.mass_update(halo.node['mass'].to_value())
        N = halo.ratio
        if len(ancestors) > 1:

            if resolution <= halo.node['mass'].to_value() <= 10 * resolution:
                no_sat.append(halo.entropy)
                no_sat_10.append(halo.entropy)

                ang_mom = (0.6950 ** 2) * (np.sqrt(
                    halo.node['angular_momentum_x'] ** 2 + halo.node['angular_momentum_y'] ** 2 + halo.node[
                        'angular_momentum_z'] ** 2))
                tot_ang_10.append(ang_mom)
                mass_10.append(halo.node['mass'])
            if 10 * resolution <= halo.node['mass'].to_value() <= 100 * resolution:
                no_sat.append(halo.entropy)
                no_sat_11.append(halo.entropy)
                ang_mom = (0.6950 ** 2) * (np.sqrt(
                    halo.node['angular_momentum_x'] ** 2 + halo.node['angular_momentum_y'] ** 2 + halo.node[
                        'angular_momentum_z'] ** 2))
                tot_ang_11.append(ang_mom)
                mass_11.append(halo.node['mass'])

            if 100 * resolution <= halo.node['mass'].to_value() <= 1000 * resolution:
                no_sat.append(halo.entropy)
                no_sat_12.append(halo.entropy)
                ang_mom = (0.6950 ** 2) * (np.sqrt(
                    halo.node['angular_momentum_x'] ** 2 + halo.node['angular_momentum_y'] ** 2 + halo.node[
                        'angular_momentum_z'] ** 2))
                tot_ang_12.append(ang_mom)
                mass_12.append(halo.node['mass'])
            if 1000 * resolution <= halo.node['mass'].to_value() <= 10000 * resolution:
                no_sat.append(halo.entropy)
                no_sat_13.append(halo.entropy)
                ang_mom = (0.6950 ** 2) * (np.sqrt(
                    halo.node['angular_momentum_x'] ** 2 + halo.node['angular_momentum_y'] ** 2 + halo.node[
                        'angular_momentum_z'] ** 2))
                tot_ang_13.append(ang_mom)
                mass_13.append(halo.node['mass'])
            if N < 0.5 and N != 0:
                no_sat_small.append(halo.entropy)

                ang_mom = (0.6950 ** 2) * (np.sqrt(
                    halo.node['angular_momentum_x'] ** 2 + halo.node['angular_momentum_y'] ** 2 + halo.node[
                        'angular_momentum_z'] ** 2))
                tot_ang_small.append(ang_mom)
            if N >= 0.5:
                no_sat_large.append(halo.entropy)

                ang_mom = (0.6950 ** 2) * (np.sqrt(
                    halo.node['angular_momentum_x'] ** 2 + halo.node['angular_momentum_y'] ** 2 + halo.node[
                        'angular_momentum_z'] ** 2))
                tot_ang_large.append(ang_mom)

    en, no = graph(no_sat)

    plt.plot(no, en / np.mean(en))

    plt.title("Number of central halos at z = 0 against tree entropy")
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Number density ($Mpc^{-3}$)")
    plt.savefig("/disk01/jhorlock/sample/centralsample.png")
    plt.clf()

    plt.scatter(no_sat_10, np.log(mass_10))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("ln M")
    plt.savefig("/disk01/jhorlock/sample/mass_10sample.png")
    plt.clf()

    plt.scatter(no_sat_10, np.log(tot_ang_10))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Total angular momentum")
    plt.savefig("/disk01/jhorlock/sample/ang_mom_10sample.png")
    plt.clf()

    plt.scatter(no_sat_11, np.log(mass_11))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("ln M")
    plt.savefig("/disk01/jhorlock/sample/mass_11sample.png")
    plt.clf()

    plt.scatter(no_sat_11, np.log(tot_ang_11))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Total angular momentum")
    plt.savefig("/disk01/jhorlock/sample/ang_mom_11sample.png")
    plt.clf()

    plt.scatter(no_sat_12, np.log(mass_12))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("ln M")
    plt.savefig("/disk01/jhorlock/sample/mass_12sample.png")
    plt.clf()

    plt.scatter(no_sat_12, np.log(tot_ang_12))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Total angular momentum")
    plt.savefig("/disk01/jhorlock/sample/ang_mom_12sample.png")
    plt.clf()

    plt.scatter(no_sat_13, np.log(mass_13))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("ln M")
    plt.savefig("/disk01/jhorlock/sample/mass_13sample.png")
    plt.clf()

    plt.scatter(no_sat_13, np.log(tot_ang_13))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Total angular momentum")
    plt.savefig("/disk01/jhorlock/sample/ang_mom_13sample.png")
    plt.clf()

    plt.scatter(no_sat_small, np.log(tot_ang_small))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("tot_ang_small")
    plt.savefig("/disk01/jhorlock/sample/tot_ang_mergesample.png")
    plt.clf()
    plt.scatter(no_sat_large, np.log(tot_ang_large))
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Total angular momentum_large")
    plt.savefig("/disk01/jhorlock/sample/tot_ang_accretesample.png")
    plt.clf()

    a, b, loc, scale = beta_fitting.fit(no_sat)
    plt.plot(np.linspace(0, 1, 100), stats.beta.pdf(np.linspace(0, 1, 100), a, b, loc, scale))
    plt.hist(no_sat, bins='auto', histtype='step', density='true')
    plt.title("Number density of central halos at z = 0 against tree entropy")
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Normalised number density")
    plt.savefig("/disk01/jhorlock/sample/beta_fit_central_fullsample.png")
    plt.clf()

    no_sat_res = []
    no_sat = []
    halo_res = []
    all_halos = []

    # loops through final generation of halo objects

    for i in generations[num_gen - 1]:

        ancestors = list(i.node.ancestors)
        # gets entropy values for all halos

        all_halos.append(i.entropy)

        # checks length of ancestors to determine central halos
        if len(ancestors) > 1:

            no_sat.append(i.entropy)
            # checks for masses above 10^11
            if i.node['mass'].to_value() >= 10 ** 11:
                no_sat_res.append(i.entropy)
        # checks for all halos with mass > 10^11
        if i.node['mass'].to_value() >= 10 ** 11:
            halo_res.append(i.entropy)

    a, b, loc, scale = beta_fitting.fit(no_sat)
    plt.plot(np.linspace(0, 1, 100), stats.beta.pdf(np.linspace(0, 1, 100), a, b, loc, scale))
    plt.hist(no_sat, bins='auto', histtype='step', density='true')
    plt.title("Number of central halos at z = 0 against tree entropy")
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Normalised number density")
    plt.savefig("/disk01/jhorlock/sample/central_testsample.png")
    plt.clf()

    a, b, loc, scale = beta_fitting.fit(no_sat_res)
    plt.plot(np.linspace(0, 1, 100), stats.beta.pdf(np.linspace(0, 1, 100), a, b, loc, scale))
    plt.hist(no_sat_res, bins='auto', histtype='step', density='true')
    plt.title("Number of central halos with masses > 10 ^11 at z = 0 against tree entropy")
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Normalised number density")
    plt.savefig("/disk01/jhorlock/sample/central_res_testsample.png")
    plt.clf()

    a, b, loc, scale = beta_fitting.fit(halo_res)
    plt.plot(np.linspace(0, 1, 100), stats.beta.pdf(np.linspace(0, 1, 100), a, b, loc, scale))
    plt.hist(halo_res, bins='auto', histtype='step', density='true')
    plt.title("Number of halos with masses > 10 ^11 at z = 0 against tree entropy")
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Normalised number density")
    plt.savefig("/disk01/jhorlock/sample/res_testsample.png")
    plt.clf()

    a, b, loc, scale = beta_fitting.fit(all_halos)
    plt.plot(np.linspace(0, 1, 100), stats.beta.pdf(np.linspace(0, 1, 100), a, b, loc, scale))
    plt.hist(all_halos, bins='auto', histtype='step', density='true')
    plt.title("Number of halos at z = 0 against tree entropy")
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Normalised number density")
    plt.savefig("/disk01/jhorlock/sample/testsample.png")
    plt.clf()

    plot_x = []
    plot_y = []
    plot_x_no_sat = []
    plot_y_no_sat = []
    z_list = []
    print(len(indices))
    print("length of last generation is")
    print(len(generations[num_gen - 1]))

    for i in indices:
        plot_list = []
        plot_list_no_sat = []
        current_list = generations[i]
        for j in current_list:
            if j.node['mass'].to_value() >= 10 * resolution:
                plot_list.append(j.entropy)
                if len(list(j.node.ancestors)) > 1:
                    plot_list_no_sat.append(j.entropy)

        en, no = graph(plot_list)
        plot_x.append(en / 100 ** 3)
        plot_y.append(no)
        en, no = graph(plot_list_no_sat)
        plot_x_no_sat.append(en / 100 ** 3)
        plot_y_no_sat.append(no)
        z_list.append(current_list[0].node['redshift'])





    plot_list = []
    plot_list_no_sat = []
    for j in generations[num_gen - 1]:
        if j.node['mass'].to_value() >= 10 * resolution:
            plot_list.append(j.entropy)
            if len(list(j.node.ancestors)) > 1:
                plot_list_no_sat.append(j.entropy)

    en, no = graph(plot_list)
    plot_x.append(en / 100 ** 3)
    plot_y.append(no)
    en, no = graph(plot_list_no_sat)
    plot_x_no_sat.append(en / 100 ** 3)
    plot_y_no_sat.append(no)
    z_list.append(0)
    a, b, loc, scale = beta_fitting.fit(plot_list)
    plt.plot(np.linspace(0, 1, 100), stats.beta.pdf(np.linspace(0, 1, 100), a, b, loc, scale))
    plt.hist(plot_list, bins='auto', histtype='step', density='true')
    plt.ylabel("Normalised Number Density")
    plt.xlabel("Tree Entropy s")
    plt.title("Number density of halos against tree entropy")
    plt.savefig("/disk01/jhorlock/sample/beta_fit_fullsample.png")
    plt.clf()
    for i in range(len(plot_x)):
        redshift = round(z_list[i], 1)
        plt.plot(plot_y[i], plot_x[i] / 100 ** 3, label="z = {}".format(redshift))

    plt.title("Number density of halos against tree entropy at different redshifts")
    plt.legend()
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Number density ($Mpc^{-3}$)")
    plt.savefig("/disk01/jhorlock/sample/spline_fullsample.png")

    plt.clf()

    for i in range(len(plot_x_no_sat)):
        redshift = round(z_list[i], 1)
        plt.plot(plot_y_no_sat[i], plot_x_no_sat[i] / 100 ** 3, label="z = {}".format(redshift))
    plt.title("Number density of central halos against tree entropy at different redshifts")
    plt.legend()
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Number density ($Mpc^{-3}$)")
    plt.savefig("/disk01/jhorlock/sample/spline_central_fullsample.png")

    plt.clf()

    average_entropy = []
    average_entropy_no_sat = []
    for i in range(len(generations)):

        calc = []
        calc_no_sat = []
        plot_list = generations[i]
        for j in plot_list:

            calc.append(j.entropy)
            if len(list(j.node.ancestors)) > 1:
                calc_no_sat.append(j.entropy)

        average_entropy.append(np.sum(calc) / len(calc))
        average_entropy_no_sat.append(np.sum(calc_no_sat) / len(calc_no_sat))

    no_gens = []
    for i in generations:
        redshift = i[0].node['redshift']
        red_in = 0
        while redshift > initial_redshift:
            red_in = red_in + 1
            redshift = i[red_in].node['redshift']
        no_gens.append(i[red_in].node['redshift'])
    plt.plot(no_gens, average_entropy)
    plt.title("Average Entropy of each generation")
    plt.xlabel("Redshift z")
     
    plt.ylabel("Tree Entropy s")
    plt.savefig("/disk01/jhorlock/sample/av_ensample.png")

    plt.clf()

    plt.plot(no_gens, average_entropy_no_sat)
    plt.title("Average Entropy of central halos each generation")
    plt.xlabel("Redshift z")
    plt.xlim(initial_redshift, 0)
    plt.ylabel("Tree Entropy s")
    plt.savefig("/disk01/jhorlock/sample/av_en_censample.png")

    plt.clf()

    redshift = round(z_list[len(plot_x) - 1], 1)
    plt.plot(plot_y[len(plot_x) - 1], plot_x[len(plot_x) - 1] / 100 ** 3, label="z = {}".format(redshift))

    plt.title("Number density of halos with different tree entropy at z = 0 ")
    plt.xlabel("Tree Entropy s")
    plt.ylabel("Number density ($Mpc^{-3}$)")
    plt.xlim(0, 1.0)
    plt.savefig("/disk01/jhorlock/sample/final_fullsample.png")


    f = open("/disk01/jhorlock/sample/out.txt", "a")
    f.write(str(less*100 / (less + other)) + " full")
    f.write(str(no_res*100 / (no_res + res)) + " resolved in full")
    f.close()


main()

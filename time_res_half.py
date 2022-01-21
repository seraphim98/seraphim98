import ytree
import math as m
import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.stats import beta as beta_fitting
from scipy import stats
sys.path.insert(1, '/disk01/jhorlock')
from Halo import Halo

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

    snew = entropy * (np.sum(masses) / halo_mass)**gamma
    if snew > 1:
        pass
    else:
        entropy = snew

    halo.entropy_update(entropy)


def one_ancestor(halo, ancestor, gamma):

    halo_mass = halo.node['mass'].to_value()
    mass = ancestor.node['mass'].to_value()
    entropy = ancestor.entropy
    snew = entropy * (mass/halo_mass) ** gamma
    if snew > 1:
        pass
    else:
        entropy = snew

    halo.entropy_update(entropy)

#/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/rockstar/trees/

def main():
    treefile = "/disk12/legacy/GVD_C700_l100n2048_SLEGAC/dm_gadget/rockstar/trees/tree_0_0_0.dat"
    tree = ytree.load(treefile)
    fn = tree.save_arbor()
    tree = ytree.load(fn)
    halos_0 = tree.select_halos('tree["tree", "redshift"] == 0 ')
    halos = list(halos_0)

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

    
    half_generations = []
    half_generations_id = []
    calc_list = []
    num_gen = len(list(halos_0[0]['prog']))

    no_res = 0
    res = 0
    resolution = 10 ** 10
    indices = []

    for i in range(int(num_gen/2)):
        half_generations.append([])
        half_generations_id.append([])
    print(len(halos))
    for i in halos:
        prog_list = list(i['prog'])
        if len(prog_list) % 2 == 1:
            index = int((num_gen - len(prog_list) - 1) / 2)
            index2 = len(prog_list) - 1
            half_generations[index].append(Halo(prog_list[index2], 0, prog_list[index2]['mass'].to_value()))
            half_generations_id[index].append(prog_list[index2]['uid'])
        if len(prog_list) % 2 == 0:
            index = int((num_gen - len(prog_list)) / 2)
            index2 = len(prog_list) - 2
            first_halo = prog_list[index2]
            half_generations[index].append(Halo(first_halo, 0, first_halo['mass'].to_value()))
            half_generations_id[index].append(prog_list[index2]['uid'])
    print("Halo loop complete")
    for i in range(len(half_generations)):

        test_list = []
        calc_gen = half_generations[i]
        if i == 0:
            initial_redshift = calc_gen[0].node['redshift']
        if i != (len(half_generations) - 1):

            for j in range(len(calc_gen)):
                halo_calc = calc_gen[j].node.descendent
                halo_calc = halo_calc.descendent
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
                        calc_ancestors = []
                        for k in halo_calc_ancestors:
                            second_ancestors = list(k.ancestors)
                            for l in second_ancestors:
                                calc_ancestors.append(l)
                        halo_calc_ancestors = calc_ancestors

                        if len(halo_calc_ancestors) == 1:
                            k = halo_calc_ancestors[0]
                            x = half_generations_id[i].index(k['uid'])
                            initial_mass = calc_gen[x].initial_mass
                            halo_object = Halo(halo_calc, initial_entropy, initial_mass)
                            halo_object.multi_step_ancestors(1)
                            one_ancestor(halo_object, calc_gen[x], gamma)
                            half_generations[i + 1].append(halo_object)
                            half_generations_id[i + 1].append(halo_object.node['uid'])
                            less = less + 1
                        else:
                            other = other + 1

                            test_list.append(halo_calc['uid'])

                            for k in halo_calc_ancestors:
                                if k['uid'] not in half_generations_id[i]:
                                    calc_list.append(Halo(k, 0, k['mass'].to_value()))
                                else:
                                    x = half_generations_id[i].index(k['uid'])
                                    calc_list.append(calc_gen[x])
                            initial_mass = calc_list[0].initial_mass
                            halo_object = Halo(halo_calc, initial_entropy, initial_mass)
                            halo_object.multi_step_ancestors(len(halo_calc_ancestors))
                            x_list = mass_calc(calc_list)
                            H = h_calc(x_list, alpha, f)
                            entropy_calc(halo_object, x_list, H, calc_list, a, b, gamma)

                            half_generations[i+1].append(halo_object)
                            half_generations_id[i + 1].append(halo_object.node['uid'])


        calc_list = []
        print("Calculations...")
        print(str((i+1) * 100/len(half_generations)) + '%')

    half_final_gen = half_generations[int(num_gen/2 - 1)]
    f = open("/disk01/jhorlock/time_res/half_text.txt", 'w')

    for i in range(len(half_final_gen)):


        half = half_final_gen[i]
        if half.num_ancestors > 1 and half.node['mass'].to_value() > 10**11:
            f.write("{0:<9f} {1:<15f}\n".format(half.entropy, half.node['uid']))


    f.close()

    f = open("/disk01/jhorlock/time_res/half_av_en.txt", 'w')
    f_no_sat = open("/disk01/jhorlock/time_res/half_av_en_no_sat.txt", 'w')



    for i in half_generations :

        calc = []
        calc_no_sat = []
        plot_list = i
        for j in plot_list:

            calc.append(j.entropy)
            if j.num_ancestors > 1:
                calc_no_sat.append(j.entropy)
        redshift = i[0].node['redshift']
        red_in = 0
        while redshift > initial_redshift:
            red_in = red_in + 1
            redshift = i[red_in].node['redshift']

        f.write("{0:<9f} {1:<15f}\n".format(np.mean(calc), redshift))
        if len(calc_no_sat) == 0:
            pass
        else:
            f_no_sat.write("{0:<9f} {1:<15f}\n".format(np.mean(calc_no_sat), redshift))
    f.close()
    f_no_sat.close()

main()

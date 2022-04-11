import GH_Code
import F_Code
import I_Code
import misc_functions as MF

import subprocess
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os
from os.path import exists
from itertools import product

def make_series(directories, all):
    file_info = {
        # "Force": ["forceY.txt", ' '],
        "Energy": ["time_vs_energyInput.txt", '\t'],
        # "rotation": ["rotation.txt", '\t'],
        # "velocity": ["velocity.txt", '\t'],
        "kinetic": ["kinetic.txt", '\t'],
        # "corners": ["corners.txt", '\t']
    }
    energy_file_generation(directories["Local_Sims"])
    # work_file_generation(directories["Local_Sims"])

    # Gathering all data into a single dataframe
    property_array = []
    W_array = []
    A_array = []
    T_array = []
    data_array = []
    borderline_array = []
    # Uncomment next line if it
    # directories["Local_Sims"] = directories["Local_Sims"] + "wet_0.4/4.0/"
    borderlines = np.loadtxt((directories["Local_Sims Base"] + "Borderline_Data.csv"), delimiter=',',
                             skiprows=1)  # loading borderline cases

    for (dirpath, dirnames, filenames) in os.walk(directories["Local_Sims"]):
        for file in filenames:
            for key in file_info:
                if file == file_info[key][0]:  # tests if file is forceY.txt etc
                    dirpath = dirpath.replace('\\', '/')
                    path_list = dirpath.split('/')
                    ref_point = path_list.index("rigid_particle")
                    W = float(path_list[ref_point + 1][4:])
                    A = float(path_list[ref_point + 2])
                    T = float(path_list[ref_point + 3][path_list[ref_point + 3].index('T') + 1:])
                    ifborderline = Borderline_Test(borderlines, W, A, T)

                    if all == True or ifborderline:
                        file_generation_individual(dirpath, force=False)
                        print(str(W) + " " + str(A) + " " + str(T) + " " + key)
                        property_array.append(key), W_array.append(W), A_array.append(A), T_array.append(
                            T), borderline_array.append(ifborderline)
                        data = (np.loadtxt(os.path.join(dirpath, file), delimiter=file_info[key][1]))
                        data_array.append(data)

    index = pd.MultiIndex.from_arrays([borderline_array, W_array, A_array, T_array, property_array],
                                      names=['Borderline', 'W', 'A', 'T', 'Property'])
    ps = pd.Series(index=index, data=data_array)
    return ps

def overall_energy_well(location_overall, wet_list):
    data_wet = {}
    for wet in wet_list:
        data = np.loadtxt(location_overall + wet + "/energyWell.txt", delimiter=" ")
        E_Overall = data[-1, 1]
        data_wet[wet] = E_Overall
    return data_wet[wet]

def energy_file_generation(Local_Sims):
    for (dirpath, dirnames, filenames) in os.walk(Local_Sims):
        for file in filenames:
            if file == "forceY.txt":
                file1 = open(dirpath + '/forceY.txt', 'r')
                Line = file1.readline(1)
                if Line == '#':
                    print("Generating Energy txt file")
                    subprocess.call("sed -i '1d' ./forceY.txt", shell=True, cwd=dirpath)
                    subprocess.call("sed -i '1d' ./forceY.txt", shell=True, cwd=dirpath)
                    subprocess.call('g++ -o work_new workDone.cpp', shell=True, cwd=dirpath + "/postprocess")
                    subprocess.call('./work_new', shell=True, cwd=dirpath + "/postprocess")

def work_file_generation(Local_Sims):
    for (dirpath, dirnames, filenames) in os.walk(Local_Sims):
        for file in filenames:
            if file == "velocity.txt":
                path_list = dirpath.split('/')
                if path_list[-1] != "postprocess":
                    subprocess.call('rm 1_ke_vs_time.txt', shell=True, cwd=dirpath)
                    subprocess.call('rm rotation.cpp', shell=True, cwd=dirpath)
                    subprocess.call('rm rotation', shell=True, cwd=dirpath)
                    subprocess.call('rm rotation.txt', shell=True, cwd=dirpath)
                    subprocess.call('rm velocity.txt', shell=True, cwd=dirpath)
                    subprocess.call('rm t_meas.cpp', shell=True, cwd=dirpath)
                    subprocess.call('rm t_meas.cpp', shell=True, cwd=dirpath)
                    subprocess.call('rm t_meas_new', shell=True, cwd=dirpath)

def file_generation_individual(dirpath, force):
    path_list = dirpath.split('/')
    # if path_list[-1] != "postprocess":
    #     subprocess.call('rm 1_ke_vs_time.txt', shell=True, cwd=dirpath)
    #     subprocess.call('rm rotation.cpp', shell=True, cwd=dirpath)
    #     subprocess.call('rm rotation', shell=True, cwd=dirpath)
    #     subprocess.call('rm rotation.txt', shell=True, cwd=dirpath)
    #     subprocess.call('rm t_meas.cpp', shell=True, cwd=dirpath)
    #     subprocess.call('rm t_meas.cpp', shell=True, cwd=dirpath)
    #     subprocess.call('rm t_meas_new', shell=True, cwd=dirpath)

    if path_list[-1] == "postprocess":
        # checks if file KE exists. If not it copies it over and generates data
        if not exists(dirpath + '/1_ke_vs_time.txt'):
            print("Generating t_meas txt file")
            subprocess.call('cp /media/cameron/Data_Storage/Simulation_Files/Misc/t_meas.cpp ' + dirpath,
                            shell=True, cwd=dirpath)
            subprocess.call('g++ -o t_meas_new t_meas.cpp', shell=True, cwd=dirpath)
            subprocess.call('./t_meas_new', shell=True, cwd=dirpath)

        # checks if rotation.cpp file exists. If not it copies it over and generates data
        if not exists(dirpath + '/rotation.cpp') or force:
            print("Generating rotation and velocity txt files")
            print(dirpath)
            subprocess.call('cp /media/cameron/Data_Storage/Simulation_Files/Misc/rotation.cpp ' + dirpath,
                            shell=True, cwd=dirpath)
            subprocess.call('g++ -o rotation rotation.cpp', shell=True, cwd=dirpath)
            subprocess.call('./rotation', shell=True, cwd=dirpath)

        # checks for issues in /1_ke_vs_time.txt file
        file1 = open(dirpath + '/1_ke_vs_time.txt', 'r')
        Line = file1.readline(30)

        if Line[-1] == '\t':
            print(Line)
            print("copying")
            subprocess.call('cp /media/cameron/Data_Storage/Simulation_Files/Misc/t_meas.cpp ' + dirpath,
                            shell=True, cwd=dirpath)
            print("Generating t_meas txt file")
            subprocess.call('g++ -o t_meas_new t_meas.cpp', shell=True, cwd=dirpath)
            subprocess.call('./t_meas_new', shell=True, cwd=dirpath)


def theta_plot():
    # Angle plots for different particle types
    subprocess.call('rm -r Theta_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    subprocess.call('mkdir Theta_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")

    for W in numerical_wet_list:
        sliced_ps = ps.xs(("rotation", W, True), level=("Property", "W", "Borderline"))
        Amplitude = []
        Average_delta_theta = []

        for index, value in sliced_ps.items():

            for j in range(len(value[0, :])):
                for i in range(len(value[:, 0])-2):

                    # shifting up if going negitive
                    if value[i, j]-value[i+1, j] > 70:
                        value[i+1:, j] = value[i+1:, j] + 94.08

                    # shifting down if going positive
                    if value[i, j]-value[i+1, j] < -70:
                        value[i+1:, j] = value[i+1:, j] - 94.08

            delta_x = value[:, 1] - value[:, 2]
            delta_y = value[:, 3] - value[:, 4]
            delta_z = value[:, 5] - value[:, 6]

            theta_xy = ((np.arctan2(delta_y, delta_x) * 180 / np.pi))
            theta_zy = ((np.arctan2(delta_y, delta_z) * 180 / np.pi))

            # theta_xy[theta_xy < 0] += 360
            # theta_zy[theta_zy < 0] += 360

            for i in range(len(theta_xy)-2):

                # shifting up if going negitive
                if theta_xy[i]-theta_xy[i+1] > 100:
                    theta_xy[i+1:] = theta_xy[i+1:] + 360

                if theta_zy[i]-theta_zy[i+1] > 100:
                    theta_zy[i+1:] = theta_zy[i+1:] + 360

                # shifting down if going positive
                if theta_xy[i]-theta_xy[i+1] < -100:
                    theta_xy[i+1:] = theta_xy[i+1:] - 360

                if theta_zy[i]-theta_zy[i+1] < -100:
                    theta_zy[i+1:] = theta_zy[i+1:] - 360


            Amplitude.append(index[0])
            Average_delta_theta.append(np.average(np.abs(np.diff(theta_zy))) + np.average(np.abs(np.diff(theta_xy))))

            name = str(index[0])
            plt.plot(value[:len(theta_xy), 0], theta_xy[:len(theta_xy)], linewidth=0.45, label=name + " xy")
            plt.plot(value[:len(theta_xy), 0], theta_zy[:len(theta_xy)], linewidth=0.45, label=name + " zy",
                     linestyle=(0, (5, 1)))

            # plt.plot(value[:,0], delta_x[:], linewidth=0.45, label=name + " x")
            # plt.plot(value[:, 0], delta_y[:], linewidth=0.45, label=name + " y")
            # plt.plot(value[:, 0], delta_z[:], linewidth=0.45, label=name + " z")

        leg = plt.legend(loc='upper left', bbox_to_anchor=(1, 0.9), title = "Amplitude (\u212B)")
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)

        plt.tight_layout(pad=3)
        plt.title("Total Angle Swept vs Time")
        plt.xlabel("Time (fs)")
        plt.xlim(left=0)
        plt.ylabel("Total Angle Swept (degrees)")
        plt.annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction")
        plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Theta_Graphs/Theta " + str(W) + ".png",
                    dpi=600), plt.close()

        plt.xlabel("Amplitude (\u212B)")
        plt.ylabel("Average rate of rotation (degrees/fs)")
        plt.title("Average rate of Rotation vs Amplitude")
        plt.scatter(Amplitude, np.array(Average_delta_theta)/100)
        plt.ylim([0, 1.2 * max(Average_delta_theta)])
        plt.annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction")
        plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Theta_Graphs/Average Theta " + str(W) + ".png",
                    dpi=600), plt.close()

    # combined average angle plot:
    for W in numerical_wet_list:
        sliced_ps = ps.xs(("rotation", W, True), level=("Property", "W", "Borderline"))
        Amplitude = []
        Average_delta_theta = []

        for index, value in sliced_ps.items():

            for j in range(len(value[0, :])):
                for i in range(len(value[:, 0])-2):

                    # shifting up if going negitive
                    if value[i, j]-value[i+1, j] > 70:
                        value[i+1:, j] = value[i+1:, j] + 94.08

                    # shifting down if going positive
                    if value[i, j]-value[i+1, j] < -70:
                        value[i+1:, j] = value[i+1:, j] - 94.08

            delta_x = value[:, 1] - value[:, 2]
            delta_y = value[:, 3] - value[:, 4]
            delta_z = value[:, 5] - value[:, 6]

            theta_xy = ((np.arctan2(delta_y, delta_x) * 180 / np.pi))
            theta_zy = ((np.arctan2(delta_y, delta_z) * 180 / np.pi))

            # theta_xy[theta_xy < 0] += 360
            # theta_zy[theta_zy < 0] += 360

            for i in range(len(theta_xy)-2):

                # shifting up if going negitive
                if theta_xy[i]-theta_xy[i+1] > 100:
                    theta_xy[i+1:] = theta_xy[i+1:] + 360

                if theta_zy[i]-theta_zy[i+1] > 100:
                    theta_zy[i+1:] = theta_zy[i+1:] + 360

                # shifting down if going positive
                if theta_xy[i]-theta_xy[i+1] < -100:
                    theta_xy[i+1:] = theta_xy[i+1:] - 360

                if theta_zy[i]-theta_zy[i+1] < -100:
                    theta_zy[i+1:] = theta_zy[i+1:] - 360


            Amplitude.append(index[0])
            Average_delta_theta.append(np.average(np.abs(np.diff(theta_zy))) + np.average(np.abs(np.diff(theta_xy))))

        plt.title("Average rate of Rotation vs Amplitude")
        plt.xlabel("Amplitude (\u212B)")
        plt.ylabel("Average rate of rotation (degrees/fs)")
        plt.plot(Amplitude, np.array(Average_delta_theta)/100)
        plt.scatter(Amplitude, np.array(Average_delta_theta)/100, label=str(W))

        leg = plt.legend(loc = 'best', title = "\u03B5 (kcal/mol)")
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)

    plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Theta_Graphs/Average Theta 2.png",
                dpi=600)
    plt.close()
    # plt.clf()

def corners_plot():
    # Angle plots for different particle types
    # subprocess.call('rm -r min_dist_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    # subprocess.call('mkdir min_dist_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")

    plot_range = [0.0, 0.2]

    for W in numerical_wet_list:
        corner_ps = ps.xs(("corners", W, True), level=("Property", "W", "Borderline"))
        force_ps = ps.xs(("Force", W, True), level=("Property", "W", "Borderline"))
        # subprocess.call('mkdir ' + str(W), shell=True,
        #                 cwd=directories["Local_Sims Base"] + "Overview_Graphs/min_dist_Graphs/")

        plt.clf()

        for index, force_value in force_ps.items():
            # Corner Data
            corners_value = corner_ps[index[0], index[1]]
            x_array = np.array(corners_value[:,0])
            min_array = np.array(corners_value[:, 1:8]).min(axis=1)
            y_array_min = min_array - corners_value[:, 9]
            y_array_avg = np.average(corners_value[:, 1:8], axis=1) - corners_value[:, 9]

            # force data
            fig, ax = plt.subplots()
            arr_len = len(force_value[:, 0])
            i_min = int(arr_len * plot_range[0])
            i_max = int(arr_len * plot_range[1])
            lns1 = ax.plot(force_value[i_min:i_max, 0] - force_value[0, 0], force_value[i_min:i_max, 1], linewidth=1, color="green", label="$F_{particle}$")

            title_name = "(\u03B5=" + str(W) + " kcal/mol, A = " + str(index[0]) + "\u212B)"
            plt.annotate(title_name, xy=(0.05, 0.9), xycoords="axes fraction", bbox = dict(boxstyle="round,pad=0.3", fc="w"))
            plt.title("Force and Particle Seperation Over Time")
            ax.set_xlabel("Time (fs)")
            ax.set_ylabel("Force (kcal/mol-\u212B)")
            ax.set_xlim(left = 0)

            ax2 = ax.twinx()
            arr_len = len(x_array)
            i_min = int(arr_len*plot_range[0])
            i_max = int(arr_len*plot_range[1])
            lns2 = ax2.plot(x_array[i_min:i_max] / 2 - x_array[0], y_array_min[i_min:i_max], linewidth=1, label="$d_{edge}$")
            lns3 = ax2.plot(x_array[i_min:i_max] / 2 - x_array[0], y_array_avg[i_min:i_max], linewidth=1, label="$d_{centre}$")

            # added these three lines
            lns = lns1 + lns2 + lns3
            labs = [l.get_label() for l in lns]
            leg = ax2.legend(lns, labs, loc='upper right')
            for legobj in leg.legendHandles:
                legobj.set_linewidth(3.0)

            ax2.set_ylabel("Seperation, d (\u212B)")

            fig.savefig(directories["Local_Sims Base"] + "Overview_Graphs/min_dist_Graphs/" + str(W) + "/min_dist " + str(index[0]) + ".png",
                        dpi=600)
            plt.clf()
            plt.close()

def energy_plot():
    # Energy plot for each wettability:
    subprocess.call('rm -r Energy_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    subprocess.call('mkdir Energy_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    for W in numerical_wet_list:
        # plot of Energy diagrams for bordeline cases
        sliced_ps = ps.xs(("Energy", W), level=("Property", "W"))
        plot_range = [0.0, 0.2]
        for index, value in sliced_ps.items():
            length = len(value[:, 0])
            low = int(length * plot_range[0])
            high = int(length * plot_range[1])

            x_array = value[low:high, 0]
            plt.plot(value[low:high, 0], value[low:high, 1], linewidth=0.45,
                         label=str(index[1]))

        plt.title("Work Done by System vs Time")

        leg = plt.legend(loc='lower right', title = "Amplitude (\u212B)")
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)

        plt.xlabel("Time (fs)")
        plt.ylabel("Work Done by System (kcal/mol)")
        plt.xlim(left=0)
        plt.plot(x_array, np.ones(len(x_array)) * energy_well[W], linewidth=0.5, linestyle='dashed')
        plt.annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction")
        plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Energy_Graphs/" + str(W) + " Energy.png",
                    dpi=600), plt.close()

def force_and_energy_plot(numerical_wet_list):
    # Energy plot for each wettability:
    # subprocess.call('rm -r Force_and_Energy_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    # subprocess.call('mkdir Force_and_Energy_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")

    # Combined Plot of normalized velocities

    for W in numerical_wet_list:
        fig, ax = plt.subplots(1, 2, facecolor='w', figsize=(10, 5))
        ax = ax.ravel()
        # plot of Energy diagrams for bordeline cases
        energy_ps = ps.xs(("Energy", W), level=("Property", "W"))
        plot_range = [0.0, 0.2]
        for index, value in energy_ps.items():
            length = len(value[:, 0])
            low = int(length * plot_range[0])
            high = int(length * plot_range[1])

            x_array = value[low:high, 0]
            ax[1].plot(value[low:high, 0]*5, value[low:high, 1], linewidth=1.0,
                         label=str(index[1]))

        plt.title("Work Done by System vs Time")

        leg = ax[1].legend(loc='lower right', title = "Amplitude (\u212B)")
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)

        ax[1].set_xlabel("Time (fs)")
        ax[1].set_ylabel("Work Done by System (kcal/mol)")
        ax[1].set_title("Work Done by System vs Time")
        ax[1].set_xlim(left=0)
        ax[1].plot(x_array*5, np.ones(len(x_array)) * energy_well[W], linewidth=1.0, linestyle='dashed')
        ax[1].annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction", fontsize=10,
                     bbox=dict(boxstyle="round,pad=0.3", fc="w"))

        # plot of Energy diagrams for bordeline cases
        sliced_ps = ps.xs(("Force", W), level=("Property", "W"))
        plot_range = [0.0, 0.2]

        for index, value in sliced_ps.items():
            length = len(value[:, 0])
            low = int(length * plot_range[0])
            high = int(length * plot_range[1])

            x_array = value[low:high, 0] - value[low, 0]
            ax[0].plot(x_array*10, value[low:high, 1], linewidth=1.0,
                     label=str(index[1]))

        ax[0].set_title("Force Acting on Particle vs Time")
        # leg = ax[0].legend(loc='upper right', title="Amplitude (\u212B)")
        # for legobj in leg.legendHandles:
        #     legobj.set_linewidth(3.0)
        ax[0].set_xlabel("Time (fs)")
        ax[0].set_xlim(left=0)
        ax[0].set_ylabel("Force Acting on Particle (kcal/mol-\u212B)")
        ax[0].plot(x_array*10, np.zeros(len(x_array)), linewidth=1.0, linestyle='dashed')
        ax[0].annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction", fontsize=10,
                     bbox=dict(boxstyle="round,pad=0.3", fc="w"))
        leg = ax[0].legend(loc='upper right', title="Amplitude (\u212B)")
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)

        ax[0].xaxis.set_major_locator(plt.MaxNLocator(5))
        ax[1].xaxis.set_major_locator(plt.MaxNLocator(5))

        fig.tight_layout()
        fig.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Force_and_Energy_Graphs/" + str(W) + "Force_and_Energy.png",
                    dpi=600), plt.clf()

def energy_well_plot(numerical_wet_list):
    # Energy plot for each wettability:
    # subprocess.call('rm -r Energy_Well_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    # subprocess.call('mkdir Energy_Well_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")

    # for W in numerical_wet_list:
    #     # plot of Energy diagrams for bordeline cases
    #     path = directories["Local_Test"] + "wet_" + str(W) + "/energyWell.txt"
    #     data = (np.loadtxt(path, delimiter=' ', skiprows=2))
    #
    #     plt.plot(data[:,0], data[:,1], linewidth=1.0, label=str(W))
    #
    #     plt.title("Work Done by System vs Time (Wettabilty = " + str(W) + "(kcal/mol)")
    #     plt.xlabel("Time (fs)")
    #     plt.ylabel("Energy Well of Particle (kcal/mol)")
    #     plt.xlim(left=0)
    #     plt.annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction", fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="w"))
    #     plt.savefig(
    #         directories["Local_Sims Base"] + "Overview_Graphs/Energy_Well_Graphs/" + str(W) + "_Energy_Well.png",
    #         dpi=600), plt.close()

    for W in numerical_wet_list:
        # plot of Energy diagrams for bordeline cases
        path = directories["Local_Test"] + "wet_" + str(W) + "/energyWell.txt"
        data = (np.loadtxt(path, delimiter=' ', skiprows=2))

        plt.plot(data[:, 0], data[:, 1], linewidth=1.0, label=str(W))

    plt.title("Adhesion Energy of Particle vs Time")
    leg = plt.legend(loc='lower right', title="\u03B5 (kcal/mol)", ncol=3)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)

    plt.xlabel("Time (fs)")
    plt.ylabel("Adhesion Energy of Particle (kcal/mol)")
    plt.ylim([0, -375])
    plt.xlim(left=0)
    plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
    # plt.text(0, 0, "Test", ha="center", va="center", rotation=0, size=15,
    #     bbox=dict(boxstyle="round,pad=0.3", lw=2))
    plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Energy_Well_Graphs/Combined_Energy_Well.png",
                    dpi=600), plt.close()

def border_line_liftoff_plot(numerical_wet_list):
    # Energy plot for each wettability:
    # subprocess.call('rm -r Amp_vs_TP_Plot', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    # subprocess.call('mkdir Amp_vs_TP_Plot', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")

    borderlines = np.loadtxt((directories["Local_Sims Base"] + "Borderline_Data.csv"), delimiter=',',
                             skiprows=1)  # loading borderline cases

    # Single Borderline Plots (Amplitude vs Time Period)
    borderlines_split = np.split(borderlines, 6, axis=0)
    for data in borderlines_split:
        W = data[0, 0]
        plt.plot(data[:, 1], data[:, 2], linewidth=1.0, label=str(W), color = wet_colour[W])
        plt.scatter(data[:, 1], data[:, 2], linewidth=1.0, color = wet_colour[W])

    plt.title("Borderline Time Period vs Amplitude")
    plt.xlabel("Amplitude (\u212B)")
    plt.ylabel("Borderline Time Period (fs)")
    leg = plt.legend(loc='upper left', title="\u03B5 (kcal/mol)")
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)

    plt.ylim(bottom=15000)
    # plt.xlim(left=4)
    plt.savefig(
        directories["Local_Sims Base"] + "Overview_Graphs/Amp_vs_TP_Plot/TP_vs_A.png",
        dpi=600), plt.close()



    # Combined Plot
    # Time period vs amplitude
    borderlines_split = np.split(borderlines, 6, axis=0)
    fig, ax = plt.subplots(1,3, figsize=(15,5))
    ax = ax.ravel()

    for data in borderlines_split:
        W = data[0, 0]

        ax[0].plot(data[:, 1], data[:, 2], linewidth=1.0, label=str(W), color = wet_colour[W])
        ax[0].scatter(data[:, 1], data[:, 2], linewidth=1.0, color = wet_colour[W])

    ax[0].set_title("$\it{TP} \; vs \; Amplitude$")
    ax[0].set_xlabel("Amplitude (\u212B)")
    ax[0].set_ylabel("$\it{TP} (fs)$")
    leg = ax[0].legend(loc='upper left', title="\u03B5 (kcal/mol)", ncol=2)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)

    # ax[0].set_ylim(bottom=0)
    # ax[0].set_xlim(right=11)
    # ax[0].annotate("\u0394TP", xy=(10.2, 57000), fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="w"))
    # ax[0].annotate(text='', xy=(10,62437), xytext=(10,52040), arrowprops=dict(arrowstyle='<->'))

    # sorted by amplitude
    #borderlines_split = np.split(borderlines[np.argsort(borderlines[:, 1])], 6, axis=0)
    borderlines_split = np.split(borderlines[np.argsort(borderlines[:, 1])], 7, axis=0)
    print(borderlines_split[0])
    for data in borderlines_split:
        data = data[np.argsort(data[:, 0])]
        A = data[0, 1]
        ax[1].plot(data[:, 0], data[:, 2], linewidth=1.0, label=str(A))
        ax[1].scatter(data[:, 0], data[:, 2], linewidth=1.0)

    ax[1].set_title("$\it{TP} \; vs \; \epsilon$")
    leg = ax[1].legend(loc = 'upper right', title="Amplitude (\u212B)", ncol = 2)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)
    ax[1].set_xlabel("\u03B5 (kcal/mol)")
    ax[1].set_ylabel("$\it{TP} \; (fs)$")

    borderlines_split = np.split(borderlines[np.argsort(borderlines[:, 1])], 7, axis=0)
    for data in borderlines_split:
        data = data[np.argsort(data[:, 0])]
        A = data[0, 1]
        ax[2].plot(data[:-1, 0], np.diff(data[:, 2]), linewidth=1.0, label=str(A))
        ax[2].scatter(data[:-1, 0], np.diff(data[:, 2]), linewidth=1.0)

    ax[2].set_title(r"$\nabla \it{TP} \; vs \; \epsilon$")
    leg = ax[2].legend(loc = 'lower right', title="Amplitude (\u212B)", ncol = 1)
    for legobj in leg.legendHandles:
        legobj.set_linewidth(3.0)
    ax[2].set_xlabel("\u03B5 (kcal/mol)")
    ax[2].set_ylabel(r"$\nabla \it{TP} (fs \, kcal/mol)$")

    fig.suptitle("Borderline Lift-off Plots")
    fig.tight_layout()
    fig.savefig(
        directories["Local_Sims Base"] + "Overview_Graphs/Amp_vs_TP_Plot/Combined_Plot.png",
        dpi=600), plt.clf()

def force_plot():
    # Force plot for each wettability:
    # subprocess.call('rm -r Force_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    # subprocess.call('mkdir Force_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    for W in numerical_wet_list:
        # plot of Energy diagrams for bordeline cases
        sliced_ps = ps.xs(("Force", W), level=("Property", "W"))
        plot_range = [0.0, 0.2]

        for index, value in sliced_ps.items():
            length = len(value[:, 0])
            low = int(length * plot_range[0])
            high = int(length * plot_range[1])

            x_array = value[low:high, 0]-value[low, 0]
            plt.plot(x_array, value[low:high, 1], linewidth=1.0,
                     label=str(index[1]))

        plt.title("Force Acting on Particle vs Time")
        leg = plt.legend(loc='upper right', title = "Amplitude (\u212B)")
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)
        plt.xlabel("Time (fs)")
        plt.xlim(left=0)
        plt.ylabel("Force Acting on Particle (kcal/mole-Angstrom)")
        plt.plot(x_array, np.zeros(len(x_array)) * energy_well[W], linewidth=0.5, linestyle='dashed')
        plt.annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction", fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="w"))
        plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Force_Graphs/" + str(W) + " Force.png",
                    dpi=600), plt.close()

def energy_scatter_plot(energy_well_array):
    # Energy plot of expected energy vs actual energy for each wettability:
    # subprocess.call('rm -r L_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    # subprocess.call('mkdir L_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")

    # for W in numerical_wet_list:
    #     energy_well = energy_well_array[W]
    #     # plot of Energy diagrams for bordeline cases
    #     energy_ps = ps.xs(("Energy", W, True), level=("Property", "W", "Borderline"))
    #     ke_ps = ps.xs(("kinetic", W, True), level=("Property", "W", "Borderline"))
    #     x_array, D_array, Att_array, theory_array, index_array = [], [], [], [], []
    #     ke_array = []
    #
    #     Total_work_done = []
    #     for index, value in energy_ps.items():
    #         x_array.append(index[0])
    #         D_array.append(np.min(value[:,1]))
    #         index_array.append([np.argmin(value[:,1]),index[0],index[1]])
    #         Att_array.append(np.min(value[:, 1]) - float(value[-1:, 1]))
    #         Total_work_done.append(float(value[-1:, 1]))
    #
    #
    #         Amplitude, Time_period = index[0], index[1]
    #         Amplitude = Amplitude * 10 ** (-10)
    #         Frequency = 1 / (Time_period * 10 ** (-15))
    #         mass = (195.084 * (10973 - 10369))*1.6605*10**(-27)
    #         KE = -(0.5 * mass * (2 * Amplitude * 2 * np.pi * Frequency) ** 2) * 1/4184 * 6.0221 * 10**23
    #
    #         theory_array.append(KE)
    #
    #     for index, value in ke_ps.items():
    #         ke_array.append(value[-1, 1])
    #
    #     ke_array = np.array(ke_array)
    #
    #     x_array = np.array(x_array)
    #     D_array = np.array(D_array)
    #     Att_array = np.array(Att_array)
    #     theory_array = np.array(theory_array)
    #
    #     print("KE array is:", ke_array)
    #
    #     y_array_KE = energy_well / (theory_array - Att_array)
    #     y_Actual = (energy_well)/ (D_array - Att_array)
    #     y_Actual_ke = (energy_well - ke_array) / (D_array - Att_array)
    #
    #     plt.scatter(x_array, y_array_KE, label="Predicted")
    #     plt.scatter(x_array, y_Actual, label="Actual - Work Done")
    #     plt.scatter(x_array, y_Actual_ke, label="Actual - Work Done + Final Ke")
    #
    #     # L* Scatter
    #     plt.plot(x_array, np.ones(len(x_array)), linewidth=0.5, linestyle='dashed')
    #     plt.title("L* vs Amplitude")
    #     leg = plt.legend(loc='best')
    #     for legobj in leg.legendHandles:
    #         legobj.set_linewidth(3.0)
    #     plt.ylim([0,max(1.25,np.min(y_array_KE))])
    #     plt.xlabel("Amplitude (\u212B)"), plt.ylabel("L*")
    #     plt.annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction", fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="w"))
    #     plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/L_Graphs/L_Star plot " + str(W) + ".png", dpi=600), plt.close()
    #
    #     # Energy-D vs KE Scatter
    #     plt.scatter(x_array, theory_array, label="Theroretical")
    #     plt.scatter(x_array, D_array, label="Actual")
    #     plt.title("Work Done vs Amplitude")
    #     leg = plt.legend(loc='lower right')
    #     for legobj in leg.legendHandles:
    #         legobj.set_linewidth(3.0)
    #     plt.xlabel("Amplitude (\u212B)"), plt.ylabel("Work Done at Point D (KCal/Mol)")
    #     plt.ylim([0, 1.2*min(np.min(theory_array), np.min(D_array))])
    #     plt.annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction", fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="w"))
    #     plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/L_Graphs/Energy-D vs KE " + str(W) + ".png", dpi=600), plt.close()
    #
    #     # KE over time plot
    #     for index, value in ke_ps.items():
    #         plt.plot(value[:,0], value[:,1], linewidth=0.5, label=str(index[0]))
    #         plt.title("Kinetic Energy of Particle vs Time")
    #         plt.xlabel("Time (fs)"), plt.ylabel("Kinetic Energy of Particle (KCal/Mol)")
    #     leg = plt.legend(loc='upper right', title = "Amplitude (\u212B)")
    #     for legobj in leg.legendHandles:
    #         legobj.set_linewidth(3.0)
    #     plt.annotate("\u03B5="+str(W) + "(kcal/mol)", xy=(0.1, 0.9), xycoords="axes fraction", fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="w"))
    #     plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/L_Graphs/KE_Plot " + str(W) + ".png", dpi=600), plt.close()

    fig_KE, ax_KE = plt.subplots(2, 4, facecolor='w', figsize=(10, 5))
    ax_KE = ax_KE.ravel()

    fig_L, ax_L = plt.subplots(2, 3, facecolor='w', figsize=(10, 5))
    ax_L = ax_L.ravel()

    fig_FKE, ax_FKE = plt.subplots(2, 3, facecolor='w', figsize=(10, 5))
    ax_FKE = ax_FKE.ravel()

    # fig_KET, ax_KET = plt.subplots(2, 4, facecolor='w', figsize=(10, 5))
    # ax_KET = ax_KET.ravel()
    average_difference = np.zeros([6,2])

    i = 0
    for W in numerical_wet_list:

        order_array = [0, 1, 2, 4, 5, 6]
        energy_well = energy_well_array[W]
        # plot of Energy diagrams for bordeline cases
        energy_ps = ps.xs(("Energy", W, True), level=("Property", "W", "Borderline"))
        ke_ps = ps.xs(("kinetic", W, True), level=("Property", "W", "Borderline"))
        x_array, D_array, Att_array, theory_array, index_array = [], [], [], [], []
        ke_array = []

        Total_work_done = []
        for index, value in energy_ps.items():
            x_array.append(index[0])
            D_array.append(np.min(value[:,1]))
            index_array.append([np.argmin(value[:,1]),index[0],index[1]])
            Att_array.append(np.min(value[:, 1]) - float(value[-1:, 1]))
            Total_work_done.append(float(value[-1:, 1]))


            Amplitude, Time_period = index[0], index[1]
            Amplitude = Amplitude * 10 ** (-10)
            Frequency = 1 / (Time_period * 10 ** (-15))
            mass = (195.084 * (10973 - 10369))*1.6605*10**(-27)
            KE = -(0.5 * mass * (2 * Amplitude * 2 * np.pi * Frequency) ** 2) * 1/4184 * 6.0221 * 10**23

            theory_array.append(KE)

        for index, value in ke_ps.items():
            ke_array.append(value[-1, 1])


        ke_array = np.array(ke_array)
        x_array = np.array(x_array)
        D_array = np.array(D_array)
        Att_array = np.array(Att_array)
        theory_array = np.array(theory_array)
        y_array_KE = energy_well / (theory_array - Att_array)
        y_Actual = (energy_well)/ (D_array - Att_array)
        y_Actual_ke = (energy_well - ke_array) / (D_array - Att_array)

        # FKE Scatter
        if i < 6:
            ax_FKE[i].scatter(x_array, ke_array)
            ax_FKE[i].set_title("\u03B5="+str(W)+" kcal/mol")
            ax_FKE[i].set_ylim([0,60])

            if i == 3 or i == 0:
                ax_FKE[i].set_ylabel(r"$KE_{Final} (kcal/mol)$")
            else:
                ax_FKE[i].set_yticklabels([])

            if i == 3 or i == 4 or i == 5:
                ax_FKE[i].set_xlabel("Amplitude (\u212B)")
            else:
                ax_FKE[i].set_xticklabels([])

        # L* Scatter
        if i < 6:
            ax_L[i].scatter(x_array, y_array_KE, label="Predicted")
            ax_L[i].scatter(x_array, y_Actual, label="Actual - Work Done")
            ax_L[i].scatter(x_array, y_Actual_ke, label="Actual - Work Done + Final Ke")
            ax_L[i].plot(x_array, np.ones(len(x_array)), linewidth=0.5, linestyle='dashed')
            ax_L[i].set_title("\u03B5="+str(W)+" kcal/mol")
            ax_L[i].set_ylim([0,max(1.25,np.min(y_array_KE))])

            if i == 3 or i == 0:
                ax_L[i].set_ylabel(r"$L^*$")
            else:
                ax_L[i].set_yticklabels([])

            if i == 3 or i == 4 or i == 5:
                ax_L[i].set_xlabel("Amplitude (\u212B)")
            else:
                ax_L[i].set_xticklabels([])

        average_difference[i, :] = [W, np.mean(np.abs(D_array - theory_array))]

        # Energy-D vs KE Scatter
        ax_KE[order_array[i]].scatter(x_array, theory_array, label="Theroretical")
        ax_KE[order_array[i]].scatter(x_array, D_array, label="Actual")
        ax_KE[order_array[i]].set_title("\u03B5="+str(W) + " kcal/mol")

        if i == 3 or i == 0:
            ax_KE[order_array[i]].set_ylabel(r"$W_D (kcal/mol)$")
        else:
            ax_KE[order_array[i]].set_yticklabels([])

        if i == 3 or i == 4 or i == 5:
            ax_KE[order_array[i]].set_xlabel("Amplitude (\u212B)")
        else:
            ax_KE[order_array[i]].set_xticklabels([])

        ax_KE[order_array[i]].set_ylim([0, -1200])

        average_difference[i, :] = [W, np.mean(np.abs(D_array - theory_array))]
        i += 1

    # leg = ax_L[4].legend(bbox_to_anchor=(3, 0.5), title="Legend", ncol = 2)
    # for legobj in leg.legendHandles:
    #     legobj.set_linewidth(3.0)

    fig_L.suptitle("$L^* \; vs \; Amplitude$")
    fig_L.savefig(directories["Local_Sims Base"] + "Overview_Graphs/L_Graphs/L_Star plot Combined.png",
                  dpi=600)

    fig_FKE.suptitle("$KE_{Final} \; vs \; Amplitude$")
    fig_FKE.savefig(directories["Local_Sims Base"] + "Overview_Graphs/L_Graphs/FKE_Scatter.png",
                  dpi=600)


    # Energy-D vs KE Scatter
    ax_KE[7].scatter(average_difference[:,0], average_difference[:,1], color = 'g')
    ax_KE[7].set_title("$ \overline{\Delta W_D} \; vs \; \epsilon$")
    ax_KE[7].set_ylim([0, 1.2*max(average_difference[:,1])])
    ax_KE[7].set_xlabel("$\epsilon\; (kcal/mol)$"), ax_KE[7].set_ylabel(r"$\overline{\Delta W_D} (kcal/mol)$")

    # moving to the right
    box = ax_KE[7].get_position()
    box.x0 = box.x0 + 0.05
    box.x1 = box.x1 + 0.05
    ax_KE[7].set_position(box)

    fig_KE.suptitle("Work Done at Point D vs Amplitude")
    fig_KE.delaxes(ax_KE[3])
    # fig_KE.tight_layout()
    fig_KE.savefig(
        directories["Local_Sims Base"] + "Overview_Graphs/L_Graphs/Energy-D vs KE Combined.png",
        dpi=600), plt.clf()

def F_adhesion_gravity():
    # Energy plot for each wettability:
    mass = np.logspace(4,15, 1000)
    F_min_array = []

    for W in numerical_wet_list:
        sliced_ps = ps.xs(("Force", W), level=("Property", "W"))
        F_min = 0
        for index, value in sliced_ps.items():
            F_min = min(np.min(value), F_min)

        F_min_array.append([F_min, W])

    F_min_array = np.array(F_min_array)
    F_grav = mass * 9.81 * 1 / (1.6605 * 10 ** (-27)) * 1 / (1 / 4184 * 6.0221 * 10 ** (23)) * 10 ** (-10)

    for F_min, W in F_min_array:
        p = (604*194.084)/8000
        F_adh = -(F_min/20*20) * (mass/p)**(2/3)

        # x_array = np.where(F_adh == your_array_here)[0][0]
        # print(x_array, closest[0], F_adh[x_array])
        # plt.plot(y_array, y_array, linewidth=0.4, label="F_adh - " + str(W))
        plt.plot(mass, F_adh, linewidth=0.7, label = "F_adh - \u03B5=" + str(W))

    plt.plot(mass, F_grav, linewidth=0.7, label="F_Gravity")
    # for W in numerical_wet_list:
    #     plt.xlim(left=0)
    #
    plt.legend(loc = 'best')
    plt.xlabel("Mass (grams/mole)")
    plt.xscale('log')
    plt.xlim(left = 10**4)

    plt.ylabel("Force (kcal/mol-\u212B)")
    plt.yscale('log')
    plt.ylim(bottom=1)
    plt.title("Force vs Mass")

    plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Base_Forces_Mass.png",
                dpi=600), plt.close()

    # length scale:
    Length = (mass/p)**(2/3)

    F_min_array = np.array(F_min_array)
    F_grav = mass * 9.81 * 1 / (1.6605 * 10 ** (-27)) * 1 / (1 / 4184 * 6.0221 * 10 ** (23)) * 10 ** (-10)

    for F_min, W in F_min_array:
        p = (604*194.084)/8000
        F_adh = -(F_min/20*20) * (mass/p)**(2/3)
        plt.plot(Length, F_adh, linewidth=0.7, label = "F_adh - " + str(W))

    plt.plot(Length, F_grav, linewidth=0.7, label="F_Gravity")
    plt.legend(loc = 'best')
    plt.xlabel("Characterisitc Length (\u212B)")
    plt.xscale('log')
    plt.xlim(left = 100)

    plt.ylabel("Force (kcal/mol-\u212B)")
    plt.yscale('log')
    plt.ylim(bottom=1)
    plt.title("Force vs Length")

    plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Base_Forces_Length.png",
                dpi=600), plt.close()

def normalized_velocity_plot(numerical_wet_list):
    # plot of Kinetic Energy diagrams for bordeline cases
    # subprocess.call('rm -r Velocity_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")
    # subprocess.call('mkdir Velocity_Graphs', shell=True, cwd=directories["Local_Sims Base"] + "Overview_Graphs/")

    # plot for each wettability:
    for W in numerical_wet_list:
        # plot of Energy diagrams for bordeline cases
        sliced_ps = ps.xs(("velocity", True, W), level=("Property", "Borderline", "W"))

        for index, value in sliced_ps.items():
            # key parameters
            A, T = float(index[0]), float(index[1])


            factor = (A/T * 2 * np.pi)
            y_data = value[1:, 1]/factor

            y_data_change = y_data[:-1] - y_data[1::1]
            liftoff_pos = np.argmax(y_data_change > 0)
            x_data = value[1:, 0] / (value[liftoff_pos, 0])

            # print(value[liftoff_pos, 0])

            plt.plot(x_data, y_data, linewidth=1.0, label=str(A))

        plt.xlim([0, 6])
        leg = plt.legend(loc='upper right', title = "Amplitude (\u212B)")
        for legobj in leg.legendHandles:
            legobj.set_linewidth(3.0)
        plt.xlabel("T*")
        plt.ylabel("V*")
        plt.title("Normalized Veloctiy Plot for \u03B5=" + str(W) + "(kcal/mol)")
        plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Velocity_Graphs/" + str(W) + " velocity_plot.png", dpi=600)
        plt.close()

    # Combined Plot of normalized velocities
    fig, ax = plt.subplots(2, 4, facecolor='w', figsize=(10,5))
    ax = ax.ravel()
    i = 0
    order_array = [0, 1, 2, 4, 5, 6]
    max_array = np.zeros([6,2])
    for W in numerical_wet_list:
        # plot of Energy diagrams for bordeline cases
        sliced_ps = ps.xs(("velocity", True, W), level=("Property", "Borderline", "W"))

        for index, value in sliced_ps.items():
            # key parameters
            A, T = float(index[0]), float(index[1])

            factor = (A / T * 2 * np.pi)
            y_data = value[1:, 1] / factor

            y_data_change = y_data[:-1] - y_data[1::1]
            liftoff_pos = np.argmax(y_data_change > 0)
            x_data = value[1:, 0] / (value[liftoff_pos, 0])

            max_array[i,:] = [W, value[liftoff_pos, 1] / factor - value[1, 1] / factor]
            # print(value[liftoff_pos, 0])

            ax[order_array[i]].plot(x_data, y_data, linewidth=1.0, label=str(A) + " - " + str(T))

        ax[order_array[i]].set_xlim([0, 6])
        ax[order_array[i]].set_ylim([-1.2, 1.2])

        if i == 3 or i == 0:
            ax[order_array[i]].set_ylabel("V*")
        else:
            ax[order_array[i]].set_yticklabels([])
            # ax[order_array[i]].set_yticks([])
            # ax[order_array[i]].tick_params(axis="y",direction="in", pad=-22)

        if i == 3 or i == 4 or i == 5:
            ax[order_array[i]].set_xlabel("T*")
        else:
            ax[order_array[i]].set_xticklabels([])

        ax[order_array[i]].set_title("\u03B5="+str(W) + "(kcal/mol)")
        i += 1

    ax[7].scatter(max_array[:,0], max_array[:,1], linewidth=1.0)
    ax[7].set_xlabel("\u03B5 (kcal/mol)")
    ax[7].set_ylabel("\u0394V*")
    ax[7].set_title("\u0394V* vs \u03B5")
    ax[7].set_ylim([0,2.2])
    fig.delaxes(ax[3])
    fig.suptitle("Normalized Velocity against Normalized Time")

    # moving to the right
    box = ax[7].get_position()
    box.x0 = box.x0 + 0.05
    box.x1 = box.x1 + 0.05
    ax[7].set_position(box)

    plt.savefig(directories["Local_Sims Base"] + "Overview_Graphs/Velocity_Graphs/multi_velocity_plot.png",
                dpi=600, bbox_inches='tight'), plt.clf()

def Borderline_Test(borderlines, W, A, T):
    for row in borderlines:
        if W == row[0] and A == row[1] and T == row[2]:
            return True
    return False

directories = {
    "Local_Test": "/home/cameron/Documents/Disertation_Work/test_dir/rigid_particle/",
    "Local_Sims Base": "D:/Simulation_Data/rigid_particle/",
    "Local_Sims": "D:/Simulation_Data/rigid_particle/",
    "Local_Long_Sims": "/home/cameron/Documents/Disertation_Work/sim_files/LongRunTemplate/LongRunTemplate/",
    "HPC_Sims": "s1731890@cirrus.epcc.ed.ac.uk:/lustre/home/sc062/s1731890/sims/",

}
wet_colour = {
    0.3: 'b',
    0.4: 'g',
    0.5: 'r',
    0.6: 'c',
    0.7: 'm',
    0.8: 'y'
}
wet_list = ["wet_0.3", "wet_0.4", "wet_0.5", "wet_0.6", "wet_0.7", "wet_0.8"]
numerical_wet_list = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
A_list = [5.0, 6.0]
energy_well = { 0.3:-128.075, 0.4:-173.012, 0.5:-217.097, 0.6:-259.48, 0.7:-306.839, 0.8:-351.421}

ps = make_series(directories, all = False)
ps.index = ps.index.swaplevel(1, 3); ps = ps.sort_index(); ps.index = ps.index.swaplevel(1, 3)
## ps.to_csv(directories["Local_Sims Base"]+"data_dump.csv")

# corners_plot()
energy_scatter_plot(energy_well)
# force_plot()
# energy_plot()
# force_and_energy_plot(numerical_wet_list)
# normalized_velocity_plot(numerical_wet_list)
# theta_plot()

# energy_well_plot(numerical_wet_list)
# border_line_liftoff_plot(numerical_wet_list)
# F_adhesion_gravity()

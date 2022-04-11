import subprocess
import numpy as np


def individual_energy_well(location_individual):
    file1 = open(location_individual + '/forceY.txt', 'r')
    Lines = file1.readlines()

    if Lines[0][0] == '#':
        subprocess.call("sed -i '1d' ./forceY.txt", shell=True, cwd=location_individual)
        subprocess.call("sed -i '1d' ./forceY.txt", shell=True, cwd=location_individual)
        subprocess.call('g++ -o work_new workDone.cpp', shell=True, cwd=location_individual + "/postprocess")
        subprocess.call('./work_new', shell=True, cwd=location_individual + "/postprocess")

    data = np.loadtxt(location_individual + "/postprocess/time_vs_energyInput.txt", delimiter="\t")
    E_individual = data[-1, 1]
    print("The individual energy: " + str(E_individual), end="\t")
    return E_individual

def individual_energy_well_array(location_individual):
    file1 = open(location_individual + '/forceY.txt', 'r')
    Lines = file1.readlines()

    if Lines[0][0] == '#':
        subprocess.call("sed -i '1d' ./forceY.txt", shell=True, cwd=location_individual)
        subprocess.call("sed -i '1d' ./forceY.txt", shell=True, cwd=location_individual)
        subprocess.call('g++ -o work_new workDone.cpp', shell=True, cwd=location_individual + "/postprocess")
        subprocess.call('./work_new', shell=True, cwd=location_individual + "/postprocess")

    data = np.loadtxt(location_individual + "/postprocess/time_vs_energyInput.txt", delimiter="\t")

    return data


def overall_energy_well(location_overall):
    data = np.loadtxt(location_overall + "energyWell.txt", delimiter=" ")
    print(location_overall)
    E_Overall = data[-1, 1]
    print(E_Overall)
    return E_Overall


def energy_well_threshold(sim, directories, threshold_frac):
    is_close = False
    E_Overall = overall_energy_well(directories["Local_Long_Sims"])
    E_individual = individual_energy_well(directories["Local_Sims"] + sim)
    print("The total energy difference:", np.abs(E_Overall - E_individual), end="\t")
    print("The frac of total energy difference:", np.abs(E_Overall / threshold_frac), end="\t")
    if np.abs(E_Overall - E_individual) < np.abs(E_Overall / threshold_frac):
        is_close = True
    return is_close, np.abs(E_Overall - E_individual), E_Overall


def Lift_off_test(sim, TimePeriod, directories, energy_well, energy_diff):
    lift_off = False
    data = np.loadtxt(directories["Local_Sims"] + sim + "/forceY.txt", delimiter=" ")

    # This tests if the rebound happened:
    Test_A = (np.average(data[int(1/6*TimePeriod):, 1]) < 0.01)
    Test_B = (np.average(data[int(1/6*TimePeriod):int(TimePeriod), 1]) < 0.01)
    print(TimePeriod, " ", len(data[:,1]))
    Test_C =  np.max(data[int(1/3*TimePeriod):, 1]) < 50

    no_rebound =  Test_A and Test_B and Test_C

    energy_data = individual_energy_well_array(directories["Local_Sims"] + sim)
    #
    #
    # print("\nmin_value:", np.min(np.abs(data[int(1/2*TimePeriod):,1])))
    # print("max_value:", np.max(data[int(1 / 2 * TimePeriod):, 1]))
    # print("min energy:", np.min(energy_data[int(1 / 4 * TimePeriod):int(TimePeriod),1]))
    # print("min location:", np.argmin(energy_data))
    # print("energy_well:", energy_well)
    # print("energy_at_D:", energy_data[int(1/4*TimePeriod), 1])
    # print("Time period:", sim)

    #This tests for zero force at end of sim:
    no_force = (data[-1, 1] == 0.0)

    if no_force == True and no_rebound == True:
        lift_off = True

    # print("the min energy: ", np.min(data[:, 1]), ", the energy well: ", energy_well)
    # if np.min(data[:, 1]) < energy_well*3:
    #     lift_off = True

    print("is there lift off: ", lift_off, end="\t")
    return lift_off


def F_main(full_sims, threshold_frac, directories):


    print("F-Main started")
    # defining a previous sim that didn't lift_off
    energy_dict = {}
    previous_value = 10000

    for sim in reversed(full_sims):
        is_close, energy_diff, energy_well = energy_well_threshold(sim, directories, threshold_frac)
        lift_off = Lift_off_test(sim, full_sims[sim][1], directories, energy_well, energy_diff)
        energy_dict[sim] = energy_diff

        if lift_off == True:
            break

        previous_value = full_sims[sim][1]

    return previous_value, full_sims[sim][1]
        # print(sim, lift_off)
        # if Lift_off_test(sim, full_sims[sim][1], directories, energy_well) != True or -energy_diff/energy_well > 2:
        #     energy_dict[sim] = 1000




    #     else:
    #         energy_dict[sim] = energy_diff
    #         is_Final = True
    #         return previous_sim, sim, is_Final
    #
    # min_value = 1000
    # current_sim = "First"
    # print(energy_dict)
    #
    # for sim in reversed(sorted(energy_dict)):
    #     # FIND A WAY TO DEAL WITH SINGLE LIFTOFFF CASES # Min dictionary approach?
    #     energy_diff = energy_dict[sim]
    #
    #     if energy_diff > min_value:
    #         # This returns the sim before min energy sim:
    #         is_Final = False #Because if it was below the threshold we would have exited already
    #         return previous_sim, current_sim, is_Final
    #
    #     if energy_diff < min_value:
    #         min_value = energy_diff
    #
    #
    #     previous_sim = current_sim
    #     current_sim = sim
    #
    # is_Final = False
    # return previous_sim, current_sim, is_Final
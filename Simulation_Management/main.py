import GH_Code
import F_Code
import I_Code
import misc_functions as MF

wet_params = {
    "TS_meas": 100000,
    "TS_init": 4100000,
    "wetability": 'wet_0.4'
}

directories = {
    "Local_Test": "/home/cameron/Documents/Disertation_Work/test_dir/rigid_particle/",
    "Local_Sims Base": "/media/cameron/Data_Storage/Simulation_Files/rigid_particle/",
    "Local_Sims": "/media/cameron/Data_Storage/Simulation_Files/rigid_particle/",
    "Local_Long_Sims": "/home/cameron/Documents/Disertation_Work/sim_files/LongRunTemplate/LongRunTemplate/",
    "HPC_Sims": "**************************",
}
#Amplitude = 4.0
for wet_params["wetability"] in ["wet_0.3", "wet_0.4", "wet_0.5", "wet_0.6", "wet_0.7", "wet_0.8"]:
    period_lower = 10000
    period_upper = 30000
    directories["Local_Long_Sims"] = directories["Local_Test"] + wet_params["wetability"] + "/"

    for Amplitude in [4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]:
        period_upper += 20000
        directories["Local_Sims"] = directories["Local_Sims Base"] + wet_params["wetability"] + "/" + str(Amplitude) + "/"

        # checks for existing sims in directory
        period_lower, period_upper = MF.Check_Existing_sims(directories["Local_Sims"], Amplitude, period_lower, period_upper)

        escape_loop = ((period_upper - period_lower) <= 1000)
        is_Final = False
        i = 0
        while is_Final == False and (not escape_loop) and i <= 8:
            i += 1
            escape_loop = ((period_upper - period_lower) <= 1000)

            # Create local sim skeletons for analysis
            dir_name_dict = GH_Code.making_sim_directories(period_lower, period_upper, Amplitude, wet_params, directories)

            # Clears HPC folder
            cmd = ["cd sims", "rm -r Temp", "mkdir Temp"]
            MF.ssh_multi(cmd)

            # Copy to HPC
            MF.ssh_scp_to_HPC('scp -r ' + directories["Local_Test"] + "Temp/ " + directories["HPC_Sims"])

            # Time check
            MF.TimeCheck()

            # Run on HPC
            MF.run_hpc(delay_time=25)

            success = False
            while not success:
                # Clearing Local Directory
                GH_Code.clearing_local_dir(directories, str(Amplitude), Wetability=wet_params["wetability"])

                # Copy finished sims to Local Sims
                MF.ssh_scp('scp -r ' + directories["HPC_Sims"] + "Temp/* " + directories["Local_Sims"], timeout=6000)

                # Tests if copy over was successful
                success = MF.Check_Sims_Finished(directories["Local_Sims"])

                if not success:
                    print("Waiting until all results downloaded:")
                    MF.sleep_timer(20, 5)

            # Tests sims for lift_off/final solution
            period_upper, period_lower = F_Code.F_main(full_sims=dir_name_dict, threshold_frac=25,
                                                        directories=directories)

            print(period_upper, period_lower)

            if period_upper - period_lower < 1000:
                period_lower = period_upper - 1000
            success = True

            print("The upper value is: ", period_upper, " and the lower value is: ", period_lower)



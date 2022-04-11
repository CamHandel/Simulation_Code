import subprocess
import numpy as np

def replace(file, location, text, value):
    # Time period inclusion
    command_string = 'sed \'s/'+text+'/' + value + '/g\' ' + file + ' > tmp.file'
    subprocess.call(command_string, shell=True, cwd=location)
    subprocess.call('mv tmp.file ' + file , shell=True, cwd=location)

def make_dir(A, time_period, base_dir, TS_meas, TS_init, wetability):
    dir_name_dict = {}
    for T in time_period:
        dir_name = 'A' + str(A) + '-T' + str(T)
        dir_name_dict[dir_name] = [A, T]

    for name in dir_name_dict:
        # making directory
        subprocess.call('mkdir ' + name, shell=True, cwd=base_dir + "Temp/")

        # copy over restart file:
        subprocess.call('cp ' + base_dir + wetability +'/restart.' + str(TS_init) + ' ' + base_dir + "Temp/" + name + '/', shell=True, cwd=base_dir)

        # copy over post processing script
        #subprocess.call('cp ' + base_dir + '/EquilSetup/postprocess/t_meas.cpp ' + base_dir + name + '/', shell=True,
        #                cwd=base_dir)

        # copy over measure specific .in file, properties_file and runscript
        subprocess.call('cp ' + base_dir + 'Measure_Files/in.meas ' + base_dir + "Temp/" +name + '/', shell=True,
                       cwd=base_dir)
        subprocess.call('cp ' + base_dir + 'Measure_Files/runscript.slurm ' + base_dir + "Temp/" +name + '/', shell=True,
                        cwd=base_dir)
        subprocess.call('cp ' + base_dir + 'Measure_Files/SimProperties ' + base_dir + "Temp/" +name + '/', shell=True,
                        cwd=base_dir)

        # copy over post process folder
        subprocess.call('cp -r ' + base_dir + 'Measure_Files/postprocess ' + base_dir + "Temp/" +name + '/', shell=True,
                        cwd=base_dir)

        # rename the measure values
        A = str(dir_name_dict[name][0])
        T = str(dir_name_dict[name][1])
        TS_meas = str(TS_meas)
        TS_init = str(TS_init)
        # This is now manual!!!!
        WET = "0.8"
        # WET = wetability[4:]

        # Time period inclusion
        meas_dir = base_dir + "Temp/" + name + '/'
        replace(file='in.meas', location=meas_dir, text='AMP', value=A)
        replace(file='in.meas', location=meas_dir, text='TIMEP', value=T)
        replace(file='in.meas', location=meas_dir, text='NTSTEPS', value=TS_meas)
        replace(file='in.meas', location=meas_dir, text='INIT', value=TS_init)
        replace(file='in.meas', location=meas_dir, text='WET', value=WET)

        # remove top of force Y


        # Post Processing: Amplitude and Time Period inclusion
        post_dir = base_dir + "Temp/" + name + '/' + 'postprocess/'
        #print("This is the post dir:"+post_dir)
        replace(file='workDone.cpp', location=post_dir, text='AMP', value=A)
        replace(file='workDone.cpp', location=post_dir, text='TIMEP', value=T)
        replace(file='workDone.cpp', location=post_dir, text='NTSTEPS', value=TS_meas)

    return dir_name_dict

def rsh_write(dir_name_dict):
    with open('run.sh', 'w') as rsh:
        for name in dir_name_dict:
            command = 'cd ' + name
            rsh.write('cd ' + name + '\n')
            rsh.write('sbatch runscript.slurm \n')
            rsh.write('cd .. \n')
    rsh.close()

def clearing_local_dir(directories, Amplitude, Wetability):
    subprocess.call('rm -r ' + str(Amplitude), shell=True, cwd=directories["Local_Sims Base"] + Wetability)
    subprocess.call('mkdir ' + str(Amplitude), shell=True, cwd=directories["Local_Sims Base"] + Wetability)

def clearing_test_dir(directories):
    subprocess.call('rm -r Temp', shell=True, cwd=directories["Local_Test"])
    subprocess.call('mkdir Temp', shell=True, cwd=directories["Local_Test"])

def making_sim_directories(time_period_lower, time_period_upper, amplitude, wet_params, directories):
    # Setting up range of sims to run
    clearing_test_dir(directories)

    n_dir = 11
    time_period = np.linspace(time_period_lower, time_period_upper, n_dir).astype(int)
    TS_meas = wet_params["TS_meas"]
    TS_init = wet_params["TS_init"]
    wetability = wet_params["wetability"]

    dir_name_dict = make_dir(amplitude, time_period, directories["Local_Test"], TS_meas, TS_init, wetability)

    # copys script to base dir
    rsh_write(dir_name_dict)
    subprocess.call('cp run.sh ' + directories["Local_Test"] + "Temp/", shell=True)

    return dir_name_dict
    # '{test,runscript.slurm}'

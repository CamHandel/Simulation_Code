import pexpect
import time
import tempfile
import os
import numpy as np

def ssh(cmd, timeout=30, bg_run=False):
    """SSH'es to a host using the supplied credentials and executes a command.
   Throws an exception if the command doesn't return 0.
   bgrun: run command in the background"""
    user = "******************"
    host = "******************"
    password = "******************"

    ssh_cmd = 'ssh %s@%s' % (user, host)
    child = pexpect.spawn(ssh_cmd, timeout=timeout)
    fout = open('mylog.txt','wb')
    child.logfile = fout
    child.expect(['[pP]assword: '])
    child.sendline(password)
    child.expect(["\$"])
    child.sendline(cmd)
    child.expect(["\$"])
    child.close()
    fout.close()

def ssh_multi(cmd, timeout=60, bg_run=False):
    """SSH'es to a host using the supplied credentials and executes a command.
   Throws an exception if the command doesn't return 0.
   bgrun: run command in the background"""
    user = "******************"
    host = "******************"
    password = "******************"
    print("Multi-Command SSH started")

    ssh_cmd = 'ssh %s@%s' % (user, host)
    child = pexpect.spawn(ssh_cmd, timeout=timeout)
    fout = open('mylog.txt','wb')
    child.logfile = fout
    child.expect(['[pP]assword: '])
    child.sendline(password)

    for sub_cmd in cmd:
        child.expect(["\$"])
        child.sendline(sub_cmd)

    child.expect(["\$"])
    child.close()
    fout.close()

def ssh_scp(copy_cmd, timeout=120, bg_run=False):
    user = "******************"
    host = "******************"
    password = "******************"

    print("Copying for SSH started")

    #ssh_cmd = 'ssh %s@%s' % (user, host)
    child = pexpect.spawn(copy_cmd, timeout=timeout)
    fout = open('scp_mylog.txt','wb')
    child.logfile = fout
    child.expect(['[pP]assword: '])
    child.sendline(password)
    #time.sleep(timeout)
    child.expect(["run.sh"])
    child.close()
    fout.close()

def ssh_scp_to_HPC(copy_cmd, timeout=4, bg_run=False):
    user = "******************"
    host = "******************"
    password = "******************"

    print("To HPC copying started")

    child = pexpect.spawn(copy_cmd, timeout=timeout*60)
    fout = open('scp_mylog.txt','wb')
    child.logfile = fout
    child.expect(['[pP]assword: '])
    child.sendline(password)
    child.expect(["run.sh"])
    time.sleep(30)
    #sleep_timer(timeout, 1)
    child.close()
    fout.close()

def F_function_error(previous_sim):
    if previous_sim == "Last":
        # this means there is no liftoff in any of the sims, expanded range needed
        raise ValueError("The scope is too small for this sim, no liftoff reached")

    elif previous_sim == "First":
        raise ValueError("The scope is too small for this sim, liftoff is achived in the first sim")

def sleep_timer(minutes, interval):
    n = minutes/interval
    print("Total time: ")
    for i in range(int(n)-1):
        time.sleep(interval*60)
        print(i*interval+interval, end=" ")
    print(int(n)*interval)

def run_hpc(delay_time):
    print("Running on HPC started")
    cmd = ["cd sims/Temp/", "chmod a+x run.sh", "./run.sh"]
    ssh_multi(cmd)
    sleep_timer(delay_time, 5)

def TimeCheck():
    if time.gmtime().tm_hour > 23 and time.gmtime().tm_min > 30:
        time.sleep(31*60)

def Check_Existing_sims(dir, Amplitude, period_lower, period_upper):
    values = []
    for (dirpath, dirnames, filenames) in os.walk(dir):
        for dirname in dirnames:
            if dirname != "postprocess" and dirname != "Graphs":
                value = int(dirname[dirname.index('T')+1:])
                values.append(value)
    array = np.array(values)
    if len(array) == 0:
        return period_lower, period_upper

    elif array.min() <= period_lower or array.max() >= period_upper:
        return period_lower, period_upper

    return array.min(), array.max()

def Check_Sims_Finished(dir):
    n_force = 0
    n_runscript = 0
    for (dirpath, dirnames, filenames) in os.walk(dir):
        for file in filenames:
            path = os.path.join(dirpath, file)
            size = os.stat(path).st_size

            if file == "forceY.txt" and size > 900000:
                n_force += 1
            if file == "runscript.slurm":
                n_runscript += 1

    print("Number of files (force and slurm): " + str(n_force) + " and " + str(n_runscript))
    if n_force == n_runscript and n_force != 0:
        return True
    else:
        return False
import shutil
import os
import random
from numpy import arange

def create_call(L, r, seed, tau, M, crit):
    str_base = "./exe"
    call_params = " {} {} {} {} {} {}".format(L,
                                              r, seed,
                                              tau, M,
                                              crit)
    return str_base + call_params

def copy_header(job_name):
    shutil.copyfile('./job.sh', job_name)

def write_job(job_name, call_string):
    with open(job_name, "a") as job:
        job.write(call_string + '\n')

def generate_calls():
    calls = []
    rs = [round(r,4) for r in arange(r_i, r_f+dr, dr)]
    seed = first_seed
    for r in rs:
        for _ in range(samples):
            calls.append(create_call(L, r, seed, tau, M, crit))
            seed += 1

    return calls


L = 50
tau = 1.1
M = 10
crit = 0.5
max_jobs = 20

########################################################

samples = 50 # num of samples for each r
r_i = 0.5 # first r
r_f = 2.0 # last r
dr = 0.1 # r increment
first_seed = 9075421428

# total number of simulations
num_simuls = samples*(r_f + dr - r_i)/dr

submits_per_job = num_simuls / max_jobs
# print("submits per job:", submits_per_job)

########################################################

if not os.path.exists("JOBS"):
    os.mkdir("JOBS")

calls = generate_calls()
random.shuffle(calls)

num_submits = 0
num_jobs = 1
job_name = "run_{}.sh".format(num_jobs)
copy_header(job_name)
for call in calls:
    write_job(job_name, call)
    num_submits += 1
    if num_submits >= submits_per_job:
        num_submits = 0
        shutil.move(job_name, "JOBS")
        num_jobs += 1
        if num_jobs > max_jobs:
            break
        job_name = "run_{}.sh".format(num_jobs)
        copy_header(job_name)
if os.path.exists(job_name):
    shutil.move(job_name, "JOBS")

########################################################

# with open("calls.log", "w") as f:
#     for call in calls:
#         f.write(call+"\n")

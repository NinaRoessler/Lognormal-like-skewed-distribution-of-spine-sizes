#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 00:17:28 2023

@author: Alex
"""
from sbi import analysis as analysis
from sbi import utils as utils
from sbi.utils.get_nn_models import posterior_nn
from sbi.inference import SNPE, prepare_for_sbi, simulate_for_sbi
import torch
import numpy as np
from numba import jit
from scipy import stats

@jit(nopython=True)
def intrinsic_noise(n_neuron,n_spine,n_step,tS,kn,gn):
    spine_sizes=np.log(np.clip(tS/n_spine+0.01*np.random.randn(n_spine,n_neuron),0,None))
    spine_sizes*=tS/np.sum(spine_sizes,0) # Heterosynaptic part
    for dt_ind in range(n_step):
        spine_sizes=np.clip(spine_sizes*(np.ones((n_spine,n_neuron))+kn*np.random.randn(n_spine,n_neuron))+gn*np.random.randn(n_spine,n_neuron),0,None) # Stochastic growth
        spine_sizes*=tS/np.sum(spine_sizes,0) # Heterosynaptic part
    return spine_sizes

@jit(nopython=True)
def extrinsic_noise(spine_sizes,n_step,tS,kn,gn,
                    steps_per_in,bck_pr,n_act,
                    threshold,tau_m,
                    tau_STDP,dep_rate,pot_rate
                    ):
    n_neuron=spine_sizes.shape[1]
    n_spine=spine_sizes.shape[0]
    
    net_E=threshold*np.random.rand(n_neuron) # (initially uniform between reset and threshold)
    pre_times=np.zeros((n_spine,n_neuron)) #% Time since last EPSP at each spine
    post_times=np.zeros(n_neuron) # % Time since last spike at each neuron
    since_spike=0
    
    for dt_ind in range(n_step):
        since_spike+=1
        net_E*=np.exp(-1/tau_m) # Membrane leak
    
        if since_spike%steps_per_in==0:
            net_E+=np.sum(spine_sizes[:n_act:],0)
            pre_times[:n_act,:]=0
            
        bc=np.random.rand(n_spine,n_neuron) # Background
        S=np.zeros((n_spine,n_neuron))
        for ii in range(n_spine):
            for jj in range(n_neuron):
                if bc[ii,jj]<bck_pr:
                    S[ii,jj]=1
                    pre_times[ii,jj]=0
                else:
                    pre_times[ii,jj]+=1  # Increment times since EPSPs
        net_E+=np.sum(spine_sizes*S,0)
        
        for jj in range(n_neuron):
            if net_E[jj]>=threshold:
                net_E[jj]=0 # % Reset after spikes
                post_times[jj]=0 # % Reset after spikes
        
       # STDP
        spine_sizes+=pot_rate*np.exp(-pre_times/tau_STDP)# % Potentiate
        for ii in range(n_spine):
            for jj in range(n_neuron):
                spine_sizes[ii,jj]*=(1-dep_rate*np.exp(-post_times[jj]/tau_STDP)) # Depress
        post_times+=1 # Increment times since spikes
            
        spine_sizes=np.clip(spine_sizes*(np.ones((n_spine,n_neuron))+kn*np.random.randn(n_spine,n_neuron))+gn*np.random.randn(n_spine,n_neuron),0,None) # Stochastic growth
        spine_sizes*=tS/np.sum(spine_sizes,0) # Heterosynaptic part
    return spine_sizes

S=intrinsic_noise(500,25,10,100,0.001,0)
E=extrinsic_noise(S,10,100,0.001,0,100,0.001,5,10,100,100,0.1,0.1)

def granule_STDP_fst(tau_STDP_rw, # STDP timescale
                          tS, # Total size
                          kn, # Kesten noise (multi)
                          gn,# Gauss noise (add)
                          pr_act, # Prop. spines receiving inputs
                          pot_rate, # Potentiation rate
                          dep_rate): # Depression rate
    dt=0.0005
    
    max_time_1=25
    max_time_2=100
    
    n_neuron=500
    n_spine=25
    tau_m_rw=20
    threshold=10
    
    back_rate=5
    input_rate_1=10
    input_rate_2=200
    
    n_step_1=round(max_time_1/dt)
    n_step_2=round(max_time_2/dt)
    steps_per_in_1=round(input_rate_1/dt)
    steps_per_in_2=round(input_rate_2/dt)
    bck_pr=dt*back_rate
    n_act=round(pr_act*n_spine)
    tau_m=tau_m_rw/(1000*dt)
    tau_STDP=tau_STDP_rw/(1000*dt)
    
    Si=intrinsic_noise(n_neuron,n_spine,n_step_1,tS,kn,gn)
    Se_10=extrinsic_noise(Si,n_step_2,tS,kn,gn,
                       steps_per_in_1,bck_pr,n_act,
                       threshold,tau_m,
                       tau_STDP,dep_rate,pot_rate
                       )
    Se_200=extrinsic_noise(Si,n_step_2,tS,kn,gn,
                       steps_per_in_2,bck_pr,n_act,
                       threshold,tau_m,
                       tau_STDP,dep_rate,pot_rate
                       )
    return Si , Se_10 , Se_200





num_dim = 7
lb = torch.ones(num_dim)
ub = torch.ones(num_dim)

# STDP timescale (log10)
lb[0] = 0
ub[0] = 5

# Total_size (log10)
lb[1] = 1
ub[1] = 3

# Kesten noise (log10)
lb[2] = -10
ub[2] = -3

# Gauss noise (linear)
lb[3] = 0
ub[3] = 0.0001

# Prop spines
lb[4] = 0
ub[4] = 1

# Pot rate
lb[5] = 0
ub[5] = 0.35

# Dep rate
lb[6] = 0
ub[6] = 0.15

prior = utils.BoxUniform(low=lb, high=ub)


def simulator_GC(prmst):
    Si , Se_10 , Se_200 = granule_STDP_fst(10**prmst[0].item(), # STDP timescale
                              10**prmst[1].item(), # Total size
                              10**prmst[2].item(), # Kesten noise (multi)
                              prmst[3].item(),# Gauss noise (add)
                              prmst[4].item(), # Prop. spines receiving inputs
                              prmst[5].item(), # Potentiation rate
                              prmst[6].item()) # Depression rate

    mean_sil=Si.mean()
    mean_data=0.3990
    Si*=mean_data/mean_sil

    mean_10=Se_10.mean()   
    mean_std=0.1372;
    Se_10*=mean_std/mean_10
    
    mean_200=Se_200.mean()  
    mean_HFS=0.1162;
    Se_200*=mean_HFS/mean_200
   
    prm_sil=stats.lognorm.fit(Si.flatten(),floc=0)
    prm_10=stats.lognorm.fit(Se_10.flatten(),floc=0)
    prm_200=stats.lognorm.fit(Se_200.flatten(),floc=0)
   
    prms_out=torch.zeros(6,)
    
    prms_out[0]=prm_sil[0]
    prms_out[1]=prm_sil[2]
    prms_out[2]=prm_10[0]
    prms_out[3]=prm_10[2]
    prms_out[4]=prm_200[0]
    prms_out[5]=prm_200[2]
    
    
    return prms_out


num_rounds = 5
# The specific observation we want to focus the inference on.
x_o = torch.zeros(6,)
x_o[0] = 0.6976
x_o[1] = 0.3368
x_o[2] = 0.6326
x_o[3] = 0.1122
x_o[4] = 0.7564
x_o[5] = 0.0902

simulator, prior = prepare_for_sbi(simulator_GC, prior)
inference = SNPE(prior=prior)

posteriors = []
proposal = prior

for _ in range(num_rounds):
    theta, x = simulate_for_sbi(
        simulator, proposal, num_simulations=4096,num_workers=32)

    # In `SNLE` and `SNRE`, you should not pass the `proposal` to `.append_simulations()`
    density_estimator = inference.append_simulations(
        theta, x, proposal=proposal
    ).train()
    posterior = inference.build_posterior(density_estimator)
    posteriors.append(posterior)
    proposal = posterior.set_default_x(x_o)

##
posterior_samples = posterior.sample((10000,), x=x_o)
_ = analysis.pairplot(
    posterior_samples, limits=[[lb[0].item(), ub[0].item()],
                                [lb[1].item(), ub[1].item()],
                                [lb[2].item(), ub[2].item()],
                                [lb[3].item(), ub[3].item()],
                                [lb[4].item(), ub[4].item()],
                                [lb[5].item(), ub[5].item()],
                                [lb[6].item(), ub[6].item()],],
    figsize=(10, 10)
    )

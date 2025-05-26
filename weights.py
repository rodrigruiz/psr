#!/usr/bin/env python

import numpy as np

import km3flux  # from km3flux.flux import Honda2015  # from km3flux.flux import flux
# from km3services.oscprob import OscProb
import ROOT

# Attention: weight_one_year actually in unit 1/s

# osc_prob_dict
oscprob_dict = {14: 1, -14: 1, 12: 0, -12: 0, 16: 2, -16: 2}


# Hondaflux
def get_honda_flux(particle_type, energy, cos_zenith):
    honda = km3flux.flux.Honda()
    f = honda.flux(2014, "Frejus", solar="min", averaged="azimuth")

    flux = f[particle_type](energy, cos_zenith)

    return flux

def set_oscillation_parameters(p, osc_dict):
    p.SetAngle(1, 2, osc_dict['th12'] * np.pi / 180)
    p.SetAngle(2, 3, osc_dict['th23'] * np.pi / 180)
    p.SetAngle(1, 3, osc_dict['th13'] * np.pi / 180)

    p.SetDm(2, osc_dict['dm21'])
    p.SetDm(3, osc_dict['dm31'])

    p.SetDelta(1, 3, osc_dict['delta_cp'] * np.pi / 180)

    return p

# OscProb
def get_osc_prob(flv_in, flv_out, energy, cos_zenith):
    # oscprob = OscProb()
    # prob = oscprob.oscillationprobabilities(flv_in, flv_out, energy, dir_z)
    # oscillation parameters from NuFit v5.0:
    # osc_pars_NO = {'th12': 33.44, 'th23': 49.2, 'th13': 8.57, 'dm21': 7.42e-5, 'dm31': 2.517e-3, 'delta_cp': 197}
    # oscillation parameters from NuFit v6.0:
    osc_pars_NO = {'th12': 33.68, 'th23': 48.5, 'th13': 8.52, 'dm21': 7.49e-5, 'dm31': 2.534e-3, 'delta_cp': 177}

    #ROOT.gSystem.Load('/sps/km3net/users/ngeissel/Code/Swim2.0/lib/libOscProb.so')
    #prem = ROOT.OscProb.PremModel('Code/Swim2.0/OscProb/PremTables/prem_15layers.txt')

    ROOT.gSystem.Load('/Swim/lib/libOscProb.so')
    prem = ROOT.OscProb.PremModel('/home/hpc/capn/capn107h/software/psr/flux_weighting/prem_15layers.txt')

    oscprob_idx_flv_in = oscprob_dict[flv_in]
    oscprob_idx_flv_out = oscprob_dict[flv_out]

    p = ROOT.OscProb.PMNS_Fast()
    p.SetAvgProbPrec(0.01)

    set_oscillation_parameters(p, osc_pars_NO)
    #print(flv_in)
    if flv_in < 0:
        p.SetIsNuBar(True)

    prob = np.zeros(cos_zenith.shape)
    for i, c in enumerate(cos_zenith):
        prem.FillPath(c)
        p.SetPath(prem.GetNuPath())
        prob[i] = p.Prob(oscprob_idx_flv_in, oscprob_idx_flv_out, energy[i])

    return prob


# oscillated weight
def get_weight(is_cc, pdg, weight_one_year, energy, dir_z):
    if is_cc:
        return get_weight_cc(pdg, weight_one_year, energy, dir_z)
    else:
        return get_weight_nc(pdg, weight_one_year, energy, dir_z)


def get_weight_cc(pdg, weight_one_year, energy, dir_z):
    if pdg == 12:
        w_osc = _get_weight_nue(energy, dir_z)
    elif pdg == -12:
        w_osc = _get_weight_anue(energy, dir_z)
    elif pdg == 14:
        w_osc = _get_weight_numu(energy, dir_z)
    elif pdg == -14:
        w_osc = _get_weight_anumu(energy, dir_z)
    elif pdg == 16:
        w_osc = _get_weight_nutau(energy, dir_z)
    elif pdg == -16:
        w_osc = _get_weight_anutau(energy, dir_z)
    else:
        w_osc = np.nan

    w = weight_one_year * w_osc
    return w


def get_weight_nc(pdg, weight_one_year, energy, dir_z):
    if pdg == 14:
        w = _get_weight_nu_nc(energy, dir_z)
    elif pdg == -14:
        w = _get_weight_anu_nc(energy, dir_z)
    else:
        w = np.nan

    return weight_one_year * w


def _get_weight_nue(energy, dir_z):
    flux_numu = get_honda_flux('numu', energy, dir_z)
    flux_nue = get_honda_flux('nue', energy, dir_z)

    p_numu_nue = get_osc_prob(14, 12, energy, dir_z)
    p_nue_nue = get_osc_prob(12, 12, energy, dir_z)

    w_osc = flux_numu * p_numu_nue + flux_nue * p_nue_nue
    return w_osc


def _get_weight_anue(energy, dir_z):
    flux_anumu = get_honda_flux('anumu', energy, dir_z)
    flux_anue = get_honda_flux('anue', energy, dir_z)

    p_anumu_anue = get_osc_prob(-14, -12, energy, dir_z)
    p_anue_anue = get_osc_prob(-12, -12, energy, dir_z)

    w_osc = flux_anumu * p_anumu_anue + flux_anue * p_anue_anue
    return w_osc


def _get_weight_numu(energy, dir_z):
    flux_numu = get_honda_flux('numu', energy, dir_z)
    flux_nue = get_honda_flux('nue', energy, dir_z)

    p_numu_numu = get_osc_prob(14, 14, energy, dir_z)
    p_nue_numu = get_osc_prob(12, 14, energy, dir_z)

    w_osc = flux_numu * p_numu_numu + flux_nue * p_nue_numu
    return w_osc


def _get_weight_anumu(energy, dir_z):
    flux_anumu = get_honda_flux('anumu', energy, dir_z)
    flux_anue = get_honda_flux('anue', energy, dir_z)

    p_anumu_anumu = get_osc_prob(-14, -14, energy, dir_z)
    p_anue_anumu = get_osc_prob(-12, -14, energy, dir_z)

    w_osc = flux_anumu * p_anumu_anumu + flux_anue * p_anue_anumu
    return w_osc


def _get_weight_nutau(energy, dir_z):
    flux_numu = get_honda_flux('numu', energy, dir_z)
    flux_nue = get_honda_flux('nue', energy, dir_z)

    p_numu_nutau = get_osc_prob(14, 16, energy, dir_z)
    p_nue_nutau = get_osc_prob(12, 16, energy, dir_z)

    w_osc = flux_numu * p_numu_nutau + flux_nue * p_nue_nutau
    return w_osc


def _get_weight_anutau(energy, dir_z):
    flux_anumu = get_honda_flux('anumu', energy, dir_z)
    flux_anue = get_honda_flux('anue', energy, dir_z)

    p_anumu_anutau = get_osc_prob(-14, -16, energy, dir_z)
    p_anue_anutau = get_osc_prob(-12, -16, energy, dir_z)

    w_osc = flux_anumu * p_anumu_anutau + flux_anue * p_anue_anutau
    return w_osc


def _get_weight_nu_nc(energy, dir_z):
    flux_numu = get_honda_flux('numu', energy, dir_z)
    flux_nue = get_honda_flux('nue', energy, dir_z)

    w = flux_numu + flux_nue
    return w


def _get_weight_anu_nc(energy, dir_z):
    flux_anumu = get_honda_flux('anumu', energy, dir_z)
    flux_anue = get_honda_flux('anue', energy, dir_z)

    w = flux_anumu + flux_anue
    return w


# atmospheric (not oscillated) weight --> use w3 for now
def get_weight_no_osc(w3, n_gen):
    # if pdg == 12:
    #     flux = get_honda_flux('nu_e', energy, dir_z)
    # elif pdg == -12:
    #     w_osc = get_honda_flux('anu_e', energy, dir_z)
    # elif pdg == 14:
    #     w_osc = get_honda_flux('nu_mu', energy, dir_z)
    # elif pdg == -14:
    #     w_osc = get_honda_flux('anu_mu', energy, dir_z)
    # elif pdg == 16:
    #     w_osc = _get_weight_nutau(energy, dir_z)
    # elif pdg == -16:
    #     w_osc = _get_weight_anutau(energy, dir_z)
    # else:
    #     w_osc = np.nan

    w = w3 / n_gen
    return w


def get_weights_pid_output(df):
    n_events = len(df['energy'])
    w_osc = np.zeros(n_events)

    energy = np.array(df['energy'])
    dir_z = np.array(df['cos_zenith_true'])  # np.array(df['dir_z'])

    w2 = np.array(df['w2'])
    # w3 = np.array(df['w3'])
    n_gen = np.array(df['n_gen'])
    w_one_sec = w2 / n_gen
    # w_atm = w3 / n_gen

    try:
        pid = np.array(df['__pidClass'])
    except KeyError:
        pid = get_pid_class(df['pdgid'], df['is_cc'])
    pdg = np.array(df['pdgid'])

    idx_nue_cc = np.logical_and(pid == 'elec_cc', pdg > 0)
    idx_anue_cc = np.logical_and(pid == 'elec_cc', pdg < 0)
    idx_numu_cc = np.logical_and(pid == 'muon_cc', pdg > 0)
    idx_anumu_cc = np.logical_and(pid == 'muon_cc', pdg < 0)
    idx_nutau_cc = np.logical_and(pid == 'tau_cc', pdg > 0)
    idx_anutau_cc = np.logical_and(pid == 'tau_cc', pdg < 0)
    idx_numu_nc = np.logical_and(pid == 'muon_nc', pdg > 0)
    idx_anumu_nc = np.logical_and(pid == 'muon_nc', pdg < 0)

    indices = [idx_nue_cc, idx_anue_cc, idx_numu_cc, idx_anumu_cc, idx_nutau_cc, idx_anutau_cc, idx_numu_nc,
               idx_anumu_nc]

    for idx in indices:
        if np.sum(idx) == 0:
            continue
        pdg_id = pdg[idx][0]
        is_cc = np.array(df['is_cc'])[idx][0]
        e = energy[idx]
        z = dir_z[idx]

        w_one = w_one_sec[idx]
        w = get_weight(is_cc, pdg_id, w_one, e, z)  # oscillated weight [1/s]
        w_osc[idx] = w  # [0, :]

    return w_osc


def get_pid_class(pdg, is_cc):
    pid = np.empty(len(pdg), dtype='<U9')
    pdg = np.array(pdg)
    is_cc = np.array(is_cc)

    is_elec_cc = np.logical_and(np.abs(pdg) == 12, is_cc == 1)
    is_muon_cc = np.logical_and(np.abs(pdg) == 14, is_cc == 1)
    is_tau_cc = np.logical_and(np.abs(pdg) == 16, is_cc == 1)
    is_muon_nc = np.logical_and(np.abs(pdg) == 14, is_cc == 0)

    pid[is_elec_cc] = 'elec_cc'
    pid[is_muon_cc] = 'muon_cc'
    pid[is_tau_cc] = 'tau_cc'
    pid[is_muon_nc] = 'muon_nc'

    return pid


def get_run_duration_event(run_ids):
    run_dict = utils.get_runduration_dict()  # run: duration [s]
    run_duration = np.array([run_dict[int(run_ids[i])] for i in range(0, len(run_ids))])  # run duration for each event

    return run_duration


def get_total_run_duration(run_ids):
    run_dict = utils.get_runduration_dict()  # run: duration [s]
    run_set = list(set(run_ids))

    total_run_dur = 0
    for run in run_set:
        total_run_dur += run_dict[int(run)]

    return total_run_dur


def get_weights_one_year(pdgid, is_cc, energy, e_min, e_max, run_id, run_duration, w_osc):

    weights_one_year = w_osc * run_duration

    pid_class = get_pid_class(pdgid, is_cc)
    pid_set = list(set(pid_class))

    # loop over all pids
    for pid in pid_set:
        is_pid = pid_class == pid
        pdgs = list(set(pdgid[is_pid]))

        # loop over nu and anti-nu
        for pdg in pdgs:
            is_pdg = pdgid == pdg
            is_type_all = np.logical_and(is_pid, is_pdg)

            e_min_max = list(set(list(zip(e_min[is_type_all], e_max[is_type_all]))))

            # loop over all different energy ranges
            for e_min_gen, e_max_gen in e_min_max:
                is_energy = np.logical_and(energy <= e_max_gen, energy > e_min_gen)
                is_type = np.logical_and(is_type_all, is_energy)

                runs_type = list(set(run_id[is_type]))

                livetime_type = 0
                for run in runs_type:
                    is_run = run_id[is_type] == run
                    duration_run = run_duration[is_type][is_run][0]

                    if not np.all(run_duration[is_type][is_run] == duration_run):
                        print('Problem with run durations')

                    livetime_type += duration_run

                if livetime_type == 0:
                    continue

                livetime_type_years = livetime_type / (60 * 60 * 24 * 365.25)
                scale = 1 / livetime_type_years

                weights_one_year[is_type] = weights_one_year[is_type] * scale

    return weights_one_year

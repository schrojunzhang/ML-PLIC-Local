#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-05

import os
import shutil
from base_scripts.base_class.soft_base import soft


class scoring_function(soft):

    def __init__(self, job_id):
        super().__init__(job_id)
        self.sf2energy_term = {
            'affiscore': ['name', 'Orientation_Score', 'Score(Heavy_Ligand_Atoms)', 'Affinity_Score',
                          'Buried_Protein_Hydrophobic_Term', 'Hydrophobic_Complementarity_Term', 'Polar_Component_Term',
                          'Number_of_Protein-Ligand_Hydrophobic_Contacts', 'Number_of_Protein-Ligand_H-bonds',
                          'Number_of_Protein-Ligand_Salt-bridges', 'Number_of_Metal-Ligand_Bonds',
                          'Number_of_Interfacial_Unsatisfied_Polar_Atoms',
                          'Number_of_Interfacial_Unsatisfied_Charged_Atoms',
                          'Buried_Carbons'],
            'affinitydG': ['name', 'affinitydG', 'hb', 'ion', 'mlig', 'hh_hp_aa'],
            'alphaHB': ['name', 'AlphaHB', 'repulsion_term', 'attraction_term', 'HB_term'],
            'ase': ['name', 'ASE'],
            'asp': ['name', 'ASP.Fitness', 'ASP.ASP', 'ASP.Map', 'ASP.DEClash', 'ASP.DEInternal'],
            'autodock': ['name', 'AutoDock4.2Score', 'qq', 'hb', 'vdw', 'dsolv', 'tors'],
            'chemgauss': ['name', 'Chemgauss4 score', 'Steric', 'Clash', 'ProDesolv', 'LigDesolv', 'LigDesolvHB', 'HB'],
            'chemplp': ['name', 'PLP.Fitness', 'PLP.PLP', 'PLP.part.hbond', 'PLP.part.metal', 'PLP.part.buried',
                        'PLP.part.nonpolar', 'PLP.part.repulsive', 'PLP.ligand.clash', 'PLP.ligand.torsion',
                        'PLP.Chemscore.Hbond', 'Chemscore.CHOScore', 'PLP.Chemscore.Metal'],
            'chemscore': ['name', 'Chemscore.Fitness', 'Chemscore.DG', 'Chemscore.ZeroCoef', 'Chemscore.Hbond',
                          'Chemscore.Metal', 'Chemscore.Lipo', 'Chemscore.Internal_Hbond', 'Chemscore.DEClash',
                          'Chemscore.DEInternal', 'Chemscore.Rot'],
            'contact': ['name', 'Contact_Score'],
            'continuous': ['name', 'Continuous_Score', 'Continuous_vdw_energy', 'Continuous_es_energy',
                           'Continuous_internal_repulsive'],
            'cyscore': ['name', 'cyscore', 'cyscore_hydrophobic', 'cyscore_Vdw', 'cyscore_HBond', 'cyscore_Ent'],
            'dsx': ['name', 'atom_pairs', 'intra_clashes', 'sas_score'],
            'galaxy': ['name', 'total', 'qq_pl', 'desolv_pl', 'vdw_pl', 'hbond_pl', 'qq_l', 'desolv_l', 'vdw_l',
                       'hbond_l', 'DrugScore', 'HMscore', 'PLP_tor'],
            'GBVIWSAdG': ['name', 'GBVIWSA_dG', 'coul', 'sol', 'vdw', 'sa'],
            'goldscore': ['name', 'Goldscore.Fitness', 'Goldscore.Internal.Vdw', 'Goldscore.Internal.Torsion',
                          'Goldscore.External.Vdw', 'Goldscore.Internal.HBond', 'Goldscore.External.HBond'],
            'grid': ['name', 'Grid_Score', 'Grid_vdw_energy', 'Grid_es_energy', 'Grid_internal_repulsive'],
            'hawkins': ['name', 'Hawkins_GBSA_Score', 'Hawkins_GBSA_vdw_energy', 'Hawkins_GBSA_es_energy',
                        'Hawkins_GBSA_gb_energy', 'Hawkins_GBSA_sa_energy', 'Hawkins_internal_repulsive'],
            'LondondG': ['name', 'total', 'hbond', 'hbondA', 'metal', 'metal_r', 'ionP', 'ionN', 'aroalk', 'repulse',
                         'flex', 'p', 'l'],
            'nn_vina': ['nn_vina_affinity', 'nn_vina_gauss_1', 'nn_vina_gauss_2', 'nn_vina_repulsion',
                        'nn_vina_hydrophobic',
                        'nn_vina_hydrogen'],
            'nn_close': ['A_MN', 'OA_SA', 'HD_N', 'N_ZN', 'A_MG', 'HD_NA', 'A_CL', 'MG_OA', 'FE_HD', 'A_OA', 'NA_ZN',
                         'A_N', 'C_OA', 'F_HD', 'C_HD', 'NA_SA', 'A_ZN', 'C_NA', 'N_N', 'MN_N', 'F_N', 'FE_OA', 'HD_I',
                         'BR_C', 'MG_NA', 'C_ZN', 'CL_MG', 'BR_OA', 'A_FE', 'CL_OA', 'CL_N', 'NA_OA', 'F_ZN', 'HD_P',
                         'CL_ZN', 'C_C', 'C_CL', 'FE_N', 'HD_S', 'HD_MG', 'C_F', 'A_NA', 'BR_HD', 'HD_OA', 'HD_MN',
                         'A_SA', 'A_F', 'HD_SA', 'A_C', 'A_A', 'F_SA', 'C_N', 'HD_ZN', 'OA_OA', 'N_SA', 'CL_FE', 'C_MN',
                         'CL_HD', 'OA_ZN', 'MN_OA', 'C_MG', 'F_OA', 'CD_OA', 'S_ZN', 'N_OA', 'C_SA', 'N_NA', 'A_HD',
                         'HD_HD', 'SA_ZN'],
            'nn_semi': ['I_N', 'OA_SA', 'FE_NA', 'HD_NA', 'A_CL', 'MG_SA', 'A_CU', 'P_SA', 'C_NA', 'MN_NA', 'F_N',
                        'HD_N', 'HD_I', 'CL_MG', 'HD_S', 'CL_MN', 'F_OA', 'HD_OA', 'F_HD', 'A_SA', 'A_BR', 'BR_HD',
                        'SA_SA', 'A_MN', 'N_ZN', 'A_MG', 'I_OA', 'C_C', 'N_S', 'N_N', 'FE_N', 'NA_SA', 'BR_N', 'MN_N',
                        'A_P', 'BR_C', 'A_FE', 'MN_P', 'CL_OA', 'CU_HD', 'MN_S', 'A_S', 'FE_OA', 'NA_ZN', 'P_ZN', 'A_F',
                        'A_C', 'A_A', 'A_N', 'HD_MN', 'A_I', 'N_SA', 'C_OA', 'MG_P', 'BR_SA', 'CU_N', 'MN_OA', 'MG_N',
                        'HD_HD', 'C_FE', 'CL_NA', 'MG_OA', 'A_OA', 'CL_ZN', 'BR_OA', 'HD_ZN', 'HD_P', 'OA_P', 'OA_S',
                        'N_P', 'A_NA', 'CL_FE', 'HD_SA', 'C_MN', 'CL_HD', 'C_MG', 'FE_HD', 'MG_S', 'NA_S', 'NA_P',
                        'FE_SA', 'P_S', 'C_HD', 'A_ZN', 'CL_P', 'S_SA', 'CL_S', 'OA_ZN', 'N_NA', 'MN_SA', 'CL_N',
                        'NA_OA', 'C_ZN', 'C_CD', 'HD_MG', 'C_F', 'C_I', 'C_CL', 'C_N', 'C_P', 'C_S', 'A_HD', 'F_SA',
                        'MG_NA', 'OA_OA', 'CL_SA', 'S_ZN', 'N_OA', 'C_SA', 'SA_ZN'],
            'nn_atm': ['A', 'C', 'CL', 'I', 'N', 'P', 'S', 'BR', 'HD', 'NA', 'F', 'OA', 'SA'],
            'nn_elec': ['I_N', 'OA_SA', 'FE_NA', 'HD_NA', 'A_CL', 'MG_SA', 'P_SA', 'C_NA', 'MN_NA', 'F_N', 'HD_N',
                        'HD_I', 'CL_MG', 'HD_S', 'CL_MN', 'F_OA', 'HD_OA', 'F_HD', 'A_SA', 'A_BR', 'BR_HD', 'SA_SA',
                        'A_MN', 'N_ZN', 'A_MG', 'I_OA', 'C_C', 'N_S', 'N_N', 'FE_N', 'NA_SA', 'BR_N', 'MN_N', 'A_P',
                        'BR_C', 'A_FE', 'MN_P', 'CL_OA', 'CU_HD', 'MN_S', 'A_S', 'FE_OA', 'NA_ZN', 'P_ZN', 'A_F', 'A_C',
                        'A_A', 'A_N', 'HD_MN', 'A_I', 'N_SA', 'C_OA', 'MG_P', 'BR_SA', 'CU_N', 'MN_OA', 'MG_N', 'HD_HD',
                        'C_FE', 'CL_NA', 'MG_OA', 'A_OA', 'CL_ZN', 'BR_OA', 'HD_ZN', 'HD_P', 'OA_P', 'OA_S', 'N_P',
                        'A_NA', 'CL_FE', 'HD_SA', 'C_MN', 'CL_HD', 'C_MG', 'FE_HD', 'MG_S', 'NA_S', 'NA_P', 'FE_SA',
                        'P_S', 'C_HD', 'A_ZN', 'CL_P', 'S_SA', 'CL_S', 'OA_ZN', 'N_NA', 'MN_SA', 'CL_N', 'NA_OA',
                        'F_ZN',
                        'C_ZN', 'HD_MG', 'C_F', 'C_I', 'C_CL', 'C_N', 'C_P', 'C_S', 'A_HD', 'F_SA', 'MG_NA', 'OA_OA',
                        'CL_SA', 'S_ZN', 'N_OA', 'C_SA', 'SA_ZN'],
            'nn_rot': ['rot_bonds'],
            'nn_active_site_flexibility': ['SIDECHAIN_OTHER', 'SIDECHAIN_ALPHA', 'BACKBONE_ALPHA', 'SIDECHAIN_BETA',
                                           'BACKBONE_BETA', 'BACKBONE_OTHER'],
            'nn_hb': ['HDONOR-LIGAND_SIDECHAIN_BETA', 'HDONOR-LIGAND_BACKBONE_OTHER', 'HDONOR-LIGAND_SIDECHAIN_ALPHA',
                      'HDONOR-RECEPTOR_SIDECHAIN_OTHER', 'HDONOR-RECEPTOR_BACKBONE_ALPHA',
                      'HDONOR-RECEPTOR_SIDECHAIN_BETA',
                      'HDONOR-RECEPTOR_SIDECHAIN_ALPHA', 'HDONOR-LIGAND_SIDECHAIN_OTHER', 'HDONOR-LIGAND_BACKBONE_BETA',
                      'HDONOR-RECEPTOR_BACKBONE_BETA', 'HDONOR-RECEPTOR_BACKBONE_OTHER',
                      'HDONOR-LIGAND_BACKBONE_ALPHA'],
            'nn_hydrophobics': ['SIDECHAIN_OTHER', 'SIDECHAIN_ALPHA', 'BACKBONE_ALPHA', 'SIDECHAIN_BETA',
                                'BACKBONE_BETA',
                                'BACKBONE_OTHER'],
            'nn_stacking': ['ALPHA', 'BETA', 'OTHER'],
            'nn_pi_cation': ['LIGAND-CHARGED_BETA', 'LIGAND-CHARGED_ALPHA', 'RECEPTOR-CHARGED_BETA',
                             'RECEPTOR-CHARGED_OTHER',
                             'RECEPTOR-CHARGED_ALPHA', 'LIGAND-CHARGED_OTHER'],
            'nn_t_shape': ['ALPHA', 'BETA', 'OTHER'],
            'nn_salt_bridges': ['ALPHA', 'BETA', 'OTHER'],
            'plp': ['name', 'TOTAL_SCORE', 'SCORE_RB_PEN', 'SCORE_NORM_HEVATOMS', 'SCORE_NORM_CRT_HEVATOMS',
                    'SCORE_NORM_WEIGHT',
                    'SCORE_NORM_CRT_WEIGHT', 'SCORE_RB_PEN_NORM_CRT_HEVATOMS', 'SCORE_NORM_CONTACT', 'PLPtotal',
                    'PLPparthbond',
                    'PLPpartsteric', 'PLPpartmetal', 'PLPpartrepulsive', 'PLPpartburpolar', 'LIG_NUM_CLASH',
                    'LIG_NUM_CONTACT',
                    'LIG_NUM_NO_CONTACT', 'CHEMPLP_CLASH2', 'TRIPOS_TORS', 'ATOMS_OUTSIDE_BINDINGSITE'],
            'plp95': ['name', 'TOTAL_SCORE', 'SCORE_RB_PEN', 'SCORE_NORM_HEVATOMS', 'SCORE_NORM_CRT_HEVATOMS',
                      'SCORE_NORM_WEIGHT', 'SCORE_NORM_CRT_WEIGHT', 'SCORE_RB_PEN_NORM_CRT_HEVATOMS',
                      'SCORE_NORM_CONTACT',
                      'PLPparthbond', 'PLPpartSteric', 'PLPpartmetal', 'LIG_NUM_CLASH', 'LIG_NUM_CONTACT',
                      'LIG_NUM_NO_CONTACT',
                      'CHEMPLP_CLASH2', 'TRIPOS_TORS', 'ATOMS_OUTSIDE_BINDINGSITE'],
            'rdock': ['name', 'rdock', 'rdock.INTER', 'rdock.INTER.CONST', 'rdock.INTER.POLAR', 'rdock.INTER.REPUL',
                      'rdock.INTER.ROT', 'rdock.INTER.VDW', 'rdock.INTER.norm', 'rdock.INTRA', 'rdock.INTRA.DIHEDRAL',
                      'rdock.INTRA.DIHEDRAL.0', 'rdock.INTRA.POLAR', 'rdock.INTRA.POLAR.0', 'rdock.INTRA.REPUL',
                      'rdock.INTRA.REPUL.0', 'rdock.INTRA.VDW', 'rdock.INTRA.VDW.0', 'rdock.INTRA.norm', 'rdock.RESTR',
                      'rdock.RESTR.CAVITY', 'rdock.RESTR.norm', 'rdock.SYSTEM', 'rdock.SYSTEM.CONST',
                      'rdock.SYSTEM.DIHEDRAL',
                      'rdock.SYSTEM.norm', 'rdock.heavy', 'rdock.norm'],
            'rdock_sol': ['name', 'rdock_sol', 'rdock_sol.INTER', 'rdock_sol.INTER.CONST', 'rdock_sol.INTER.POLAR',
                          'rdock_sol.INTER.ROT',
                          'rdock_sol.INTER.SOLV', 'rdock_sol.INTER.VDW', 'rdock_sol.INTER.norm', 'rdock_sol.INTRA',
                          'rdock_sol.INTRA.DIHEDRAL',
                          'rdock_sol.INTRA.DIHEDRAL.0', 'rdock_sol.INTRA.POLAR', 'rdock_sol.INTRA.POLAR.0',
                          'rdock_sol.INTRA.REPUL',
                          'rdock_sol.INTRA.REPUL.0', 'rdock_sol.INTRA.SOLV', 'rdock_sol.INTRA.SOLV.lig_0',
                          'rdock_sol.INTRA.VDW',
                          'rdock_sol.INTRA.VDW.0', 'rdock_sol.INTRA.norm', 'rdock_sol.RESTR', 'rdock_sol.RESTR.CAVITY',
                          'rdock_sol.RESTR.norm',
                          'rdock_sol.SYSTEM', 'rdock_sol.SYSTEM.CONST', 'rdock_sol.SYSTEM.DIHEDRAL',
                          'rdock_sol.SYSTEM.SOLV', 'rdock_sol.SYSTEM.norm',
                          'rdock_sol.heavy', 'rdock_sol.norm'],
            'rfscore_credo': ['name', 'covalent-B0', 'covalent-B1', 'covalent-B2', 'covalent-B3', 'covalent-B4',
                              'covalent-B5',
                              'vdw_clash-B0', 'vdw_clash-B1', 'vdw_clash-B2', 'vdw_clash-B3', 'vdw_clash-B4',
                              'vdw_clash-B5',
                              'vdw-B0', 'vdw-B1', 'vdw-B2', 'vdw-B3', 'vdw-B4', 'vdw-B5', 'proximal-B0', 'proximal-B1',
                              'proximal-B2',
                              'proximal-B3', 'proximal-B4', 'proximal-B5', 'hbond-B0', 'hbond-B1', 'hbond-B2',
                              'hbond-B3',
                              'hbond-B4', 'hbond-B5', 'weak_hbond-B0', 'weak_hbond-B1', 'weak_hbond-B2',
                              'weak_hbond-B3',
                              'weak_hbond-B4', 'weak_hbond-B5', 'xbond-B0', 'xbond-B1', 'xbond-B2', 'xbond-B3',
                              'xbond-B4',
                              'xbond-B5', 'ionic-B0', 'ionic-B1', 'ionic-B2', 'ionic-B3', 'ionic-B4', 'ionic-B5',
                              'metal_complex-B0',
                              'metal_complex-B1', 'metal_complex-B2', 'metal_complex-B3', 'metal_complex-B4',
                              'metal_complex-B5',
                              'aromatic-B0', 'aromatic-B1', 'aromatic-B2', 'aromatic-B3', 'aromatic-B4', 'aromatic-B5',
                              'hydrophobic-B0', 'hydrophobic-B1', 'hydrophobic-B2', 'hydrophobic-B3', 'hydrophobic-B4',
                              'hydrophobic-B5',
                              'carbonyl-B0', 'carbonyl-B1', 'carbonyl-B2', 'carbonyl-B3', 'carbonyl-B4', 'carbonyl-B5'],
            'rfscore_element': ['name', 'C.C-B0', 'C.C-B1', 'C.C-B2', 'C.C-B3', 'C.C-B4', 'C.C-B5', 'C.N-B0', 'C.N-B1',
                                'C.N-B2', 'C.N-B3', 'C.N-B4', 'C.N-B5', 'C.O-B0', 'C.O-B1', 'C.O-B2', 'C.O-B3',
                                'C.O-B4', 'C.O-B5', 'C.F-B0', 'C.F-B1', 'C.F-B2', 'C.F-B3', 'C.F-B4', 'C.F-B5',
                                'C.P-B0', 'C.P-B1', 'C.P-B2', 'C.P-B3', 'C.P-B4', 'C.P-B5', 'C.S-B0', 'C.S-B1',
                                'C.S-B2', 'C.S-B3', 'C.S-B4', 'C.S-B5', 'C.Cl-B0', 'C.Cl-B1', 'C.Cl-B2', 'C.Cl-B3',
                                'C.Cl-B4', 'C.Cl-B5', 'C.Br-B0', 'C.Br-B1', 'C.Br-B2', 'C.Br-B3', 'C.Br-B4', 'C.Br-B5',
                                'C.I-B0', 'C.I-B1', 'C.I-B2', 'C.I-B3', 'C.I-B4', 'C.I-B5', 'N.C-B0', 'N.C-B1',
                                'N.C-B2', 'N.C-B3', 'N.C-B4', 'N.C-B5', 'N.N-B0', 'N.N-B1', 'N.N-B2', 'N.N-B3',
                                'N.N-B4', 'N.N-B5', 'N.O-B0', 'N.O-B1', 'N.O-B2', 'N.O-B3', 'N.O-B4', 'N.O-B5',
                                'N.F-B0', 'N.F-B1', 'N.F-B2', 'N.F-B3', 'N.F-B4', 'N.F-B5', 'N.P-B0', 'N.P-B1',
                                'N.P-B2', 'N.P-B3', 'N.P-B4', 'N.P-B5', 'N.S-B0', 'N.S-B1', 'N.S-B2', 'N.S-B3',
                                'N.S-B4', 'N.S-B5', 'N.Cl-B0', 'N.Cl-B1', 'N.Cl-B2', 'N.Cl-B3', 'N.Cl-B4', 'N.Cl-B5',
                                'N.Br-B0', 'N.Br-B1', 'N.Br-B2', 'N.Br-B3', 'N.Br-B4', 'N.Br-B5', 'N.I-B0', 'N.I-B1',
                                'N.I-B2', 'N.I-B3', 'N.I-B4', 'N.I-B5', 'O.C-B0', 'O.C-B1', 'O.C-B2', 'O.C-B3',
                                'O.C-B4', 'O.C-B5', 'O.N-B0', 'O.N-B1', 'O.N-B2', 'O.N-B3', 'O.N-B4', 'O.N-B5',
                                'O.O-B0', 'O.O-B1', 'O.O-B2', 'O.O-B3', 'O.O-B4', 'O.O-B5', 'O.F-B0', 'O.F-B1',
                                'O.F-B2', 'O.F-B3', 'O.F-B4', 'O.F-B5', 'O.P-B0', 'O.P-B1', 'O.P-B2', 'O.P-B3',
                                'O.P-B4', 'O.P-B5', 'O.S-B0', 'O.S-B1', 'O.S-B2', 'O.S-B3', 'O.S-B4', 'O.S-B5',
                                'O.Cl-B0', 'O.Cl-B1', 'O.Cl-B2', 'O.Cl-B3', 'O.Cl-B4', 'O.Cl-B5', 'O.Br-B0', 'O.Br-B1',
                                'O.Br-B2', 'O.Br-B3', 'O.Br-B4', 'O.Br-B5', 'O.I-B0', 'O.I-B1', 'O.I-B2', 'O.I-B3',
                                'O.I-B4', 'O.I-B5', 'F.C-B0', 'F.C-B1', 'F.C-B2', 'F.C-B3', 'F.C-B4', 'F.C-B5',
                                'F.N-B0', 'F.N-B1', 'F.N-B2', 'F.N-B3', 'F.N-B4', 'F.N-B5', 'F.O-B0', 'F.O-B1',
                                'F.O-B2', 'F.O-B3', 'F.O-B4', 'F.O-B5', 'F.F-B0', 'F.F-B1', 'F.F-B2', 'F.F-B3',
                                'F.F-B4', 'F.F-B5', 'F.P-B0', 'F.P-B1', 'F.P-B2', 'F.P-B3', 'F.P-B4', 'F.P-B5',
                                'F.S-B0', 'F.S-B1', 'F.S-B2', 'F.S-B3', 'F.S-B4', 'F.S-B5', 'F.Cl-B0', 'F.Cl-B1',
                                'F.Cl-B2', 'F.Cl-B3', 'F.Cl-B4', 'F.Cl-B5', 'F.Br-B0', 'F.Br-B1', 'F.Br-B2', 'F.Br-B3',
                                'F.Br-B4', 'F.Br-B5', 'F.I-B0', 'F.I-B1', 'F.I-B2', 'F.I-B3', 'F.I-B4', 'F.I-B5',
                                'P.C-B0', 'P.C-B1', 'P.C-B2', 'P.C-B3', 'P.C-B4', 'P.C-B5', 'P.N-B0', 'P.N-B1',
                                'P.N-B2', 'P.N-B3', 'P.N-B4', 'P.N-B5', 'P.O-B0', 'P.O-B1', 'P.O-B2', 'P.O-B3',
                                'P.O-B4', 'P.O-B5', 'P.F-B0', 'P.F-B1', 'P.F-B2', 'P.F-B3', 'P.F-B4', 'P.F-B5',
                                'P.P-B0', 'P.P-B1', 'P.P-B2', 'P.P-B3', 'P.P-B4', 'P.P-B5', 'P.S-B0', 'P.S-B1',
                                'P.S-B2', 'P.S-B3', 'P.S-B4', 'P.S-B5', 'P.Cl-B0', 'P.Cl-B1', 'P.Cl-B2', 'P.Cl-B3',
                                'P.Cl-B4', 'P.Cl-B5', 'P.Br-B0', 'P.Br-B1', 'P.Br-B2', 'P.Br-B3', 'P.Br-B4', 'P.Br-B5',
                                'P.I-B0', 'P.I-B1', 'P.I-B2', 'P.I-B3', 'P.I-B4', 'P.I-B5', 'S.C-B0', 'S.C-B1',
                                'S.C-B2', 'S.C-B3', 'S.C-B4', 'S.C-B5', 'S.N-B0', 'S.N-B1', 'S.N-B2', 'S.N-B3',
                                'S.N-B4', 'S.N-B5', 'S.O-B0', 'S.O-B1', 'S.O-B2', 'S.O-B3', 'S.O-B4', 'S.O-B5',
                                'S.F-B0', 'S.F-B1', 'S.F-B2', 'S.F-B3', 'S.F-B4', 'S.F-B5', 'S.P-B0', 'S.P-B1',
                                'S.P-B2', 'S.P-B3', 'S.P-B4', 'S.P-B5', 'S.S-B0', 'S.S-B1', 'S.S-B2', 'S.S-B3',
                                'S.S-B4', 'S.S-B5', 'S.Cl-B0', 'S.Cl-B1', 'S.Cl-B2', 'S.Cl-B3', 'S.Cl-B4', 'S.Cl-B5',
                                'S.Br-B0', 'S.Br-B1', 'S.Br-B2', 'S.Br-B3', 'S.Br-B4', 'S.Br-B5', 'S.I-B0', 'S.I-B1',
                                'S.I-B2', 'S.I-B3', 'S.I-B4', 'S.I-B5', 'Cl.C-B0', 'Cl.C-B1', 'Cl.C-B2', 'Cl.C-B3',
                                'Cl.C-B4', 'Cl.C-B5', 'Cl.N-B0', 'Cl.N-B1', 'Cl.N-B2', 'Cl.N-B3', 'Cl.N-B4', 'Cl.N-B5',
                                'Cl.O-B0', 'Cl.O-B1', 'Cl.O-B2', 'Cl.O-B3', 'Cl.O-B4', 'Cl.O-B5', 'Cl.F-B0', 'Cl.F-B1',
                                'Cl.F-B2', 'Cl.F-B3', 'Cl.F-B4', 'Cl.F-B5', 'Cl.P-B0', 'Cl.P-B1', 'Cl.P-B2', 'Cl.P-B3',
                                'Cl.P-B4', 'Cl.P-B5', 'Cl.S-B0', 'Cl.S-B1', 'Cl.S-B2', 'Cl.S-B3', 'Cl.S-B4', 'Cl.S-B5',
                                'Cl.Cl-B0', 'Cl.Cl-B1', 'Cl.Cl-B2', 'Cl.Cl-B3', 'Cl.Cl-B4', 'Cl.Cl-B5', 'Cl.Br-B0',
                                'Cl.Br-B1', 'Cl.Br-B2', 'Cl.Br-B3', 'Cl.Br-B4', 'Cl.Br-B5', 'Cl.I-B0', 'Cl.I-B1',
                                'Cl.I-B2', 'Cl.I-B3', 'Cl.I-B4', 'Cl.I-B5', 'Br.C-B0', 'Br.C-B1', 'Br.C-B2', 'Br.C-B3',
                                'Br.C-B4', 'Br.C-B5', 'Br.N-B0', 'Br.N-B1', 'Br.N-B2', 'Br.N-B3', 'Br.N-B4', 'Br.N-B5',
                                'Br.O-B0', 'Br.O-B1', 'Br.O-B2', 'Br.O-B3', 'Br.O-B4', 'Br.O-B5', 'Br.F-B0', 'Br.F-B1',
                                'Br.F-B2', 'Br.F-B3', 'Br.F-B4', 'Br.F-B5', 'Br.P-B0', 'Br.P-B1', 'Br.P-B2', 'Br.P-B3',
                                'Br.P-B4', 'Br.P-B5', 'Br.S-B0', 'Br.S-B1', 'Br.S-B2', 'Br.S-B3', 'Br.S-B4', 'Br.S-B5',
                                'Br.Cl-B0', 'Br.Cl-B1', 'Br.Cl-B2', 'Br.Cl-B3', 'Br.Cl-B4', 'Br.Cl-B5', 'Br.Br-B0',
                                'Br.Br-B1', 'Br.Br-B2', 'Br.Br-B3', 'Br.Br-B4', 'Br.Br-B5', 'Br.I-B0', 'Br.I-B1',
                                'Br.I-B2', 'Br.I-B3', 'Br.I-B4', 'Br.I-B5', 'I.C-B0', 'I.C-B1', 'I.C-B2', 'I.C-B3',
                                'I.C-B4', 'I.C-B5', 'I.N-B0', 'I.N-B1', 'I.N-B2', 'I.N-B3', 'I.N-B4', 'I.N-B5',
                                'I.O-B0', 'I.O-B1', 'I.O-B2', 'I.O-B3', 'I.O-B4', 'I.O-B5', 'I.F-B0', 'I.F-B1',
                                'I.F-B2', 'I.F-B3', 'I.F-B4', 'I.F-B5', 'I.P-B0', 'I.P-B1', 'I.P-B2', 'I.P-B3',
                                'I.P-B4', 'I.P-B5', 'I.S-B0', 'I.S-B1', 'I.S-B2', 'I.S-B3', 'I.S-B4', 'I.S-B5',
                                'I.Cl-B0', 'I.Cl-B1', 'I.Cl-B2', 'I.Cl-B3', 'I.Cl-B4', 'I.Cl-B5', 'I.Br-B0', 'I.Br-B1',
                                'I.Br-B2', 'I.Br-B3', 'I.Br-B4', 'I.Br-B5', 'I.I-B0', 'I.I-B1', 'I.I-B2', 'I.I-B3',
                                'I.I-B4', 'I.I-B5'],
            'rfscore_sybyl': ['name', 'Br.Br-B0', 'Br.Br-B1', 'Br.Br-B2', 'Br.Br-B3', 'Br.Br-B4',
                              'Br.Br-B5', 'Br.C+-B0', 'Br.C+-B1', 'Br.C+-B2', 'Br.C+-B3',
                              'Br.C+-B4', 'Br.C+-B5', 'Br.C1-B0', 'Br.C1-B1', 'Br.C1-B2',
                              'Br.C1-B3', 'Br.C1-B4', 'Br.C1-B5', 'Br.C2-B0', 'Br.C2-B1',
                              'Br.C2-B2', 'Br.C2-B3', 'Br.C2-B4', 'Br.C2-B5', 'Br.C3-B0',
                              'Br.C3-B1', 'Br.C3-B2', 'Br.C3-B3', 'Br.C3-B4', 'Br.C3-B5',
                              'Br.Cac-B0', 'Br.Cac-B1', 'Br.Cac-B2', 'Br.Cac-B3', 'Br.Cac-B4',
                              'Br.Cac-B5', 'Br.Car-B0', 'Br.Car-B1', 'Br.Car-B2', 'Br.Car-B3',
                              'Br.Car-B4', 'Br.Car-B5', 'Br.Cl-B0', 'Br.Cl-B1', 'Br.Cl-B2',
                              'Br.Cl-B3', 'Br.Cl-B4', 'Br.Cl-B5', 'Br.F-B0', 'Br.F-B1',
                              'Br.F-B2', 'Br.F-B3', 'Br.F-B4', 'Br.F-B5', 'Br.I-B0', 'Br.I-B1',
                              'Br.I-B2', 'Br.I-B3', 'Br.I-B4', 'Br.I-B5', 'Br.N1-B0', 'Br.N1-B1',
                              'Br.N1-B2', 'Br.N1-B3', 'Br.N1-B4', 'Br.N1-B5', 'Br.N2-B0',
                              'Br.N2-B1', 'Br.N2-B2', 'Br.N2-B3', 'Br.N2-B4', 'Br.N2-B5',
                              'Br.N3-B0', 'Br.N3-B1', 'Br.N3-B2', 'Br.N3-B3', 'Br.N3-B4',
                              'Br.N3-B5', 'Br.N3+-B0', 'Br.N3+-B1', 'Br.N3+-B2', 'Br.N3+-B3',
                              'Br.N3+-B4', 'Br.N3+-B5', 'Br.Nam-B0', 'Br.Nam-B1', 'Br.Nam-B2',
                              'Br.Nam-B3', 'Br.Nam-B4', 'Br.Nam-B5', 'Br.Nar-B0', 'Br.Nar-B1',
                              'Br.Nar-B2', 'Br.Nar-B3', 'Br.Nar-B4', 'Br.Nar-B5', 'Br.Ng+-B0',
                              'Br.Ng+-B1', 'Br.Ng+-B2', 'Br.Ng+-B3', 'Br.Ng+-B4', 'Br.Ng+-B5',
                              'Br.Npl-B0', 'Br.Npl-B1', 'Br.Npl-B2', 'Br.Npl-B3', 'Br.Npl-B4',
                              'Br.Npl-B5', 'Br.Ntr-B0', 'Br.Ntr-B1', 'Br.Ntr-B2', 'Br.Ntr-B3',
                              'Br.Ntr-B4', 'Br.Ntr-B5', 'Br.O--B0', 'Br.O--B1', 'Br.O--B2',
                              'Br.O--B3', 'Br.O--B4', 'Br.O--B5', 'Br.O.co2-B0', 'Br.O.co2-B1',
                              'Br.O.co2-B2', 'Br.O.co2-B3', 'Br.O.co2-B4', 'Br.O.co2-B5',
                              'Br.O2-B0', 'Br.O2-B1', 'Br.O2-B2', 'Br.O2-B3', 'Br.O2-B4',
                              'Br.O2-B5', 'Br.O3-B0', 'Br.O3-B1', 'Br.O3-B2', 'Br.O3-B3',
                              'Br.O3-B4', 'Br.O3-B5', 'Br.P-B0', 'Br.P-B1', 'Br.P-B2', 'Br.P-B3',
                              'Br.P-B4', 'Br.P-B5', 'Br.Pac-B0', 'Br.Pac-B1', 'Br.Pac-B2',
                              'Br.Pac-B3', 'Br.Pac-B4', 'Br.Pac-B5', 'Br.S2-B0', 'Br.S2-B1',
                              'Br.S2-B2', 'Br.S2-B3', 'Br.S2-B4', 'Br.S2-B5', 'Br.S3-B0',
                              'Br.S3-B1', 'Br.S3-B2', 'Br.S3-B3', 'Br.S3-B4', 'Br.S3-B5',
                              'Br.So2-B0', 'Br.So2-B1', 'Br.So2-B2', 'Br.So2-B3', 'Br.So2-B4',
                              'Br.So2-B5', 'Br.Sox-B0', 'Br.Sox-B1', 'Br.Sox-B2', 'Br.Sox-B3',
                              'Br.Sox-B4', 'Br.Sox-B5', 'C+.Br-B0', 'C+.Br-B1', 'C+.Br-B2',
                              'C+.Br-B3', 'C+.Br-B4', 'C+.Br-B5', 'C+.C+-B0', 'C+.C+-B1',
                              'C+.C+-B2', 'C+.C+-B3', 'C+.C+-B4', 'C+.C+-B5', 'C+.C1-B0',
                              'C+.C1-B1', 'C+.C1-B2', 'C+.C1-B3', 'C+.C1-B4', 'C+.C1-B5',
                              'C+.C2-B0', 'C+.C2-B1', 'C+.C2-B2', 'C+.C2-B3', 'C+.C2-B4',
                              'C+.C2-B5', 'C+.C3-B0', 'C+.C3-B1', 'C+.C3-B2', 'C+.C3-B3',
                              'C+.C3-B4', 'C+.C3-B5', 'C+.Cac-B0', 'C+.Cac-B1', 'C+.Cac-B2',
                              'C+.Cac-B3', 'C+.Cac-B4', 'C+.Cac-B5', 'C+.Car-B0', 'C+.Car-B1',
                              'C+.Car-B2', 'C+.Car-B3', 'C+.Car-B4', 'C+.Car-B5', 'C+.Cl-B0',
                              'C+.Cl-B1', 'C+.Cl-B2', 'C+.Cl-B3', 'C+.Cl-B4', 'C+.Cl-B5',
                              'C+.F-B0', 'C+.F-B1', 'C+.F-B2', 'C+.F-B3', 'C+.F-B4', 'C+.F-B5',
                              'C+.I-B0', 'C+.I-B1', 'C+.I-B2', 'C+.I-B3', 'C+.I-B4', 'C+.I-B5',
                              'C+.N1-B0', 'C+.N1-B1', 'C+.N1-B2', 'C+.N1-B3', 'C+.N1-B4',
                              'C+.N1-B5', 'C+.N2-B0', 'C+.N2-B1', 'C+.N2-B2', 'C+.N2-B3',
                              'C+.N2-B4', 'C+.N2-B5', 'C+.N3-B0', 'C+.N3-B1', 'C+.N3-B2',
                              'C+.N3-B3', 'C+.N3-B4', 'C+.N3-B5', 'C+.N3+-B0', 'C+.N3+-B1',
                              'C+.N3+-B2', 'C+.N3+-B3', 'C+.N3+-B4', 'C+.N3+-B5', 'C+.Nam-B0',
                              'C+.Nam-B1', 'C+.Nam-B2', 'C+.Nam-B3', 'C+.Nam-B4', 'C+.Nam-B5',
                              'C+.Nar-B0', 'C+.Nar-B1', 'C+.Nar-B2', 'C+.Nar-B3', 'C+.Nar-B4',
                              'C+.Nar-B5', 'C+.Ng+-B0', 'C+.Ng+-B1', 'C+.Ng+-B2', 'C+.Ng+-B3',
                              'C+.Ng+-B4', 'C+.Ng+-B5', 'C+.Npl-B0', 'C+.Npl-B1', 'C+.Npl-B2',
                              'C+.Npl-B3', 'C+.Npl-B4', 'C+.Npl-B5', 'C+.Ntr-B0', 'C+.Ntr-B1',
                              'C+.Ntr-B2', 'C+.Ntr-B3', 'C+.Ntr-B4', 'C+.Ntr-B5', 'C+.O--B0',
                              'C+.O--B1', 'C+.O--B2', 'C+.O--B3', 'C+.O--B4', 'C+.O--B5',
                              'C+.O.co2-B0', 'C+.O.co2-B1', 'C+.O.co2-B2', 'C+.O.co2-B3',
                              'C+.O.co2-B4', 'C+.O.co2-B5', 'C+.O2-B0', 'C+.O2-B1', 'C+.O2-B2',
                              'C+.O2-B3', 'C+.O2-B4', 'C+.O2-B5', 'C+.O3-B0', 'C+.O3-B1',
                              'C+.O3-B2', 'C+.O3-B3', 'C+.O3-B4', 'C+.O3-B5', 'C+.P-B0',
                              'C+.P-B1', 'C+.P-B2', 'C+.P-B3', 'C+.P-B4', 'C+.P-B5', 'C+.Pac-B0',
                              'C+.Pac-B1', 'C+.Pac-B2', 'C+.Pac-B3', 'C+.Pac-B4', 'C+.Pac-B5',
                              'C+.S2-B0', 'C+.S2-B1', 'C+.S2-B2', 'C+.S2-B3', 'C+.S2-B4',
                              'C+.S2-B5', 'C+.S3-B0', 'C+.S3-B1', 'C+.S3-B2', 'C+.S3-B3',
                              'C+.S3-B4', 'C+.S3-B5', 'C+.So2-B0', 'C+.So2-B1', 'C+.So2-B2',
                              'C+.So2-B3', 'C+.So2-B4', 'C+.So2-B5', 'C+.Sox-B0', 'C+.Sox-B1',
                              'C+.Sox-B2', 'C+.Sox-B3', 'C+.Sox-B4', 'C+.Sox-B5', 'C1.Br-B0',
                              'C1.Br-B1', 'C1.Br-B2', 'C1.Br-B3', 'C1.Br-B4', 'C1.Br-B5',
                              'C1.C+-B0', 'C1.C+-B1', 'C1.C+-B2', 'C1.C+-B3', 'C1.C+-B4',
                              'C1.C+-B5', 'C1.C1-B0', 'C1.C1-B1', 'C1.C1-B2', 'C1.C1-B3',
                              'C1.C1-B4', 'C1.C1-B5', 'C1.C2-B0', 'C1.C2-B1', 'C1.C2-B2',
                              'C1.C2-B3', 'C1.C2-B4', 'C1.C2-B5', 'C1.C3-B0', 'C1.C3-B1',
                              'C1.C3-B2', 'C1.C3-B3', 'C1.C3-B4', 'C1.C3-B5', 'C1.Cac-B0',
                              'C1.Cac-B1', 'C1.Cac-B2', 'C1.Cac-B3', 'C1.Cac-B4', 'C1.Cac-B5',
                              'C1.Car-B0', 'C1.Car-B1', 'C1.Car-B2', 'C1.Car-B3', 'C1.Car-B4',
                              'C1.Car-B5', 'C1.Cl-B0', 'C1.Cl-B1', 'C1.Cl-B2', 'C1.Cl-B3',
                              'C1.Cl-B4', 'C1.Cl-B5', 'C1.F-B0', 'C1.F-B1', 'C1.F-B2', 'C1.F-B3',
                              'C1.F-B4', 'C1.F-B5', 'C1.I-B0', 'C1.I-B1', 'C1.I-B2', 'C1.I-B3',
                              'C1.I-B4', 'C1.I-B5', 'C1.N1-B0', 'C1.N1-B1', 'C1.N1-B2',
                              'C1.N1-B3', 'C1.N1-B4', 'C1.N1-B5', 'C1.N2-B0', 'C1.N2-B1',
                              'C1.N2-B2', 'C1.N2-B3', 'C1.N2-B4', 'C1.N2-B5', 'C1.N3-B0',
                              'C1.N3-B1', 'C1.N3-B2', 'C1.N3-B3', 'C1.N3-B4', 'C1.N3-B5',
                              'C1.N3+-B0', 'C1.N3+-B1', 'C1.N3+-B2', 'C1.N3+-B3', 'C1.N3+-B4',
                              'C1.N3+-B5', 'C1.Nam-B0', 'C1.Nam-B1', 'C1.Nam-B2', 'C1.Nam-B3',
                              'C1.Nam-B4', 'C1.Nam-B5', 'C1.Nar-B0', 'C1.Nar-B1', 'C1.Nar-B2',
                              'C1.Nar-B3', 'C1.Nar-B4', 'C1.Nar-B5', 'C1.Ng+-B0', 'C1.Ng+-B1',
                              'C1.Ng+-B2', 'C1.Ng+-B3', 'C1.Ng+-B4', 'C1.Ng+-B5', 'C1.Npl-B0',
                              'C1.Npl-B1', 'C1.Npl-B2', 'C1.Npl-B3', 'C1.Npl-B4', 'C1.Npl-B5',
                              'C1.Ntr-B0', 'C1.Ntr-B1', 'C1.Ntr-B2', 'C1.Ntr-B3', 'C1.Ntr-B4',
                              'C1.Ntr-B5', 'C1.O--B0', 'C1.O--B1', 'C1.O--B2', 'C1.O--B3',
                              'C1.O--B4', 'C1.O--B5', 'C1.O.co2-B0', 'C1.O.co2-B1',
                              'C1.O.co2-B2', 'C1.O.co2-B3', 'C1.O.co2-B4', 'C1.O.co2-B5',
                              'C1.O2-B0', 'C1.O2-B1', 'C1.O2-B2', 'C1.O2-B3', 'C1.O2-B4',
                              'C1.O2-B5', 'C1.O3-B0', 'C1.O3-B1', 'C1.O3-B2', 'C1.O3-B3',
                              'C1.O3-B4', 'C1.O3-B5', 'C1.P-B0', 'C1.P-B1', 'C1.P-B2', 'C1.P-B3',
                              'C1.P-B4', 'C1.P-B5', 'C1.Pac-B0', 'C1.Pac-B1', 'C1.Pac-B2',
                              'C1.Pac-B3', 'C1.Pac-B4', 'C1.Pac-B5', 'C1.S2-B0', 'C1.S2-B1',
                              'C1.S2-B2', 'C1.S2-B3', 'C1.S2-B4', 'C1.S2-B5', 'C1.S3-B0',
                              'C1.S3-B1', 'C1.S3-B2', 'C1.S3-B3', 'C1.S3-B4', 'C1.S3-B5',
                              'C1.So2-B0', 'C1.So2-B1', 'C1.So2-B2', 'C1.So2-B3', 'C1.So2-B4',
                              'C1.So2-B5', 'C1.Sox-B0', 'C1.Sox-B1', 'C1.Sox-B2', 'C1.Sox-B3',
                              'C1.Sox-B4', 'C1.Sox-B5', 'C2.Br-B0', 'C2.Br-B1', 'C2.Br-B2',
                              'C2.Br-B3', 'C2.Br-B4', 'C2.Br-B5', 'C2.C+-B0', 'C2.C+-B1',
                              'C2.C+-B2', 'C2.C+-B3', 'C2.C+-B4', 'C2.C+-B5', 'C2.C1-B0',
                              'C2.C1-B1', 'C2.C1-B2', 'C2.C1-B3', 'C2.C1-B4', 'C2.C1-B5',
                              'C2.C2-B0', 'C2.C2-B1', 'C2.C2-B2', 'C2.C2-B3', 'C2.C2-B4',
                              'C2.C2-B5', 'C2.C3-B0', 'C2.C3-B1', 'C2.C3-B2', 'C2.C3-B3',
                              'C2.C3-B4', 'C2.C3-B5', 'C2.Cac-B0', 'C2.Cac-B1', 'C2.Cac-B2',
                              'C2.Cac-B3', 'C2.Cac-B4', 'C2.Cac-B5', 'C2.Car-B0', 'C2.Car-B1',
                              'C2.Car-B2', 'C2.Car-B3', 'C2.Car-B4', 'C2.Car-B5', 'C2.Cl-B0',
                              'C2.Cl-B1', 'C2.Cl-B2', 'C2.Cl-B3', 'C2.Cl-B4', 'C2.Cl-B5',
                              'C2.F-B0', 'C2.F-B1', 'C2.F-B2', 'C2.F-B3', 'C2.F-B4', 'C2.F-B5',
                              'C2.I-B0', 'C2.I-B1', 'C2.I-B2', 'C2.I-B3', 'C2.I-B4', 'C2.I-B5',
                              'C2.N1-B0', 'C2.N1-B1', 'C2.N1-B2', 'C2.N1-B3', 'C2.N1-B4',
                              'C2.N1-B5', 'C2.N2-B0', 'C2.N2-B1', 'C2.N2-B2', 'C2.N2-B3',
                              'C2.N2-B4', 'C2.N2-B5', 'C2.N3-B0', 'C2.N3-B1', 'C2.N3-B2',
                              'C2.N3-B3', 'C2.N3-B4', 'C2.N3-B5', 'C2.N3+-B0', 'C2.N3+-B1',
                              'C2.N3+-B2', 'C2.N3+-B3', 'C2.N3+-B4', 'C2.N3+-B5', 'C2.Nam-B0',
                              'C2.Nam-B1', 'C2.Nam-B2', 'C2.Nam-B3', 'C2.Nam-B4', 'C2.Nam-B5',
                              'C2.Nar-B0', 'C2.Nar-B1', 'C2.Nar-B2', 'C2.Nar-B3', 'C2.Nar-B4',
                              'C2.Nar-B5', 'C2.Ng+-B0', 'C2.Ng+-B1', 'C2.Ng+-B2', 'C2.Ng+-B3',
                              'C2.Ng+-B4', 'C2.Ng+-B5', 'C2.Npl-B0', 'C2.Npl-B1', 'C2.Npl-B2',
                              'C2.Npl-B3', 'C2.Npl-B4', 'C2.Npl-B5', 'C2.Ntr-B0', 'C2.Ntr-B1',
                              'C2.Ntr-B2', 'C2.Ntr-B3', 'C2.Ntr-B4', 'C2.Ntr-B5', 'C2.O--B0',
                              'C2.O--B1', 'C2.O--B2', 'C2.O--B3', 'C2.O--B4', 'C2.O--B5',
                              'C2.O.co2-B0', 'C2.O.co2-B1', 'C2.O.co2-B2', 'C2.O.co2-B3',
                              'C2.O.co2-B4', 'C2.O.co2-B5', 'C2.O2-B0', 'C2.O2-B1', 'C2.O2-B2',
                              'C2.O2-B3', 'C2.O2-B4', 'C2.O2-B5', 'C2.O3-B0', 'C2.O3-B1',
                              'C2.O3-B2', 'C2.O3-B3', 'C2.O3-B4', 'C2.O3-B5', 'C2.P-B0',
                              'C2.P-B1', 'C2.P-B2', 'C2.P-B3', 'C2.P-B4', 'C2.P-B5', 'C2.Pac-B0',
                              'C2.Pac-B1', 'C2.Pac-B2', 'C2.Pac-B3', 'C2.Pac-B4', 'C2.Pac-B5',
                              'C2.S2-B0', 'C2.S2-B1', 'C2.S2-B2', 'C2.S2-B3', 'C2.S2-B4',
                              'C2.S2-B5', 'C2.S3-B0', 'C2.S3-B1', 'C2.S3-B2', 'C2.S3-B3',
                              'C2.S3-B4', 'C2.S3-B5', 'C2.So2-B0', 'C2.So2-B1', 'C2.So2-B2',
                              'C2.So2-B3', 'C2.So2-B4', 'C2.So2-B5', 'C2.Sox-B0', 'C2.Sox-B1',
                              'C2.Sox-B2', 'C2.Sox-B3', 'C2.Sox-B4', 'C2.Sox-B5', 'C3.Br-B0',
                              'C3.Br-B1', 'C3.Br-B2', 'C3.Br-B3', 'C3.Br-B4', 'C3.Br-B5',
                              'C3.C+-B0', 'C3.C+-B1', 'C3.C+-B2', 'C3.C+-B3', 'C3.C+-B4',
                              'C3.C+-B5', 'C3.C1-B0', 'C3.C1-B1', 'C3.C1-B2', 'C3.C1-B3',
                              'C3.C1-B4', 'C3.C1-B5', 'C3.C2-B0', 'C3.C2-B1', 'C3.C2-B2',
                              'C3.C2-B3', 'C3.C2-B4', 'C3.C2-B5', 'C3.C3-B0', 'C3.C3-B1',
                              'C3.C3-B2', 'C3.C3-B3', 'C3.C3-B4', 'C3.C3-B5', 'C3.Cac-B0',
                              'C3.Cac-B1', 'C3.Cac-B2', 'C3.Cac-B3', 'C3.Cac-B4', 'C3.Cac-B5',
                              'C3.Car-B0', 'C3.Car-B1', 'C3.Car-B2', 'C3.Car-B3', 'C3.Car-B4',
                              'C3.Car-B5', 'C3.Cl-B0', 'C3.Cl-B1', 'C3.Cl-B2', 'C3.Cl-B3',
                              'C3.Cl-B4', 'C3.Cl-B5', 'C3.F-B0', 'C3.F-B1', 'C3.F-B2', 'C3.F-B3',
                              'C3.F-B4', 'C3.F-B5', 'C3.I-B0', 'C3.I-B1', 'C3.I-B2', 'C3.I-B3',
                              'C3.I-B4', 'C3.I-B5', 'C3.N1-B0', 'C3.N1-B1', 'C3.N1-B2',
                              'C3.N1-B3', 'C3.N1-B4', 'C3.N1-B5', 'C3.N2-B0', 'C3.N2-B1',
                              'C3.N2-B2', 'C3.N2-B3', 'C3.N2-B4', 'C3.N2-B5', 'C3.N3-B0',
                              'C3.N3-B1', 'C3.N3-B2', 'C3.N3-B3', 'C3.N3-B4', 'C3.N3-B5',
                              'C3.N3+-B0', 'C3.N3+-B1', 'C3.N3+-B2', 'C3.N3+-B3', 'C3.N3+-B4',
                              'C3.N3+-B5', 'C3.Nam-B0', 'C3.Nam-B1', 'C3.Nam-B2', 'C3.Nam-B3',
                              'C3.Nam-B4', 'C3.Nam-B5', 'C3.Nar-B0', 'C3.Nar-B1', 'C3.Nar-B2',
                              'C3.Nar-B3', 'C3.Nar-B4', 'C3.Nar-B5', 'C3.Ng+-B0', 'C3.Ng+-B1',
                              'C3.Ng+-B2', 'C3.Ng+-B3', 'C3.Ng+-B4', 'C3.Ng+-B5', 'C3.Npl-B0',
                              'C3.Npl-B1', 'C3.Npl-B2', 'C3.Npl-B3', 'C3.Npl-B4', 'C3.Npl-B5',
                              'C3.Ntr-B0', 'C3.Ntr-B1', 'C3.Ntr-B2', 'C3.Ntr-B3', 'C3.Ntr-B4',
                              'C3.Ntr-B5', 'C3.O--B0', 'C3.O--B1', 'C3.O--B2', 'C3.O--B3',
                              'C3.O--B4', 'C3.O--B5', 'C3.O.co2-B0', 'C3.O.co2-B1',
                              'C3.O.co2-B2', 'C3.O.co2-B3', 'C3.O.co2-B4', 'C3.O.co2-B5',
                              'C3.O2-B0', 'C3.O2-B1', 'C3.O2-B2', 'C3.O2-B3', 'C3.O2-B4',
                              'C3.O2-B5', 'C3.O3-B0', 'C3.O3-B1', 'C3.O3-B2', 'C3.O3-B3',
                              'C3.O3-B4', 'C3.O3-B5', 'C3.P-B0', 'C3.P-B1', 'C3.P-B2', 'C3.P-B3',
                              'C3.P-B4', 'C3.P-B5', 'C3.Pac-B0', 'C3.Pac-B1', 'C3.Pac-B2',
                              'C3.Pac-B3', 'C3.Pac-B4', 'C3.Pac-B5', 'C3.S2-B0', 'C3.S2-B1',
                              'C3.S2-B2', 'C3.S2-B3', 'C3.S2-B4', 'C3.S2-B5', 'C3.S3-B0',
                              'C3.S3-B1', 'C3.S3-B2', 'C3.S3-B3', 'C3.S3-B4', 'C3.S3-B5',
                              'C3.So2-B0', 'C3.So2-B1', 'C3.So2-B2', 'C3.So2-B3', 'C3.So2-B4',
                              'C3.So2-B5', 'C3.Sox-B0', 'C3.Sox-B1', 'C3.Sox-B2', 'C3.Sox-B3',
                              'C3.Sox-B4', 'C3.Sox-B5', 'Cac.Br-B0', 'Cac.Br-B1', 'Cac.Br-B2',
                              'Cac.Br-B3', 'Cac.Br-B4', 'Cac.Br-B5', 'Cac.C+-B0', 'Cac.C+-B1',
                              'Cac.C+-B2', 'Cac.C+-B3', 'Cac.C+-B4', 'Cac.C+-B5', 'Cac.C1-B0',
                              'Cac.C1-B1', 'Cac.C1-B2', 'Cac.C1-B3', 'Cac.C1-B4', 'Cac.C1-B5',
                              'Cac.C2-B0', 'Cac.C2-B1', 'Cac.C2-B2', 'Cac.C2-B3', 'Cac.C2-B4',
                              'Cac.C2-B5', 'Cac.C3-B0', 'Cac.C3-B1', 'Cac.C3-B2', 'Cac.C3-B3',
                              'Cac.C3-B4', 'Cac.C3-B5', 'Cac.Cac-B0', 'Cac.Cac-B1', 'Cac.Cac-B2',
                              'Cac.Cac-B3', 'Cac.Cac-B4', 'Cac.Cac-B5', 'Cac.Car-B0',
                              'Cac.Car-B1', 'Cac.Car-B2', 'Cac.Car-B3', 'Cac.Car-B4',
                              'Cac.Car-B5', 'Cac.Cl-B0', 'Cac.Cl-B1', 'Cac.Cl-B2', 'Cac.Cl-B3',
                              'Cac.Cl-B4', 'Cac.Cl-B5', 'Cac.F-B0', 'Cac.F-B1', 'Cac.F-B2',
                              'Cac.F-B3', 'Cac.F-B4', 'Cac.F-B5', 'Cac.I-B0', 'Cac.I-B1',
                              'Cac.I-B2', 'Cac.I-B3', 'Cac.I-B4', 'Cac.I-B5', 'Cac.N1-B0',
                              'Cac.N1-B1', 'Cac.N1-B2', 'Cac.N1-B3', 'Cac.N1-B4', 'Cac.N1-B5',
                              'Cac.N2-B0', 'Cac.N2-B1', 'Cac.N2-B2', 'Cac.N2-B3', 'Cac.N2-B4',
                              'Cac.N2-B5', 'Cac.N3-B0', 'Cac.N3-B1', 'Cac.N3-B2', 'Cac.N3-B3',
                              'Cac.N3-B4', 'Cac.N3-B5', 'Cac.N3+-B0', 'Cac.N3+-B1', 'Cac.N3+-B2',
                              'Cac.N3+-B3', 'Cac.N3+-B4', 'Cac.N3+-B5', 'Cac.Nam-B0',
                              'Cac.Nam-B1', 'Cac.Nam-B2', 'Cac.Nam-B3', 'Cac.Nam-B4',
                              'Cac.Nam-B5', 'Cac.Nar-B0', 'Cac.Nar-B1', 'Cac.Nar-B2',
                              'Cac.Nar-B3', 'Cac.Nar-B4', 'Cac.Nar-B5', 'Cac.Ng+-B0',
                              'Cac.Ng+-B1', 'Cac.Ng+-B2', 'Cac.Ng+-B3', 'Cac.Ng+-B4',
                              'Cac.Ng+-B5', 'Cac.Npl-B0', 'Cac.Npl-B1', 'Cac.Npl-B2',
                              'Cac.Npl-B3', 'Cac.Npl-B4', 'Cac.Npl-B5', 'Cac.Ntr-B0',
                              'Cac.Ntr-B1', 'Cac.Ntr-B2', 'Cac.Ntr-B3', 'Cac.Ntr-B4',
                              'Cac.Ntr-B5', 'Cac.O--B0', 'Cac.O--B1', 'Cac.O--B2', 'Cac.O--B3',
                              'Cac.O--B4', 'Cac.O--B5', 'Cac.O.co2-B0', 'Cac.O.co2-B1',
                              'Cac.O.co2-B2', 'Cac.O.co2-B3', 'Cac.O.co2-B4', 'Cac.O.co2-B5',
                              'Cac.O2-B0', 'Cac.O2-B1', 'Cac.O2-B2', 'Cac.O2-B3', 'Cac.O2-B4',
                              'Cac.O2-B5', 'Cac.O3-B0', 'Cac.O3-B1', 'Cac.O3-B2', 'Cac.O3-B3',
                              'Cac.O3-B4', 'Cac.O3-B5', 'Cac.P-B0', 'Cac.P-B1', 'Cac.P-B2',
                              'Cac.P-B3', 'Cac.P-B4', 'Cac.P-B5', 'Cac.Pac-B0', 'Cac.Pac-B1',
                              'Cac.Pac-B2', 'Cac.Pac-B3', 'Cac.Pac-B4', 'Cac.Pac-B5',
                              'Cac.S2-B0', 'Cac.S2-B1', 'Cac.S2-B2', 'Cac.S2-B3', 'Cac.S2-B4',
                              'Cac.S2-B5', 'Cac.S3-B0', 'Cac.S3-B1', 'Cac.S3-B2', 'Cac.S3-B3',
                              'Cac.S3-B4', 'Cac.S3-B5', 'Cac.So2-B0', 'Cac.So2-B1', 'Cac.So2-B2',
                              'Cac.So2-B3', 'Cac.So2-B4', 'Cac.So2-B5', 'Cac.Sox-B0',
                              'Cac.Sox-B1', 'Cac.Sox-B2', 'Cac.Sox-B3', 'Cac.Sox-B4',
                              'Cac.Sox-B5', 'Car.Br-B0', 'Car.Br-B1', 'Car.Br-B2', 'Car.Br-B3',
                              'Car.Br-B4', 'Car.Br-B5', 'Car.C+-B0', 'Car.C+-B1', 'Car.C+-B2',
                              'Car.C+-B3', 'Car.C+-B4', 'Car.C+-B5', 'Car.C1-B0', 'Car.C1-B1',
                              'Car.C1-B2', 'Car.C1-B3', 'Car.C1-B4', 'Car.C1-B5', 'Car.C2-B0',
                              'Car.C2-B1', 'Car.C2-B2', 'Car.C2-B3', 'Car.C2-B4', 'Car.C2-B5',
                              'Car.C3-B0', 'Car.C3-B1', 'Car.C3-B2', 'Car.C3-B3', 'Car.C3-B4',
                              'Car.C3-B5', 'Car.Cac-B0', 'Car.Cac-B1', 'Car.Cac-B2',
                              'Car.Cac-B3', 'Car.Cac-B4', 'Car.Cac-B5', 'Car.Car-B0',
                              'Car.Car-B1', 'Car.Car-B2', 'Car.Car-B3', 'Car.Car-B4',
                              'Car.Car-B5', 'Car.Cl-B0', 'Car.Cl-B1', 'Car.Cl-B2', 'Car.Cl-B3',
                              'Car.Cl-B4', 'Car.Cl-B5', 'Car.F-B0', 'Car.F-B1', 'Car.F-B2',
                              'Car.F-B3', 'Car.F-B4', 'Car.F-B5', 'Car.I-B0', 'Car.I-B1',
                              'Car.I-B2', 'Car.I-B3', 'Car.I-B4', 'Car.I-B5', 'Car.N1-B0',
                              'Car.N1-B1', 'Car.N1-B2', 'Car.N1-B3', 'Car.N1-B4', 'Car.N1-B5',
                              'Car.N2-B0', 'Car.N2-B1', 'Car.N2-B2', 'Car.N2-B3', 'Car.N2-B4',
                              'Car.N2-B5', 'Car.N3-B0', 'Car.N3-B1', 'Car.N3-B2', 'Car.N3-B3',
                              'Car.N3-B4', 'Car.N3-B5', 'Car.N3+-B0', 'Car.N3+-B1', 'Car.N3+-B2',
                              'Car.N3+-B3', 'Car.N3+-B4', 'Car.N3+-B5', 'Car.Nam-B0',
                              'Car.Nam-B1', 'Car.Nam-B2', 'Car.Nam-B3', 'Car.Nam-B4',
                              'Car.Nam-B5', 'Car.Nar-B0', 'Car.Nar-B1', 'Car.Nar-B2',
                              'Car.Nar-B3', 'Car.Nar-B4', 'Car.Nar-B5', 'Car.Ng+-B0',
                              'Car.Ng+-B1', 'Car.Ng+-B2', 'Car.Ng+-B3', 'Car.Ng+-B4',
                              'Car.Ng+-B5', 'Car.Npl-B0', 'Car.Npl-B1', 'Car.Npl-B2',
                              'Car.Npl-B3', 'Car.Npl-B4', 'Car.Npl-B5', 'Car.Ntr-B0',
                              'Car.Ntr-B1', 'Car.Ntr-B2', 'Car.Ntr-B3', 'Car.Ntr-B4',
                              'Car.Ntr-B5', 'Car.O--B0', 'Car.O--B1', 'Car.O--B2', 'Car.O--B3',
                              'Car.O--B4', 'Car.O--B5', 'Car.O.co2-B0', 'Car.O.co2-B1',
                              'Car.O.co2-B2', 'Car.O.co2-B3', 'Car.O.co2-B4', 'Car.O.co2-B5',
                              'Car.O2-B0', 'Car.O2-B1', 'Car.O2-B2', 'Car.O2-B3', 'Car.O2-B4',
                              'Car.O2-B5', 'Car.O3-B0', 'Car.O3-B1', 'Car.O3-B2', 'Car.O3-B3',
                              'Car.O3-B4', 'Car.O3-B5', 'Car.P-B0', 'Car.P-B1', 'Car.P-B2',
                              'Car.P-B3', 'Car.P-B4', 'Car.P-B5', 'Car.Pac-B0', 'Car.Pac-B1',
                              'Car.Pac-B2', 'Car.Pac-B3', 'Car.Pac-B4', 'Car.Pac-B5',
                              'Car.S2-B0', 'Car.S2-B1', 'Car.S2-B2', 'Car.S2-B3', 'Car.S2-B4',
                              'Car.S2-B5', 'Car.S3-B0', 'Car.S3-B1', 'Car.S3-B2', 'Car.S3-B3',
                              'Car.S3-B4', 'Car.S3-B5', 'Car.So2-B0', 'Car.So2-B1', 'Car.So2-B2',
                              'Car.So2-B3', 'Car.So2-B4', 'Car.So2-B5', 'Car.Sox-B0',
                              'Car.Sox-B1', 'Car.Sox-B2', 'Car.Sox-B3', 'Car.Sox-B4',
                              'Car.Sox-B5', 'Cl.Br-B0', 'Cl.Br-B1', 'Cl.Br-B2', 'Cl.Br-B3',
                              'Cl.Br-B4', 'Cl.Br-B5', 'Cl.C+-B0', 'Cl.C+-B1', 'Cl.C+-B2',
                              'Cl.C+-B3', 'Cl.C+-B4', 'Cl.C+-B5', 'Cl.C1-B0', 'Cl.C1-B1',
                              'Cl.C1-B2', 'Cl.C1-B3', 'Cl.C1-B4', 'Cl.C1-B5', 'Cl.C2-B0',
                              'Cl.C2-B1', 'Cl.C2-B2', 'Cl.C2-B3', 'Cl.C2-B4', 'Cl.C2-B5',
                              'Cl.C3-B0', 'Cl.C3-B1', 'Cl.C3-B2', 'Cl.C3-B3', 'Cl.C3-B4',
                              'Cl.C3-B5', 'Cl.Cac-B0', 'Cl.Cac-B1', 'Cl.Cac-B2', 'Cl.Cac-B3',
                              'Cl.Cac-B4', 'Cl.Cac-B5', 'Cl.Car-B0', 'Cl.Car-B1', 'Cl.Car-B2',
                              'Cl.Car-B3', 'Cl.Car-B4', 'Cl.Car-B5', 'Cl.Cl-B0', 'Cl.Cl-B1',
                              'Cl.Cl-B2', 'Cl.Cl-B3', 'Cl.Cl-B4', 'Cl.Cl-B5', 'Cl.F-B0',
                              'Cl.F-B1', 'Cl.F-B2', 'Cl.F-B3', 'Cl.F-B4', 'Cl.F-B5', 'Cl.I-B0',
                              'Cl.I-B1', 'Cl.I-B2', 'Cl.I-B3', 'Cl.I-B4', 'Cl.I-B5', 'Cl.N1-B0',
                              'Cl.N1-B1', 'Cl.N1-B2', 'Cl.N1-B3', 'Cl.N1-B4', 'Cl.N1-B5',
                              'Cl.N2-B0', 'Cl.N2-B1', 'Cl.N2-B2', 'Cl.N2-B3', 'Cl.N2-B4',
                              'Cl.N2-B5', 'Cl.N3-B0', 'Cl.N3-B1', 'Cl.N3-B2', 'Cl.N3-B3',
                              'Cl.N3-B4', 'Cl.N3-B5', 'Cl.N3+-B0', 'Cl.N3+-B1', 'Cl.N3+-B2',
                              'Cl.N3+-B3', 'Cl.N3+-B4', 'Cl.N3+-B5', 'Cl.Nam-B0', 'Cl.Nam-B1',
                              'Cl.Nam-B2', 'Cl.Nam-B3', 'Cl.Nam-B4', 'Cl.Nam-B5', 'Cl.Nar-B0',
                              'Cl.Nar-B1', 'Cl.Nar-B2', 'Cl.Nar-B3', 'Cl.Nar-B4', 'Cl.Nar-B5',
                              'Cl.Ng+-B0', 'Cl.Ng+-B1', 'Cl.Ng+-B2', 'Cl.Ng+-B3', 'Cl.Ng+-B4',
                              'Cl.Ng+-B5', 'Cl.Npl-B0', 'Cl.Npl-B1', 'Cl.Npl-B2', 'Cl.Npl-B3',
                              'Cl.Npl-B4', 'Cl.Npl-B5', 'Cl.Ntr-B0', 'Cl.Ntr-B1', 'Cl.Ntr-B2',
                              'Cl.Ntr-B3', 'Cl.Ntr-B4', 'Cl.Ntr-B5', 'Cl.O--B0', 'Cl.O--B1',
                              'Cl.O--B2', 'Cl.O--B3', 'Cl.O--B4', 'Cl.O--B5', 'Cl.O.co2-B0',
                              'Cl.O.co2-B1', 'Cl.O.co2-B2', 'Cl.O.co2-B3', 'Cl.O.co2-B4',
                              'Cl.O.co2-B5', 'Cl.O2-B0', 'Cl.O2-B1', 'Cl.O2-B2', 'Cl.O2-B3',
                              'Cl.O2-B4', 'Cl.O2-B5', 'Cl.O3-B0', 'Cl.O3-B1', 'Cl.O3-B2',
                              'Cl.O3-B3', 'Cl.O3-B4', 'Cl.O3-B5', 'Cl.P-B0', 'Cl.P-B1',
                              'Cl.P-B2', 'Cl.P-B3', 'Cl.P-B4', 'Cl.P-B5', 'Cl.Pac-B0',
                              'Cl.Pac-B1', 'Cl.Pac-B2', 'Cl.Pac-B3', 'Cl.Pac-B4', 'Cl.Pac-B5',
                              'Cl.S2-B0', 'Cl.S2-B1', 'Cl.S2-B2', 'Cl.S2-B3', 'Cl.S2-B4',
                              'Cl.S2-B5', 'Cl.S3-B0', 'Cl.S3-B1', 'Cl.S3-B2', 'Cl.S3-B3',
                              'Cl.S3-B4', 'Cl.S3-B5', 'Cl.So2-B0', 'Cl.So2-B1', 'Cl.So2-B2',
                              'Cl.So2-B3', 'Cl.So2-B4', 'Cl.So2-B5', 'Cl.Sox-B0', 'Cl.Sox-B1',
                              'Cl.Sox-B2', 'Cl.Sox-B3', 'Cl.Sox-B4', 'Cl.Sox-B5', 'F.Br-B0',
                              'F.Br-B1', 'F.Br-B2', 'F.Br-B3', 'F.Br-B4', 'F.Br-B5', 'F.C+-B0',
                              'F.C+-B1', 'F.C+-B2', 'F.C+-B3', 'F.C+-B4', 'F.C+-B5', 'F.C1-B0',
                              'F.C1-B1', 'F.C1-B2', 'F.C1-B3', 'F.C1-B4', 'F.C1-B5', 'F.C2-B0',
                              'F.C2-B1', 'F.C2-B2', 'F.C2-B3', 'F.C2-B4', 'F.C2-B5', 'F.C3-B0',
                              'F.C3-B1', 'F.C3-B2', 'F.C3-B3', 'F.C3-B4', 'F.C3-B5', 'F.Cac-B0',
                              'F.Cac-B1', 'F.Cac-B2', 'F.Cac-B3', 'F.Cac-B4', 'F.Cac-B5',
                              'F.Car-B0', 'F.Car-B1', 'F.Car-B2', 'F.Car-B3', 'F.Car-B4',
                              'F.Car-B5', 'F.Cl-B0', 'F.Cl-B1', 'F.Cl-B2', 'F.Cl-B3', 'F.Cl-B4',
                              'F.Cl-B5', 'F.F-B0', 'F.F-B1', 'F.F-B2', 'F.F-B3', 'F.F-B4',
                              'F.F-B5', 'F.I-B0', 'F.I-B1', 'F.I-B2', 'F.I-B3', 'F.I-B4',
                              'F.I-B5', 'F.N1-B0', 'F.N1-B1', 'F.N1-B2', 'F.N1-B3', 'F.N1-B4',
                              'F.N1-B5', 'F.N2-B0', 'F.N2-B1', 'F.N2-B2', 'F.N2-B3', 'F.N2-B4',
                              'F.N2-B5', 'F.N3-B0', 'F.N3-B1', 'F.N3-B2', 'F.N3-B3', 'F.N3-B4',
                              'F.N3-B5', 'F.N3+-B0', 'F.N3+-B1', 'F.N3+-B2', 'F.N3+-B3',
                              'F.N3+-B4', 'F.N3+-B5', 'F.Nam-B0', 'F.Nam-B1', 'F.Nam-B2',
                              'F.Nam-B3', 'F.Nam-B4', 'F.Nam-B5', 'F.Nar-B0', 'F.Nar-B1',
                              'F.Nar-B2', 'F.Nar-B3', 'F.Nar-B4', 'F.Nar-B5', 'F.Ng+-B0',
                              'F.Ng+-B1', 'F.Ng+-B2', 'F.Ng+-B3', 'F.Ng+-B4', 'F.Ng+-B5',
                              'F.Npl-B0', 'F.Npl-B1', 'F.Npl-B2', 'F.Npl-B3', 'F.Npl-B4',
                              'F.Npl-B5', 'F.Ntr-B0', 'F.Ntr-B1', 'F.Ntr-B2', 'F.Ntr-B3',
                              'F.Ntr-B4', 'F.Ntr-B5', 'F.O--B0', 'F.O--B1', 'F.O--B2', 'F.O--B3',
                              'F.O--B4', 'F.O--B5', 'F.O.co2-B0', 'F.O.co2-B1', 'F.O.co2-B2',
                              'F.O.co2-B3', 'F.O.co2-B4', 'F.O.co2-B5', 'F.O2-B0', 'F.O2-B1',
                              'F.O2-B2', 'F.O2-B3', 'F.O2-B4', 'F.O2-B5', 'F.O3-B0', 'F.O3-B1',
                              'F.O3-B2', 'F.O3-B3', 'F.O3-B4', 'F.O3-B5', 'F.P-B0', 'F.P-B1',
                              'F.P-B2', 'F.P-B3', 'F.P-B4', 'F.P-B5', 'F.Pac-B0', 'F.Pac-B1',
                              'F.Pac-B2', 'F.Pac-B3', 'F.Pac-B4', 'F.Pac-B5', 'F.S2-B0',
                              'F.S2-B1', 'F.S2-B2', 'F.S2-B3', 'F.S2-B4', 'F.S2-B5', 'F.S3-B0',
                              'F.S3-B1', 'F.S3-B2', 'F.S3-B3', 'F.S3-B4', 'F.S3-B5', 'F.So2-B0',
                              'F.So2-B1', 'F.So2-B2', 'F.So2-B3', 'F.So2-B4', 'F.So2-B5',
                              'F.Sox-B0', 'F.Sox-B1', 'F.Sox-B2', 'F.Sox-B3', 'F.Sox-B4',
                              'F.Sox-B5', 'I.Br-B0', 'I.Br-B1', 'I.Br-B2', 'I.Br-B3', 'I.Br-B4',
                              'I.Br-B5', 'I.C+-B0', 'I.C+-B1', 'I.C+-B2', 'I.C+-B3', 'I.C+-B4',
                              'I.C+-B5', 'I.C1-B0', 'I.C1-B1', 'I.C1-B2', 'I.C1-B3', 'I.C1-B4',
                              'I.C1-B5', 'I.C2-B0', 'I.C2-B1', 'I.C2-B2', 'I.C2-B3', 'I.C2-B4',
                              'I.C2-B5', 'I.C3-B0', 'I.C3-B1', 'I.C3-B2', 'I.C3-B3', 'I.C3-B4',
                              'I.C3-B5', 'I.Cac-B0', 'I.Cac-B1', 'I.Cac-B2', 'I.Cac-B3',
                              'I.Cac-B4', 'I.Cac-B5', 'I.Car-B0', 'I.Car-B1', 'I.Car-B2',
                              'I.Car-B3', 'I.Car-B4', 'I.Car-B5', 'I.Cl-B0', 'I.Cl-B1',
                              'I.Cl-B2', 'I.Cl-B3', 'I.Cl-B4', 'I.Cl-B5', 'I.F-B0', 'I.F-B1',
                              'I.F-B2', 'I.F-B3', 'I.F-B4', 'I.F-B5', 'I.I-B0', 'I.I-B1',
                              'I.I-B2', 'I.I-B3', 'I.I-B4', 'I.I-B5', 'I.N1-B0', 'I.N1-B1',
                              'I.N1-B2', 'I.N1-B3', 'I.N1-B4', 'I.N1-B5', 'I.N2-B0', 'I.N2-B1',
                              'I.N2-B2', 'I.N2-B3', 'I.N2-B4', 'I.N2-B5', 'I.N3-B0', 'I.N3-B1',
                              'I.N3-B2', 'I.N3-B3', 'I.N3-B4', 'I.N3-B5', 'I.N3+-B0', 'I.N3+-B1',
                              'I.N3+-B2', 'I.N3+-B3', 'I.N3+-B4', 'I.N3+-B5', 'I.Nam-B0',
                              'I.Nam-B1', 'I.Nam-B2', 'I.Nam-B3', 'I.Nam-B4', 'I.Nam-B5',
                              'I.Nar-B0', 'I.Nar-B1', 'I.Nar-B2', 'I.Nar-B3', 'I.Nar-B4',
                              'I.Nar-B5', 'I.Ng+-B0', 'I.Ng+-B1', 'I.Ng+-B2', 'I.Ng+-B3',
                              'I.Ng+-B4', 'I.Ng+-B5', 'I.Npl-B0', 'I.Npl-B1', 'I.Npl-B2',
                              'I.Npl-B3', 'I.Npl-B4', 'I.Npl-B5', 'I.Ntr-B0', 'I.Ntr-B1',
                              'I.Ntr-B2', 'I.Ntr-B3', 'I.Ntr-B4', 'I.Ntr-B5', 'I.O--B0',
                              'I.O--B1', 'I.O--B2', 'I.O--B3', 'I.O--B4', 'I.O--B5',
                              'I.O.co2-B0', 'I.O.co2-B1', 'I.O.co2-B2', 'I.O.co2-B3',
                              'I.O.co2-B4', 'I.O.co2-B5', 'I.O2-B0', 'I.O2-B1', 'I.O2-B2',
                              'I.O2-B3', 'I.O2-B4', 'I.O2-B5', 'I.O3-B0', 'I.O3-B1', 'I.O3-B2',
                              'I.O3-B3', 'I.O3-B4', 'I.O3-B5', 'I.P-B0', 'I.P-B1', 'I.P-B2',
                              'I.P-B3', 'I.P-B4', 'I.P-B5', 'I.Pac-B0', 'I.Pac-B1', 'I.Pac-B2',
                              'I.Pac-B3', 'I.Pac-B4', 'I.Pac-B5', 'I.S2-B0', 'I.S2-B1',
                              'I.S2-B2', 'I.S2-B3', 'I.S2-B4', 'I.S2-B5', 'I.S3-B0', 'I.S3-B1',
                              'I.S3-B2', 'I.S3-B3', 'I.S3-B4', 'I.S3-B5', 'I.So2-B0', 'I.So2-B1',
                              'I.So2-B2', 'I.So2-B3', 'I.So2-B4', 'I.So2-B5', 'I.Sox-B0',
                              'I.Sox-B1', 'I.Sox-B2', 'I.Sox-B3', 'I.Sox-B4', 'I.Sox-B5',
                              'N1.Br-B0', 'N1.Br-B1', 'N1.Br-B2', 'N1.Br-B3', 'N1.Br-B4',
                              'N1.Br-B5', 'N1.C+-B0', 'N1.C+-B1', 'N1.C+-B2', 'N1.C+-B3',
                              'N1.C+-B4', 'N1.C+-B5', 'N1.C1-B0', 'N1.C1-B1', 'N1.C1-B2',
                              'N1.C1-B3', 'N1.C1-B4', 'N1.C1-B5', 'N1.C2-B0', 'N1.C2-B1',
                              'N1.C2-B2', 'N1.C2-B3', 'N1.C2-B4', 'N1.C2-B5', 'N1.C3-B0',
                              'N1.C3-B1', 'N1.C3-B2', 'N1.C3-B3', 'N1.C3-B4', 'N1.C3-B5',
                              'N1.Cac-B0', 'N1.Cac-B1', 'N1.Cac-B2', 'N1.Cac-B3', 'N1.Cac-B4',
                              'N1.Cac-B5', 'N1.Car-B0', 'N1.Car-B1', 'N1.Car-B2', 'N1.Car-B3',
                              'N1.Car-B4', 'N1.Car-B5', 'N1.Cl-B0', 'N1.Cl-B1', 'N1.Cl-B2',
                              'N1.Cl-B3', 'N1.Cl-B4', 'N1.Cl-B5', 'N1.F-B0', 'N1.F-B1',
                              'N1.F-B2', 'N1.F-B3', 'N1.F-B4', 'N1.F-B5', 'N1.I-B0', 'N1.I-B1',
                              'N1.I-B2', 'N1.I-B3', 'N1.I-B4', 'N1.I-B5', 'N1.N1-B0', 'N1.N1-B1',
                              'N1.N1-B2', 'N1.N1-B3', 'N1.N1-B4', 'N1.N1-B5', 'N1.N2-B0',
                              'N1.N2-B1', 'N1.N2-B2', 'N1.N2-B3', 'N1.N2-B4', 'N1.N2-B5',
                              'N1.N3-B0', 'N1.N3-B1', 'N1.N3-B2', 'N1.N3-B3', 'N1.N3-B4',
                              'N1.N3-B5', 'N1.N3+-B0', 'N1.N3+-B1', 'N1.N3+-B2', 'N1.N3+-B3',
                              'N1.N3+-B4', 'N1.N3+-B5', 'N1.Nam-B0', 'N1.Nam-B1', 'N1.Nam-B2',
                              'N1.Nam-B3', 'N1.Nam-B4', 'N1.Nam-B5', 'N1.Nar-B0', 'N1.Nar-B1',
                              'N1.Nar-B2', 'N1.Nar-B3', 'N1.Nar-B4', 'N1.Nar-B5', 'N1.Ng+-B0',
                              'N1.Ng+-B1', 'N1.Ng+-B2', 'N1.Ng+-B3', 'N1.Ng+-B4', 'N1.Ng+-B5',
                              'N1.Npl-B0', 'N1.Npl-B1', 'N1.Npl-B2', 'N1.Npl-B3', 'N1.Npl-B4',
                              'N1.Npl-B5', 'N1.Ntr-B0', 'N1.Ntr-B1', 'N1.Ntr-B2', 'N1.Ntr-B3',
                              'N1.Ntr-B4', 'N1.Ntr-B5', 'N1.O--B0', 'N1.O--B1', 'N1.O--B2',
                              'N1.O--B3', 'N1.O--B4', 'N1.O--B5', 'N1.O.co2-B0', 'N1.O.co2-B1',
                              'N1.O.co2-B2', 'N1.O.co2-B3', 'N1.O.co2-B4', 'N1.O.co2-B5',
                              'N1.O2-B0', 'N1.O2-B1', 'N1.O2-B2', 'N1.O2-B3', 'N1.O2-B4',
                              'N1.O2-B5', 'N1.O3-B0', 'N1.O3-B1', 'N1.O3-B2', 'N1.O3-B3',
                              'N1.O3-B4', 'N1.O3-B5', 'N1.P-B0', 'N1.P-B1', 'N1.P-B2', 'N1.P-B3',
                              'N1.P-B4', 'N1.P-B5', 'N1.Pac-B0', 'N1.Pac-B1', 'N1.Pac-B2',
                              'N1.Pac-B3', 'N1.Pac-B4', 'N1.Pac-B5', 'N1.S2-B0', 'N1.S2-B1',
                              'N1.S2-B2', 'N1.S2-B3', 'N1.S2-B4', 'N1.S2-B5', 'N1.S3-B0',
                              'N1.S3-B1', 'N1.S3-B2', 'N1.S3-B3', 'N1.S3-B4', 'N1.S3-B5',
                              'N1.So2-B0', 'N1.So2-B1', 'N1.So2-B2', 'N1.So2-B3', 'N1.So2-B4',
                              'N1.So2-B5', 'N1.Sox-B0', 'N1.Sox-B1', 'N1.Sox-B2', 'N1.Sox-B3',
                              'N1.Sox-B4', 'N1.Sox-B5', 'N2.Br-B0', 'N2.Br-B1', 'N2.Br-B2',
                              'N2.Br-B3', 'N2.Br-B4', 'N2.Br-B5', 'N2.C+-B0', 'N2.C+-B1',
                              'N2.C+-B2', 'N2.C+-B3', 'N2.C+-B4', 'N2.C+-B5', 'N2.C1-B0',
                              'N2.C1-B1', 'N2.C1-B2', 'N2.C1-B3', 'N2.C1-B4', 'N2.C1-B5',
                              'N2.C2-B0', 'N2.C2-B1', 'N2.C2-B2', 'N2.C2-B3', 'N2.C2-B4',
                              'N2.C2-B5', 'N2.C3-B0', 'N2.C3-B1', 'N2.C3-B2', 'N2.C3-B3',
                              'N2.C3-B4', 'N2.C3-B5', 'N2.Cac-B0', 'N2.Cac-B1', 'N2.Cac-B2',
                              'N2.Cac-B3', 'N2.Cac-B4', 'N2.Cac-B5', 'N2.Car-B0', 'N2.Car-B1',
                              'N2.Car-B2', 'N2.Car-B3', 'N2.Car-B4', 'N2.Car-B5', 'N2.Cl-B0',
                              'N2.Cl-B1', 'N2.Cl-B2', 'N2.Cl-B3', 'N2.Cl-B4', 'N2.Cl-B5',
                              'N2.F-B0', 'N2.F-B1', 'N2.F-B2', 'N2.F-B3', 'N2.F-B4', 'N2.F-B5',
                              'N2.I-B0', 'N2.I-B1', 'N2.I-B2', 'N2.I-B3', 'N2.I-B4', 'N2.I-B5',
                              'N2.N1-B0', 'N2.N1-B1', 'N2.N1-B2', 'N2.N1-B3', 'N2.N1-B4',
                              'N2.N1-B5', 'N2.N2-B0', 'N2.N2-B1', 'N2.N2-B2', 'N2.N2-B3',
                              'N2.N2-B4', 'N2.N2-B5', 'N2.N3-B0', 'N2.N3-B1', 'N2.N3-B2',
                              'N2.N3-B3', 'N2.N3-B4', 'N2.N3-B5', 'N2.N3+-B0', 'N2.N3+-B1',
                              'N2.N3+-B2', 'N2.N3+-B3', 'N2.N3+-B4', 'N2.N3+-B5', 'N2.Nam-B0',
                              'N2.Nam-B1', 'N2.Nam-B2', 'N2.Nam-B3', 'N2.Nam-B4', 'N2.Nam-B5',
                              'N2.Nar-B0', 'N2.Nar-B1', 'N2.Nar-B2', 'N2.Nar-B3', 'N2.Nar-B4',
                              'N2.Nar-B5', 'N2.Ng+-B0', 'N2.Ng+-B1', 'N2.Ng+-B2', 'N2.Ng+-B3',
                              'N2.Ng+-B4', 'N2.Ng+-B5', 'N2.Npl-B0', 'N2.Npl-B1', 'N2.Npl-B2',
                              'N2.Npl-B3', 'N2.Npl-B4', 'N2.Npl-B5', 'N2.Ntr-B0', 'N2.Ntr-B1',
                              'N2.Ntr-B2', 'N2.Ntr-B3', 'N2.Ntr-B4', 'N2.Ntr-B5', 'N2.O--B0',
                              'N2.O--B1', 'N2.O--B2', 'N2.O--B3', 'N2.O--B4', 'N2.O--B5',
                              'N2.O.co2-B0', 'N2.O.co2-B1', 'N2.O.co2-B2', 'N2.O.co2-B3',
                              'N2.O.co2-B4', 'N2.O.co2-B5', 'N2.O2-B0', 'N2.O2-B1', 'N2.O2-B2',
                              'N2.O2-B3', 'N2.O2-B4', 'N2.O2-B5', 'N2.O3-B0', 'N2.O3-B1',
                              'N2.O3-B2', 'N2.O3-B3', 'N2.O3-B4', 'N2.O3-B5', 'N2.P-B0',
                              'N2.P-B1', 'N2.P-B2', 'N2.P-B3', 'N2.P-B4', 'N2.P-B5', 'N2.Pac-B0',
                              'N2.Pac-B1', 'N2.Pac-B2', 'N2.Pac-B3', 'N2.Pac-B4', 'N2.Pac-B5',
                              'N2.S2-B0', 'N2.S2-B1', 'N2.S2-B2', 'N2.S2-B3', 'N2.S2-B4',
                              'N2.S2-B5', 'N2.S3-B0', 'N2.S3-B1', 'N2.S3-B2', 'N2.S3-B3',
                              'N2.S3-B4', 'N2.S3-B5', 'N2.So2-B0', 'N2.So2-B1', 'N2.So2-B2',
                              'N2.So2-B3', 'N2.So2-B4', 'N2.So2-B5', 'N2.Sox-B0', 'N2.Sox-B1',
                              'N2.Sox-B2', 'N2.Sox-B3', 'N2.Sox-B4', 'N2.Sox-B5', 'N3.Br-B0',
                              'N3.Br-B1', 'N3.Br-B2', 'N3.Br-B3', 'N3.Br-B4', 'N3.Br-B5',
                              'N3.C+-B0', 'N3.C+-B1', 'N3.C+-B2', 'N3.C+-B3', 'N3.C+-B4',
                              'N3.C+-B5', 'N3.C1-B0', 'N3.C1-B1', 'N3.C1-B2', 'N3.C1-B3',
                              'N3.C1-B4', 'N3.C1-B5', 'N3.C2-B0', 'N3.C2-B1', 'N3.C2-B2',
                              'N3.C2-B3', 'N3.C2-B4', 'N3.C2-B5', 'N3.C3-B0', 'N3.C3-B1',
                              'N3.C3-B2', 'N3.C3-B3', 'N3.C3-B4', 'N3.C3-B5', 'N3.Cac-B0',
                              'N3.Cac-B1', 'N3.Cac-B2', 'N3.Cac-B3', 'N3.Cac-B4', 'N3.Cac-B5',
                              'N3.Car-B0', 'N3.Car-B1', 'N3.Car-B2', 'N3.Car-B3', 'N3.Car-B4',
                              'N3.Car-B5', 'N3.Cl-B0', 'N3.Cl-B1', 'N3.Cl-B2', 'N3.Cl-B3',
                              'N3.Cl-B4', 'N3.Cl-B5', 'N3.F-B0', 'N3.F-B1', 'N3.F-B2', 'N3.F-B3',
                              'N3.F-B4', 'N3.F-B5', 'N3.I-B0', 'N3.I-B1', 'N3.I-B2', 'N3.I-B3',
                              'N3.I-B4', 'N3.I-B5', 'N3.N1-B0', 'N3.N1-B1', 'N3.N1-B2',
                              'N3.N1-B3', 'N3.N1-B4', 'N3.N1-B5', 'N3.N2-B0', 'N3.N2-B1',
                              'N3.N2-B2', 'N3.N2-B3', 'N3.N2-B4', 'N3.N2-B5', 'N3.N3-B0',
                              'N3.N3-B1', 'N3.N3-B2', 'N3.N3-B3', 'N3.N3-B4', 'N3.N3-B5',
                              'N3.N3+-B0', 'N3.N3+-B1', 'N3.N3+-B2', 'N3.N3+-B3', 'N3.N3+-B4',
                              'N3.N3+-B5', 'N3.Nam-B0', 'N3.Nam-B1', 'N3.Nam-B2', 'N3.Nam-B3',
                              'N3.Nam-B4', 'N3.Nam-B5', 'N3.Nar-B0', 'N3.Nar-B1', 'N3.Nar-B2',
                              'N3.Nar-B3', 'N3.Nar-B4', 'N3.Nar-B5', 'N3.Ng+-B0', 'N3.Ng+-B1',
                              'N3.Ng+-B2', 'N3.Ng+-B3', 'N3.Ng+-B4', 'N3.Ng+-B5', 'N3.Npl-B0',
                              'N3.Npl-B1', 'N3.Npl-B2', 'N3.Npl-B3', 'N3.Npl-B4', 'N3.Npl-B5',
                              'N3.Ntr-B0', 'N3.Ntr-B1', 'N3.Ntr-B2', 'N3.Ntr-B3', 'N3.Ntr-B4',
                              'N3.Ntr-B5', 'N3.O--B0', 'N3.O--B1', 'N3.O--B2', 'N3.O--B3',
                              'N3.O--B4', 'N3.O--B5', 'N3.O.co2-B0', 'N3.O.co2-B1',
                              'N3.O.co2-B2', 'N3.O.co2-B3', 'N3.O.co2-B4', 'N3.O.co2-B5',
                              'N3.O2-B0', 'N3.O2-B1', 'N3.O2-B2', 'N3.O2-B3', 'N3.O2-B4',
                              'N3.O2-B5', 'N3.O3-B0', 'N3.O3-B1', 'N3.O3-B2', 'N3.O3-B3',
                              'N3.O3-B4', 'N3.O3-B5', 'N3.P-B0', 'N3.P-B1', 'N3.P-B2', 'N3.P-B3',
                              'N3.P-B4', 'N3.P-B5', 'N3.Pac-B0', 'N3.Pac-B1', 'N3.Pac-B2',
                              'N3.Pac-B3', 'N3.Pac-B4', 'N3.Pac-B5', 'N3.S2-B0', 'N3.S2-B1',
                              'N3.S2-B2', 'N3.S2-B3', 'N3.S2-B4', 'N3.S2-B5', 'N3.S3-B0',
                              'N3.S3-B1', 'N3.S3-B2', 'N3.S3-B3', 'N3.S3-B4', 'N3.S3-B5',
                              'N3.So2-B0', 'N3.So2-B1', 'N3.So2-B2', 'N3.So2-B3', 'N3.So2-B4',
                              'N3.So2-B5', 'N3.Sox-B0', 'N3.Sox-B1', 'N3.Sox-B2', 'N3.Sox-B3',
                              'N3.Sox-B4', 'N3.Sox-B5', 'N3+.Br-B0', 'N3+.Br-B1', 'N3+.Br-B2',
                              'N3+.Br-B3', 'N3+.Br-B4', 'N3+.Br-B5', 'N3+.C+-B0', 'N3+.C+-B1',
                              'N3+.C+-B2', 'N3+.C+-B3', 'N3+.C+-B4', 'N3+.C+-B5', 'N3+.C1-B0',
                              'N3+.C1-B1', 'N3+.C1-B2', 'N3+.C1-B3', 'N3+.C1-B4', 'N3+.C1-B5',
                              'N3+.C2-B0', 'N3+.C2-B1', 'N3+.C2-B2', 'N3+.C2-B3', 'N3+.C2-B4',
                              'N3+.C2-B5', 'N3+.C3-B0', 'N3+.C3-B1', 'N3+.C3-B2', 'N3+.C3-B3',
                              'N3+.C3-B4', 'N3+.C3-B5', 'N3+.Cac-B0', 'N3+.Cac-B1', 'N3+.Cac-B2',
                              'N3+.Cac-B3', 'N3+.Cac-B4', 'N3+.Cac-B5', 'N3+.Car-B0',
                              'N3+.Car-B1', 'N3+.Car-B2', 'N3+.Car-B3', 'N3+.Car-B4',
                              'N3+.Car-B5', 'N3+.Cl-B0', 'N3+.Cl-B1', 'N3+.Cl-B2', 'N3+.Cl-B3',
                              'N3+.Cl-B4', 'N3+.Cl-B5', 'N3+.F-B0', 'N3+.F-B1', 'N3+.F-B2',
                              'N3+.F-B3', 'N3+.F-B4', 'N3+.F-B5', 'N3+.I-B0', 'N3+.I-B1',
                              'N3+.I-B2', 'N3+.I-B3', 'N3+.I-B4', 'N3+.I-B5', 'N3+.N1-B0',
                              'N3+.N1-B1', 'N3+.N1-B2', 'N3+.N1-B3', 'N3+.N1-B4', 'N3+.N1-B5',
                              'N3+.N2-B0', 'N3+.N2-B1', 'N3+.N2-B2', 'N3+.N2-B3', 'N3+.N2-B4',
                              'N3+.N2-B5', 'N3+.N3-B0', 'N3+.N3-B1', 'N3+.N3-B2', 'N3+.N3-B3',
                              'N3+.N3-B4', 'N3+.N3-B5', 'N3+.N3+-B0', 'N3+.N3+-B1', 'N3+.N3+-B2',
                              'N3+.N3+-B3', 'N3+.N3+-B4', 'N3+.N3+-B5', 'N3+.Nam-B0',
                              'N3+.Nam-B1', 'N3+.Nam-B2', 'N3+.Nam-B3', 'N3+.Nam-B4',
                              'N3+.Nam-B5', 'N3+.Nar-B0', 'N3+.Nar-B1', 'N3+.Nar-B2',
                              'N3+.Nar-B3', 'N3+.Nar-B4', 'N3+.Nar-B5', 'N3+.Ng+-B0',
                              'N3+.Ng+-B1', 'N3+.Ng+-B2', 'N3+.Ng+-B3', 'N3+.Ng+-B4',
                              'N3+.Ng+-B5', 'N3+.Npl-B0', 'N3+.Npl-B1', 'N3+.Npl-B2',
                              'N3+.Npl-B3', 'N3+.Npl-B4', 'N3+.Npl-B5', 'N3+.Ntr-B0',
                              'N3+.Ntr-B1', 'N3+.Ntr-B2', 'N3+.Ntr-B3', 'N3+.Ntr-B4',
                              'N3+.Ntr-B5', 'N3+.O--B0', 'N3+.O--B1', 'N3+.O--B2', 'N3+.O--B3',
                              'N3+.O--B4', 'N3+.O--B5', 'N3+.O.co2-B0', 'N3+.O.co2-B1',
                              'N3+.O.co2-B2', 'N3+.O.co2-B3', 'N3+.O.co2-B4', 'N3+.O.co2-B5',
                              'N3+.O2-B0', 'N3+.O2-B1', 'N3+.O2-B2', 'N3+.O2-B3', 'N3+.O2-B4',
                              'N3+.O2-B5', 'N3+.O3-B0', 'N3+.O3-B1', 'N3+.O3-B2', 'N3+.O3-B3',
                              'N3+.O3-B4', 'N3+.O3-B5', 'N3+.P-B0', 'N3+.P-B1', 'N3+.P-B2',
                              'N3+.P-B3', 'N3+.P-B4', 'N3+.P-B5', 'N3+.Pac-B0', 'N3+.Pac-B1',
                              'N3+.Pac-B2', 'N3+.Pac-B3', 'N3+.Pac-B4', 'N3+.Pac-B5',
                              'N3+.S2-B0', 'N3+.S2-B1', 'N3+.S2-B2', 'N3+.S2-B3', 'N3+.S2-B4',
                              'N3+.S2-B5', 'N3+.S3-B0', 'N3+.S3-B1', 'N3+.S3-B2', 'N3+.S3-B3',
                              'N3+.S3-B4', 'N3+.S3-B5', 'N3+.So2-B0', 'N3+.So2-B1', 'N3+.So2-B2',
                              'N3+.So2-B3', 'N3+.So2-B4', 'N3+.So2-B5', 'N3+.Sox-B0',
                              'N3+.Sox-B1', 'N3+.Sox-B2', 'N3+.Sox-B3', 'N3+.Sox-B4',
                              'N3+.Sox-B5', 'Nam.Br-B0', 'Nam.Br-B1', 'Nam.Br-B2', 'Nam.Br-B3',
                              'Nam.Br-B4', 'Nam.Br-B5', 'Nam.C+-B0', 'Nam.C+-B1', 'Nam.C+-B2',
                              'Nam.C+-B3', 'Nam.C+-B4', 'Nam.C+-B5', 'Nam.C1-B0', 'Nam.C1-B1',
                              'Nam.C1-B2', 'Nam.C1-B3', 'Nam.C1-B4', 'Nam.C1-B5', 'Nam.C2-B0',
                              'Nam.C2-B1', 'Nam.C2-B2', 'Nam.C2-B3', 'Nam.C2-B4', 'Nam.C2-B5',
                              'Nam.C3-B0', 'Nam.C3-B1', 'Nam.C3-B2', 'Nam.C3-B3', 'Nam.C3-B4',
                              'Nam.C3-B5', 'Nam.Cac-B0', 'Nam.Cac-B1', 'Nam.Cac-B2',
                              'Nam.Cac-B3', 'Nam.Cac-B4', 'Nam.Cac-B5', 'Nam.Car-B0',
                              'Nam.Car-B1', 'Nam.Car-B2', 'Nam.Car-B3', 'Nam.Car-B4',
                              'Nam.Car-B5', 'Nam.Cl-B0', 'Nam.Cl-B1', 'Nam.Cl-B2', 'Nam.Cl-B3',
                              'Nam.Cl-B4', 'Nam.Cl-B5', 'Nam.F-B0', 'Nam.F-B1', 'Nam.F-B2',
                              'Nam.F-B3', 'Nam.F-B4', 'Nam.F-B5', 'Nam.I-B0', 'Nam.I-B1',
                              'Nam.I-B2', 'Nam.I-B3', 'Nam.I-B4', 'Nam.I-B5', 'Nam.N1-B0',
                              'Nam.N1-B1', 'Nam.N1-B2', 'Nam.N1-B3', 'Nam.N1-B4', 'Nam.N1-B5',
                              'Nam.N2-B0', 'Nam.N2-B1', 'Nam.N2-B2', 'Nam.N2-B3', 'Nam.N2-B4',
                              'Nam.N2-B5', 'Nam.N3-B0', 'Nam.N3-B1', 'Nam.N3-B2', 'Nam.N3-B3',
                              'Nam.N3-B4', 'Nam.N3-B5', 'Nam.N3+-B0', 'Nam.N3+-B1', 'Nam.N3+-B2',
                              'Nam.N3+-B3', 'Nam.N3+-B4', 'Nam.N3+-B5', 'Nam.Nam-B0',
                              'Nam.Nam-B1', 'Nam.Nam-B2', 'Nam.Nam-B3', 'Nam.Nam-B4',
                              'Nam.Nam-B5', 'Nam.Nar-B0', 'Nam.Nar-B1', 'Nam.Nar-B2',
                              'Nam.Nar-B3', 'Nam.Nar-B4', 'Nam.Nar-B5', 'Nam.Ng+-B0',
                              'Nam.Ng+-B1', 'Nam.Ng+-B2', 'Nam.Ng+-B3', 'Nam.Ng+-B4',
                              'Nam.Ng+-B5', 'Nam.Npl-B0', 'Nam.Npl-B1', 'Nam.Npl-B2',
                              'Nam.Npl-B3', 'Nam.Npl-B4', 'Nam.Npl-B5', 'Nam.Ntr-B0',
                              'Nam.Ntr-B1', 'Nam.Ntr-B2', 'Nam.Ntr-B3', 'Nam.Ntr-B4',
                              'Nam.Ntr-B5', 'Nam.O--B0', 'Nam.O--B1', 'Nam.O--B2', 'Nam.O--B3',
                              'Nam.O--B4', 'Nam.O--B5', 'Nam.O.co2-B0', 'Nam.O.co2-B1',
                              'Nam.O.co2-B2', 'Nam.O.co2-B3', 'Nam.O.co2-B4', 'Nam.O.co2-B5',
                              'Nam.O2-B0', 'Nam.O2-B1', 'Nam.O2-B2', 'Nam.O2-B3', 'Nam.O2-B4',
                              'Nam.O2-B5', 'Nam.O3-B0', 'Nam.O3-B1', 'Nam.O3-B2', 'Nam.O3-B3',
                              'Nam.O3-B4', 'Nam.O3-B5', 'Nam.P-B0', 'Nam.P-B1', 'Nam.P-B2',
                              'Nam.P-B3', 'Nam.P-B4', 'Nam.P-B5', 'Nam.Pac-B0', 'Nam.Pac-B1',
                              'Nam.Pac-B2', 'Nam.Pac-B3', 'Nam.Pac-B4', 'Nam.Pac-B5',
                              'Nam.S2-B0', 'Nam.S2-B1', 'Nam.S2-B2', 'Nam.S2-B3', 'Nam.S2-B4',
                              'Nam.S2-B5', 'Nam.S3-B0', 'Nam.S3-B1', 'Nam.S3-B2', 'Nam.S3-B3',
                              'Nam.S3-B4', 'Nam.S3-B5', 'Nam.So2-B0', 'Nam.So2-B1', 'Nam.So2-B2',
                              'Nam.So2-B3', 'Nam.So2-B4', 'Nam.So2-B5', 'Nam.Sox-B0',
                              'Nam.Sox-B1', 'Nam.Sox-B2', 'Nam.Sox-B3', 'Nam.Sox-B4',
                              'Nam.Sox-B5', 'Nar.Br-B0', 'Nar.Br-B1', 'Nar.Br-B2', 'Nar.Br-B3',
                              'Nar.Br-B4', 'Nar.Br-B5', 'Nar.C+-B0', 'Nar.C+-B1', 'Nar.C+-B2',
                              'Nar.C+-B3', 'Nar.C+-B4', 'Nar.C+-B5', 'Nar.C1-B0', 'Nar.C1-B1',
                              'Nar.C1-B2', 'Nar.C1-B3', 'Nar.C1-B4', 'Nar.C1-B5', 'Nar.C2-B0',
                              'Nar.C2-B1', 'Nar.C2-B2', 'Nar.C2-B3', 'Nar.C2-B4', 'Nar.C2-B5',
                              'Nar.C3-B0', 'Nar.C3-B1', 'Nar.C3-B2', 'Nar.C3-B3', 'Nar.C3-B4',
                              'Nar.C3-B5', 'Nar.Cac-B0', 'Nar.Cac-B1', 'Nar.Cac-B2',
                              'Nar.Cac-B3', 'Nar.Cac-B4', 'Nar.Cac-B5', 'Nar.Car-B0',
                              'Nar.Car-B1', 'Nar.Car-B2', 'Nar.Car-B3', 'Nar.Car-B4',
                              'Nar.Car-B5', 'Nar.Cl-B0', 'Nar.Cl-B1', 'Nar.Cl-B2', 'Nar.Cl-B3',
                              'Nar.Cl-B4', 'Nar.Cl-B5', 'Nar.F-B0', 'Nar.F-B1', 'Nar.F-B2',
                              'Nar.F-B3', 'Nar.F-B4', 'Nar.F-B5', 'Nar.I-B0', 'Nar.I-B1',
                              'Nar.I-B2', 'Nar.I-B3', 'Nar.I-B4', 'Nar.I-B5', 'Nar.N1-B0',
                              'Nar.N1-B1', 'Nar.N1-B2', 'Nar.N1-B3', 'Nar.N1-B4', 'Nar.N1-B5',
                              'Nar.N2-B0', 'Nar.N2-B1', 'Nar.N2-B2', 'Nar.N2-B3', 'Nar.N2-B4',
                              'Nar.N2-B5', 'Nar.N3-B0', 'Nar.N3-B1', 'Nar.N3-B2', 'Nar.N3-B3',
                              'Nar.N3-B4', 'Nar.N3-B5', 'Nar.N3+-B0', 'Nar.N3+-B1', 'Nar.N3+-B2',
                              'Nar.N3+-B3', 'Nar.N3+-B4', 'Nar.N3+-B5', 'Nar.Nam-B0',
                              'Nar.Nam-B1', 'Nar.Nam-B2', 'Nar.Nam-B3', 'Nar.Nam-B4',
                              'Nar.Nam-B5', 'Nar.Nar-B0', 'Nar.Nar-B1', 'Nar.Nar-B2',
                              'Nar.Nar-B3', 'Nar.Nar-B4', 'Nar.Nar-B5', 'Nar.Ng+-B0',
                              'Nar.Ng+-B1', 'Nar.Ng+-B2', 'Nar.Ng+-B3', 'Nar.Ng+-B4',
                              'Nar.Ng+-B5', 'Nar.Npl-B0', 'Nar.Npl-B1', 'Nar.Npl-B2',
                              'Nar.Npl-B3', 'Nar.Npl-B4', 'Nar.Npl-B5', 'Nar.Ntr-B0',
                              'Nar.Ntr-B1', 'Nar.Ntr-B2', 'Nar.Ntr-B3', 'Nar.Ntr-B4',
                              'Nar.Ntr-B5', 'Nar.O--B0', 'Nar.O--B1', 'Nar.O--B2', 'Nar.O--B3',
                              'Nar.O--B4', 'Nar.O--B5', 'Nar.O.co2-B0', 'Nar.O.co2-B1',
                              'Nar.O.co2-B2', 'Nar.O.co2-B3', 'Nar.O.co2-B4', 'Nar.O.co2-B5',
                              'Nar.O2-B0', 'Nar.O2-B1', 'Nar.O2-B2', 'Nar.O2-B3', 'Nar.O2-B4',
                              'Nar.O2-B5', 'Nar.O3-B0', 'Nar.O3-B1', 'Nar.O3-B2', 'Nar.O3-B3',
                              'Nar.O3-B4', 'Nar.O3-B5', 'Nar.P-B0', 'Nar.P-B1', 'Nar.P-B2',
                              'Nar.P-B3', 'Nar.P-B4', 'Nar.P-B5', 'Nar.Pac-B0', 'Nar.Pac-B1',
                              'Nar.Pac-B2', 'Nar.Pac-B3', 'Nar.Pac-B4', 'Nar.Pac-B5',
                              'Nar.S2-B0', 'Nar.S2-B1', 'Nar.S2-B2', 'Nar.S2-B3', 'Nar.S2-B4',
                              'Nar.S2-B5', 'Nar.S3-B0', 'Nar.S3-B1', 'Nar.S3-B2', 'Nar.S3-B3',
                              'Nar.S3-B4', 'Nar.S3-B5', 'Nar.So2-B0', 'Nar.So2-B1', 'Nar.So2-B2',
                              'Nar.So2-B3', 'Nar.So2-B4', 'Nar.So2-B5', 'Nar.Sox-B0',
                              'Nar.Sox-B1', 'Nar.Sox-B2', 'Nar.Sox-B3', 'Nar.Sox-B4',
                              'Nar.Sox-B5', 'Ng+.Br-B0', 'Ng+.Br-B1', 'Ng+.Br-B2', 'Ng+.Br-B3',
                              'Ng+.Br-B4', 'Ng+.Br-B5', 'Ng+.C+-B0', 'Ng+.C+-B1', 'Ng+.C+-B2',
                              'Ng+.C+-B3', 'Ng+.C+-B4', 'Ng+.C+-B5', 'Ng+.C1-B0', 'Ng+.C1-B1',
                              'Ng+.C1-B2', 'Ng+.C1-B3', 'Ng+.C1-B4', 'Ng+.C1-B5', 'Ng+.C2-B0',
                              'Ng+.C2-B1', 'Ng+.C2-B2', 'Ng+.C2-B3', 'Ng+.C2-B4', 'Ng+.C2-B5',
                              'Ng+.C3-B0', 'Ng+.C3-B1', 'Ng+.C3-B2', 'Ng+.C3-B3', 'Ng+.C3-B4',
                              'Ng+.C3-B5', 'Ng+.Cac-B0', 'Ng+.Cac-B1', 'Ng+.Cac-B2',
                              'Ng+.Cac-B3', 'Ng+.Cac-B4', 'Ng+.Cac-B5', 'Ng+.Car-B0',
                              'Ng+.Car-B1', 'Ng+.Car-B2', 'Ng+.Car-B3', 'Ng+.Car-B4',
                              'Ng+.Car-B5', 'Ng+.Cl-B0', 'Ng+.Cl-B1', 'Ng+.Cl-B2', 'Ng+.Cl-B3',
                              'Ng+.Cl-B4', 'Ng+.Cl-B5', 'Ng+.F-B0', 'Ng+.F-B1', 'Ng+.F-B2',
                              'Ng+.F-B3', 'Ng+.F-B4', 'Ng+.F-B5', 'Ng+.I-B0', 'Ng+.I-B1',
                              'Ng+.I-B2', 'Ng+.I-B3', 'Ng+.I-B4', 'Ng+.I-B5', 'Ng+.N1-B0',
                              'Ng+.N1-B1', 'Ng+.N1-B2', 'Ng+.N1-B3', 'Ng+.N1-B4', 'Ng+.N1-B5',
                              'Ng+.N2-B0', 'Ng+.N2-B1', 'Ng+.N2-B2', 'Ng+.N2-B3', 'Ng+.N2-B4',
                              'Ng+.N2-B5', 'Ng+.N3-B0', 'Ng+.N3-B1', 'Ng+.N3-B2', 'Ng+.N3-B3',
                              'Ng+.N3-B4', 'Ng+.N3-B5', 'Ng+.N3+-B0', 'Ng+.N3+-B1', 'Ng+.N3+-B2',
                              'Ng+.N3+-B3', 'Ng+.N3+-B4', 'Ng+.N3+-B5', 'Ng+.Nam-B0',
                              'Ng+.Nam-B1', 'Ng+.Nam-B2', 'Ng+.Nam-B3', 'Ng+.Nam-B4',
                              'Ng+.Nam-B5', 'Ng+.Nar-B0', 'Ng+.Nar-B1', 'Ng+.Nar-B2',
                              'Ng+.Nar-B3', 'Ng+.Nar-B4', 'Ng+.Nar-B5', 'Ng+.Ng+-B0',
                              'Ng+.Ng+-B1', 'Ng+.Ng+-B2', 'Ng+.Ng+-B3', 'Ng+.Ng+-B4',
                              'Ng+.Ng+-B5', 'Ng+.Npl-B0', 'Ng+.Npl-B1', 'Ng+.Npl-B2',
                              'Ng+.Npl-B3', 'Ng+.Npl-B4', 'Ng+.Npl-B5', 'Ng+.Ntr-B0',
                              'Ng+.Ntr-B1', 'Ng+.Ntr-B2', 'Ng+.Ntr-B3', 'Ng+.Ntr-B4',
                              'Ng+.Ntr-B5', 'Ng+.O--B0', 'Ng+.O--B1', 'Ng+.O--B2', 'Ng+.O--B3',
                              'Ng+.O--B4', 'Ng+.O--B5', 'Ng+.O.co2-B0', 'Ng+.O.co2-B1',
                              'Ng+.O.co2-B2', 'Ng+.O.co2-B3', 'Ng+.O.co2-B4', 'Ng+.O.co2-B5',
                              'Ng+.O2-B0', 'Ng+.O2-B1', 'Ng+.O2-B2', 'Ng+.O2-B3', 'Ng+.O2-B4',
                              'Ng+.O2-B5', 'Ng+.O3-B0', 'Ng+.O3-B1', 'Ng+.O3-B2', 'Ng+.O3-B3',
                              'Ng+.O3-B4', 'Ng+.O3-B5', 'Ng+.P-B0', 'Ng+.P-B1', 'Ng+.P-B2',
                              'Ng+.P-B3', 'Ng+.P-B4', 'Ng+.P-B5', 'Ng+.Pac-B0', 'Ng+.Pac-B1',
                              'Ng+.Pac-B2', 'Ng+.Pac-B3', 'Ng+.Pac-B4', 'Ng+.Pac-B5',
                              'Ng+.S2-B0', 'Ng+.S2-B1', 'Ng+.S2-B2', 'Ng+.S2-B3', 'Ng+.S2-B4',
                              'Ng+.S2-B5', 'Ng+.S3-B0', 'Ng+.S3-B1', 'Ng+.S3-B2', 'Ng+.S3-B3',
                              'Ng+.S3-B4', 'Ng+.S3-B5', 'Ng+.So2-B0', 'Ng+.So2-B1', 'Ng+.So2-B2',
                              'Ng+.So2-B3', 'Ng+.So2-B4', 'Ng+.So2-B5', 'Ng+.Sox-B0',
                              'Ng+.Sox-B1', 'Ng+.Sox-B2', 'Ng+.Sox-B3', 'Ng+.Sox-B4',
                              'Ng+.Sox-B5', 'Npl.Br-B0', 'Npl.Br-B1', 'Npl.Br-B2', 'Npl.Br-B3',
                              'Npl.Br-B4', 'Npl.Br-B5', 'Npl.C+-B0', 'Npl.C+-B1', 'Npl.C+-B2',
                              'Npl.C+-B3', 'Npl.C+-B4', 'Npl.C+-B5', 'Npl.C1-B0', 'Npl.C1-B1',
                              'Npl.C1-B2', 'Npl.C1-B3', 'Npl.C1-B4', 'Npl.C1-B5', 'Npl.C2-B0',
                              'Npl.C2-B1', 'Npl.C2-B2', 'Npl.C2-B3', 'Npl.C2-B4', 'Npl.C2-B5',
                              'Npl.C3-B0', 'Npl.C3-B1', 'Npl.C3-B2', 'Npl.C3-B3', 'Npl.C3-B4',
                              'Npl.C3-B5', 'Npl.Cac-B0', 'Npl.Cac-B1', 'Npl.Cac-B2',
                              'Npl.Cac-B3', 'Npl.Cac-B4', 'Npl.Cac-B5', 'Npl.Car-B0',
                              'Npl.Car-B1', 'Npl.Car-B2', 'Npl.Car-B3', 'Npl.Car-B4',
                              'Npl.Car-B5', 'Npl.Cl-B0', 'Npl.Cl-B1', 'Npl.Cl-B2', 'Npl.Cl-B3',
                              'Npl.Cl-B4', 'Npl.Cl-B5', 'Npl.F-B0', 'Npl.F-B1', 'Npl.F-B2',
                              'Npl.F-B3', 'Npl.F-B4', 'Npl.F-B5', 'Npl.I-B0', 'Npl.I-B1',
                              'Npl.I-B2', 'Npl.I-B3', 'Npl.I-B4', 'Npl.I-B5', 'Npl.N1-B0',
                              'Npl.N1-B1', 'Npl.N1-B2', 'Npl.N1-B3', 'Npl.N1-B4', 'Npl.N1-B5',
                              'Npl.N2-B0', 'Npl.N2-B1', 'Npl.N2-B2', 'Npl.N2-B3', 'Npl.N2-B4',
                              'Npl.N2-B5', 'Npl.N3-B0', 'Npl.N3-B1', 'Npl.N3-B2', 'Npl.N3-B3',
                              'Npl.N3-B4', 'Npl.N3-B5', 'Npl.N3+-B0', 'Npl.N3+-B1', 'Npl.N3+-B2',
                              'Npl.N3+-B3', 'Npl.N3+-B4', 'Npl.N3+-B5', 'Npl.Nam-B0',
                              'Npl.Nam-B1', 'Npl.Nam-B2', 'Npl.Nam-B3', 'Npl.Nam-B4',
                              'Npl.Nam-B5', 'Npl.Nar-B0', 'Npl.Nar-B1', 'Npl.Nar-B2',
                              'Npl.Nar-B3', 'Npl.Nar-B4', 'Npl.Nar-B5', 'Npl.Ng+-B0',
                              'Npl.Ng+-B1', 'Npl.Ng+-B2', 'Npl.Ng+-B3', 'Npl.Ng+-B4',
                              'Npl.Ng+-B5', 'Npl.Npl-B0', 'Npl.Npl-B1', 'Npl.Npl-B2',
                              'Npl.Npl-B3', 'Npl.Npl-B4', 'Npl.Npl-B5', 'Npl.Ntr-B0',
                              'Npl.Ntr-B1', 'Npl.Ntr-B2', 'Npl.Ntr-B3', 'Npl.Ntr-B4',
                              'Npl.Ntr-B5', 'Npl.O--B0', 'Npl.O--B1', 'Npl.O--B2', 'Npl.O--B3',
                              'Npl.O--B4', 'Npl.O--B5', 'Npl.O.co2-B0', 'Npl.O.co2-B1',
                              'Npl.O.co2-B2', 'Npl.O.co2-B3', 'Npl.O.co2-B4', 'Npl.O.co2-B5',
                              'Npl.O2-B0', 'Npl.O2-B1', 'Npl.O2-B2', 'Npl.O2-B3', 'Npl.O2-B4',
                              'Npl.O2-B5', 'Npl.O3-B0', 'Npl.O3-B1', 'Npl.O3-B2', 'Npl.O3-B3',
                              'Npl.O3-B4', 'Npl.O3-B5', 'Npl.P-B0', 'Npl.P-B1', 'Npl.P-B2',
                              'Npl.P-B3', 'Npl.P-B4', 'Npl.P-B5', 'Npl.Pac-B0', 'Npl.Pac-B1',
                              'Npl.Pac-B2', 'Npl.Pac-B3', 'Npl.Pac-B4', 'Npl.Pac-B5',
                              'Npl.S2-B0', 'Npl.S2-B1', 'Npl.S2-B2', 'Npl.S2-B3', 'Npl.S2-B4',
                              'Npl.S2-B5', 'Npl.S3-B0', 'Npl.S3-B1', 'Npl.S3-B2', 'Npl.S3-B3',
                              'Npl.S3-B4', 'Npl.S3-B5', 'Npl.So2-B0', 'Npl.So2-B1', 'Npl.So2-B2',
                              'Npl.So2-B3', 'Npl.So2-B4', 'Npl.So2-B5', 'Npl.Sox-B0',
                              'Npl.Sox-B1', 'Npl.Sox-B2', 'Npl.Sox-B3', 'Npl.Sox-B4',
                              'Npl.Sox-B5', 'Ntr.Br-B0', 'Ntr.Br-B1', 'Ntr.Br-B2', 'Ntr.Br-B3',
                              'Ntr.Br-B4', 'Ntr.Br-B5', 'Ntr.C+-B0', 'Ntr.C+-B1', 'Ntr.C+-B2',
                              'Ntr.C+-B3', 'Ntr.C+-B4', 'Ntr.C+-B5', 'Ntr.C1-B0', 'Ntr.C1-B1',
                              'Ntr.C1-B2', 'Ntr.C1-B3', 'Ntr.C1-B4', 'Ntr.C1-B5', 'Ntr.C2-B0',
                              'Ntr.C2-B1', 'Ntr.C2-B2', 'Ntr.C2-B3', 'Ntr.C2-B4', 'Ntr.C2-B5',
                              'Ntr.C3-B0', 'Ntr.C3-B1', 'Ntr.C3-B2', 'Ntr.C3-B3', 'Ntr.C3-B4',
                              'Ntr.C3-B5', 'Ntr.Cac-B0', 'Ntr.Cac-B1', 'Ntr.Cac-B2',
                              'Ntr.Cac-B3', 'Ntr.Cac-B4', 'Ntr.Cac-B5', 'Ntr.Car-B0',
                              'Ntr.Car-B1', 'Ntr.Car-B2', 'Ntr.Car-B3', 'Ntr.Car-B4',
                              'Ntr.Car-B5', 'Ntr.Cl-B0', 'Ntr.Cl-B1', 'Ntr.Cl-B2', 'Ntr.Cl-B3',
                              'Ntr.Cl-B4', 'Ntr.Cl-B5', 'Ntr.F-B0', 'Ntr.F-B1', 'Ntr.F-B2',
                              'Ntr.F-B3', 'Ntr.F-B4', 'Ntr.F-B5', 'Ntr.I-B0', 'Ntr.I-B1',
                              'Ntr.I-B2', 'Ntr.I-B3', 'Ntr.I-B4', 'Ntr.I-B5', 'Ntr.N1-B0',
                              'Ntr.N1-B1', 'Ntr.N1-B2', 'Ntr.N1-B3', 'Ntr.N1-B4', 'Ntr.N1-B5',
                              'Ntr.N2-B0', 'Ntr.N2-B1', 'Ntr.N2-B2', 'Ntr.N2-B3', 'Ntr.N2-B4',
                              'Ntr.N2-B5', 'Ntr.N3-B0', 'Ntr.N3-B1', 'Ntr.N3-B2', 'Ntr.N3-B3',
                              'Ntr.N3-B4', 'Ntr.N3-B5', 'Ntr.N3+-B0', 'Ntr.N3+-B1', 'Ntr.N3+-B2',
                              'Ntr.N3+-B3', 'Ntr.N3+-B4', 'Ntr.N3+-B5', 'Ntr.Nam-B0',
                              'Ntr.Nam-B1', 'Ntr.Nam-B2', 'Ntr.Nam-B3', 'Ntr.Nam-B4',
                              'Ntr.Nam-B5', 'Ntr.Nar-B0', 'Ntr.Nar-B1', 'Ntr.Nar-B2',
                              'Ntr.Nar-B3', 'Ntr.Nar-B4', 'Ntr.Nar-B5', 'Ntr.Ng+-B0',
                              'Ntr.Ng+-B1', 'Ntr.Ng+-B2', 'Ntr.Ng+-B3', 'Ntr.Ng+-B4',
                              'Ntr.Ng+-B5', 'Ntr.Npl-B0', 'Ntr.Npl-B1', 'Ntr.Npl-B2',
                              'Ntr.Npl-B3', 'Ntr.Npl-B4', 'Ntr.Npl-B5', 'Ntr.Ntr-B0',
                              'Ntr.Ntr-B1', 'Ntr.Ntr-B2', 'Ntr.Ntr-B3', 'Ntr.Ntr-B4',
                              'Ntr.Ntr-B5', 'Ntr.O--B0', 'Ntr.O--B1', 'Ntr.O--B2', 'Ntr.O--B3',
                              'Ntr.O--B4', 'Ntr.O--B5', 'Ntr.O.co2-B0', 'Ntr.O.co2-B1',
                              'Ntr.O.co2-B2', 'Ntr.O.co2-B3', 'Ntr.O.co2-B4', 'Ntr.O.co2-B5',
                              'Ntr.O2-B0', 'Ntr.O2-B1', 'Ntr.O2-B2', 'Ntr.O2-B3', 'Ntr.O2-B4',
                              'Ntr.O2-B5', 'Ntr.O3-B0', 'Ntr.O3-B1', 'Ntr.O3-B2', 'Ntr.O3-B3',
                              'Ntr.O3-B4', 'Ntr.O3-B5', 'Ntr.P-B0', 'Ntr.P-B1', 'Ntr.P-B2',
                              'Ntr.P-B3', 'Ntr.P-B4', 'Ntr.P-B5', 'Ntr.Pac-B0', 'Ntr.Pac-B1',
                              'Ntr.Pac-B2', 'Ntr.Pac-B3', 'Ntr.Pac-B4', 'Ntr.Pac-B5',
                              'Ntr.S2-B0', 'Ntr.S2-B1', 'Ntr.S2-B2', 'Ntr.S2-B3', 'Ntr.S2-B4',
                              'Ntr.S2-B5', 'Ntr.S3-B0', 'Ntr.S3-B1', 'Ntr.S3-B2', 'Ntr.S3-B3',
                              'Ntr.S3-B4', 'Ntr.S3-B5', 'Ntr.So2-B0', 'Ntr.So2-B1', 'Ntr.So2-B2',
                              'Ntr.So2-B3', 'Ntr.So2-B4', 'Ntr.So2-B5', 'Ntr.Sox-B0',
                              'Ntr.Sox-B1', 'Ntr.Sox-B2', 'Ntr.Sox-B3', 'Ntr.Sox-B4',
                              'Ntr.Sox-B5', 'O-.Br-B0', 'O-.Br-B1', 'O-.Br-B2', 'O-.Br-B3',
                              'O-.Br-B4', 'O-.Br-B5', 'O-.C+-B0', 'O-.C+-B1', 'O-.C+-B2',
                              'O-.C+-B3', 'O-.C+-B4', 'O-.C+-B5', 'O-.C1-B0', 'O-.C1-B1',
                              'O-.C1-B2', 'O-.C1-B3', 'O-.C1-B4', 'O-.C1-B5', 'O-.C2-B0',
                              'O-.C2-B1', 'O-.C2-B2', 'O-.C2-B3', 'O-.C2-B4', 'O-.C2-B5',
                              'O-.C3-B0', 'O-.C3-B1', 'O-.C3-B2', 'O-.C3-B3', 'O-.C3-B4',
                              'O-.C3-B5', 'O-.Cac-B0', 'O-.Cac-B1', 'O-.Cac-B2', 'O-.Cac-B3',
                              'O-.Cac-B4', 'O-.Cac-B5', 'O-.Car-B0', 'O-.Car-B1', 'O-.Car-B2',
                              'O-.Car-B3', 'O-.Car-B4', 'O-.Car-B5', 'O-.Cl-B0', 'O-.Cl-B1',
                              'O-.Cl-B2', 'O-.Cl-B3', 'O-.Cl-B4', 'O-.Cl-B5', 'O-.F-B0',
                              'O-.F-B1', 'O-.F-B2', 'O-.F-B3', 'O-.F-B4', 'O-.F-B5', 'O-.I-B0',
                              'O-.I-B1', 'O-.I-B2', 'O-.I-B3', 'O-.I-B4', 'O-.I-B5', 'O-.N1-B0',
                              'O-.N1-B1', 'O-.N1-B2', 'O-.N1-B3', 'O-.N1-B4', 'O-.N1-B5',
                              'O-.N2-B0', 'O-.N2-B1', 'O-.N2-B2', 'O-.N2-B3', 'O-.N2-B4',
                              'O-.N2-B5', 'O-.N3-B0', 'O-.N3-B1', 'O-.N3-B2', 'O-.N3-B3',
                              'O-.N3-B4', 'O-.N3-B5', 'O-.N3+-B0', 'O-.N3+-B1', 'O-.N3+-B2',
                              'O-.N3+-B3', 'O-.N3+-B4', 'O-.N3+-B5', 'O-.Nam-B0', 'O-.Nam-B1',
                              'O-.Nam-B2', 'O-.Nam-B3', 'O-.Nam-B4', 'O-.Nam-B5', 'O-.Nar-B0',
                              'O-.Nar-B1', 'O-.Nar-B2', 'O-.Nar-B3', 'O-.Nar-B4', 'O-.Nar-B5',
                              'O-.Ng+-B0', 'O-.Ng+-B1', 'O-.Ng+-B2', 'O-.Ng+-B3', 'O-.Ng+-B4',
                              'O-.Ng+-B5', 'O-.Npl-B0', 'O-.Npl-B1', 'O-.Npl-B2', 'O-.Npl-B3',
                              'O-.Npl-B4', 'O-.Npl-B5', 'O-.Ntr-B0', 'O-.Ntr-B1', 'O-.Ntr-B2',
                              'O-.Ntr-B3', 'O-.Ntr-B4', 'O-.Ntr-B5', 'O-.O--B0', 'O-.O--B1',
                              'O-.O--B2', 'O-.O--B3', 'O-.O--B4', 'O-.O--B5', 'O-.O.co2-B0',
                              'O-.O.co2-B1', 'O-.O.co2-B2', 'O-.O.co2-B3', 'O-.O.co2-B4',
                              'O-.O.co2-B5', 'O-.O2-B0', 'O-.O2-B1', 'O-.O2-B2', 'O-.O2-B3',
                              'O-.O2-B4', 'O-.O2-B5', 'O-.O3-B0', 'O-.O3-B1', 'O-.O3-B2',
                              'O-.O3-B3', 'O-.O3-B4', 'O-.O3-B5', 'O-.P-B0', 'O-.P-B1',
                              'O-.P-B2', 'O-.P-B3', 'O-.P-B4', 'O-.P-B5', 'O-.Pac-B0',
                              'O-.Pac-B1', 'O-.Pac-B2', 'O-.Pac-B3', 'O-.Pac-B4', 'O-.Pac-B5',
                              'O-.S2-B0', 'O-.S2-B1', 'O-.S2-B2', 'O-.S2-B3', 'O-.S2-B4',
                              'O-.S2-B5', 'O-.S3-B0', 'O-.S3-B1', 'O-.S3-B2', 'O-.S3-B3',
                              'O-.S3-B4', 'O-.S3-B5', 'O-.So2-B0', 'O-.So2-B1', 'O-.So2-B2',
                              'O-.So2-B3', 'O-.So2-B4', 'O-.So2-B5', 'O-.Sox-B0', 'O-.Sox-B1',
                              'O-.Sox-B2', 'O-.Sox-B3', 'O-.Sox-B4', 'O-.Sox-B5', 'O.co2.Br-B0',
                              'O.co2.Br-B1', 'O.co2.Br-B2', 'O.co2.Br-B3', 'O.co2.Br-B4',
                              'O.co2.Br-B5', 'O.co2.C+-B0', 'O.co2.C+-B1', 'O.co2.C+-B2',
                              'O.co2.C+-B3', 'O.co2.C+-B4', 'O.co2.C+-B5', 'O.co2.C1-B0',
                              'O.co2.C1-B1', 'O.co2.C1-B2', 'O.co2.C1-B3', 'O.co2.C1-B4',
                              'O.co2.C1-B5', 'O.co2.C2-B0', 'O.co2.C2-B1', 'O.co2.C2-B2',
                              'O.co2.C2-B3', 'O.co2.C2-B4', 'O.co2.C2-B5', 'O.co2.C3-B0',
                              'O.co2.C3-B1', 'O.co2.C3-B2', 'O.co2.C3-B3', 'O.co2.C3-B4',
                              'O.co2.C3-B5', 'O.co2.Cac-B0', 'O.co2.Cac-B1', 'O.co2.Cac-B2',
                              'O.co2.Cac-B3', 'O.co2.Cac-B4', 'O.co2.Cac-B5', 'O.co2.Car-B0',
                              'O.co2.Car-B1', 'O.co2.Car-B2', 'O.co2.Car-B3', 'O.co2.Car-B4',
                              'O.co2.Car-B5', 'O.co2.Cl-B0', 'O.co2.Cl-B1', 'O.co2.Cl-B2',
                              'O.co2.Cl-B3', 'O.co2.Cl-B4', 'O.co2.Cl-B5', 'O.co2.F-B0',
                              'O.co2.F-B1', 'O.co2.F-B2', 'O.co2.F-B3', 'O.co2.F-B4',
                              'O.co2.F-B5', 'O.co2.I-B0', 'O.co2.I-B1', 'O.co2.I-B2',
                              'O.co2.I-B3', 'O.co2.I-B4', 'O.co2.I-B5', 'O.co2.N1-B0',
                              'O.co2.N1-B1', 'O.co2.N1-B2', 'O.co2.N1-B3', 'O.co2.N1-B4',
                              'O.co2.N1-B5', 'O.co2.N2-B0', 'O.co2.N2-B1', 'O.co2.N2-B2',
                              'O.co2.N2-B3', 'O.co2.N2-B4', 'O.co2.N2-B5', 'O.co2.N3-B0',
                              'O.co2.N3-B1', 'O.co2.N3-B2', 'O.co2.N3-B3', 'O.co2.N3-B4',
                              'O.co2.N3-B5', 'O.co2.N3+-B0', 'O.co2.N3+-B1', 'O.co2.N3+-B2',
                              'O.co2.N3+-B3', 'O.co2.N3+-B4', 'O.co2.N3+-B5', 'O.co2.Nam-B0',
                              'O.co2.Nam-B1', 'O.co2.Nam-B2', 'O.co2.Nam-B3', 'O.co2.Nam-B4',
                              'O.co2.Nam-B5', 'O.co2.Nar-B0', 'O.co2.Nar-B1', 'O.co2.Nar-B2',
                              'O.co2.Nar-B3', 'O.co2.Nar-B4', 'O.co2.Nar-B5', 'O.co2.Ng+-B0',
                              'O.co2.Ng+-B1', 'O.co2.Ng+-B2', 'O.co2.Ng+-B3', 'O.co2.Ng+-B4',
                              'O.co2.Ng+-B5', 'O.co2.Npl-B0', 'O.co2.Npl-B1', 'O.co2.Npl-B2',
                              'O.co2.Npl-B3', 'O.co2.Npl-B4', 'O.co2.Npl-B5', 'O.co2.Ntr-B0',
                              'O.co2.Ntr-B1', 'O.co2.Ntr-B2', 'O.co2.Ntr-B3', 'O.co2.Ntr-B4',
                              'O.co2.Ntr-B5', 'O.co2.O--B0', 'O.co2.O--B1', 'O.co2.O--B2',
                              'O.co2.O--B3', 'O.co2.O--B4', 'O.co2.O--B5', 'O.co2.O.co2-B0',
                              'O.co2.O.co2-B1', 'O.co2.O.co2-B2', 'O.co2.O.co2-B3',
                              'O.co2.O.co2-B4', 'O.co2.O.co2-B5', 'O.co2.O2-B0', 'O.co2.O2-B1',
                              'O.co2.O2-B2', 'O.co2.O2-B3', 'O.co2.O2-B4', 'O.co2.O2-B5',
                              'O.co2.O3-B0', 'O.co2.O3-B1', 'O.co2.O3-B2', 'O.co2.O3-B3',
                              'O.co2.O3-B4', 'O.co2.O3-B5', 'O.co2.P-B0', 'O.co2.P-B1',
                              'O.co2.P-B2', 'O.co2.P-B3', 'O.co2.P-B4', 'O.co2.P-B5',
                              'O.co2.Pac-B0', 'O.co2.Pac-B1', 'O.co2.Pac-B2', 'O.co2.Pac-B3',
                              'O.co2.Pac-B4', 'O.co2.Pac-B5', 'O.co2.S2-B0', 'O.co2.S2-B1',
                              'O.co2.S2-B2', 'O.co2.S2-B3', 'O.co2.S2-B4', 'O.co2.S2-B5',
                              'O.co2.S3-B0', 'O.co2.S3-B1', 'O.co2.S3-B2', 'O.co2.S3-B3',
                              'O.co2.S3-B4', 'O.co2.S3-B5', 'O.co2.So2-B0', 'O.co2.So2-B1',
                              'O.co2.So2-B2', 'O.co2.So2-B3', 'O.co2.So2-B4', 'O.co2.So2-B5',
                              'O.co2.Sox-B0', 'O.co2.Sox-B1', 'O.co2.Sox-B2', 'O.co2.Sox-B3',
                              'O.co2.Sox-B4', 'O.co2.Sox-B5', 'O2.Br-B0', 'O2.Br-B1', 'O2.Br-B2',
                              'O2.Br-B3', 'O2.Br-B4', 'O2.Br-B5', 'O2.C+-B0', 'O2.C+-B1',
                              'O2.C+-B2', 'O2.C+-B3', 'O2.C+-B4', 'O2.C+-B5', 'O2.C1-B0',
                              'O2.C1-B1', 'O2.C1-B2', 'O2.C1-B3', 'O2.C1-B4', 'O2.C1-B5',
                              'O2.C2-B0', 'O2.C2-B1', 'O2.C2-B2', 'O2.C2-B3', 'O2.C2-B4',
                              'O2.C2-B5', 'O2.C3-B0', 'O2.C3-B1', 'O2.C3-B2', 'O2.C3-B3',
                              'O2.C3-B4', 'O2.C3-B5', 'O2.Cac-B0', 'O2.Cac-B1', 'O2.Cac-B2',
                              'O2.Cac-B3', 'O2.Cac-B4', 'O2.Cac-B5', 'O2.Car-B0', 'O2.Car-B1',
                              'O2.Car-B2', 'O2.Car-B3', 'O2.Car-B4', 'O2.Car-B5', 'O2.Cl-B0',
                              'O2.Cl-B1', 'O2.Cl-B2', 'O2.Cl-B3', 'O2.Cl-B4', 'O2.Cl-B5',
                              'O2.F-B0', 'O2.F-B1', 'O2.F-B2', 'O2.F-B3', 'O2.F-B4', 'O2.F-B5',
                              'O2.I-B0', 'O2.I-B1', 'O2.I-B2', 'O2.I-B3', 'O2.I-B4', 'O2.I-B5',
                              'O2.N1-B0', 'O2.N1-B1', 'O2.N1-B2', 'O2.N1-B3', 'O2.N1-B4',
                              'O2.N1-B5', 'O2.N2-B0', 'O2.N2-B1', 'O2.N2-B2', 'O2.N2-B3',
                              'O2.N2-B4', 'O2.N2-B5', 'O2.N3-B0', 'O2.N3-B1', 'O2.N3-B2',
                              'O2.N3-B3', 'O2.N3-B4', 'O2.N3-B5', 'O2.N3+-B0', 'O2.N3+-B1',
                              'O2.N3+-B2', 'O2.N3+-B3', 'O2.N3+-B4', 'O2.N3+-B5', 'O2.Nam-B0',
                              'O2.Nam-B1', 'O2.Nam-B2', 'O2.Nam-B3', 'O2.Nam-B4', 'O2.Nam-B5',
                              'O2.Nar-B0', 'O2.Nar-B1', 'O2.Nar-B2', 'O2.Nar-B3', 'O2.Nar-B4',
                              'O2.Nar-B5', 'O2.Ng+-B0', 'O2.Ng+-B1', 'O2.Ng+-B2', 'O2.Ng+-B3',
                              'O2.Ng+-B4', 'O2.Ng+-B5', 'O2.Npl-B0', 'O2.Npl-B1', 'O2.Npl-B2',
                              'O2.Npl-B3', 'O2.Npl-B4', 'O2.Npl-B5', 'O2.Ntr-B0', 'O2.Ntr-B1',
                              'O2.Ntr-B2', 'O2.Ntr-B3', 'O2.Ntr-B4', 'O2.Ntr-B5', 'O2.O--B0',
                              'O2.O--B1', 'O2.O--B2', 'O2.O--B3', 'O2.O--B4', 'O2.O--B5',
                              'O2.O.co2-B0', 'O2.O.co2-B1', 'O2.O.co2-B2', 'O2.O.co2-B3',
                              'O2.O.co2-B4', 'O2.O.co2-B5', 'O2.O2-B0', 'O2.O2-B1', 'O2.O2-B2',
                              'O2.O2-B3', 'O2.O2-B4', 'O2.O2-B5', 'O2.O3-B0', 'O2.O3-B1',
                              'O2.O3-B2', 'O2.O3-B3', 'O2.O3-B4', 'O2.O3-B5', 'O2.P-B0',
                              'O2.P-B1', 'O2.P-B2', 'O2.P-B3', 'O2.P-B4', 'O2.P-B5', 'O2.Pac-B0',
                              'O2.Pac-B1', 'O2.Pac-B2', 'O2.Pac-B3', 'O2.Pac-B4', 'O2.Pac-B5',
                              'O2.S2-B0', 'O2.S2-B1', 'O2.S2-B2', 'O2.S2-B3', 'O2.S2-B4',
                              'O2.S2-B5', 'O2.S3-B0', 'O2.S3-B1', 'O2.S3-B2', 'O2.S3-B3',
                              'O2.S3-B4', 'O2.S3-B5', 'O2.So2-B0', 'O2.So2-B1', 'O2.So2-B2',
                              'O2.So2-B3', 'O2.So2-B4', 'O2.So2-B5', 'O2.Sox-B0', 'O2.Sox-B1',
                              'O2.Sox-B2', 'O2.Sox-B3', 'O2.Sox-B4', 'O2.Sox-B5', 'O3.Br-B0',
                              'O3.Br-B1', 'O3.Br-B2', 'O3.Br-B3', 'O3.Br-B4', 'O3.Br-B5',
                              'O3.C+-B0', 'O3.C+-B1', 'O3.C+-B2', 'O3.C+-B3', 'O3.C+-B4',
                              'O3.C+-B5', 'O3.C1-B0', 'O3.C1-B1', 'O3.C1-B2', 'O3.C1-B3',
                              'O3.C1-B4', 'O3.C1-B5', 'O3.C2-B0', 'O3.C2-B1', 'O3.C2-B2',
                              'O3.C2-B3', 'O3.C2-B4', 'O3.C2-B5', 'O3.C3-B0', 'O3.C3-B1',
                              'O3.C3-B2', 'O3.C3-B3', 'O3.C3-B4', 'O3.C3-B5', 'O3.Cac-B0',
                              'O3.Cac-B1', 'O3.Cac-B2', 'O3.Cac-B3', 'O3.Cac-B4', 'O3.Cac-B5',
                              'O3.Car-B0', 'O3.Car-B1', 'O3.Car-B2', 'O3.Car-B3', 'O3.Car-B4',
                              'O3.Car-B5', 'O3.Cl-B0', 'O3.Cl-B1', 'O3.Cl-B2', 'O3.Cl-B3',
                              'O3.Cl-B4', 'O3.Cl-B5', 'O3.F-B0', 'O3.F-B1', 'O3.F-B2', 'O3.F-B3',
                              'O3.F-B4', 'O3.F-B5', 'O3.I-B0', 'O3.I-B1', 'O3.I-B2', 'O3.I-B3',
                              'O3.I-B4', 'O3.I-B5', 'O3.N1-B0', 'O3.N1-B1', 'O3.N1-B2',
                              'O3.N1-B3', 'O3.N1-B4', 'O3.N1-B5', 'O3.N2-B0', 'O3.N2-B1',
                              'O3.N2-B2', 'O3.N2-B3', 'O3.N2-B4', 'O3.N2-B5', 'O3.N3-B0',
                              'O3.N3-B1', 'O3.N3-B2', 'O3.N3-B3', 'O3.N3-B4', 'O3.N3-B5',
                              'O3.N3+-B0', 'O3.N3+-B1', 'O3.N3+-B2', 'O3.N3+-B3', 'O3.N3+-B4',
                              'O3.N3+-B5', 'O3.Nam-B0', 'O3.Nam-B1', 'O3.Nam-B2', 'O3.Nam-B3',
                              'O3.Nam-B4', 'O3.Nam-B5', 'O3.Nar-B0', 'O3.Nar-B1', 'O3.Nar-B2',
                              'O3.Nar-B3', 'O3.Nar-B4', 'O3.Nar-B5', 'O3.Ng+-B0', 'O3.Ng+-B1',
                              'O3.Ng+-B2', 'O3.Ng+-B3', 'O3.Ng+-B4', 'O3.Ng+-B5', 'O3.Npl-B0',
                              'O3.Npl-B1', 'O3.Npl-B2', 'O3.Npl-B3', 'O3.Npl-B4', 'O3.Npl-B5',
                              'O3.Ntr-B0', 'O3.Ntr-B1', 'O3.Ntr-B2', 'O3.Ntr-B3', 'O3.Ntr-B4',
                              'O3.Ntr-B5', 'O3.O--B0', 'O3.O--B1', 'O3.O--B2', 'O3.O--B3',
                              'O3.O--B4', 'O3.O--B5', 'O3.O.co2-B0', 'O3.O.co2-B1',
                              'O3.O.co2-B2', 'O3.O.co2-B3', 'O3.O.co2-B4', 'O3.O.co2-B5',
                              'O3.O2-B0', 'O3.O2-B1', 'O3.O2-B2', 'O3.O2-B3', 'O3.O2-B4',
                              'O3.O2-B5', 'O3.O3-B0', 'O3.O3-B1', 'O3.O3-B2', 'O3.O3-B3',
                              'O3.O3-B4', 'O3.O3-B5', 'O3.P-B0', 'O3.P-B1', 'O3.P-B2', 'O3.P-B3',
                              'O3.P-B4', 'O3.P-B5', 'O3.Pac-B0', 'O3.Pac-B1', 'O3.Pac-B2',
                              'O3.Pac-B3', 'O3.Pac-B4', 'O3.Pac-B5', 'O3.S2-B0', 'O3.S2-B1',
                              'O3.S2-B2', 'O3.S2-B3', 'O3.S2-B4', 'O3.S2-B5', 'O3.S3-B0',
                              'O3.S3-B1', 'O3.S3-B2', 'O3.S3-B3', 'O3.S3-B4', 'O3.S3-B5',
                              'O3.So2-B0', 'O3.So2-B1', 'O3.So2-B2', 'O3.So2-B3', 'O3.So2-B4',
                              'O3.So2-B5', 'O3.Sox-B0', 'O3.Sox-B1', 'O3.Sox-B2', 'O3.Sox-B3',
                              'O3.Sox-B4', 'O3.Sox-B5', 'P.Br-B0', 'P.Br-B1', 'P.Br-B2',
                              'P.Br-B3', 'P.Br-B4', 'P.Br-B5', 'P.C+-B0', 'P.C+-B1', 'P.C+-B2',
                              'P.C+-B3', 'P.C+-B4', 'P.C+-B5', 'P.C1-B0', 'P.C1-B1', 'P.C1-B2',
                              'P.C1-B3', 'P.C1-B4', 'P.C1-B5', 'P.C2-B0', 'P.C2-B1', 'P.C2-B2',
                              'P.C2-B3', 'P.C2-B4', 'P.C2-B5', 'P.C3-B0', 'P.C3-B1', 'P.C3-B2',
                              'P.C3-B3', 'P.C3-B4', 'P.C3-B5', 'P.Cac-B0', 'P.Cac-B1',
                              'P.Cac-B2', 'P.Cac-B3', 'P.Cac-B4', 'P.Cac-B5', 'P.Car-B0',
                              'P.Car-B1', 'P.Car-B2', 'P.Car-B3', 'P.Car-B4', 'P.Car-B5',
                              'P.Cl-B0', 'P.Cl-B1', 'P.Cl-B2', 'P.Cl-B3', 'P.Cl-B4', 'P.Cl-B5',
                              'P.F-B0', 'P.F-B1', 'P.F-B2', 'P.F-B3', 'P.F-B4', 'P.F-B5',
                              'P.I-B0', 'P.I-B1', 'P.I-B2', 'P.I-B3', 'P.I-B4', 'P.I-B5',
                              'P.N1-B0', 'P.N1-B1', 'P.N1-B2', 'P.N1-B3', 'P.N1-B4', 'P.N1-B5',
                              'P.N2-B0', 'P.N2-B1', 'P.N2-B2', 'P.N2-B3', 'P.N2-B4', 'P.N2-B5',
                              'P.N3-B0', 'P.N3-B1', 'P.N3-B2', 'P.N3-B3', 'P.N3-B4', 'P.N3-B5',
                              'P.N3+-B0', 'P.N3+-B1', 'P.N3+-B2', 'P.N3+-B3', 'P.N3+-B4',
                              'P.N3+-B5', 'P.Nam-B0', 'P.Nam-B1', 'P.Nam-B2', 'P.Nam-B3',
                              'P.Nam-B4', 'P.Nam-B5', 'P.Nar-B0', 'P.Nar-B1', 'P.Nar-B2',
                              'P.Nar-B3', 'P.Nar-B4', 'P.Nar-B5', 'P.Ng+-B0', 'P.Ng+-B1',
                              'P.Ng+-B2', 'P.Ng+-B3', 'P.Ng+-B4', 'P.Ng+-B5', 'P.Npl-B0',
                              'P.Npl-B1', 'P.Npl-B2', 'P.Npl-B3', 'P.Npl-B4', 'P.Npl-B5',
                              'P.Ntr-B0', 'P.Ntr-B1', 'P.Ntr-B2', 'P.Ntr-B3', 'P.Ntr-B4',
                              'P.Ntr-B5', 'P.O--B0', 'P.O--B1', 'P.O--B2', 'P.O--B3', 'P.O--B4',
                              'P.O--B5', 'P.O.co2-B0', 'P.O.co2-B1', 'P.O.co2-B2', 'P.O.co2-B3',
                              'P.O.co2-B4', 'P.O.co2-B5', 'P.O2-B0', 'P.O2-B1', 'P.O2-B2',
                              'P.O2-B3', 'P.O2-B4', 'P.O2-B5', 'P.O3-B0', 'P.O3-B1', 'P.O3-B2',
                              'P.O3-B3', 'P.O3-B4', 'P.O3-B5', 'P.P-B0', 'P.P-B1', 'P.P-B2',
                              'P.P-B3', 'P.P-B4', 'P.P-B5', 'P.Pac-B0', 'P.Pac-B1', 'P.Pac-B2',
                              'P.Pac-B3', 'P.Pac-B4', 'P.Pac-B5', 'P.S2-B0', 'P.S2-B1',
                              'P.S2-B2', 'P.S2-B3', 'P.S2-B4', 'P.S2-B5', 'P.S3-B0', 'P.S3-B1',
                              'P.S3-B2', 'P.S3-B3', 'P.S3-B4', 'P.S3-B5', 'P.So2-B0', 'P.So2-B1',
                              'P.So2-B2', 'P.So2-B3', 'P.So2-B4', 'P.So2-B5', 'P.Sox-B0',
                              'P.Sox-B1', 'P.Sox-B2', 'P.Sox-B3', 'P.Sox-B4', 'P.Sox-B5',
                              'Pac.Br-B0', 'Pac.Br-B1', 'Pac.Br-B2', 'Pac.Br-B3', 'Pac.Br-B4',
                              'Pac.Br-B5', 'Pac.C+-B0', 'Pac.C+-B1', 'Pac.C+-B2', 'Pac.C+-B3',
                              'Pac.C+-B4', 'Pac.C+-B5', 'Pac.C1-B0', 'Pac.C1-B1', 'Pac.C1-B2',
                              'Pac.C1-B3', 'Pac.C1-B4', 'Pac.C1-B5', 'Pac.C2-B0', 'Pac.C2-B1',
                              'Pac.C2-B2', 'Pac.C2-B3', 'Pac.C2-B4', 'Pac.C2-B5', 'Pac.C3-B0',
                              'Pac.C3-B1', 'Pac.C3-B2', 'Pac.C3-B3', 'Pac.C3-B4', 'Pac.C3-B5',
                              'Pac.Cac-B0', 'Pac.Cac-B1', 'Pac.Cac-B2', 'Pac.Cac-B3',
                              'Pac.Cac-B4', 'Pac.Cac-B5', 'Pac.Car-B0', 'Pac.Car-B1',
                              'Pac.Car-B2', 'Pac.Car-B3', 'Pac.Car-B4', 'Pac.Car-B5',
                              'Pac.Cl-B0', 'Pac.Cl-B1', 'Pac.Cl-B2', 'Pac.Cl-B3', 'Pac.Cl-B4',
                              'Pac.Cl-B5', 'Pac.F-B0', 'Pac.F-B1', 'Pac.F-B2', 'Pac.F-B3',
                              'Pac.F-B4', 'Pac.F-B5', 'Pac.I-B0', 'Pac.I-B1', 'Pac.I-B2',
                              'Pac.I-B3', 'Pac.I-B4', 'Pac.I-B5', 'Pac.N1-B0', 'Pac.N1-B1',
                              'Pac.N1-B2', 'Pac.N1-B3', 'Pac.N1-B4', 'Pac.N1-B5', 'Pac.N2-B0',
                              'Pac.N2-B1', 'Pac.N2-B2', 'Pac.N2-B3', 'Pac.N2-B4', 'Pac.N2-B5',
                              'Pac.N3-B0', 'Pac.N3-B1', 'Pac.N3-B2', 'Pac.N3-B3', 'Pac.N3-B4',
                              'Pac.N3-B5', 'Pac.N3+-B0', 'Pac.N3+-B1', 'Pac.N3+-B2',
                              'Pac.N3+-B3', 'Pac.N3+-B4', 'Pac.N3+-B5', 'Pac.Nam-B0',
                              'Pac.Nam-B1', 'Pac.Nam-B2', 'Pac.Nam-B3', 'Pac.Nam-B4',
                              'Pac.Nam-B5', 'Pac.Nar-B0', 'Pac.Nar-B1', 'Pac.Nar-B2',
                              'Pac.Nar-B3', 'Pac.Nar-B4', 'Pac.Nar-B5', 'Pac.Ng+-B0',
                              'Pac.Ng+-B1', 'Pac.Ng+-B2', 'Pac.Ng+-B3', 'Pac.Ng+-B4',
                              'Pac.Ng+-B5', 'Pac.Npl-B0', 'Pac.Npl-B1', 'Pac.Npl-B2',
                              'Pac.Npl-B3', 'Pac.Npl-B4', 'Pac.Npl-B5', 'Pac.Ntr-B0',
                              'Pac.Ntr-B1', 'Pac.Ntr-B2', 'Pac.Ntr-B3', 'Pac.Ntr-B4',
                              'Pac.Ntr-B5', 'Pac.O--B0', 'Pac.O--B1', 'Pac.O--B2', 'Pac.O--B3',
                              'Pac.O--B4', 'Pac.O--B5', 'Pac.O.co2-B0', 'Pac.O.co2-B1',
                              'Pac.O.co2-B2', 'Pac.O.co2-B3', 'Pac.O.co2-B4', 'Pac.O.co2-B5',
                              'Pac.O2-B0', 'Pac.O2-B1', 'Pac.O2-B2', 'Pac.O2-B3', 'Pac.O2-B4',
                              'Pac.O2-B5', 'Pac.O3-B0', 'Pac.O3-B1', 'Pac.O3-B2', 'Pac.O3-B3',
                              'Pac.O3-B4', 'Pac.O3-B5', 'Pac.P-B0', 'Pac.P-B1', 'Pac.P-B2',
                              'Pac.P-B3', 'Pac.P-B4', 'Pac.P-B5', 'Pac.Pac-B0', 'Pac.Pac-B1',
                              'Pac.Pac-B2', 'Pac.Pac-B3', 'Pac.Pac-B4', 'Pac.Pac-B5',
                              'Pac.S2-B0', 'Pac.S2-B1', 'Pac.S2-B2', 'Pac.S2-B3', 'Pac.S2-B4',
                              'Pac.S2-B5', 'Pac.S3-B0', 'Pac.S3-B1', 'Pac.S3-B2', 'Pac.S3-B3',
                              'Pac.S3-B4', 'Pac.S3-B5', 'Pac.So2-B0', 'Pac.So2-B1', 'Pac.So2-B2',
                              'Pac.So2-B3', 'Pac.So2-B4', 'Pac.So2-B5', 'Pac.Sox-B0',
                              'Pac.Sox-B1', 'Pac.Sox-B2', 'Pac.Sox-B3', 'Pac.Sox-B4',
                              'Pac.Sox-B5', 'S2.Br-B0', 'S2.Br-B1', 'S2.Br-B2', 'S2.Br-B3',
                              'S2.Br-B4', 'S2.Br-B5', 'S2.C+-B0', 'S2.C+-B1', 'S2.C+-B2',
                              'S2.C+-B3', 'S2.C+-B4', 'S2.C+-B5', 'S2.C1-B0', 'S2.C1-B1',
                              'S2.C1-B2', 'S2.C1-B3', 'S2.C1-B4', 'S2.C1-B5', 'S2.C2-B0',
                              'S2.C2-B1', 'S2.C2-B2', 'S2.C2-B3', 'S2.C2-B4', 'S2.C2-B5',
                              'S2.C3-B0', 'S2.C3-B1', 'S2.C3-B2', 'S2.C3-B3', 'S2.C3-B4',
                              'S2.C3-B5', 'S2.Cac-B0', 'S2.Cac-B1', 'S2.Cac-B2', 'S2.Cac-B3',
                              'S2.Cac-B4', 'S2.Cac-B5', 'S2.Car-B0', 'S2.Car-B1', 'S2.Car-B2',
                              'S2.Car-B3', 'S2.Car-B4', 'S2.Car-B5', 'S2.Cl-B0', 'S2.Cl-B1',
                              'S2.Cl-B2', 'S2.Cl-B3', 'S2.Cl-B4', 'S2.Cl-B5', 'S2.F-B0',
                              'S2.F-B1', 'S2.F-B2', 'S2.F-B3', 'S2.F-B4', 'S2.F-B5', 'S2.I-B0',
                              'S2.I-B1', 'S2.I-B2', 'S2.I-B3', 'S2.I-B4', 'S2.I-B5', 'S2.N1-B0',
                              'S2.N1-B1', 'S2.N1-B2', 'S2.N1-B3', 'S2.N1-B4', 'S2.N1-B5',
                              'S2.N2-B0', 'S2.N2-B1', 'S2.N2-B2', 'S2.N2-B3', 'S2.N2-B4',
                              'S2.N2-B5', 'S2.N3-B0', 'S2.N3-B1', 'S2.N3-B2', 'S2.N3-B3',
                              'S2.N3-B4', 'S2.N3-B5', 'S2.N3+-B0', 'S2.N3+-B1', 'S2.N3+-B2',
                              'S2.N3+-B3', 'S2.N3+-B4', 'S2.N3+-B5', 'S2.Nam-B0', 'S2.Nam-B1',
                              'S2.Nam-B2', 'S2.Nam-B3', 'S2.Nam-B4', 'S2.Nam-B5', 'S2.Nar-B0',
                              'S2.Nar-B1', 'S2.Nar-B2', 'S2.Nar-B3', 'S2.Nar-B4', 'S2.Nar-B5',
                              'S2.Ng+-B0', 'S2.Ng+-B1', 'S2.Ng+-B2', 'S2.Ng+-B3', 'S2.Ng+-B4',
                              'S2.Ng+-B5', 'S2.Npl-B0', 'S2.Npl-B1', 'S2.Npl-B2', 'S2.Npl-B3',
                              'S2.Npl-B4', 'S2.Npl-B5', 'S2.Ntr-B0', 'S2.Ntr-B1', 'S2.Ntr-B2',
                              'S2.Ntr-B3', 'S2.Ntr-B4', 'S2.Ntr-B5', 'S2.O--B0', 'S2.O--B1',
                              'S2.O--B2', 'S2.O--B3', 'S2.O--B4', 'S2.O--B5', 'S2.O.co2-B0',
                              'S2.O.co2-B1', 'S2.O.co2-B2', 'S2.O.co2-B3', 'S2.O.co2-B4',
                              'S2.O.co2-B5', 'S2.O2-B0', 'S2.O2-B1', 'S2.O2-B2', 'S2.O2-B3',
                              'S2.O2-B4', 'S2.O2-B5', 'S2.O3-B0', 'S2.O3-B1', 'S2.O3-B2',
                              'S2.O3-B3', 'S2.O3-B4', 'S2.O3-B5', 'S2.P-B0', 'S2.P-B1',
                              'S2.P-B2', 'S2.P-B3', 'S2.P-B4', 'S2.P-B5', 'S2.Pac-B0',
                              'S2.Pac-B1', 'S2.Pac-B2', 'S2.Pac-B3', 'S2.Pac-B4', 'S2.Pac-B5',
                              'S2.S2-B0', 'S2.S2-B1', 'S2.S2-B2', 'S2.S2-B3', 'S2.S2-B4',
                              'S2.S2-B5', 'S2.S3-B0', 'S2.S3-B1', 'S2.S3-B2', 'S2.S3-B3',
                              'S2.S3-B4', 'S2.S3-B5', 'S2.So2-B0', 'S2.So2-B1', 'S2.So2-B2',
                              'S2.So2-B3', 'S2.So2-B4', 'S2.So2-B5', 'S2.Sox-B0', 'S2.Sox-B1',
                              'S2.Sox-B2', 'S2.Sox-B3', 'S2.Sox-B4', 'S2.Sox-B5', 'S3.Br-B0',
                              'S3.Br-B1', 'S3.Br-B2', 'S3.Br-B3', 'S3.Br-B4', 'S3.Br-B5',
                              'S3.C+-B0', 'S3.C+-B1', 'S3.C+-B2', 'S3.C+-B3', 'S3.C+-B4',
                              'S3.C+-B5', 'S3.C1-B0', 'S3.C1-B1', 'S3.C1-B2', 'S3.C1-B3',
                              'S3.C1-B4', 'S3.C1-B5', 'S3.C2-B0', 'S3.C2-B1', 'S3.C2-B2',
                              'S3.C2-B3', 'S3.C2-B4', 'S3.C2-B5', 'S3.C3-B0', 'S3.C3-B1',
                              'S3.C3-B2', 'S3.C3-B3', 'S3.C3-B4', 'S3.C3-B5', 'S3.Cac-B0',
                              'S3.Cac-B1', 'S3.Cac-B2', 'S3.Cac-B3', 'S3.Cac-B4', 'S3.Cac-B5',
                              'S3.Car-B0', 'S3.Car-B1', 'S3.Car-B2', 'S3.Car-B3', 'S3.Car-B4',
                              'S3.Car-B5', 'S3.Cl-B0', 'S3.Cl-B1', 'S3.Cl-B2', 'S3.Cl-B3',
                              'S3.Cl-B4', 'S3.Cl-B5', 'S3.F-B0', 'S3.F-B1', 'S3.F-B2', 'S3.F-B3',
                              'S3.F-B4', 'S3.F-B5', 'S3.I-B0', 'S3.I-B1', 'S3.I-B2', 'S3.I-B3',
                              'S3.I-B4', 'S3.I-B5', 'S3.N1-B0', 'S3.N1-B1', 'S3.N1-B2',
                              'S3.N1-B3', 'S3.N1-B4', 'S3.N1-B5', 'S3.N2-B0', 'S3.N2-B1',
                              'S3.N2-B2', 'S3.N2-B3', 'S3.N2-B4', 'S3.N2-B5', 'S3.N3-B0',
                              'S3.N3-B1', 'S3.N3-B2', 'S3.N3-B3', 'S3.N3-B4', 'S3.N3-B5',
                              'S3.N3+-B0', 'S3.N3+-B1', 'S3.N3+-B2', 'S3.N3+-B3', 'S3.N3+-B4',
                              'S3.N3+-B5', 'S3.Nam-B0', 'S3.Nam-B1', 'S3.Nam-B2', 'S3.Nam-B3',
                              'S3.Nam-B4', 'S3.Nam-B5', 'S3.Nar-B0', 'S3.Nar-B1', 'S3.Nar-B2',
                              'S3.Nar-B3', 'S3.Nar-B4', 'S3.Nar-B5', 'S3.Ng+-B0', 'S3.Ng+-B1',
                              'S3.Ng+-B2', 'S3.Ng+-B3', 'S3.Ng+-B4', 'S3.Ng+-B5', 'S3.Npl-B0',
                              'S3.Npl-B1', 'S3.Npl-B2', 'S3.Npl-B3', 'S3.Npl-B4', 'S3.Npl-B5',
                              'S3.Ntr-B0', 'S3.Ntr-B1', 'S3.Ntr-B2', 'S3.Ntr-B3', 'S3.Ntr-B4',
                              'S3.Ntr-B5', 'S3.O--B0', 'S3.O--B1', 'S3.O--B2', 'S3.O--B3',
                              'S3.O--B4', 'S3.O--B5', 'S3.O.co2-B0', 'S3.O.co2-B1',
                              'S3.O.co2-B2', 'S3.O.co2-B3', 'S3.O.co2-B4', 'S3.O.co2-B5',
                              'S3.O2-B0', 'S3.O2-B1', 'S3.O2-B2', 'S3.O2-B3', 'S3.O2-B4',
                              'S3.O2-B5', 'S3.O3-B0', 'S3.O3-B1', 'S3.O3-B2', 'S3.O3-B3',
                              'S3.O3-B4', 'S3.O3-B5', 'S3.P-B0', 'S3.P-B1', 'S3.P-B2', 'S3.P-B3',
                              'S3.P-B4', 'S3.P-B5', 'S3.Pac-B0', 'S3.Pac-B1', 'S3.Pac-B2',
                              'S3.Pac-B3', 'S3.Pac-B4', 'S3.Pac-B5', 'S3.S2-B0', 'S3.S2-B1',
                              'S3.S2-B2', 'S3.S2-B3', 'S3.S2-B4', 'S3.S2-B5', 'S3.S3-B0',
                              'S3.S3-B1', 'S3.S3-B2', 'S3.S3-B3', 'S3.S3-B4', 'S3.S3-B5',
                              'S3.So2-B0', 'S3.So2-B1', 'S3.So2-B2', 'S3.So2-B3', 'S3.So2-B4',
                              'S3.So2-B5', 'S3.Sox-B0', 'S3.Sox-B1', 'S3.Sox-B2', 'S3.Sox-B3',
                              'S3.Sox-B4', 'S3.Sox-B5', 'So2.Br-B0', 'So2.Br-B1', 'So2.Br-B2',
                              'So2.Br-B3', 'So2.Br-B4', 'So2.Br-B5', 'So2.C+-B0', 'So2.C+-B1',
                              'So2.C+-B2', 'So2.C+-B3', 'So2.C+-B4', 'So2.C+-B5', 'So2.C1-B0',
                              'So2.C1-B1', 'So2.C1-B2', 'So2.C1-B3', 'So2.C1-B4', 'So2.C1-B5',
                              'So2.C2-B0', 'So2.C2-B1', 'So2.C2-B2', 'So2.C2-B3', 'So2.C2-B4',
                              'So2.C2-B5', 'So2.C3-B0', 'So2.C3-B1', 'So2.C3-B2', 'So2.C3-B3',
                              'So2.C3-B4', 'So2.C3-B5', 'So2.Cac-B0', 'So2.Cac-B1', 'So2.Cac-B2',
                              'So2.Cac-B3', 'So2.Cac-B4', 'So2.Cac-B5', 'So2.Car-B0',
                              'So2.Car-B1', 'So2.Car-B2', 'So2.Car-B3', 'So2.Car-B4',
                              'So2.Car-B5', 'So2.Cl-B0', 'So2.Cl-B1', 'So2.Cl-B2', 'So2.Cl-B3',
                              'So2.Cl-B4', 'So2.Cl-B5', 'So2.F-B0', 'So2.F-B1', 'So2.F-B2',
                              'So2.F-B3', 'So2.F-B4', 'So2.F-B5', 'So2.I-B0', 'So2.I-B1',
                              'So2.I-B2', 'So2.I-B3', 'So2.I-B4', 'So2.I-B5', 'So2.N1-B0',
                              'So2.N1-B1', 'So2.N1-B2', 'So2.N1-B3', 'So2.N1-B4', 'So2.N1-B5',
                              'So2.N2-B0', 'So2.N2-B1', 'So2.N2-B2', 'So2.N2-B3', 'So2.N2-B4',
                              'So2.N2-B5', 'So2.N3-B0', 'So2.N3-B1', 'So2.N3-B2', 'So2.N3-B3',
                              'So2.N3-B4', 'So2.N3-B5', 'So2.N3+-B0', 'So2.N3+-B1', 'So2.N3+-B2',
                              'So2.N3+-B3', 'So2.N3+-B4', 'So2.N3+-B5', 'So2.Nam-B0',
                              'So2.Nam-B1', 'So2.Nam-B2', 'So2.Nam-B3', 'So2.Nam-B4',
                              'So2.Nam-B5', 'So2.Nar-B0', 'So2.Nar-B1', 'So2.Nar-B2',
                              'So2.Nar-B3', 'So2.Nar-B4', 'So2.Nar-B5', 'So2.Ng+-B0',
                              'So2.Ng+-B1', 'So2.Ng+-B2', 'So2.Ng+-B3', 'So2.Ng+-B4',
                              'So2.Ng+-B5', 'So2.Npl-B0', 'So2.Npl-B1', 'So2.Npl-B2',
                              'So2.Npl-B3', 'So2.Npl-B4', 'So2.Npl-B5', 'So2.Ntr-B0',
                              'So2.Ntr-B1', 'So2.Ntr-B2', 'So2.Ntr-B3', 'So2.Ntr-B4',
                              'So2.Ntr-B5', 'So2.O--B0', 'So2.O--B1', 'So2.O--B2', 'So2.O--B3',
                              'So2.O--B4', 'So2.O--B5', 'So2.O.co2-B0', 'So2.O.co2-B1',
                              'So2.O.co2-B2', 'So2.O.co2-B3', 'So2.O.co2-B4', 'So2.O.co2-B5',
                              'So2.O2-B0', 'So2.O2-B1', 'So2.O2-B2', 'So2.O2-B3', 'So2.O2-B4',
                              'So2.O2-B5', 'So2.O3-B0', 'So2.O3-B1', 'So2.O3-B2', 'So2.O3-B3',
                              'So2.O3-B4', 'So2.O3-B5', 'So2.P-B0', 'So2.P-B1', 'So2.P-B2',
                              'So2.P-B3', 'So2.P-B4', 'So2.P-B5', 'So2.Pac-B0', 'So2.Pac-B1',
                              'So2.Pac-B2', 'So2.Pac-B3', 'So2.Pac-B4', 'So2.Pac-B5',
                              'So2.S2-B0', 'So2.S2-B1', 'So2.S2-B2', 'So2.S2-B3', 'So2.S2-B4',
                              'So2.S2-B5', 'So2.S3-B0', 'So2.S3-B1', 'So2.S3-B2', 'So2.S3-B3',
                              'So2.S3-B4', 'So2.S3-B5', 'So2.So2-B0', 'So2.So2-B1', 'So2.So2-B2',
                              'So2.So2-B3', 'So2.So2-B4', 'So2.So2-B5', 'So2.Sox-B0',
                              'So2.Sox-B1', 'So2.Sox-B2', 'So2.Sox-B3', 'So2.Sox-B4',
                              'So2.Sox-B5', 'Sox.Br-B0', 'Sox.Br-B1', 'Sox.Br-B2', 'Sox.Br-B3',
                              'Sox.Br-B4', 'Sox.Br-B5', 'Sox.C+-B0', 'Sox.C+-B1', 'Sox.C+-B2',
                              'Sox.C+-B3', 'Sox.C+-B4', 'Sox.C+-B5', 'Sox.C1-B0', 'Sox.C1-B1',
                              'Sox.C1-B2', 'Sox.C1-B3', 'Sox.C1-B4', 'Sox.C1-B5', 'Sox.C2-B0',
                              'Sox.C2-B1', 'Sox.C2-B2', 'Sox.C2-B3', 'Sox.C2-B4', 'Sox.C2-B5',
                              'Sox.C3-B0', 'Sox.C3-B1', 'Sox.C3-B2', 'Sox.C3-B3', 'Sox.C3-B4',
                              'Sox.C3-B5', 'Sox.Cac-B0', 'Sox.Cac-B1', 'Sox.Cac-B2',
                              'Sox.Cac-B3', 'Sox.Cac-B4', 'Sox.Cac-B5', 'Sox.Car-B0',
                              'Sox.Car-B1', 'Sox.Car-B2', 'Sox.Car-B3', 'Sox.Car-B4',
                              'Sox.Car-B5', 'Sox.Cl-B0', 'Sox.Cl-B1', 'Sox.Cl-B2', 'Sox.Cl-B3',
                              'Sox.Cl-B4', 'Sox.Cl-B5', 'Sox.F-B0', 'Sox.F-B1', 'Sox.F-B2',
                              'Sox.F-B3', 'Sox.F-B4', 'Sox.F-B5', 'Sox.I-B0', 'Sox.I-B1',
                              'Sox.I-B2', 'Sox.I-B3', 'Sox.I-B4', 'Sox.I-B5', 'Sox.N1-B0',
                              'Sox.N1-B1', 'Sox.N1-B2', 'Sox.N1-B3', 'Sox.N1-B4', 'Sox.N1-B5',
                              'Sox.N2-B0', 'Sox.N2-B1', 'Sox.N2-B2', 'Sox.N2-B3', 'Sox.N2-B4',
                              'Sox.N2-B5', 'Sox.N3-B0', 'Sox.N3-B1', 'Sox.N3-B2', 'Sox.N3-B3',
                              'Sox.N3-B4', 'Sox.N3-B5', 'Sox.N3+-B0', 'Sox.N3+-B1', 'Sox.N3+-B2',
                              'Sox.N3+-B3', 'Sox.N3+-B4', 'Sox.N3+-B5', 'Sox.Nam-B0',
                              'Sox.Nam-B1', 'Sox.Nam-B2', 'Sox.Nam-B3', 'Sox.Nam-B4',
                              'Sox.Nam-B5', 'Sox.Nar-B0', 'Sox.Nar-B1', 'Sox.Nar-B2',
                              'Sox.Nar-B3', 'Sox.Nar-B4', 'Sox.Nar-B5', 'Sox.Ng+-B0',
                              'Sox.Ng+-B1', 'Sox.Ng+-B2', 'Sox.Ng+-B3', 'Sox.Ng+-B4',
                              'Sox.Ng+-B5', 'Sox.Npl-B0', 'Sox.Npl-B1', 'Sox.Npl-B2',
                              'Sox.Npl-B3', 'Sox.Npl-B4', 'Sox.Npl-B5', 'Sox.Ntr-B0',
                              'Sox.Ntr-B1', 'Sox.Ntr-B2', 'Sox.Ntr-B3', 'Sox.Ntr-B4',
                              'Sox.Ntr-B5', 'Sox.O--B0', 'Sox.O--B1', 'Sox.O--B2', 'Sox.O--B3',
                              'Sox.O--B4', 'Sox.O--B5', 'Sox.O.co2-B0', 'Sox.O.co2-B1',
                              'Sox.O.co2-B2', 'Sox.O.co2-B3', 'Sox.O.co2-B4', 'Sox.O.co2-B5',
                              'Sox.O2-B0', 'Sox.O2-B1', 'Sox.O2-B2', 'Sox.O2-B3', 'Sox.O2-B4',
                              'Sox.O2-B5', 'Sox.O3-B0', 'Sox.O3-B1', 'Sox.O3-B2', 'Sox.O3-B3',
                              'Sox.O3-B4', 'Sox.O3-B5', 'Sox.P-B0', 'Sox.P-B1', 'Sox.P-B2',
                              'Sox.P-B3', 'Sox.P-B4', 'Sox.P-B5', 'Sox.Pac-B0', 'Sox.Pac-B1',
                              'Sox.Pac-B2', 'Sox.Pac-B3', 'Sox.Pac-B4', 'Sox.Pac-B5',
                              'Sox.S2-B0', 'Sox.S2-B1', 'Sox.S2-B2', 'Sox.S2-B3', 'Sox.S2-B4',
                              'Sox.S2-B5', 'Sox.S3-B0', 'Sox.S3-B1', 'Sox.S3-B2', 'Sox.S3-B3',
                              'Sox.S3-B4', 'Sox.S3-B5', 'Sox.So2-B0', 'Sox.So2-B1', 'Sox.So2-B2',
                              'Sox.So2-B3', 'Sox.So2-B4', 'Sox.So2-B5', 'Sox.Sox-B0',
                              'Sox.Sox-B1', 'Sox.Sox-B2', 'Sox.Sox-B3', 'Sox.Sox-B4',
                              'Sox.Sox-B5'],
            'sasa': ['name', 'SASA_Score', 'SASA_com_lig_sasa_tot', 'SASA_lig_sasa_tot', 'SASA_phobic_lig_exposed',
                     'SASA_com_lig_sasa_phobic', 'SASA_lig_sasa_phobic', 'SASA_philic_lig_exposed',
                     'SASA_com_lig_sasa_philic', 'SASA_lig_sasa_philic', 'SASA_other_lig_exposed',
                     'SASA_com_lig_sasa_other',
                     'SASA_lig_sasa_other', 'SASA_internal_repulsive'],
            'smina': ['name', 'gauss(o=0,_w=0.3,_c=8)', 'gauss(o=0.5,_w=0.3,_c=8)', 'gauss(o=1,_w=0.3,_c=8)',
             'gauss(o=1.5,_w=0.3,_c=8)', 'gauss(o=2,_w=0.3,_c=8)', 'gauss(o=2.5,_w=0.3,_c=8)', 'gauss(o=0,_w=0.5,_c=8)'
                , 'gauss(o=1,_w=0.5,_c=8)', 'gauss(o=2,_w=0.5,_c=8)', 'gauss(o=0,_w=0.7,_c=8)',
             'gauss(o=1,_w=0.7,_c=8)',
             'gauss(o=2,_w=0.7,_c=8)', 'gauss(o=0,_w=0.9,_c=8)', 'gauss(o=1,_w=0.9,_c=8)', 'gauss(o=2,_w=0.9,_c=8)',
             'gauss(o=3,_w=0.9,_c=8)', 'gauss(o=0,_w=1.5,_c=8)', 'gauss(o=1,_w=1.5,_c=8)', 'gauss(o=2,_w=1.5,_c=8)',
             'gauss(o=3,_w=1.5,_c=8)', 'gauss(o=4,_w=1.5,_c=8)', 'gauss(o=0,_w=2,_c=8)', 'gauss(o=1,_w=2,_c=8)',
             'gauss(o=2,_w=2,_c=8)', 'gauss(o=3,_w=2,_c=8)', 'gauss(o=4,_w=2,_c=8)', 'gauss(o=0,_w=3,_c=8)',
             'gauss(o=1,_w=3,_c=8)', 'gauss(o=2,_w=3,_c=8)', 'gauss(o=3,_w=3,_c=8)', 'gauss(o=4,_w=3,_c=8)',
             'repulsion(o=0.4,_c=8)', 'repulsion(o=0.2,_c=8)', 'repulsion(o=0,_c=8)', 'repulsion(o=-0.2,_c=8)',
             'repulsion(o=-0.4,_c=8)', 'repulsion(o=-0.6,_c=8)', 'repulsion(o=-0.8,_c=8)', 'repulsion(o=-1,_c=8)',
             'hydrophobic(g=0.5,_b=1.5,_c=8)', 'hydrophobic(g=0.5,_b=1,_c=8)', 'hydrophobic(g=0.5,_b=2,_c=8)',
             'hydrophobic(g=0.5,_b=3,_c=8)', 'non_hydrophobic(g=0.5,_b=1.5,_c=8)', 'vdw(i=4,_j=8,_s=0,_^=100,_c=8)',
             'vdw(i=6,_j=12,_s=1,_^=100,_c=8)', 'non_dir_h_bond(g=-0.7,_b=0,_c=8)', 'non_dir_h_bond(g=-0.7,_b=0.2,_c=8)'
                , 'non_dir_h_bond(g=-0.7,_b=0.5,_c=8)', 'non_dir_h_bond(g=-1,_b=0,_c=8)',
             'non_dir_h_bond(g=-1,_b=0.2,_c=8)',
             'non_dir_h_bond(g=-1,_b=0.5,_c=8)', 'non_dir_h_bond(g=-1.3,_b=0,_c=8)',
             'non_dir_h_bond(g=-1.3,_b=0.2,_c=8)',
             'non_dir_h_bond(g=-1.3,_b=0.5,_c=8)', 'non_dir_anti_h_bond_quadratic(o=0,_c=8)',
             'non_dir_anti_h_bond_quadratic(o=0.5,_c=8)', 'non_dir_anti_h_bond_quadratic(o=1,_c=8)',
             'donor_donor_quadratic(o=0,_c=8)', 'donor_donor_quadratic(o=0.5,_c=8)', 'donor_donor_quadratic(o=1,_c=8)',
             'acceptor_acceptor_quadratic(o=0,_c=8)', 'acceptor_acceptor_quadratic(o=0.5,_c=8)',
             'acceptor_acceptor_quadratic(o=1,_c=8)', 'non_dir_h_bond_lj(o=-0.7,_^=100,_c=8)',
             'non_dir_h_bond_lj(o=-1,_^=100,_c=8)', 'non_dir_h_bond_lj(o=-1.3,_^=100,_c=8)',
             'ad4_solvation(d-sigma=3.6,_s/q=0.01097,_c=8)', 'ad4_solvation(d-sigma=3.6,_s/q=0,_c=8)',
             'electrostatic(i=1,_^=100,_c=8)', 'electrostatic(i=2,_^=100,_c=8)', 'num_tors_div', 'num_tors_div_simple',
             'num_heavy_atoms_div', 'num_heavy_atoms', 'num_tors_add', 'num_tors_sqr', 'num_tors_sqrt',
             'num_hydrophobic_atoms', 'ligand_length', 'num_ligands'],
            'smog2016': ['name', 'total', 'SMoG2016_KBP2016', 'SMoG2016_LJP', 'SMoG2016_Rotor', 'SMoG2016_lnMass'],
            'surflex': ['name', 'surflex_rot_bond', 'surflex_score', 'surflex_crash', 'surflex_self-crash', 'surflex_polar'],
            'sp': ['name', 'rotatable_bonds', 'dockingscore', 'ecoul', 'eff_state_penalty', 'einternal', 'emodel',
                   'energy',
                   'erotb', 'esite', 'evdw', 'gscore', 'hbond', 'ligand_efficiency', 'ligand_efficiency_ln',
                   'ligand_efficiency_sa',
                   'lipo', 'metal', 'rewards'],
            'xp': ['name', 'XP_rotatable_bonds', 'XP_ClBr', 'XP_Electro', 'XP_ExposPenal', 'XP_GScore', 'XP_HBPenal',
                   'XP_HBond',
                   'XP_LipophilicEvdW', 'XP_LowMW', 'XP_Penalties', 'XP_PhobEn', 'XP_PhobEnHB', 'XP_PhobEnPairHB',
                   'XP_PiCat', 'XP_PiStack', 'XP_RotPenal', 'XP_Sitemap', 'XP_Zpotr', 'XP_dockingscore', 'XP_ecoul',
                   'XP_eff_state_penalty',
                   'XP_einternal', 'XP_emodel', 'XP_energy', 'XP_evdw', 'XP_gscore', 'XP_ligand_efficiency',
                   'XP_ligand_efficiency_ln', 'XP_ligand_efficiency_sa'],
            'vina': ['name', 'vina_score', 'gauss1', 'gauss2', 'repulsion', 'hydrophobic', 'HB', 'rbond'],
            'xscore': ['name', 'VDW', 'HB', 'HP', 'HM', 'HS', 'RT', 'XSCORE']
        }


    def cal_affiscore(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/affiscore'.format(path_job)  # 目标路径
        # 蛋白文件
        pre_protein_file = '{}/protein_file.pdb'.format(dst_path)
        # 小分子预处理
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)
        lig_temp = '{}/{}_temp.mol2'.format(dst_path, lig_name)
        # 小分子预处理
        self.slide_lig_pre(ligand, lig_temp)
        # 定义log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 计算
        cmdline = 'export SLIDE_DIR=/home/xujun/Soft/Score_Function/SLIDE/SLIDE-master &&'
        cmdline += '${SLIDE_DIR}/bin/slide_score -p %s -l %s > %s' % (pre_protein_file, lig_temp, log_file)
        os.system(cmdline)

    def cal_autodock(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/autodock'.format(path_job)  # 目标路径
        # 蛋白文件
        pre_protein_file = '{}/protein_file.pdbqt'.format(dst_path)
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)
        lig_pdbqt = '{}/{}.pdbqt'.format(dst_path, lig_name)
        # 定义log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 准备小分子并计算
        cmdline = 'module purge &&'
        cmdline += 'module load autodock &&'
        cmdline += 'prepare_ligand4.py -l {} -o {} -A hydrogens &&'.format(ligand, lig_pdbqt)
        cmdline += 'compute_AutoDock41_score.py -r {} -l {} -o {}'.format(pre_protein_file, lig_pdbqt, log_file)
        os.system(cmdline)

    def cal_chemgauss(self, ligand):
        # 还原路径
        path_job = self.path_job # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/chemgauss'.format(path_job)  # 目标路径
        # 蛋白文件
        pre_protein_file = '{}/protein_file.oeb'.format(dst_path)
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)
        lig_oeb = '{}/{}.oeb'.format(dst_path, lig_name)
        # 定义log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 准备小分子
        cmd = 'module load openeye &&convert.py {} {} '.format(ligand, lig_oeb)
        os.system(cmd)
        # 计算
        cmd = 'cd {}&&export PATH=$PATH:/opt/openeye/applications-2018.11.3/bin &&fred -receptor {} -dbase {} ' \
              '-dock_resolution High -save_component_scores -score_file {} -hitlist_size 1'.format(
            dst_path, pre_protein_file, lig_oeb, log_file)
        os.system(cmd)

    def cal_cyscore(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/cyscore'.format(path_job)  # 目标路径
        # 蛋白文件
        protein_file = '{}/protein_file.pdb'.format(path_files)
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)
        # 定义log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 计算
        cmd = 'cd {}&&export PATH=${{PATH}}:/home/xujun/Soft/Score_Function/Cyscore/CyscoreV2.0/bin &&' \
              'Cyscore {} {} > {}'.format(dst_path, protein_file, ligand, log_file)
        os.system(cmd)

    def cal_dsx(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig # 分子存放路径
        dst_path = '{}/dsx'.format(path_job)  # 目标路径
        # 蛋白文件
        protein_file = self.protein
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)
        # 循环计算三种能量项
        for i in [0, 2, 3]:
            # 定义log文件
            log_file = '{}/{}_{}.txt'.format(dst_path, lig_name, i)
            # 定义打分项
            score_type = ''.join(['-T%s 1 ' % i] + ['-T%s 0 ' % j for j in range(5) if j != i])
            # 计算
            cmd = 'cd {}&&/home/xujun/Soft/Score_Function/dsx/dsx090_and_hotspotsx061_linux/linux64/dsx_linux_64.lnx ' \
                  '-D /home/xujun/Soft/Score_Function/dsx/dsx090_and_hotspotsx061_linux/pdb_pot_0511 ' \
                  '-P {} -L {} -F {} '.format(dst_path, protein_file, ligand, log_file) + score_type
            os.system(cmd)

    def cal_galaxy(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/galaxy'.format(path_job)  # 目标路径
        # 蛋白文件
        pre_protein_file = '{}/protein_file.pdb'.format(dst_path)
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)
        dst_lig_path = '{}/{}'.format(dst_path, lig_name)
        # 创建小分子文件夹
        os.mkdir(dst_lig_path)
        # 计算  要进入到对应文件夹内，log文件才会生成在那个文件夹中
        cmd = 'cd {}&&module load anaconda2&&/home/xujun/Soft/Score_Function/Galaxy/GalaxyDock_BP2/script/calc_energy.py' \
              ' -d /home/xujun/Soft/Score_Function/Galaxy/GalaxyDock_BP2 -p {} -l {}'.format(dst_lig_path,
                                                                                             pre_protein_file, ligand)
        os.system(cmd)

    def cal_nnscore(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/nnscore'.format(path_job)  # 目标路径
        # 蛋白文件
        pre_protein_file = '{}/protein_file.pdbqt'.format(dst_path)
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)  # 源配体
        pre_ligand = '{}/{}.pdbqt'.format(dst_path, lig_name)  # 准备后配体
        # log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 准备小分子
        self.vina_lig_prep(ligand, pre_ligand)  # 加H，转换为pdbqtg格式
        if not os.path.exists(pre_ligand):
            return None 
        with open(pre_ligand, 'r')as f:  # 读取文件
            lines = f.readlines()
        new_lines = []
        for line in lines:
            if line.startswith('ATOM'):
                new_lines.append(line[:23] + '   ' + line[26:])
            else:
                new_lines.append(line)
        new_lig = ''.join(new_lines)  # 合并字符串
        with open(pre_ligand, 'w')as f:  # 生成新文件
            f.write(new_lig)
        # 计算
        # 辅助py
        cal_nn = '{}/cal_nn.py'.format(self.help_path)
        # 指令
        cmd = 'export PYTHONPATH=/home/xujun/Soft/Score_Function/NNscore:${{PYTHONPATH}}&& timeout 60 python2 {} {} {} {}'.format(
            cal_nn, pre_protein_file, pre_ligand, log_file)
        # 执行
        os.system(cmd)

    def cal_smina(self, ligand, x, y, z):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/smina'.format(path_job)  # 目标路径
        # 蛋白文件
        pre_protein_file = '{}/protein_file.pdbqt'.format(dst_path)
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)  # 源配体
        pre_ligand = '{}/{}.pdbqt'.format(dst_path, lig_name)  # 准备后配体
        # log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 准备小分子
        self.vina_lig_prep(ligand, pre_ligand)  # 加H，转换为pdbqtg格式
        # 计算 
        # 指令 
        cmd = 'timeout 60 /home/xujun/Soft/Score_Function/smina/smina --receptor {} --ligand {} ' \
                   '--center_x {} --center_y {} --center_z {} --size_x 18.75 --size_y 18.75 --size_z 18.75 --log {} ' \
                  '--custom_scoring /home/xujun/Soft/Score_Function/smina/total.score '\
                   '--score_only --cpu 1'.format(pre_protein_file, pre_ligand, x, y, z, log_file)
        # 执行
        os.system(cmd)

    def cal_smog2016(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/smog2016'.format(path_job)  # 目标路径
        # 蛋白文件
        protein_file = self.protein
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)  # 源配体
        # log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 计算
        # 指令
        cmd = 'module load openbabel&&cd /home/xujun/Soft/Score_Function/SMoGxxx &&chmod +x SMoG2016.exe &&' \
              '/home/xujun/Soft/Score_Function/SMoGxxx/SMoG2016.exe {} {} 0 > {}'.format(protein_file, ligand, log_file)
        # 执行
        os.system(cmd)

    def cal_surflex(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/surflex'.format(path_job)  # 目标路径
        # 蛋白文件
        pre_protein = '{}/protein_file.mol2'.format(dst_path)
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)  # 源配体
        # 分子文件夹
        dst_lig_path = '{}/{}'.format(dst_path, lig_name)
        os.mkdir(dst_lig_path)
        # list文件
        ligand_list = '{}/list'.format(dst_lig_path)
        with open(ligand_list, 'w') as f:  # 创建list文件
            f.write('{}'.format(ligand))
        # log文件
        log_file = '{}/{}/loggeom'.format(dst_path, lig_name)
        # 计算
        # 指令
        cmd = 'cd {}&&module purge &&module load tripos/sybylx2.1.1 &&surflex-dock.exe +misc_premin score_list ./list  ../p1-protomol.mol2 {} {}'.format(dst_lig_path, pre_protein, log_file)
        # 执行
        os.system(cmd)
        # 删除多余文件
        os.remove(ligand_list)  # list文件  76B
        os.remove('{}/scores'.format(dst_lig_path))  # scores  2KB
        # os.remove('{}/p1-protomol.mol2'.format(dst_path))  # p1-protomol.mol2  27KB
        os.remove('{}/loggeom-results.mol2'.format(dst_lig_path))  # loggeom-results.mol2 6KB

    def cal_vina(self, ligand, x, y, z):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/vina'.format(path_job)  # 目标路径
        # 蛋白文件
        pre_protein_file = '{}/protein_file.pdbqt'.format(dst_path)
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)  # 源配体
        pre_ligand = '{}/{}.pdbqt'.format(dst_path, lig_name)  # 准备后配体
        # log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 准备小分子
        self.vina_lig_prep(ligand, pre_ligand)  # 加H，转换为pdbqtg格式
        # 计算
        # 指令
        cmd = 'module purge&&module load vina&&vina --receptor {} --ligand {} --center_x {} --center_y {} --center_z {} ' \
              '--size_x 18.75 --size_y 18.75 --size_z 18.75 --out out.pdbqt --log {} --score_only --cpu 1'.format(pre_protein_file, pre_ligand, x, y, z, log_file)
        # 执行
        os.system(cmd)

    def cal_xscore(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/xscore'.format(path_job)  # 目标路径
        # 蛋白文件
        protein_file = self.protein
        # 定义小分子
        lig_name = ligand.split('.')[0]
        ligand = '{}/{}'.format(path_files, ligand)
        # 定义log文件
        log_file = '{}/{}.txt'.format(dst_path, lig_name)
        # 定义无用文件
        table_file = '{}/{}.table'.format(dst_path, lig_name)
        mdb_file = '{}/{}.mdb '.format(dst_path, lig_name)
        # 创建配置文件
        config_file = '{}/{}.input'.format(dst_path, lig_name)  # 配置文件
        with open(config_file, 'w') as f:  # 创建
            f.write(self.xscore_score.format(protein_file, ligand, table_file, log_file, mdb_file))
        # 计算
        cmd = 'cd /home/xujun/Soft/Score_Function/xscorelinux/bin&&./xscore {}'.format(config_file)
        os.system(cmd)
        # 删除无用文件
        os.remove(table_file)
        os.remove(config_file)

    # moe 软件  ########################################################################################################

    def cal_affinitydG(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/affinitydG'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 小分子预处理
        ligand = '{}/{}'.format(path_files, ligand)
        # 小分子预处理
        self.moe_lig_pre(dst_lig_path, ligand)
        # 蛋白预处理
        self.moe_pro_pre(path_files, dst_lig_path, ligand)
        # 定义配置文件
        dock_batch = '{}/dock_batch.svl'.format(dst_path)  # 对接文件
        score_batch = '{}/score_batch.svl'.format(dst_path)  # 打分文件
        # 创建配置文件
        if not os.path.exists(dock_batch):  # 若不存在则创建对接文件
            with open(dock_batch, 'w') as f:
                f.write(self.affinitydG_dock)
        if not os.path.exists(score_batch):  # 若不存在则创建打分文件
            with open(score_batch, 'w') as f:
                f.write(self.affinitydG_score)
        # 计算
        self.moe_excute(dst_lig_path, dock_batch, score_batch)

    def cal_alphaHB(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/alphaHB'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 小分子预处理
        ligand = '{}/{}'.format(path_files, ligand)
        # 小分子预处理
        self.moe_lig_pre(dst_lig_path, ligand)
        # 蛋白预处理
        self.moe_pro_pre(path_files, dst_lig_path, ligand)
        # 定义配置文件
        dock_batch = '{}/dock_batch.svl'.format(dst_path)  # 对接文件
        score_batch = '{}/score_batch.svl'.format(dst_path)  # 打分文件
        # 创建配置文件
        if not os.path.exists(dock_batch):  # 若不存在则创建对接文件
            with open(dock_batch, 'w') as f:
                f.write(self.alphaHB_dock)
        if not os.path.exists(score_batch):  # 若不存在则创建打分文件
            with open(score_batch, 'w') as f:
                f.write(self.alphaHB_score)
        # 计算
        self.moe_excute(dst_lig_path, dock_batch, score_batch)

    def cal_ase(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/ase'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 小分子预处理
        ligand = '{}/{}'.format(path_files, ligand)
        # 小分子预处理
        self.moe_lig_pre(dst_lig_path, ligand)
        # 蛋白预处理
        self.moe_pro_pre(path_files, dst_lig_path, ligand)
        # 定义配置文件
        dock_batch = '{}/dock_batch.svl'.format(dst_path)  # 对接文件
        # 创建配置文件
        if not os.path.exists(dock_batch):  # 若不存在则创建对接文件
            with open(dock_batch, 'w') as f:
                f.write(self.ase_dock)
        # 计算
        cmd = 'cd {}&&module load moe&&moebatch -run {} -lig ./l.mdb'.format(dst_lig_path, dock_batch)
        # 执行指令
        os.system(cmd)
        # 从mdb获取结果
        self.moe_mdb2result(dst_path, dst_lig_path)

    def cal_GBVIWSAdG(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/GBVIWSAdG'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 小分子预处理
        ligand = '{}/{}'.format(path_files, ligand)
        # 小分子预处理
        self.moe_lig_pre(dst_lig_path, ligand)
        # 蛋白预处理
        self.moe_pro_pre(path_files, dst_lig_path, ligand)
        # 定义配置文件
        dock_batch = '{}/dock_batch.svl'.format(dst_path)  # 对接文件
        score_batch = '{}/score_batch.svl'.format(dst_path)  # 打分文件
        # 创建配置文件
        if not os.path.exists(dock_batch):  # 若不存在则创建对接文件
            with open(dock_batch, 'w') as f:
                f.write(self.GBVIWSAdG_dock)
        if not os.path.exists(score_batch):  # 若不存在则创建打分文件
            with open(score_batch, 'w') as f:
                f.write(self.GBVIWSAdG_score)
        # 计算
        self.moe_excute(dst_lig_path, dock_batch, score_batch)

    def cal_LondondG(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/LondondG'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 小分子预处理
        ligand = '{}/{}'.format(path_files, ligand)
        # 小分子预处理
        self.moe_lig_pre(dst_lig_path, ligand)
        # 蛋白预处理
        self.moe_pro_pre(path_files, dst_lig_path, ligand)
        # 定义配置文件
        dock_batch = '{}/dock_batch.svl'.format(dst_path)  # 对接文件
        score_batch = '{}/score_batch.svl'.format(dst_path)  # 打分文件
        # 创建配置文件
        if not os.path.exists(dock_batch):  # 若不存在则创建对接文件
            with open(dock_batch, 'w') as f:
                f.write(self.LondondG_dock)
        if not os.path.exists(score_batch):  # 若不存在则创建打分文件
            with open(score_batch, 'w') as f:
                f.write(self.LondondG_score)
        # 计算
        self.moe_excute(dst_lig_path, dock_batch, score_batch)

    # gold #############################################################################################################

    def cal_asp(self, ligand, x, y, z):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/asp'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 定义文件
        protein_file = self.protein  # 蛋白
        pre_protein_file = '{}/protein_file.pdb'.format(dst_path)  # 准备后的蛋白
        ligand_file = '{}/{}'.format(path_files, ligand)
        result_sdf = '{}/tem.sdf'.format(dst_lig_path)  # 生成的结果sdf文件
        result_csv = '{}/log.csv'.format(dst_lig_path)  # 由sdf转成的结果数据
        # 蛋白预处理
        if not os.path.exists(pre_protein_file):
            self.gold_pro_pre(protein_file, pre_protein_file)
        # 创建小分子文件夹
        os.mkdir(dst_lig_path)
        # 定义配置文件
        config_file = '{}/score.conf'.format(dst_lig_path)  # 打分文件
        con = self.gold_config.format(x, y, z, ligand_file, result_sdf, 'asp', pre_protein_file)  # 文件内容
        # 创建配置文件
        with open(config_file, 'w')as f:
            f.write(con)
        # 计算
        self.gold_excute(dst_lig_path, config_file)
        # 将sdf转为数据csv
        cmd = 'module load openeye&&convert.py {} {}'.format(result_sdf, result_csv)
        os.system(cmd)

    def cal_chemplp(self, ligand, x, y, z):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/chemplp'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 定义文件
        protein_file = self.protein  # 蛋白
        pre_protein_file = '{}/protein_file.pdb'.format(dst_path)  # 准备后的蛋白
        ligand_file = '{}/{}'.format(path_files, ligand)
        result_sdf = '{}/tem.sdf'.format(dst_lig_path)  # 生成的结果sdf文件
        result_csv = '{}/log.csv'.format(dst_lig_path)  # 由sdf转成的结果数据
        # 蛋白预处理
        if not os.path.exists(pre_protein_file):
            self.gold_pro_pre(protein_file, pre_protein_file)
        # 创建小分子文件夹
        os.mkdir(dst_lig_path)
        # 定义配置文件
        config_file = '{}/score.conf'.format(dst_lig_path)  # 打分文件
        con = self.gold_config.format(x, y, z, ligand_file, result_sdf, 'plp', pre_protein_file)  # 文件内容
        # 创建配置文件
        with open(config_file, 'w')as f:
            f.write(con)
        # 计算
        self.gold_excute(dst_lig_path, config_file)
        # 将sdf转为数据csv
        cmd = 'module load openeye&&convert.py {} {}'.format(result_sdf, result_csv)
        os.system(cmd)

    def cal_chemscore(self, ligand, x, y, z):
        # 还原路径
        path_job = self.path_job # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/chemscore'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 定义文件
        protein_file = self.protein  # 蛋白
        pre_protein_file = '{}/protein_file.pdb'.format(dst_path)  # 准备后的蛋白
        ligand_file = '{}/{}'.format(path_files, ligand)
        result_sdf = '{}/tem.sdf'.format(dst_lig_path)  # 生成的结果sdf文件
        result_csv = '{}/log.csv'.format(dst_lig_path)  # 由sdf转成的结果数据
        # 蛋白预处理
        if not os.path.exists(pre_protein_file):
            self.gold_pro_pre(protein_file, pre_protein_file)
        # 创建小分子文件夹
        os.mkdir(dst_lig_path)
        # 定义配置文件
        config_file = '{}/score.conf'.format(dst_lig_path)  # 打分文件
        con = self.gold_config.format(x, y, z, ligand_file, result_sdf, 'plp', pre_protein_file)  # 文件内容
        # 创建配置文件
        with open(config_file, 'w')as f:
            f.write(con)
        # 计算
        self.gold_excute(dst_lig_path, config_file)
        # 将sdf转为数据csv
        cmd = 'module load openeye&&convert.py {} {}'.format(result_sdf, result_csv)
        os.system(cmd)

    def cal_goldscore(self, ligand, x, y, z):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/goldscore'.format(path_job)  # 目标路径
        dst_lig_path = '{}/{}'.format(dst_path, ligand.split('.')[0])  # 单独小分子文件
        # 定义文件
        protein_file = self.protein  # 蛋白
        pre_protein_file = '{}/protein_file.pdb'.format(dst_path)  # 准备后的蛋白
        ligand_file = '{}/{}'.format(path_files, ligand)
        result_sdf = '{}/tem.sdf'.format(dst_lig_path)  # 生成的结果sdf文件
        result_csv = '{}/log.csv'.format(dst_lig_path)  # 由sdf转成的结果数据
        # 蛋白预处理
        if not os.path.exists(pre_protein_file):
            self.gold_pro_pre(protein_file, pre_protein_file)
        # 创建小分子文件夹
        os.mkdir(dst_lig_path)
        # 定义配置文件
        config_file = '{}/score.conf'.format(dst_lig_path)  # 打分文件
        con = self.gold_config.format(x, y, z, ligand_file, result_sdf, 'plp', pre_protein_file)  # 文件内容
        # 创建配置文件
        with open(config_file, 'w')as f:
            f.write(con)
        # 计算
        self.gold_excute(dst_lig_path, config_file)
        # 将sdf转为数据csv
        cmd = 'module load openeye&&convert.py {} {}'.format(result_sdf, result_csv)
        os.system(cmd)

    # dock #############################################################################################################

    def cal_contact(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig # 分子存放路径
        dst_path = '{}/contact'.format(path_job)  # 目标路径
        # 小分子预处理
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        # 小分子文件
        lig_temp = '{}/{}.sdf'.format(dst_path, lig_name)
        dst_ligand = '{}/{}.mol2'.format(dst_path, lig_name)
        # 小分子转换格式
        self.openeye_convert(ligand_file, lig_temp)  # mol2 to sdf
        # 小分子加电荷
        self.dock_lig_addcharge(lig_temp, dst_path, lig_name)
        # 计算
        config_file = '{}/{}.in'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.dock_contact.format(dst_ligand, lig_name))
        # 执行计算
        self.dock_excute(dst_path, lig_name)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove(lig_temp)  # 删除SDF文件

    def cal_continuous(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/continuous'.format(path_job)  # 目标路径
        # 小分子预处理
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        # 小分子文件
        lig_temp = '{}/{}.sdf'.format(dst_path, lig_name)
        dst_ligand = '{}/{}.mol2'.format(dst_path, lig_name)
        # 小分子转换格式
        self.openeye_convert(ligand_file, lig_temp)  # mol2 to sdf
        # 小分子加电荷
        self.dock_lig_addcharge(lig_temp, dst_path, lig_name)
        # 计算
        config_file = '{}/{}.in'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.dock_continuous.format(dst_ligand, lig_name))
        # 执行计算
        self.dock_excute(dst_path, lig_name)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove(lig_temp)  # 删除SDF文件

    def cal_grid(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/grid'.format(path_job)  # 目标路径
        # 小分子预处理
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        # 小分子文件
        lig_temp = '{}/{}.sdf'.format(dst_path, lig_name)
        dst_ligand = '{}/{}.mol2'.format(dst_path, lig_name)
        # 小分子转换格式
        self.openeye_convert(ligand_file, lig_temp)  # mol2 to sdf
        # 小分子加电荷
        self.dock_lig_addcharge(lig_temp, dst_path, lig_name)
        # 计算
        config_file = '{}/{}.in'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.dock_grid.format(dst_ligand, lig_name))
        # 执行计算
        self.dock_excute(dst_path, lig_name)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove(lig_temp)  # 删除SDF文件

    def cal_hawkins(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/hawkins'.format(path_job)  # 目标路径
        # 小分子预处理
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        # 小分子文件
        lig_temp = '{}/{}.sdf'.format(dst_path, lig_name)
        dst_ligand = '{}/{}.mol2'.format(dst_path, lig_name)
        # 小分子转换格式
        self.openeye_convert(ligand_file, lig_temp)  # mol2 to sdf
        # 小分子加电荷
        self.dock_lig_addcharge(lig_temp, dst_path, lig_name)
        # 计算
        config_file = '{}/{}.in'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.dock_hawkins.format(dst_ligand, lig_name))
        # 执行计算
        self.dock_excute(dst_path, lig_name)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove(lig_temp)  # 删除SDF文件

    def cal_sasa(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/sasa'.format(path_job)  # 目标路径
        # 小分子预处理
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        # 小分子文件
        lig_temp = '{}/{}.sdf'.format(dst_path, lig_name)
        dst_ligand = '{}/{}.mol2'.format(dst_path, lig_name)
        # 小分子转换格式
        self.openeye_convert(ligand_file, lig_temp)  # mol2 to sdf
        # 小分子加电荷
        self.dock_lig_addcharge(lig_temp, dst_path, lig_name)
        # 计算
        config_file = '{}/{}.in'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.dock_sasa.format(dst_ligand, lig_name))
        # 执行计算
        self.dock_excute(dst_path, lig_name)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove(lig_temp)  # 删除SDF文件

    # glide ############################################################################################################

    def cal_sp(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/sp'.format(path_job)  # 目标路径
        # 格点文件
        grid_zip = '{}/glide-grid.zip'.format(dst_path)
        # 小分子预处理
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        # 计算
        config_file = '{}/{}.in'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.glide_score.format(grid_zip, ligand_file, 'SP'))
        # 执行计算
        self.glide_excute(dst_path, config_file, lig_name)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove('{}/{}_raw.maegz'.format(dst_path, lig_name))  # 删除结果文件

    def cal_xp(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/xp'.format(path_job)  # 目标路径
        # 格点文件
        grid_zip = '{}/glide-grid.zip'.format(dst_path)
        # 小分子预处理
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        # 计算
        config_file = '{}/{}.in'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.glide_score.format(grid_zip, ligand_file, 'XP') + '\nWRITE_XP_DESC   True')
        # 执行计算
        self.glide_excute(dst_path, config_file, lig_name)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove('{}/{}_raw.maegz'.format(dst_path, lig_name))  # 删除结果文件

    # plants ###########################################################################################################

    def cal_plp(self, ligand, x, y, z):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/plp'.format(path_job)  # 目标路径
        # 小分子与蛋白
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        pre_protein_file = '{}/protein_file.mol2'.format(dst_path)
        dst_lig_path = '{}/{}'.format(dst_path, lig_name)  # 小分子文件夹
        # 创建文件夹  # 无需创建分子文件夹，软件计算时会创建
        # 计算
        config_file = '{}/{}.config'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.plants_score.format('plp', pre_protein_file, ligand_file, dst_lig_path, x, y, z))
        # 执行计算
        self.plants_excute(config_file)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件

    def cal_plp95(self, ligand, x, y, z):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/plp95'.format(path_job)  # 目标路径
        # 小分子与蛋白
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)
        pre_protein_file = '{}/protein_file.mol2'.format(dst_path)
        dst_lig_path = '{}/{}'.format(dst_path, lig_name)  # 小分子文件夹
        # 创建文件夹  # 无需创建分子文件夹，软件计算时会创建
        # 计算
        config_file = '{}/{}.config'.format(dst_path, lig_name)  # 定义文件路径
        # 创建文件
        with open(config_file, 'w') as f:
            f.write(self.plants_score.format('plp95', pre_protein_file, ligand_file, dst_lig_path, x, y, z))
        # 执行计算
        self.plants_excute(config_file)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件

    # rdock ############################################################################################################

    def cal_rdock(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/rdock'.format(path_job)  # 目标路径
        # 小分子与蛋白
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)  # 源小分子
        pre_protein_file = '{}/protein_file.mol2'.format(dst_path)  # 准备后的蛋白
        pre_ligand = '{}/{}.sdf'.format(dst_path, lig_name)  # 准备后的小分子
        dockprm = '/opt/rDock/2013.1/data/scripts/minimise.prm'  # 参数文件
        config_file = '{}/{}.prm'.format(dst_path, lig_name)  # 配置文件
        # 转换分子格式
        self.openeye_convert(ligand_file, pre_ligand)
        # 创建配置文件
        with open(config_file, 'w') as f:
            f.write(self.rdock_score.format(pre_protein_file, pre_ligand))
        # 执行计算
        self.rdock_excute(lig_name, dst_path, config_file, pre_ligand, dockprm)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove(pre_ligand)  # 删除小分子 4KB
        os.remove('{}/{}.as'.format(dst_path, lig_name))  # 删除中间文件 4MB
        os.remove('{}/{}_cav1.grd'.format(dst_path, lig_name))  # 删除中间文件 4MB

    def cal_rdock_sol(self, ligand):
        # 还原路径
        path_job = self.path_job  # 任务文件夹路径
        path_files = self.path_for_lig  # 分子存放路径
        dst_path = '{}/rdock_sol'.format(path_job)  # 目标路径
        # 小分子与蛋白
        lig_name = ligand.split('.')[0]
        ligand_file = '{}/{}'.format(path_files, ligand)  # 源小分子
        pre_protein_file = '{}/protein_file.mol2'.format(dst_path)  # 准备后的蛋白
        pre_ligand = '{}/{}.sdf'.format(dst_path, lig_name)  # 准备后的小分子
        dockprm = '/opt/rDock/2013.1/data/scripts/minimise_solv.prm'  # 参数文件
        config_file = '{}/{}.prm'.format(dst_path, lig_name)  # 配置文件
        # 转换分子格式
        self.openeye_convert(ligand_file, pre_ligand)
        # 创建配置文件
        with open(config_file, 'w') as f:
            f.write(self.rdock_score.format(pre_protein_file, pre_ligand))
        # 执行计算
        self.rdock_excute(lig_name, dst_path, config_file, pre_ligand, dockprm)
        # 删除多余文件
        os.remove(config_file)  # 删除配置文件
        os.remove(pre_ligand)  # 删除小分子 4KB
        os.remove('{}/{}.as'.format(dst_path, lig_name))  # 删除中间文件 4MB
        os.remove('{}/{}_cav1.grd'.format(dst_path, lig_name))  # 删除中间文件 4MB

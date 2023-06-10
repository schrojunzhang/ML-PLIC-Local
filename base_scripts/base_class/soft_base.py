#!usr/bin/env python
# -*- coding:utf-8 -*-
# author = zhang xujun
# time = 2020-07-09

import os
import glob
import shutil
from base_scripts.base_class.global_base import job_control


class soft(job_control):

    def __init__(self, job_id):
        super().__init__(job_id)
        self.affinitydG_dock = '''#svl
function DockAtoms, DockFile;
function DockMDBwAtoms, DockMDBwFile;
local function RunDock[outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file]
    if isnull outf then
    outf = 'dock.mdb';
    endif;
    local psys = SystemPush [];

    Open './complex_prep.moe';

    local c = Chains[];
    local c0 = cAtoms c;
    local l0 = length c0;
    local site = c0(l0);
    local i = 0;
    local rec = [];
    local lig;
    while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop

    local opt = [
    outrmsd: 1,
    sel_ent_only_rec: 0,
    sel_ent_only: 0,
    wall: [ '', 0, [ 0, 0, 0 ], [ 1000000, 1000000, 1000000 ], 0 ],
    csearch: 1,
    confGenMethod: 'Rotate Bonds',
    placement: 'None',
    placement_opt: [  ],
    scoring: 'None',
    scoring_opt: [  ],
    dup_placement: 1,
    maxpose: 1,
    refine: 'Rigid Receptor',
    refine_opt: [ fixrec : 'Fix',	ed_map : 'Fo',	ed_f : 'Simulated',	ed_phi : 'Simulated',	ed_f2 : 'Simulated',	ed_path : '',	ed_res : 2.5,	ed_sfdata : [ [  ], [  ], [  ], 'Simulated', 'Simulated', 'Fo' ],	ed_surflevelD : 3,	cutoff : 6,	wholeres : 1,	mmgbvi : 1,	packsidechains : 1,	rigidlig : 0,	tether : 10,	gtest : 0.01,	maxit : 500,	OverrideSetup : 1,	k_potl : 100,	roffset : 0.4 ],
    rescoring: 'Affinity dG',
    rescoring_opt: [ hbond : -0.66,	ionic : 1,	mlig : -1,	cHH : -0.01235,	cHP : 0.02497,	cXX : -0.00834 ],
    dup_refine: 1,
    remaxpose: 1,
    descexpr: '',
    receptor_mfield: '',
    ligand_mfield: 'mol',
    rxnFile: '',
    rxsite: [  ],
    edsupport: 1,
    ed_data: [ ed_dockpath : '' ],
    check_pose_geom: [  ],
    multiLigand: 0,
    need_dmat: 1,
    gen_plif: 1,
    tempDB: '/tmp/1iep_ligand_prep98896.mdb',
    ph4: ph4file,
    ligmdbname: ligmdbname,
    recmdbname: recmdbname,
    BatchFile: 'dock_batch.svl'
    ];
    pot_Load '$MOE/lib/Amber10EHT.ff.gz';

    pot_Setup [
    strEnable: 1,
    angEnable: 1,
    stbEnable: 1,
    oopEnable: 1,
    torEnable: 1,
    vdwEnable: 1,
    eleEnable: 1,
    solEnable: 0,
    resEnable: 1,
    strWeight: 1,
    angWeight: 1,
    stbWeight: 1,
    oopWeight: 1,
    torWeight: 1,
    vdwWeight: 1,
    eleWeight: 1,
    solWeight: 1,
    resWeight: 1,
    cutoffEnable: 1,
    cutoffOn: 8,
    cutoffOff: 10,
    eleDist: 2,
    vdwScale14: 0.5,
    vdwBuffer1: 0,
    vdwBuffer2: 0,
    eleScale14: 0.833333,
    eleDielectric: 1,
    eleBuffer: 0,
    solDielectric: 80,
    solDielectricOffset: 0,
    state0: 1,
    state1: 0,
    state2: 1,
    threadCount: 0
    ];
    if isnull ligmdb and isnull recmdb then
    DockAtoms [rec, site, lig, outf, opt];
    elseif length ligmdb and isnull recmdb then
    DockFile [rec, site, ligmdb, outf, opt];
    elseif isnull ligmdb and length recmdb then
    DockMDBwAtoms [recmdb, site, lig, outf, opt];
    elseif length ligmdb and length recmdb then
    DockMDBwFile [recmdb, site, ligmdb, outf, opt];
    endif
    write ['Docking finished at {}.', asctime []];
    SystemPop psys;
endfunction;

global argv;
function ArgvPull;

local function main []
    ArgvReset ArgvExpand argv;
    local [recmdb, ligmdb, ph4file, outf] = ArgvPull [
    ['-rec', '-lig','-ph4','-o'],
    1
    ];

    local ligmdbname, recmdbname;

    local delfile1 = [];
    local delfile2 = [];
    local delfile3 = [];

    local handle;
    if isnull outf then outf = 'dock.mdb'; endif
    if not isnull ligmdb then
    handle = _db_Open [ligmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read ligand mdb file {}', ligmdb];
    endif
    db_Close handle;
    ligmdbname = ligmdb;
    endif

    if not isnull recmdb then
    handle = _db_Open [recmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read receptor mdb file {}', recmdb];
    endif
    db_Close handle;
    recmdbname = recmdb;
    endif

    if not isnull ph4file then
    if not (ftype ph4file == 'file') then
        exit twrite ['File {} not found', ph4file];
    endif
    endif

    local oldttl = task_title -1;
    if second task_wfork [master: 'none'] == 'child' then
    task_settitle [-1, twrite ['Dock: [{}] {}',
        second app token wordsplit[ string oldttl, "' "], outf]];
    RunDock [outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file];
    exit [];
    endif

    fdelete delfile1;
    fdelete delfile2;
    fdelete delfile2;
endfunction
'''
        self.affinitydG_score = '''#svl
// Implementation of the Affinity dG scoring function, which is relatively
// straightforward. This is a good boilerplate example for developing other
// scoring functions.

function dock_ReceptorOpen;
function dock_LigandOpen;
function dock_LigandClose;
function dock_ReceptorClose;

local function dock_score_Affinity_dG [cmd, arg, opt]
    const DEFAULTS = [
    hbond:	-0.66,		// hbond coefficient
    ionic:	1.00,		// ionic contact
    mlig:	-1.00,		// metal ligation
    cHH:	-0.01235,	//  1 [H][H]
    cHP:	0.02497,	//  2 [H][P]
    cXX:	-0.00834	//  3 [*][*]
    ];

    const PANEL = [
    Mbox : [
        columns: 2, columnMajor:1,
        Text : [
        title: 'Hydrogen Bond:', name: 'hbond', type: 'real',
        shortcut: ['-0.5','-0.6','-0.7','-0.8','-0.9','-1.0'],
        bubbleHelp: 'The ideal hydrogen bond affinity.'
        ],
        Text : [
        title: 'Ionic Contact:', name: 'ionic', type: 'real',
        shortcut: ['0.6','0.7','0.8','0.9','1.0','1.1','1.2','1.3'],
        bubbleHelp: 'The ionic Coulomb coefficient.'
        ],
        Text : [
        title: 'Metal Ligation:', name: 'mlig', type: 'real',
        shortcut: ['-0.7','-0.8','-0.9','-1.0','-1.1','-1.2','-1.3'],
        bubbleHelp: 'The metal ligation affinity.'
        ],
        Text : [
        title: 'Hydrophobic Contact:', name: 'cHH', type: 'real',
        bubbleHelp: 'Hydrophobic-Hydrophobic contact energy.'
        ],
        Text : [
        title: 'Hydrophobic-Polar:', name: 'cHP', type: 'real',
        bubbleHelp: 'Hydrophobic-Polar contact energy.'
        ],
        Text : [
        title: 'Atom-Atom:', name: 'cXX', type: 'real',
        bubbleHelp: 'General atom contact energy.'
        ]
    ]
    ];

    const CS_HBCUTOFF = 3.5;		// proximity cutoffs
    const CS_MCUTOFF = 3.0;
    const CS_CUTOFF = 7.5;
    const CS_ICUTOFF = 7.5;

    function atype obj;
    local Qindex = *obj.Qindex;
    local at = rep [0, length Qindex];			// X
    at | *obj.m_don[Qindex] or *obj.m_acc[Qindex] = 1;	// P
    at | *obj.m_gre[Qindex] or *obj.m_hyd[Qindex] = 2;	// H
    return at;
    endfunction

    if cmd === 'ID' then
    return 'Affinity dG';

    elseif cmd === 'configpanelwidgets' then
    return PANEL;

    elseif cmd === 'configpanelevent' then	// arg = [val,trig]
    return 0;

    elseif cmd === 'configvalues' then		// arg = val
    return tagcat [arg, DEFAULTS];

    elseif cmd === 'openReceptor' then		// arg = rec
    dock_score_Affinity_dG ['closeReceptor', arg, opt];

    local rQindex = *arg.Qindex;
    local rpos = apt get [mol_aPos *arg.mol, [rQindex]];
    local rrad = *arg.radius[rQindex];

    local rhb = bitor [
        select [1, 0, *arg.Qmask and *arg.m_don],
        select [2, 0, *arg.Qmask and *arg.m_acc]
    ];
    local rhbpos = apt mget [mol_aPos *arg.mol, [rhb]];
    local rmpos = apt mget [mol_aPos *arg.mol, [*arg.m_tmetal]];
    local ipos = apt mget [mol_aPos *arg.mol, [*arg.fc <> 0]];

    *arg.Affinity_dG = [
        rrad:	rrad,
        prox:	prox_open [CS_CUTOFF, rpos, CS_CUTOFF],
        hbprox:	prox_open [CS_HBCUTOFF, rhbpos, CS_HBCUTOFF],
        rhb:	pack rhb,
        mprox:	prox_open [CS_MCUTOFF, rmpos, CS_MCUTOFF],
        iprox:	prox_open [CS_ICUTOFF, ipos, CS_ICUTOFF],
        rfc:	*arg.fc | *arg.fc <> 0,
        rtype:	atype arg
    ];
    return;

    elseif cmd === 'closeReceptor' then		// arg = rec
    prox_close *arg.Affinity_dG.prox;
    prox_close *arg.Affinity_dG.hbprox;
    prox_close *arg.Affinity_dG.mprox;
    prox_close *arg.Affinity_dG.iprox;
    *arg.Affinity_dG = [];
    return;

    elseif cmd === 'openLigand' then		// arg = lig
    dock_score_Affinity_dG ['closeLigand', arg, opt];

        // If specified, we'll restrict the various constants to the
        // used atoms (provided by the caller)

    local uQmask = *arg.Qmask, uQindex = *arg.Qindex;
    local ufc    = *arg.fc,    um_mlig = *arg.m_mlig;
    if length *arg.ligData.useAtoms then
        local m = not *arg.ligData.useAtoms;
        uQmask  | m = 0;
        ufc     | m = 0;
        um_mlig | m = 0;

        uQindex = uQindex | indexof [uQindex, x_pack *arg.ligData.useAtoms];
    endif

    local lhb = bitor [
        select [2, 0, uQmask and *arg.m_don],
        select [1, 0, uQmask and *arg.m_acc]
    ];
    *arg.Affinity_dG = [
        lQindex:    uQindex,
        lrad:	*arg.radius[uQindex],		// radius
        ltype:	atype arg,			// X/H/P
        xilig:	x_pack (ufc <> 0),		// ion indices
        lfc:	pack ufc,
        xmlig:	x_pack um_mlig,			// metal ligation
        xhlig:	x_pack lhb,			// hb indices
        lhb:	pack lhb
    ];
    return;

    elseif cmd === 'closeLigand' then		// arg = lig
    *arg.Affinity_dG = [];
    return;

    elseif not (cmd === 'score') then		// arg = [ligpos,lig,rec]
    return;
    endif

    // perform the scoring of the ligpos using the lig,rec,opt

    local [ligpos, lig, rec] = arg;
    local rctx = *rec.Affinity_dG;
    local lctx = *lig.Affinity_dG;

    // calculate the hydrogen bond term

    local lhbpos = apt get [ligpos, [lctx.xhlig]];
    local [seg,idx,r2] = prox_find [rctx.hbprox, lhbpos, 0];
    local hbond = bitor [stretch [lctx.lhb, seg], rctx.rhb[idx]];
    local hbpot = maxE[0, minE[1,
    1 - abs (sqrt r2 -2.85) * inv(CS_HBCUTOFF-2.85)
    ]];
    local HB = add ((hbond == 3) * hbpot);

    // calculate the metal ligation term

    local lmligpos = apt get [ligpos, [lctx.xmlig]];
    [seg,idx,r2] = prox_find [rctx.mprox, lmligpos, 0];
    local mpot = maxE[0, minE[1,
    1 - abs (sqrt r2 - 2.2) * inv(CS_MCUTOFF - 2.2)
    ]];
    local MLIG = add mpot;

    // calculate the ionic part

    local iligpos = apt get [ligpos, [lctx.xilig]];
    [seg,idx,r2] = prox_find [rctx.iprox, iligpos, 0];
    r2 = maxE [1, sqrt r2];
    local ipot = (
      maxE[0, minE[1, (CS_ICUTOFF - r2) * inv CS_ICUTOFF]]
    * (rctx.rfc[idx] * stretch [lctx.lfc, seg])
    );
    local ION = add ipot;

    // calculate the coefficients

    const CMAP = [ 3, 3, 3, 3, 3, 2, 3, 2, 1 ];
    local coeff = [ opt.cHH + opt.cXX, opt.cHP + opt.cXX, opt.cXX ][CMAP];

    local lQindex = lctx.lQindex;		// ligand data
    local lpos = apt get [ligpos, [lQindex]];
    local lrad = lctx.lrad;
    local ltype = lctx.ltype;

    [seg,idx,r2] = prox_find [rctx.prox, lpos, 0];
    local r = sqrt r2;

    local lidx = stretch [x_id lrad, seg];
    local rsum = lrad[lidx] + rctx.rrad[idx];

    local cindex = inc iadd [ (imul [ltype,3])[lidx], rctx.rtype[idx] ];

    local t = maxE [0, minE [1, (CS_CUTOFF - r + rsum) * inv CS_CUTOFF]];

    t = t * coeff[cindex];


    // add in the repulsion term

    local REPSCALE = 6.0;
    local i = x_pack (r < rsum);

    local v = maxE[ 0, minE[ 1, 1 - r[i]/rsum[i] ]];
    t[i] = t[i] + (REPSCALE * cube v * (6 * sqr v - 15 * v + 10));

    local affinitydG = opt.hbond * HB + opt.ionic * ION + opt.mlig * MLIG + add t;

    return [affinitydG, HB, ION, MLIG, add t];
endfunction

local function main []
    local comen = Open './complex_prep.moe';
    local c = Chains[];
    local c0 = cAtoms c;
    local l0 = length c0;
    local lig = c0(l0);
    local i = 0;
    local rec = [];
    while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop
    //Close comen;

    local ligen = db_ReadColumn [ 'dock.mdb', 'mol' ];
    local cl = mol_Create ligen(1);	
    local cl0 = cAtoms cl;
    lig = cl0(1);		

    local rec_dock = dock_ReceptorOpen [rec, [aPos lig],[],[]];
    local lig_dock = dock_LigandOpen [lig, []];
    dock_score_Affinity_dG ['openReceptor', rec_dock, []];
    dock_score_Affinity_dG ['openLigand', lig_dock, []];	
    local affinitydG = dock_score_Affinity_dG ['score', [aPos lig,lig_dock,rec_dock], [ hbond : -0.66,	ionic : 1,	mlig : -1,	cHH : -0.01235,	cHP : 0.02497,	cXX : -0.00834 ]];
    dock_score_Affinity_dG ['closeReceptor', rec_dock, []];
    dock_score_Affinity_dG ['closeLigand', lig_dock, []];	
    dock_ReceptorClose rec_dock;
    dock_LigandClose lig_dock;
    //print affinitydG;
    fwrite [ 'log.txt', '{},{},{},{},{}', affinitydG[1], affinitydG[2], affinitydG[3], affinitydG[4], affinitydG[5]]; 
endfunction
'''
        self.alphaHB_dock = '''#svl
function DockAtoms, DockFile;
function DockMDBwAtoms, DockMDBwFile;
local function RunDock[outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file]
    if isnull outf then
    outf = 'dock.mdb';
    endif;
    local psys = SystemPush [];


    Open './complex_prep.moe';

    local c = Chains[];
    local c0 = cAtoms c;
    local l0 = length c0;
    local site = c0(l0);
    local i = 0;
    local rec = [];
    local lig;
    while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop

    local opt = [
    outrmsd: 1,
    sel_ent_only_rec: 0,
    sel_ent_only: 0,
    wall: [ '', 0, [ 0, 0, 0 ], [ 1000000, 1000000, 1000000 ], 0 ],
    csearch: 1,
    confGenMethod: 'Rotate Bonds',
    placement: 'None',
    placement_opt: [  ],
    scoring: 'None',
    scoring_opt: [  ],
    dup_placement: 1,
    maxpose: 1,
    refine: 'Rigid Receptor',
    refine_opt: [ fixrec : 'Fix',	ed_map : 'Fo',	ed_f : 'Simulated',	ed_phi : 'Simulated',	ed_f2 : 'Simulated',	ed_path : '',	ed_res : 2.5,	ed_sfdata : [ [  ], [  ], [  ], 'Simulated', 'Simulated', 'Fo' ],	ed_surflevelD : 3,	cutoff : 6,	wholeres : 1,	mmgbvi : 1,	packsidechains : 1,	rigidlig : 0,	tether : 10,	gtest : 0.01,	maxit : 500,	OverrideSetup : 1,	k_potl : 100,	roffset : 0.4 ],
    rescoring: 'Alpha HB',
    rescoring_opt: [ was : 10,	whbpf : 4 ],
    dup_refine: 1,
    remaxpose: 5,
    descexpr: '',
    receptor_mfield: '',
    ligand_mfield: 'mol',
    rxnFile: '',
    rxsite: [  ],
    edsupport: 1,
    ed_data: [ ed_dockpath : '' ],
    check_pose_geom: [  ],
    multiLigand: 0,
    need_dmat: 1,
    gen_plif: 1,
    tempDB: '/tmp/1iep_ligand_prep98896.mdb',
    ph4: ph4file,
    ligmdbname: ligmdbname,
    recmdbname: recmdbname,
    BatchFile: 'dock_batch.svl'
    ];
    pot_Load '$MOE/lib/Amber10EHT.ff.gz';

    pot_Setup [
    strEnable: 1,
    angEnable: 1,
    stbEnable: 1,
    oopEnable: 1,
    torEnable: 1,
    vdwEnable: 1,
    eleEnable: 1,
    solEnable: 0,
    resEnable: 1,
    strWeight: 1,
    angWeight: 1,
    stbWeight: 1,
    oopWeight: 1,
    torWeight: 1,
    vdwWeight: 1,
    eleWeight: 1,
    solWeight: 1,
    resWeight: 1,
    cutoffEnable: 1,
    cutoffOn: 8,
    cutoffOff: 10,
    eleDist: 2,
    vdwScale14: 0.5,
    vdwBuffer1: 0,
    vdwBuffer2: 0,
    eleScale14: 0.833333,
    eleDielectric: 1,
    eleBuffer: 0,
    solDielectric: 80,
    solDielectricOffset: 0,
    state0: 1,
    state1: 0,
    state2: 1,
    threadCount: 0
    ];
    if isnull ligmdb and isnull recmdb then
    DockAtoms [rec, site, lig, outf, opt];
    elseif length ligmdb and isnull recmdb then
    DockFile [rec, site, ligmdb, outf, opt];
    elseif isnull ligmdb and length recmdb then
    DockMDBwAtoms [recmdb, site, lig, outf, opt];
    elseif length ligmdb and length recmdb then
    DockMDBwFile [recmdb, site, ligmdb, outf, opt];
    endif
    write ['Docking finished at {}.', asctime []];
    SystemPop psys;
endfunction;

global argv;
function ArgvPull;

local function main []
    ArgvReset ArgvExpand argv;
    local [recmdb, ligmdb, ph4file, outf] = ArgvPull [
    ['-rec', '-lig','-ph4','-o'],
    1
    ];

    local ligmdbname, recmdbname;

    local delfile1 = [];
    local delfile2 = [];
    local delfile3 = [];

    local handle;
    if isnull outf then outf = 'dock.mdb'; endif
    if not isnull ligmdb then
    handle = _db_Open [ligmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read ligand mdb file {}', ligmdb];
    endif
    db_Close handle;
    ligmdbname = ligmdb;
    endif

    if not isnull recmdb then
    handle = _db_Open [recmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read receptor mdb file {}', recmdb];
    endif
    db_Close handle;
    recmdbname = recmdb;
    endif

    if not isnull ph4file then
    if not (ftype ph4file == 'file') then
        exit twrite ['File {} not found', ph4file];
    endif
    endif

    local oldttl = task_title -1;
    if second task_wfork [master: 'none'] == 'child' then
    task_settitle [-1, twrite ['Dock: [{}] {}',
        second app token wordsplit[ string oldttl, "' "], outf]];
    RunDock [outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file];
    exit [];
    endif

    fdelete delfile1;
    fdelete delfile2;
    fdelete delfile2;
endfunction
'''
        self.alphaHB_score = '''#svl
function ph4_QueryMatch_Open, ph4_QueryMatch_Close;
function ph4_QueryMatch_Qdata;

// ---- Temporary ph4 functions ---

local function LoosenPH4 [ph4data, amount]
    if isnull ph4data then return []; endif
    if amount <= 0.0 then return ph4data; endif
    ph4data.F.r = ph4data.F.r + amount;
    local m = stretch [ph4data.V.ebits,ph4data.V.size] == 1;
    local y = ph4data.VS.r;
    y | m = apt max [(y | m) - amount, 0.0];
    y | not m = (y | not m) + amount;
    ph4data.VS.r = y;
    return ph4data;
endfunction

// ---- End of temporary ph4 functions ----

// ===========================================================================
// -------------- HYDROGEN BONDING PROJECTED FEATURES SCORE ------------------
// ===========================================================================

// Transition metal as donors added on 2005.8.25.
// Coefficient for this term is assumed to be the same as other donors.

function pro_StandardRes;

// cent is the position of the projecting atom.
// nbrs is / are the positions of cent's neighbours.
// Returns the projected point in the direction opposite to the vectors
// nbrs - cent.

local function Projecting1 [cent, nbrs, dist]
    local ofs = vnormalize app add vnormalize (cent - nbrs);
    return cent + dist * ofs;
endfunction

// cent is the position of the projecting atom.
// nbrs are the positions of cent's TWO neighbours.
//
// It is the caller's responsibility to ensure nbrs has length 2.
//
// Returns two points projected from cent in the direction opposite to nbrs.
// The points are out of the plane of cent and nbrs.

local function Projecting2t [cent, nbrs, d_inline, d_oop]
    local vs = cent - nbrs;
    local ofs = d_inline * vnormalize app add vnormalize vs;
    local nn = d_oop * vnormalize vcross tr vs;
    return [cent + ofs + nn , cent + ofs - nn];
endfunction

// cent is the position of the projecting atom.
// nbr1 is the position of the only neighbour of cent.
// nbr2 are the positions of the TWO neighbours of nbr1 (excl. cent).
//
// It is the caller's responsibility to ensure nbrs has length 2.
//
// Returns two points projected from cent in the direction opposite to nbrs.
// The points are in the plane of cent, nbr1, nbr2.

local function Projecting2p [cent, nbr1, nbr2, d_inline, d_ppd]
    local ofs = d_inline * vnormalize (cent - nbr1);
    local nn = vcross tr (cent - nbr2);
    nn = d_ppd * vnormalize vcross [nn, ofs];
    return [cent + ofs + nn, cent + ofs - nn];
endfunction

function ph4_aDonor, ph4_aTautomerDonor, ph4_aPlanar;
function ph4_aAcceptor, ph4_aTautomerAcceptor;

// Generates projected features iven a list of projecting atoms pjatoms.

local function GenProjFeatIn [pjatoms, type1, type2]
    const DHB = 2.9;
    const D1 = DHB / 2.0;
    const D2 = D1 * sqrt 3.0;

    local src = [];
    local pftype = [];
    local pfpos = [];
    local noprojlist = [];

    local currat, p1, p2;
    for currat in pjatoms loop
	local nnb = aHeavyValence currat;
	local cent = aPos currat;
	local b1 = cat aBonds currat;
	b1 = b1 | aAtomicNumber b1 > 1;
	if nnb == 2 and ph4_aPlanar currat then		// Type 1

	    p1 = Projecting1 [cent, aPos b1, DHB];
	    src = append [src, currat];
	    pftype = append [pftype, type1];
	    pfpos = append [pfpos, p1];

	elseif nnb == 1 and ph4_aPlanar currat then	// Type 2

	    local b2 = cat aBonds b1;
	    b2 = diff [b2, currat];
	    if length b2 <> 2 then
		write ['*** Unusual topology near \t{t:} {t:} : {t:}',
		    rName aResidue b1, rUID aResidue b1, aName b1
		];
		continue;
	    endif
	    [p1, p2] = Projecting2p [cent, aPos b1, aPos b2, D1, D2];
	    src = append [src, currat];
	    pftype = append [pftype, type2];
	    pfpos = append [pfpos, p1];
	    src = append [src, currat];
	    pftype = append [pftype, type2];
	    pfpos = append [pfpos, p2];

	else

	    noprojlist = append [noprojlist, currat];

	endif    
    endloop
    return [src, pftype, tr pfpos, noprojlist];
endfunction

// Patch to fix imidazoles. In general, histidines are assigned
// a positive charge. Hence the N's will not be detected as possible
// acceptors.

local function hisnmask atoms =
    (rName aResidue atoms == 'HIS') and m_join [aName atoms, ['ND1','NE2']]
;

local function cyssmask atoms =
    (rName aResidue atoms == 'CYS') and (aName atoms == 'SG')
;

local function PO4SO4mask atoms =
    sm_Match ['O~P(~O)~O', atoms] or sm_Match ['O~S(~O)~O', atoms]
;

local function SulfonamideN atoms =
    sm_Match ['[NQ1]~S(~[OQ1])~[OQ1]', atoms]
;

// For docking applications, use HBPF_InitLigand instead.

local function NHAtomTyping ligat
    ligat = ligat | aAtomicNumber ligat > 1;
    return iadd[		// 1 for don, 2 for acc
	(ph4_aDonor ligat or ph4_aTautomerDonor ligat)
	and not PO4SO4mask ligat,
	imul [2, ph4_aAcceptor ligat or ph4_aTautomerAcceptor ligat
	    or SulfonamideN ligat
	]
    ];
endfunction

// No need to "Close Ligand".

local function HBPF_InitLigand lig
    if not isnull *lig.HBPF_nhatype then return; endif

    local ligch = mol_Create *lig.mol;
    *lig.HBPF_nhatype = NHAtomTyping cat cAtoms ligch;
    oDestroy ligch;
endfunction

// Generate projecting features and atom self features from the acceptors
// and donors amongst atoms. Features outside of the box (plow, phigh)
// are discarded. See code for precise description of accs and dons.
// Projected features that has an atom witin 1.5 AA of it are discarded.
//
// Returns [src, type, pos] where
//    src is a list of the source atoms,
//    type is a list of the type of features,
//    pos is the positions of the features.
//
// Type A1, A2 are projected features from acceptors,
//	D1, D2				   donors,
//	B1, B2				   dual character atoms,
//	AN, DN, BN are non-projecting features (atom selfs).

local function GenProjFeatLoop [atoms, plow, phigh]
    local srcall=[], pftypeall=[], pfposall=[[],[],[]];

    local pos = aPos atoms;
    atoms = atoms | andE ((plow-3.0) < pos) and andE (pos < (phigh+3.0));

    local proxk = prox_open [3.0, aPos atoms, 0.0];

    local dons = atoms | (ph4_aDonor atoms or ph4_aTautomerDonor atoms)
	and not PO4SO4mask atoms
    ;
    local accs = atoms | ph4_aAcceptor atoms or ph4_aTautomerAcceptor atoms;
    accs = cat [accs, atoms | sm_Match ['[#G7Q0]', atoms]];
    local both = join [dons, accs];

    both = uniq cat [both, atoms | hisnmask atoms];

    dons = diff [dons, both];
    accs = diff [accs, both];

    local [srcb, pftypeb, pfposb, nplb] = GenProjFeatIn [both, 'B1', 'B2'];
    local m = andE (plow < pfposb) and andE (pfposb < phigh);
    srcb = srcb | m;
    pftypeb = pftypeb | m;
    pfposb = pfposb || [m];

    if anytrue m then
	m = not first prox_find [proxk, pfposb, 1.5];
	srcb = srcb | m;
	pftypeb = pftypeb | m;
	pfposb = pfposb || [m];
    endif

    srcall = cat [srcall, srcb, nplb];
    pftypeall = cat [pftypeall, pftypeb, rep['BN', length nplb]];
    pfposall = apt cat [pfposall, pfposb, aPos nplb];

    local [srcd, pftyped, pfposd, npld] = GenProjFeatIn [dons, 'D1', 'D2'];
    m = andE (plow < pfposd) and andE (pfposd < phigh);
    srcd = srcd | m;
    pftyped = pftyped | m;
    pfposd = pfposd || [m];

    if anytrue m then
	m = not first prox_find [proxk, pfposd, 1.5];
	srcd = srcd | m;
	pftyped = pftyped | m;
	pfposd = pfposd || [m];
    endif

    srcall = cat [srcall, srcd, npld];
    pftypeall = cat [pftypeall, pftyped, rep['DN', length npld]];
    pfposall = apt cat [pfposall, pfposd, aPos npld];

    local [srca, pftypea, pfposa, npla] = GenProjFeatIn [accs, 'A1', 'A2'];
    m = andE (plow < pfposa) and andE (pfposa < phigh);
    srca = srca | m;
    pftypea = pftypea | m;
    pfposa = pfposa || [m];

    if anytrue m then
	m = not first prox_find [proxk, pfposa, 1.5];
	srca = srca | m;
	pftypea = pftypea | m;
	pfposa = pfposa || [m];
    endif

    srcall = cat [srcall, srca, npla];
    pftypeall = cat [pftypeall, pftypea, rep['AN', length npla]];
    pfposall = apt cat [pfposall, pfposa, aPos npla];

    prox_close proxk;

    return [srcall, pftypeall, pfposall];
endfunction

// For docking applications, use HBPF_OpenReceptor instead.

local function InitProjFeat [recat, plow, phigh]
    local [featsrc, feattype, featpos] = GenProjFeatLoop [recat,
	plow, phigh
    ];

    local PFOptr = dvar_open [];
    *PFOptr.feat_accnb_pos = featpos || [
	indexof [feattype, ['A1','A2']] and not (aName featsrc == 'O' and
	indexof [rName aResidue featsrc, pro_StandardRes []])
    ];
    *PFOptr.feat_acc_pos = featpos || [indexof [feattype, ['A1','A2']]];
    *PFOptr.feat_don_pos = featpos || [indexof [feattype, ['D1','D2']]];
    *PFOptr.feat_both_pos = featpos || [indexof [feattype, ['B1','B2']]];
    *PFOptr.feat_aa_pos = featpos || [feattype == 'AN'];
    *PFOptr.feat_da_pos = featpos || [feattype == 'DN'];
    *PFOptr.feat_ba_pos = featpos || [feattype == 'BN'];
#if 1
    *PFOptr.feat_met_pos = aPos (recat | sm_Match ['[#T]', recat]);
#else
    *PFOptr.feat_da_pos = apt cat [*PFOptr.feat_da_pos,
	aPos (recat | sm_Match ['[#T]', recat])
    ];
#endif
    return PFOptr;
endfunction

local function HBPF_OpenReceptor rec
    if isfalse dvar_ref *rec.HBPFObj then
	local snew = SystemOpen [];
	local scurr = SystemCurrent snew;
	local recch = mol_Create *rec.mol;
	local recat = cat cAtoms recch;
	*rec.HBPFObj = InitProjFeat [recat, *rec.sitemin, *rec.sitemax];
	oDestroy recch;
	SystemCurrent scurr;
	SystemClose snew;
    else
	*rec.HBPFObj = dvar_open *rec.HBPFObj;
    endif
endfunction

// Simply use dvar_close. Arrays will automatically be freed.
// For docking applications, use HBPF_CloseReceptor instead.
#if 0
global function CloseProjFeat PFOptr
    *PFOptr.feat_accnb_pos = [];
    *PFOptr.feat_acc_pos = [];
    *PFOptr.feat_don_pos = [];
    *PFOptr.feat_both_pos = [];
    *PFOptr.feat_aa_pos = [];
    *PFOptr.feat_da_pos = [];
    *PFOptr.feat_ba_pos = [];
    *PFOptr.feat_da_pos = [];
    *PFOptr.feat_met_pos = [];
    dvar_close PFOptr;
endfunction
#endif

local function HBPF_CloseReceptor rec = dvar_close *rec.HBPFObj;

local function ProjFeatScore [ligpos, ligtype, PFObj]
#if 0
    if l_length ligpos <> length ligtype then
	exit 'Wrong input for HBPF_Score';
    endif
#endif

    local proxk = prox_open [3.5, ligpos, 0.0];
    local seg, idx, d2, occtype, fitvec;

	// all acceptor projecting

    [seg, idx, d2] = prox_find [proxk, *PFObj.feat_acc_pos, 1.5];
    occtype = split [get [ligtype, idx], pack seg];
    fitvec = app anytrue bitand [1, occtype];
    local nmat_apj = add fitvec;
    local nnm_apj = add not fitvec;

	// donor projecting

    [seg, idx, d2] = prox_find [proxk, *PFObj.feat_don_pos, 1.5];
    occtype = split [get [ligtype, idx], pack seg];
    fitvec = app anytrue bitand [2, occtype];
    local nmat_dpj = add fitvec;
    local nnm_dpj = add not fitvec;

	// dual projecting

    [seg, idx, d2] = prox_find [proxk, *PFObj.feat_both_pos, 1.5];
    occtype = split [get [ligtype, idx], pack seg];
    fitvec = app anytrue occtype;
    local nmat_bpj = add fitvec;
    local nnm_bpj = add not fitvec;

	// acceptor point

    idx = second prox_find [proxk, *PFObj.feat_aa_pos, 3.5];
    fitvec = bitand [1, get [ligtype, idx]];
    local nmat_aself = add fitvec;
    local nnm_aself = add not fitvec;

	// donor point

    idx = second prox_find [proxk, *PFObj.feat_da_pos, 3.5];
    fitvec = bitand [2, get [ligtype, idx]];
    local nmat_dself = length pack fitvec;
    local nnm_dself = add not fitvec;

	// dual point

    idx = second prox_find [proxk, *PFObj.feat_ba_pos, 3.5];
    occtype = get [ligtype, idx];
    local nmat_bself = length pack occtype;
    local nnm_bself = add not occtype;

	// metal

    idx = second prox_find [proxk, *PFObj.feat_met_pos, 3.5];
    fitvec = bitand [2, get [ligtype, idx]];
    local nmat_met = length pack fitvec;
    local nnm_met = add not fitvec;

    prox_close proxk;

    return (
	2.0 * (nmat_apj + nmat_dpj + nmat_bpj)
	- (nnm_apj + nnm_dpj + nnm_bpj)
	+ (nmat_aself + nmat_dself + nmat_bself)
	- (nnm_aself + nnm_dself + nnm_bself)
	+ 6.0 * nmat_met - 3.0 * nnm_met
    );
endfunction

local function HBPF_Score [ligpos, lig, rec]
    return ProjFeatScore [ligpos, *lig.HBPF_nhatype, *rec.HBPFObj];
endfunction

function dock_AlphaScoreOpen, dock_AlphaScoreClose, dock_AlphaScore;

global function dock_score_AlphaHB [cmd, arg, opt]
    const DEFAULTS = [
	was: 10.0,
	whbpf: 4.0
    ];

    const PANEL = [
	Mbox : [
	    Text : [
		title: 'Alpha Weight', name: 'was', type: 'real',
		bubbleHelp: 'Weight for the alpha score.'
	    ],
	    Text : [
		title: 'Hydrogen Bond Weight', name: 'whbpf', type: 'real',
		bubbleHelp: 'Weight for hydrogen bond projected feature score.'
	    ]
	]
    ];

    if cmd === 'ID' then
	return 'Alpha HB';

    elseif cmd === 'configpanelwidgets' then
	return PANEL;

    elseif cmd === 'configpanelevent' then
	return 0;

    elseif cmd === 'configvalues' then
	arg = tagcat [arg, DEFAULTS];
	return arg;

    elseif cmd === 'openReceptor' then
	HBPF_OpenReceptor arg;
	return;

    elseif cmd === 'closeReceptor' then
	HBPF_CloseReceptor arg;
	return;

    elseif cmd === 'openLigand' then
	HBPF_InitLigand arg;
	return;

    elseif cmd === 'closeLigand' then
	return;

    elseif not (cmd === 'score') then
//	exit twrite ['[Alpha HB scoring]: unknown command {}', cmd];
	return;

    endif

    opt = tagcat [opt, DEFAULTS];

    local [ligpos, lig, rec] = arg;
    local asctx  = dock_AlphaScoreOpen [rec, lig];
    local ascore = dock_AlphaScore [ligpos, lig, rec, asctx];

//  local ligposnh = apt get [ligpos, [*lig.Qindex]];
    local ligposnh = apt get [ligpos, [asctx.lQindex]];
    local h        = - HBPF_Score [ligposnh, lig, rec];

    dock_AlphaScoreClose asctx;
    return opt.was * ascore + opt.whbpf * h;
endfunction


// dock_AlphaScoreOpen score
const MAXRADIUS = 2.2;
const AS_RCUTOFF = 2 * MAXRADIUS;
const AS_ACUTOFF = 3.0;

function dock_ReceptorOpen;
function dock_LigandOpen;
function dock_LigandClose;
function dock_ReceptorClose;
function dock_score_AlphaHB;


global function dock_AlphaScore2 [ligpos, lig, rec, ctx]
    local rpos = ctx.rpos;
    local lpos = apt get [ligpos, [ctx.lQindex]];

    local [seg,idx,r2] = prox_find [ctx.rprox, lpos, ctx.lrad];
    local rsum = stretch [ctx.lrad,seg] + ctx.rrad[idx];

	// add in the repulsion term

    local REPSCALE = 1.0;
    local v = maxE[ 0, minE[ 1, 1 - sqrt r2 / rsum ]];
    local t = REPSCALE * add (cube v * (6 * sqr v - 15 * v + 10));

	// add in alpha attraction term

    local [aseg,aidx,ar2] = prox_find [ctx.aprox, lpos, 0];
//  local a = -0.04 * 1.85 * 1.85 * add exp (-0.5 * ar2);
    local a = -0.20 * 1.85 * 1.85 * add exp (-0.5 * s_min [ar2, pack aseg]);

//  local a = -0.2 * add (1 - sqrt ar2 * inv AS_ACUTOFF);

    return [t + a, t, a];
endfunction


local function main []
	Open './complex_prep.moe';
	local c = Chains[];
	local c0 = cAtoms c;
	local l0 = length c0;
	local lig = c0(l0);
	local i = 0;
	local rec = [];
	while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop

	//Close comen;

	local ligen = db_ReadColumn [ 'dock.mdb', 'mol' ];
	local cl = mol_Create ligen(1);	
	local cl0 = cAtoms cl;
	lig = cl0(1);	

	local rec_dock = dock_ReceptorOpen [rec, [aPos lig],[],[]];
	local lig_dock = dock_LigandOpen [lig, []];
	local asctx = dock_AlphaScoreOpen [rec_dock, lig_dock];
	local alphascore = dock_AlphaScore2 [aPos lig, lig_dock, rec_dock, asctx];	
	dock_score_AlphaHB ['openReceptor', rec_dock, []];
	dock_score_AlphaHB ['openLigand', lig_dock, []];
	local alphaHB_score = dock_score_AlphaHB ['score', [aPos lig,lig_dock,rec_dock], [ was : 10,	whbpf : 4 ]];
	local repulsion_term = alphascore[2];
	local attraction_term = alphascore[3];
	local ligposnh = apt get [aPos lig, [asctx.lQindex]];
	local HB_term = - HBPF_Score [ligposnh, lig_dock, rec_dock];
	dock_AlphaScoreClose asctx;
	dock_score_AlphaHB ['closeReceptor', rec_dock, []];
	dock_score_AlphaHB ['closeLigand', lig_dock, []];	
	dock_ReceptorClose rec_dock;
	dock_LigandClose lig_dock;
	//print [alphaHB_score, repulsion_term, attraction_term, HB_term];
	fwrite [ 'log.txt', '{},{},{},{}', alphaHB_score, repulsion_term, attraction_term, HB_term]; 
endfunction
'''
        self.ase_dock = '''#svl
function DockAtoms, DockFile;
function DockMDBwAtoms, DockMDBwFile;
local function RunDock[outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file]
    if isnull outf then
    outf = 'dock.mdb';
    endif;
    local psys = SystemPush [];

    Open './complex_prep.moe';

    local c = Chains[];
    local c0 = cAtoms c;
    local l0 = length c0;
    local site = c0(l0);
    local i = 0;
    local rec = [];
    local lig;
    while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop

    local opt = [
    outrmsd: 1,
    sel_ent_only_rec: 0,
    sel_ent_only: 0,
    wall: [ '', 0, [ 0, 0, 0 ], [ 1000000, 1000000, 1000000 ], 0 ],
    csearch: 1,
    confGenMethod: 'Rotate Bonds',
    placement: 'None',
    placement_opt: [  ],
    scoring: 'None',
    scoring_opt: [  ],
    dup_placement: 1,
    maxpose: 1,
    refine: 'Rigid Receptor',
    refine_opt: [ fixrec : 'Fix',	ed_map : 'Fo',	ed_f : 'Simulated',	ed_phi : 'Simulated',	ed_f2 : 'Simulated',	ed_path : '',	ed_res : 2.5,	ed_sfdata : [ [  ], [  ], [  ], 'Simulated', 'Simulated', 'Fo' ],	ed_surflevelD : 3,	cutoff : 6,	wholeres : 1,	mmgbvi : 1,	packsidechains : 1,	rigidlig : 0,	tether : 10,	gtest : 0.01,	maxit : 500,	OverrideSetup : 1,	k_potl : 100,	roffset : 0.4 ],
    rescoring: 'ASE',
    rescoring_opt: [ weight : 0.035 ],
    dup_refine: 1,
    remaxpose: 1,
    descexpr: '',
    receptor_mfield: '',
    ligand_mfield: 'mol',
    rxnFile: '',
    rxsite: [  ],
    edsupport: 1,
    ed_data: [ ed_dockpath : '' ],
    check_pose_geom: [  ],
    multiLigand: 0,
    need_dmat: 1,
    gen_plif: 1,
    tempDB: '/tmp/1iep_ligand_prep98896.mdb',
    ph4: ph4file,
    ligmdbname: ligmdbname,
    recmdbname: recmdbname,
    BatchFile: 'dock_batch.svl'
    ];
    pot_Load '$MOE/lib/Amber10EHT.ff.gz';

    pot_Setup [
    strEnable: 1,
    angEnable: 1,
    stbEnable: 1,
    oopEnable: 1,
    torEnable: 1,
    vdwEnable: 1,
    eleEnable: 1,
    solEnable: 0,
    resEnable: 1,
    strWeight: 1,
    angWeight: 1,
    stbWeight: 1,
    oopWeight: 1,
    torWeight: 1,
    vdwWeight: 1,
    eleWeight: 1,
    solWeight: 1,
    resWeight: 1,
    cutoffEnable: 1,
    cutoffOn: 8,
    cutoffOff: 10,
    eleDist: 2,
    vdwScale14: 0.5,
    vdwBuffer1: 0,
    vdwBuffer2: 0,
    eleScale14: 0.833333,
    eleDielectric: 1,
    eleBuffer: 0,
    solDielectric: 80,
    solDielectricOffset: 0,
    state0: 1,
    state1: 0,
    state2: 1,
    threadCount: 0
    ];
    if isnull ligmdb and isnull recmdb then
    DockAtoms [rec, site, lig, outf, opt];
    elseif length ligmdb and isnull recmdb then
    DockFile [rec, site, ligmdb, outf, opt];
    elseif isnull ligmdb and length recmdb then
    DockMDBwAtoms [recmdb, site, lig, outf, opt];
    elseif length ligmdb and length recmdb then
    DockMDBwFile [recmdb, site, ligmdb, outf, opt];
    endif
    write ['Docking finished at {}.', asctime []];
    SystemPop psys;
endfunction;

global argv;
function ArgvPull;

local function main []
    ArgvReset ArgvExpand argv;
    local [recmdb, ligmdb, ph4file, outf] = ArgvPull [
    ['-rec', '-lig','-ph4','-o'],
    1
    ];

    local ligmdbname, recmdbname;
    local delfile1 = [];
    local delfile2 = [];
    local delfile3 = [];

    local handle;
    if isnull outf then outf = 'dock.mdb'; endif
    if not isnull ligmdb then
    handle = _db_Open [ligmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read ligand mdb file {}', ligmdb];
    endif
    db_Close handle;
    ligmdbname = ligmdb;
    endif

    if not isnull recmdb then
    handle = _db_Open [recmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read receptor mdb file {}', recmdb];
    endif
    db_Close handle;
    recmdbname = recmdb;
    endif

    if not isnull ph4file then
    if not (ftype ph4file == 'file') then
        exit twrite ['File {} not found', ph4file];
    endif
    endif

    local oldttl = task_title -1;
    if second task_wfork [master: 'none'] == 'child' then
    task_settitle [-1, twrite ['Dock: [{}] {}',
        second app token wordsplit[ string oldttl, "' "], outf]];
    RunDock [outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file];
    exit [];
    endif

    fdelete delfile1;
    fdelete delfile2;
    fdelete delfile2;
endfunction
'''
        self.GBVIWSAdG_dock = '''#svl
function DockAtoms, DockFile;
function DockMDBwAtoms, DockMDBwFile;
local function RunDock[outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file]
    if isnull outf then
    outf = 'dock.mdb';
    endif;
    local psys = SystemPush [];


    Open './complex_prep.moe';

    local c = Chains[];
    local c0 = cAtoms c;
    local l0 = length c0;
    local site = c0(l0);
    local i = 0;
    local rec = [];
    local lig;
    while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop

    local opt = [
    outrmsd: 1,
    sel_ent_only_rec: 0,
    sel_ent_only: 0,
    wall: [ '', 0, [ 0, 0, 0 ], [ 1000000, 1000000, 1000000 ], 0 ],
    csearch: 1,
    confGenMethod: 'Rotate Bonds',
    placement: 'None',
    placement_opt: [  ],
    scoring: 'None',
    scoring_opt: [  ],
    dup_placement: 1,
    maxpose: 1,
    refine: 'Rigid Receptor',
    refine_opt: [ fixrec : 'Fix',	ed_map : 'Fo',	ed_f : 'Simulated',	ed_phi : 'Simulated',	ed_f2 : 'Simulated',	ed_path : '',	ed_res : 2.5,	ed_sfdata : [ [  ], [  ], [  ], 'Simulated', 'Simulated', 'Fo' ],	ed_surflevelD : 3,	cutoff : 6,	wholeres : 1,	mmgbvi : 1,	packsidechains : 1,	rigidlig : 0,	tether : 10,	gtest : 0.01,	maxit : 500,	OverrideSetup : 1,	k_potl : 100,	roffset : 0.4 ],
    rescoring: 'GBVI/WSA dG',
    rescoring_opt: [  ],
    dup_refine: 1,
    remaxpose: 1,
    descexpr: '',
    receptor_mfield: '',
    ligand_mfield: 'mol',
    rxnFile: '',
    rxsite: [  ],
    edsupport: 1,
    ed_data: [ ed_dockpath : '' ],
    check_pose_geom: [  ],
    multiLigand: 0,
    need_dmat: 1,
    gen_plif: 1,
    tempDB: '/tmp/1iep_ligand_prep98896.mdb',
    ph4: ph4file,
    ligmdbname: ligmdbname,
    recmdbname: recmdbname,
    BatchFile: 'dock_batch.svl'
    ];
    pot_Load '$MOE/lib/Amber10EHT.ff.gz';

    pot_Setup [
    strEnable: 1,
    angEnable: 1,
    stbEnable: 1,
    oopEnable: 1,
    torEnable: 1,
    vdwEnable: 1,
    eleEnable: 1,
    solEnable: 0,
    resEnable: 1,
    strWeight: 1,
    angWeight: 1,
    stbWeight: 1,
    oopWeight: 1,
    torWeight: 1,
    vdwWeight: 1,
    eleWeight: 1,
    solWeight: 1,
    resWeight: 1,
    cutoffEnable: 1,
    cutoffOn: 8,
    cutoffOff: 10,
    eleDist: 2,
    vdwScale14: 0.5,
    vdwBuffer1: 0,
    vdwBuffer2: 0,
    eleScale14: 0.833333,
    eleDielectric: 1,
    eleBuffer: 0,
    solDielectric: 80,
    solDielectricOffset: 0,
    state0: 1,
    state1: 0,
    state2: 1,
    threadCount: 0
    ];
    if isnull ligmdb and isnull recmdb then
    DockAtoms [rec, site, lig, outf, opt];
    elseif length ligmdb and isnull recmdb then
    DockFile [rec, site, ligmdb, outf, opt];
    elseif isnull ligmdb and length recmdb then
    DockMDBwAtoms [recmdb, site, lig, outf, opt];
    elseif length ligmdb and length recmdb then
    DockMDBwFile [recmdb, site, ligmdb, outf, opt];
    endif
    write ['Docking finished at {}.', asctime []];
    SystemPop psys;
endfunction;

global argv;
function ArgvPull;

local function main []
    ArgvReset ArgvExpand argv;
    local [recmdb, ligmdb, ph4file, outf] = ArgvPull [
    ['-rec', '-lig','-ph4','-o'],
    1
    ];

    local ligmdbname, recmdbname;

    local delfile1 = [];
    local delfile2 = [];
    local delfile3 = [];

    local handle;
    if isnull outf then outf = 'dock.mdb'; endif
    if not isnull ligmdb then
    handle = _db_Open [ligmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read ligand mdb file {}', ligmdb];
    endif
    db_Close handle;
    ligmdbname = ligmdb;
    endif

    if not isnull recmdb then
    handle = _db_Open [recmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read receptor mdb file {}', recmdb];
    endif
    db_Close handle;
    recmdbname = recmdb;
    endif

    if not isnull ph4file then
    if not (ftype ph4file == 'file') then
        exit twrite ['File {} not found', ph4file];
    endif
    endif

    local oldttl = task_title -1;
    if second task_wfork [master: 'none'] == 'child' then
    task_settitle [-1, twrite ['Dock: [{}] {}',
        second app token wordsplit[ string oldttl, "' "], outf]];
    RunDock [outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file];
    exit [];
    endif

    fdelete delfile1;
    fdelete delfile2;
    fdelete delfile2;
endfunction
'''
        self.GBVIWSAdG_score = '''#svl
//	dock_gbviwsa.svl			GBVI/WSA scoring function
//
//	13-nov-2015 (em) Add "Amber*" check
//	14-nov-2011 (cc) Add PFROSST FF parameters for scoring
//	30-jun-2011 (cc) Add AMBER FF trained parameters for scoring function
//	03-dec-2010 (cc) created
//
// COPYRIGHT (C) 2006-2017 CHEMICAL COMPUTING GROUP ULC ("CCG").
// ALL RIGHTS RESERVED.
//
// PERMISSION TO USE, COPY, MODIFY AND DISTRIBUTE THIS SOFTWARE IS HEREBY
// GRANTED PROVIDED THAT: (1) UNMODIFIED OR FUNCTIONALLY EQUIVALENT SOFTWARE
// DERIVED FROM THIS SOFTWARE MUST CONTAIN THIS NOTICE; (2) ALL CODE DERIVED
// FROM THIS SOFTWARE MUST ACKNOWLEDGE THE AUTHOR(S) AND INSTITUTION(S); (3)
// THE NAMES OF THE AUTHOR(S) AND INSTITUTION(S) NOT BE USED IN ADVERTISING
// OR PUBLICITY PERTAINING TO THIS SOFTWARE WITHOUT SPECIFIC WRITTEN PRIOR
// PERMISSION; (4) ALL CODE DERIVED FROM THIS SOFTWARE BE EXECUTED WITH THE
// MOLECULAR OPERATING ENVIRONMENT LICENSED FROM CCG.
//
// CCG DISCLAIMS ALL WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING
// ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS, AND IN NO EVENT
// SHALL CCG BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL DAMAGES OR
// ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER
// IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT
// OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//
#set title	'Dock GBVI/WSA Scoring Function'
#set version	'2015.11'
#set class	'MOE:simulation'

function SphereSurfaceArea;
function _Atoms;
function dock_ReceptorOpen, dock_ReceptorClose;
function dock_LigandOpen, dock_LigandClose;

// The energy of the receptor should be calculated when the receptor
// opened. This can not be currently done because of a bug
// when calculating the energy of a HEME group in the
// presence of a inert ligand. To enable chaching of energy data
// set CACHEDATA = 1.

const CACHEDATA = 0;

// All Constants are in the following format
// [Values using MMFF, AMBER]

const VDWCONST	= [0.12720,0.09150];
const ELECONST	= VDWCONST*(2/3);
const SACONST	= [0.01857,0.01244];
const INTERCEPT	= [-1.64247,-2.09251];
const DOTDEN	= 1;

// pocket computes the pocket score

local function clamp [x, lo, hi] = maxE [lo, minE [hi, x]];

local function pocket [pos, rad, MAXD]

    const R_WATER = 1.4;		// "water" probe defines interior
    const R_AVG = 1.2;		// "average" atom radius
    const R_EXT = 8.0;		// "ligand" probe defines exterior
    const W_EXTERIOR = sqr (R_AVG + R_EXT) - sqr (R_AVG + R_WATER);
    const M = 1024;			// round to 1/1024 of an Angstrom
    local [faces, links] = ialpha_complex3 [M*pos, sqr(M*rad + M*R_WATER)];

    local size = faces(1);				// #of vertices per face
    local weight = inv sqr M * faces(2);		// face weight
    local center = inv M * drop [faces, 2];		// face center
    local [f, g, attached] = links;			// face links

    // Extract simplices (index into faces) and vertices (index into faces)
    // Keep only links between smx1 and their vertices

    local smx = x_pack (size == max size);		// simplex id
    local smx1 = smx | weight[smx] >= 0;		// noninterior
    local smx2 = smx | weight[smx] > W_EXTERIOR;	// exterior

    local N = iadd (size <= 1);				// number of vertices
    local vtx1 = uniq f[x_join [g,smx1]];
    local vtx2 = uniq f[x_join [g,smx2]];

    [f,g] = [f,g] || [leE [1,f,N] and indexof [g, smx1]];

    // A_vs: for each vertex, list all simplices that contain the vertex
    // A_sv: for each smx1, list all vertices contained in the simplex
    // ff: for each face, if simplex from smx1, list all vertices it
    // contains A: list all vertices that share a simplex with the vertex

    local A_vs = apt get [[g], indicesof [igen N, f]];
    local A_sv = apt get [[f], indicesof [smx1, g]];
    local ff = put [[[]], smx1, A_sv];
    local A = app uniq app cat apt get [[ff], A_vs];

    // D: distance to the exterior
    // L: vertex-vertex distance

    MAXD = clamp [MAXD, 0.1, 50]; // De-facto infinity (veeery deeep)
    local D = rep [MAXD, N];	// Start with Inf everywhere
    D[pack vtx1] = MAXD-0.01;	// Non-interior vertices (almost Inf)
    D[pack vtx2] = 0;		// Exterior vertices
    local i, L = [];
    for i = N, 1, -1 loop
    local pos_i = apt peek [center, i];		// this vertex
    local pos_adj = apt get [center,[A(i)]];	// adjacent vertices
    L(i) = norm (pos_adj - pos_i);
    endloop

    // Calculate the weighted graph distance to the closest ext.vertex

    loop
    local new_D = minE [D, app min (apt get [[D], A] + L)];
    if D === new_D then break; endif
    D = new_D;
    endloop
    return inv MAXD * D;	// 0=exterior, 1=very deep
endfunction

const POTSETUP = [
    strEnable:0, angEnable:0, stbEnable:0, oopEnable:0, torEnable:0,
    vdwEnable:1, eleEnable:1, solEnable:1, resEnable:0,

    strWeight:1, angWeight:1, stbWeight:1, oopWeight:1,
    torWeight:1, vdwWeight:1, eleWeight:1, solWeight:1,

    cutoffEnable:1,
    cutoffOn:8,
    cutoffOff:10,

#   if 0			// !!! leave them as defaults
    vdwScale14:1,
    vdwBuffer1:0.07,		// !!! only valid for MMFF
    vdwBuffer2:0.12,
    eleScale14:0.75,
    eleBuffer:0.05,
#   endif

    eleDist:0,
    eleDielectric:1,
    solDielectric:80

];

// FUNCTION: 	Get Coulombic, Solvation and vdw terms.
// INPUT:	ff=1, MMFF, ff=2, AMBER
// RETURN:	[ Coulomb energy + solvation energy + van der Waals energy ]

local function GetEnergy ff 
    return first Potential [
    dX:0,
    setup:tagcat[
        [
        vdwWeight:VDWCONST(ff),
        eleWeight:ELECONST(ff),
        solWeight:ELECONST(ff)
        ],
        POTSETUP
    ]
    ];
endfunction

local function GetEnergy_d ff
    local setup = tagcat[
        [
            vdwWeight:VDWCONST(ff),
            eleWeight:ELECONST(ff),
            solWeight:ELECONST(ff)
        ],
        POTSETUP
    ];
    local old_setup = pot_Setup setup;
    local atoms = Atoms [];
    //local Ecoul = ELECONST(ff) * (first pot_eleEnergy atoms);
    //local Esol = ELECONST(ff) * (first pot_solEnergy atoms);
    //local Evdw = VDWCONST(ff) * (first pot_vdwEnergy atoms);
    local Ecoul = (first pot_eleEnergy atoms);
    local Esol = (first pot_solEnergy atoms);
    local Evdw = (first pot_vdwEnergy atoms);
    pot_Setup old_setup;
    return [Ecoul, Esol, Evdw];
endfunction

// closeReceptor and openReceptor handle the scoring function operators

local function openReceptor [rec, opt]
    local nskey = SystemOpen [];
    local prev_skey = SystemCurrent nskey;

    *rec.mmgbviwsa_atoms = cat oAtoms mol_Create *rec.mol;
    *rec.mmgbviwsa_prox = prox_open [
    3 + max [ 0.1, max aRadius *rec.mmgbviwsa_atoms],
    aPos *rec.mmgbviwsa_atoms,
    aRadius *rec.mmgbviwsa_atoms
    ];

    *rec.mmgbviwsa_saweight = pocket [
    aPos *rec.mmgbviwsa_atoms,
    aRadius *rec.mmgbviwsa_atoms,
    50
    ];

#   if CACHEDATA
    pot_SetCharges [];
    *rec.mmgbviwsa_Erec = GetEnergy opt;
#   endif
    *rec.mmgbviwsa_nskey = nskey;
    SystemCurrent prev_skey;
endfunction

local function closeReceptor [rec, opt]
    *rec.mmgbviwsa_atoms = [];
    *rec.mmgbviwsa_sa = [];
    *rec.mmgbviwsa_saweight = [];
    *rec.mmgbviwsa_Erec = [];
    prox_close *rec.mmgbviwsa_prox;
    *rec.mmgbviwsa_prox = [];
    SystemClose *rec.mmgbviwsa_nskey;
endfunction

// ------------------------------ SPECIAL ENTRY POINTS -----------------------

// ScoreGBVIWSA_dG takes receptor atoms and ligand atoms and returns an 
// estimate of the binding affinity using the GBVIWSA dG. This is a one-shot 
// calculation not intended for speed since the receptor is initialized each 
// time
// FUNCTION : 	Calculate the GBVIWSA score between recatoms and ligatoms
// INPUT: 	recatoms = receptor atoms
//		ligatoms = ligand atoms
//		ff = options (1 = MMFF parameters,2 = AMBER Parameters)
//		rec = if docking cached receptor energy terms.
// RETURN:	GBVIWSA score.
// NOTE:	This function must be used within a private system.

local function ScoreGBVIWSA_dG [recatoms, ligatoms, ff, rec]

    // Get initial inert state and set all atoms to inert. This accounts
    // for cases when there are additional active atoms in the system then
    // the ones being passed.

    local ainert = aInert Atoms [];
    aSetInert [Atoms [], 1];
    aSetInert [[recatoms,ligatoms], 0];

    // Get receptor terms. If rec is null then calculate terms otherwise
    // use cashed values.

    local Erec, SArec, weight, prox;
    prox = *rec.mmgbviwsa_prox;

    local [seg, at, d2] = prox_find [prox, aPos ligatoms, aRadius ligatoms];
    at = uniq at;

#   if CACHEDATA
    Erec = *rec.mmgbviwsa_Erec;
#   else
    aSetInert [[recatoms,ligatoms],[0,1]];
    Erec = GetEnergy ff;
    local [Erec_coul, Erec_sol, Erec_vdw] = GetEnergy_d ff;
#   endif

    weight = *rec.mmgbviwsa_saweight;
    SArec = SphereSurfaceArea [aPos recatoms[at], aRadius recatoms[at], DOTDEN];

    // Get Complex Energy Terms

    aSetInert [[recatoms,ligatoms], 0];
    local comatoms = cat [recatoms[at],ligatoms];
    local Ecom = GetEnergy ff;
    local [Ecom_coul, Ecom_sol, Ecom_vdw] = GetEnergy_d ff;
    local SAcom = SphereSurfaceArea [aPos comatoms, aRadius comatoms, DOTDEN];

    // Get ligand only energy terms

    aSetInert [recatoms,1];
    local Elig = GetEnergy ff;
    local [Elig_coul, Elig_sol, Elig_vdw] = GetEnergy_d ff;
    local SAlig = SphereSurfaceArea [aPos ligatoms, aRadius ligatoms, DOTDEN];

    // Since we are not using atoms any more reset to initial inert state

    aSetInert [Atoms [], ainert ];

    // if wsa == 2 then use weighted surface area

    local SA;
    local mask = m_join [comatoms,recatoms[at]];
    SAcom = add ((SAcom|mask) * weight[at]) + add (SAcom| not mask);
    SArec = add (SArec * weight[at]);
    SA = SAcom - SArec - (add SAlig);

    return [Ecom - Erec - Elig + SACONST(ff)*SA +INTERCEPT(ff),
            INTERCEPT(ff),
            Ecom_coul - Erec_coul - Elig_coul,
            Ecom_sol  - Erec_sol  - Elig_sol,
            Ecom_vdw  - Erec_vdw  - Elig_vdw,
            SACONST(ff)*SA
            ];
endfunction

// FUNCTION:	Main entry point for GBSA scoring function.
// INPUT:	cmd = what do you want the function to do.
//		arg = arguments for the cmd.
//		opt = cmd options
// Return:	GBSA score

global function dock_score_GBVIWSA_dG [cmd, arg, opt]

    local ligpos = [], lig = [], rec=[], score;

    if cmd === 'ID' then
    return 'GBVI/WSA dG';

    elseif cmd === 'configpanelwidgets' then
    return [];

    elseif cmd === 'configpanelevent' then	// arg = [val,trig]
    return 0;

    elseif cmd === 'configvalues' then		// arg = val
    return [];

    elseif cmd === 'openReceptor' then		// arg = rec
    openReceptor [arg, opt];
    return;

    elseif cmd === 'closeReceptor' then		// arg = rec
    closeReceptor [arg,opt];
    return;

    elseif cmd === 'openLigand' then		// arg = lig
    return;

    elseif cmd === 'closeLigand' then		// arg = lig
    return;

    elseif cmd === 'scoreLigX' then		// arg = [ recatoms, ligatoms]
    local old_prio = task_prio 0;
    local [recatoms,ligatoms] = arg;

    rec = dvar_open [];
    lig = dvar_open [];
    *rec.mol = mol_Extract recatoms;
    *lig.mol = mol_Extract ligatoms;

    dock_score_GBVIWSA_dG ['openReceptor', rec, opt];
    dock_score_GBVIWSA_dG ['openLigand', lig, opt];

    score = dock_score_GBVIWSA_dG [
        'score',
        [aPos ligatoms, lig, rec]
    ];

    dock_score_GBVIWSA_dG ['closeReceptor', rec, opt];
    dock_score_GBVIWSA_dG ['closeLigand', lig, opt];
    dvar_close rec;
    dvar_close lig;

    task_prio old_prio;
    return score;

    elseif not (cmd === 'score') then		// arg = [ligpos,lig,rec]
    return;
    endif

    [ligpos, lig, rec] = arg;

    local prev_skey = SystemCurrent *rec.mmgbviwsa_nskey;

    local ff=1;
    if length findmatch [['Amber*', 'AMBER*','PFROSST'],(pot_Info []).title]
    then
    ff = 2;
    endif

    // If it is the same ligand but different ligand coordinate use 
    // previously create atoms keys. If not the same ligand destroy
    // atom keys and create new ligand

    if neL [*rec.mmgbviwsa_ligmol,*lig.mol] then
    oDestroy oChains *rec.mmgbviwsa_ligatoms;
    *rec.mmgbviwsa_ligmol = *lig.mol;
    *rec.mmgbviwsa_ligatoms = cat oAtoms mol_Create *lig.mol;
    endif

    aSetPos [ *rec.mmgbviwsa_ligatoms, ligpos ];

    // As ScoreGBVIWSA_dG assumes it can muck with the inert state,
    // we will set the stub atoms Inert and exclude them from the
    // ligand passed in.  This will prevent it from making them uninert,
    // and should exclude them from FF calculations.  We'll restore the
    // original inert state upon completion.

    local uAtoms = *rec.mmgbviwsa_ligatoms, origInert = aInert uAtoms;
    if length *lig.ligData.useAtoms then
    uAtoms = uAtoms | *lig.ligData.useAtoms;
    aSetInert [*rec.mmgbviwsa_ligatoms, not *lig.ligData.useAtoms];
    endif

    score = ScoreGBVIWSA_dG [
    *rec.mmgbviwsa_atoms,
//	*rec.mmgbviwsa_ligatoms,
    uAtoms,
    ff,
    rec
    ];

    aSetInert [*rec.mmgbviwsa_ligatoms, origInert];
    SystemCurrent prev_skey;
    return score;
endfunction

local function main []
    Open './complex_prep.moe';
    local c = Chains[];
    local c0 = cAtoms c;
    local l0 = length c0;
    local lig = c0(l0);
    local i = 0;
    local rec = [];
    while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop

    //Close comen;

    local ligen = db_ReadColumn [ 'dock.mdb', 'mol' ];
    local cl = mol_Create ligen(1);	
    local cl0 = cAtoms cl;
    lig = cl0(1);		

    local rec_dock = dock_ReceptorOpen [rec, [aPos lig],[],[]];
    local lig_dock = dock_LigandOpen [lig, []];
    dock_score_GBVIWSA_dG ['openReceptor', rec_dock, []];
    dock_score_GBVIWSA_dG ['openLigand', lig_dock, []];	
    local GBVIWSA_dG = dock_score_GBVIWSA_dG ['score', [aPos lig,lig_dock,rec_dock], [  ]];
    dock_score_GBVIWSA_dG ['closeReceptor', rec_dock, []];
    dock_score_GBVIWSA_dG ['closeLigand', lig_dock, []];	
    dock_ReceptorClose rec_dock;
    dock_LigandClose lig_dock;
    //print GBVIWSA_dG;
    fwrite [ 'log.txt', '{},{},{},{},{},{}', GBVIWSA_dG[1], GBVIWSA_dG[2], GBVIWSA_dG[3], GBVIWSA_dG[4], GBVIWSA_dG[5], GBVIWSA_dG[6]]; 
endfunction
'''
        self.LondondG_dock = '''#svl
function DockAtoms, DockFile;
function DockMDBwAtoms, DockMDBwFile;
local function RunDock[outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file]
    if isnull outf then
    outf = 'dock.mdb';
    endif;
    local psys = SystemPush [];

    Open './complex_prep.moe';

    local c = Chains[];
    local c0 = cAtoms c;
    local l0 = length c0;
    local site = c0(l0);
    local i = 0;
    local rec = [];
    local lig;
    while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop

    local opt = [
    outrmsd: 1,
    sel_ent_only_rec: 0,
    sel_ent_only: 0,
    wall: [ '', 0, [ 0, 0, 0 ], [ 1000000, 1000000, 1000000 ], 0 ],
    csearch: 1,
    confGenMethod: 'Rotate Bonds',
    placement: 'None',
    placement_opt: [  ],
    scoring: 'None',
    scoring_opt: [  ],
    dup_placement: 1,
    maxpose: 1,
    refine: 'Rigid Receptor',
    refine_opt: [ fixrec : 'Fix',	ed_map : 'Fo',	ed_f : 'Simulated',	ed_phi : 'Simulated',	ed_f2 : 'Simulated',	ed_path : '',	ed_res : 2.5,	ed_sfdata : [ [  ], [  ], [  ], 'Simulated', 'Simulated', 'Fo' ],	ed_surflevelD : 3,	cutoff : 6,	wholeres : 1,	mmgbvi : 1,	packsidechains : 1,	rigidlig : 0,	tether : 10,	gtest : 0.01,	maxit : 500,	OverrideSetup : 1,	k_potl : 100,	roffset : 0.4 ],
    rescoring: 'London dG',
    rescoring_opt: [  ],
    dup_refine: 1,
    remaxpose: 1,
    descexpr: '',
    receptor_mfield: '',
    ligand_mfield: 'mol',
    rxnFile: '',
    rxsite: [  ],
    edsupport: 1,
    ed_data: [ ed_dockpath : '' ],
    check_pose_geom: [  ],
    multiLigand: 0,
    need_dmat: 1,
    gen_plif: 1,
    tempDB: '/tmp/1iep_ligand_prep98896.mdb',
    ph4: ph4file,
    ligmdbname: ligmdbname,
    recmdbname: recmdbname,
    BatchFile: 'dock_batch.svl'
    ];
    pot_Load '$MOE/lib/Amber10EHT.ff.gz';

    pot_Setup [
    strEnable: 1,
    angEnable: 1,
    stbEnable: 1,
    oopEnable: 1,
    torEnable: 1,
    vdwEnable: 1,
    eleEnable: 1,
    solEnable: 0,
    resEnable: 1,
    strWeight: 1,
    angWeight: 1,
    stbWeight: 1,
    oopWeight: 1,
    torWeight: 1,
    vdwWeight: 1,
    eleWeight: 1,
    solWeight: 1,
    resWeight: 1,
    cutoffEnable: 1,
    cutoffOn: 8,
    cutoffOff: 10,
    eleDist: 2,
    vdwScale14: 0.5,
    vdwBuffer1: 0,
    vdwBuffer2: 0,
    eleScale14: 0.833333,
    eleDielectric: 1,
    eleBuffer: 0,
    solDielectric: 80,
    solDielectricOffset: 0,
    state0: 1,
    state1: 0,
    state2: 1,
    threadCount: 0
    ];
    if isnull ligmdb and isnull recmdb then
    DockAtoms [rec, site, lig, outf, opt];
    elseif length ligmdb and isnull recmdb then
    DockFile [rec, site, ligmdb, outf, opt];
    elseif isnull ligmdb and length recmdb then
    DockMDBwAtoms [recmdb, site, lig, outf, opt];
    elseif length ligmdb and length recmdb then
    DockMDBwFile [recmdb, site, ligmdb, outf, opt];
    endif
    write ['Docking finished at {}.', asctime []];
    SystemPop psys;
endfunction;

global argv;
function ArgvPull;

local function main []
    ArgvReset ArgvExpand argv;
    local [recmdb, ligmdb, ph4file, outf] = ArgvPull [
    ['-rec', '-lig','-ph4','-o'],
    1
    ];

    local ligmdbname, recmdbname;

    local delfile1 = [];
    local delfile2 = [];
    local delfile3 = [];

    local handle;
    if isnull outf then outf = 'dock.mdb'; endif
    if not isnull ligmdb then
    handle = _db_Open [ligmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read ligand mdb file {}', ligmdb];
    endif
    db_Close handle;
    ligmdbname = ligmdb;
    endif

    if not isnull recmdb then
    handle = _db_Open [recmdb, 'read'];
    if handle == 0 then
        exit twrite ['Cannot read receptor mdb file {}', recmdb];
    endif
    db_Close handle;
    recmdbname = recmdb;
    endif

    if not isnull ph4file then
    if not (ftype ph4file == 'file') then
        exit twrite ['File {} not found', ph4file];
    endif
    endif

    local oldttl = task_title -1;
    if second task_wfork [master: 'none'] == 'child' then
    task_settitle [-1, twrite ['Dock: [{}] {}',
        second app token wordsplit[ string oldttl, "' "], outf]];
    RunDock [outf, recmdb, recmdbname, ligmdb, ligmdbname, ph4file];
    exit [];
    endif

    fdelete delfile1;
    fdelete delfile2;
    fdelete delfile2;
endfunction
'''
        self.LondondG_score = '''#svl
function dock_ReceptorOpen, dock_ReceptorClose;
function dock_LigandOpen, dock_LigandClose;
function pboltz_AtomParameters;
function ph4_aPlanar, ph4_aMetalLigand;

const CUTOFF    = 8.0;			// main interaction cutoff distanced
const MCUTOFF   = 3.0;			// maximum metal contact distance 
const AROCUTOFF = 5.0;			// aro-alk cutoff distance

// ------------------------------ IONIC TERM ----------------------------------
// !!! NOT USED !!! EXPERIMENTAL

local function FormalCharge atoms
    const FCHARGE = tr [
    [  0.0, '[#1]'				],
    [  0.0, '[CX4]'				],

#if 1
    [  0.0, '[NX3!i][C+0]=[N+!$(*[-*])]'		],	// NCN+ res
    [  0.0, '[N+!$(*[-*])]=[C+0][NX3!i]'		],
    [  1.0, '[C+0](=[N+!$(*[-*])])[NX3!i]'		],
    [  1.0, '[C+1]([NX3!i])[NX3!i]'			],

//	[  0.0, '[n;i1+!$(*[-*])]:[c+0]:[n;i2]'		],	// ncn+ res
//	[  0.0, '[n;i2]:[c+0]:[n;i1+!$(*[-*])]'		],
//	[  1.0, '[c+0](:[n;i2]):[n;i1+!$(*[-*])]'	],

    [  0.0, '[N;X3!i][c+0]([n;i1+!$(*[-*])])'	],
    [  0.0, '[n;i1+!$(*[-*])][c+0][N;X3!i]'		],
    [  1.0, '[c+0]([n;i1+!$(*[-*])])[N;X3!i]'	],

    [  0.0, '[OX1][C+0]=[OX1]'			],	// -CO2
    [  0.0, '[OX1]=[C+0][OX1]'			],
    [ -1.0, '[C+0](=[OX1])-[OX1]'			],

    [ -2.0, '[SX4]([OX1])([OX1])([OX1])[OX1]'	],	// -SO4
    [ -1.0, '[SX4]([OX1])([OX1])[OX1]'		],	// -SO3
    [  0.0, '[OX1][SX4]([OX1])[OX1]'		],

    [ -1.0, '[PX4]([OX1])([OX1])[OX1]'		],	// -PO3/-PO4
    [ -1.0, '[PX4]([OX1])[OX1]'			],	// -PO2
    [  0.0, '[OX1][PX4][OX1]'			],

    [  1.0, '[NX4!i!$(*~[-])]'			],	// N+
    [  1.0, '[#7+;iX3!$(*~[-])]'			],	// n+

#else
    [ -0.9, '[OX1]-[C+0]=[OX1]'		],		// -COO
    [ -0.9, '[OX1]=[C+0]-[OX1]'		],
    [  0.8, '[C+0](=[OX1])-[OX1]'		],

    [  1.0, '[NX4!i]'			],		// N+

//	[ -0.5,	'[OX1]=[C+0]'			],		// >C=O
//	[  0.5,	'[C+0][OX1]'			],		// >C=O

    [ -0.9, '[OX1][SX4]([OX1])([OX1])[OX1]'		],	// -SO4
    [  1.6, '[SX4]([OX1])([OX1])([OX1])[OX1]'	],
    [ -0.8, '[OX1][SX4]([OX1])[OX1]'		],	// -SO3
    [  1.4, '[SX4]([OX1])([OX1])[OX1]'		],
    [ -0.6, '[OX1][SX4][OX1]'			],	// -SO2
    [  1.2, '[SX4]([OX1])[OX1]'			],

    [ -1.0, '[OX1][PX4]([OX1])[OX1]'		],	// -PO3
    [  1.0, '[PX4]([OX1])([OX1])[OX1]'		],
    [ -1.0, '[OX1][PX4][OX1]'			],	// -PO2
    [  1.0, '[PX4]([OX1])[OX1]'			],
    [ -1.0, '[OX1][PX4]'				],	// -PO
    [  1.0, '[PX4][OX1]'				],

    [  1/8, '[N+!$(*[-*])]=[C+0]([NX3!i])[NX3!i]'	 ],	// N2CN+ res
    [  1/8, '[NX3!i][C+0](=[N+!$(*[-*])])[NX3!i]'	 ],
    [  5/8, '[C+0](=[N+!$(*[-*])])([NX3!i])[NX3!i]' ],

    [  1/5, '[NX3!i][C+0]=[N+!$(*[-*])]'		],	// NCN+ res
    [  1/5, '[N+!$(*[-*])]=[C+0][NX3!i]'		],
    [  3/5, '[C+0](=[N+!$(*[-*])])[NX3!i]'		],

    [  1/5, '[n;i1+!$(*[-*])]:[c+0]:[n;i2]'		],	// ncn+ res
    [  1/5, '[n;i2]:[c+0]:[n;i1+!$(*[-*])]'		],
    [  3/5, '[c+0](:[n;i2]):[n;i1+!$(*[-*])]'	],
#endif

    [  1.0, '[#G1X0]'			],
#if 0
    [  2.0, '[#G2X0]'			],
#else
    [  2.0, '[#G2,#T;X0]'			],
#endif
    [ -1.0, '[#G7X0]'			],

    [  0.0, '[#T]'				],
    [  1.0, '[#G1D0]'			],
#if 0
    [  2.0, '[#G2D0]'			],
#else
    [  2.0, '[#G2,#T;D0]'			],
#endif

    [ -1.0, '[-*!$(*~[+*])]'		],
    [  1.0, '[+*!$(*~[-*])]'		],

    [  0.0, '*'				]
    ];

    return FCHARGE(1)[sm_Indexof[atoms, FCHARGE(2)]];
endfunction

#ifnbif invxerfx
    local function invxerfx x
    local y = rep [2 / sqrt PI, length x];	// lim{x->0} = 2 / sqrt PI
    local idx = x_pack (abs x > 1e-8);
    x = x[idx];
    y[idx] = inv x * erf x;
    return y;
    endfunction
#endif

local function ReactionField [pos, q]
    const L = 6.0;
    local prox = prox_open [L, pos, L];

    local [seg, xB, r2] = prox_find [prox, pos, 0];
    local xA = stretch [igen l_length pos, seg];
    [xA,xB,r2] = [xA, xB, r2] || [r2 > 0.25];
    prox_close prox;

#if 0
    local r = sqrt r2;
    local eps = 1 + 60 * (1 - exp (-0.1 * r));
    local qq = q[xA] * q[xB];
    return add (qq * inv (r * eps));
#else
    local qq = q[xA] * q[xB];
    return add (qq * 0.25 / sqr (1 + sqrt r2));
#endif
endfunction

local function IonicTerm [rpos, rq, lpos, lq]
    const CUTOFF = 8.0;

    lpos = lpos || [lq <> 0];
    lq = lq | lq <> 0;

    rpos = rpos || [rq <> 0];
    rq = rq | rq <> 0;

    local prox = prox_open [CUTOFF, rpos, CUTOFF];
    local [seg,xB,r2] = prox_find [prox, lpos, 0];
    local xA = stretch [igen l_length lpos, seg];
    [xA,xB,r2] = [xA, xB, r2] || [r2 > sqr 0.25];
    prox_close prox;

    return COULOMB_SCALE * 0.25 * add (
    lq[xA] * rq[xB] * (inv sqr (1 + sqrt r2) - inv sqr (1 + CUTOFF))
    );

    return (
      ReactionField [apt cat [rpos,lpos], cat[rq,lq]]
    - ReactionField [rpos, rq]
    - ReactionField [lpos, lq]
    );
endfunction

// =========================== SCORING FUNCTION ===============================

// SolvationRadius accepts atoms and returns the solvation radius.  We
// take these to be the pboltz (OPLS-AA) radii except that C+ atoms in NCN+
// structures are given Car radii

local function SolvationRadius atoms
    local solR = first pboltz_AtomParameters atoms;
#if 1
    const CNN = tr [
    [ '[C+0](=[N+!$(*[i])])[N;X3!i]',	1	],
    [ '[C+1]([N;X3!i])[N;X3!i]',		1	],
    [ '*',					0	]
    ];
    solR | CNN(2)[sm_Indexof[atoms, CNN(1)]] = 1.992;
    return 0.5 + solR;
#else
    solR = solR * pow[2,1/6];
    return solR;
#endif
endfunction

// GetUdir gets the 'u' direction from a mol: the opposite of the
// average normalize heavy neighbor direction

local function GetUdir [mol, pos, idxlist]
    local el = mol(4)(MOL_ATOM_EL);
    local xbond = mol(4)(MOL_ATOM_BONDS);
    local n = l_length mol(4);

    local udir = rep [[], length idxlist];
    local idx;

    for idx = 1, length idxlist loop
    local i = idxlist(idx);
    local pi = apt peek [pos, i];
    local nbr = xbond(i) | not indexof [el[xbond(i)],['H','LP']];
    local pn = vnormalize (apt get [pos, [nbr]] - pi);
    udir(idx) = app add pn * invz l_length pn;
    endloop

    return - vnormalize tr udir;
endfunction

local function GetZdir [mol, pos, idxlist]
    local el = mol(4)(MOL_ATOM_EL);
    local xbond = mol(4)(MOL_ATOM_BONDS);
    local n = l_length mol(4);

    local zdir = rep [[], length idxlist];
    local j;

    for j = 1, length idxlist loop
    local i = idxlist(j);
    local nbr = xbond(i) | not indexof [el[xbond(i)],['H','LP']];

    if length nbr == 1 then
        i = nbr;
        nbr = xbond(i) | not indexof [el[xbond(i)],['H','LP']];
    endif
    if length nbr <= 1 or length nbr >= 4 then
        zdir(j) = [0,0,0];
        continue;
    endif

    local pi = apt peek [pos, i];
    local pn = vnormalize (apt get [pos, [nbr]] - pi);

    if length nbr == 2 then
        zdir(j) = vnormalize vcross tr pn;
    else
        zdir(j) = vnormalize app add vnormalize vcross [pn, app rotl pn];
    endif
    endloop

    return tr zdir;
endfunction

// -------------------------- FLEXIBILITY MEASURE -----------------------------

// AutomorphismCount counts the number of automorphisms of a molecule

local function logAutomorphismCount atoms
    local xbond = BondGraph atoms;
    local prio = aPrioZQH atoms;
    return graph_automorphism_logcount [xbond, prio];
endfunction

local function FlexibilityScore atoms
    atoms = atoms | aAtomicNumber atoms > 1;

    const PIATOM = [				// planar atoms
        '[i]'
    ,	'[N!i!X4]C=[#G6]'
    ,	'[N!i!X4]C=N'
    ,	'[O!i!X4]C=[#G6]'
    ,	'[C+][N;X3!i]'
    ];
    local atype = aHeavyValence atoms;
    atype | sm_Indexof [atoms, PIATOM] = 1;

    local A = stretch [atoms, aBondCount atoms];
    local B = cat aBonds atoms;

    [A,B] = [A,B] || [A < B];
    [A,B] = [A,B] || [not bInRing[A,B] and bOrder [A,B] == 1];
    [A,B] = [A,B] || [aHeavyValence A > 1 and aHeavyValence B > 1];

    function smi2 [pat1,pat2] = (
       sm_Match [pat1, A] and sm_Match [pat2, B]
    or sm_Match [pat2, A] and sm_Match [pat1, B]
    );

    [A,B] = [A,B] || not [
       smi2 [ '[N;!i]', '[C+0]=[#G6]'	]	// (thio)amide
    or smi2 [ '[N;!i]', '[C+0]=[N+]'	]	// NCN+ resonance
    or smi2 [ '[N;!i]', '[c+0]:[n+Q2]'	]	// NCN+ resonance
    or smi2 [ '[N;!i]', '[C+1][N;!i]'	]	// NCN+ resonance
    or smi2 [ '[O;X2]', '[C+0]=[#G6]'	]	// -O- in ester
    ];

    const PARAM = cat [
        [ 0.5226,	0.5070,		0.4606,		0.7870	]
    ,   [ 0.5070,	0.1382,		0.1614,		0.1986	]
    ,   [ 0.4606,	0.1614,		0.2000,		0.2000	]
    ,   [ 0.7870,	0.1986,		0.2000,		0.2000	]
    ];

    local xA = indexof [A, atoms];
    local xB = indexof [B, atoms];

    // We will pack [xA, xB] so as to remove zeros.  These can result
    // from (bogus) reactions creating atoms with more than 4 bonds that
    // are not considered unimportant.  We don't wish to beep...

    local m  = andE [xA, xB];
    [xA, xB] = [xA, xB] || [m];

#   if 0
    local data = [];
    local tA = totok atype[xA], tB = totok atype[xB];
    local stype = tok_cat [ minE [tA, tB], '-', maxE [tA, tB] ];
    if length stype then
        local [idx, mask] = sam stype;
        data = tagcat [data, tag [(stype[idx]|mask), mtoc mask]];
    endif
    return data;
#   endif

    return (
           add PARAM[4 * dec atype[xA] + atype[xB]]
    - 0.3361 * logAutomorphismCount atoms
    );
endfunction

// --------------------------- METAL LIGATION --------------------------------

// MetalLigandRadius is used to determine the covalent radii for atoms used
// for the metal ligation 9-6 potential.

local function MetalLigandRadius el
    local rad = el_COV_Radius el;
    const TABLE = tr [
    [ 'O',		0.77 ]
    ,	[ 'N',		0.75 ]
    ,	[ 'S',		1.00 ]
    ,	[ 'P',		1.00 ]
    ,	[ 'Zn',		1.38 ]	// Transition
    ,	[ 'Mn',		1.38 ]
    ,	[ 'Fe',		1.38 ]
    ,	[ 'Ni',		1.38 ]
    ,	[ 'Cu',		1.38 ]
    ,	[ 'Cr',		1.38 ]
    ,	[ 'Na',		1.53 ]	// Group I
    ,	[ 'K',		1.76 ]
    ,	[ 'Mg',		1.43 ]	// Group II
    ,	[ 'Ca',		1.81 ]
    ];
    local idx = indexof [el, TABLE(1)];
    rad | idx = TABLE(2)[pack idx];
    return rad;
endfunction

// The metal ligation parameter table encodes ligand preferred bond angles
// with its 'u' direction (opposite to bonds).

const MLIG_PARAM = tr [				// !!! VERIFY
//			ang	aisd
//			---	----
    [ '[Q0]',		  0,	  0	]

,   [ '[Q1]#*',		  0,	1/15	]
,   [ '[Q1]=*',		 60,	1/15	]	// S=, O=, P=
,   [ '[Q1]*=[#G6]',	 60,	1/15	]	// O-C=O
,   [ '[Q2;i]',		  0,	1/15	]	// O-C=O

,   [ '[OQ1][i]',	 60,	1/15	]
,   [ '[OQ1]',		 85,	1/15	]
,   [ '[SQ1]',		 90,	1/15	]

,   [ '[Q2]',		 35,	1/15	]
,   [ '[Q3]',		  0,	1/15	]
,   [ '[Q1]',		 70,	1/15	]

,   [ '*',		  0,	  0	]
];

// The metal ligation term measures a number of things using heavy atoms:
//
//	1) distance potential between ligator and metal
//	2) angular environment of ligator
//	3) bond angles about the metal center (sd3 hybridization)

local function MetalContactScore [lpos, lig, rec, opt]
    const kT = KBOLTZ * 300, inv_kT = 1/kT;

    // obtain the list of candidate metal ligations

    local xL = *lig.LdG_Mxml;
    local [seg, xM, r] = prox_find [*rec.LdG_Mprox, apt get [lpos,[xL]], 0];
    if not l_length [xM, r] then return [0,0]; endif

    xL = stretch [xL, seg];
    r = sqrt r;

    // calculate the log of the distance score using a 9-6 potential
    // which approximates the contact statistics radial distribution
    // !!!! WHAT CONSTANT TO USE FOR 2.485

    local r0 = *rec.LdG_Mrad[xM] + *lig.LdG_Mrad[xL];
    local x = r0 / maxE [r0, r], x3 = cube x, x6 = sqr x3;
    local logS = ( -inv_kT * 2.485 * (x6 * (2 * x3 - 3) + 1) );
    local bscore = exp logS;

    // make sure that the ligating atom makes a reasonable C-O..M angle
    // there shouldn't be very many contacts so we can be slightly

    local xR = *rec.LdG_Mxmet[xM];		// convert to mol indices
    local rmol = *rec.mol, rpos = mol_aPos rmol;
    local mpos = apt get [rpos, [xR]];
    local mligpos = apt get [lpos, [xL]];

    local mpx = *lig.LdG_Mpx[xL];
    local ang0 = MLIG_PARAM(2)[mpx], aisd = MLIG_PARAM(3)[mpx];

    local udir = GetUdir [*lig.mol, lpos, xL];
    local hvec = mpos - mligpos;
    local ang = (180/PI) * acos maxE[-1, minE[1, add (hvec * udir) * invz r]];

    logS = logS - 0.5 * sqr ( (ang - ang0) * aisd );

//    return [add exp logS, 0];

    // evaluate the valbond angle energy for the metal evaluated
    // add sd3 hybridization (include ligand - M - receptor interactions)

    local function sd3Energy cos_a		// valbond (s=1, p=0, d=3)
    local f = 0.25 * (1 + 3/2 * (3 * sqr cos_a - 1));
    return 0.5 * (1+sqrt(5*3)) * (
        1 - sqrt 0.5 * sqrt(1 + sqrt(1 - sqr f))
    );
    endfunction

    function MetalEnergy [mpos, rnbr, lnbr, lS]
    local i, j;
    local E = 0;

    for i = 1, length lnbr loop
        local u = vnormalize (apt peek [lpos, lnbr(i)] - mpos);
        local v = vnormalize (apt get [rpos, [rnbr]] - mpos);
        E = E + lS(i) * add sd3Energy add (u * v);
        v = vnormalize (apt get [lpos, [drop[lnbr,i]]] - mpos);
        E = E + add ( sd3Energy add (u * v) * sqrt (lS(i) * drop[lS,i]) );
    endloop

    return E;
    endfunction

    local i, j;
    local angEnergy = 0;

    local [Midx,Mmask] = sam xR, Mdeg = mtoc Mmask;
    local Mlnbr = split [xL[Midx], Mdeg];
    local Mlsco = split [(bscore)[Midx], Mdeg];
    local Mbond = *rec.LdG_Mxnbr[xM];

    for i = 1, length xR loop
    angEnergy = angEnergy + MetalEnergy [
         apt peek [mpos, i], Mbond(i), Mlnbr(i), Mlsco(i)
    ];
    endloop

    return [add exp logS, angEnergy];
endfunction

// ------------------------------ HYDROGEN BONDS -----------------------------

// The hydrogen bond parameter table is designed to approximate the
// contact statistics distributions for the various atom types.  The
// distance distribution is a 12-10 potential, the udir angle is gaussian
// and the out of plane is exponential.

// !!! SULFUR =S !!!!

const HBPARAM = tr [
//			a0	aisd	zisd	r0    risd
//			--      ----    ----    ----  --------
    [ '[OX1]=*',	60,	1/10,	1/35,	2.72, inv 0.25, 2.0	]
,   [ '[OX1]*=[OX1]',	60,	1/10,	1/35,	2.72, inv 0.25, 2.0	]
,   [ '[OQ1][i]',	63,	1/ 7,	0,52.65, inv 0.25, 2.5	]

,   [ '[OQ1]',		65,	1/12,	0,	2.70, inv 0.25, 2.5	]
,   [ '[OQ2]',		30,	1/35,	0,	2.72, inv 0.30, 1.4	]
,   [ '[OQ0]',		 0,	   0,	0,	2.70, inv 0.25, 2.1	]

,   [ '[NQ1X2]#*',	 0,	1/10,	0,	2.83, inv 0.25, 1.6	]
,   [ '[NQ1X1]#*',	60,	1/15,	0,	2.83, inv 0.25, 1.6	]
,   [ '[NQ1X1]=*=*',	60,	1/15,	0,	2.83, inv 0.25, 1.6	]

,   [ '[N;iQ1]',	60,	1/13,	1/20,	2.83, inv 0.25, 1.6	]
,   [ '[N;Q1][i]',	60,	1/13,	1/20,	2.83, inv 0.25, 1.6	]

,   [ '[#7;iQ2]',	 0,	1/15,	1/15,	2.80, inv 0.20, 2.0	]
,   [ '[NQ2][i]',	 0,	1/15,	1/15,	2.80, inv 0.20, 2.0	]

,   [ '[N!iQ3]',	 0,	1/25,	0,	2.75, inv 0.22, 1.8	]
,   [ '[N!iQ2]',	55,	1/25,	0,	2.75, inv 0.22, 1.8	]
,   [ '[N!iQ1]',	72,	1/14,	0,	2.75, inv 0.22, 1.8	]

,   [ '[#M]',		 0,	   0,	0,	2.70, inv 0.25, 2.5	]
];

// HBondScore scores the hydrogen bonds given the contact list:
//
//	rpos	: atom positions of receptor
//	lQpos	: heavy atom positions of ligand
//	xL	: index of contact in lQpos
//	xR	: index of contact in receptor (all atom)
//	r	: distance of contact

local function HBondScore [rec, lig, rpos, lQpos, xL, xR, r, scale]
    local mask, idx;

    const HCUTOFF = 4.5;

    // examine only the H-bond candidates within the cutoff region
    // and only those that possibly are a don-acc pair

    [xL, xR, r, scale] = [xL, xR, r, scale] || [r < HCUTOFF];
    [xL, xR, r, scale] = [xL, xR, r, scale] || [
    bitand [*rec.LdG_hbcode[xR], *lig.LdG_Qhbcode[xL]]
    ];
    if not l_length [xL, xR, r] then return [0,0]; endif

    // set hvec to the R -> L vector between the contact atoms
    // then get the parameters from the HBPARAM table

    local hvec = (apt get [lQpos, [xL]] - apt get [rpos, [xR]]);
    local l_px = *lig.LdG_Qhbpx[xL], r_px = *rec.LdG_hbpx[xR];

    // score the distance with a gaussian based upon the mean
    // separation distance between the atoms

    local r0 = 0.5 * (HBPARAM(5)[r_px] + HBPARAM(5)[l_px]);
    local s0 = 0.5 * (HBPARAM(7)[r_px] + HBPARAM(7)[l_px]);	// !!! 7

    // 80 radius maximum ? !!!
    local x = r0 / maxE[r0, r], x3 = cube x, x9 = cube x3;
    local logS = ( -1/(KBOLTZ*300) * s0 * (x9 * (5*x3 - 6*x) + 1) );

    // score the udir direction with respect to the receptor
    // u direction (opposite to bonds) divided by two (geometric mean)
    // geometric mean'd with the ligand u direction score

    local ang0 = HBPARAM(2)[r_px], aisd = HBPARAM(3)[r_px];
    local ang = (180/PI) * acos maxE [-1, minE [1,
    add (hvec * apt get [*rec.LdG_udir, [xR]]) * invz r
    ]];
    logS = logS + 0.5 * ( -0.5 * sqr ( (ang - ang0) * aisd ) ); 

    // !!! MAKE FASTER / DON'T USE HEAVY MOL
    local lmol = *lig.LdG_Qmol;
    local l_udir = GetUdir [lmol, lQpos, xL];

    ang0 = HBPARAM(2)[l_px]; aisd = HBPARAM(3)[l_px];
    ang = (180/PI) * acos maxE [-1, minE [1, - add (hvec * l_udir) * invz r ]];
    logS = logS + 0.5 * ( -0.5 * sqr ( (ang - ang0) * aisd ) ); 

    // perform the out of plane test with respect to the
    // receptor and ligand (geometric mean)

    local oop = (180/PI) * asin abs maxE [-1, minE [1,
    add (hvec * apt get [*rec.LdG_zdir, [xR]]) * invz r
    ]];
    local itheta = HBPARAM(4)[r_px];
    logS = logS - 0.5 * (oop * itheta);

    local l_zdir = GetZdir [lmol, lQpos, xL];
    oop = (180/PI) * asin abs maxE[-1, minE[1, add (hvec * l_zdir) * invz r]];
    itheta = HBPARAM(4)[l_px];
    logS = logS - 0.5 * (oop * itheta);

#if 0
    return [add (scale * exp logS), 0];
#else
    local hbAngle = 0;

    // h-bond network score
    mask = logS > log 0.1;
    [xL,xR] = [xL,xR] || [mask];

    hbAngle = add log mtoc (rotrpoke [xL, 0] <> xL);

#if 0
    xL = uniq xL;
    xR = uniq xR;

    local mligpos = apt get [lQpos, [xL]];
    local mpos = apt get [mol_aPos *rec.mol, [xR]];

print *rec.mol(4)(MOL_ATOM_EL)[xR];
print *lig.mol(4)(MOL_ATOM_EL)[*lig.Qindex][xL];

    local cpos = apt cat [mpos, mligpos];
    local prox = prox_open [3.0, cpos, 3.0];
    local [deg, xbond, r2] = prox_find [prox, cpos, 0];
    prox_close prox;

    xbond = split [xbond, deg];
    xbond = xbond || (xbond <> x_id xbond);

    local clist = cat graph_scycle_list xbond;
    clist = clist | app length clist < 7;

print ['cycles', clist];

//print tr [MLIG_PARAM(1)[mpx], ang, ang0, aisd];
#endif


    return [add (scale * exp logS), hbAngle];
#endif
endfunction

// ------------------------------- DISPERSION ---------------------------------

const MIN_SUMI = inv cube 50;		// min value for sumI

local function EffectiveRadius sumI = (
    cbrt inv maxE [MIN_SUMI, sumI]
);

const HAND = 1;

const XP_PARAM = tr [
    [ '?',	'[#1]',				0	],

#if 1
    [ 'NO2',	'[N+]([OX1])=[OX1]',		0.8	],	// nitro
    [ 'NO2',	'[OX1][N+]=[OX1]',		0.8	],

//    [ 'CO2',	'[CX3]([OX1])=[OX1]',		0.8	],	// carboxylate
//    [ 'CO2',	'[OX1][CX3]=[OX1]',		0.8	],

    [ 'CH3',	'[CX4H3]',			0.5	],
    [ 'CH2',	'[CX4H2]',			0.5	],
    [ 'CH1',	'[CX4H1]',			0.5	],
    [ 'Car',	'[c]',				0.5	],
    [ 'C',	'[#6]',				0.5	],

    [ 'N',	'[#7]',				0.5	],
    [ 'O',	'[#8]',				0.5	],
    [ 'P',	'[#15]',			0.5	],
    [ 'S',	'[#16]',			0.5	],
    [ 'F',	'[#9]',				0.5	],
    [ 'Cl',	'[#17]',			0.5	],
    [ 'BrI',	'[#35]',			0.5	],
    [ 'BrI',	'[#53]',			0.5	],
    [ 'M',	'[#M]',			       -2.5	],

#elseif HAND
    [ 'NO2',	'[N+]([OX1])=[OX1]',		0.8	],	// nitro
    [ 'NO2',	'[OX1][N+]=[OX1]',		0.8	],

    [ 'Chyd',	'[CQ4]',			1.0	],
    [ 'Chyd',	'[#6Q3;i]',			1.0	],
    [ 'Chyd',	'[#6Q2]#*',			1.0	],
    [ 'Chyd',	'[#6Q2](=*)=*',			1.0	],

    [ 'CX4',	'[CX4]C=O',		        0.5	],
    [ 'CX4',	'[CX4][N+]',		        0.5	],
    [ 'CX4',	'[CX4][S+*][OX1]',	        0.5	],
    [ 'CX4',	'[CX4][P+][OX1]',	        0.5	],
    [ 'CX4',	'[CX4][NX3!i][C+1]',	        0.5	],
    [ 'CX4',	'[CX4][NX3!i][C+0]=[N+]',       0.5	],

    [ 'Chyd',	'[CX4]',			1.0	],
    [ 'C',	'[#6]',				0.5	],	// pi carbon

    [ 'Nb',	'[#7;iQ3]',		       -0.5	],
    [ 'Nb',	'[N;Q3!i][i]',		       -0.5	],
    [ 'Nb',	'[N;Q2](=*)=*',		       -0.5	],
    [ 'Nb',	'[N;Q2]#*',		       -0.5	],
    [ 'NH',	'[n;X2]:a:[n;i2H]',	       -1.3	],
    [ 'NH',	'[n;X2]:[n;i2H]',	       -1.3	],
    [ 'NH',	'[n;i2H]:a:[n;X2]',	       -1.3	],
    [ 'NH',	'[n;i2H]:a:[n;i1H]',	       -1.3	],
    [ 'NH',	'[#7!H0]',		       -1.8	],
    [ 'N',	'[#7]',		      	       -1.0	],

//    [ 'O=S',	'[OX1][S+*]',		        0.5	],
//    [ 'O=P',	'[OX1][P+][OX1]',	        0.5	],
    [ 'Oar',	'o',				0.6	],	// furan O
    [ 'O=',	'[OX1]',		       -1.0	],
    [ 'OH',	'[O!H0]',		       -1.3	],
    [ 'O',	'[#8]',			       -1.3 	],

    [ 'Sar',	's',				0.6	],
    [ 'SO2',	'[S+*][OX1]',		        1.5	],
    [ 'SH',	'[SX2!H0]',		       -0.5	],
    [ 'S',	'[#16]',			1.5	],

    [ 'P',	'[#15]',		       -1.0	],

    [ 'BrI',  	'[#53]',		        1.5	],	// halogens
    [ 'BrI',  	'[#35]',		        1.5	],
    [ 'Cl',  	'[#17]',			0.3	],
    [ 'F',  	'[#9]',				0.7	],

//    [ 'TM',	'[#T]',			       -3.0	],	// metals
    [ 'M',	'[#M]',			       -2.5	],
#else
    [ 'CX4',	'[CX4]',			0.6861	],
    [ 'C',	'[#6]',				0.4205	],	// pi carbon

    [ 'Nar',	'[n]',			       -0.1078	],
    [ 'Nar',	'[NQ3X3!i]C=[#G6]',	       -0.1078	],
    [ 'Nar',	'[NQ3X3!i]([i])[i]',	       -0.1078	],

    [ 'N=',	'[N+!$(*[-])]=[C+0][N;X3!i]',   0.4183	],	// !!! N= sign
    [ 'N=',	'[N;X3!i][C+0]=[N+!$(*[-])]',   0.4183	],
    [ 'N=',	'[N;X3!i][C+1][NX3!i]',	        0.4183	],
    [ 'N=',	'[N;X3!i][c+0]:[n+!$(*[-])]',   0.4183	],
    [ 'N=',	'[#7;i]',		        0.4183	],
    [ 'N=',	'[N;i]',		        0.4183	],
    [ 'N=',	'[NQ3!i!X4][i]',	        0.4183	],
    [ 'Ni',	'[N][i]',		       -1.4803	],
    [ 'N',	'[#7]',			       -1.3563	],

    [ 'O2',	'[OX1]=C[OX1]',		       -1.0147	],
    [ 'O2',	'[OX1]C=[OX1]',		       -1.0147	],
    [ 'O2',	'[OX1][SX4]([OX1])[OX1]',      -1.0147	],
    [ 'O2',	'[OX1][PX4][OX1]',	       -1.0147	],
    [ 'O2',	'[O-!$(*[+*])]',	       -1.0147	],
//    [ 'O=',	'[OX1]',		        0.8448	],	// !!! O= sign
    [ 'O=',	'[OX1]',		       -1.0147	],	// !!! O= sign
    [ 'OH',	'[#8;!H0]',		       -1.3421	],
    [ 'O',	'[#8]',			       -0.9966	],

    [ 'P',	'[#15]',			3.3223	],
    [ 'S',	'[#16]',			1.7216	],
    [ 'X',  	'[#G7]',			0.6774	],	// halogens
    [ 'M',	'[#T]',				5.7941	],	// metals
#endif

    [ '?',	'*',				0	]	// default
];

// IntegrationRadius returns the effective integration radius of an atom
// by taking into accout its bonded neighbors.  We require
//
//	solR		the atom radius
//	catxbond	the cat'd bond graph
//	bondlen		the cat'd lengths of the bonds
//	seg		the degree of each atom

local function IntegrationRadius [solR, catxbond, bondlen, seg]
    local r1 = stretch [solR, seg], r2 = solR[catxbond];
    local d = maxE [abs (r1 - r2), minE [bondlen, (r1 + r2)]];

    local h1 = (sqr r2 - sqr (r1 - d)) / (d + d);
    local h2 = (sqr r1 - sqr (r2 - d)) / (d + d);

    local ov = s_add [ 0.5 * (sqr h1 * (3*r1-h1) + sqr h2 * (3*r2-h2)), seg ];
    local V = (PI/3) * maxE [0, 4 * cube solR - ov];

    return 0.95 * cbrt (V * (3/(4*PI)));
endfunction

// Heavy_solS calculates the heavy atom integration radius for a given
// collection of atoms

local function Heavy_solS [mol, solR, Qmask]
    local n = l_length append [mol(4),solR];
    local pos = mol_aPos mol;

    solR = resize [solR, n];
    mol(4) = apt resize [mol(4), n];
    pos = apt resize [pos, n];

    local xbond = graph_mget [mol(4)(MOL_ATOM_BONDS), Qmask];
    solR = solR | Qmask;
    pos = pos || [Qmask];

    local deg = app length xbond;
    local xA = stretch [x_id xbond, deg], xB = cat xbond;
    local blen = norm (apt get [pos, [xA]] - apt get [pos, [xB]]);
    local solS = IntegrationRadius [ solR, xB, blen, deg ];

    return mput [rep [0, n], Qmask, solS];
endfunction

// xAromaticRings returns the indices of the aromatic rings of a
// collection of atoms

local function xAromaticRings atoms
    local rings = cat [
    sm_MatchAtoms [ 'a1:a:a:a:a:a:1', atoms ],
    sm_MatchAtoms [ 'a1:a:a:a:a:1', atoms ]
    ];
    rings = rings | app length rings;
    rings = rings | m_uniq app sort rings;
    rings = split [ indexof [cat rings, atoms], app length rings ];
    rings = rings | app andE rings;
    return rings;
endfunction

// closeReceptor and openReceptor handle the scoring function operators

local function closeReceptor [rec, opt]
    local t;

    if isfalse *rec.LdG_init then return; endif
    if (*rec.LdG_init = *rec.LdG_init - 1) then return; endif

    prox_close *rec.LdG_Mprox;		// metal data
    *rec.LdG_Mrad = [];
    *rec.LdG_Mxmet = [];

    prox_close *rec.LdG_prox;
    prox_close *rec.LdG_qprox;
    prox_close *rec.LdG_aroprox;
    prox_close *rec.LdG_alkprox;

    for t in uniq (tags *rec | m_findmatch ['LdG_*', tags *rec]) loop
    *rec.(t) = [];
    endloop
endfunction

local function openReceptor [rec, opt]
    if add *rec.LdG_init then			// already initialized?
    *rec.LdG_init = *rec.LdG_init + 1;
    return;
    endif

    local mol = *rec.mol, pos = mol_aPos mol;
    local Qmask = *rec.Qmask, Qindex = *rec.Qindex;

    local el = mol(4)(MOL_ATOM_EL);

    *rec.LdG_zdir = GetZdir [mol, pos, igen l_length pos];
    *rec.LdG_udir = GetUdir [mol, pos, igen l_length pos];

    // create the terms for metal ligation.  We locate the metal
    // atoms in the receptor and set up a prox for fast searching;
    // we use a prox to locate each atom's receptor ligators

    local xmetal = x_pack (el_Metal el and el_spValence el <= 3);

    *rec.LdG_Mprox = prox_open [MCUTOFF, apt get [pos, [xmetal]], MCUTOFF];
    *rec.LdG_Mrad  = MetalLigandRadius el[xmetal];
    *rec.LdG_Mxmet = xmetal;

    local prox = prox_open [MCUTOFF, apt get [pos, [Qindex]], MCUTOFF];
    local [seg, idx] = prox_find [prox, apt get [pos,[xmetal]], 0];
    prox_close prox;

    *rec.LdG_Mxnbr = apt diff [split [Qindex[idx], seg], xmetal];

    // !!!! GET THIS OUT OF HERE AND INTO DOCK_U (RECEPTOR STD)

    local chains = [];
    local atoms = *rec.__atoms;

    local old_prio = task_prio 0;

    local flag = length atoms <> mol_aCount mol;
    if flag then
    local snew = SystemOpen [];
    local scurr = SystemCurrent snew;

    chains = mol_Create mol;
    atoms = cat cAtoms chains;
    endif

    local solR = SolvationRadius atoms;

    *rec.LdG_type = sm_Indexof [atoms, XP_PARAM(2)];		// !!!
    *rec.LdG_hbpx = sm_Indexof [atoms, HBPARAM(1)];		// !!!

    local rq = FormalCharge atoms;

    // InitAroAlkTerm initializes a prox used to find the
    // aromatic ring-center alkane term

    function InitAroAlkTerm []
    local xring = xAromaticRings atoms;
    local rsize = app length xring;
    local center = (
        apt s_add [aPos atoms[cat xring], [rsize]] * [invz rsize]
    );
    *rec.LdG_aroprox = prox_open [AROCUTOFF, center, AROCUTOFF];

    local alkpos = aPos (atoms | sm_Match ['[CX4!H0]', atoms]);
    *rec.LdG_alkprox = prox_open [AROCUTOFF, alkpos, AROCUTOFF];
    endfunction

    InitAroAlkTerm [];

    if flag then
    oDestroy chains;

    SystemCurrent scurr;
    SystemClose snew;
    endif

    task_prio old_prio;
    // !!!

    // create a prox that contains the heavy atoms of the
    // receptor (since we ignore light atoms); then create
    // a prox for calculating the ionic term

    local rpos = apt get [pos, [Qindex]];
    *rec.LdG_prox = prox_open [CUTOFF, rpos, CUTOFF];

    local qpos = pos || [rq <> 0];
    rq = rq | rq <> 0;

    *rec.LdG_qprox = prox_open [CUTOFF, qpos, CUTOFF];
    *rec.LdG_q = rq;

    // compute the solvation radius and the integration radius
    // and initialize the born radius accumulators to zero (they
    // will be filled in on demand based on ligand positions)

    *rec.LdG_solS = Heavy_solS [mol, solR, Qmask];
    *rec.LdG_solR = solR;

    *rec.LdG_effR = zero Qmask;		// uncomplexed born radius
    *rec.LdG_sumI = zero Qmask;		// uncomplexed born integration

    local rmetal = not *rec.m_tmetal and el_Metal *rec.mol(4)(MOL_ATOM_EL);
    local hbcode = bitor [
    select [0x1, 0, *rec.m_don or rmetal],
    select [0x2, 0, *rec.m_acc]
    ];
    *rec.LdG_hbcode = mput [hbcode, not *rec.LdG_hbpx, 0];

    *rec.LdG_init = 1;
endfunction

// closeLigand and openLigand handle the scoring function operators

local function closeLigand [lig, opt]
    local t;
    for t in uniq (tags *lig | m_findmatch ['LdG_*', tags *lig]) loop
    *lig.(t) = [];
    endloop
endfunction

local function openLigand [lig, opt]
    local mol = *lig.mol, el = mol(4)(MOL_ATOM_EL), m;
    local Qmask = *lig.Qmask, Qindex = *lig.Qindex;

    // "Turn off" unimportant atoms: remove them from Qmask and Qindex

    if length *lig.ligData.useAtoms then
    Qmask | not *lig.ligData.useAtoms = 0;

    m = indexof [Qindex, x_pack *lig.ligData.useAtoms];
    Qindex = Qindex | m;
    endif

    *lig.LdG_Qmol = mol_aMask [mol, Qmask];

    // compute the flexibility measure and atom types

    local atoms = *lig.__atoms;
    local chains = [], psys = [];

    local old_prio = task_prio 0;

    if length atoms <> mol_aCount mol then
    psys = SystemPush [];
    chains = mol_Create mol;
    atoms = cat cAtoms chains;
    endif

//  *lig.LdG_flex  = FlexibilityScore atoms;
    *lig.LdG_flex  = FlexibilityScore (atoms | Qmask); // remove unimportant
    *lig.LdG_Qtype = sm_Indexof [atoms[Qindex], XP_PARAM(2)];
    *lig.LdG_Qhbpx = sm_Indexof [atoms[Qindex], HBPARAM(1)];

    local solR = SolvationRadius atoms;

    // Mrad(i) is the metal ligation radius for the i'th atom
    // Mxml are the indices of the metal ligators in the molecule

    *lig.LdG_Mrad = MetalLigandRadius el;
    *lig.LdG_Mxml = x_pack(ph4_aMetalLigand atoms and aAtomicNumber atoms > 1);
    *lig.LdG_Mpx  = sm_Indexof [atoms, MLIG_PARAM(1)];

    *lig.LdG_q = FormalCharge atoms;
    *lig.LdG_xq = x_pack (*lig.LdG_q <> 0);

    *lig.LdG_xalk = x_pack sm_Match ['[CX4!H0]', atoms ];

    local xring = xAromaticRings atoms;
    *lig.LdG_xaroring = cat xring;
    *lig.LdG_xaroring_size = app length xring;

    if length psys then SystemPop psys; endif
    task_prio old_prio;

    // compute the solvation radius and the integration radius
    // by converting the solvation radii using the bonded neighbors

    *lig.LdG_QsolS = (Heavy_solS [mol, solR, Qmask])[Qindex];
    *lig.LdG_QsolR = solR[Qindex];

    local hbcode = bitor [
    select [0x2, 0, *lig.m_don[Qindex]],	// opposite of receptor
    select [0x1, 0, *lig.m_acc[Qindex]]
    ];
    *lig.LdG_Qhbcode = mput [hbcode, not *lig.LdG_Qhbpx, 0];
endfunction

// this is the main entry point for the scoring function in the docking system

local function dock_score_London_dG [cmd, arg, opt]
    const DEFAULTS = [
    train:	0
    ];

    if cmd === 'ID' then
    return 'London dG';

    elseif cmd === 'configpanelwidgets' then
    return [];

    elseif cmd === 'configpanelevent' then	// arg = [val,trig]
    return 0;

    elseif cmd === 'configvalues' then		// arg = val
    return [];

    elseif cmd === 'openReceptor' then		// arg = rec
    openReceptor [arg, opt];
    return;

    elseif cmd === 'closeReceptor' then		// arg = rec
    closeReceptor [arg, opt];
    return;

    elseif cmd === 'openLigand' then		// arg = lig
    openLigand [arg, opt];
    return;

    elseif cmd === 'closeLigand' then		// arg = lig
    closeLigand [arg, opt];
    return;

    elseif not (cmd === 'score') then		// arg = [ligpos,lig,rec]
    return;
    endif

    // ExposureDelta takes integration summations in the complex
    // and isolated and returns the delta

    local function ExposureDelta [solR, sumIc, sumI] = (
    cube solR * (maxE [MIN_SUMI, sumIc] - maxE [MIN_SUMI, sumI])
    );

    local i, seg, idx, r2, mask;

    // perform the scoring of the ligpos using the lig,rec,opt

    local [ligpos, lig, rec] = arg;

    local rpos = mol_aPos *rec.mol;		// receptor data
    local rQindex = *rec.Qindex;		// heavy atom indices
    local r_solR = *rec.LdG_solR;
    local r_solS = *rec.LdG_solS;
    local r_sumI = *rec.LdG_sumI;
    local r_effR = *rec.LdG_effR;

//  local lQmask  = *lig.Qmask;			// mask for heavy atoms (unused)
    local lQindex = *lig.Qindex;		// indices of heavy atoms
    if length *lig.ligData.useAtoms then
    local m   = indexof [*lig.Qindex, x_pack *lig.ligData.useAtoms];
    lQindex   = lQindex | m;
    endif
    local lQpos   = apt get [ligpos, [lQindex]];

    local l_QsolR = *lig.LdG_QsolR;
    local l_QsolS = *lig.LdG_QsolS;

    // determine the interactions between the ligand the receptor

    local [xRseg, xR, xD] = prox_find [*rec.LdG_prox, lQpos, 0];
    local xL = stretch [x_id lQindex, xRseg];
    xR = rQindex[xR];
    xD = sqrt xD;

    local [metScore,metAngle] = MetalContactScore [ligpos, lig, rec, opt];
    local [hbScore,hbAngle] = HBondScore [rec,lig,rpos,lQpos,xL,xR,xD,1.0];

    // determine the ionic interactions

    function IonContactScore [rec, lig, ligpos]
//	return [0,0];
    local rq = *rec.LdG_q, lq = *lig.LdG_q;
    local lqpos = apt get [ligpos, [*lig.LdG_xq]];
    local [seg,xB,r2] = prox_find [*rec.LdG_qprox, lqpos, 0];
    local xA = stretch [*lig.LdG_xq, seg];
    r2 = maxE [r2, sqr 1.0];
    local qq = lq[xA] * rq[xB];
    local E = COULOMB_SCALE * 0.25 * (
        qq * (inv sqr (1 + sqrt r2) - inv sqr (1 + CUTOFF))
//	    qq * (inv (1 + sqrt r2) - inv (1 + CUTOFF))
    );
    return [add (E | qq > 0), add (E | qq < 0)];
    endfunction

    local [ionScoreP,ionScoreN] = IonContactScore [rec, lig, ligpos];

    // calculate the aro-alkane term

    local ar2 = third prox_find [
    *rec.LdG_aroprox, apt get [ligpos, [*lig.LdG_xalk]], 0
    ];
    local aroalk = add sqr (1 - sqrt ar2 / AROCUTOFF);

    if length *lig.LdG_xaroring_size then
    local larocenter = apt s_add [
        apt get [ ligpos, [*lig.LdG_xaroring] ],
        [ *lig.LdG_xaroring_size ]
    ] * [ invz *lig.LdG_xaroring_size ];
    ar2 = third prox_find [ *rec.LdG_alkprox, larocenter, 0 ];
    aroalk = aroalk + add sqr (1 - sqrt ar2 / AROCUTOFF);
    endif

    // calculate a repulsive term

#if 0
    local REPSCALE = 1.0;
    local rsum = *lig.radius[lQindex][xL] + *rec.radius[xR];

    i = x_pack (xD < rsum);
    local v = maxE[ 0, minE[ 1, 1 - xD[i]/rsum[i] ]];
    local repulse = REPSCALE * add (cube v * (6 * sqr v - 15 * v + 10));
#elseif 0
    local rsum = 0.7 * (*lig.radius[lQindex][xL] + *rec.radius[xR]);
    i = x_pack (xD < rsum);

    local repulse = 5.0 * add (1 - xD[i] / rsum[i]);
#else
    const REPSCALE = 70.0;
    local rsum = 0.9 * (*lig.radius[lQindex][xL] + *rec.radius[xR]);
    i = x_pack (xD < rsum);

    local repulse = REPSCALE * add cube (1 - maxE [0.8, xD[i] / rsum[i]]);
#endif

    // locate the interacting atoms between the receptor and the
    // ligand (convert receptor contacts to atom indices)
    //
    // if the receptor atoms have not been initialized (for the
    // born calculation) then do so now; this avoids calculating
    // the born stuff for the whole receptor in the beginning
    //
    // r_effR = uncomplexed receptor atom born radius
    // r_sumI = uncomplexed receptor initial integral sums

    if orE (mask = r_effR[xR] <= 0) then
    for i in uniq (xR | mask) loop
        [[], idx, r2] = prox_find [ *rec.LdG_prox, apt peek [rpos, i], 0 ];
        [idx, r2] = [idx, r2] || [r2 > sqr 0.5];
        r_sumI(i) = cube inv r_solR(i) - add first _pot_BornIntegral [
        sqrt r2, r_solR(i), r_solS[rQindex[idx]]
        ];
    endloop

    r_effR = EffectiveRadius r_sumI;
    *rec.LdG_sumI = r_sumI;
    *rec.LdG_effR = r_effR;
    endif

    // compute the ligand uncomplexed data and the interaction data
    // l_effR = uncomplexed heavy atom born radius
    // l_sumI = uncomplexed initial integral sums

    local l_sumI = zero l_QsolR;
    local lprox = prox_open [CUTOFF, lQpos, CUTOFF];

    for i = 1, length lQindex loop
    [[], idx, r2] = prox_find [ lprox, apt peek [lQpos, i], 0 ];
    [idx, r2] = [idx, r2] || [r2 > sqr 0.5];
    l_sumI(i) = cube inv l_QsolR(i) - add first _pot_BornIntegral [
        sqrt r2, l_QsolR(i), l_QsolS[idx]
    ];
    endloop

    prox_close lprox;

    // compute the born radius of the complexed ligand by adding in
    // the interaction summations from the receptor
    // l_sumI_c = complexed ligand heavy integral sums

    local l_sumI_c = l_sumI - s_add [
    first _pot_BornIntegral [xD, l_QsolR[xL], r_solS[xR]],
    xRseg
    ];
    local l_delta = ExposureDelta [ l_QsolR, l_sumI_c, l_sumI ];
    local l_type = *lig.LdG_Qtype;

    // compute the bron radius of the complexed receptor by adding in
    // the interaction summations from the ligand
    // r_sumI_c = complexed receptor integral sums

    [idx, mask] = sam xR;
    [xL,xR,xD] = apt get [[xL,xR,xD], [idx]];

    local r_sumI_c = r_sumI[xR|mask] - s_add [
    first _pot_BornIntegral [xD, r_solR[xR], l_QsolS[xL]],
    mtoc mask
    ];

    xR = xR | mask;
    local r_delta = ExposureDelta [ r_solR[xR], r_sumI_c, r_sumI[xR] ];
    local r_type = *rec.LdG_type[xR];

    // if we are in training mode form the tagged vector of descriptors
    // that will be used to write to a training database (for QSAR)

    if istrue opt.train then
    local delta = cat [r_delta, l_delta];
    local dtype = XP_PARAM(1)[cat [r_type, l_type]];

    [idx,mask] = sam dtype;			// merge results
    delta = s_add [delta[idx], mtoc mask];
    dtype = dtype[idx] | mask;

    local utype = sort uniq XP_PARAM(1);
    local udata = put [zero utype, indexof [dtype, utype], delta];

    [utype,udata] = apt cat [
        [utype, udata],
        ['hbond', hbScore],
        ['hbondA', hbAngle],
        ['metal', metScore],
        ['metal_r', metAngle],
        ['ionP', ionScoreP],
        ['ionN', ionScoreN],
        ['aroalk', aroalk],
        ['repulse', repulse / REPSCALE],
        ['flex', *lig.LdG_flex]
    ];
    utype = tok_cat ['LdG_', utype];
    return tag [utype,udata];
    endif

    local score = (
    - 3.16
    + 0.2113 * *lig.LdG_flex
#if HAND
    - 2.0 * hbScore
    - 2.5 * metScore
    + 3.5 * metAngle
    + 0.0 * hbAngle
    - 1.0 * aroalk
#else
    - 1.4382 * hbScore
    - 3.2190 * metScore
    + 2.5065 * metAngle
#endif
    + 0.0000 * ionScoreP
    + 0.5000 * ionScoreN
    + add ((cat XP_PARAM(3))[r_type] * r_delta)
    + add ((cat XP_PARAM(3))[l_type] * l_delta)
    );

    //return score + repulse;
    return [score + repulse, hbScore, hbAngle, metScore, metAngle, ionScoreP, ionScoreN, aroalk, repulse, *lig.LdG_flex, add ((cat XP_PARAM(3))[r_type] * r_delta), add ((cat XP_PARAM(3))[l_type] * l_delta)];
endfunction



local function main []
    Open './complex_prep.moe';
    local c = Chains[];
    local c0 = cAtoms c;
    local l0 = length c0;
    local lig = c0(l0);
    local i = 0;
    local rec = [];
    while (i = inc i) < l0 loop
        rec = cat [rec, c0(i)];
    endloop
    //Close comen;

    local ligen = db_ReadColumn [ 'dock.mdb', 'mol' ];
    local cl = mol_Create ligen(1);	
    local cl0 = cAtoms cl;
    lig = cl0(1);	

    local rec_dock = dock_ReceptorOpen [rec, [aPos lig],[],[]];
    local lig_dock = dock_LigandOpen [lig, []];
    dock_score_London_dG ['openReceptor', rec_dock, []];
    dock_score_London_dG ['openLigand', lig_dock, []];	
    local LondondG = dock_score_London_dG ['score', [aPos lig,lig_dock,rec_dock], [ hbond : -0.66,	ionic : 1,	mlig : -1,	cHH : -0.01235,	cHP : 0.02497,	cXX : -0.00834 ]];
    dock_score_London_dG ['closeReceptor', rec_dock, []];
    dock_score_London_dG ['closeLigand', lig_dock, []];	
    dock_ReceptorClose rec_dock;
    dock_LigandClose lig_dock;
    //print LondondG;
    fwrite [ 'log.txt', '{},{},{},{},{},{},{},{},{},{},{},{}', LondondG[1], LondondG[2], LondondG[3], LondondG[4], LondondG[5], LondondG[6], LondondG[7], LondondG[8], LondondG[9], LondondG[10], LondondG[11], LondondG[12]]; 
endfunction
'''
        self.gold_config = '''  GOLD CONFIGURATION FILE

  AUTOMATIC SETTINGS
autoscale = 1

  POPULATION
popsiz = auto
select_pressure = auto
n_islands = auto
maxops = auto
niche_siz = auto

  GENETIC OPERATORS
pt_crosswt = auto
allele_mutatewt = auto
migratewt = auto

  FLOOD FILL
radius = 12
origin = {} {} {}
do_cavity = 0
floodfill_atom_no = 0
cavity_file = 
floodfill_center = point

  DATA FILES
ligand_data_file {} 10
param_file = DEFAULT
set_ligand_atom_types = 1
set_protein_atom_types = 0
directory = .
tordist_file = DEFAULT
make_subdirs = 0
save_lone_pairs = 1
fit_points_file = fit_pts.mol2      
read_fitpts = 0

  FLAGS
internal_ligand_h_bonds = 0
flip_free_corners = 0
match_ring_templates = 0
flip_amide_bonds = 0
flip_planar_n = 1 flip_ring_NRR flip_ring_NHR
flip_pyramidal_n = 0
rotate_carboxylic_oh = flip
use_tordist = 1
postprocess_bonds = 1
rotatable_bond_override_file = DEFAULT
solvate_all = 1

  TERMINATION
early_termination = 1
n_top_solutions = 3
rms_tolerance = 1.5

  CONSTRAINTS
force_constraints = 0

  COVALENT BONDING
covalent = 0

  SAVE OPTIONS
save_score_in_file = 1 unweighted
save_protein_torsions = 1
concatenated_output = {}    
clean_up_option delete_all_solutions
clean_up_option delete_redundant_log_files
clean_up_option delete_all_initialised_ligands
clean_up_option delete_empty_directories
clean_up_option delete_rank_file
clean_up_option delete_all_log_files
output_file_format = MACCS

  FITNESS FUNCTION SETTINGS
initial_virtual_pt_match_max = 3
relative_ligand_energy = 1
gold_fitfunc_path = {}
start_vdw_linear_cutoff = 6
score_param_file = DEFAULT

  RUN TYPE
run_flag = RESCORE retrieve

  PROTEIN DATA
protein_datafile = {}  
'''
        self.dock_contact = '''conformer_search_type 	rigid
use_internal_energy	yes
internal_energy_rep_exp	12
internal_energy_cutoff	100.0
ligand_atom_file	{}
limit_max_ligands	no
skip_molecule	no
read_mol_solvation	no
calculate_rmsd	no
use_database_filter	no
orient_ligand	no
bump_filter	no
score_molecules	yes
contact_score_primary	yes
contact_score_secondary	no 
contact_score_cutoff_distance	4.5
contact_score_clash_overlap	0.75
contact_score_clash_penalty	50
contact_score_grid_prefix	grid
grid_score_primary	no
grid_score_secondary	no
multigrid_score_secondary	no
dock3.5_score_secondary	no
continuous_score_secondary	no
footprint_similarity_score_secondary	no
pharmacophore_score_secondary	no
descriptor_score_secondary	no
gbsa_zou_score_secondary	no
gbsa_hawkins_score_secondary	no
SASA_score_secondary	no
amber_score_secondary	no
minimize_ligand	no
atom_model	all
vdw_defn_file	/opt/dock/6.8/parameters/vdw_AMBER_parm99.defn
flex_defn_file	/opt/dock/6.8/parameters/flex.defn
flex_drive_file	/opt/dock/6.8/parameters/flex_drive.tbl
chem_defn_file /opt/dock/6.8/parameters/chem.defn
pharmacophore_defn_file /opt/dock/6.8/parameters/ph4.defn
ligand_outfile_prefix	{}
write_orientations	no
num_scored_conformers	1
write_conformations	no
rank_ligands	no
'''
        self.dock_continuous = '''conformer_search_type 	rigid
use_internal_energy	yes
internal_energy_rep_exp	12
internal_energy_cutoff	100.0
ligand_atom_file	{}
limit_max_ligands	no
skip_molecule	no
read_mol_solvation	no
calculate_rmsd	no
use_database_filter	no
orient_ligand	no
bump_filter	no
score_molecules	yes
contact_score_primary	no
contact_score_secondary	no
grid_score_primary	no
grid_score_secondary	no
multigrid_score_primary	no
multigrid_score_secondary	no
dock3.5_score_primary	no
dock3.5_score_secondary	no
continuous_score_primary	yes
continuous_score_secondary	no
cont_score_rec_filename	./protein_99sb.mol2
cont_score_att_exp	6
cont_score_rep_exp	12
cont_score_rep_rad_scale	1.0
cont_score_use_dist_dep_dielectric	yes
cont_score_dielectric	4.0
cont_score_vdw_scale	1.0
cont_score_es_scale		1.0
footprint_similarity_score_secondary	no
pharmacophore_score_secondary	no
descriptor_score_secondary	no
gbsa_zou_score_secondary	no
gbsa_hawkins_score_secondary	no
SASA_score_secondary	no
amber_score_secondary	no
minimize_ligand	no
atom_model	all
vdw_defn_file	/opt/dock/6.8/parameters/vdw_AMBER_parm99.defn
flex_defn_file	/opt/dock/6.8/parameters/flex.defn
flex_drive_file	/opt/dock/6.8/parameters/flex_drive.tbl
chem_defn_file /opt/dock/6.8/parameters/chem.defn
pharmacophore_defn_file /opt/dock/6.8/parameters/ph4.defn
ligand_outfile_prefix	{}
write_orientations	no
num_scored_conformers	1
write_conformations	no
rank_ligands	no
'''
        self.dock_grid = '''conformer_search_type 	rigid
use_internal_energy	yes
internal_energy_rep_exp	12
internal_energy_cutoff	100.0
ligand_atom_file	{}
limit_max_ligands	no
skip_molecule	no
read_mol_solvation	no
calculate_rmsd	no
use_database_filter	no
orient_ligand	no
bump_filter	no
score_molecules	yes
contact_score_primary	no
contact_score_secondary	no
grid_score_primary	yes
grid_score_secondary	no
grid_score_rep_rad_scale 	1
grid_score_vdw_scale	1
grid_score_es_scale	1
grid_score_grid_prefix 	grid
multigrid_score_secondary	no
dock3.5_score_secondary	no
continuous_score_secondary	no
footprint_similarity_score_secondary	no
pharmacophore_score_secondary	no
descriptor_score_secondary	no
gbsa_zou_score_secondary	no
gbsa_hawkins_score_secondary	no
SASA_score_secondary	no
amber_score_secondary	no
minimize_ligand	no
atom_model	all
vdw_defn_file	/opt/dock/6.8/parameters/vdw_AMBER_parm99.defn
flex_defn_file	/opt/dock/6.8/parameters/flex.defn
flex_drive_file	/opt/dock/6.8/parameters/flex_drive.tbl
chem_defn_file /opt/dock/6.8/parameters/chem.defn
pharmacophore_defn_file /opt/dock/6.8/parameters/ph4.defn
ligand_outfile_prefix	{}
write_orientations	no
num_scored_conformers	1
write_conformations	no
rank_ligands	no
'''
        self.dock_hawkins = '''conformer_search_type 	rigid
use_internal_energy	yes
internal_energy_rep_exp	12
internal_energy_cutoff	100.0
ligand_atom_file	{}
limit_max_ligands	no
skip_molecule	no
read_mol_solvation	no
calculate_rmsd	no
use_database_filter	no
orient_ligand	no
bump_filter	no
score_molecules	yes
contact_score_primary	no
contact_score_secondary	no
grid_score_primary	no
grid_score_secondary	no
multigrid_score_primary	no
multigrid_score_secondary	no
dock3.5_score_primary	no
dock3.5_score_secondary	no
continuous_score_primary	no
continuous_score_secondary	no
footprint_similarity_score_primary	no
footprint_similarity_score_secondary	no
pharmacophore_score_primary	no
pharmacophore_score_secondary	no
descriptor_score_primary	no
descriptor_score_secondary	no
gbsa_zou_score_primary	no
gbsa_zou_score_secondary	no
gbsa_hawkins_score_primary	yes
gbsa_hawkins_score_secondary	no
gbsa_hawkins_score_rec_filename	./protein_99sb.mol2
gbsa_hawkins_score_solvent_dielectric	78.5
gbsa_hawkins_use_salt_screen	no
gbsa_hawkins_score_gb_offset	0.09
gbsa_hawkins_score_cont_vdw_and_es	no
gbsa_hawkins_score_vdw_att_exp	6
gbsa_hawkins_score_vdw_rep_exp	12
grid_score_rep_rad_scale	1.0
gbsa_hawkins_score_grid_prefix	grid
SASA_score_secondary	no
amber_score_secondary	no
minimize_ligand	no
atom_model	all
vdw_defn_file	/opt/dock/6.8/parameters/vdw_AMBER_parm99.defn
flex_defn_file	/opt/dock/6.8/parameters/flex.defn
flex_drive_file	/opt/dock/6.8/parameters/flex_drive.tbl
chem_defn_file /opt/dock/6.8/parameters/chem.defn
pharmacophore_defn_file /opt/dock/6.8/parameters/ph4.defn
ligand_outfile_prefix	{}
write_orientations	no
num_scored_conformers	1
write_conformations	no
rank_ligands	no
'''
        self.dock_sasa = '''conformer_search_type 	rigid
use_internal_energy	yes
internal_energy_rep_exp	12
internal_energy_cutoff	100.0
ligand_atom_file	{}
limit_max_ligands	no
skip_molecule	no
read_mol_solvation	no
calculate_rmsd	no
use_database_filter	no
orient_ligand	no
bump_filter	no
score_molecules	yes
contact_score_primary	no
contact_score_secondary	no
grid_score_primary	no
grid_score_secondary	no
multigrid_score_primary	no
multigrid_score_secondary	no
dock3.5_score_primary	no
dock3.5_score_secondary	no
continuous_score_primary	no
continuous_score_secondary	no
footprint_similarity_score_primary	no
footprint_similarity_score_secondary	no
pharmacophore_score_primary	no
pharmacophore_score_secondary	no
descriptor_score_primary	no
descriptor_score_secondary	no
gbsa_zou_score_primary	no
gbsa_zou_score_secondary	no
gbsa_hawkins_score_primary	no
gbsa_hawkins_score_secondary	no
SASA_score_primary	yes
SASA_score_secondary	no
SASA_score_rec_filename	./protein_99sb.mol2
amber_score_primary	no
amber_score_secondary	no
minimize_ligand	no
atom_model	all
vdw_defn_file	/opt/dock/6.8/parameters/vdw_AMBER_parm99.defn
flex_defn_file	/opt/dock/6.8/parameters/flex.defn
flex_drive_file	/opt/dock/6.8/parameters/flex_drive.tbl
chem_defn_file /opt/dock/6.8/parameters/chem.defn
pharmacophore_defn_file /opt/dock/6.8/parameters/ph4.defn
ligand_outfile_prefix	{}
write_orientations	no
num_scored_conformers	1
write_conformations	no
rank_ligands	no
'''
        self.glide_grid = '''GRID_CENTER   {:.10f}, {:.10f}, {:.10f}
GRIDFILE   glide-grid.zip
INNERBOX   10, 10, 10
OUTERBOX   30, 30, 30
RECEP_FILE  {}
RECEP_VSCALE   1.0
RECEP_CCUT   0.25
'''
        self.glide_score = '''GRIDFILE   {}
LIGANDFILE   {}
POSES_PER_LIG   1
POSE_OUTTYPE   ligandlib
##POSTDOCK   False
NOSORT TRUE
DOCKING_METHOD   mininplace
PRECISION   {}
'''
        self.plants_score = '''# scoring function and search settings
scoring_function {}
rescore_mode simplex

# input
protein_file {}
ligand_file  {}

# output
output_dir {}

# write single mol2 files (e.g. for RMSD calculation)
write_multi_mol2 0

# binding site definition
bindingsite_center {} {} {}
bindingsite_radius 20


# cluster algorithm
cluster_structures 10
cluster_rmsd 2.0
'''
        self.rdock_score = '''RBT_PARAMETER_FILE_V1.00
TITLE rdock

RECEPTOR_FILE {}
### RECEPTOR_FLEX 3.0

##################################################################
### CAVITY DEFINITION: REFERENCE LIGAND METHOD
##################################################################
SECTION MAPPER
    SITE_MAPPER RbtLigandSiteMapper
    REF_MOL {}
	RADIUS 10.0
	SMALL_SPHERE 1.5
##	LARGE_SPHERE 4.0
	MAX_CAVITIES 1
	MIN_VOLUME 100
 	VOL_INCR 0.0
	GRIDSTEP 0.5
END_SECTION

#################################
#CAVITY RESTRAINT PENALTY
#################################
SECTION CAVITY
    SCORING_FUNCTION RbtCavityGridSF
    WEIGHT 1.0
END_SECTION	
'''
        self.xscore_score = '''
#######################################################################
#                            XTOOL/SCORE                             # 
######################################################################
###
FUNCTION	SCORE
###
### set up input and output files ------------------------------------
###
#
RECEPTOR_PDB_FILE    {}
#prepare protein
FixPDB           
#非必需
#REFERENCE_MOL2_FILE  ./1ppc_ligand.mol2
#COFACTOR_MOL2_FILE  none               
LIGAND_MOL2_FILE     {}
#prepare ligand
FixMol2
#
OUTPUT_TABLE_FILE    {}
OUTPUT_LOG_FILE      {}
###
### how many top hits to extract from the LIGAND_MOL2_FILE?
###
NUMBER_OF_HITS       1 
HITS_DIRECTORY       {}
###
### want to include atomic binding scores in the resulting Mol2 files?
###
SHOW_ATOM_BIND_SCORE	YES		[YES/NO]
###
### set up scoring functions -----------------------------------------
###
APPLY_HPSCORE         YES             	[YES/NO]
    HPSCORE_CVDW  0.004 
    HPSCORE_CHB   0.053
    HPSCORE_CHP   0.011
    HPSCORE_CRT  -0.061
    HPSCORE_C0    3.448
APPLY_HMSCORE         YES             	[YES/NO]
    HMSCORE_CVDW  0.004
    HMSCORE_CHB   0.094
    HMSCORE_CHM   0.394
    HMSCORE_CRT  -0.099
    HMSCORE_C0    3.585
APPLY_HSSCORE         YES 	  	[YES/NO]
    HSSCORE_CVDW  0.004
    HSSCORE_CHB   0.069
    HSSCORE_CHS   0.004
    HSSCORE_CRT  -0.092
    HSSCORE_C0    3.349
###
### set up chemical rules for pre-pipline ligand molecules ---------
###（类药性五规则）
APPLY_CHEMICAL_RULES    NO            [YES/NO]	
    MAXIMAL_MOLECULAR_WEIGHT      600.0
    MINIMAL_MOLECULAR_WEIGHT      200.0
    MAXIMAL_LOGP                  6.00
    MINIMAL_LOGP                  1.00
    MAXIMAL_HB_ATOM               8 
    MINIMAL_HB_ATOM               2 
###

###
#以上文件为 input_parameter_file，设置其中的输入文件及输出文件路径，选择使用的打分函数
#本项目只需按照上述参数执行即可
#XSCORE使用命令：
#xscore input_parameter_file
                    '''

    def get_xyz(self, crystal_file):
        # 根据共晶配体确定坐标
        x = os.popen(
            "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $3}' | awk '{x+=$1} END {print x/(NR-2)}'" % crystal_file).read()
        y = os.popen(
            "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $4}' | awk '{y+=$1} END {print y/(NR-2)}'" % crystal_file).read()
        z = os.popen(
            "cat %s | sed -n '/@<TRIPOS>ATOM/,/@<TRIPOS>BOND/'p | awk '{print $5}' | awk '{z+=$1} END {print z/(NR-2)}'" % crystal_file).read()
        xyz = [float(x.strip()), float(y.strip()), float(z.strip())]
        return xyz

    def get_xyz_pdbqt(self, crystal_file):
        # 根据配体来
        x = os.popen(
            "cat %s | awk '{if ($1==\"ATOM\") print $(NF-6)}' | awk '{x+=$1} END {print x/(NR)}'" % crystal_file).read()
        y = os.popen(
            "cat %s | awk '{if ($1==\"ATOM\") print $(NF-5)}' | awk '{y+=$1} END {print y/(NR)}'" % crystal_file).read()
        z = os.popen(
            "cat %s | awk '{if ($1==\"ATOM\") print $(NF-4)}' | awk '{z+=$1} END {print z/(NR)}'" % crystal_file).read()
        xyz = [float(x.strip()), float(y.strip()), float(z.strip())]
        return xyz

    def openeye_convert(self, src_file, dst_file):
        # 指令
        cmd = 'module load openeye&&convert.py {} {}'.format(src_file, dst_file)
        # 执行
        os.system(cmd)

    def vina_lig_prep(self, src_file, dst_file):
        # 指令
        cmd = 'module purge&&module load vina&&timeout 60 prepare_ligand4.py -l {} -o {} -A hydrogens'.format(src_file, dst_file)
        # 执行
        os.system(cmd) 

    def surflex_getpocket(self, ligand, pre_protein, dst_path):
        # 指令
        cmd = 'cd {}&&module purge&&module load tripos/sybylx2.1.1 &&surflex-dock.exe proto {} {} p1' \
              '&&rm -rf p1-*.pdb'.format(dst_path, ligand, pre_protein)
        # 执行
        os.system(cmd)

    def slide_lig_pre(self, src_lig, dst_lig):
        with open(src_lig, 'r') as f:
            lig_data = f.read()
        if '@<TRIPOS>DICT' in lig_data:
            cmdline = 'cat %s | sed \'/@<TRIPOS>DICT/,/@<TRIPOS>ATOM/c @<TRIPOS>ATOM\' > %s ' % (src_lig, dst_lig)
            os.system(cmdline)
        else:
            cmdline = 'cp %s %s' % (src_lig, dst_lig)
            os.system(cmdline)
        newline = []
        with open(dst_lig, 'r') as f:
            lines = f.readlines()
        for line in lines:
            if len(line.split()) >= 3:
                if line.split()[1] == 'OXT':
                    newline.append(line[:47] + 'O.co2' + line[52:])
                else:
                    newline.append(line)
            else:
                newline.append(line)
        with open(dst_lig, 'w') as f:
            f.write(''.join(newline))

    def moe_lig_pre(self, dst_lig_path, src_lig):
        # 创建文件夹
        os.mkdir(dst_lig_path)
        # 定义配置文件
        con = '''#svl
ligmdb = db_Open [ 'l.mdb','create' ];
db_ImportMOL2 ['{}', 'l.mdb', 'mol'];
//db_ImportSD ['l.mdb', '../%s.sdf', 'mol', [], [], [], [append:0,file_field:0,no_field:1]];
[fieldnames, fieldtypes] = db_Fields 'l.mdb';
l = length fieldnames;
newfields = split [fieldnames, [1,l-1]];
db_DeleteField ['l.mdb', newfields(2)];
db_Close ligmdb;'''.format(src_lig)
        # 生成配置文件
        svl = '{}/ligprep.svl'.format(dst_lig_path)
        with open(svl, 'w') as f:
            f.write(con)
        # 执行配置文件
        try:
            cmd = 'cd {}&&module load moe&&moebatch -script ligprep.svl'.format(dst_lig_path)
            os.system(cmd)
        except:
            pass

    def moe_pro_pre(self, path_files, dst_lig_path, src_lig):
        # 定义蛋白
        protein = '{}/protein_file.pdb'.format(path_files)
        # 定义配置文件
        con = '''#svl
rec0 = Open '{}';
lig0 = Open '{}';
c = Chains [];
issues = StructurePreparation [cmd:'cli', batch_protonate3d:1];
SaveAs 'complex_prep.moe';
Close [];'''.format(protein, src_lig)
        # 生成配置文件
        svl = '{}/prot_prep.svl'.format(dst_lig_path)
        with open(svl, 'w') as f:
            f.write(con)
        # 执行配置文件
        try:
            cmd = 'cd {}&&module load moe&&moebatch -script prot_prep.svl'.format(dst_lig_path)
            os.system(cmd)
        except:
            pass

    def moe_mdb2result(self, dst_path, dst_lig_path):
        # 写配置文件
        con = '''#svl
ligmdb = db_Open [ 'dock.mdb','read' ];
db_ExportSD ['dock.mdb', 'dock_ase.sd', [], []];
db_Close ligmdb;
        '''
        # 第一次生成配置文件
        svl = '{}/result_extract.svl'.format(dst_path)
        if not os.path.exists(svl):
            with open(svl, 'w') as f:
                f.write(con)
        # 执行命令
        try:
            cmd = 'cd {}&&module load moe&&moebatch -script {}'.format(dst_lig_path, svl)
            os.system(cmd)
        except:
            pass

    def moe_excute(self, dst_lig_path, dock_batch, score_batch):
        # 命令行
        cmd = 'cd {}&&module load moe&&moebatch -run {} -lig ./l.mdb&&moebatch -run {}'.format(dst_lig_path, dock_batch,
                                                                                               score_batch)
        # 执行命令行
        os.system(cmd)

    def gold_pro_pre(self, src_protein, dst_protein):
        # 蛋白准备脚本
        pre_py = '{}/gold_pro_pre.py'.format(self.help_path)
        # 指令
        cmd = 'module load ccdc&&export CCDC_PYTHON_API_NO_QAPPLICATION=10&&' \
              'python {} {} {}'.format(pre_py, src_protein, dst_protein)
        # 执行指令
        os.system(cmd)

    def gold_excute(self, dst_lig_path, config_file):
        # 指令  加上export CCDC_PYTHON_API_NO_QAPPLICATION=10这条指令才能在节点上运行，否则报错
        cmd = 'cd {}&&module load ccdc&&export CCDC_PYTHON_API_NO_QAPPLICATION=10&&' \
              'gold_auto {}'.format(dst_lig_path, config_file)
        # 执行
        os.system(cmd)

    def dock_pro_pre(self, pre_protein_file, dst_path):
        # 加电荷脚本内容
        con = '''from chimera import runCommand
runCommand("open 1 {0}")
runCommand("addh")
runCommand("addcharge all chargeModel 99sb method am1")
runCommand("write format mol2 atomTypes sybyl 1 {1}/protein_99sb.mol2")
runCommand("delete #1@H@H?@H??@H???")
runCommand("write format pdb 1 {1}/protein_noH.pdb")
runCommand("close 1") 
        '''.format(pre_protein_file, dst_path)
        # 加电荷脚本
        add_py = '{}/prot_addcharge.py'.format(dst_path)
        # 创建脚本
        with open(add_py, 'w') as f:
            f.write(con)
        # 命令行 先后生成了protein_99sb.mol2、protein_noH.pdb
        cmd = 'cd {}&&module load chimera&&chimera --nogui --nostatus --script prot_addcharge.py&&rm -rf ' \
              'prot_addcharge.py'.format(dst_path)
        # 执行
        os.system(cmd)

    def dock_active_site(self, dst_path, crystal_ligand):
        # 创建INSPH文件（后面get_sph中会用到）
        con = '''protein_noH.dms\nR\nX\n0.0\n4.0\n1.4\nprotein_noH.sph'''  # 文件内容
        INSPH = '{}/INSPH'.format(dst_path)  # 定义文件
        with open(INSPH, 'w') as f:  # 创建文件
            f.write(con)
        # 指令  生成protein_noH.dms, protein_noH.sph
        cmd = 'cd {}&&module load dock&&dms ./protein_noH.pdb -n -w 1.4 -o ./protein_noH.dms' \
              '&&sphgen_cpp&&sphere_selector protein_noH.sph {} 10.0'.format(dst_path, crystal_ligand)
        # 执行
        os.system(cmd)

    def dock_lig_addcharge(self, temp_ligand, dst_path, lig_name):
        # 加电荷准备脚本
        pre_py = '{}/dock_pybel.py'.format(self.help_path)
        # 指令
        cmd = "export PYTHONPATH=/opt/openbabel/2.4.1/lib64/python2.6/site-packages:" \
              "${{PYTHONPATH}}&&module load openbabel&&python2 {} {} {} {}".format(pre_py, temp_ligand, dst_path,
                                                                                   lig_name)
        # 执行指令
        os.system(cmd)

    def dock_generate_grid(self, dst_path):
        # 定义showbox文件
        box_pbd = '{}/box.pdb'.format(dst_path)  # 待生成的蛋白
        showbox = '{}/showbox.in'.format(dst_path)  # 配置文件
        # 若文件不存在则创建
        if not os.path.exists(showbox):
            con = '''Y\n8.0\nselected_spheres.sph\n1\n{}'''.format(box_pbd)
            with open(showbox, 'w') as f:
                f.write(con)
        # 定义格点配置文件
        con = '''compute_grids yes
grid_spacing 0.3
output_molecule no
contact_score yes
contact_cutoff_distance	4.5
energy_score yes
energy_cutoff_distance 10
atom_model a
attractive_exponent 6
repulsive_exponent 12
distance_dielectric yes
dielectric_factor 4.0
bump_filter yes
bump_overlap 0.75
receptor_file {}/protein_99sb.mol2
box_file {}
vdw_definition_file /opt/dock/6.8/parameters/vdw_AMBER_parm99.defn
score_grid_prefix grid
        '''.format(dst_path, box_pbd)  # 文件内容
        grid_file = '{}/grid.in'.format(dst_path)  # 文件路径
        if not os.path.exists(grid_file):  # 若不存在则创建
            with open(grid_file, 'w') as f:
                f.write(con)
        # 指令  生成grid.cnt, grid.bmp ，用于后续的打分
        cmd = 'cd {}&&module load dock&&showbox < showbox.in&&grid -i grid.in > gridinfo.out'.format(dst_path)
        # 执行
        os.system(cmd)

    def dock_excute(self, dst_path, lig_name):
        # 定义 指令
        cmd = 'cd {}&&module load dock&&dock6 -i ./{}.in'.format(dst_path, lig_name)
        # 执行
        os.system(cmd)

    def glide_generate_grid(self, pre_protein_file, x, y, z, dst_path):
        # 格点配置文件
        grid_file = '{}/grid.in'.format(dst_path)  # 格点文件
        con = self.glide_grid.format(x, y, z, pre_protein_file)  # 文件内容
        # 生成配置文件
        with open(grid_file, 'w') as f:
            f.write(con)
        # 生成格点文件
        cmd = 'cd {}&&module load schrodinger&&glide {} -NOJOBID'.format(dst_path, grid_file)
        os.system(cmd)

    def glide_excute(self, dst_path, config_file, lig_name):
        # 指令
        cmd = 'cd {0}&&module load schrodinger&&glide {1} -NOJOBID' \
              '&&canvasConvert -imae {2}_raw.maegz -ocsv {2}.csv'.format(dst_path, config_file, lig_name)
        # 执行
        os.system(cmd)

    def plants_excute(self, config_file):
        # 指令
        cmd = 'module load plants &&/opt/plants/1.2/plants --mode rescore {}'.format(config_file)
        # 执行
        os.system(cmd)

    def rdock_excute(self, lig_name, dst_lig_path, config_file, pre_ligand, dockprm):
        # 指令
        cmd = 'cd {0}&&module load rDock&&rbcavity -was -d -r {1}' \
              '&&rbdock -i {2} -o {3} -r {1} -p {4} -n 1 &&sdreport -l {3}.sd > {3}.txt'.format \
            (dst_lig_path, config_file, pre_ligand, lig_name, dockprm)
        # 执行
        os.system(cmd)

    def rfscore_excute(self, path_job, dst_path, score_type):
        # 辅助py
        cal_rf = '{}/cal_rf.py'.format(self.help_path)
        # 指令
        cmd = 'module load anaconda2&&'
        cmd += 'module load openbabel&&'
        cmd += 'export PYTHONPATH=${PYTHONPATH}:/home/xujun/My_project:&&'
        cmd += 'export LD_LIBRARY_PATH=/opt/anaconda2/5.3.0/lib:${LD_LIBRARY_PATH}&&'
        cmd += 'python {} {} {} {}'.format(cal_rf, path_job, dst_path, score_type)
        # 执行
        os.system(cmd)

    def oddt_excute(self, protein, ligand, csv_file, fp_type):
        # 辅助py
        cal_ifp = '{}/cal_ifp.py'.format(self.help_path)
        # 指令
        cmd = '{} {} {} {} {} {}'.format(self.ifp_module, cal_ifp, protein, ligand, csv_file, fp_type)
        # 执行
        os.system(cmd)

    def rdkit_excute(self, pre_ligand, csv_file, fp_type):
        # 辅助py
        cal_fp = '{}/cal_fp.py'.format(self.help_path)
        # 指令
        cmd = '{} {} {} {} {}'.format(self.ifp_module, cal_fp, pre_ligand, csv_file, fp_type)
        # 执行
        os.system(cmd)

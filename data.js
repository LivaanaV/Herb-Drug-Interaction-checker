/* HDI Engine v2 — seed dataset
 *
 * PROVENANCE IS PART OF THE DATA MODEL.
 * Every numeric parameter carries `src`, one of:
 *   'trial'        a human DDI crossover study (none seeded yet — curate from PubMed)
 *   'consensus'    a standard, widely reported literature value; verify before publication
 *   'illustrative' a placeholder of the right order of magnitude; MUST be replaced
 *
 * The engine propagates this. Any prediction touching an 'illustrative' value is
 * marked provisional and capped at LOW confidence, whatever the arithmetic says.
 */
(function (root) {
  'use strict';

  // ---------------------------------------------------------------- proteins
  // kdeg in h^-1. Required by every time-dependent inactivation term.
  const PROTEINS = {
    CYP3A4_H:  { symbol: 'CYP3A4', uniprot: 'P08684', chembl: 'CHEMBL340',  organ: 'hepatic',    kind: 'enzyme',      kdeg: 0.0193, kdegSrc: 'consensus' },
    CYP3A4_G:  { symbol: 'CYP3A4', uniprot: 'P08684', chembl: 'CHEMBL340',  organ: 'intestinal', kind: 'enzyme',      kdeg: 0.0300, kdegSrc: 'consensus' },
    CYP2C9_H:  { symbol: 'CYP2C9', uniprot: 'P11712', chembl: 'CHEMBL3397',  organ: 'hepatic',    kind: 'enzyme',      kdeg: 0.0068, kdegSrc: 'illustrative' },
    CYP2C19_H: { symbol: 'CYP2C19', uniprot: 'P33261', chembl: 'CHEMBL3622', organ: 'hepatic',    kind: 'enzyme',      kdeg: 0.0257, kdegSrc: 'illustrative' },
    CYP2D6_H:  { symbol: 'CYP2D6', uniprot: 'P10635', chembl: 'CHEMBL289',  organ: 'hepatic',    kind: 'enzyme',      kdeg: 0.0106, kdegSrc: 'illustrative' },
    CYP1A2_H:  { symbol: 'CYP1A2', uniprot: 'P05177', chembl: 'CHEMBL3356',  organ: 'hepatic',    kind: 'enzyme',      kdeg: 0.0183, kdegSrc: 'illustrative' },
    CYP2B6_H:  { symbol: 'CYP2B6', uniprot: 'P20813', chembl: 'CHEMBL4729',  organ: 'hepatic',    kind: 'enzyme',      kdeg: 0.0217, kdegSrc: 'illustrative' },
    UGT1A1_H:  { symbol: 'UGT1A1', uniprot: 'P22309', chembl: null,  organ: 'hepatic',    kind: 'enzyme',      kdeg: 0.0200, kdegSrc: 'illustrative' },
    ABCB1_G:   { symbol: 'ABCB1', uniprot: 'P08183', chembl: null,   organ: 'intestinal', kind: 'transporter', note: 'P-glycoprotein efflux' },
    ABCG2_G:   { symbol: 'ABCG2', uniprot: null, chembl: null,   organ: 'intestinal', kind: 'transporter', note: 'BCRP efflux' },
    SLCO1B1_H: { symbol: 'SLCO1B1', uniprot: 'Q9Y6L6', chembl: null, organ: 'hepatic',    kind: 'transporter', note: 'OATP1B1 uptake' }
  };

  // ------------------------------------------------------------ constituents
  // pk.cmax in µM at pk.dose mg. faFg = fraction absorbed x fraction escaping
  // gut metabolism, for the CONSTITUENT. Known faFg promotes a constituent from
  // an interval (tier B) to a point estimate (tier A).
  const CONSTITUENTS = {
    curcumin:            { name: 'Curcumin',              mw: 368.4, solubility_uM: 9.5, solubility_note: 'Reported aqueous solubility spans 11 ng/mL (0.03 uM, PubChem) to 3.5 ug/mL (9.5 uM). The higher, more conservative figure is used. Curcumin is BCS class IV.', pk: null,
      pk_note: 'Vareed 2008 (PMID 18559556) gave 10 g and 12 g doses. Free curcumin was below the 50 ng/mL detection limit at all but one of 14 time points, so no free-plasma Cmax is supported. Only glucuronide and sulfate conjugates were quantified. The earlier 0.011 uM value is removed because its citation does not report it. Curcumin is gut-only here (tier C).' },
    demethoxycurcumin:   { name: 'Demethoxycurcumin',     mw: 338.4, solubility_uM: 19.4, solubility_note: '6.56e-3 mg/mL aqueous.', pk: null },
    bisdemethoxycurcumin:{ name: 'Bisdemethoxycurcumin',  mw: 308.3, solubility_uM: 48.7, solubility_note: '0.015 mg/mL aqueous; the least insoluble of the three curcuminoids.', pk: null },
    piperine:            { name: 'Piperine',              mw: 285.3, pk: { cmax: 0.18,  fu: 0.03,  dose: 20,   faFg: null, src: 'illustrative' } },
    withaferin_a:        { name: 'Withaferin A',          mw: 470.6, pk: null },
    withanolide_a:       { name: 'Withanolide A',         mw: 470.6, pk: null },
    hyperforin:          { name: 'Hyperforin',            mw: 536.8, pk: { cmax: 0.28,  fu: 0.005, dose: 14.8, faFg: null, src: 'sourced', pmid: '9684946', fuSrc: 'illustrative',
      note: 'Biber et al. 1998 (Pharmacopsychiatry 31 Suppl 1:36-43, PMID 9684946): about 150 ng/mL (about 280 nM) at 3.5 h after a 300 mg extract representing 14.8 mg hyperforin; half-life 9 h; kinetics linear to 600 mg extract. Fa*Fg and fu are not in the abstract, so both are unsourced (fa*Fg null, fu illustrative).' } },
    hypericin:           { name: 'Hypericin',             mw: 504.4, pk: { cmax: 0.01,  fu: 0.005, dose: 900,  faFg: null, src: 'illustrative' } },
    gingerol_6:          { name: '6-Gingerol',            mw: 294.4, pk: null },
    shogaol_6:           { name: '6-Shogaol',             mw: 276.4, pk: null },
    allicin:             { name: 'Allicin',               mw: 162.3, pk: null },
    diallyl_disulfide:   { name: 'Diallyl disulfide',     mw: 146.3, pk: null },
    silybin:             { name: 'Silybin',               mw: 482.4, pk: { cmax: 0.28,  fu: 0.05,  dose: 420,  faFg: null, src: 'partial', pmid: 'PMC3876803', fuSrc: 'illustrative',
      note: 'Steady-state silybin A Cmax 134.7 ng/mL (0.28 uM) on silymarin 140 mg three times daily for 28 days (PMC3876803, silybin A only; silybin B reported separately). Dose is 420 mg silymarin/day, so the scaling mass is silymarin, not silybin, which makes the scaling approximate. Absolute oral bioavailability (0.45-0.95%) is not Fa*Fg (point 30), so Fa*Fg is null. fu is unsourced.' } },
    glycyrrhizin:        { name: 'Glycyrrhizin',          mw: 822.9, pk: { cmax: 2.00,  fu: 0.02,  dose: 3000, faFg: null, src: 'illustrative' } },
    liquiritigenin:      { name: 'Liquiritigenin',        mw: 256.3, pk: null },
    berberine:           { name: 'Berberine',             mw: 336.4, pk: { cmax: 0.004, fu: 0.05,  dose: 1000, faFg: null, src: 'illustrative' } },
    tinosporaside:       { name: 'Tinosporaside',         mw: 450.4, pk: null },
    bacoside_a:          { name: 'Bacoside A',            mw: 929.1, pk: null },
    gallic_acid:         { name: 'Gallic acid',           mw: 170.1, pk: { cmax: 1.50,  fu: 0.30,  dose: 3000, faFg: null, src: 'illustrative' } },
    ellagic_acid:        { name: 'Ellagic acid',          mw: 302.2, pk: { cmax: 0.10,  fu: 0.10,  dose: 3000, faFg: null, src: 'illustrative' } },
    eugenol:             { name: 'Eugenol',               mw: 164.2, pk: null },
    ursolic_acid:        { name: 'Ursolic acid',          mw: 456.7, pk: null },
    trigonelline:        { name: 'Trigonelline',          mw: 137.1, pk: { cmax: 3.00,  fu: 0.60,  dose: 5000, faFg: null, src: 'illustrative' } },
    diosgenin:           { name: 'Diosgenin',             mw: 414.6, pk: null }
  };

  // ------------------------------------------------------------ preparations
  // The query entity. Species + plant part + preparation type + daily dose.
  // `pct` is % w/w of the preparation. `marker` means a pharmacopoeial assay
  // constituent, which passes the gate even without a quantified abundance.
  const PREPARATIONS = [
    { id: 'turmeric-std', species: 'Curcuma longa', common: 'Turmeric', vernacular: 'Haldi / Haridra',
      part: 'rhizome', prep: 'standardised curcuminoid extract', dose: 1500, doseSrc: 'consensus',
      constituents: [
        { cid: 'curcumin',             pct: 75,  src: 'uncited',    marker: true },
        { cid: 'demethoxycurcumin',    pct: 15,  src: 'uncited',    marker: true },
        { cid: 'bisdemethoxycurcumin', pct: 5,   src: 'uncited',    marker: true }
      ] },
    { id: 'turmeric-powder', species: 'Curcuma longa', common: 'Turmeric', vernacular: 'Haldi / Haridra',
      part: 'rhizome', prep: 'dried powder', dose: 3000, doseSrc: 'illustrative',
      constituents: [
        { cid: 'curcumin',             pct: 3.0, src: 'illustrative', marker: true },
        { cid: 'demethoxycurcumin',    pct: 0.6, src: 'illustrative', marker: false },
        { cid: 'bisdemethoxycurcumin', pct: 0.2, src: 'illustrative', marker: false }
      ] },
    { id: 'blackpepper-fruit', species: 'Piper nigrum', common: 'Black pepper', vernacular: 'Kali mirch / Maricha',
      part: 'fruit', prep: 'dried powder', dose: 1000, doseSrc: 'illustrative',
      constituents: [ { cid: 'piperine', pct: 5, src: 'uncited', marker: true } ] },
    { id: 'piperine-std', species: 'Piper nigrum', common: 'Black pepper', vernacular: 'Kali mirch / Maricha',
      part: 'fruit', prep: 'standardised piperine extract', dose: 20, doseSrc: 'consensus',
      constituents: [ { cid: 'piperine', pct: 95, src: 'uncited', marker: true } ] },
    { id: 'ashwagandha-root', species: 'Withania somnifera', common: 'Ashwagandha', vernacular: 'Asgandh',
      part: 'root', prep: 'hydroalcoholic extract', dose: 600, doseSrc: 'consensus',
      constituents: [
        { cid: 'withaferin_a',  pct: 0.3, src: 'illustrative', marker: true },
        { cid: 'withanolide_a', pct: 0.5, src: 'illustrative', marker: true }
      ] },
    { id: 'stjohnswort-std', species: 'Hypericum perforatum', common: "St John's wort", vernacular: 'reference inducer',
      part: 'aerial parts', prep: 'standardised extract', dose: 900, doseSrc: 'consensus',
      constituents: [
        { cid: 'hyperforin', pct: 3.0, src: 'uncited', marker: true },
        { cid: 'hypericin',  pct: 0.3, src: 'uncited', marker: true }
      ] },
    { id: 'ginger-rhizome', species: 'Zingiber officinale', common: 'Ginger', vernacular: 'Adrak / Shunthi',
      part: 'rhizome', prep: 'dried powder', dose: 2000, doseSrc: 'illustrative',
      constituents: [
        { cid: 'gingerol_6', pct: 1.5, src: 'illustrative', marker: true },
        { cid: 'shogaol_6',  pct: 0.3, src: 'illustrative', marker: false }
      ] },
    { id: 'garlic-bulb', species: 'Allium sativum', common: 'Garlic', vernacular: 'Lahsun / Rasona',
      part: 'bulb', prep: 'dried powder', dose: 1000, doseSrc: 'illustrative',
      constituents: [
        { cid: 'allicin',           pct: 0.4, src: 'illustrative', marker: true },
        { cid: 'diallyl_disulfide', pct: 0.2, src: 'illustrative', marker: false }
      ] },
    { id: 'milkthistle-std', species: 'Silybum marianum', common: 'Milk thistle', vernacular: 'reference inhibitor',
      part: 'seed', prep: 'standardised silymarin extract', dose: 420, doseSrc: 'consensus',
      constituents: [ { cid: 'silybin', pct: 30, src: 'uncited', marker: true } ] },
    { id: 'licorice-root', species: 'Glycyrrhiza glabra', common: 'Licorice', vernacular: 'Mulethi / Yashtimadhu',
      part: 'root', prep: 'dried powder', dose: 3000, doseSrc: 'illustrative',
      constituents: [
        { cid: 'glycyrrhizin',   pct: 4.0, src: 'uncited',    marker: true },
        { cid: 'liquiritigenin', pct: 0.3, src: 'illustrative', marker: false }
      ] },
    { id: 'guduchi-stem', species: 'Tinospora cordifolia', common: 'Guduchi', vernacular: 'Giloy',
      part: 'stem', prep: 'aqueous extract', dose: 1000, doseSrc: 'illustrative',
      constituents: [
        { cid: 'berberine',     pct: 0.1, src: 'illustrative', marker: true },
        { cid: 'tinosporaside', pct: 0.2, src: 'illustrative', marker: false }
      ] },
    { id: 'brahmi-std', species: 'Bacopa monnieri', common: 'Brahmi', vernacular: 'Jalnaveri',
      part: 'whole plant', prep: 'standardised extract', dose: 300, doseSrc: 'consensus',
      constituents: [ { cid: 'bacoside_a', pct: 50, src: 'uncited', marker: true } ] },
    { id: 'amla-fruit', species: 'Phyllanthus emblica', common: 'Amla', vernacular: 'Amalaki',
      part: 'fruit', prep: 'dried powder', dose: 3000, doseSrc: 'illustrative',
      constituents: [
        { cid: 'gallic_acid',  pct: 1.2, src: 'illustrative', marker: true },
        { cid: 'ellagic_acid', pct: 0.5, src: 'illustrative', marker: false }
      ] },
    { id: 'tulsi-leaf', species: 'Ocimum tenuiflorum', common: 'Tulsi', vernacular: 'Holy basil',
      part: 'leaf', prep: 'dried powder', dose: 2000, doseSrc: 'illustrative',
      constituents: [
        { cid: 'eugenol',      pct: 0.5, src: 'illustrative', marker: true },
        { cid: 'ursolic_acid', pct: 0.4, src: 'illustrative', marker: false }
      ] },
    { id: 'fenugreek-seed', species: 'Trigonella foenum-graecum', common: 'Fenugreek', vernacular: 'Methi',
      part: 'seed', prep: 'dried powder', dose: 5000, doseSrc: 'illustrative',
      constituents: [
        { cid: 'trigonelline', pct: 0.4, src: 'illustrative', marker: true },
        { cid: 'diosgenin',    pct: 0.5, src: 'illustrative', marker: false }
      ] }
  ];

  // -------------------------------------------------------------------- drugs
  // fg  = fraction of the dose escaping gut metabolism
  // fm  = fraction of systemic clearance through that enzyme
  // nti = graded, with the fold-change at which the drug becomes clinically unsafe
  const DRUGS = [
    { id: 'midazolam', inn: 'Midazolam', brands: ['Fulsed', 'Mezolam'], note: 'FDA index sensitive substrate, CYP3A4',
      moieties: [ { name: 'midazolam', fg: 0.57, fgSrc: 'sourced', fgRef: 'Paine et al. 1996, Clin Pharmacol Ther 60:14-24, doi:10.1016/s0009-9236(96)90162-9 (mean intestinal extraction 0.43)', gutEnzyme: 'CYP3A4_G',
        clearance: [ { protein: 'CYP3A4_H', fm: 0.93, src: 'uncited' } ], transporters: [], nti: null } ] },

    { id: 'warfarin', inn: 'Warfarin', brands: ['Warf', 'Sofarin', 'Uniwarfin'], note: 'Enantiomers cleared by different enzymes',
      moieties: [
        { name: 'S-warfarin', fg: 0.93, fgSrc: 'uncited', gutEnzyme: null,
          clearance: [ { protein: 'CYP2C9_H', fm: 0.85, src: 'uncited' } ], transporters: [],
          nti: { grade: 'high', threshold: 1.25 } },
        { name: 'R-warfarin', fg: 0.93, fgSrc: 'uncited', gutEnzyme: null,
          clearance: [ { protein: 'CYP1A2_H', fm: 0.40, src: 'uncited' }, { protein: 'CYP3A4_H', fm: 0.40, src: 'uncited' } ],
          transporters: [], nti: { grade: 'high', threshold: 1.25 } }
      ] },

    { id: 'digoxin', inn: 'Digoxin', brands: ['Lanoxin', 'Dixin'], note: 'Index P-gp probe; renally cleared, no CYP pathway',
      moieties: [ { name: 'digoxin', fg: 1.0, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [], transporters: ['ABCB1_G'], nti: { grade: 'high', threshold: 1.30 } } ] },

    { id: 'omeprazole', inn: 'Omeprazole', brands: ['Omez', 'Ocid', 'Omecip'], note: 'FDA index substrate, CYP2C19',
      moieties: [ { name: 'omeprazole', fg: 0.90, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [ { protein: 'CYP2C19_H', fm: 0.80, src: 'uncited' }, { protein: 'CYP3A4_H', fm: 0.20, src: 'uncited' } ],
        transporters: [], nti: null } ] },

    { id: 'dextromethorphan', inn: 'Dextromethorphan', brands: ['Benadryl DR', 'Tossex'], note: 'FDA index probe, CYP2D6',
      moieties: [ { name: 'dextromethorphan', fg: 0.60, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [ { protein: 'CYP2D6_H', fm: 0.90, src: 'uncited' } ], transporters: [], nti: null } ] },

    { id: 'theophylline', inn: 'Theophylline', brands: ['Deriphyllin', 'Theobid'], note: 'NTI, CYP1A2',
      moieties: [ { name: 'theophylline', fg: 1.0, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [ { protein: 'CYP1A2_H', fm: 0.80, src: 'uncited' } ], transporters: [],
        nti: { grade: 'high', threshold: 1.30 } } ] },

    { id: 'caffeine', inn: 'Caffeine', brands: ['Caffeine citrate'], note: 'CYP1A2 probe used in most human cocktail studies',
      moieties: [ { name: 'caffeine', fg: 1.0, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [ { protein: 'CYP1A2_H', fm: 0.95, src: 'uncited' } ], transporters: [], nti: null } ] },

    { id: 'simvastatin', inn: 'Simvastatin', brands: ['Simvotin', 'Simlup'], note: 'Very low Fg, so intestinal inhibition dominates',
      moieties: [ { name: 'simvastatin', fg: 0.34, fgSrc: 'uncited', gutEnzyme: 'CYP3A4_G',
        clearance: [ { protein: 'CYP3A4_H', fm: 0.80, src: 'uncited' } ], transporters: ['SLCO1B1_H'], nti: null } ] },

    { id: 'atorvastatin', inn: 'Atorvastatin', brands: ['Atorva', 'Storvas', 'Lipicure'], note: 'Enzyme and uptake-transporter substrate',
      moieties: [ { name: 'atorvastatin', fg: 0.24, fgSrc: 'uncited', gutEnzyme: 'CYP3A4_G',
        clearance: [ { protein: 'CYP3A4_H', fm: 0.68, src: 'uncited' } ], transporters: ['SLCO1B1_H'], nti: null } ] },

    { id: 'rosuvastatin', inn: 'Rosuvastatin', brands: ['Rosuvas', 'Rozavel'], note: 'Transporter-mediated only; minimal CYP',
      moieties: [ { name: 'rosuvastatin', fg: 0.90, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [], transporters: ['ABCG2_G', 'SLCO1B1_H'], nti: null } ] },

    { id: 'clopidogrel', inn: 'Clopidogrel', brands: ['Clopilet', 'Deplatt', 'Plavix'], note: 'Prodrug — loss of activation is the hazard',
      isProdrug: true, criticalEfficacy: true,
      moieties: [ { name: 'clopidogrel active metabolite', fg: 0.50, fgSrc: 'illustrative', gutEnzyme: null,
        clearance: [], activation: [ { protein: 'CYP2C19_H', fAct: 0.45, src: 'uncited' }, { protein: 'CYP3A4_H', fAct: 0.35, src: 'uncited' } ],
        transporters: [], nti: null } ] },

    { id: 'tacrolimus', inn: 'Tacrolimus', brands: ['Pangraf', 'Tacrograf'], note: 'NTI, very low Fg, graft rejection on loss of exposure',
      moieties: [ { name: 'tacrolimus', fg: 0.14, fgSrc: 'uncited', gutEnzyme: 'CYP3A4_G',
        clearance: [ { protein: 'CYP3A4_H', fm: 0.95, src: 'uncited' } ], transporters: ['ABCB1_G'],
        nti: { grade: 'high', threshold: 1.20 } } ] },

    { id: 'cyclosporine', inn: 'Cyclosporine', brands: ['Panimun Bioral', 'Sandimmun'], note: 'NTI, enzyme and transporter substrate',
      moieties: [ { name: 'cyclosporine', fg: 0.44, fgSrc: 'uncited', gutEnzyme: 'CYP3A4_G',
        clearance: [ { protein: 'CYP3A4_H', fm: 0.80, src: 'uncited' } ], transporters: ['ABCB1_G'],
        nti: { grade: 'high', threshold: 1.20 } } ] },

    { id: 'phenytoin', inn: 'Phenytoin', brands: ['Eptoin', 'Dilantin'], note: 'NTI with non-linear kinetics — small clearance changes amplify',
      moieties: [ { name: 'phenytoin', fg: 1.0, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [ { protein: 'CYP2C9_H', fm: 0.80, src: 'uncited' }, { protein: 'CYP2C19_H', fm: 0.15, src: 'uncited' } ],
        transporters: [], nti: { grade: 'very high', threshold: 1.15 } } ] },

    { id: 'alprazolam', inn: 'Alprazolam', brands: ['Alprax', 'Restyl'], note: 'CYP3A4 substrate used in human interaction studies',
      moieties: [ { name: 'alprazolam', fg: 0.90, fgSrc: 'uncited', gutEnzyme: 'CYP3A4_G',
        clearance: [ { protein: 'CYP3A4_H', fm: 0.80, src: 'uncited' } ], transporters: [], nti: null } ] },

    { id: 'tolbutamide', inn: 'Tolbutamide', brands: ['Rastinon'], note: 'CYP2C9 probe',
      moieties: [ { name: 'tolbutamide', fg: 1.0, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [ { protein: 'CYP2C9_H', fm: 0.85, src: 'uncited' } ], transporters: [], nti: null } ] },

    { id: 'bupropion', inn: 'Bupropion', brands: ['Bupron SR', 'Smoquit'], note: 'CYP2B6 probe',
      moieties: [ { name: 'bupropion', fg: 0.90, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [ { protein: 'CYP2B6_H', fm: 0.85, src: 'uncited' } ], transporters: [], nti: null } ] },

    { id: 'raltegravir', inn: 'Raltegravir', brands: ['Isentress'], note: 'UGT1A1 probe — the phase II pathway',
      moieties: [ { name: 'raltegravir', fg: 1.0, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [ { protein: 'UGT1A1_H', fm: 0.70, src: 'uncited' } ], transporters: [], nti: null } ] },

    { id: 'metformin', inn: 'Metformin', brands: ['Glycomet', 'Glucophage', 'Obimet'], note: 'Negative control — renal clearance, no hepatic oxidation',
      moieties: [ { name: 'metformin', fg: 1.0, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [], transporters: [], nti: null } ] },

    { id: 'lithium', inn: 'Lithium carbonate', brands: ['Lithosun', 'Licab'], note: 'Negative control — renal clearance',
      moieties: [ { name: 'lithium', fg: 1.0, fgSrc: 'uncited', gutEnzyme: null,
        clearance: [], transporters: [], nti: { grade: 'high', threshold: 1.25 } } ] }
  ];

  // ------------------------------------------------------- typed bioactivity
  // Reversible inhibition. ki in µM.
  const REVERSIBLE = [
    // --- Volak 2008, PMID 18480186 — curcuminoids and piperine, human liver microsomes
    { cid: 'curcumin',          protein: 'CYP3A4_H',  ki: 11.0, ic50: 25.3, ki_basis: 'Ki of curcuminoid extract (Volak 2008), applied to curcumin', mode: 'competitive',    probe: 'triazolam 1-OH',  system: 'HLM',    species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'curcumin',          protein: 'CYP3A4_G',  ki: 11.0, ic50: 25.3, ki_basis: 'Ki of curcuminoid extract (Volak 2008); hepatic value assumed for intestinal enzyme', mode: 'competitive', probe: 'triazolam 1-OH', system: 'HLM', species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'curcumin',          protein: 'UGT1A1_H',  ki: 1.1,  ic50: 2.2,  ki_basis: 'IC50/2; acetaminophen glucuronidation is not UGT1A1-selective', mode: 'n/a', probe: 'acetaminophen glucuronidation', system: 'intact LS180 cells', species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'demethoxycurcumin', protein: 'CYP3A4_H',  ki: 21.8, ic50: 43.5, ki_basis: 'IC50/2', mode: 'competitive',    probe: 'triazolam 1-OH',  system: 'HLM',    species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'demethoxycurcumin', protein: 'CYP3A4_G',  ki: 21.8, ic50: 43.5, ki_basis: 'IC50/2, hepatic value assumed for intestinal enzyme', mode: 'competitive', probe: 'triazolam 1-OH', system: 'HLM', species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'demethoxycurcumin', protein: 'CYP2C9_H',  ki: 4.4,  ic50: 8.8,  ki_basis: 'IC50/2', mode: 'competitive',    probe: 'flurbiprofen OH', system: 'HLM',    species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'demethoxycurcumin', protein: 'UGT1A1_H',  ki: 7.0,  ic50: 13.9, ki_basis: 'IC50/2; probe not UGT1A1-selective', mode: 'n/a', probe: 'acetaminophen glucuronidation', system: 'intact LS180 cells', species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'bisdemethoxycurcumin', protein: 'CYP3A4_H', ki: 20.9, ic50: 41.7, ki_basis: 'IC50/2', mode: 'competitive',  probe: 'triazolam 1-OH',  system: 'HLM',    species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'bisdemethoxycurcumin', protein: 'CYP3A4_G', ki: 20.9, ic50: 41.7, ki_basis: 'IC50/2, hepatic value assumed for intestinal enzyme', mode: 'competitive', probe: 'triazolam 1-OH', system: 'HLM', species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'bisdemethoxycurcumin', protein: 'UGT1A1_H', ki: 4.9, ic50: 9.8, ki_basis: 'IC50/2; probe not UGT1A1-selective', mode: 'n/a', probe: 'acetaminophen glucuronidation', system: 'intact LS180 cells', species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'piperine',          protein: 'CYP3A4_H',  ki: 5.4,  ic50: 5.5,  ki_basis: 'measured', mode: 'noncompetitive', probe: 'triazolam 1-OH', system: 'HLM',   species: 'human', src: 'sourced', pmid: '18480186' },
    { cid: 'piperine',          protein: 'CYP3A4_G',  ki: 5.4,  ic50: 5.5,  ki_basis: 'measured, hepatic value assumed for intestinal enzyme', mode: 'noncompetitive', probe: 'triazolam 1-OH', system: 'HLM', species: 'human', src: 'sourced', pmid: '18480186' },
    { cid: 'piperine',          protein: 'CYP2C9_H',  ki: 20.4, ic50: 40.7, ki_basis: 'IC50/2', mode: 'competitive',    probe: 'flurbiprofen OH', system: 'HLM',    species: 'human', src: 'partial', pmid: '18480186' },
    { cid: 'piperine',          protein: 'CYP1A2_H',  ki: 14.9, ic50: 29.8, ki_basis: 'IC50/2', mode: 'competitive',    probe: 'phenacetin O-de', system: 'HLM',    species: 'human', src: 'partial', pmid: '18480186' },
    // Volak 2008 found NO piperine inhibition of UGT or SULT up to 50 uM, and no
    // piperine inhibition of CYP2C19, CYP2D6, CYP2E1 or CYP2B6 (all IC50 > 50 uM).
    // The placeholder piperine-UGT1A1 record that used to sit here was wrong and is gone.

    // --- still placeholders, awaiting curation
    { cid: 'curcumin',          protein: 'ABCB1_G',   ki: 10,  mode: 'n/a',         probe: 'digoxin flux',    system: 'Caco-2',      species: 'human', src: 'illustrative' },
    { cid: 'curcumin',          protein: 'SLCO1B1_H', ki: 2,   mode: 'n/a',         probe: 'estrone-3S',      system: 'HEK-OATP1B1', species: 'human', src: 'illustrative' },
    { cid: 'piperine',          protein: 'ABCB1_G',   ki: 15,  mode: 'n/a',         probe: 'digoxin flux',    system: 'Caco-2',      species: 'human', src: 'illustrative' },
    { cid: 'withaferin_a',      protein: 'CYP3A4_H',  ki: 20,  mode: 'competitive', probe: 'testosterone 6B', system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'withanolide_a',     protein: 'CYP2C19_H', ki: 30,  mode: 'competitive', probe: 'S-mephenytoin',   system: 'HLM',         species: 'human', src: 'illustrative' },
    // Komoroski 2004 (DMD 32:512): hypericin had NO effect on CYP3A4, CYP2C9,
    // CYP1A2 or CYP2D6 in human hepatocytes. The placeholder hypericin-CYP3A4
    // record previously here was wrong and is deleted.
    { cid: 'hypericin',         protein: 'ABCB1_G',   ki: 20,  mode: 'n/a',         probe: 'digoxin flux',    system: 'Caco-2',      species: 'human', src: 'illustrative' },
    { cid: 'silybin',           protein: 'UGT1A1_H',  ki: 1.4, ic50: 1.4, ki_basis: 'IC50 used as Ki proxy (Sridar 2004 reports IC50 only)', mode: 'competitive', probe: '7-hydroxy-4-trifluoromethylcoumarin glucuronidation', system: 'recombinant UGT1A1', species: 'human', src: 'partial', pmid: '15155549' },
    { cid: 'silybin',           protein: 'CYP2C9_H',  ki: 18.0, ic50: 18, ki_basis: 'IC50 used as Ki proxy (JPET 2010, silybin A, human liver microsomes, IC50 18 uM)', mode: 'competitive', probe: '(S)-warfarin 7-hydroxylation', system: 'HLM', species: 'human', src: 'partial', source_id: 'PMC2835426' },
    // Sridar 2004 also reports UGT1A6 28, UGT1A9 20, UGT2B7 92, UGT2B15 75 uM —
    // silybin is 14- to 20-fold selective for UGT1A1 over UGT1A9 and UGT1A6.
    // Silybin inhibition of CYP2C9 and CYP3A4 is mechanism-based, so it is held
    // in the TDI table below rather than as a reversible record.
    { cid: 'silybin',           protein: 'SLCO1B1_H', ki: 1.0, mode: 'n/a',         probe: 'estrone-3S',      system: 'HEK-OATP1B1', species: 'human', src: 'illustrative' },
    { cid: 'glycyrrhizin',      protein: 'CYP3A4_H',  ki: 40,  mode: 'competitive', probe: 'midazolam 1-OH',  system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'liquiritigenin',    protein: 'CYP2C9_H',  ki: 8,   mode: 'competitive', probe: 'diclofenac 4-OH', system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'berberine',         protein: 'CYP2D6_H',  ki: 4,   mode: 'competitive', probe: 'dextrorphan',     system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'berberine',         protein: 'CYP3A4_H',  ki: 15,  mode: 'competitive', probe: 'midazolam 1-OH',  system: 'HLM',         species: 'human', src: 'illustrative' },
    // Ginger constituents. The placeholders here were several-fold too potent.
    // No pure-compound CYP2C9 value was found for 6-gingerol, so that record is
    // gone rather than guessed; the CYP2C9 signal attributed to ginger comes
    // from whole extract and from 6-shogaol.
    { cid: 'gingerol_6',        protein: 'CYP3A4_H',  ki: 54,   ic50: 108,   ki_basis: 'IC50/2', mode: 'competitive', probe: 'CYP3A4 marker', system: 'HLM', species: 'human', src: 'partial', note: 'Not reproduced from the abstracts reviewed. Nature APS 2013 reports 8-gingerol values, not these, and no CYP1A2 result.' },
    { cid: 'gingerol_6',        protein: 'CYP3A4_G',  ki: 54,   ic50: 108,   ki_basis: 'IC50/2, hepatic value assumed for intestinal enzyme', mode: 'competitive', probe: 'CYP3A4 marker', system: 'HLM', species: 'human', src: 'partial', note: 'Not reproduced from the abstracts reviewed. Nature APS 2013 reports 8-gingerol values, not these, and no CYP1A2 result.' },
    { cid: 'gingerol_6',        protein: 'CYP2C19_H', ki: 18,   ic50: 36,    ki_basis: 'IC50/2', mode: 'competitive', probe: 'S-mephenytoin', system: 'HLM', species: 'human', src: 'partial', note: 'Not reproduced from the abstracts reviewed. Nature APS 2013 reports 8-gingerol values, not these, and no CYP1A2 result.' },
    { cid: 'gingerol_6',        protein: 'CYP1A2_H',  ki: 25.5, ic50: 51,    ki_basis: 'IC50/2', mode: 'competitive', probe: 'phenacetin O-de', system: 'HLM', species: 'human', src: 'partial', note: 'Not reproduced from the abstracts reviewed. Nature APS 2013 reports 8-gingerol values, not these, and no CYP1A2 result.' },
    { cid: 'shogaol_6',         protein: 'CYP2C9_H',  ki: 29.02, ic50: 29.20, ki_basis: 'measured Ki, noncompetitive (Kim 2016)', mode: 'noncompetitive', probe: 'CYP2C9 marker', system: 'HLM', species: 'human', src: 'sourced', source_id: 'JAKO201608450940050' },
    { cid: 'shogaol_6',         protein: 'CYP1A2_H',  ki: 6.33,  ic50: 20.68, ki_basis: 'measured Ki, competitive (Kim 2016)', mode: 'competitive', probe: 'phenacetin O-de', system: 'HLM', species: 'human', src: 'sourced', source_id: 'JAKO201608450940050' },
    { cid: 'shogaol_6',         protein: 'CYP2C19_H', ki: 19.26, ic50: 18.78, ki_basis: 'measured Ki, noncompetitive (Kim 2016)', mode: 'noncompetitive', probe: 'S-mephenytoin', system: 'HLM', species: 'human', src: 'sourced', source_id: 'JAKO201608450940050' },
    { cid: 'allicin',           protein: 'CYP3A4_H',  ki: 45,  mode: 'competitive', probe: 'midazolam 1-OH',  system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'diallyl_disulfide', protein: 'CYP2C19_H', ki: 20,  mode: 'competitive', probe: 'S-mephenytoin',   system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'gallic_acid',       protein: 'CYP2C9_H',  ki: 50,  mode: 'competitive', probe: 'diclofenac 4-OH', system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'ellagic_acid',      protein: 'CYP1A2_H',  ki: 6,   mode: 'competitive', probe: 'phenacetin O-de', system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'eugenol',           protein: 'CYP2C9_H',  ki: 25,  mode: 'competitive', probe: 'diclofenac 4-OH', system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'ursolic_acid',      protein: 'CYP1A2_H',  ki: 18,  mode: 'competitive', probe: 'phenacetin O-de', system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'trigonelline',      protein: 'CYP3A4_H',  ki: 80,  mode: 'competitive', probe: 'midazolam 1-OH',  system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'diosgenin',         protein: 'CYP3A4_H',  ki: 40,  mode: 'competitive', probe: 'midazolam 1-OH',  system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'bacoside_a',        protein: 'CYP3A4_H',  ki: 35,  mode: 'competitive', probe: 'midazolam 1-OH',  system: 'HLM',         species: 'human', src: 'illustrative' },
    { cid: 'tinosporaside',     protein: 'CYP3A4_H',  ki: 60,  mode: 'competitive', probe: 'midazolam 1-OH',  system: 'HLM',         species: 'human', src: 'illustrative' }
  ];

  // Time-dependent inactivation. KI in µM, kinact in min^-1.
  const TDI = [
    // Volak 2008 (PMID 18480186): preincubation did NOT enhance piperine inhibition of
    // CYP3A, so the placeholder piperine TDI records previously here were wrong and removed.
    { cid: 'hyperforin', protein: 'CYP3A4_H', KI: 8,  kinact: 0.030, system: 'HLM', species: 'human', src: 'illustrative' },
    // Sridar 2004, PMID 15155549 — silybin is a mechanism-based inactivator.
    { cid: 'silybin',    protein: 'CYP3A4_H', KI: 32, kinact: 0.060, system: 'recombinant CYP3A4', probe: 'BFC O-debenzylation', species: 'human', src: 'partial', pmid: '15155549',
      note: 'The same paper reports KI 166 uM / kinact 0.08 per min with testosterone 6-beta-hydroxylation as the probe. Probe-dependent; the BFC value is used here.' },
    { cid: 'silybin',    protein: 'CYP2C9_H', KI: 5,  kinact: 0.140, system: 'purified CYP2C9', probe: '7-ethoxy-4-trifluoromethylcoumarin O-deethylation', species: 'human', src: 'partial', pmid: '15155549', note: 'Purified recombinant enzyme, not liver tissue. This is the likely source of the tolbutamide false positive (point 69). The abstract confirms inactivation but not these values.' }
  ];

  // Induction. emax as fold over vehicle, ec50 in µM.
  const INDUCTION = [
    { cid: 'hyperforin', protein: 'CYP3A4_H',  emax: 7, ec50: 0.30, system: 'primary human hepatocytes', readout: 'activity', src: 'partial', pmid: '16203790',
      ec50_basis: 'ESTIMATED. Induction was dose-dependent over 0.1-1.5 uM reaching 2.6- to 7-fold; the paper reports no EC50, so 0.30 uM is read off that range and is not a published value.' },
    { cid: 'hyperforin', protein: 'CYP3A4_G',  emax: 8,  ec50: 0.30, system: 'LS180',                     readout: 'mRNA', src: 'illustrative' },
    { cid: 'hyperforin', protein: 'CYP2C9_H',  emax: 4,  ec50: 0.50, system: 'primary human hepatocytes', readout: 'mRNA, protein and activity', src: 'partial', pmid: '15100169',
      note: 'Komoroski 2004 confirms hyperforin induces CYP2C9 but publishes no Emax or EC50; both numbers here are estimates.' },
    { cid: 'hyperforin', protein: 'CYP2C19_H', emax: 4,  ec50: 0.50, system: 'primary human hepatocytes', readout: 'mRNA', src: 'illustrative' },
    { cid: 'hyperforin', protein: 'CYP2B6_H',  emax: 6,  ec50: 0.40, system: 'primary human hepatocytes', readout: 'mRNA', src: 'illustrative' },
    { cid: 'hyperforin', protein: 'ABCB1_G',   emax: 3,  ec50: 0.60, system: 'LS180',                     readout: 'mRNA', src: 'illustrative' }
  ];

  // Human crossover DDI studies. Deliberately empty: HIGH confidence is
  // unreachable until these are curated with PMIDs and extracted AUC ratios.
  const CLINICAL = [
    // Kawaguchi-Suzuki 2014, PMID 25028567. Nine healthy volunteers, four probe
    // drugs before and after a standardised milk thistle extract. No significant
    // change in any pathway. This is a STRATUM S2 HARD NEGATIVE: in vitro signal
    // exists (Sridar 2004) yet the clinical effect is null, which is exactly the
    // case the benchmark must contain to mean anything.
    { prep: 'milkthistle-std', drug: 'caffeine',        design: 'open-label, before/after', n: 9, aucRatio: 1.00, outcome: 'no significant change', stratum: 'S2', pmid: '25028567' },
    { prep: 'milkthistle-std', drug: 'tolbutamide',     design: 'open-label, before/after', n: 9, aucRatio: 1.00, outcome: 'no significant change', stratum: 'S2', pmid: '25028567' },
    { prep: 'milkthistle-std', drug: 'dextromethorphan',design: 'open-label, before/after', n: 9, aucRatio: 1.00, outcome: 'no significant change', stratum: 'S2', pmid: '25028567' },
    { prep: 'milkthistle-std', drug: 'midazolam',       design: 'open-label, before/after', n: 9, aucRatio: 1.00, outcome: 'no significant change', stratum: 'S2', pmid: '25028567' },
    // Markowitz 2003, JAMA. 14-day St John's wort course, healthy volunteers.
    // A two-fold DECREASE in alprazolam AUC, i.e. AUCR about 0.5. STRATUM S1.
    { prep: 'stjohnswort-std', drug: 'alprazolam',      design: 'before/after, healthy volunteers', n: 12, aucRatio: 0.50, outcome: 'CYP3A4 induction, alprazolam AUC halved', stratum: 'S1', doi: '10.1001/jama.290.11.1500', note: 'Markowitz et al., JAMA 2003;290:1500-1504. Abstract reports a 2-fold decrease in alprazolam AUC. The PMID previously recorded (14722142) could not be confirmed.' },
    // Silibinin proof-of-concept (Gurley et al., CPT Pharmacometrics Syst Pharmacol 2014;3:e107, PMID 24670388).
    // High-dose silibinin 1,650 mg/day. Probe AUCs rose 9% (midazolam) and 13% (S-warfarin). This is a
    // DOSE MISMATCH: the model uses the standardised extract dose, so the prediction is not like-for-like.
    { prep: 'milkthistle-std', drug: 'midazolam', dose: 1650, design: 'proof-of-concept, healthy volunteers', n: null, aucRatio: 1.09, outcome: 'AUC +9% (silibinin 1,650 mg/day)', stratum: 'S2', pmid: '24670388', note: 'n not extracted from the abstract.' },
    { prep: 'milkthistle-std', drug: 'warfarin', dose: 1650, design: 'proof-of-concept, healthy volunteers', n: null, aucRatio: 1.13, outcome: 'S-warfarin AUC +13% (silibinin 1,650 mg/day)', stratum: 'S2', pmid: '24670388', note: 'Drug is warfarin (S-enantiomer effect). n not extracted from the abstract.' }
  ];

  // Official FDA classifications (public domain). Separate from numeric PK parameters.
  const FDA_SUBSTRATES = {
    "midazolam": {
      "classification": "Sensitive index substrate (clinical)",
      "fda_enzyme": "CYP3A",
      "in_vitro_marker": "CYP3A4/5",
      "evidence": "clinical"
    },
    "warfarin": {
      "classification": "Moderately sensitive index substrate, listed as S-warfarin",
      "fda_enzyme": "CYP2C9",
      "in_vitro_marker": "S-warfarin 7-hydroxylation",
      "evidence": "clinical"
    },
    "digoxin": {
      "classification": "In vitro transporter substrate",
      "fda_enzyme": "P-gp (ABCB1)",
      "evidence": "in_vitro"
    },
    "omeprazole": {
      "classification": "Sensitive index substrate (clinical)",
      "fda_enzyme": "CYP2C19",
      "in_vitro_marker": "CYP1A2 inducer",
      "evidence": "clinical"
    },
    "dextromethorphan": {
      "classification": "Sensitive index substrate (clinical)",
      "fda_enzyme": "CYP2D6",
      "in_vitro_marker": "CYP2D6",
      "evidence": "clinical"
    },
    "theophylline": {
      "classification": "Not listed",
      "fda_enzyme": null,
      "evidence": null
    },
    "caffeine": {
      "classification": "Sensitive index substrate (clinical)",
      "fda_enzyme": "CYP1A2",
      "evidence": "clinical"
    },
    "simvastatin": {
      "classification": "Not listed on this table",
      "fda_enzyme": null,
      "evidence": null,
      "note": "appears under CYP3A and OATP1B1/3 substrates on the FDA examples page"
    },
    "atorvastatin": {
      "classification": "Not listed on this table",
      "fda_enzyme": null,
      "evidence": null,
      "note": "appears under CYP3A and OATP1B1/3 substrates on the FDA examples page"
    },
    "rosuvastatin": {
      "classification": "In vitro transporter substrate",
      "fda_enzyme": "BCRP (ABCG2); OATP1B1/OATP1B3",
      "evidence": "in_vitro"
    },
    "clopidogrel": {
      "classification": "Not a substrate; in vitro selective CYP2B6 inhibitor and clinical moderate CYP2C8 index inhibitor",
      "fda_enzyme": "CYP2B6, CYP2C8 (as inhibitor)",
      "evidence": "inhibitor"
    },
    "tacrolimus": {
      "classification": "Not listed on this table",
      "fda_enzyme": null,
      "evidence": null
    },
    "cyclosporine": {
      "classification": "Not a substrate; in vitro transporter inhibitor",
      "fda_enzyme": "P-gp; OATP1B1/OATP1B3 (as inhibitor)",
      "evidence": "inhibitor"
    },
    "phenytoin": {
      "classification": "Not a substrate; clinical strong index inducer",
      "fda_enzyme": "CYP3A (as inducer)",
      "evidence": "inducer"
    },
    "alprazolam": {
      "classification": "Not listed on this table",
      "fda_enzyme": null,
      "evidence": null,
      "note": "appears under CYP3A substrates on the FDA examples page"
    },
    "tolbutamide": {
      "classification": "Moderately sensitive index substrate (clinical)",
      "fda_enzyme": "CYP2C9",
      "evidence": "clinical"
    },
    "bupropion": {
      "classification": "In vitro marker reaction only",
      "fda_enzyme": "CYP2B6 (bupropion hydroxylation)",
      "evidence": "in_vitro"
    },
    "raltegravir": {
      "classification": "Not listed on this table",
      "fda_enzyme": null,
      "evidence": null
    },
    "metformin": {
      "classification": "In vitro transporter substrate",
      "fda_enzyme": "MATE1/MATE2-K and OCT2",
      "evidence": "in_vitro"
    },
    "lithium": {
      "classification": "Not listed on this table",
      "fda_enzyme": null,
      "evidence": null
    }
  };

  const DATA = {
    version: 'seed-0.1.0',
    PROTEINS, CONSTITUENTS, PREPARATIONS, DRUGS, FDA_SUBSTRATES,
    REVERSIBLE, TDI, INDUCTION, CLINICAL
  };

  root.HDI_DATA = DATA;
  if (typeof module !== 'undefined' && module.exports) module.exports = DATA;
})(typeof globalThis !== 'undefined' ? globalThis : this);

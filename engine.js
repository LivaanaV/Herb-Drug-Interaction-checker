/* HDI Engine v2 — pure prediction functions.
 *
 * No I/O, no DOM, no database. Takes the seed dataset plus a query, returns a
 * prediction and the audit record that reproduces it. This is what lets the
 * validation harness run the engine without an application around it.
 *
 * Model: FDA 2020 / ICH M12 basic static screen, then the mechanistic static
 * model. Transporter edges are reported as co-signals and are NOT folded into
 * the AUCR, because in vitro efflux and uptake inhibition translate poorly.
 */
(function (root) {
  'use strict';

  const ENGINE_VERSION = '2.0.0';

  // Physiological constants for the hepatic inlet term.
  const QH = 97;      // hepatic blood flow, L/h
  const KA = 1.8;     // first-order absorption rate, h^-1
  const RB = 1.0;     // blood-to-plasma ratio
  const GUT_VOLUME = 0.25; // L, the regulatory basic-model convention
  // Bile salts, food and the extract matrix raise the apparent solubility of a
  // lipophilic compound in intestinal fluid well above its value in water.
  // 10x is a conventional screening allowance, not a measurement.
  const BILE_SOLUBILISATION = 10;
  const IND_D = 1.0;  // induction scaling factor

  // Guideline cut-offs for the basic static screen.
  const CUTOFF = { hepatic: 1.02, gut: 11, induction: 0.8 };

  // Provenance vocabulary, matching the data's `src` field. Lower rank is weaker.
  // A prediction is only as strong as its weakest contributing parameter (point 48).
  // 'uncited' is a number with no source behind it. It ranks with 'illustrative'.
  const SRC_RANK = { sourced: 3, consensus: 2, partial: 1, uncited: 0, illustrative: 0 };
  function weakestSource(list) {
    if (!list.length) return 'illustrative';
    return list.reduce((w, s) => (SRC_RANK[s] < SRC_RANK[w] ? s : w), 'sourced');
  }

  // ------------------------------------------------------------------ helpers
  function pctToUmolInGut(doseMg, pct, mw) {
    const mg = doseMg * (pct / 100);
    const mmol = mg / mw;                    // mg / (g/mol) = mmol
    return (mmol / GUT_VOLUME) * 1000;       // mM -> µM
  }
  function doseUmol(doseMg, pct, mw) {
    return (doseMg * (pct / 100) / mw) * 1000;
  }

  // --------------------------------------------------------- 1. resolution
  function resolve(data, prepId, drugId) {
    const prep = data.PREPARATIONS.find(p => p.id === prepId);
    const drug = data.DRUGS.find(d => d.id === drugId);
    if (!prep || !drug) throw new Error('unresolved query');
    return { prep, drug };
  }

  // ------------------------------------------------- 2. the abundance gate
  function gate(data, prep) {
    const admitted = [], excluded = [];
    prep.constituents.forEach(link => {
      const c = data.CONSTITUENTS[link.cid];
      const quantified = typeof link.pct === 'number' && link.pct > 0;
      if (quantified || link.marker) {
        admitted.push({ link, c });
      } else {
        excluded.push({ link, c, reason: 'no quantified abundance and not a pharmacopoeial marker' });
      }
    });
    return { admitted, excluded };
  }

  // ------------------------------------------------- 3. exposure tiering
  // Returns per constituent: tier, gut concentration, and a hepatic unbound
  // interval. Tier A needs a curated faFg; Cmax alone gives an interval (B).
  function exposureFor(entry, prep) {
    const { link, c } = entry;
    // dose/250 mL assumes the whole dose dissolves. For a poorly soluble
    // compound that is wrong by orders of magnitude, and since the gut term
    // dominates every prediction, an uncapped value drives false positives.
    const igutNominal = pctToUmolInGut(prep.dose, link.pct, c.mw);
    const solubilityCap = typeof c.solubility_uM === 'number'
      ? c.solubility_uM * BILE_SOLUBILISATION : null;
    const igut = solubilityCap === null ? igutNominal : Math.min(igutNominal, solubilityCap);
    const solubilityLimited = solubilityCap !== null && igutNominal > solubilityCap;
    const dUmol = doseUmol(prep.dose, link.pct, c.mw);
    const inletMax = (KA * dUmol) / (QH * RB);   // Fa·Fg = 1 upper bound, µM

    if (!c.pk) {
      return { tier: 'C', igut, igutNominal, solubilityLimited, ihLow: 0, ihHigh: 0, hepatic: false,
               reason: 'No constituent pharmacokinetics on file. Gut lumen concentration is computable from the dose; hepatic exposure is not.' };
    }
    // Cmax is measured at one dose. Scale it linearly to the dose actually
    // queried, using the constituent mass rather than the preparation mass.
    // Linearity is an assumption; it is recorded in the audit.
    const constituentMg = prep.dose * (link.pct / 100);
    const scale = c.pk.dose ? constituentMg / c.pk.dose : 1;
    const fu = c.pk.fu, cmax = c.pk.cmax * scale;
    if (typeof c.pk.faFg === 'number') {
      const ih = fu * (cmax + (c.pk.faFg * KA * dUmol) / (QH * RB));
      return { tier: 'A', igut, igutNominal, solubilityLimited, ihLow: ih, ihHigh: ih, hepatic: true, doseScale: scale,
               reason: 'Human pharmacokinetics with a curated absorbed fraction, so the hepatic inlet concentration is a point estimate.' };
    }
    return { tier: 'B', igut, igutNominal, solubilityLimited, ihLow: fu * cmax, ihHigh: fu * (cmax + inletMax), hepatic: true, doseScale: scale,
             reason: 'Human Cmax is known but the absorbed fraction is not, so hepatic exposure is an interval between systemic unbound and a fully-absorbed inlet bound.' };
  }

  // --------------------------------------- 4. typed bioactivity retrieval
  // Time-dependent inhibition is used only when measured in a hepatic system.
  // Recombinant and purified-enzyme inactivation constants are not transferable
  // to hepatic prediction (point 4): the same paper reports very different KI for
  // different probes and enzyme sources, and the inactivation rate depends on the
  // system. Excluded records are listed in the audit.
  const HEPATIC_SYSTEM = /HLM|microsom|hepatocyte|liver/i;
  function edgesFor(data, cids) {
    const set = new Set(cids);
    return {
      rev: data.REVERSIBLE.filter(r => set.has(r.cid)),
      tdi: data.TDI.filter(r => set.has(r.cid) && HEPATIC_SYSTEM.test(r.system || '')),
      ind: data.INDUCTION.filter(r => set.has(r.cid))
    };
  }

  // ----------------------------- 5/6. compartment terms and the static model
  // conc: map cid -> concentration in µM for THIS compartment.
  function compartmentTerms(data, proteinId, edges, conc) {
    const protein = data.PROTEINS[proteinId];
    let sumRatio = 0, sumInact = 0, bestInd = 0;
    const used = [];

    edges.rev.filter(e => e.protein === proteinId).forEach(e => {
      const I = conc[e.cid] || 0;
      if (I <= 0) return;
      const ratio = I / e.ki;
      sumRatio += ratio;
      used.push({ kind: 'reversible', cid: e.cid, I, param: e.ki, paramName: 'Ki', ratio, system: e.system, probe: e.probe, src: e.src, pmid: e.pmid, kiBasis: e.ki_basis });
    });

    edges.tdi.filter(e => e.protein === proteinId).forEach(e => {
      const I = conc[e.cid] || 0;
      if (I <= 0) return;
      const kinactH = e.kinact * 60;
      const contrib = (kinactH * I) / (e.KI + I);
      sumInact += contrib;
      used.push({ kind: 'tdi', cid: e.cid, I, param: e.KI, paramName: 'KI', ratio: contrib, system: e.system, src: e.src, pmid: e.pmid });
    });

    edges.ind.filter(e => e.protein === proteinId).forEach(e => {
      const I = conc[e.cid] || 0;
      if (I <= 0) return;
      const term = (IND_D * e.emax * I) / (e.ec50 + I);
      if (term > bestInd) bestInd = term;
      used.push({ kind: 'induction', cid: e.cid, I, param: e.ec50, paramName: 'EC50', ratio: term, system: e.system, src: e.src, pmid: e.pmid });
    });

    const kdeg = protein && protein.kdeg ? protein.kdeg : null;
    const A = 1 / (1 + sumRatio);
    const B = kdeg === null ? 1 : kdeg / (kdeg + sumInact);
    const C = 1 + bestInd;
    return { A, B, C, product: A * B * C, used, sumRatio, sumInact, bestInd };
  }

  function screenPasses(terms, compartment) {
    const r1 = 1 + terms.sumRatio;
    const r3 = 1 / (1 + terms.bestInd);
    const inhibitionCut = compartment === 'gut' ? CUTOFF.gut : CUTOFF.hepatic;
    return {
      r1, r3,
      inhibition: r1 >= inhibitionCut || terms.sumInact > 0,
      induction: r3 <= CUTOFF.induction,
      pass: r1 >= inhibitionCut || terms.sumInact > 0 || r3 <= CUTOFF.induction
    };
  }

  function ratioForMoiety(data, moiety, edges, concH, concG, isProdrug) {
    const detail = { hepatic: [], gut: null };

    if (isProdrug) {
      let amr = 0, accounted = 0;
      (moiety.activation || []).forEach(path => {
        const t = compartmentTerms(data, path.protein, edges, concH);
        amr += t.product * path.fAct;
        accounted += path.fAct;
        detail.hepatic.push({ protein: path.protein, fm: path.fAct, terms: t, role: 'activation' });
      });
      amr += (1 - accounted);
      return { ratio: amr, kind: 'AMR', detail };
    }

    let hepatic = 1;
    (moiety.clearance || []).forEach(cl => {
      const t = compartmentTerms(data, cl.protein, edges, concH);
      hepatic *= (t.product * cl.fm + (1 - cl.fm));
      detail.hepatic.push({ protein: cl.protein, fm: cl.fm, terms: t, role: 'clearance' });
    });

    let gut = 1;
    if (moiety.gutEnzyme) {
      const t = compartmentTerms(data, moiety.gutEnzyme, edges, concG);
      gut = t.product * (1 - moiety.fg) + moiety.fg;
      detail.gut = { protein: moiety.gutEnzyme, fg: moiety.fg, terms: t };
    }

    const denom = hepatic * gut;
    return { ratio: denom === 0 ? Infinity : 1 / denom, kind: 'AUCR', detail, hepaticBracket: hepatic, gutBracket: gut };
  }

  // ------------------------------------------------- 7. severity, confidence
  function severityFor(ratio, moiety, drug) {
    const nti = !!moiety.nti;
    const critical = nti || !!drug.criticalEfficacy;
    if (drug.isProdrug) {
      if (ratio <= 0.5) return critical ? 'MAJOR' : 'MODERATE';
      if (ratio <= 0.8) return critical ? 'MODERATE' : 'MINOR';
      return 'NONE';
    }
    if (ratio >= 5) return 'MAJOR';
    if (ratio >= 2) return nti ? 'MAJOR' : 'MODERATE';
    if (ratio >= 1.25) return nti ? 'MODERATE' : 'MINOR';
    if (ratio <= 0.2) return critical ? 'MAJOR' : 'MODERATE';
    if (ratio <= 0.5) return 'MODERATE';
    if (ratio <= 0.8) return critical ? 'MODERATE' : 'MINOR';
    return 'NONE';
  }

  function directionFor(ratio, drug) {
    if (drug.isProdrug) return ratio < 0.95 ? 'reduced activation' : 'no meaningful change';
    if (ratio >= 1.25) return 'accumulation';
    if (ratio <= 0.8) return 'loss of exposure';
    return 'no meaningful change';
  }

  // Confidence levels (point 49):
  //   HIGH     exposure tier A, and every contributing parameter is sourced with a PMID.
  //            Not reachable at present: fu is uncited for every constituent.
  //   MODERATE exposure tier A, but some contributing parameter is consensus-only.
  //   LOW      anything else: tier B/C/D, an illustrative or partial parameter, a
  //            transporter-only or non-human signal, or a discordant clinical study.
  // A clinical study does NOT raise confidence. The validation corpus was seen while
  // the parameters were being curated, so it cannot also confirm them (points 43-44).
  // Clinical agreement is shown beside the prediction instead.
  function confidenceFor(ctx) {
    let level = ctx.tier === 'A' ? 'MODERATE' : 'LOW';
    const caps = [];
    if (ctx.tier === 'A' && ctx.weakestSrc === 'sourced') level = 'HIGH';
    if (SRC_RANK[ctx.weakestSrc] < 2) {
      level = 'LOW'; caps.push('A contributing parameter is uncited or partly cited (provenance: ' + ctx.weakestSrc + ').');
    }
    if (ctx.uncappedGut) {
      level = 'LOW'; caps.push('A contributing constituent has no solubility on record, so its gut concentration is not capped and may be overstated.');
    }
    if (ctx.clinicalCheck && ctx.clinicalCheck.status === 'discordant') {
      level = 'LOW'; caps.push('A clinical study of this preparation and drug disagrees with the predicted direction or size.');
    }
    if (ctx.isProdrug && level === 'HIGH') { level = 'MODERATE'; caps.push('Prodrug activation is formation-rate-limited; a static model cannot capture it.'); }
    if (ctx.transporterOnly) { level = 'LOW'; caps.push('Transporter-mediated signal only. In vitro efflux and uptake inhibition predict clinical exposure change poorly.'); }
    if (ctx.nonHuman) { level = 'LOW'; caps.push('A non-human in vitro system contributed to this prediction.'); }
    if (ctx.provisional && level !== 'LOW') { level = 'LOW'; caps.push('At least one contributing parameter is an uncurated placeholder.'); }
    return { level, caps };
  }

  // Compares the governing prediction with any clinical study of the same
  // preparation and drug. Within twofold and on the same side of the no-effect
  // band counts as concordant. The result is displayed, not used to raise confidence.
  function clinicalCheckFor(studies, predicted) {
    if (!studies.length) return { status: 'none' };
    const s = studies[0];
    const obsInt = s.aucRatio <= 0.8 || s.aucRatio >= 1.25;
    const predInt = predicted <= 0.8 || predicted >= 1.25;
    const fold = Math.max(predicted, s.aucRatio) / Math.min(predicted, s.aucRatio);
    const concordant = obsInt === predInt && (!obsInt || fold <= 2);
    return { status: concordant ? 'concordant' : 'discordant', observed: s.aucRatio, predicted, fold, pmid: s.pmid, stratum: s.stratum };
  }

  // ----------------------------------------------------------- orchestration
  function predict(data, prepId, drugId, opts) {
    opts = opts || {};
    const { prep, drug } = resolve(data, prepId, drugId);
    const dose = typeof opts.dose === 'number' && opts.dose > 0 ? opts.dose : prep.dose;
    const prepAtDose = Object.assign({}, prep, { dose });

    const gated = gate(data, prepAtDose);
    const exposures = gated.admitted.map(e => Object.assign({ cid: e.link.cid, name: e.c.name }, exposureFor(e, prepAtDose), { link: e.link, c: e.c }));
    const cids = exposures.map(e => e.cid);
    const edges = edgesFor(data, cids);

    const activeCids = new Set([].concat(
      edges.rev.map(e => e.cid), edges.tdi.map(e => e.cid), edges.ind.map(e => e.cid)));
    const contributing = exposures.filter(e => activeCids.has(e.cid));

    const tierRank = { A: 4, B: 3, C: 2, D: 1 };
    const bestTier = contributing.reduce((best, e) => tierRank[e.tier] > tierRank[best] ? e.tier : best, 'D');

    // Provisional is reported with its reasons, so curation has a worklist rather
    // than a yes/no. Only inputs that actually feed the reported ratio count:
    // an unsourced transporter edge is a co-signal caveat, not a reason to
    // distrust the number.
    const ratioProteins = new Set();
    drug.moieties.forEach(m => {
      (m.clearance || []).forEach(c => ratioProteins.add(c.protein));
      (m.activation || []).forEach(a => ratioProteins.add(a.protein));
      if (m.gutEnzyme) ratioProteins.add(m.gutEnzyme);
    });
    const provisionalReasons = [];
    [].concat(edges.rev, edges.tdi, edges.ind).forEach(e => {
      if (activeCids.has(e.cid) && ratioProteins.has(e.protein) && e.src === 'illustrative')
        provisionalReasons.push({ kind: 'kinetics', cid: e.cid, protein: e.protein,
          what: (data.CONSTITUENTS[e.cid] || {}).name + ' against ' + (data.PROTEINS[e.protein] || {}).symbol });
    });
    contributing.forEach(e => {
      if (e.c.pk && e.c.pk.src === 'illustrative')
        provisionalReasons.push({ kind: 'exposure', cid: e.cid,
          what: e.name + ' pharmacokinetics (Cmax, fu, absorbed fraction)' });
      if (!e.c.pk)
        provisionalReasons.push({ kind: 'exposure-missing', cid: e.cid,
          what: e.name + ' has no human pharmacokinetics on file' });
      if (e.link.src === 'illustrative')
        provisionalReasons.push({ kind: 'abundance', cid: e.cid,
          what: e.name + ' abundance in this preparation' });
    });
    if (prepAtDose.doseSrc === 'illustrative')
      provisionalReasons.push({ kind: 'dose', what: 'typical daily dose of this preparation' });
    const provisional = provisionalReasons.length > 0;

    const audit = {
      engineVersion: ENGINE_VERSION, dataVersion: data.version,
      query: { preparation: prep.id, dose, drug: drug.id },
      gate: { admitted: cids, excluded: gated.excluded.map(e => e.link.cid) },
      exposure: contributing.map(e => ({ cid: e.cid, tier: e.tier, igut_uM: e.igut, ih_uM: [e.ihLow, e.ihHigh] })),
      constants: { QH, KA, RB, GUT_VOLUME, cutoffs: CUTOFF }
    };

    // transporter co-signals, never folded into the ratio
    const transporterSignals = [];
    drug.moieties.forEach(m => (m.transporters || []).forEach(tp => {
      const rel = edges.rev.filter(e => e.protein === tp).concat(edges.ind.filter(e => e.protein === tp));
      if (!rel.length) return;
      const organ = data.PROTEINS[tp].organ;
      rel.forEach(e => {
        const ex = contributing.find(x => x.cid === e.cid);
        if (!ex) return;
        const I = organ === 'intestinal' ? ex.igut : ex.ihHigh;
        const r1 = 1 + I / (e.ki || e.ec50);
        const cut = organ === 'intestinal' ? CUTOFF.gut : CUTOFF.hepatic;
        if (r1 >= cut) transporterSignals.push({ protein: tp, symbol: data.PROTEINS[tp].symbol, cid: e.cid, name: data.CONSTITUENTS[e.cid].name, r1, system: e.system, src: e.src });
      });
    }));

    // A drug with no panel pathway at all is the negative-control case. Say so
    // explicitly rather than reporting a ratio of exactly 1, which reads as a
    // computed finding when nothing was computed.
    const hasPathway = drug.moieties.some(m =>
      (m.clearance || []).length || (m.activation || []).length || m.gutEnzyme || (m.transporters || []).length);
    if (!hasPathway) {
      return {
        outcome: 'no-pathway', engineVersion: ENGINE_VERSION, prep: prepAtDose, drug,
        tier: null, exposures, gated, edges, transporterSignals: [], provisional, provisionalReasons, audit,
        message: drug.inn + ' is cleared by a route outside the modelled panel, so no constituent of this preparation can alter its exposure through any mechanism this engine represents.'
      };
    }

    if (bestTier === 'D' || contributing.length === 0) {
      return {
        outcome: 'hypothesis', engineVersion: ENGINE_VERSION, prep: prepAtDose, drug,
        tier: 'D', exposures, gated, edges, transporterSignals, provisional, provisionalReasons, audit,
        message: contributing.length === 0
          ? 'No constituent passing the gate has a recorded interaction with any protein this drug depends on.'
          : 'Exposure for every contributing constituent is unanchored, so no exposure ratio is computed.'
      };
    }

    const results = drug.moieties.map(moiety => {
      const concGut = {}, concHLow = {}, concHHigh = {};
      contributing.forEach(e => {
        concGut[e.cid] = e.igut;
        concHLow[e.cid] = e.hepatic ? e.ihLow : 0;
        concHHigh[e.cid] = e.hepatic ? e.ihHigh : 0;
      });

      const low  = ratioForMoiety(data, moiety, edges, concHLow,  concGut, drug.isProdrug);
      const high = ratioForMoiety(data, moiety, edges, concHHigh, concGut, drug.isProdrug);

      // For accumulation the worst case is the larger ratio; for loss of
      // exposure and prodrug activation it is the smaller one.
      const candidates = [low.ratio, high.ratio];
      const worst = drug.isProdrug ? Math.min.apply(null, candidates)
        : (Math.max.apply(null, candidates) >= 1 ? Math.max.apply(null, candidates) : Math.min.apply(null, candidates));
      const other = candidates.find(v => v !== worst);

      // Point 5: if the two exposure ends land in different severity bands, the
      // severity is not determined by the data. It is reported as INDETERMINATE
      // rather than as the worst end.
      const sevEnds = [severityFor(low.ratio, moiety, drug), severityFor(high.ratio, moiety, drug)];
      const screenGut = low.detail.gut ? screenPasses(low.detail.gut.terms, 'gut') : null;
      const screenHep = low.detail.hepatic.map(h => ({ protein: h.protein, fm: h.fm, role: h.role, screen: screenPasses(high.detail.hepatic.find(x => x.protein === h.protein).terms, 'hepatic') }));

      return {
        moiety, kind: low.kind, ratio: worst,
        interval: [Math.min(low.ratio, high.ratio), Math.max(low.ratio, high.ratio)],
        pointEstimate: bestTier === 'A',
        detailLow: low, detailHigh: high, screenGut, screenHep,
        severity: sevEnds[0] === sevEnds[1] ? sevEnds[0] : 'INDETERMINATE',
        severityRank: sevEnds.sort((x, y) => ({ NONE: 0, MINOR: 1, MODERATE: 2, MAJOR: 3 }[x] - { NONE: 0, MINOR: 1, MODERATE: 2, MAJOR: 3 }[y])).slice(-1)[0],
        severityRange: sevEnds,
        direction: directionFor(worst, drug)
      };
    });

    const order = { MAJOR: 4, MODERATE: 3, MINOR: 2, NONE: 1 };
    const governing = results.reduce((a, b) => order[b.severityRank] > order[a.severityRank] ? b : a);

    const anyEnzymeEdge = results.some(r => r.detailHigh.detail.hepatic.length || r.detailHigh.detail.gut);
    const nonHuman = [].concat(edges.rev, edges.tdi, edges.ind)
      .some(e => activeCids.has(e.cid) && e.species && e.species !== 'human');

    const clinicalCheck = clinicalCheckFor(
      data.CLINICAL.filter(s => s.prep === prep.id && s.drug === drug.id), governing.ratio);

    // Provenance of every parameter that feeds the ratio. fu is tagged illustrative
    // unless a sourced fu is recorded (point 36): the pk-level src covers Cmax and
    // Fa*Fg only.
    const contribSrc = [];
    contributing.forEach(e => {
      if (!e.c.pk) { contribSrc.push('illustrative'); return; }
      contribSrc.push(e.c.pk.fuSrc || 'illustrative');
      contribSrc.push(e.c.pk.src || 'illustrative');
    });
    [].concat(edges.rev, edges.tdi, edges.ind)
      .filter(e => activeCids.has(e.cid)).forEach(e => contribSrc.push(e.src || 'illustrative'));
    // Point 47: fm and Fg of the governing moiety are parameters too.
    if (governing.moiety) {
      if (governing.moiety.fgSrc) contribSrc.push(governing.moiety.fgSrc);
      (governing.moiety.clearance || []).forEach(c => contribSrc.push(c.src || 'uncited'));
      (governing.moiety.activation || []).forEach(a => contribSrc.push(a.src || 'uncited'));
    }
    const weakestSrc = weakestSource(contribSrc);

    const confidence = confidenceFor({
      tier: bestTier, weakestSrc, clinicalCheck,
      uncappedGut: contributing.some(e => typeof e.c.solubility_uM !== 'number' && (e.igut || 0) > 1),
      isProdrug: !!drug.isProdrug,
      transporterOnly: !anyEnzymeEdge && transporterSignals.length > 0,
      nonHuman, provisional
    });

    audit.result = { kind: governing.kind, ratio: governing.ratio, severity: governing.severity, confidence: confidence.level, tier: bestTier };

    return {
      outcome: 'scored', engineVersion: ENGINE_VERSION, prep: prepAtDose, drug,
      tier: bestTier, exposures, contributing, gated, edges,
      results, governing, confidence, clinicalCheck, weakestSrc, transporterSignals, provisional, provisionalReasons, audit
    };
  }

  const ENGINE = { ENGINE_VERSION, BILE_SOLUBILISATION, predict, gate, exposureFor, compartmentTerms, ratioForMoiety, severityFor, confidenceFor, pctToUmolInGut, CUTOFF };
  root.HDI_ENGINE = ENGINE;
  if (typeof module !== 'undefined' && module.exports) module.exports = ENGINE;
})(typeof globalThis !== 'undefined' ? globalThis : this);

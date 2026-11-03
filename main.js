/* =========================================================================
 *  Buckley-Leverett — analytical solution
 *  Units: x, L [cm]; A [cm²]; q_inj [cc/min]; dt [h]; μ [cP]
 * ========================================================================= */

"use strict";

const DEFAULTS = {
  N: 130,          // points to plot
  mu_w: 1.0,       // water viscosity [cP]
  mu_o: 1.1,       // oil viscosity [cP]
  Swc: 0.0,        // connate water saturation
  Sor: 0.1,        // residual oil saturation
  krw_end: 1.0,    // krw endpoint (Sw = 1 - Sor)
  kro_end: 1.0,    // kro endpoint (Sw = Swc)
  nw: 2.0,         // Corey exponent for water
  no: 2.0,         // Corey exponent for oil
  por: 0.22,       // porosity
  A: 11.34,        // cross-section area [cm²]
  L: 18,           // core length [cm]
  q_inj: 0.25,     // [cc/min]
  dt: 0.000625,    // [h]
  tempos: "1600, 4800",
};

const COLORS = ["#1f77b4", "#2ca02c", "#e6a700", "#d62728", "#17becf", "#9467bd", "#111111", "#8c564b"];

/* -------------------------------------------------------------------------
 *  Relative permeability (Corey) and its derivatives with respect to Sw
 * ------------------------------------------------------------------------- */
function krel(Sw, p) {
  const D = 1 - p.Swc - p.Sor;
  const SeRaw = (Sw - p.Swc) / D;
  const Se = Math.max(0, Math.min(1, SeRaw));
  const inside = SeRaw > 0 && SeRaw < 1;

  const krw = p.krw_end * Math.pow(Se, p.nw);
  const kro = p.kro_end * Math.pow(1 - Se, p.no);
  const dkrw = inside ? p.krw_end * p.nw * Math.pow(Se, p.nw - 1) / D : 0;
  const dkro = inside ? -p.kro_end * p.no * Math.pow(1 - Se, p.no - 1) / D : 0;

  return { krw, kro, dkrw, dkro };
}

/* -------------------------------------------------------------------------
 *  Fractional flow fw = λw / (λw + λo) and analytical derivative dfw/dSw
 * ------------------------------------------------------------------------- */
function fracFlow(Sw, p) {
  const k = krel(Sw, p);
  const lw = k.krw / p.mu_w;
  const lo = k.kro / p.mu_o;
  const dlw = k.dkrw / p.mu_w;
  const dlo = k.dkro / p.mu_o;
  const lt = lw + lo;
  if (lt <= 0) return { fw: 0, dfw: 0 };
  const fw = lw / lt;
  const dfw = (dlw * lo - lw * dlo) / (lt * lt);
  return { fw, dfw: Number.isFinite(dfw) ? dfw : 0 };
}

/* -------------------------------------------------------------------------
 *  Welge tangent: Sf maximizes the secant fw(S)/(S - Swi), Swi = Swc.
 *  At the optimum, fw(Sf)/(Sf - Swi) = dfw/dSw(Sf).
 * ------------------------------------------------------------------------- */
function welge(p) {
  const Swi = p.Swc;
  const Smax = 1 - p.Sor;
  const secant = (S) => fracFlow(S, p).fw / (S - Swi);

  // coarse search
  const M = 4000;
  const h = (Smax - Swi) / M;
  let iBest = 1;
  let best = -Infinity;
  for (let i = 1; i <= M; i++) {
    const v = secant(Swi + i * h);
    if (v > best) { best = v; iBest = i; }
  }

  // golden-section refinement in the neighboring interval
  let a = Swi + Math.max(iBest - 1, 1e-6) * h;
  let b = Swi + Math.min(iBest + 1, M) * h;
  const gr = (Math.sqrt(5) - 1) / 2;
  let c = b - gr * (b - a);
  let d = a + gr * (b - a);
  for (let it = 0; it < 100; it++) {
    if (secant(c) > secant(d)) b = d; else a = c;
    c = b - gr * (b - a);
    d = a + gr * (b - a);
  }
  const Sf = iBest === M ? Smax : 0.5 * (a + b);
  const ff = fracFlow(Sf, p);
  const slope = ff.fw / (Sf - Swi);   // = dimensionless shock velocity

  return {
    Swi,
    Sf,
    fwSf: ff.fw,
    slope,
    SwAvgBT: Swi + 1 / slope,         // average saturation behind the front (Welge)
    pviBT: 1 / slope,                 // pore volumes injected at breakthrough
  };
}

/* -------------------------------------------------------------------------
 *  Sw(x) profile at a given time: rarefaction x(S) = xD·f'(S) for S ∈ [Sf, 1-Sor],
 *  shock at xf = xD·f'(Sf) dropping to Swi, then Swi up to L.
 * ------------------------------------------------------------------------- */
function profile(p, w, Vinj, nPts) {
  const xD = Vinj / (p.por * p.A);
  const Smax = 1 - p.Sor;
  const xs = [0];
  const ss = [Smax];

  for (let i = 0; i <= nPts; i++) {
    const S = Smax - (Smax - w.Sf) * (i / nPts);
    const dfw = S === w.Sf ? w.slope : fracFlow(S, p).dfw;
    xs.push(xD * dfw);
    ss.push(S);
  }
  const xf = xD * w.slope;
  xs.push(xf, Math.max(xf, p.L));
  ss.push(w.Swi, w.Swi);

  return { ...clipAtL(xs, ss, p.L), xf };
}

function clipAtL(xs, ss, L) {
  const x = [xs[0]];
  const s = [ss[0]];
  for (let i = 1; i < xs.length; i++) {
    if (xs[i] <= L) {
      x.push(xs[i]); s.push(ss[i]);
      continue;
    }
    const x0 = xs[i - 1], s0 = ss[i - 1];
    if (x0 < L) {
      const t = (L - x0) / (xs[i] - x0);
      x.push(L); s.push(s0 + t * (ss[i] - s0));
    }
    break;
  }
  return { x, s };
}

/* -------------------------------------------------------------------------
 *  Full computation
 * ------------------------------------------------------------------------- */
function linspace(a, b, n) {
  if (n < 2) return [a];
  return Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
}

function solve(p) {
  const Sw = linspace(0, 1, p.N);
  const kr = Sw.map((S) => krel(S, p));
  const ff = Sw.map((S) => fracFlow(S, p));
  const w = welge(p);
  const Vp = p.por * p.A * p.L;

  const times = p.tempos.map((t, k) => {
    const Vinj = p.q_inj * t * p.dt * 60;   // [cc]
    const xD = Vinj / (p.por * p.A);
    return {
      steps: t,
      hours: t * p.dt,
      Vinj,
      pvi: Vinj / Vp,
      color: COLORS[k % COLORS.length],
      xRaw: ff.map((f) => xD * f.dfw),       // multivalued BL curve
      prof: profile(p, w, Vinj, Math.max(p.N, 200)),
    };
  });

  const tBTmin = w.pviBT * Vp / p.q_inj;
  return { Sw, kr, ff, w, Vp, times, tBT: { min: tBTmin, h: tBTmin / 60, steps: tBTmin / 60 / p.dt } };
}

/* =========================================================================
 *  User interface
 * ========================================================================= */
const NUM_FIELDS = Object.keys(DEFAULTS).filter((k) => k !== "tempos");

function readInputs() {
  const errors = [];
  const p = {};
  document.querySelectorAll("input.invalid").forEach((el) => el.classList.remove("invalid"));
  const bad = (id, msg) => { document.getElementById(id).classList.add("invalid"); errors.push(msg); };

  for (const id of NUM_FIELDS) {
    const v = parseFloat(document.getElementById(id).value);
    if (!Number.isFinite(v)) bad(id, `${id}: invalid value`);
    p[id] = v;
  }
  p.N = Math.round(p.N);
  p.tempos = document.getElementById("tempos").value
    .split(/[,;\s]+/).filter(Boolean).map(Number);

  if (!p.tempos.length || p.tempos.some((t) => !Number.isFinite(t) || t < 0)) bad("tempos", "times: invalid list");
  if (p.N < 10) bad("N", "N must be ≥ 10");
  for (const id of ["mu_w", "mu_o", "krw_end", "kro_end", "nw", "no", "por", "A", "L", "q_inj", "dt"]) {
    if (p[id] <= 0) bad(id, `${id} must be > 0`);
  }
  if (p.por > 1) bad("por", "porosity must be ≤ 1");
  if (p.Swc < 0) bad("Swc", "Swc must be ≥ 0");
  if (p.Sor < 0) bad("Sor", "Sor must be ≥ 0");
  if (p.Swc + p.Sor >= 1) { bad("Swc", "Swc + Sor must be < 1"); bad("Sor", ""); }

  return { p, errors: errors.filter(Boolean) };
}

const fmt = (v, d = 4) => (Number.isFinite(v) ? Number(v.toPrecision(d)).toString() : "—");

function renderSummary(r, p) {
  const { w, tBT } = r;
  const items = [
    ["S<sub>w</sub> at the front (S<sub>wf</sub>)", fmt(w.Sf)],
    ["f<sub>w</sub>(S<sub>wf</sub>)", fmt(w.fwSf)],
    ["df<sub>w</sub>/dS<sub>w</sub> (S<sub>wf</sub>)", fmt(w.slope)],
    ["S̄<sub>w</sub> behind the front", fmt(w.SwAvgBT)],
    ["Pore volume V<sub>p</sub> [cc]", fmt(r.Vp)],
    ["PVI at breakthrough", fmt(w.pviBT)],
    ["Breakthrough time [h]", fmt(tBT.h)],
    ["Breakthrough time [steps]", fmt(tBT.steps, 6)],
  ];
  document.getElementById("summary").innerHTML =
    items.map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join("");

  document.querySelector("#times-table tbody").innerHTML = r.times.map((t) => {
    const xf = t.prof.xf > p.L ? `> L (${fmt(t.prof.xf)})` : fmt(t.prof.xf);
    return `<tr>
      <td><span class="swatch" style="background:${t.color}"></span>${t.steps}</td>
      <td>${fmt(t.hours)}</td><td>${fmt(t.Vinj)}</td><td>${fmt(t.pvi)}</td><td>${xf}</td>
    </tr>`;
  }).join("");
}

const PLOT_CONFIG = { responsive: true, displaylogo: false, toImageButtonOptions: { format: "png", scale: 2 } };

function baseLayout(title, xTitle, yTitle, extra = {}) {
  return {
    title: { text: title, font: { size: 15 } },
    margin: { l: 60, r: 20, t: 46, b: 50 },
    xaxis: { title: { text: xTitle }, showgrid: true, zeroline: false, ...(extra.xaxis || {}) },
    yaxis: { title: { text: yTitle }, showgrid: true, zeroline: false, ...(extra.yaxis || {}) },
    legend: { bgcolor: "rgba(255,255,255,0.8)", ...(extra.legend || {}) },
    hovermode: "closest",
    paper_bgcolor: "#ffffff",
    plot_bgcolor: "#ffffff",
  };
}

function line(x, y, name, color, lineStyle = {}) {
  return { x, y, name, type: "scatter", mode: "lines", line: { color, width: 2, ...lineStyle } };
}

function renderPlots(r, p) {
  const { Sw, kr, ff, w } = r;

  // Saturation profile Sw(x)
  Plotly.react("plot-profile", r.times.map((t) =>
    line(t.prof.x, t.prof.s, `t = ${t.steps}`, t.color)
  ), baseLayout("Buckley-Leverett Solution", "x [cm]", "S<sub>w</sub>", {
    xaxis: { range: [0, p.L], dtick: p.L <= 40 ? 1 : undefined },
    yaxis: { range: [0, 1] },
    legend: { x: 0.99, xanchor: "right", y: 0.02, yanchor: "bottom" },
  }), PLOT_CONFIG);

  // Relative permeability
  Plotly.react("plot-kr", [
    line(Sw, kr.map((k) => k.krw), "k<sub>rw</sub>", "#1f77b4"),
    line(Sw, kr.map((k) => k.kro), "k<sub>ro</sub>", "#d62728"),
  ], baseLayout("Relative Permeability Curve", "S<sub>w</sub>", "k<sub>r</sub>", {
    xaxis: { range: [0, 1] },
    legend: { x: 0.99, xanchor: "right", y: 0.5 },
  }), PLOT_CONFIG);

  // Fractional flow + Welge tangent
  const tanEnd = Math.min(1, w.Swi + 1 / w.slope);
  Plotly.react("plot-fw", [
    line(Sw, ff.map((f) => f.fw), "f<sub>w</sub>", "#1f77b4"),
    line([w.Swi, tanEnd], [0, (tanEnd - w.Swi) * w.slope], "Welge tangent", "#888888", { dash: "dash", width: 1.5 }),
    { x: [w.Sf], y: [w.fwSf], name: "S<sub>wf</sub>", type: "scatter", mode: "markers", marker: { color: "#d62728", size: 9 } },
  ], baseLayout("Fractional Flow Curve", "S<sub>w</sub>", "f<sub>w</sub>(S<sub>w</sub>)", {
    xaxis: { range: [0, 1] }, yaxis: { range: [0, 1.02] },
    legend: { x: 0.01, y: 0.99 },
  }), PLOT_CONFIG);

  // Fractional flow derivative
  Plotly.react("plot-dfw", [
    line(Sw, ff.map((f) => f.dfw), "df<sub>w</sub>/dS<sub>w</sub>", "#1f77b4"),
    { x: [w.Sf], y: [w.slope], name: "S<sub>wf</sub>", type: "scatter", mode: "markers", marker: { color: "#d62728", size: 9 } },
  ], baseLayout("Fractional Flow Curve Derivative", "S<sub>w</sub>", "df<sub>w</sub>/dS<sub>w</sub>", {
    xaxis: { range: [0, 1] },
    legend: { x: 0.01, y: 0.99 },
  }), PLOT_CONFIG);

  // Multivalued BL curve x(Sw) (before the shock correction)
  Plotly.react("plot-bl", r.times.map((t) =>
    line(t.xRaw, Sw, `t = ${t.steps}`, t.color)
  ), baseLayout("Buckley-Leverett Curve (multivalued)", "x [cm]", "S<sub>w</sub>", {
    xaxis: { range: [0, p.L] }, yaxis: { range: [0, 1] },
    legend: { x: 0.99, xanchor: "right", y: 0.02, yanchor: "bottom" },
  }), PLOT_CONFIG);
}

let lastResult = null;

function run() {
  const { p, errors } = readInputs();
  const errBox = document.getElementById("error");
  if (errors.length) {
    errBox.innerHTML = errors.join("<br>");
    errBox.hidden = false;
    return;
  }
  errBox.hidden = true;
  const r = solve(p);
  lastResult = r;
  renderSummary(r, p);
  renderPlots(r, p);
}

function exportCSV() {
  if (!lastResult) return;
  const lines = [];
  for (const t of lastResult.times) {
    lines.push(`# t = ${t.steps} steps (${t.hours} h), PVI = ${t.pvi}`);
    lines.push("x_cm,Sw");
    t.prof.x.forEach((x, i) => lines.push(`${x},${t.prof.s[i]}`));
    lines.push("");
  }
  const blob = new Blob([lines.join("\n")], { type: "text/csv" });
  const a = document.createElement("a");
  a.href = URL.createObjectURL(blob);
  a.download = "buckley_leverett_profiles.csv";
  a.click();
  URL.revokeObjectURL(a.href);
}

function resetDefaults() {
  for (const [id, v] of Object.entries(DEFAULTS)) document.getElementById(id).value = v;
  run();
}

if (typeof document !== "undefined") {
  document.addEventListener("DOMContentLoaded", () => {
    const form = document.getElementById("input-form");
    form.addEventListener("submit", (e) => { e.preventDefault(); run(); });
    let timer;
    form.addEventListener("input", () => { clearTimeout(timer); timer = setTimeout(run, 300); });
    document.getElementById("btn-reset").addEventListener("click", resetDefaults);
    document.getElementById("btn-csv").addEventListener("click", exportCSV);
    run();
  });
}

if (typeof module !== "undefined") module.exports = { krel, fracFlow, welge, profile, solve, DEFAULTS };

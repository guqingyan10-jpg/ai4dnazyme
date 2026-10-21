"use strict";

// Offline display only. Data is always inserted as text, never as HTML.
const MAX_FILE_BYTES = 50 * 1024 * 1024;
const MAX_TABLE_ROWS = 2000;
const state = { result: null, selected: null, filtered: [] };
const $ = (id) => document.getElementById(id);
const finite = (value) => typeof value === "number" && Number.isFinite(value);
const clean = (value) => value === null || value === undefined || value === "" ? "Not available" : String(value);
const number = (value, digits = 3) => finite(value) ? value.toFixed(digits) : "Not available";
const probability = (value) => finite(value) && value >= 0 && value <= 1 ? `${(100 * value).toFixed(2)}%` : "Not available";
const integer = (value) => finite(value) && Number.isInteger(value) && value >= 0 ? String(value) : "Not evaluated";
const field = (object, ...names) => {
  for (const name of names) if (object[name] !== undefined && object[name] !== null) return object[name];
  return null;
};

function node(tag, text, className) {
  const element = document.createElement(tag);
  if (text !== undefined) element.textContent = String(text);
  if (className) element.className = className;
  return element;
}

function validProbability(value) { return finite(value) && value >= 0 && value <= 1; }
function screened(candidate) {
  return typeof candidate.off_target_status === "string" &&
    !/not|pending|unknown|unscanned|未|待/i.test(candidate.off_target_status) &&
    finite(field(candidate, "other_gene_exact_compatible", "other_gene_exact_ag"));
}
function statusLabel(value) { return clean(value).replaceAll("_", " "); }
function candidateKey(candidate, index) { return clean(candidate.candidate_id || `candidate_${index + 1}`); }

function validateResult(result) {
  if (!result || typeof result !== "object" || Array.isArray(result)) throw new Error("The JSON root must be a result object.");
  if (!result.target || typeof result.target !== "object") throw new Error("The target object is missing.");
  if (!Array.isArray(result.candidates)) throw new Error("The candidates array is missing.");
  if (result.candidates.length > 100000) throw new Error("More than 100,000 candidates. Export a smaller review file.");
  if (result.candidates.some((item) => !item || typeof item !== "object" || Array.isArray(item))) throw new Error("The candidates array contains a non-object entry.");
  if (result.offtarget_hits !== undefined && !Array.isArray(result.offtarget_hits)) throw new Error("offtarget_hits must be an array.");
  return result;
}

async function loadFile(event) {
  const file = event.target.files[0];
  if (!file) return;
  if (file.size > MAX_FILE_BYTES) { $("file-status").textContent = "File exceeds 50 MB. Export a candidate-summary JSON file."; return; }
  try {
    loadResult(JSON.parse(await file.text()), file.name);
  } catch (error) { $("file-status").textContent = `Could not open file: ${error.message}`; }
}

function loadResult(result, name) {
  const parsed = validateResult(result);
  // Keep the previous working result visible if a later file fails to load.
  state.result = parsed;
  state.selected = parsed.candidates[0] || null;
  $("file-status").textContent = `${name} · ${parsed.candidates.length} candidate designs`;
  $("empty-state").hidden = true;
  $("results").hidden = false;
  $("search").value = "";
  $("sort").value = "original";
  $("analysis-kind").textContent = name === "BIRC5 saved example" ? "Example · tiny reference" : "Saved analysis";
  renderSummary(); renderTable(); renderDetail();
}

function renderSummary() {
  const { target, candidates } = state.result;
  const metadata = state.result.metadata || {};
  const sequence = typeof target.sequence === "string" ? target.sequence : "";
  $("target-label").textContent = [target.gene, target.accession].filter(Boolean).join(" · ") || "Accession not specified";
  $("target-length").textContent = `${sequence ? sequence.length : clean(target.length_nt)} nt · 1-based coordinates`;
  $("candidate-count").textContent = String(candidates.length);
  $("accessible-count").textContent = `${candidates.filter((item) => validProbability(item.p_unpaired_target)).length} / ${candidates.length}`;
  $("screen-count").textContent = `${candidates.filter(screened).length} / ${candidates.length}`;
  const reference = metadata.offtarget || metadata.reference || metadata.transcriptome || {};
  $("reference-label").textContent = clean(field(metadata, "scope_label", "reference_label", "reference_fasta", "reference_path") || field(reference, "fasta_path", "label", "path", "sha256"));
  $("metadata-json").textContent = JSON.stringify(metadata, null, 2);
  $("run-warnings").replaceChildren();
  const warnings = [...(Array.isArray(metadata.warnings) ? metadata.warnings : []), ...(Array.isArray(reference.warnings) ? reference.warnings : [])];
  const access = metadata.accessibility || {};
  if (Array.isArray(access.warnings)) warnings.push(...access.warnings);
  if (typeof access.warning === "string") warnings.push(access.warning);
  if (candidates.some((item) => !screened(item))) warnings.push("Some designs were not screened against the reference. Unscreened designs cannot be treated as off-target-free; ranking is limited to the evaluated set.");
  if (candidates.some((item) => !validProbability(item.p_unpaired_target))) warnings.push("Accessibility is unavailable for some designs. Compare scores within the same analysis stage.");
  if (!warnings.length) warnings.push("Reference identity, sequence digests, folding parameters, and scoring scope are recorded below.");
  warnings.forEach((warning) => $("run-warnings").append(node("li", clean(warning))));
}

function renderTable() {
  const query = $("search").value.trim().toUpperCase();
  state.filtered = state.result.candidates.filter((item) => [item.candidate_id, item.target_rna, item.start_1based, item.cut_after_1based].some((value) => clean(value).toUpperCase().includes(query)));
  const mode = $("sort").value;
  if (mode !== "original") {
    const key = { score: "ranking_score", accessibility: "p_unpaired_target", coordinate: "cut_after_1based" }[mode];
    state.filtered.sort((a, b) => {
      if (!finite(a[key]) && !finite(b[key])) return 0;
      if (!finite(a[key])) return 1;
      if (!finite(b[key])) return -1;
      return mode === "coordinate" ? a[key] - b[key] : b[key] - a[key];
    });
  }
  $("candidate-body").replaceChildren();
  $("table-status").textContent = `${state.filtered.length} designs match · showing ${Math.min(state.filtered.length, MAX_TABLE_ROWS)}. Select a row to inspect the design.`;
  $("download-filtered").disabled = state.filtered.length === 0;
  state.filtered.slice(0, MAX_TABLE_ROWS).forEach((candidate, index) => {
    const row = node("tr");
    row.tabIndex = 0;
    if (candidate === state.selected) row.className = "selected";
    const cut = finite(candidate.cut_after_1based) ? `${candidate.cut_after_1based}|${candidate.cut_after_1based + 1}` : "Not provided";
    const values = [candidateKey(candidate, index), cut, clean(candidate.target_rna), `${clean(candidate.left_length)} + ${clean(candidate.right_length)}`, probability(candidate.p_unpaired_target), number(candidate.ranking_score), integer(field(candidate,"other_gene_exact_compatible","other_gene_exact_ag"))];
    values.forEach((value, index) => row.append(node("td", value, index === 2 ? "mono" : "")));
    const status = node("td");
    const evidence = screened(candidate) ? (validProbability(candidate.p_unpaired_target) ? "Screened + folding" : "Screened") : (validProbability(candidate.p_unpaired_target) ? "Folding only" : "Sequence only");
    const badge = node("span", evidence, screened(candidate) ? "pill" : "pill warning");
    badge.setAttribute("title", statusLabel(candidate.ranking_stage || candidate.off_target_status));
    status.append(badge);
    row.append(status);
    const select = () => { state.selected = candidate; renderTable(); renderDetail(); };
    row.addEventListener("click", select);
    row.addEventListener("keydown", (event) => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); select(); } });
    $("candidate-body").append(row);
  });
  renderMap();
}

function svgElement(tag, attributes, text) {
  const item = document.createElementNS("http://www.w3.org/2000/svg", tag);
  Object.entries(attributes || {}).forEach(([key, value]) => item.setAttribute(key, String(value)));
  if (text !== undefined) item.textContent = String(text);
  return item;
}

function renderMap() {
  const sequence = state.result.target.sequence;
  const sequenceLength = typeof sequence === "string" ? sequence.length : Number(state.result.target.length_nt);
  const length = finite(sequenceLength) && sequenceLength > 0 ? sequenceLength : Math.max(1, ...state.result.candidates.map((candidate) => finite(candidate.end_1based) ? candidate.end_1based : 0));
  // Match the actual panel width so SVG labels remain readable on mobile.
  const width = Math.max(300, Math.round($("position-map").clientWidth || 1100));
  const height = 235, left = 43, right = 17, top = 20, bottom = 54;
  const svg = svgElement("svg", { viewBox: `0 0 ${width} ${height}`, role: "img", "aria-label": "Cut-site coordinates and recognition-window unpaired probability; grey indicates unavailable values" });
  const x = (position) => left + (width - left - right) * (position - 1) / Math.max(1, length - 1);
  const y = (p) => top + (height - top - bottom) * (1 - p);
  for (const p of [0, .25, .5, .75, 1]) {
    svg.append(svgElement("line", { x1: left, y1: y(p), x2: width - right, y2: y(p), stroke: "#e0e8eb", "stroke-width": 1 }));
    svg.append(svgElement("text", { x: left - 10, y: y(p) + 4, "text-anchor": "end", fill: "#657980", "font-size": 11 }, `${100 * p}%`));
  }
  [1, Math.round(length / 4), Math.round(length / 2), Math.round(3 * length / 4), length].filter((position, index, array) => array.indexOf(position) === index).forEach((position) => {
    svg.append(svgElement("text", { x: x(position), y: height - 30, "text-anchor": "middle", fill: "#657980", "font-size": 11 }, position));
  });
  svg.append(svgElement("text", { x: width / 2, y: height - 8, "text-anchor": "middle", fill: "#657980", "font-size": 11 }, "Target mRNA position (1-based)"));
  // Sampling is only for a very large visual map, not for scoring/exporting.
  const stride = Math.max(1, Math.ceil(state.filtered.length / 10000));
  state.filtered.filter((candidate, index) => candidate === state.selected || index % stride === 0).forEach((candidate) => {
    if (!finite(candidate.cut_after_1based) || candidate.cut_after_1based < 1 || candidate.cut_after_1based > length) return;
    const selected = candidate === state.selected;
    const valid = validProbability(candidate.p_unpaired_target);
    const point = svgElement("circle", { cx: x(candidate.cut_after_1based), cy: valid ? y(candidate.p_unpaired_target) : height - bottom + 11, r: selected ? 6 : 3, fill: selected ? "#c75a34" : valid ? "#168879" : "#a9b6bb", opacity: selected ? 1 : .55, tabindex: 0, role: "button", "aria-label": `${clean(candidate.candidate_id)}, accessibility ${probability(candidate.p_unpaired_target)}` });
    point.append(svgElement("title", {}, `${clean(candidate.candidate_id)} · ${candidate.cut_after_1based}|${candidate.cut_after_1based + 1} · ${probability(candidate.p_unpaired_target)}`));
    const select = () => { state.selected = candidate; renderTable(); renderDetail(); };
    point.addEventListener("click", select);
    point.addEventListener("keydown", (event) => { if (event.key === "Enter" || event.key === " ") { event.preventDefault(); select(); } });
    svg.append(point);
  });
  $("position-map").replaceChildren(svg);
  $("map-legend").textContent = stride > 1 ? `Map sampled at 1 in ${stride} designs; the table and export retain all records. Grey: unavailable.` : "Y axis: joint unpaired probability. Grey: unavailable; orange: selected design.";
}

function addValues(id, rows) {
  const list = $(id); list.replaceChildren();
  rows.forEach(([label, value]) => { list.append(node("dt", label), node("dd", value)); });
}

function renderContext(candidate) {
  const container = $("sequence-context"); container.replaceChildren();
  const sequence = state.result.target.sequence;
  const { start_1based: start, end_1based: end, cut_after_1based: cut } = candidate;
  if (typeof sequence !== "string" || !Number.isInteger(start) || !Number.isInteger(end) || !Number.isInteger(cut) || start < 1 || end > sequence.length || end < start || cut < start || cut >= end) {
    container.append(node("span", clean(candidate.target_rna)));
    return;
  }
  const before = Math.max(0, start - 1 - 12), after = Math.min(sequence.length, end + 12);
  container.append(node("span", `${before + 1} … `, "flank"), node("span", sequence.slice(before, start - 1), "flank"));
  container.append(node("span", sequence.slice(start - 1, cut - 1), "arm"));
  container.append(node("span", sequence[cut - 1], "motif"), node("span", "|", "cut"), node("span", sequence[cut], "motif"));
  container.append(node("span", sequence.slice(cut + 1, end), "arm"));
  container.append(node("span", sequence.slice(end, after), "flank"), node("span", ` … ${after}`, "flank"));
}

function renderDetail() {
  const candidate = state.selected;
  $("download-selected").disabled = !candidate;
  $("download-fasta").disabled = !candidate;
  if (!candidate) {
    $("detail-title").textContent = "No candidate designs"; $("hit-status").textContent = "This result contains no candidate designs.";
    ["detail-summary", "sequence-context", "sequence-values", "feature-values", "advanced-feature-values", "candidate-warnings", "candidate-json", "hit-body"].forEach((id) => $(id).replaceChildren());
    ["selected-name", "selected-site", "selected-accessibility", "selected-hits", "preview-left-arm", "preview-right-arm", "hits-mm0", "hits-mm1", "hits-mm2"].forEach((id) => { $(id).textContent = "—"; });
    return;
  }
  $("selected-name").textContent = clean(candidate.candidate_id);
  $("selected-site").textContent = `Window ${clean(candidate.start_1based)}–${clean(candidate.end_1based)} nt · arms ${clean(candidate.left_length)} + ${clean(candidate.right_length)} nt`;
  $("selected-accessibility").textContent = probability(candidate.p_unpaired_target);
  $("selected-hits").textContent = screened(candidate) ? integer(field(candidate,"other_gene_exact_compatible","other_gene_exact_ag")) : "Not evaluated";
  $("preview-left-arm").textContent = clean(candidate.left_arm_dna);
  $("preview-right-arm").textContent = clean(candidate.right_arm_dna);
  [["hits-mm0","other_gene_exact_compatible","other_gene_exact_ag"],["hits-mm1","other_gene_1mm_compatible","other_gene_1mm_ag"],["hits-mm2","other_gene_2mm_compatible","other_gene_2mm_ag"]].forEach(([id,modern,legacy]) => {
    $(id).textContent = screened(candidate) ? integer(field(candidate, modern, legacy)) : "Not evaluated";
  });
  $("detail-title").textContent = clean(candidate.candidate_id);
  $("detail-summary").textContent = `Window ${clean(candidate.start_1based)}–${clean(candidate.end_1based)} nt · ${clean(candidate.left_length)} + ${clean(candidate.right_length)} nt arms`;
  renderContext(candidate);
  addValues("sequence-values", [["Target RNA, 5′→3′", clean(candidate.target_rna)], ["Left DNA arm, 5′→3′", clean(candidate.left_arm_dna)], ["Right DNA arm, 5′→3′", clean(candidate.right_arm_dna)], ["Left partzyme", candidate.left_partzyme ? clean(candidate.left_partzyme) : "Arms only; no scaffold specified"], ["Right partzyme", candidate.right_partzyme ? clean(candidate.right_partzyme) : "Arms only; no scaffold specified"]]);
  const diagnostic = (state.result.structure_diagnostics || {})[candidate.candidate_id] || null;
  const featureRows = [["Cut-motif rule", clean(candidate.motif_pattern || "AG")], ["Window accessibility", probability(candidate.p_unpaired_target)], ["Left-arm accessibility", probability(candidate.p_unpaired_left_arm)], ["Right-arm accessibility", probability(candidate.p_unpaired_right_arm)], ["Opening-energy proxy, kcal/mol", number(candidate.opening_dg_kcal_mol)], ["Hybridization model", clean(candidate.binding_model)], ["Left / right ΔG, kcal/mol", `${number(candidate.left_binding_dg_kcal_mol)} / ${number(candidate.right_binding_dg_kcal_mol)}`], ["Arm reference Tm, °C", `${number(candidate.left_tm_c, 1)} / ${number(candidate.right_tm_c, 1)}`], ["Arm-length trade-off score", number(candidate.length_tradeoff_score)], ["Trade-off factors", clean(candidate.length_tradeoff_reasons)], ["Self / cross-complementary run", `${integer(candidate.max_self_complementary_run)} / ${integer(candidate.max_cross_complementary_run)}`], ["Other-gene compatible hits: 0 / 1 / 2 mismatches", `${integer(field(candidate,"other_gene_exact_compatible","other_gene_exact_ag"))} / ${integer(field(candidate,"other_gene_1mm_compatible","other_gene_1mm_ag"))} / ${integer(field(candidate,"other_gene_2mm_compatible","other_gene_2mm_ag"))}`], ["Unannotated hits", integer(candidate.unknown_gene_hits)], ["3D diagnostic / scoring", diagnostic ? `${clean(diagnostic.mapping_status)} / Not used in ranking` : clean(candidate.structure_status || candidate.three_d_status)], ["3D target coverage", diagnostic ? probability(diagnostic.target_coordinate_coverage_fraction) : "Not evaluated"], ["Target pLDDT (declared source)", diagnostic ? number(diagnostic.target_mean_plddt, 1) : "Not evaluated"]];
  const primaryFeatures = new Set(["Cut-motif rule", "Window accessibility", "Left / right ΔG, kcal/mol",
    "Arm reference Tm, °C", "Arm-length trade-off score", "Self / cross-complementary run"]);
  addValues("feature-values", featureRows.filter(([label]) => primaryFeatures.has(label)));
  addValues("advanced-feature-values", featureRows.filter(([label]) => !primaryFeatures.has(label)));
  const warnings = [];
  if (Array.isArray(candidate.warnings)) warnings.push(...candidate.warnings);
  if (Array.isArray(candidate.binding_warnings)) warnings.push(...candidate.binding_warnings);
  if (diagnostic && Array.isArray(diagnostic.warnings)) warnings.push(...diagnostic.warnings);
  if (!screened(candidate)) warnings.push("Reference screening is incomplete or hit counts are unknown for this design.");
  if (!validProbability(candidate.p_unpaired_target)) warnings.push("Full-window accessibility is unavailable. Check the target sequence and the associated lunp file.");
  if (!candidate.left_partzyme || !candidate.right_partzyme) warnings.push("This design contains substrate-binding arms only. Full partzymes require a specified catalytic core, split position, and trigger scaffold.");
  warnings.push("The dinucleotide rule is a design constraint. Cleavage specificity and rates require independent validation.");
  warnings.push("Tm describes each isolated DNA–RNA arm under model reference conditions, not the assembled MNAzyme.");
  $("candidate-warnings").replaceChildren(...warnings.map((warning) => node("li", clean(warning))));
  $("candidate-json").textContent = JSON.stringify(diagnostic ? { candidate, structure_diagnostic: diagnostic } : candidate, null, 2);
  renderHits(candidate);
}

function renderHits(candidate) {
  const hits = (state.result.offtarget_hits || []).filter((hit) => hit && hit.candidate_id === candidate.candidate_id);
  $("hit-body").replaceChildren();
  $("hit-status").textContent = !screened(candidate) ? "This design has not been reference-screened. No displayed records does not mean zero off-target matches." : `${hits.length} records saved · showing ${Math.min(hits.length, MAX_TABLE_ROWS)}. Total hit counts can exceed the saved detail limit.`;
  hits.slice(0, MAX_TABLE_ROWS).forEach((hit) => {
    const row = node("tr");
    [field(hit, "transcript_accession", "accession", "transcript_id", "transcript_header"), [field(hit, "gene", "gene_id"), field(hit, "gene_classification", "category")].filter((value) => value !== null).join(" / ") || null, field(hit, "start_1based", "position_1based"), field(hit, "arm_mismatches", "mismatches", "mismatch_count"), field(hit, "cut_motif", "ag_intact", "motif"), field(hit, "window_rna", "matched_rna", "rna", "sequence", "target_rna")].forEach((value, index) => row.append(node("td", clean(value), index === 5 ? "mono" : "")));
    $("hit-body").append(row);
  });
}

function csvValue(value) {
  let text = value === null || value === undefined ? "" : typeof value === "object" ? JSON.stringify(value) : String(value);
  // Prevent spreadsheet formula execution on data downloaded from an untrusted JSON.
  if (/^(?:\s*[=+\-@]|[\t\r\n])/.test(text) && !finite(value)) text = `'${text}`;
  return `"${text.replaceAll('"', '""')}"`;
}
function downloadSelected() {
  if (!state.selected) return;
  const fields = Object.keys(state.selected);
  const content = "\ufeff" + fields.map(csvValue).join(",") + "\r\n" + fields.map((key) => csvValue(state.selected[key])).join(",") + "\r\n";
  const blob = new Blob([content], { type: "text/csv;charset=utf-8" });
  const url = URL.createObjectURL(blob);
  const link = node("a"); link.href = url; link.download = "selected_candidate.csv";
  document.body.append(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

function saveText(content, mime, filename) {
  const url = URL.createObjectURL(new Blob([content], { type: mime }));
  const link = node("a"); link.href = url; link.download = filename;
  document.body.append(link); link.click(); link.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}
function filteredCsv() {
  const fields = [...new Set(state.filtered.flatMap((candidate) => Object.keys(candidate)))];
  return "\ufeff" + fields.map(csvValue).join(",") + "\r\n" +
    state.filtered.map((candidate) => fields.map((key) => csvValue(candidate[key])).join(",")).join("\r\n") + "\r\n";
}
function selectedFasta() {
  if (!state.selected) return "";
  const candidate = state.selected;
  const id = clean(candidate.candidate_id).replace(/[^A-Za-z0-9_.-]/g, "_");
  return [["target_RNA",candidate.target_rna],["left_DNA_arm",candidate.left_arm_dna],["right_DNA_arm",candidate.right_arm_dna],
    ["left_partzyme",candidate.left_partzyme],["right_partzyme",candidate.right_partzyme]]
    .filter(([, sequence]) => typeof sequence === "string" && /^[ACGTUacgtu]+$/.test(sequence))
    .map(([label, sequence]) => ">" + id + "_" + label + "\n" + sequence + "\n").join("");
}
function loadDemo() {
  if (window.ZYME_DEMO) { loadResult(window.ZYME_DEMO, "BIRC5 saved example"); return; }
  const button = $("load-demo");
  button.disabled = true;
  $("file-status").textContent = "Opening the saved BIRC5 example…";
  const script = document.createElement("script");
  script.src = "demo-data.js";
  script.onload = () => {
    button.disabled = false;
    try {
      if (!window.ZYME_DEMO) throw new Error("The example data did not load.");
      loadResult(window.ZYME_DEMO, "BIRC5 saved example");
    } catch (error) { $("file-status").textContent = error.message; }
  };
  script.onerror = () => {
    button.disabled = false;
    $("file-status").textContent = "No bundled example is available here. Open your result.json instead.";
  };
  document.body.append(script);
}

$("result-file").addEventListener("change", loadFile);
$("search").addEventListener("input", () => { if (state.result) renderTable(); });
$("sort").addEventListener("change", () => { if (state.result) renderTable(); });
$("download-selected").addEventListener("click", downloadSelected);
$("download-filtered").addEventListener("click", () => { if (state.filtered.length) saveText(filteredCsv(), "text/csv;charset=utf-8", "filtered_candidates.csv"); });
$("download-fasta").addEventListener("click", () => { if (state.selected) saveText(selectedFasta(), "text/plain;charset=utf-8", "selected_design.fasta"); });
$("load-demo").addEventListener("click", loadDemo);
if (typeof document.querySelectorAll === "function") {
  document.querySelectorAll(".side-nav a[href^='#']").forEach((link) => {
    link.addEventListener("click", () => {
      document.querySelectorAll(".side-nav .active").forEach((item) => item.classList.remove("active"));
      link.classList.add("active");
    });
  });
}
if (window.ZYME_RESULT) {
  try { loadResult(window.ZYME_RESULT, window.ZYME_IS_EXAMPLE ? "BIRC5 saved example" : "result.json"); }
  catch (error) { $("file-status").textContent = `Bundled results could not be loaded: ${error.message}. Choose result.json manually.`; }
}

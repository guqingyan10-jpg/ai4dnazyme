"use strict";

// Reproducible, dependency-free UI logic checks. This minimal DOM does not
// replace the separate headless-browser visual QA of layout and real events.
// Run: node tests/test_ui.cjs  (from any working directory)
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const assert = require("node:assert/strict");

class Element {
  constructor(tag = "div") {
    this.tagName = tag; this.children = []; this.listeners = {};
    this.value = ""; this.files = []; this.hidden = false; this._text = "";
  }
  set textContent(value) { this._text = String(value); this.children = []; }
  get textContent() { return this._text + this.children.map((child) => child.textContent || "").join(""); }
  append(...items) { this.children.push(...items); }
  replaceChildren(...items) { this._text = ""; this.children = items; }
  setAttribute(name, value) { this[name] = value; }
  addEventListener(name, callback) { this.listeners[name] = callback; }
  remove() {}
  click() { if (this.listeners.click) this.listeners.click(); }
}

const elements = {};
const doc = {
  body: new Element("body"),
  getElementById: (id) => elements[id] ??= new Element(),
  createElement: (tag) => new Element(tag),
  createElementNS: (_namespace, tag) => new Element(tag),
};
const candidate = {
  candidate_id: "AG8_L7_R7", start_1based: 1, end_1based: 16, cut_after_1based: 8,
  target_rna: "AUUGCAAAGGAAACCA", left_arm_dna: "TTGCAAT", right_arm_dna: "TGGTTTC",
  left_length: 7, right_length: 7, p_unpaired_target: 0,
  p_unpaired_left_arm: .29865, p_unpaired_right_arm: .34357, ranking_score: 45,
  ranking_stage: "screened_with_accessibility", off_target_status: "screened",
  other_gene_exact_ag: 0, other_gene_1mm_ag: 1, other_gene_2mm_ag: 2,
  left_tm_c: 4.7, right_tm_c: 5.1,
};
const missing = {
  ...candidate, candidate_id: "missing", p_unpaired_target: null, ranking_score: 20,
  off_target_status: "not_screened", ranking_stage: "sequence_only_not_screened",
  other_gene_exact_ag: null, other_gene_1mm_ag: null, other_gene_2mm_ag: null,
};
const result = {
  metadata: { warnings: [], scope_label: "small fixture, not human", offtarget: { warnings: ["fixture warning"] } },
  target: { accession: "NM_FIXTURE.1<img src=x onerror=alert(1)>", gene: "TEST", sequence: candidate.target_rna },
  candidates: [candidate, missing],
  offtarget_hits: [{
    candidate_id: candidate.candidate_id, transcript_accession: "NM_OTHER.1",
    gene: "OTHER", gene_classification: "other_gene", start_1based: 3,
    arm_mismatches: 1, cut_motif: "AG", window_rna: "AUUGCUAAGGAAACCA",
    target_rna: candidate.target_rna,
  }],
  structure_diagnostics: { [candidate.candidate_id]: {
    mapping_status: "explicit_residue_mapping", target_coordinate_coverage_fraction: .5,
    target_mean_plddt: null, warnings: ["PDB is not a cleavage predictor"],
  } },
};
const context = vm.createContext({ document: doc, window: { ZYME_RESULT: result }, console, Blob, URL, setTimeout });
const source = fs.readFileSync(path.resolve(__dirname, "../ui/app.js"), "utf8");
const html = fs.readFileSync(path.resolve(__dirname, "../ui/index.html"), "utf8");
const structureHtml = fs.readFileSync(path.resolve(__dirname, "../ui/structure.html"), "utf8");
const stylesheet = fs.readFileSync(path.resolve(__dirname, "../ui/styles.css"), "utf8");
vm.runInContext(source, context);
let passed = 0;
const check = (callback) => { callback(); passed += 1; };

check(() => assert.ok(html.includes('<html lang="en">')));
check(() => assert.ok(structureHtml.includes('<html lang="en">')));
check(() => assert.ok(!/[\u3400-\u9fff]/u.test(html + structureHtml)));
check(() => assert.ok(!/https?:\/\//i.test(html + structureHtml + stylesheet)));
check(() => assert.ok(html.includes('<details class="notes"><summary>Design notes and interpretation')));
check(() => assert.ok(structureHtml.includes('<details><summary>Highlight a region')));
check(() => assert.ok(stylesheet.includes('"Segoe UI",system-ui')));
check(() => assert.equal(elements["candidate-count"].textContent, "2"));
check(() => assert.equal(elements["screen-count"].textContent, "1 / 2"));
check(() => assert.equal(elements["accessible-count"].textContent, "1 / 2"));
check(() => assert.equal(elements["candidate-body"].children[0].children[4].textContent, "0.00%"));
check(() => assert.equal(elements["candidate-body"].children[1].children[4].textContent, "Not available"));
check(() => assert.equal(elements["candidate-body"].children[1].children[6].textContent, "Not evaluated"));
check(() => assert.equal(elements["hit-body"].children[0].children[0].textContent, "NM_OTHER.1"));
check(() => assert.equal(elements["hit-body"].children[0].children[5].textContent, "AUUGCUAAGGAAACCA"));
check(() => assert.ok(elements["target-label"].textContent.includes("<img")));
check(() => assert.ok(elements["sequence-context"].textContent.includes("A|G")));
check(() => assert.equal(elements["reference-label"].textContent, "small fixture, not human"));
check(() => assert.ok(elements["run-warnings"].textContent.includes("fixture warning")));
check(() => assert.ok(elements["advanced-feature-values"].textContent.includes("50.00%")));
check(() => assert.ok(elements["candidate-warnings"].textContent.includes("PDB is not a cleavage predictor")));
check(() => assert.ok(!source.includes(".innerHTML")));
check(() => assert.equal(vm.runInContext("csvValue('=2+2')", context), '"\'=2+2"'));
check(() => assert.equal(vm.runInContext("csvValue('-SUM(1,2)')", context), '"\'-SUM(1,2)"'));
check(() => assert.equal(vm.runInContext("csvValue(' =2+2')", context), '"\' =2+2"'));
check(() => assert.equal(vm.runInContext("csvValue('\\n=2+2')", context), '"\'\n=2+2"'));
check(() => assert.equal(vm.runInContext("csvValue(-5)", context), '"-5"'));
check(() => assert.equal(vm.runInContext("csvValue(null)", context), '""'));
check(() => assert.equal(vm.runInContext("csvValue('a\"b')", context), '"a""b"'));

elements["candidate-body"].children[1].listeners.click();
check(() => assert.equal(elements["detail-title"].textContent, "missing"));
check(() => assert.ok(elements["hit-status"].textContent.includes("does not mean zero off-target matches")));
elements.search.value = "AG8";
elements.search.listeners.input();
check(() => assert.equal(elements["candidate-body"].children.length, 1));
check(() => assert.throws(() => vm.runInContext("validateResult({target:{},candidates:[null]})", context)));

vm.runInContext('loadResult({target:{sequence:"ACGU"},candidates:[]},"empty")', context);
check(() => assert.equal(elements["detail-title"].textContent, "No candidate designs"));
check(() => assert.equal(elements["sequence-context"].textContent, ""));
check(() => assert.equal(elements["feature-values"].textContent, ""));
check(() => assert.equal(elements["selected-name"].textContent, "—"));
check(() => assert.equal(elements["hits-mm0"].textContent, "—"));
check(() => assert.equal(elements["download-fasta"].disabled, true));
vm.runInContext("loadResult(" + JSON.stringify(result) + ", 'fixture')", context);
check(() => assert.equal(elements["selected-hits"].textContent, "0"));
check(() => assert.equal(elements["hits-mm1"].textContent, "1"));
check(() => assert.ok(vm.runInContext("selectedFasta()", context).includes(">AG8_L7_R7_left_DNA_arm\nTTGCAAT\n")));
check(() => assert.ok(vm.runInContext("filteredCsv()", context).includes("AG8_L7_R7")));
elements["candidate-body"].children[1].listeners.click();
check(() => assert.equal(elements["selected-hits"].textContent, "Not evaluated"));
check(() => assert.equal(elements["hits-mm2"].textContent, "Not evaluated"));
check(() => assert.ok(html.includes('id="load-demo"')));
check(() => assert.ok(html.includes('href="guide.html"')));
check(() => assert.ok(html.includes("connect-src 'none'")));
console.log(`UI regression: ${passed} assertions passed (English/offline layout, null vs zero, data scope, actual matched sequence, 3D mapping, safe text/CSV, selection/filter and empty reset).`);

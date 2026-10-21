"use strict";
// Offline independent parser / mapping QA. No browser/network/library install.
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const vm = require("node:vm");
const io = require("../ui/structure.js");

function pdbAtom({serial=1, atom="P", alt="", base="A", chain="C", resid=1,
                  insertion="", x="1.000", y="2.000", z="3.000", b="85.00"}={}) {
  const chars = Array(80).fill(" ");
  const put = (start,width,value,right=false) => {
    const text = (right ? String(value).padStart(width) : String(value).padEnd(width)).slice(0,width);
    for (let i=0;i<width;i++) chars[start+i]=text[i];
  };
  put(0,6,"ATOM");put(6,5,serial,true);put(12,4,atom);put(16,1,alt);
  put(17,3,base,true);put(21,1,chain);put(22,4,resid,true);put(26,1,insertion);
  put(30,8,x,true);put(38,8,y,true);put(46,8,z,true);put(54,6,"1.00",true);
  put(60,6,b,true);put(76,2,"P",true);
  return chars.join("");
}
const cols = ["group_PDB","id","type_symbol","label_atom_id","label_alt_id","label_comp_id",
              "label_asym_id","label_seq_id","auth_asym_id","auth_seq_id","pdbx_PDB_ins_code",
              "Cartn_x","Cartn_y","Cartn_z","B_iso_or_equiv","pdbx_PDB_model_num"];
const cif = rows => "data_fixture\nloop_\n" + cols.map(c=>"_atom_site."+c).join("\n") + "\n" + rows.join("\n") + "\n#\n";
const row = (updates={}) => {
  const values={group_PDB:"ATOM",id:"1",type_symbol:"P",label_atom_id:"P",label_alt_id:".",
    label_comp_id:"A",label_asym_id:"C",label_seq_id:"1",auth_asym_id:"Q",auth_seq_id:"101",
    pdbx_PDB_ins_code:"?",Cartn_x:"0",Cartn_y:"2",Cartn_z:"3",B_iso_or_equiv:"85",pdbx_PDB_model_num:"1",...updates};
  return cols.map(c=>values[c]).join(" ");
};
const failures=[];let passed=0;
async function check(name,callback) {
  try {await callback();passed++;console.log("PASS: "+name);}
  catch(error){failures.push({name,error});console.error("FAIL: "+name+" — "+error.message);}
}
function rejectOrDrop(callback) {
  let atoms;
  try {atoms=callback();} catch(error){return;}
  assert.equal(atoms.length,0,"missing/invalid coordinates must not be synthesized as zero");
}

class Element {
  constructor() {this.listeners={};this.children=[];this.value="";this.checked=false;this.files=[];this._text="";this.style={};}
  set textContent(value){this._text=String(value);}
  get textContent(){return this._text;}
  addEventListener(event,callback){this.listeners[event]=callback;}
  append(...items){this.children.push(...items);}
  replaceChildren(...items){this.children=items;}
  click(){}
}
function mockViewer() {
  const elements={};const painted=[];
  const ctx={fillStyle:"",setTransform(){},fillRect(){},fillText(){},beginPath(){},moveTo(){},lineTo(){},stroke(){},arc(){},fill(){painted.push(this.fillStyle);}};
  const document={getElementById(id){return elements[id]??=new Element();},createElement(){return new Element();},createTextNode(text){return {textContent:text};}};
  const canvas=document.getElementById("structure-canvas");canvas.getContext=()=>ctx;canvas.clientWidth=900;canvas.clientHeight=570;
  document.getElementById("representation").value="atoms";
  document.getElementById("colour-mode").value="chain";
  const context=vm.createContext({document,console,window:{devicePixelRatio:1,addEventListener(){}},module:{exports:{}}});
  vm.runInContext(fs.readFileSync(path.resolve(__dirname,"../ui/structure.js"),"utf8"),context);
  return {elements,painted,async load(id,text,name){const element=document.getElementById(id);const file={name,size:text.length,text:async()=>text};element.files=[file];await element.listeners.change({target:element});}};
}

(async()=>{
  await check("PDB first model only and canonical altloc",()=>{
    const text="MODEL        7\n"+pdbAtom({x:"7.000"})+"\n"+pdbAtom({serial:2,alt:"B",x:"8.000"})+
      "\n"+pdbAtom({serial:3,alt:"A",atom:"C1'",x:"9.000"})+"\nENDMDL\nMODEL        8\n"+pdbAtom({x:"80.000"})+"\nENDMDL\n";
    const atoms=io.parsePDB(text);assert.equal(atoms.length,2);assert.deepEqual(atoms.map(a=>a.x),[7,9]);
  });
  await check("PDB blank and truncated coordinates never become fabricated origins",()=>{
    rejectOrDrop(()=>io.parsePDB(pdbAtom({x:""})));
    rejectOrDrop(()=>io.parsePDB(pdbAtom().slice(0,30)));
  });
  await check("PDB missing confidence is unknown, not zero pLDDT",()=>{
    assert.ok(!Number.isFinite(io.parsePDB(pdbAtom({b:""}))[0].b));
  });
  await check("mmCIF first encountered model, label identifiers and quoted apostrophe",()=>{
    const atoms=io.parseCIF(cif([row({pdbx_PDB_model_num:"7",label_atom_id:'"C1\'"'}),
      row({id:"2",pdbx_PDB_model_num:"8",Cartn_x:"80"}),
      row({id:"3",pdbx_PDB_model_num:"7",label_alt_id:"B"}),
      row({id:"4",pdbx_PDB_model_num:"7",label_alt_id:"A",label_seq_id:"2",label_comp_id:"G"})]));
    assert.equal(atoms.length,2);assert.equal(atoms[0].name,"C1'");assert.equal(atoms[0].x,0);
    assert.equal(atoms[0].chain,"C");assert.equal(atoms[0].resid,"1");assert.equal(atoms[0].model,"7");
  });
  await check("mmCIF auth fallback and missing B-factor",()=>{
    const atom=io.parseCIF(cif([row({label_asym_id:".",label_seq_id:"?",B_iso_or_equiv:"?"})]))[0];
    assert.equal(atom.chain,"Q");assert.equal(atom.resid,"101");assert.ok(Number.isNaN(atom.b));
  });
  await check("mmCIF missing/nonfinite coordinates rejected or dropped",()=>{
    for(const value of [".","?",'""',"NaN","Infinity","nonnumeric"])
      rejectOrDrop(()=>io.parseCIF(cif([row({Cartn_y:value})])));
  });
  await check("CIF comments, multiline quoted metadata, quoted keywords do not corrupt loops",()=>{
    const text="data_fixture\n_struct.title\n;quoted title\nloop_ _not_a_tag\n;\n# comment\n"+
      cif([row({label_atom_id:'"_quoted_atom"'})]).replace("data_fixture\n","");
    const atoms=io.parseCIF(text);assert.equal(atoms.length,1);assert.equal(atoms[0].name,"_quoted_atom");
    assert.throws(()=>io.cifTokens("data_bad\n'title not closed"));
    assert.throws(()=>io.parseCIF("data_bad\nloop_\n_atom_site.id\n_atom_site.Cartn_x\n1\n#\n"));
  });
  await check("automatic mapping never claims success while using wrong residue numbers",async()=>{
    const viewer=mockViewer();
    await viewer.load("coordinate-file",pdbAtom({resid:101})+"\n"+pdbAtom({serial:2,resid:102,base:"G"}),"offset.pdb");
    viewer.painted.length=0;
    const mapping={RNA_chain:"C",RNA_chain_site_start_1based:1,RNA_chain_site_end_1based:2,
      RNA_chain_transcript_start_1based:401,RNA_chain_transcript_end_1based:402,
      chains:[{id:"C",type:"rna",role:"RNA_substrate_context",sequence:"AG"}]};
    await viewer.load("mapping-file",JSON.stringify(mapping),"mapping.json");
    const status=viewer.elements["structure-status"].textContent;
    assert.ok(!status.includes("sequence verified")||viewer.painted.filter(c=>c==="#ff6286").length===2,
              "correct sequence alone cannot prove residue-number mapping; either explicitly map or reject");
  });
  await check("sequence mismatch and insertion codes never pass automatic mapping",async()=>{
    for(const [text,expected] of [
      [pdbAtom({resid:1})+"\n"+pdbAtom({serial:2,resid:2,base:"G"}),"AC"],
      [pdbAtom({resid:1})+"\n"+pdbAtom({serial:2,resid:1,insertion:"A",base:"G"}),"AG"]]) {
      const viewer=mockViewer();await viewer.load("coordinate-file",text,"ambiguous_mapping.pdb");
      const mapping={RNA_chain:"C",RNA_chain_site_start_1based:1,RNA_chain_site_end_1based:2,
        RNA_chain_transcript_start_1based:401,RNA_chain_transcript_end_1based:402,
        chains:[{id:"C",type:"rna",sequence:expected}]};
      await viewer.load("mapping-file",JSON.stringify(mapping),"mapping.json");
      assert.ok(!viewer.elements["structure-status"].textContent.includes("sequence verified"));
    }
  });
  await check("confidence colours require affirmative pLDDT-source confirmation",async()=>{
    const viewer=mockViewer();await viewer.load("coordinate-file",pdbAtom({chain:"Q",b:"85.00"}),"experimental.pdb");
    viewer.elements["colour-mode"].value="confidence";viewer.painted.length=0;
    viewer.elements["colour-mode"].listeners.change();
    assert.equal(viewer.painted.at(-1),"#57b8ff","experimental B-factor not treated as pLDDT automatically");
    viewer.elements["plddt-confirm"].checked=true;viewer.painted.length=0;
    viewer.elements["plddt-confirm"].listeners.change();
    assert.equal(viewer.painted.at(-1),"#72c7f5","explicitly declared 85 pLDDT may use confidence palette");
  });
  console.log(`Independent structure QA: ${passed} passed, ${failures.length} failed.`);
  if(failures.length)process.exitCode=1;
})();

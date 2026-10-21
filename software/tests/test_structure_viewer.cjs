"use strict";
const assert=require("node:assert/strict");
const io=require("../ui/structure.js");
const pdb="ATOM      1  P     A A   1      10.000  20.000  30.000  1.00 85.00           P  \n"+
          "ATOM      2  P     G A   2      11.000  21.000  31.000  1.00 75.00           P  \nEND\n";
let passed=0;const check=f=>{f();passed++;};
check(()=>assert.equal(io.parsePDB(pdb).length,2));
check(()=>assert.equal(io.summarize(io.parsePDB(pdb))[0].sequence,"AG"));
check(()=>assert.equal(io.parsePDB(pdb)[0].b,85));
const cif=`data_demo
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 1 P P A C 1 1.0 2.0 3.0 90 1
ATOM 2 C "C1'" A C 1 1.1 2.1 3.1 88 1
ATOM 3 P P G C 2 2.0 3.0 4.0 60 1
ATOM 4 P P G C 2 4.0 5.0 6.0 50 2
#
`;
check(()=>assert.equal(io.parseCIF(cif).length,3));
check(()=>assert.equal(io.parseCIF(cif)[1].name,"C1'"));
check(()=>assert.equal(io.summarize(io.parseCIF(cif))[0].sequence,"AG"));
check(()=>assert.equal(io.parseCIF(cif)[0].chain,"C"));
check(()=>assert.equal(io.parseCIF(cif)[0].resid,"1"));
check(()=>assert.throws(()=>io.parse("no coordinates","bad.pdb")));
check(()=>assert.throws(()=>io.parseCIF(cif.replace("ATOM 1 P P A C 1 1.0 2.0 3.0 90 1","ATOM 1 P P A C 1 1.0 2.0"))));
check(()=>assert.throws(()=>io.cifTokens("data_x\n;not closed")));
check(()=>assert.equal(io.parse(cif,"x.cif").length,3));
console.log(`Structure viewer: ${passed} parser assertions passed; not topology or biological validation.`);

"use strict";
// Original dependency-free perspective projection of supplied coordinates.
// No network requests, modeled coordinates, inferred chemistry or topology.
const StructureIO = (() => {
  const bases = {A:"A",C:"C",G:"G",U:"U",T:"T",DA:"A",DC:"C",DG:"G",DT:"T",DU:"U"};
  function parsePDB(text) {
    const atoms = [];
    let sawModel = false;
    for (const line of text.split(/\r?\n/)) {
      if (line.startsWith("MODEL")) { if (sawModel) break; sawModel = true; }
      if (line.startsWith("ENDMDL")) break;
      if (!/^(ATOM  |HETATM)/.test(line)) continue;
      const alternate = line.slice(16,17).trim();
      if (alternate && alternate !== "A") continue;
      const xyz = [line.slice(30,38),line.slice(38,46),line.slice(46,54)];
      if (line.length < 54 || xyz.some(value=>!value.trim() || !Number.isFinite(Number(value))))
        throw new Error("Missing or invalid coordinates in a PDB atom record.");
      const atom = {name:line.slice(12,16).trim(),resname:line.slice(17,20).trim(),
        chain:line.slice(21,22).trim()||"_",resid:line.slice(22,26).trim(),insertion:line.slice(26,27).trim(),
        x:Number(line.slice(30,38)),y:Number(line.slice(38,46)),z:Number(line.slice(46,54)),
        b:line.slice(60,66).trim()?Number(line.slice(60,66)):NaN,element:line.slice(76,78).trim(),model:"1"};
      if ([atom.x,atom.y,atom.z].every(Number.isFinite)) atoms.push(atom);
    }
    return atoms;
  }
  function cifTokens(text) {
    const out = [];
    let i = 0;
    while (i < text.length) {
      if (/\s/.test(text[i])) {i++;continue;}
      if (text[i] === "#") {while(i<text.length && text[i] !== "\n") i++;continue;}
      if (text[i] === ";" && (i===0 || text[i-1] === "\n")) {
        const start=++i; let end=text.indexOf("\n;",i);
        if (end<0) throw new Error("Unterminated mmCIF multiline value");
        out.push({v:text.slice(start,end),quoted:true});i=end+2;continue;
      }
      if (text[i] === "'" || text[i] === '"') {
        const quote=text[i++],start=i;
        while(i<text.length && !(text[i]===quote && (i+1===text.length || /\s/.test(text[i+1])))) i++;
        if(i===text.length) throw new Error("Unterminated quoted mmCIF value");
        out.push({v:text.slice(start,i++),quoted:true});continue;
      }
      const start=i;while(i<text.length&&!/\s/.test(text[i]))i++;
      out.push({v:text.slice(start,i),quoted:false});
    }
    return out;
  }
  function parseCIF(text) {
    const tokens=cifTokens(text),atoms=[];let selectedModel=null;
    const keyword=(t)=>!t.quoted&&(/^(loop_|stop_|data_|save_|global_)/i.test(t.v)||t.v.startsWith("_"));
    for(let i=0;i<tokens.length;) {
      if(tokens[i].v!=="loop_"||tokens[i].quoted){i++;continue;}
      i++;const cols=[];while(i<tokens.length&&!tokens[i].quoted&&tokens[i].v.startsWith("_"))cols.push(tokens[i++].v);
      if(!cols.length)throw new Error("mmCIF loop has no column names");
      const atomLoop=cols.some(c=>c.startsWith("_atom_site."));
      const indices=Object.fromEntries(cols.map((c,n)=>[c.replace(/^_atom_site\./,""),n]));
      while(i<tokens.length&&!keyword(tokens[i])) {
        if(i+cols.length>tokens.length||tokens.slice(i,i+cols.length).some(keyword))throw new Error("Incomplete mmCIF loop row");
        const row=tokens.slice(i,i+cols.length).map(t=>t.v);i+=cols.length;
        if(!atomLoop)continue;
        const get=(...keys)=>{for(const key of keys){const value=row[indices[key]];if(value!==undefined&&value!=="."&&value!=="?")return value;}return "";};
        const model=get("pdbx_PDB_model_num")||"1";if(selectedModel===null)selectedModel=model;
        if(model!==selectedModel)continue;
        const alt=get("label_alt_id");if(alt&&alt!=="A")continue;
        const atom={name:get("label_atom_id","auth_atom_id"),resname:get("label_comp_id","auth_comp_id"),
          chain:get("label_asym_id","auth_asym_id")||"_",resid:get("label_seq_id","auth_seq_id"),
          insertion:get("pdbx_PDB_ins_code"),x:Number(get("Cartn_x")),y:Number(get("Cartn_y")),z:Number(get("Cartn_z")),
          b:get("B_iso_or_equiv")?Number(get("B_iso_or_equiv")):NaN,element:get("type_symbol"),model};
        if(!get("Cartn_x")||!get("Cartn_y")||!get("Cartn_z"))throw new Error("Missing mmCIF atom coordinates");
        if([atom.x,atom.y,atom.z].every(Number.isFinite))atoms.push(atom);
      }
    }
    return atoms;
  }
  function summarize(atoms) {
    const chains=new Map();
    for(const a of atoms){if(!chains.has(a.chain))chains.set(a.chain,new Map());
      const key=a.resid+":"+a.insertion;if(!chains.get(a.chain).has(key))chains.get(a.chain).set(key,a);}
    return [...chains].map(([chain,residues])=>({chain,residues:residues.size,
      contiguous_one_based:[...residues.values()].every((a,i)=>a.resid===String(i+1)&&!a.insertion),
      sequence:[...residues.values()].map(a=>bases[a.resname]||"?").join("")}));
  }
  function parse(text,name) {
    const atoms=/\.(cif|mmcif)$/i.test(name)||/^\s*data_/m.test(text)?parseCIF(text):parsePDB(text);
    if(!atoms.length)throw new Error("No valid coordinates found. Check the file format and contents.");
    if(atoms.length>300000)throw new Error("More than 300,000 atoms. Provide a smaller structure for preview.");
    return atoms;
  }
  return {parsePDB,parseCIF,parse,summarize,cifTokens};
})();
if(typeof module!=="undefined"&&module.exports)module.exports=StructureIO;

if(typeof document!=="undefined") {
  const $=id=>document.getElementById(id),canvas=$("structure-canvas"),ctx=canvas.getContext("2d");
  const colours=["#57b8ff","#ffb86b","#60d5b6","#d39bfa","#e8db6a","#ff849f"];
  const state={atoms:[],chains:[],visible:new Set(),mapping:null,angleX:.25,angleY:-.35,zoom:1,center:[0,0,0],radius:1,
    highlight:{chain:"C",start:1,end:16},filename:"",points:[],links:[]};
  function status(text){$("structure-status").textContent=text;}
  function applyMapping(){
    if(!state.mapping||!state.atoms.length)return;
    const mapping=state.mapping,chain=mapping.RNA_chain;
    const expected=(mapping.chains||[]).find(c=>c.id===chain),actual=state.chains.find(c=>c.chain===chain);
    if(!expected||!actual||actual.sequence!==expected.sequence){status("The mapping sequence does not match the structure. Verify the chain, then set a region manually.");return;}
    if(!actual.contiguous_one_based){status("The chain sequence matches, but residue numbering is offset or contains insertion codes. Set the region using the structure's residue numbers.");return;}
    $("highlight-chain").value=chain;$("highlight-start").value=mapping.RNA_chain_site_start_1based;
    $("highlight-end").value=mapping.RNA_chain_site_end_1based;highlight();
    status(`Chain ${chain} sequence verified · transcript positions ${mapping.RNA_chain_transcript_start_1based+mapping.RNA_chain_site_start_1based-1}–${mapping.RNA_chain_transcript_start_1based+mapping.RNA_chain_site_end_1based-1} highlighted.`);
  }
  function prepare(){
    const points=[],links=[];
    if($("representation").value==="atoms") {
      const heavy=state.atoms.filter(a=>!/^H|D$/.test(a.element||a.name)&&state.visible.has(a.chain));
      const stride=Math.max(1,Math.ceil(heavy.length/15000));
      heavy.forEach((a,i)=>{if(i%stride===0)points.push(a);});
    } else {
      const residues=new Map();
      for(const a of state.atoms){if(!state.visible.has(a.chain))continue;const key=a.chain+":"+a.resid+":"+a.insertion;
        if(!residues.has(key))residues.set(key,[]);residues.get(key).push(a);}
      const previous=new Map();
      for(const items of residues.values()) {
        const a=items.find(a=>a.name==="P")||items.find(a=>a.name==="C4'"||a.name==="C4*")||items.find(a=>a.name==="CA")||items[0];
        points.push(a);const before=previous.get(a.chain);
        if(before&&Number(a.resid)===Number(before.resid)+1&&Math.hypot(a.x-before.x,a.y-before.y,a.z-before.z)<12)links.push([before,a]);
        previous.set(a.chain,a);
      }
    }
    state.points=points;state.links=links;draw();
  }
  function confidenceColor(b){return b>=90?"#367fff":b>=70?"#72c7f5":b>=50?"#ffe277":"#ef9257";}
  function colour(a){
    if(a.chain===state.highlight.chain&&Number(a.resid)>=state.highlight.start&&Number(a.resid)<=state.highlight.end)return "#ff6286";
    if($("colour-mode").value==="confidence"&&$("plddt-confirm").checked&&Number.isFinite(a.b)&&a.b>=0&&a.b<=100)return confidenceColor(a.b);
    return colours[state.chains.findIndex(c=>c.chain===a.chain)%colours.length]||"#ffffff";
  }
  function project(a,w,h){
    let x=a.x-state.center[0],y=a.y-state.center[1],z=a.z-state.center[2];
    const cy=Math.cos(state.angleY),sy=Math.sin(state.angleY),cx=Math.cos(state.angleX),sx=Math.sin(state.angleX);
    const xx=x*cy+z*sy,zz=-x*sy+z*cy,yy=y*cx-zz*sx,depth=y*sx+zz*cx;
    const camera=state.radius*4,scale=Math.min(w,h)*.37*state.zoom/state.radius, perspective=camera/Math.max(camera*.25,camera+depth);
    return{x:w/2+xx*scale*perspective,y:h/2-yy*scale*perspective,z:depth,r:Math.max(1.4,Math.min(6,2.6*perspective))};
  }
  function draw(){
    const dpr=window.devicePixelRatio||1,w=canvas.clientWidth||900,h=canvas.clientHeight||570;
    if(canvas.width!==Math.round(w*dpr)||canvas.height!==Math.round(h*dpr)){canvas.width=Math.round(w*dpr);canvas.height=Math.round(h*dpr);}
    ctx.setTransform(dpr,0,0,dpr,0,0);ctx.fillStyle="#102a34";ctx.fillRect(0,0,w,h);
    if(!state.atoms.length){ctx.fillStyle="#aac9d2";ctx.font="16px sans-serif";ctx.textAlign="center";ctx.fillText("Open a PDB or mmCIF file",w/2,h/2);return;}
    ctx.lineWidth=2;for(const [a,b]of state.links){const p=project(a,w,h),q=project(b,w,h);ctx.strokeStyle=colour(a);ctx.beginPath();ctx.moveTo(p.x,p.y);ctx.lineTo(q.x,q.y);ctx.stroke();}
    const projected=state.points.map(a=>({a,...project(a,w,h)})).sort((a,b)=>b.z-a.z);
    for(const p of projected){ctx.fillStyle=colour(p.a);ctx.beginPath();ctx.arc(p.x,p.y,p.r,0,2*Math.PI);ctx.fill();}
    ctx.fillStyle="#b2cdd4";ctx.font="12px sans-serif";ctx.textAlign="left";
    ctx.fillText(`${state.filename} · ${state.points.length} display points · coordinates in Å`,18,h-18);
  }
  function reset(){state.angleX=.25;state.angleY=-.35;state.zoom=1;draw();}
  function highlight(){const start=Number($("highlight-start").value),end=Number($("highlight-end").value);
    if(!Number.isInteger(start)||!Number.isInteger(end)||start<1||end<start){status("Enter a valid positive-integer residue range.");return;}
    state.highlight={chain:$("highlight-chain").value.trim(),start,end};draw();}
  async function loadCoordinates(event){const file=event.target.files[0];if(!file)return;
    if(file.size>80*1024*1024){status("File exceeds 80 MB. Crop the structure before preview.");return;}
    try {const atoms=StructureIO.parse(await file.text(),file.name);state.atoms=atoms;state.filename=file.name;state.chains=StructureIO.summarize(atoms);state.visible=new Set(state.chains.map(c=>c.chain));
      state.center=[0,1,2].map(axis=>atoms.reduce((sum,a)=>sum+[a.x,a.y,a.z][axis],0)/atoms.length);
      state.radius=Math.max(1,...atoms.map(a=>Math.hypot(a.x-state.center[0],a.y-state.center[1],a.z-state.center[2])).filter((_,i)=>i%Math.max(1,Math.floor(atoms.length/10000))===0));
      $("chain-list").replaceChildren();state.chains.forEach((c,i)=>{const label=document.createElement("label");label.className="chain-row";
        const box=document.createElement("input");box.type="checkbox";box.checked=true;box.addEventListener("change",()=>{if(box.checked)state.visible.add(c.chain);else state.visible.delete(c.chain);prepare();});
        const swatch=document.createElement("span");swatch.className="swatch";swatch.style.background=colours[i%colours.length];label.append(box,swatch,document.createTextNode(`${c.chain} · ${c.residues} residues`));$("chain-list").append(label);});
      $("coordinate-info").textContent=`${atoms.length.toLocaleString()} atoms\n${state.chains.length} chains\nFirst model / preferred alternate location`;
      status(`${file.name} loaded · pink marks the selected region.`);reset();prepare();applyMapping();
    }catch(error){status(`Could not open file: ${error.message}`);}}
  $("coordinate-file").addEventListener("change",loadCoordinates);
  $("mapping-file").addEventListener("change",async event=>{const file=event.target.files[0];if(!file)return;
    try{if(file.size>1024*1024)throw new Error("Mapping file exceeds 1 MB");state.mapping=JSON.parse(await file.text());applyMapping();}catch(error){status(`Could not load mapping: ${error.message}`);}});
  $("representation").addEventListener("change",prepare);$("colour-mode").addEventListener("change",()=>{if($("colour-mode").value==="confidence"&&!$("plddt-confirm").checked)status("Confirm that B-factor stores pLDDT. Chain colours are used until confirmed.");draw();});
  $("plddt-confirm").addEventListener("change",draw);$("apply-highlight").addEventListener("click",highlight);$("reset-camera").addEventListener("click",reset);
  let pointer=null;canvas.addEventListener("pointerdown",e=>{pointer={x:e.clientX,y:e.clientY};canvas.setPointerCapture(e.pointerId);});
  canvas.addEventListener("pointermove",e=>{if(!pointer)return;state.angleY+=(e.clientX-pointer.x)*.008;state.angleX+=(e.clientY-pointer.y)*.008;pointer={x:e.clientX,y:e.clientY};draw();});
  canvas.addEventListener("pointerup",()=>pointer=null);canvas.addEventListener("wheel",e=>{e.preventDefault();state.zoom=Math.max(.2,Math.min(8,state.zoom*Math.exp(-e.deltaY*.001)));draw();},{passive:false});
  $("save-image").addEventListener("click",()=>{const link=document.createElement("a");link.href=canvas.toDataURL("image/png");link.download="structure_coordinate_preview.png";link.click();});
  window.addEventListener("resize",draw);draw();
}

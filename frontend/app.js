// app.js — όλη η συμπεριφορά (JavaScript) της εφαρμογής.
// Μιλάει με το backend μέσω fetch() και δείχνει τα αποτελέσματα στη σελίδα.

let charts = {};       // κρατάμε τα Chart.js αντικείμενα για να τα κάνουμε destroy πριν ξανασχεδιάσουμε
let sessionId = "";    // το id του εκπαιδευμένου μοντέλου ΑΥΤΟΥ του χρήστη (Tab 3)
const MAX_MB = 5;      // όριο μεγέθους αρχείου (πρέπει να ταιριάζει με το backend)

// ---------- Βοηθητικές ----------
async function postJSON(url, body){
  const res = await fetch(url, {method:"POST",
    headers:{"Content-Type":"application/json"}, body:JSON.stringify(body)});
  return res.json();
}
// Ελέγχει το μέγεθος αρχείου πριν το ανεβάσουμε (καλύτερο UX από το να το στείλουμε τζάμπα).
function fileTooBig(file){
  if(file && file.size > MAX_MB*1024*1024){
    return `Το αρχείο είναι πολύ μεγάλο (${(file.size/1048576).toFixed(1)} MB, όριο ${MAX_MB} MB).`;
  }
  return "";
}
function esc(s){ return String(s).replace(/[&<>]/g, c => ({'&':'&amp;','<':'&lt;','>':'&gt;'}[c])); }
function errBox(m){ return `<div class="err">⚠️ ${esc(m)}</div>`; }
function loading(el){ el.innerHTML = `<p style="margin-top:14px;color:#6b7280"><span class="spinner"></span>Επεξεργασία...</p>`; }
function badge(txt){ return `<span class="badge ${txt}">${txt}</span>`; }

function table(rows, badgeCols={}){
  if(!rows || !rows.length) return "<p class='sub'>Κανένα αποτέλεσμα.</p>";
  const cols = Object.keys(rows[0]);
  let h = "<table><thead><tr>"+cols.map(c=>`<th>${esc(c)}</th>`).join("")+"</tr></thead><tbody>";
  h += rows.map(r=>"<tr>"+cols.map(c=>{
    let v = r[c] === null ? "—" : r[c];
    if(badgeCols[c]) return `<td>${badge(v)}</td>`;
    return `<td>${esc(v)}</td>`;
  }).join("")+"</tr>").join("");
  return h+"</tbody></table>";
}
function makeChart(id, cfg){
  if(charts[id]) charts[id].destroy();
  charts[id] = new Chart(document.getElementById(id), cfg);
}

// ---------- Tabs ----------
document.querySelectorAll(".tab").forEach(b=>{
  b.onclick=()=>{
    document.querySelectorAll(".tab").forEach(x=>x.classList.remove("active"));
    document.querySelectorAll(".panel").forEach(p=>p.classList.remove("active"));
    b.classList.add("active");
    document.getElementById(b.dataset.tab).classList.add("active");
  };
});

// ---------- Drag & drop upload ----------
function wireDrop(dropId, inputId, nameId){
  const drop=document.getElementById(dropId), input=document.getElementById(inputId),
        name=document.getElementById(nameId);
  drop.onclick=()=>input.click();
  input.onchange=()=>{ if(input.files[0]) name.textContent="✓ "+input.files[0].name; };
  ["dragover","dragenter"].forEach(e=>drop.addEventListener(e,ev=>{ev.preventDefault();drop.classList.add("over");}));
  ["dragleave","drop"].forEach(e=>drop.addEventListener(e,ev=>{ev.preventDefault();drop.classList.remove("over");}));
  drop.addEventListener("drop",ev=>{ input.files=ev.dataTransfer.files;
    if(input.files[0]) name.textContent="✓ "+input.files[0].name; });
}
wireDrop("p1drop","p1file","p1fname");
wireDrop("p3drop","p3file","p3fname");

// ============================================================
//  TAB 1 — Ανάλυση
// ============================================================
async function analyze(){
  const out=document.getElementById("p1out");
  const file=document.getElementById("p1file").files[0];
  const big=fileTooBig(file); if(big){ out.innerHTML=errBox(big); return; }
  loading(out);
  const fasta=file ? await file.text() : document.getElementById("p1text").value;
  const d=await postJSON("/api/analyze",{fasta});
  if(d.error){ out.innerHTML=errBox(d.error); return; }

  const rep=d.report;
  const f=rep[0];
  let html = `<div class="stat-row">
      <div class="stat"><div class="k">Ακολουθίες</div><div class="v">${rep.length}</div></div>
      <div class="stat"><div class="k">Μήκος (1η)</div><div class="v">${f.length}</div></div>
      <div class="stat"><div class="k">GC (1η)</div><div class="v">${(f.gc_fraction*100).toFixed(1)}%</div></div>
      <div class="stat"><div class="k">Μεγ. ORF</div><div class="v">${f.longest_orf_aa}</div></div>
    </div>`;
  html += `<div class="chartbox"><canvas id="gcChart" height="110"></canvas></div>`;
  html += `<div class="chartbox"><canvas id="baseChart" height="110"></canvas></div>`;
  // ΕΠΙΠΕΔΟ 2: reverse complement + πρωτεΐνη της 1ης ακολουθίας
  if(d.reverse_complement){
    html += `<h3>Reverse complement (1η ακολουθία, έως 120 βάσεις)</h3>
             <div class="seqbox">${esc(d.reverse_complement)}</div>`;
  }
  if(d.protein){
    html += `<h3>Πρωτεΐνη (μετάφραση έως το 1ο stop)</h3>
             <div class="seqbox">${esc(d.protein) || "—"}</div>`;
  }
  html += `<h3>Πλήρης πίνακας</h3>`+table(rep);
  html += `<h3>ORFs στην «${esc(d.first_id)}»</h3>`+table(d.orfs);
  out.innerHTML=html;

  const ids=rep.map(r=>r.id);
  makeChart("gcChart",{type:"bar",
    data:{labels:ids,datasets:[{label:"GC fraction",data:rep.map(r=>r.gc_fraction),
      backgroundColor:"#4f46e5",borderRadius:6}]},
    options:{plugins:{title:{display:true,text:"GC content ανά ακολουθία"},legend:{display:false}},
      scales:{y:{beginAtZero:true,max:1}}}});
  makeChart("baseChart",{type:"bar",
    data:{labels:ids,datasets:[
      {label:"A",data:rep.map(r=>r.a),backgroundColor:"#60a5fa"},
      {label:"T",data:rep.map(r=>r.t),backgroundColor:"#a78bfa"},
      {label:"G",data:rep.map(r=>r.g),backgroundColor:"#f59e0b"},
      {label:"C",data:rep.map(r=>r.c),backgroundColor:"#34d399"}]},
    options:{plugins:{title:{display:true,text:"Σύσταση βάσεων ανά ακολουθία"}},
      scales:{x:{stacked:true},y:{stacked:true}}}});
}

// ============================================================
//  TAB 2 — Μεταλλάξεις
// ============================================================
async function mutations(){
  const out=document.getElementById("p2out"); loading(out);
  const d=await postJSON("/api/mutations",{
    reference:document.getElementById("p2ref").value,
    variant:document.getElementById("p2var").value,
    align:document.getElementById("p2align").checked});
  if(d.error){ out.innerHTML=errBox(d.error); return; }

  // ΕΠΙΠΕΔΟ 2: alignment mode (άνισα μήκη / indels)
  if(d.mode==="align"){
    out.innerHTML =
      `<div class="stat-row">
         <div class="stat"><div class="k">Mismatches</div><div class="v">${d.n_mismatch}</div></div>
         <div class="stat"><div class="k">Gaps (indels)</div><div class="v">${d.n_gap}</div></div>
       </div>
       <h3>Ευθυγράμμιση</h3>
       <div class="aln">${alignHTML(d.aligned_ref,d.aligned_var)}</div>`+
      table(d.diffs,{kind:true});
    return;
  }
  // SNP mode
  if(d.hamming===0){ out.innerHTML=`<div class="ok">✓ Καμία διαφορά — οι ακολουθίες είναι ίδιες.</div>`; return; }
  out.innerHTML =
    `<div class="stat-row"><div class="stat"><div class="k">Μεταλλάξεις</div>
       <div class="v">${d.hamming}</div></div></div>`+
    table(d.snps,{ts_tv:true,effect:true});
}
// Χρωματίζει την ευθυγράμμιση: κόκκινο mismatch, κίτρινο gap.
function alignHTML(a,b){
  let l1="",l2="",l3="";
  for(let i=0;i<a.length;i++){
    const x=a[i],y=b[i];
    if(x===y){ l1+=x; l2+=" "; l3+=y; }
    else if(x==="-"||y==="-"){ l1+=`<span class="gp">${x}</span>`; l2+=`<span class="gp">-</span>`; l3+=`<span class="gp">${y}</span>`; }
    else { l1+=`<span class="mm">${x}</span>`; l2+=`<span class="mm">*</span>`; l3+=`<span class="mm">${y}</span>`; }
  }
  return `ref: ${l1}\n     ${l2}\nvar: ${l3}`;
}

// ============================================================
//  TAB 3 — Ταξινόμηση
// ============================================================
async function train(){
  const out=document.getElementById("p3trainout");
  const file=document.getElementById("p3file").files[0];
  if(!file){ out.innerHTML=errBox("Διάλεξε πρώτα ένα αρχείο FASTA."); return; }
  const big=fileTooBig(file); if(big){ out.innerHTML=errBox(big); return; }
  loading(out);
  const d=await postJSON("/api/train",{fasta:await file.text()});
  if(d.error){ out.innerHTML=errBox(d.error); return; }
  sessionId = d.session_id;   // αποθηκεύουμε το id ΤΟΥ δικού μας μοντέλου
  out.innerHTML=`<div class="ok">✓ Εκπαιδεύτηκε με <b>${d.n_sequences}</b> ακολουθίες ·
    κατηγορίες: ${d.classes.map(c=>`<b>${esc(c)}</b>`).join(", ")}</div>`;
}
async function predict(){
  const out=document.getElementById("p3predout");
  if(!sessionId){ out.innerHTML=errBox("Κάνε πρώτα εκπαίδευση (Βήμα 1)."); return; }
  loading(out);
  const d=await postJSON("/api/predict",
    {sequence:document.getElementById("p3seq").value, session_id:sessionId});
  if(d.error){ out.innerHTML=errBox(d.error); return; }
  let bars="";
  for(const [k,v] of Object.entries(d.probabilities)){
    bars+=`<div class="prob"><div class="lab"><span>Κατηγορία ${esc(k)}</span><span>${(v*100).toFixed(0)}%</span></div>
      <div class="bar"><div class="fill" style="width:${v*100}%"></div></div></div>`;
  }
  out.innerHTML=`<div class="predbox"><div class="cap">Προβλεπόμενη κατηγορία</div>
    <div class="big">${esc(d.label)}</div></div>${bars}`;
}

// ============================================================
//  ΕΠΙΠΕΔΟ 2 — NCBI fetch & Ιστορικό
// ============================================================
async function fetchNcbi(){
  const acc=document.getElementById("p1ncbi").value.trim();
  const out=document.getElementById("p1out");
  if(!acc){ out.innerHTML=errBox("Δώσε ένα NCBI accession ID (π.χ. NM_000410)."); return; }
  loading(out);
  const d=await postJSON("/api/fetch_ncbi",{accession:acc});
  if(d.error){ out.innerHTML=errBox(d.error); return; }
  document.getElementById("p1text").value=d.fasta;   // γεμίζει το πεδίο FASTA
  out.innerHTML=`<div class="ok">✓ Κατέβηκε το «${esc(acc)}». Πάτα «Ανάλυσε».</div>`;
}

async function loadHistory(){
  const out=document.getElementById("p4out"); loading(out);
  const res=await fetch("/api/history"); const d=await res.json();
  if(!d.history || !d.history.length){ out.innerHTML="<p class='sub'>Δεν υπάρχει ιστορικό ακόμη.</p>"; return; }
  out.innerHTML=table(d.history);
}
// Φόρτωσε το ιστορικό αυτόματα όταν ανοίγει το tab.
document.querySelector('[data-tab="p4"]').addEventListener("click", loadHistory);

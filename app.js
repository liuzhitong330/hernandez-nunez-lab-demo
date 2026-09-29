(function () {
"use strict";
const $=id=>document.getElementById(id), NS="http://www.w3.org/2000/svg";
document.querySelectorAll('select,button').forEach(e=>e.disabled=true);
const color={ink:"#222",muted:"#999",accent:"#1f7a8c",selected:"#c07a2a"};
const S={fish:0,cell:0,event:0,budget:10,score:"nutrient",yaw:.45,pitch:.58,records:[],projected:[],cache:new Map(),seq:0};
const finite=x=>x!==null && Number.isFinite(x);
const fmt=(n,d=1)=>finite(n)?n.toFixed(d):"not estimable";
const pct=n=>finite(n)?(n>0?"+":"")+n.toFixed(1)+"%":"not estimable";
const text=(id,s)=>$(id).textContent=s;
const elem=(tag,attrs={},s)=>{const e=document.createElementNS(NS,tag);Object.entries(attrs).forEach(([k,v])=>e.setAttribute(k,v));if(s!==undefined)e.textContent=s;return e;};
const f=()=>S.records[S.fish], c=()=>f().cells[S.cell];
const metric=(cell,key=S.score)=>cell.metrics[key+"_error_reduction_pct"];
const shortlisted=r=>r.cells.filter(x=>x.discovery_rank<=S.budget);
const good=(cell,key="nutrient")=>finite(metric(cell,key)) && metric(cell,key)>0 && cell.metrics.full_r>0;
function axes(svg,xLabel,yLabel,{w=620,h=300,left=58,right=18,top=18,bottom=46}={}){
 svg.replaceChildren();const g=elem("g");svg.append(g);g.append(elem("path",{d:"M"+left+" "+top+"V"+(h-bottom)+"H"+(w-right),fill:"none",stroke:"#aaa"}));
 g.append(elem("text",{x:(left+w-right)/2,y:h-5,"text-anchor":"middle","font-size":12,fill:"#555"},xLabel));
 g.append(elem("text",{x:left,y:12,"font-size":12,fill:"#555"},yLabel));
 return {g,left,right,top,bottom,w,h};
}
function range(values){const a=values.filter(finite);if(!a.length)return [-1,1];let lo=Math.min(...a),hi=Math.max(...a);const pad=(hi-lo||1)*.08;return [lo-pad,hi+pad];}
function pathSeries(g,time,values,idx,x,y,stroke,dash){let d="",pen=false;idx.forEach(i=>{if(!finite(values[i])){pen=false;return;}d+=(pen?"L":"M")+x(time[i]).toFixed(2)+","+y(values[i]).toFixed(2);pen=true;});g.append(elem("path",{d,fill:"none",stroke,"stroke-width":1.5,...(dash?{"stroke-dasharray":"5 4"}:{})}));}
function labelEvent(e,r=f()){return (e.kind==="nutrient"?"Position "+r.target_position:e.kind==="control"?"Reference position "+r.reference_position:"Other UV")+" · "+e.time.toFixed(1)+" s · "+(e.role==="holdout"?"reserved":e.role==="discovery"?"discovery":e.role);}
function options(select,values){select.replaceChildren();values.forEach(([v,l])=>{const e=document.createElement("option");e.value=v;e.textContent=l;select.append(e);});}
function setupRecording(){
 const r=f();const order=r.cells.map((v,i)=>[i,v]).sort((a,b)=>a[1].discovery_rank-b[1].discovery_rank);
 options($("cell"),order.map(([i,v])=>[i,v.id+" · rank "+v.discovery_rank]));S.cell=order[0][0];$("cell").value=S.cell;
 options($("event"),r.events.filter(e=>e.complete).map(e=>[r.events.indexOf(e),labelEvent(e)]));
 S.event=Math.max(0,r.events.findIndex(e=>e.role==="holdout"&&e.kind==="nutrient"));$("event").value=S.event;
 text("protocol-note",r.protocol.name+". Fixed source-index sampling; "+r.counts.sampled_cells_before_camera_baseline_qc+" sampled ROIs, "+r.counts.analyzed_cells+" analyzed after QC. "+r.protocol.analysis_binning+". Response scoring starts 2 s after a pulse and stops at the next pulse or 55 s. Training excludes reserved windows and a 55 s buffer. Fixed regularization (ridge 2, smoothness 20), training-only scaling. "+r.protocol.rank_definition.replaceAll('nutrient','position-code')+". Source metadata do not verify organ identity for each position code; comparisons remain code-specific.");
 text("map-caption",r.label+" · "+r.counts.analyzed_cells+" sampled ROIs. Acquisition voxel indices; Z is expanded ×4 for display, not a physical-distance scale. No atlas or cell-type labels inferred.");
 renderAll();
}
function project(){
 const coords=f().cells.map(x=>x.xyz),mins=[0,1,2].map(k=>Math.min(...coords.map(a=>a[k]))),maxs=[0,1,2].map(k=>Math.max(...coords.map(a=>a[k])));
 const center=mins.map((x,k)=>(x+maxs[k])/2),span=Math.max(maxs[0]-mins[0],maxs[1]-mins[1],4*(maxs[2]-mins[2]))||1;
 return coords.map((v,i)=>{const x=(v[0]-center[0])/span*470,y=(v[1]-center[1])/span*470,z=(v[2]-center[2])/span*470*4;
 const a=x*Math.cos(S.yaw)+z*Math.sin(S.yaw),b=-x*Math.sin(S.yaw)+z*Math.cos(S.yaw);
 return {x:300+a,y:165+y*Math.cos(S.pitch)-b*Math.sin(S.pitch),z:y*Math.sin(S.pitch)+b*Math.cos(S.pitch),i};});
}
function drawBrain(){
 const canvas=$("brain"),ctx=canvas.getContext("2d");ctx.setTransform(2,0,0,2,0,0);ctx.clearRect(0,0,600,330);
 S.projected=project();const points=[...S.projected].sort((a,b)=>a.z-b.z);
 points.forEach(p=>{const selected=p.i===S.cell,top=f().cells[p.i].discovery_rank<=S.budget;ctx.beginPath();ctx.arc(p.x,p.y,selected?5:top?3:1.6,0,Math.PI*2);ctx.fillStyle=selected?color.selected:top?color.accent:"#a5b0b477";ctx.fill();});
 const p=S.projected[S.cell];ctx.strokeStyle=color.selected;ctx.beginPath();ctx.arc(p.x,p.y,8,0,Math.PI*2);ctx.stroke();
 ctx.fillStyle=color.ink;ctx.font="12px system-ui";ctx.fillText("Cell "+c().id+" · ("+c().xyz.map(n=>n.toFixed(1)).join(", ")+")",12,20);
 ctx.fillStyle="#777";ctx.font="11px system-ui";ctx.fillText("X / Y / Z: acquisition indices • depth display ×4",12,318);
}
async function getTrace(){
 const r=f(),cell=c(),seq=++S.seq;
 $("trace").replaceChildren();$("signals").replaceChildren();text("trial-note","Loading this cell's measured trial…");
 const key=r.recording_id+"/"+cell.shard;
 try {
  if(!S.cache.has(key)){const response=await fetch(cell.shard+"?v=2");if(!response.ok)throw new Error("Trace file unavailable");const decoded=cell.shard.endsWith('.gz')?await new Response(response.body.pipeThrough(new DecompressionStream('gzip'))).json():await response.json();S.cache.set(key,decoded);}
  if(seq!==S.seq)return;
  const data=S.cache.get(key)[String(cell.id)]; if(!data)throw new Error("Cell trace absent");
  drawTrace(data);
 }catch(e){if(seq===S.seq)text("trial-note","Trace unavailable: "+e.message+". Other analysis remains visible; reload to retry.");}
}
function drawTrace(data){
 const r=f(),ev=r.events[S.event],base=$("baseline").value;
 const low=ev.time-10,high=Math.min(ev.time+55,r.time.at(-1)),idx=r.time.map((v,i)=>i).filter(i=>r.time[i]>=low&&r.time[i]<=high);
 const all=[data.trace,...Object.values(data.predictions)].flatMap(arr=>idx.map(i=>arr[i])),[lo,hi]=range(all);
 const a=axes($("trace"),"Seconds from selected pulse","ΔF/F · source-derived",{h:250}),{g}=a;
 const x=t=>58+(t-low)/(high-low)*544,y=v=>204-(v-lo)/(hi-lo)*180;
 g.append(elem("rect",{x:x(ev.time),y:20,width:Math.max(0,x(ev.evaluation_end)-x(ev.time)),height:184,fill:"#1f7a8c0c"}));
 [lo,(lo+hi)/2,hi].forEach(v=>g.append(elem("text",{x:51,y:y(v)+4,"text-anchor":"end","font-size":11,fill:"#777"},fmt(v,2))));
 [-10,0,20,40].forEach(v=>{const t=ev.time+v;if(t>=low&&t<=high){g.append(elem("text",{x:x(t),y:225,"text-anchor":"middle","font-size":11,fill:"#777"},v));}});
 r.events.filter(e=>e.time>=low&&e.time<=high).forEach(e=>g.append(elem("line",{x1:x(e.time),x2:x(e.time),y1:20,y2:204,stroke:"#aaa","stroke-dasharray":"2 4"})));
 pathSeries(g,r.time,data.predictions[base],idx,x,y,"#999",true);pathSeries(g,r.time,data.predictions.FULL,idx,x,y,"#222");pathSeries(g,r.time,data.trace,idx,x,y,color.accent);
 const sg=$("signals");sg.replaceChildren();const gg=elem("g");sg.append(gg);const [ml,mh]=range(idx.map(i=>r.motor[i]));
 gg.append(elem("text",{x:5,y:15,"font-size":11,fill:"#555"},"Motor proxy"));pathSeries(gg,r.time,r.motor,idx,x,v=>48-(v-ml)/(mh-ml)*25,"#666");
 gg.append(elem("text",{x:5,y:74,"font-size":11,fill:"#555"},"Visual input"));pathSeries(gg,r.time,r.visual,idx,x,v=>95-v*24,color.accent);
 const score=c().trial_metrics.find(z=>z.event_id===ev.id);
 text("trial-note",labelEvent(ev)+". Shaded region ends at "+(ev.evaluation_end-ev.time).toFixed(1)+" s; first 2 s excluded from scoring. "+(score?"Full-model r "+fmt(score.full_r,2)+"; selected-site error reduction "+pct(score.nutrient_error_reduction_pct)+".":"Discovery/example window, not a reserved score.")+" "+(r.visual_available?"Visual input varies in this recording.":"No visual modulation in this recording; no visual-effect claim."));
}
function renderCell(){
 const r=f(),cell=c(),m=cell.metrics;
 text("metric-one",pct(m.nutrient_error_reduction_pct));text("metric-two",pct(m.multimodal_error_reduction_pct));text("metric-three",fmt(m.full_r,2));
 text("sensitivity","Slow-drift sensitivity: selected-position gain without drift adjustment "+pct(cell.original_no_drift.nutrient_error_reduction_pct)+"; with the common cubic trend "+pct(m.nutrient_error_reduction_pct)+". "+(cell.alternative_reference?"Alternate reference code "+cell.alternative_reference.reference_position+": "+pct(cell.alternative_reference.metrics.nutrient_error_reduction_pct)+", full r "+fmt(cell.alternative_reference.metrics.full_r,2)+". Same ROI, no validation-based reranking.":"No alternate-reference analysis for this recording."));
 const referenceStable=!cell.alternative_reference || (cell.alternative_reference.metrics.nutrient_error_reduction_pct>0 && cell.alternative_reference.metrics.full_r>0);
 text("readout","Cell "+cell.id+" was discovery rank "+cell.discovery_rank+". On reserved selected-position events, the mean per-trial error reduction from adding its position code is "+pct(m.nutrient_error_reduction_pct)+". Mean per-trial full-model Pearson r is "+fmt(m.full_r,2)+". "+(!referenceStable?"The alternate-reference check fails: this candidate is sensitive to comparison choice.":good(cell)?"Both primary descriptive checks are positive; independent re-recording is still needed.":"This does not support immediate prioritization: error gain and response shape must both be inspected.")+" Position-to-organ identity is not verified.");
 const body=$("trial-table").querySelector("tbody");body.replaceChildren();
 cell.trial_metrics.forEach(v=>{const e=r.events.find(e=>e.id===v.event_id);const tr=document.createElement("tr");[labelEvent(e),pct(v.nutrient_error_reduction_pct),pct(v.multimodal_error_reduction_pct),fmt(v.full_r,2)].forEach(t=>{const td=document.createElement("td");td.textContent=t;tr.append(td);});body.append(tr);});
 text("plan","First resolve the acquisition metadata: confirm what position codes "+r.reference_position+" and "+r.target_position+" illuminated and whether light actually reached those sites. "+(!referenceStable?"For this ROI, the priority is an interleaved, order-balanced recording: the earlier and later reference blocks lead to different conclusions. ":"")+"Establish anatomical registration or molecular identity before identifying a corresponding candidate in another animal; native coordinates alone are insufficient. Re-record with matched UV controls and a motor measurement. "+(good(cell)&&referenceStable?"A repeatable position-linked response would motivate a controlled perturbation, not establish a pathway by itself.":"Prioritize replication/QC over perturbation; compare response shape as well as average error.")+" Separate nutrient-dependent responses from site/light effects with suitable chemical controls. Vary movement/visual conditions independently where feasible; do not interpret motor adjustment as proof of a direct visceral pathway.");
 drawBrain();renderScatter();getTrace();
}
function renderScatter(){
 const r=f(),a=axes($("validation"),"Discovery position-code error reduction (%)","Reserved "+(S.score==="nutrient"?"position-code":"motor / available visual")+" gain (%)"),{g}=a;
 // symlog keeps large negative failures visible without clipping the main cloud.
 const trans=x=>Math.sign(x)*Math.log1p(Math.abs(x)),inv=x=>Math.sign(x)*Math.expm1(Math.abs(x));
 const points=r.cells.filter(v=>finite(v.discovery_score)&&finite(metric(v)));
 const [xl,xh]=range(points.map(v=>trans(v.discovery_score))),[yl,yh]=range(points.map(v=>trans(metric(v))));
 const x=v=>58+(trans(v)-xl)/(xh-xl)*544,y=v=>254-(trans(v)-yl)/(yh-yl)*236;
 [0,10,50,-10,-50,-200].forEach(v=>{if(trans(v)>=xl&&trans(v)<=xh){g.append(elem("text",{x:x(v),y:275,"text-anchor":"middle","font-size":10,fill:"#777"},v));}if(trans(v)>=yl&&trans(v)<=yh)g.append(elem("text",{x:51,y:y(v)+4,"text-anchor":"end","font-size":10,fill:"#777"},v));});
 if(yl<=0&&yh>=0)g.append(elem("line",{x1:58,x2:602,y1:y(0),y2:y(0),stroke:"#aaa","stroke-dasharray":"3 4"}));
 points.forEach(v=>{const selected=v.id===c().id,top=v.discovery_rank<=S.budget;const dot=elem("circle",{cx:x(v.discovery_score),cy:y(metric(v)),r:selected?5:top?3.2:1.7,fill:selected?color.selected:top?color.accent:"#abb4b488",class:"plot-point"});
 const title=elem("title",{},"Cell "+v.id+" · discovery "+pct(v.discovery_score)+" · reserved "+pct(metric(v)));dot.append(title);
 dot.addEventListener("click",()=>{S.cell=r.cells.findIndex(z=>z.id===v.id);$("cell").value=S.cell;renderCell();});g.append(dot);});
 const list=shortlisted(r),n=list.filter(v=>good(v,S.score)).length;
 text("shortlist-note",n+" / "+list.length+" discovery-ranked candidates have both positive reserved "+(S.score==="nutrient"?"position-code":"motor / available visual")+" error reduction and positive full-model correlation. These are descriptive checks, not validated biological targets. Axes use a signed log transform so large failures remain visible. Ranking stays fixed to the discovery position-code score in both comparisons.");
 $("cohort").querySelector("thead th:last-child").textContent="Shortlist positive on reserved "+(S.score==="nutrient"?"position-code":"motor / visual")+" checks";
 const body=$("cohort").querySelector("tbody");body.replaceChildren();
 S.records.forEach(rr=>{const list=shortlisted(rr),tr=document.createElement("tr");[rr.label,rr.counts.analyzed_cells,rr.counts.validation_nutrient_events+" / "+rr.counts.validation_control_events,list.filter(v=>good(v,S.score)).length+" / "+list.length].forEach(t=>{const td=document.createElement("td");td.textContent=t;tr.append(td);});body.append(tr);});
}
function renderAll(){renderCell();}
function download(name,content,mime){const a=document.createElement("a"),url=URL.createObjectURL(new Blob([content],{type:mime}));a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);}
$("fish").addEventListener("change",()=>{S.fish=+$("fish").value;setupRecording();});
$("cell").addEventListener("change",()=>{S.cell=+$("cell").value;renderCell();});
$("event").addEventListener("change",()=>{S.event=+$("event").value;getTrace();});
$("baseline").addEventListener("change",getTrace);
$("budget").addEventListener("change",()=>{S.budget=+$("budget").value;renderAll();});

$("score").addEventListener("change",()=>{S.score=$("score").value;renderAll();});
$("reset-view").addEventListener("click",()=>{S.yaw=.45;S.pitch=.58;drawBrain();});
$("view-xy").addEventListener("click",()=>{S.yaw=0;S.pitch=0;drawBrain();});
$("view-xz").addEventListener("click",()=>{S.yaw=0;S.pitch=Math.PI/2;drawBrain();});
let drag=null;
$("brain").addEventListener("pointerdown",e=>{drag={x:e.clientX,y:e.clientY,total:0};$("brain").setPointerCapture(e.pointerId);});
$("brain").addEventListener("pointermove",e=>{if(!drag)return;const dx=e.clientX-drag.x,dy=e.clientY-drag.y;drag.total+=Math.abs(dx)+Math.abs(dy);S.yaw+=dx*.009;S.pitch+=dy*.009;drag.x=e.clientX;drag.y=e.clientY;drawBrain();});
$("brain").addEventListener("pointerup",e=>{if(drag&&drag.total<6){const b=$("brain").getBoundingClientRect(),x=(e.clientX-b.left)/b.width*600,y=(e.clientY-b.top)/b.height*330;const p=S.projected.reduce((a,b)=>Math.hypot(a.x-x,a.y-y)<Math.hypot(b.x-x,b.y-y)?a:b);if(Math.hypot(p.x-x,p.y-y)<15){S.cell=p.i;$("cell").value=S.cell;renderCell();}}drag=null;});
$("brain").addEventListener("pointercancel",()=>{drag=null;});
$("brain").addEventListener("keydown",e=>{if(["ArrowLeft","ArrowRight","ArrowUp","ArrowDown"].includes(e.key)){e.preventDefault();S.yaw+=(e.key==="ArrowRight"?.12:e.key==="ArrowLeft"?-.12:0);S.pitch+=(e.key==="ArrowDown"?.12:e.key==="ArrowUp"?-.12:0);drawBrain();}});
$("download-cell").addEventListener("click",()=>download(f().recording_id+"-cell-"+c().id+".json",JSON.stringify({recording:f().recording_id,source:"10.25378/janelia.31967640.v1",cell:c(),protocol:f().protocol,source_coordinate_units:"native voxel indices; not physical or registered",target_position:f().target_position,reference_position:f().reference_position,trial_support:f().counts,source_manifest:"repro/"+f().recording_id+"_manifest.json",alternative_reference_protocol:f().alternative_reference??null,mapping_status:"unverified position-to-organ mapping; no causal or nutrient-specific identity",next_experiment:$("plan").textContent},null,2),"application/json"));
$("download-shortlist").addEventListener("click",()=>{
 const rows=[["recording","source_cell_id","discovery_rank","x_voxel","y_voxel","z_voxel","selected_position_code","reference_code","discovery_gain_pct","reserved_site_gain_pct","reserved_multimodal_gain_pct","reserved_full_r","no_drift_gain_pct","alternate_reference_code","alternate_site_gain_pct","alternate_full_r","reference_sensitive","mapping_status"]];
 S.records.forEach(r=>shortlisted(r).sort((a,b)=>a.discovery_rank-b.discovery_rank).forEach(v=>rows.push([r.recording_id,v.id,v.discovery_rank,...v.xyz,r.target_position,r.reference_position,v.discovery_score,v.metrics.nutrient_error_reduction_pct,v.metrics.multimodal_error_reduction_pct,v.metrics.full_r,v.original_no_drift.nutrient_error_reduction_pct,v.alternative_reference?.reference_position??"",v.alternative_reference?.metrics.nutrient_error_reduction_pct??"",v.alternative_reference?.metrics.full_r??"",v.alternative_reference?(v.alternative_reference.metrics.nutrient_error_reduction_pct<=0||v.alternative_reference.metrics.full_r<=0):"not tested","unverified organ mapping"])));
 download("discovery-shortlist.csv",rows.map(row=>row.map(v=>'"'+String(v??"").replaceAll('"','""')+'"').join(",")).join("\n"),"text/csv");});
(async()=>{try{S.records=await Promise.all(window.DEMO_RECORDS.map(async url=>{const r=await fetch(url+"?v=2");if(!r.ok)throw new Error("Analysis file unavailable");return r.json();}));options($("fish"),S.records.map((r,i)=>[i,r.label]));$("explorer").hidden=false;text("load-status",S.records.reduce((n,r)=>n+r.cells.length,0).toLocaleString()+" source-index sampled, QC-passed ROIs across "+S.records.length+" recordings. Organ mapping remains unverified; all position labels are literal acquisition codes.");document.querySelectorAll("select,button").forEach(e=>e.disabled=false);setupRecording();}catch(e){text("load-status","The analysis could not load: "+e.message+". Reload or use the linked repository; no example results are substituted.");$("load-status").classList.add("error");}})();
})();

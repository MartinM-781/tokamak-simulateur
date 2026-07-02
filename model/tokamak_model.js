/* Modèle physique tokamak v5 (0D, unités arbitraires, t en ms) — extrait tel
   quel de web/tokamak_v5.html. Le code entre MODEL-BEGIN et MODEL-END est la
   référence unique de la physique : pas de DOM, pas d'horloge, pas d'aléa
   caché (toute la stochasticité passe par mulberry32/makeGauss).
   Double export : window.TokamakModel en navigateur ET module.exports en Node
   — pas de module ES, les pages doivent rester ouvrables en file://. */
(function (root, factory) {
  if (typeof module === 'object' && typeof module.exports === 'object') {
    module.exports = factory();
  } else {
    root.TokamakModel = factory();
  }
})(typeof self !== 'undefined' ? self : this, function () {
/*MODEL-BEGIN*/
var MP={CR:0.0008,WSAT:0.35,W32SAT:0.25,WSEED:0.004,W32SEED:0.003,
  CPL:10,D32:-0.3,LBOOST:1.5,TAUV:20,LOCKFRAC:0.15,TAULOCK:1.5,
  TAUE:40,DEG:1.2,SAWP:24,SAWDROP:0.92,SAWTE:0.55,
  TQTAU:0.8,TQTE:0.03,CQTE:0.15,SPIKE0:0.08,SPTAU:1.2,CQBASE:6,CQK:50,ENDIP:0.02,
  DRS:0.306,MIRA:0.4,TMAX:2500};
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;var t=Math.imul(a^a>>>15,1|a);
  t=t+Math.imul(t^t>>>7,61|t)^t;return ((t^t>>>14)>>>0)/4294967296;};}
function makeGauss(rng){var spare=null;return function(){
  if(spare!==null){var s=spare;spare=null;return s;}
  var u=0,v=0;while(u===0)u=rng();while(v===0)v=rng();
  var m=Math.sqrt(-2*Math.log(u)),th=2*Math.PI*v;
  spare=m*Math.sin(th);return m*Math.cos(th);};}
function newState(P){
  return {t:0,W21:MP.WSEED,W32:MP.W32SEED,Om:2*Math.PI*P.f0,Om0:2*Math.PI*P.f0,
    locked:false,tq:false,cq:false,ended:false,
    Te:1,Ip:1,spike:0,sawT:0,rotPh:0,K:0,phaseId:0,
    nM:0,nT:0,nI:0,
    tLock:-1,tTQ:-1,tCQ:-1,tEnd:-1};
}
function stepModel(S,P,g,dt){
  if(S.ended){
    S.t+=dt;
    S.Te+=(0-S.Te)*dt/20;S.Ip+=(0-S.Ip)*dt/20;
    var sg0=(P.noise/100)*0.2431;
    S.nM=0.97*S.nM+sg0*0.8*g();S.nT=0.97*S.nT+sg0*g();S.nI=0.97*S.nI+sg0*0.6*g();
    return {mir:S.nM,te:S.Te+S.nT,ip:S.Ip+S.nI};
  }
  var eta=Math.pow(Math.max(S.Te,0.02),-1.5);
  var dEff=P.d0+(S.locked?MP.LBOOST:0);
  S.W21+=MP.CR*eta*dEff*(1-S.W21/MP.WSAT)*dt;
  if(S.W21<MP.WSEED)S.W21=MP.WSEED;
  if(S.W21>MP.WSAT)S.W21=MP.WSAT;
  var d32=MP.D32+MP.CPL*S.W21;
  S.W32+=MP.CR*eta*d32*(1-S.W32/MP.W32SAT)*dt;
  if(S.W32<MP.W32SEED)S.W32=MP.W32SEED;
  if(S.W32>MP.W32SAT)S.W32=MP.W32SAT;
  if(!S.locked){
    S.Om+=((S.Om0-S.Om)/MP.TAUV-P.cb*S.W21*S.W21*S.Om)*dt;
    if(S.Om<MP.LOCKFRAC*S.Om0){S.locked=true;S.tLock=S.t;}
  }else{S.Om+=(0-S.Om)*dt/MP.TAULOCK;}
  S.K=(S.W21+S.W32)/MP.DRS;
  if(!S.tq&&S.K>=1){S.tq=true;S.tTQ=S.t;S.spike=MP.SPIKE0;}
  if(S.tq){S.Te+=(MP.TQTE-S.Te)*dt/MP.TQTAU;}
  else{
    var teq=Math.max(0.1,1-MP.DEG*(S.W21+S.W32));
    S.Te+=(teq-S.Te)*dt/MP.TAUE;
    S.sawT+=dt;
    if(S.sawT>MP.SAWP&&S.Te>MP.SAWTE){S.sawT=0;S.Te*=MP.SAWDROP;}
  }
  if(S.tq&&!S.cq&&S.Te<MP.CQTE){S.cq=true;S.tCQ=S.t;}
  S.spike*=Math.exp(-dt/MP.SPTAU);
  if(S.cq){
    var tcq=MP.CQBASE+MP.CQK*Math.pow(Math.max(S.Te,0.02),1.5);
    S.Ip+=(-S.Ip/tcq)*dt;
    if(S.Ip<MP.ENDIP){S.Ip=MP.ENDIP;S.ended=true;S.tEnd=S.t;}
  }
  S.rotPh+=S.Om*dt;
  var mir=0;
  if(!S.ended){
    mir=MP.MIRA*Math.pow(S.W21/0.05,2)*(S.Om/S.Om0)*Math.sin(S.rotPh)
       +0.35*MP.MIRA*Math.pow(S.W32/0.05,2)*(S.Om/S.Om0)*Math.sin(1.5*S.rotPh+1.1);
    if(mir>3)mir=3;if(mir<-3)mir=-3;
  }
  var sg=(P.noise/100)*0.2431;
  S.nM=0.97*S.nM+sg*0.8*g();
  S.nT=0.97*S.nT+sg*g();
  S.nI=0.97*S.nI+sg*0.6*g();
  if(S.t<S.tLock||S.tLock<0){if(S.W21>=0.03)S.phaseId=1;else S.phaseId=0;}
  if(S.locked)S.phaseId=2;
  if(S.tq)S.phaseId=3;
  if(S.cq)S.phaseId=4;
  if(S.ended)S.phaseId=5;
  S.t+=dt;
  return {mir:mir+S.nM,te:S.Te+S.nT,ip:S.Ip+S.spike+S.nI};
}
/*MODEL-END*/
  return {
    MP: MP,
    mulberry32: mulberry32,
    makeGauss: makeGauss,
    newState: newState,
    stepModel: stepModel
  };
});

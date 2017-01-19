function Ln(x) { return Math.log(x) } function Exp(x) { return Math.exp(x)}
function xlx(x) { return x*Ln(x+1e-20) }
function Abs(x) { return Math.abs(x) }
function Sqrt(x) { return Math.sqrt(x) }
function Cos(x) { return Math.cos(x) }
function Pow(x,n) { return Math.pow(x,n) }
function fEnt( x ) { return x * Ln(x) / Ln(2) }
// Global variables
var Pi=3.141592653589793; Pi2=2*Pi; LnPi2 = Ln(Pi2); PiD2=Pi/2
var Cell_A; var Cell_B; var Cell_C; var Cell_D; var N
var Cell_r1; var Cell_r2; var Cell_c1; var Cell_c2; var t
var Ex_A; var Ex_B; var Ex_C; var Ex_D; var Sav_A; var Sav_B; var Sav_C; var Sav_D
var cs; var od; var rr; var kp; var fc; var mcr; var sn; var sp; var pp; var np
var arr; var rrr; var plr; var nlr; var dor; var yj; var nnd
var dp; var nn; var nmi; var cc; var ca; var cp; var yq; var ets
function CalcTots(form) {
    Cell_A = eval(form.Cell_A.value)
    Cell_B = eval(form.Cell_B.value)
    Cell_C = eval(form.Cell_C.value)
    Cell_D = eval(form.Cell_D.value)
    Cell_r1 = Cell_A+Cell_B
    Cell_r2 = Cell_C+Cell_D
    Cell_c1 = Cell_A+Cell_C
    Cell_c2 = Cell_B+Cell_D
    // Update row and column totals
    if(typeof Cell_A !== "undefined" & typeof Cell_B !== "undefined")
        document.getElementById("r1").innerText = Cell_r1 + " = r1"
    if(typeof Cell_C !== "undefined" & typeof Cell_D !== "undefined")
        document.getElementById("r2").innerText = Cell_r2 + " = r2"
    if(typeof Cell_A !== "undefined" & typeof Cell_C !== "undefined")
        document.getElementById("c1").innerText = Cell_c1 + " = c1"
    if(typeof Cell_B !== "undefined" & typeof Cell_D !== "undefined")
        document.getElementById("c2").innerText = Cell_c2 + " = c2"
    if((typeof Cell_A !== "undefined") & (typeof Cell_B !== "undefined") & (typeof Cell_C !== "undefined") & (typeof Cell_D !== "undefined")) {
        t = Cell_A + Cell_B + Cell_C + Cell_D
        document.getElementById("t" ).innerText = t + " = t"
        }
}
function CalcStats(form) {
    var LoSlop = Cell_A; if(Cell_D<Cell_A) { LoSlop = Cell_D }
    var HiSlop = Cell_B; if(Cell_C<Cell_B) { HiSlop = Cell_C }
    var LnProb1 = LnFact(Cell_r1) + LnFact(Cell_r2) + LnFact(Cell_c1) + LnFact(Cell_c2) - LnFact(t)
    var SingleP = Exp( LnProb1 - LnFact(Cell_A) - LnFact(Cell_B) - LnFact(Cell_C) - LnFact(Cell_D) )
    var FisherP=0; var LeftP=0; var RightP=0; var RosnerP=0; var SumCheck=0
    var k = Cell_A - LoSlop
    while( k<=Cell_A+HiSlop ) {
        var P = Exp( LnProb1 - LnFact(k) - LnFact(Cell_r1-k) - LnFact(Cell_c1-k) - LnFact(k+Cell_r2-Cell_c1) )
        SumCheck = SumCheck + P
        if( k<=Cell_A ) { LeftP = LeftP + P }
        if( k>=Cell_A ) { RightP = RightP + P }
        if( P<=(SingleP+1e-12) ) { FisherP = FisherP + P }
        k = k + 1
    }
    form.LeftP.value = Fmt(LeftP)
    form.RightP.value = Fmt(RightP)
    form.FisherP.value = Fmt(FisherP)
    form.SingleP.value = Fmt(SingleP)
    form.SumCheck.value = "" + SumCheck
    RosnerP = 0.5
    if( LeftP<RosnerP ) { RosnerP = LeftP }
    if( RightP<RosnerP ) { RosnerP = RightP }
    RosnerP = 2*RosnerP
    form.RosnerP.value = Fmt(RosnerP)
    Ex_A = Cell_r1*Cell_c1/t; Sav_A = Ex_A
    Ex_B = Cell_r1*Cell_c2/t; Sav_B = Ex_B
    Ex_C = Cell_r2*Cell_c1/t; Sav_C = Ex_C
    Ex_D = Cell_r2*Cell_c2/t; Sav_D = Ex_D
    // Insert Expected frequencies
    document.getElementById("ex_a").innerText  = Fmt(Ex_A);
    document.getElementById("ex_b").innerText  = Fmt(Ex_B);
    document.getElementById("ex_c").innerText  = Fmt(Ex_C);
    document.getElementById("ex_d").innerText  = Fmt(Ex_D);
    document.getElementById("ex_ab").innerText = Fmt(Ex_A+Ex_B);
    document.getElementById("ex_cd").innerText = Fmt(Ex_C+Ex_D);
    document.getElementById("ex_ac").innerText = Fmt(Ex_A+Ex_C);
    document.getElementById("ex_bd").innerText = Fmt(Ex_B+Ex_D);
    // Insert Proportions
    document.getElementById("pr_a").innerText  = Fmt(Cell_A/t)
    document.getElementById("pr_b").innerText  = Fmt(Cell_B/t)
    document.getElementById("pr_c").innerText  = Fmt(Cell_C/t)
    document.getElementById("pr_d").innerText  = Fmt(Cell_D/t)
    document.getElementById("pr_r1").innerText = Fmt((Cell_A+Cell_B)/t)
    document.getElementById("pr_r2").innerText = Fmt((Cell_C+Cell_D)/t)
    document.getElementById("pr_c1").innerText = Fmt((Cell_A+Cell_C)/t)
    document.getElementById("pr_c2").innerText = Fmt((Cell_B+Cell_D)/t)
    // Insert Col Proportions
    document.getElementById("col_pr_a").innerText  = Fmt(Cell_A/(Cell_A+Cell_C))
    document.getElementById("col_pr_b").innerText  = Fmt(Cell_B/(Cell_B+Cell_D))
    document.getElementById("col_pr_c").innerText  = Fmt(Cell_C/(Cell_A+Cell_C))
    document.getElementById("col_pr_d").innerText  = Fmt(Cell_D/(Cell_B+Cell_D))
    document.getElementById("col_pr_c1").innerText = Fmt((Cell_A+Cell_C)/t)
    document.getElementById("col_pr_c2").innerText = Fmt((Cell_B+Cell_D)/t)
    // Insert Row Proportions
    document.getElementById("row_pr_a").innerText  = Fmt(Cell_A/(Cell_A+Cell_B))
    document.getElementById("row_pr_b").innerText  = Fmt(Cell_B/(Cell_A+Cell_B))
    document.getElementById("row_pr_c").innerText  = Fmt(Cell_C/(Cell_C+Cell_D))
    document.getElementById("row_pr_d").innerText  = Fmt(Cell_D/(Cell_C+Cell_D))
    document.getElementById("row_pr_r1").innerText = Fmt((Cell_A+Cell_B)/t)
    document.getElementById("row_pr_r2").innerText = Fmt((Cell_C+Cell_D)/t)
    cs=csq(Cell_A,Ex_A,.5)+csq(Cell_B,Ex_B,.5)+csq(Cell_C,Ex_C,.5)+csq(Cell_D,Ex_D,.5)
    form.csyc.value = Fmt(cs)
    form.csyc_p.value = Fmt(Csp(cs))
    od=(Cell_A/Cell_B)/(Cell_C/Cell_D); form.od.value=Fmt(od)
    rr=(Cell_A/Cell_r1)/(Cell_C/Cell_r2); form.rr.value=Fmt(rr)
    kp=2*(Cell_A*Cell_D-Cell_B*Cell_C)/((Cell_B+Cell_A)*(Cell_B+Cell_D)+(Cell_A+Cell_C)*(Cell_D+Cell_C)); form.kp.value=Fmt(kp)
    fc=(Cell_A+Cell_D)/t; form.fc.value=Fmt(fc); form.mcr.value=Fmt(1-fc)
    sn=Cell_A/Cell_c1; form.sn.value=Fmt(sn)
    sp=Cell_D/Cell_c2; form.sp.value=Fmt(sp)
    pv=(Cell_A+Cell_C)/t; form.pv.value=Fmt(pv)     // observed (sample) prevalence
    pp=Cell_A/Cell_r1; form.pp.value=Fmt(pp)
    pf=parseFloat(form.prev.value); // user given prevalence
    ppc=(sn*pf)/( (sn*pf)+(1-sp)*(1-pf) ); form.ppc.value=Fmt(ppc)  // PP corrected for user specified prevalence
    np=Cell_D/Cell_r2; form.np.value=Fmt(np)
    npc=(sp*(1-pf))/((1-sn)*pf+sp*(1-pf)); form.npc.value=Fmt(npc)  // NP corrected for user specified prevalence
    plr=sn/(1-sp); form.plr.value=Fmt(plr)
    nlr=(1-sn)/sp; form.nlr.value=Fmt(nlr)
    dor=(sn/(1-sn))/((1-sp)/sp); form.dor.value=Fmt(dor)
    eor=(sn/(1-sn))/(sp/(1-sp)); form.eor.value=Fmt(eor)
    yj=sn+sp-1; form.yj.value=Fmt(yj)
    nnd=1/yj; form.nnd.value=Fmt(nnd)
    nnm=1/(1-fc); form.nnm.value=Fmt(nnm)
    dp=Cell_A/Cell_r1-Cell_C/Cell_r2; form.dp.value=Fmt(dp)
    arr=-dp; form.arr.value=Fmt(arr)
    rrr=arr/(Cell_C/Cell_r2); form.rrr.value=Fmt(rrr)
    nmi = 1 - ( xlx(Cell_B+Cell_A) + xlx(Cell_D+Cell_C) - xlx(Cell_B) - xlx(Cell_A) - xlx(Cell_D) - xlx(Cell_C) ) / ( xlx(t) - xlx(Cell_c2) - xlx(Cell_c1) ); form.nmi.value=Fmt(nmi)
    csny=csq(Cell_B,Ex_B,0)+csq(Cell_A,Ex_A,0)+csq(Cell_D,Ex_D,0)+csq(Cell_C,Ex_C,0)
    form.csny.value = Fmt(csny)
    form.csny_p.value = Fmt(Csp(csny))
    csmh=(t-1)*(Cell_A*Cell_D-Cell_B*Cell_C)*(Cell_A*Cell_D-Cell_B*Cell_C)/(Cell_r1*Cell_r2*Cell_c1*Cell_c2)
    form.csmh.value = Fmt(csmh)
    form.csmh_p.value = Fmt(Csp(csmh))
    cc=Sqrt(csny/(csny+t)); form.cc.value=Fmt(cc)
    ca=cc*Sqrt(2); form.ca.value=Fmt(ca)
    rtet=Cos(Pi/(1+Sqrt(od))); form.rtet.value = Fmt(rtet)
    cp=(Cell_A*Cell_D-Cell_B*Cell_C)/Sqrt(Cell_r1*Cell_r2*Cell_c2*Cell_c1); form.cp.value=Fmt(cp)
    yq=(Cell_A*Cell_D-Cell_B*Cell_C)/(Cell_B*Cell_C+Cell_A*Cell_D); form.yq.value=Fmt(yq)
    ets = (Cell_A - Sav_A) / (Cell_A + Cell_B + Cell_C - Sav_A); form.ets.value=Fmt(ets)
    EntR = - ( fEnt(Cell_r1/t) + fEnt(Cell_r2/t) ); form.EntR.value=Fmt(EntR)
    EntC = - ( fEnt(Cell_c1/t) + fEnt(Cell_c2/t) ); form.EntC.value=Fmt(EntC)
    EntRC = - ( fEnt(Cell_A/t) + fEnt(Cell_B/t) + fEnt(Cell_C/t) + fEnt(Cell_D/t) ); form.EntRC.value=Fmt(EntRC)
    EntIB = EntR + EntC - EntRC; form.EntIB.value=Fmt(EntIB)
    EntIA = EntRC - EntR; form.EntIA.value=Fmt(EntIA)
    EntIC = EntRC - EntC; form.EntIC.value=Fmt(EntIC)
    EntSim = EntIB / ( EntIA + EntIB + EntIC ); form.EntSim.value=Fmt(EntSim)
    EntDif = 1 - EntSim; form.EntDif.value=Fmt(EntDif)
    var pcrit = (100-form.ConfLevel.value)/100
    var del=LoSlop
    Ex_B=Cell_B+LoSlop; Ex_A=Cell_A-LoSlop; Ex_D=Cell_D-LoSlop; Ex_C=Cell_C+LoSlop; var pval=0
    while(del>0.000001) {
        del=del/2
        if(pval<pcrit) {
            Ex_B=Ex_B-del
        } else {
            Ex_B=Ex_B+del
        }
        Ex_A=Cell_r1-Ex_B; Ex_D=Cell_c2-Ex_B; Ex_C=Cell_r2-Ex_D
        pval=Csp(csq(Cell_B,Ex_B,0.5)+csq(Cell_A,Ex_A,0.5)+csq(Cell_D,Ex_D,0.5)+csq(Cell_C,Ex_C,0.5))
    }
    form.Low_A.value=Fmt(Ex_A); form.Low_B.value=Fmt(Ex_B); form.Low_C.value=Fmt(Ex_C); form.Low_D.value=Fmt(Ex_D);
    od=(Ex_A/Ex_B)/(Ex_C/Ex_D); form.od_lo.value=Fmt(od)
    rr=(Ex_A/Cell_r1)/(Ex_C/Cell_r2); form.rr_lo.value=Fmt(rr)
    kp=2*(Ex_A*Ex_D-Ex_B*Ex_C)/((Ex_B+Ex_A)*(Ex_B+Ex_D)+(Ex_A+Ex_C)*(Ex_D+Ex_C)); form.kp_lo.value=Fmt(kp)
    fc=(Ex_A+Ex_D)/(Ex_B+Ex_A+Ex_D+Ex_C); form.fc_lo.value=Fmt(fc); form.mcr_hi.value=Fmt(1-fc)
    sn=Ex_A/Cell_c1; form.sn_lo.value=Fmt(sn)
    sp=Ex_D/Cell_c2; form.sp_lo.value=Fmt(sp)
    plr=sn/(1-sp); form.plr_lo.value=Fmt(plr)
    nlr=(1-sn)/sp; form.nlr_hi.value=Fmt(nlr)
    dor=(sn/(1-sn))/((1-sp)/sp); form.dor_lo.value=Fmt(dor)
    eor=(sn/(1-sn))/(sp/(1-sp)); form.eor_lo.value=Fmt(eor)
    yj=sn+sp-1; form.yj_lo.value=Fmt(yj)
    nnd=1/yj; form.nnd_hi.value=Fmt(nnd)
    nnm=1/(1-fc); form.nnm_lo.value=Fmt(nnm)
    form.pv_lo.value=Fmt(ciw(Cell_c1,t, pcrit,0))
    form.pv_hi.value=Fmt(ciw(Cell_c1,t, pcrit,1))
    pp=Ex_A/Cell_r1; form.pp_lo.value=Fmt(pp)
    form.ppc_lo.value=Fmt(ciw(ppc*Cell_r1,Cell_r1, pcrit,0))
    form.ppc_hi.value=Fmt(ciw(ppc*Cell_r1,Cell_r1, pcrit,1))
    np=Ex_D/Cell_r2; form.np_lo.value=Fmt(np)
    dplo=Ex_A/Cell_r1-Ex_C/Cell_r2; form.dp_lo.value=Fmt(dplo)
    arr=-dplo; form.arr_hi.value=Fmt(arr)
    rrr=arr/(Ex_C/Cell_r2); form.rrr_hi.value=Fmt(rrr)
    nmi = 1 - ( xlx(Ex_B+Ex_A) + xlx(Ex_D+Ex_C) - xlx(Ex_B) - xlx(Ex_A) - xlx(Ex_D) - xlx(Ex_C) ) / ( xlx(t) - xlx(Cell_c2) - xlx(Cell_c1) ); form.nmi_lo.value=Fmt(nmi)
    csny=csq(Ex_B,Sav_B,0)+csq(Ex_A,Sav_A,0)+csq(Ex_D,Sav_D,0)+csq(Ex_C,Sav_C,0)
    cc=Sqrt(csny/(csny+t)); form.cc_lo.value=Fmt(cc)
    ca=cc*Sqrt(2); form.ca_lo.value=Fmt(ca)
    rtet=Cos(Pi/(1+Sqrt(od))); form.rtet_lo.value = Fmt(rtet)
    cp=(Ex_A*Ex_D-Ex_B*Ex_C)/Sqrt(Cell_r1*Cell_r2*Cell_c2*Cell_c1); form.cp_lo.value=Fmt(cp)
    yq=(Ex_A*Ex_D-Ex_B*Ex_C)/(Ex_B*Ex_C+Ex_A*Ex_D); form.yq_lo.value=Fmt(yq)
    ets = (Ex_A - Sav_A) / (Ex_A + Ex_B + Ex_C - Sav_A); form.ets_lo.value=Fmt(ets)
    EntR = - ( fEnt(Cell_r1/t) + fEnt(Cell_r2/t) ); form.EntR_lo.value=Fmt(EntR)
    EntC = - ( fEnt(Cell_c1/t) + fEnt(Cell_c2/t) ); form.EntC_lo.value=Fmt(EntC)
    EntRC = - ( fEnt(Ex_A/t) + fEnt(Ex_B/t) + fEnt(Ex_C/t) + fEnt(Ex_D/t) ); form.EntRC_lo.value=Fmt(EntRC)
    EntIB = EntR + EntC - EntRC; form.EntIB_hi.value=Fmt(EntIB)
    EntIA = EntRC - EntR; form.EntIA_lo.value=Fmt(EntIA)
    EntIC = EntRC - EntC; form.EntIC_lo.value=Fmt(EntIC)
    EntSim = EntIB / ( EntIA + EntIB + EntIC ); form.EntSim_hi.value=Fmt(EntSim)
    EntDif = 1 - EntSim; form.EntDif_lo.value=Fmt(EntDif)
    del=HiSlop
    Ex_B=Cell_B-HiSlop; Ex_A=Cell_A+HiSlop; Ex_D=Cell_D+HiSlop; Ex_C=Cell_C-HiSlop; var pval=0
    while(del>0.000001) {
        del=del/2
        if(pval<pcrit) {
            Ex_B=Ex_B+del
        } else {
            Ex_B=Ex_B-del
        }
        Ex_A=Cell_r1-Ex_B; Ex_D=Cell_c2-Ex_B; Ex_C=Cell_r2-Ex_D
        pval=Csp(csq(Cell_B,Ex_B,0.5)+csq(Cell_A,Ex_A,0.5)+csq(Cell_D,Ex_D,0.5)+csq(Cell_C,Ex_C,0.5))
    }
    form.High_A.value=Fmt(Ex_A); form.High_B.value=Fmt(Ex_B); form.High_C.value=Fmt(Ex_C); form.High_D.value=Fmt(Ex_D);
    od=(Ex_A/Ex_B)/(Ex_C/Ex_D); form.od_hi.value=Fmt(od)
    rr=(Ex_A/Cell_r1)/(Ex_C/Cell_r2); form.rr_hi.value=Fmt(rr)
    kp=2*(Ex_A*Ex_D-Ex_B*Ex_C)/((Ex_B+Ex_A)*(Ex_B+Ex_D)+(Ex_A+Ex_C)*(Ex_D+Ex_C)); form.kp_hi.value=Fmt(kp)
    fc=(Ex_A+Ex_D)/(Ex_B+Ex_A+Ex_D+Ex_C); form.fc_hi.value=Fmt(fc); form.mcr_lo.value=Fmt(1-fc)
    sn=Ex_A/Cell_c1; form.sn_hi.value=Fmt(sn)
    sp=Ex_D/Cell_c2; form.sp_hi.value=Fmt(sp)
    plr=sn/(1-sp); form.plr_hi.value=Fmt(plr)
    nlr=(1-sn)/sp; form.nlr_lo.value=Fmt(nlr)
    dor=(sn/(1-sn))/((1-sp)/sp); form.dor_hi.value=Fmt(dor)
    eor=(sn/(1-sn))/(sp/(1-sp)); form.eor_hi.value=Fmt(eor)
    yj=sn+sp-1; form.yj_hi.value=Fmt(yj)
    nnd=1/yj; form.nnd_lo.value=Fmt(nnd)
    nnm=1/(1-fc); form.nnm_hi.value=Fmt(nnm)
    pp=Ex_A/Cell_r1; form.pp_hi.value=Fmt(pp)
    np=Ex_D/Cell_r2; form.np_hi.value=Fmt(np)
    dphi=Ex_A/Cell_r1-Ex_C/Cell_r2; form.dp_hi.value=Fmt(dphi)
    arr=-dphi; form.arr_lo.value=Fmt(arr)
    rrr=arr/(Ex_C/Cell_r2); form.rrr_lo.value=Fmt(rrr)
    nmi = 1 - ( xlx(Ex_B+Ex_A) + xlx(Ex_D+Ex_C) - xlx(Ex_B) - xlx(Ex_A) - xlx(Ex_D) - xlx(Ex_C) ) / ( xlx(t) - xlx(Cell_c2) - xlx(Cell_c1) ); form.nmi_hi.value=Fmt(nmi)
    csny=csq(Ex_B,Sav_B,0)+csq(Ex_A,Sav_A,0)+csq(Ex_D,Sav_D,0)+csq(Ex_C,Sav_C,0)
    cc=Sqrt(csny/(csny+t)); form.cc_hi.value=Fmt(cc)
    ca=cc*Sqrt(2); form.ca_hi.value=Fmt(ca)
    rtet=Cos(Pi/(1+Sqrt(od))); form.rtet_hi.value = Fmt(rtet)
    cp=(Ex_A*Ex_D-Ex_B*Ex_C)/Sqrt(Cell_r1*Cell_r2*Cell_c2*Cell_c1); form.cp_hi.value=Fmt(cp)
    yq=(Ex_A*Ex_D-Ex_B*Ex_C)/(Ex_B*Ex_C+Ex_A*Ex_D); form.yq_hi.value=Fmt(yq)
    ets = (Ex_A - Sav_A) / (Ex_A + Ex_B + Ex_C - Sav_A); form.ets_hi.value=Fmt(ets)
    EntR = - ( fEnt(Cell_r1/t) + fEnt(Cell_r2/t) ); form.EntR_hi.value=Fmt(EntR)
    EntC = - ( fEnt(Cell_c1/t) + fEnt(Cell_c2/t) ); form.EntC_hi.value=Fmt(EntC)
    EntRC = - ( fEnt(Ex_A/t) + fEnt(Ex_B/t) + fEnt(Ex_C/t) + fEnt(Ex_D/t) ); form.EntRC_hi.value=Fmt(EntRC)
    EntIB = EntR + EntC - EntRC; form.EntIB_lo.value=Fmt(EntIB)
    EntIA = EntRC - EntR; form.EntIA_hi.value=Fmt(EntIA)
    EntIC = EntRC - EntC; form.EntIC_hi.value=Fmt(EntIC)
    EntSim = EntIB / ( EntIA + EntIB + EntIC ); form.EntSim_lo.value=Fmt(EntSim)
    EntDif = 1 - EntSim; form.EntDif_hi.value=Fmt(EntDif)
    if(dp==0) { form.nn.value="Infinite" } else { form.nn.value=Fmt(Math.abs(1/dp)) }
    form.nn_lo.value="Unknown"; form.nn_hi.value="Unknown"
    if(dplo<0 & dphi<0) { form.nn_lo.value=Fmt(-1/dplo); form.nn_hi.value=Fmt(-1/dphi) }
    if(dplo<0 & dphi==0) { form.nn_Lo.value=Fmt(-1/dplo); form.nn_hi.value="Infinite" }
    if(dplo<0 & dphi>0) { form.nn_lo.value=Fmt(1/Math.max(Math.abs(dplo),Math.abs(dphi))); form.nn_hi.value="Infinite" }
    if(dplo==0 & dphi>0) { form.nn_lo.value=Fmt(1/dphi); form.nn_hi.value="Infinite" }
    if(dplo>0 & dphi>0) { form.nn_lo.value=Fmt(1/dphi); form.nn_hi.value=Fmt(1/dplo) }
    if(Math.abs(form.prev.value-pv)<.001){  // if prevalences is alike, use same CI for predictive and adjusted
        form.ppc_lo.value=form.pp_lo.value
        form.ppc_hi.value=form.pp_hi.value
        form.npc_lo.value=form.np_lo.value
        form.npc_hi.value=form.np_hi.value
    } else {
        form.ppc_lo.value=Fmt(ciw(ppc*Cell_r1,Cell_r1, pcrit,0))
        form.ppc_hi.value=Fmt(ciw(ppc*Cell_r1,Cell_r1, pcrit,1))
        form.npc_lo.value=Fmt(ciw(npc*Cell_r2,Cell_r2, pcrit,0))
        form.npc_hi.value=Fmt(ciw(npc*Cell_r2,Cell_r2, pcrit,1))
    }
    var RIOC = rioc(Cell_A, Cell_D, Cell_c1, Cell_r1, t, pcrit)
    form.RIOC.value=Fmt(RIOC[0])
    form.RIOC_lo.value=Fmt(RIOC[1])
    form.RIOC_hi.value=Fmt(RIOC[2])
}

function saveCSV(form){
    var data = [
            ["Observed Contingency Table"],
                ["", "Outcome Occurred","Outcome did not Occur","Totals"],
                ["Risk Factor Present or Dx Test Positive", Cell_A, Cell_B, Cell_A + Cell_B],
                ["Risk Factor Absentor Dx Test Negative",   Cell_C, Cell_D, Cell_C + Cell_D],
                ["Totals",                                  Cell_A + Cell_C, Cell_B + Cell_D, t],
            [],
            ["Confidence Level", form.ConfLevel.value + "%"],
            [],
            ["Expected Frequencies"],
                [Sav_A, Sav_B],
                [Sav_C, Sav_D],
            [],
            ["Cell Proportions"],
                [Cell_A/t, Cell_B/t, (Cell_A+Cell_B)/t],
                [Cell_C/t, Cell_D/t, (Cell_C+Cell_D)/t],
                [(Cell_A+Cell_C)/t, (Cell_B+Cell_D)/t],
            [],
            ["Row Proportions"],
                [Cell_A/(Cell_A+Cell_B), Cell_B/(Cell_A+Cell_B), (Cell_A+Cell_B)/t],
                [Cell_C/(Cell_C+Cell_D), Cell_D/(Cell_C+Cell_D), (Cell_C+Cell_D)/t],
            [],
            ["Column Proportions"],
                [Cell_A/(Cell_A+Cell_C), Cell_B/(Cell_B+Cell_D)],
                [Cell_C/(Cell_A+Cell_C), Cell_D/(Cell_B+Cell_D)],
                [(Cell_A+Cell_C)/t, (Cell_B+Cell_D)/t],
            [],
            ["Chi-Square Tests"],
                ["Type of Test", "Chi Square", "d.f.", "p-value"],
                ["Pearson Uncorrected", form.csny.value, form.dfny.value, form.csny_p.value],
                ["Yates Corrected",     form.csyc.value, form.dfyc.value, form.csyc_p.value],
                ["Mantel-Haenszel",     form.csmh.value, form.dfmh.value, form.csmh_p.value],
            [],
            ["Fisher Exact Test"],
                ["Type of comparison (Alternate Hypothesis)", "p-value"],
                ["Two-tailed (to test if the Odds Ratio is significantly different from 1). If you don't know which Fisher Exact p-value to use then use this one", form.FisherP.value],
                ["Left-tailed (to test if the Odds Ratio is significantly less than 1)", form.LeftP.value],
                ["Right-tailed (to test if the Odds Ratio is significantly greater than 1)", form.RightP.value],
                ["Two-tailed p-value calculated as described in Rosner's book: (2 times whichever is smallest: left-tail or right-tail or 0.5) It tends to agree closely with Yates Chi-Square p-value", form.RosnerP.value],
                ["Probability of getting exactly the observed table. (This is not really a p-value; don't use this as a significance test)", form.SingleP.value],
                ["Verification of computational accuracy (This number should be very close to 1.0. The closer the better)", form.SumCheck.value],
            [],
            ["Quantities derived from a 2-by-2 table"],
                ["Statistics", "Value", "Low 95% CI", "High 95% CI"],
                ["Odds Ratio (OR) = (a/b)/(c/d)",                                                                       form.od.value,      form.od_lo.value,     form.od_hi.value],
                ["Relative Risk (RR) = (a/r1)/(c/r2)",                                                                  form.rr.value,      form.rr_lo.value,     form.rr_hi.value],
                ["Kappa",                                                                                               form.kp.value,      form.kp_lo.value,     form.kp_hi.value],
                ["Overall Fraction Correct = (a+d)/t ; (often referred to simply as 'Accuracy')",                       form.fc.value,      form.fc_lo.value,     form.fc_hi.value],
                ["Mis-classification Rate = 1 - Overall Fraction Correct",                                              form.mcr.value,     form.mcr_lo.value,    form.mcr_hi.value],
                ["Sensitivity = a/c1; (use exact Binomial confidence intervals instead of these)",                      form.sn.value,      form.sn_lo.value,     form.sn_hi.value],
                ["Specificity = d/c2; (use exact Binomial confidence intervals instead of these)",                      form.sp.value,      form.sp_lo.value,     form.sp_hi.value],
                ["Prevalence (estimated from sample)",                                                                  form.pv.value,      form.pv_lo.value,     form.pv_hi.value],
                ["Positive Predictive Value (PPV) = a/r1; (use exact Binomial confidence intervals instead of these)",  form.pp.value,      form.pp_lo.value,     form.pp_hi.value],
                ["Adjusted PPV (prevalence: " + form.prev.value*100 + "%)",                                             form.ppc.value,     form.ppc_lo.value,    form.ppc_hi.value],
                ["Negative Predictive Value (NPV) = d/r2; (use exact Binomial confidence intervals instead of these)",  form.np.value,      form.np_lo.value,     form.np_hi.value],
                ["Adjusted NPV (prevalence: " + form.prev.value*100 + "%)",                                             form.npc.value,     form.npc_lo.value,    form.npc_hi.value],
                ["Difference in Proportions (DP) = a/r1 - c/r2",                                                        form.dp.value,      form.dp_lo.value,     form.dp_hi.value],
                ["Number Needed to Treat (NNT) = 1 / absolute value of DP which = 1 / absolute value of ARR",           form.nn.value,      form.nn_lo.value,     form.nn_hi.value],
                ["Absolute Risk Reduction (ARR) = c/r2 - a/r1 which = - DP",                                            form.arr.value,     form.arr_lo.value,    form.arr_hi.value],
                ["Relative Risk Reduction (RRR) = ARR/(c/r2)",                                                          form.rrr.value,     form.rrr_lo.value,    form.rrr_hi.value],
                ["Positive Likelihood Ratio (+LR) = Sensitivity / (1 - Specificity)",                                   form.plr.value,     form.plr_lo.value,    form.plr_hi.value],
                ["Negative Likelihood Ratio (-LR) = (1 - Sensitivity) / Specificity",                                   form.nlr.value,     form.nlr_lo.value,    form.nlr_hi.value],
                ["Diagnostic Odds Ratio = (Sensitivity/(1-Sensitivity))/((1-Specificity)/Specificity)",                 form.dor.value,     form.dor_lo.value,    form.dor_hi.value],
                ["Error Odds Ratio = (Sensitivity/(1-Sensitivity))/(Specificity/(1-Specificity))",                      form.eor.value,     form.eor_lo.value,    form.eor_hi.value],
                ["Youden's J = Sensitivity + Specificity - 1",                                                          form.yj.value,      form.yj_lo.value,     form.yj_hi.value],
                ["Number Needed to Diagnose (NND) = 1 / (Sensitivity - (1 - Specificity) ) = 1 / (Youden's J)",         form.nnd.value,     form.nnd_lo.value,    form.nnd_hi.value],
                ["Number Needed to Mis-diagnose (NNM) = 1 / ( 1 - Accuracy )",                                          form.nnm.value,     form.nnm_lo.value,    form.nnm_hi.value],
                ["Forbes' NMI Index",                                                                                   form.nmi.value,     form.nmi_lo.value,    form.nmi_hi.value],
                ["Contingency Coefficient",                                                                             form.cc.value,      form.cc_lo.value,     form.cc_hi.value],
                ["Adjusted Contingency Coefficient",                                                                    form.ca.value,      form.ca_lo.value,     form.ca_hi.value],
                ["Tetrachoric (terachoric) Correlation Coefficient = Cos( Pi / (1 + Sqrt( OR ) ) )",                    form.rtet.value,    form.rtet_lo.value,   form.rtet_hi.value],
                ["Phi Coefficient (Cramer's Phi and = Cohen's w Index for 2x2 table)",                                  form.cp.value,      form.cp_lo.value,     form.cp_hi.value],
                ["Yule's Q = (a*d-b*c)/(a*d+b*c) = (OR - 1) / (OR + 1)",                                                form.yq.value,      form.yq_lo.value,     form.yq_hi.value],
                ["Equitable Threat Score = (a-e)/(a+b+c-e) where e = r1*c1/t",                                          form.ets.value,     form.ets_lo.value,    form.ets_hi.value],
                ["Entropy H(r) = - ( (r1/t)log2(r1/t) + (r2/t)log2(r2/t) )",                                            form.EntR.value,    form.EntR_lo.value,   form.EntR_hi.value],
                ["Entropy H(c) = - ( (c1/t)log2(c1/t) + (c2/t)log2(c2/t) )",                                            form.EntC.value,    form.EntC_lo.value,   form.EntC_hi.value],
                ["Entropy H(r;c) = - ( (a/t)log2(a/t) + (b/t)log2(b/t) + (c/t)log2(c/t) + (d/t)log2(d/t) )",            form.EntRC.value,   form.EntRC_lo.value,  form.EntRC_hi.value],
                ["Information shared by descriptors r and c: B = H(r) + H(c) - H(r;c)",                                 form.EntIB.value,   form.EntIB_lo.value,  form.EntIB_hi.value],
                ["A = H(r;c) - H(r)",                                                                                   form.EntIA.value,   form.EntIA_lo.value,  form.EntIA_hi.value],
                ["C = H(r;c) - H(c)",                                                                                   form.EntIC.value,   form.EntIC_lo.value,  form.EntIC_hi.value],
                ["Similarity of descriptors r and c: S(r;c) = B / (A + B + C)",                                         form.EntSim.value,  form.EntSim_lo.value, form.EntSim_hi.value],
                ["Distance between r and c: D(r;c) = (A + C) / (A + B + C)",                                            form.EntDif.value,  form.EntDif_lo.value, form.EntDif_hi.value],
                ["Relative Improvement Over Chance (RIOC)",                                                             form.RIOC.value,    form.RIOC_lo.value,   form.RIOC_hi.value],
            [],
            ["Lower limiting table"],
                [form.Low_A.value, form.Low_B.value],
                [form.Low_C.value, form.Low_D.value],
            [],
            ["Upper limiting table"],
                [form.High_A.value, form.High_B.value],
                [form.High_C.value, form.High_D.value]
            ];
    var csvContent = "data:text/csv;charset=utf-8";
    data.forEach(function(infoArray, index){
       dataString = infoArray.join(",");
       csvContent += index < data.length ? dataString+ "\n" : dataString;
    }); 

    var encodedUri = encodeURI(csvContent);
    var link = document.createElement("a");
    link.setAttribute("href", encodedUri);
    link.setAttribute("download", "ctab2x2.csv");
    document.body.appendChild(link); // Required for FF
    link.click(); // This will download the data file named "my_data.csv".
}
function csq(o,e,y) {
    if(e==0) { return 0 }
    var x=Abs(o-e)-y; if(x<0) { return 0 }
    return x*x/e
}
function Csp(x) {
    return ChiSq(x,1)
}
function ChiSq(x,n) {
    if(x>1000 | n>1000) { var q=Norm((Pow(x/n,1/3)+2/(9*n)-1)/Sqrt(2/(9*n)))/2; if (x>n) {return q}{return 1-q} }
    var p=Exp(-0.5*x); if((n%2)==1) { p=p*Sqrt(2*x/Pi) }
    var k=n; while(k>=2) { p=p*x/k; k=k-2 }
    var t=p; var Cell_B=n; while(t>1e-15*p) { Cell_B=Cell_B+2; t=t*x/Cell_B; p=p+t }
    return 1-p
}
function Norm(z) {
    var q=z*z
    if(Abs(z)>7) {return (1-1/q+3/(q*q))*Exp(-q/2)/(Abs(z)*Sqrt(PiD2))} {return ChiSq(q,1) }
}
function Fmt(x) {
    var v
    if(x>=0) { v=''+(x+0.0005) } else { v=''+(x-0.0005) }
    return v.substring(0,v.indexOf('.')+4)
}
function LnFact(z) {
    if(z<2) { return 0 }
    if(z<17) { f=z; while(z>2) { z=z-1; f=f*z }; return Ln(f) }
    return (z+0.5)*Ln(z) - z + LnPi2/2 + 1/(12*z) - 1/(360*Pow(z,3)) + 1/(1260*Pow(z,5)) - 1/(1680*Pow(z,7))
}
function CalcFromDiagnostics(form) {
    var prev = form.prev.value
    var sens = form.sens.value
    var spec = form.spec.value
    var N    = form.tssz.value
    if(N=="")
        N=1000;
    form.Cell_A.value=N*sens*prev;
    form.Cell_B.value=(N*(1-prev))-(N*(1-prev)*spec)
    form.Cell_C.value=(N*prev)-(N*prev*sens)
    form.Cell_D.value=N*(1-prev)*spec
}
function ciw(r,n,p,lv){
// Conf. Int. for proportion. Ref: Statistics with Confidence DG Altman et al. 2.ed; Single sampleCI, Wilson method; p. 46
// r: observed positive, n: sample size, p: alfa, lv: CI level (lower:0, upper:1)
    var z = normsInv(1-p/2, 0, 1)
    var q = 1-r/n
    var a = 2*r+z*z
    var b = z*Math.sqrt(z*z+4*r*q)
    var c = 2*(n+z*z)
    var lv; var ci;
    if(lv==0)
        ci=(a-b)/c
    else
        ci=(a+b)/c
    return ci
}
function rioc(a, d, c, r, t, p){
// Relative Improvement Over Chance (RIOC) and Phi as Measures of Predictive Efficiency and Strength of Association in 2 × 2 TablesDavid P. Farrington and Rolf Loeber
// Journal of Quantitative Criminology, Vol. 5, No. 3 (September 1989), pp. 201-213
    var x = ( t*(a+d)-(r*c+(t-r)*(t-c)) )
    var y = ( Math.min(t*(c+t-r), t*(r+t-c))-(r*c+(t-r)*(t-c)) )
    var rioc = x/y
    var sd = Math.sqrt((r*(t-c))/(t*c*(t-r)))
    var z = normsInv(1-p/2, 0, 1)
    return [rioc, rioc-z*sd, rioc+z*sd]
}


$(document).ready(function(){

    function foo(){
        // Show prevalence in adjusted predictive values fields
        $("td#ppc").html("<font style='color:red;'>*</font>Adjusted PPV (user set prevalence: " +  $("input[name=prev]").val()*100 + "%)");
        $("td#npc").html("<font style='color:red;'>*</font>Adjusted NPV (user set prevalence: " +  $("input[name=prev]").val()*100 + "%)");
        // Show conf. level in table head
        $("th#ConfLevel").html($("input[name=ConfLevel]").val() + "% CI");
    }

    foo();
    $("input[name=Cell_A]").focus();

    $("input#compute").click(function(){
        $("input#Cell_D").focus().blur();   // Hack: needed to update the changed form
        foo();
        $("input#csv").show();
        $("div.hidden-tables").show();
        CalcStats(this.form);
    });

    $("input#diagnostics").click(function(){
        CalcFromDiagnostics(this.form);
        $("input#Cell_D").focus().blur();   // Hack: needed to update the changed form
        foo();
        $("input#csv").show();
        $("table.hidden-tables").show();
        CalcStats(this.form);
    });

    $("input#csv").click(function(){
        saveCSV(this.form);
        $(this).hide();
    });

    $("td#showHelpDialog").click(function(){
        $("#dialog").dialog({
            width: "400px",
            collision: "flipfit",
            show: true, hide: true,
            resizable: false,
            close: function( event, ui ) { $("#help").blur(); }
        });
    });

    var source = [
        "99",
        "97.5",
        "95",
        "90",
        "75",
        "50",
    ];
    
    // Create a jqxComboBox
    $("#jqxcombobox").jqxComboBox({ source: source, selectedIndex: 2, width: '60px', height: '25px' });

});
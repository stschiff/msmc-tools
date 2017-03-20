(* ::Package:: *)

LoadMSMC[fn_,cross_:False,\[Mu]_:1.25 10^-8,gen_:29]:=Module[{t,timesStart,timesEnd,vals,lbl},
t=Import[fn,"Table","HeaderLines"->1];
timesStart=t[[All,2]]/\[Mu] gen;
timesEnd=t[[All,3]]/\[Mu] gen;
If[cross,
	vals=2t[[All,5]]/(t[[All,4]]+t[[All,6]]);
	Transpose@Dataset[<|"timeStart"->timesStart,"timeEnd"->timesEnd,"ccr"->vals|>],
	vals=(1/t[[All,4;;]])/(2\[Mu]);
	Transpose@Dataset[<|"timeStart"->timesStart,"timeEnd"->timesEnd,
		Sequence@Table[ToString@StringForm["popsizes_``",i]->vals[[All,i]],{i,Length@vals[[1]]}]|>]
]]


LoadMSMCcombined[fnLeft_,fnRight_,fnBoth_,\[Mu]_:1.25 10^-8,gen_:29]:=Module[
	{mLeft,mRight,mBoth,interpLeft,interpRight,row,tLeft,tRight,f},
	mLeft=LoadMSMC[fnLeft,False,\[Mu],gen];
	mRight=LoadMSMC[fnRight,False,\[Mu],gen];
	mBoth=LoadMSMC[fnBoth,False,\[Mu],gen];
	interpLeft=MakeInterp[mLeft];
	interpRight=MakeInterp[mRight];
	f[row_] := Module[{tValues,denom,val},
		tValues=Subdivide[row["timeStart"],row["timeEnd"],10];
		denom = Mean[1/interpLeft[#]&/@tValues]+Mean[1/interpRight[#]&/@tValues];
		val = (2/row[[3]])/denom;
		<|"timeStart"->row["timeStart"],"timeEnd"->row["timeEnd"],"ccr"->val|>
	];
	mBoth[;;-2,f]
]


LoadMSMCloop[fn_,row_:-1, cross_:False,\[Mu]_:1.25 10^-8,gen_:29]:=Module[{d,last,timeBoundaries,lambdaVec,popVec,i},
d=Import[fn,"Table"];
last=d[[row]];
timeBoundaries=ImportString[last[[3]],"Table","FieldSeparators"->","][[1]]/\[Mu] gen;
lambdaVec=ImportString[last[[4]],"Table","FieldSeparators"->","][[1]];
popVec=(1/lambdaVec)/(2\[Mu]);
If[cross,
Dataset@Table[
<|"timeStart"->timeBoundaries[[i]],"timeEnd"->timeBoundaries[[i+1]],"ccr"->2popVec[[3*(i-1)+2]]/(popVec[[3*(i-1)+1]]+popVec[[3*i]])|>,
{i,Length@timeBoundaries-1}
],
Dataset@Table[
<|"timeStart"->timeBoundaries[[i]],"timeEnd"->timeBoundaries[[i+1]],"popsize"->popVec[[i]]|>,
{i,Length@timeBoundaries-1}
]
]
]


MakeInterp[dat_]:=Module[
	{f,minT,maxT,interpolationTable,vLeft,vRight},
	interpolationTable=Normal@dat[All,{(#timeStart+#timeEnd)/2,#[[3]]}&];
	f=Interpolation[interpolationTable[[;;-2]],InterpolationOrder->1];
	minT=interpolationTable[[1,1]];
	maxT=interpolationTable[[-2,1]];
	vLeft = dat[1,3];
	vRight = dat[-1,3];
	Piecewise[{{vLeft,#<=minT},{vRight,#>=maxT}},f[#]]&
]


FindCCpoint[interp_,val_]:=Block[{t},Max[0,t/.FindMinimum[(interp[t]-val)^2,t][[2]]]]


GetEstimates[interp_]:={FindCCpoint[#,0.5],FindCCpoint[#,0.25],FindCCpoint[#,0.75]}&@interp//Quiet

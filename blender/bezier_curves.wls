#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Title:: *)
(*Connection between Bezier curves and polynomial curves*)


(* ::Text:: *)
(*Basis functions for Bezier curve*)


bernsteinBasis[i_,j_]:=Simplify[PiecewiseExpand@BernsteinBasis[i,j,t],Assumptions->0<t<1];


(* ::Text:: *)
(*Transform parameter t interval to \[Xi]\[Element][-1,1]*)


var\[Xi]=t->1/2(\[Xi]+1);


(* ::Text:: *)
(*A third order Bezier curve is based on 4 points*)


points=Table[{x[i],y[i]},{i,4}];
polynomx=Simplify[points[[All,1]].Table[bernsteinBasis[3,i],{i,0,3}]/.var\[Xi]]
polynomy=Simplify[points[[All,2]].Table[bernsteinBasis[3,i],{i,0,3}]/.var\[Xi]]
fPoly={polynomx,polynomy};


(* ::Text:: *)
(*Connection between the 2nd and 3rd Bezier points and the tangents of the curve -> The vector from the first to the second (and third to fourth) Bezier point is 2/3 the length of the tangent*)


tangent=D[fPoly,\[Xi]];
Simplify[tangent/.\[Xi]->-1]
Simplify[tangent/.\[Xi]->1]


(* ::Text:: *)
(*Plot the curve for sample points*)


varPoints=Thread[Flatten[points]->{0,0,1,1,6,-1,5,2}];
pltBezier=Graphics[{BezierCurve[points/.varPoints],Green,Line[points/.varPoints],Red,Point[points/.varPoints]}];
pltPolynomial=ParametricPlot[fPoly/.varPoints,{\[Xi],-1,1},PlotStyle->{Thick,Dashed,Red}];

Show[pltBezier,pltPolynomial]

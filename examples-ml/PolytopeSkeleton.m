(* 
   PolytopeSkeleton.m,
   by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version March 3, 1999.
   -- This Needs View3D package ---
*)

Visible[ineq_Integer, {a_,b_}, viewpoint_]:=
	If[N[a[[ineq]].viewpoint] >= N[b[[ineq]]], True, False];
		
Visible[{a_,b_}, viewpoint_]:=
	Map[Visible[#,{a,b},viewpoint]&, Range[Length[b]]];
		
VisibleSkeleton[vertices_,ecd_, ead_, {a_, b_},viewpoint_]:=
	Block[{grtemp, edgelist, i, visibility, viseges, invisedges, visibles,invisibles},
		edgelist={};
		Do[	edgelist=Join[edgelist,Map[{i, #}&, ead[[i]]]],
			{i,1,Length[ead]}
		];edgelist=Select[edgelist,#[[1]]<#[[2]]&];
		visibility=Flatten[Position[Visible[{a, b},viewpoint],True]];
		visedges=Select[edgelist, 
		(Intersection[visibility,ecd[[#[[1]]]],ecd[[#[[2]]]] ] != {}) &];
		invisedges=Complement[edgelist,visedges];
		visibles=(vertices[[#]]&) /@ visedges;
		invisibles=(vertices[[#]]&) /@ invisedges;
		grtemp=Join[{Thickness[ 0.0085]},
			Line /@ visibles, 
			{Thickness[ 0.0035], Dashing[{0.01,0.012}]}, Line /@ invisibles];
		grtemp
	];
		
VisibleSkeleton2[vertices_,ecd_, ead_, flist_,{a_, b_},viewpoint_]:=
	Block[{grtemp, edgelist, i, edgelist0, visibility, viseges, invisedges, visibles,invisibles},
		edgelist0={};
		Do[	edgelist0=Join[edgelist0,Map[{i, #}&, ead[[i]]]],
			{i,1,Length[ead]}
		];edgelist0=Select[edgelist0,#[[1]]<#[[2]]&];
		edgelist=Select[edgelist0,(Intersection[ecd[[#[[1]]]],ecd[[#[[2]]]],flist] != {})&]; 
		visibility=Flatten[Position[Visible[{a, b},viewpoint],True]];
		visedges=Select[edgelist, 
		(Intersection[visibility,ecd[[#[[1]]]],ecd[[#[[2]]]] ] != {}) &];
		invisedges=Complement[edgelist,visedges];
		visibles=(vertices[[#]]&) /@ visedges;
		invisibles=(vertices[[#]]&) /@ invisedges;
		grtemp=Join[{Thickness[ 0.0085]},
			Line /@ visibles, 
			{Thickness[ 0.0035], Dashing[{0.01,0.012}]}, Line /@ invisibles];
		grtemp
	];

getMmaViewPoint[vp_,gr_]:=
	ViewPointFromUser[vp,
		FullOptions[gr,PlotRange],
		FullOptions[gr,BoxRatios]
	];

getUserViewPoint[vp_,gr_]:=
	ViewPointToUser[vp,
		FullOptions[gr,PlotRange],
		FullOptions[gr,BoxRatios]
	];
	
circle[alpha_, r_:2, h_:2]:=
	{r*Cos[alpha], r*Sin[alpha], h}

(* end of PolytopeSkeleton.m *)
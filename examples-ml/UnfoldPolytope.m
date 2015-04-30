(*

        Unfolding  Convex Polytopes

             Version 1.0 Beta
              February, 1992

           Copyright (c) 1992 by
              Makoto Namiki

This package contains Mathematica implementations
of unfolding 3-dimsional polytope.

This package is copyright 1992 by Makoto Namiki.
This package may be copied in its entirety for 
nonprofit purposes only. Sale, other than for
the direct cost of the media, is prohibited.  
This copyright notice must accompany all copies.

The authors make no representations, express or 
implied, with respond to this documentation, of 
the software it describes and contains, including
without limitations, any implied warranties of 
mechantability or fitness for a particular purpose,
all of which are expressly disclaimed.  The authors
shall in no event be liable for any indirect,
incidental, or consequential damages.

This beta release is designed to run under 
Version 1.2 & 2.0 of Mathematica. Any comments, 
bug reports, or requests to get on the 
UnfoldPolytope mailing list should be 
forwarded to:

  Makoto Namiki
  Department of Information Science,
  Tokyo Institute of Technology, Tokyo
  2-12-1 Oh-okayama, Meguro-ku
  Tokyo 145, Japan

  +81-3-3726-1111 (Ex. 4131)
  namiki@titisna.is.titech.ac.jp
 
*)

BeginPackage["UnfoldPolytope`"]

UnfoldPolytope::usage="UnfoldPolytope[facets] gives graphics of unfoldings of a polytope which is given as a list of facets.";

OrderFacets::usage="OrderFacets[facets] returns a correctly order lists of facets.";

ProperSubsetQ::usage="ProperSubsetQ[t,s] returns True if t is a proper subset of s.";

MakeEdgesFromFacets::usage="MakeEdgesFromFacets[faces] returns a cordinats of lower faces and a adjacency of faces.";

MakeFacetsFromZerosets::usage="MakeFacetsFromZerosets[v,z] returns a list of facets consisting of cordinates of vertices, where v denotes a list of vertices and z denotes a list of zerosets.";
 
MakeAdjTable::usage="MakeAdjTable[n,l] or MakeAdjTable[n,l,costs] return an adjacency table where l denotes a list of unordered pairs consist of {1,2,...,n} and a costs denotes an adjacency cost.";

BreadthFirstSearch::usage="BreadthFirstSearch[adj_list,n] gives a incidence matrix representing a tree of a graph obtained by breadth first search from a vertex n, where this graph is given as a incidence matrix of adj_list."; 

DepthFirstSearch::usage="DepththFirstSearch[adj_list,n] gives a incidence matrix representing a tree of a graph obtained by depth first search from a vertex n, where this graph is given as a incidence matrix of adj_list."; 

IntersectionOfList::usage="IntersectionOfList[{l_1,l_2,...l_n}] returns a Intersection[l_1,l_2,...,l_n].";

SelectLeaf::usage="SelectLeaf[adj_Mat] returns {i,j} such that node i is a child of node j for adjacent matrix adj_Mat.";

ProjectVertex::usage="ProjectVetex[v1,v2,v3,v4] returns coordinates which is a projection of v4 onto affine space {v1,v2,v3} around the line {v1,v2}.";

ProjectFaces::usage="ProjectFaces[child,parent,faces,vervect] argument child is a list of faces, parent is a face, and faces is a list of all faces. This procedure project faces of child onto a face of parent by using ProjectVertices.";

IntersectOfFaces::usage="IntersectOfFaces[f1_List,f2_List] returns veticeswhich commonly belong to f1 and f2.";

VectEQ::usage="VectEQ[v1,v2] returns True if v1 == v2 otherwise False";

Unfold::usage="Unfold[sptree,tree,maxf] open one of a leaf of Tree of Facet.";

Zonotope::usage="Zonotope[vlist] returns a list extreme points of zonotope which is defined by vlist, a list of vectors."; 

VerticalVector::usage="VerticalVector[maxf_List] returns a list of Vertical vector of each face.";
 
OrderFace::usage="OrderFace[facet]: reorder a facet correctly.
									where argument facet consits of {v_1,v_2,..,v_n}
					 				,v_i is a cordinat of a vetex.";
                  
Begin["`private`"]

UnfoldPolytope[facets_List]:=
	Block[{odfacets,edg,faAdj,vertices,veAdj,
	t,tree,cotree,sptree,i,tr,vervec},
    odfacets = facets;
    {edg,faAdj} = MakeEdgesFromFacets[odfacets];
    vertices = Union[Flatten[facets,1]];
	veAdj = Flatten[{Position[vertices,#[[1]]],
                Position[vertices,#[[2]]]}]& /@ edg;
	t = BreadthFirstSearch[MakeAdjTable[Length[vertices],veAdj],1];
	tree = Flatten[Position[veAdj,#]& /@ Position[t,1]];
	cotree = Complement[Table[i,{i,1,Length[edg]}],tree];
	sptree = MakeAdjTable[Length[odfacets],faAdj[[cotree]]];
	tr = Table[{i},{i,Length[facets]}];
	vervec = VerticalVector[odfacets];
	Show[Graphics3D[Polygon /@ odfacets],Boxed->False];
	Do[{sptree,tr,odfacets,vervec} = Unfold[sptree,tr,N[odfacets],vervec];
  	 Show[Graphics3D[(Polygon /@ odfacets)],Boxed->False];
	,{i,1,Length[odfacets]-1}];
 	Show[Graphics3D[(Polygon /@ odfacets)],Boxed->False,ViewPoint->vervec[[1]][[1]]];
	];


OrderFacets[facets_List]:= 
	Block[{i,edges,facetsAdj,f,out,cyc,cand,cand2cyc,ff},
		{edges,facetsAdj} = MakeEdgesFromFacets[facets];
		out = {};
		For[i=1,i<=Length[facets],i++,
			f = edges[[Transpose[Position[facetsAdj,i]][[1]]]];
			cand = Drop[f,1];
			cyc = {First[f]};
			ff = {};
			While[(cand2cyc = Select[cand,(Length[Intersection[Last[cyc],#]]==1)&]) != {},
				AppendTo[ff,Intersection[Last[cyc],cand2cyc[[1]]][[1]]];
				AppendTo[cyc,cand2cyc[[1]]];
				cand = Complement[cand,{cand2cyc[[1]]}];
			];
			AppendTo[ff,Intersection[Last[cyc],First[cyc]][[1]]];
			AppendTo[out,ff];
		];
		Return[out];
	];
	
NonNegativeVectorQ[v_List]:=
	Block[{i},
		For[i=1, i<= Length[v], ++i,
			If [v[[i]] < 0, Return[False]];
		];
		Return[True];
	];
		

MakeEdgesFromFacets[l_List] := 
	Block[{i,j,lf={},adj={}},
		For[i=1,i<Length[l],i++,
			For[j=i+1,j<=Length[l],j++,
				If[Length[Intersection[l[[i]],l[[j]]]] == 2,
					AppendTo[lf,Intersection[l[[i]],l[[j]]]];
					AppendTo[adj,{i,j}];
				];
			];
		];
		Return[{lf,adj}];
	];
	
MakeFacetsFromZerosets[v_List,z_List]:=
	Block[{i,j,flist,maxf,flag,out},
		flist = Transpose[Position[z,#]][[1]]& /@ Union[Flatten[z]];
		maxf = {};
		For[i=1,i<=Length[flist],i++,
			flag = 0;
			For[j=1,j<=Length[flist],j++,
				If[ProperSubsetQ[flist[[i]],flist[[j]]],
					flag = 1,
					If[SubsetQ[flist[[i]],flist[[j]]] && i<j,
						flag = 1;
					];
				];
			];
			If[flag == 0,
				AppendTo[maxf,flist[[i]]];
			];
		];
		(* Return[maxf]; *)
		out = v[[#]]& /@ (MakeOrder[#,z]& /@ maxf);
		Return[out];
	];
	
MakeOrder[ver_List,zerova_List]:=
	Block[{candidates, cycle ,candidate2cycle},
		candidates = Drop[ver,1];
		cycle = { First[ver] };
		(* when one is list, comparing needs other is list *)
		While[{} != (candidate2cycle = Select[candidates,AdjQ[Last[cycle],#,zerova] &]),
			AppendTo[cycle,candidate2cycle[[1]]];
			candidates = Complement[candidates,candidate2cycle[[{1}]]];
		];
		While[{} != (candidate2cycle = Select[candidates,AdjQ[First[cycle],#,zerova] &]),
			PrependTo[cycle,candidate2cycle[[1]]];
			candidates = Complement[candidates,candidate2cycle[[{1}]]];
		];
		cycle
	];


AdjQ[i_Integer,j_Integer,l_List]:=
    Block[{subface},
	subface=Intersection[l[[i]],l[[j]]];
	Length[Select[l,SubsetQ[subface,#]&]]==2];
		
ProperSubsetQ[s_List,t_List]:=
	Length[s]==Length[Intersection[s,t]] && Length[s]<Length[t];

SubsetQ[s_List,t_List]:=
	Length[s]==Length[Intersection[s,t]];
	
MakeAdjTable[n_Integer,l_List]:=
	Block[{t=Table[Table[0,{n}],{n}],i},
		For[i=1,i<=Length[l],i++,
			t[[l[[i,1]],l[[i,2]]]] = 1;
			t[[l[[i,2]],l[[i,1]]]] = 1;
		];
		Return[t];
	];
	
MakeAdjTable[n_Integer,l_List,c_List]:=
	Block[{t=Table[Table[0,{n}],{n}],i},
		For[i=1,i<=Length[l],i++,
			t[[l[[i,1]],l[[i,2]]]] = c[[i]];
			t[[l[[i,2]],l[[i,1]]]] = c[[i]];
		];
		Return[t];
	];

BreadthFirstSearch[adj_?MatrixQ,n_Integer]:=
	Block[{out={},search,find={},output,
			i,obj,cd,l=Length[adj]},
		search = {n};
		While[search != {},
			obj = search[[1]];
			search = Drop[search,1];
			cd = {};
			For[i=1, i<=l, i++,
				If[adj[[obj,i]] != 0, AppendTo[cd,i]];
			];
			cd = Complement[cd,search];
			cd = Complement[cd,find];
			For[i=1, i<=Length[cd],i++,
				AppendTo[search,cd[[i]]];
				AppendTo[out,{obj,cd[[i]]}];
				AppendTo[out,{cd[[i]],obj}];
			];
			AppendTo[find,obj];
		];
		Return[MakeAdjTable[l,out]];
	];
	
DepthFirstSearch[adj_?MatrixQ,n_Integer]:=
	Block[{out={},search,find={},output,
			i,obj,cd,l=Length[adj]},
		search = {n};
		While[search != {},
			obj = search[[1]];
			search = Drop[search,1];
			cd = {};
			For[i=1, i<=l, i++,
				If[adj[[obj,i]] != 0, AppendTo[cd,i]];
			];
			cd = Complement[cd,search];
			cd = Complement[cd,find];
			search = Flatten[Append[cd,search]];
			For[i=1, i<=Length[cd],i++,
				AppendTo[out,{obj,cd[[i]]}];
				AppendTo[out,{cd[[i]],obj}];
			];
			AppendTo[find,obj];
		];
		Return[MakeAdjTable[l,out]];
	];
	

IntersectionOfList[l_List]:=
	Block[{},
		If[Length[l]==1,Return[l[[1]]],
			Return[Intersection[l[[1]],
				   IntersectionOfList[Drop[l,{1}]]]];
		];
	];
    
Maximals[l_List]:=
    Block[{i,maximals={}},
	Do[If[MaximalQ[i,l],
	    AppendTo[maximals,l[[i]]]
	    ], {i,Length[l]}
	];
	maximals]

MaximalQ[i_Integer,l_List]:=
    Block[{candidate},
	candidate=l[[i]];
	Length[Union[Select[l,ProperSubsetQ[candidate,#]&]]]==0]

SubsetQ[s_List,t_List]:=
    Length[s]==Length[Intersection[s,t]];

Is[s_List]:=
	Block[{},
		If [Length[s] == 1,
			Return[s[[1]]],
			Return[Intersection[s[[1]],Is[Drop[s,1]]]]
		]; 
	];

SelectLeaf[g_List]:=
	Block[{size,i,j,deg,child,par,cp},
	size = Length[g];
	For[i=2,i<=size,i++,
		deg = 0;
		child = i;
		For[j=1,j<=size,j++,
			If[g[[i]][[j]] > 0, 
				deg = deg + 1;
			 	par = j;
			];
		];
		If[deg == 1,
			Return[{child,par}];
		];
	];
	]

ProjectFaces[child_List,par_Integer,face_List,vervec_List]:=
	Block[{v1,v2,nch,i,fl},
	{v1,v2} = IntersectOfFaces[face[[child[[1]]]],face[[par]]];
	For[i=1, i<=Length[face[[par]]], ++i,
		If[!VectEQ[face[[par]][[i]],v1] && !VectEQ[face[[par]][[i]],v2],
			v3 = face[[par]][[i]];
			Break[];
		];
	];
	fl = {};
	For[i=1, i<=Length[child], i++, 
		AppendTo[fl,ProjectVertex[v1,v2,v3,#,vervec]& /@ face[[child[[i]]]]];
	];
	Return[fl];
	]

IntersectOfFaces[f1_List,f2_List]:=
	Block[{i,j,lvect},
	lvect = {};
	For[i=1,i<=Length[f1],i++,
		For[j=1,j<=Length[f2],j++,
			If[VectEQ[f1[[i]],f2[[j]]],AppendTo[lvect,f1[[i]]]];
		];
	];
	Return[{lvect[[1]],lvect[[2]]}];
	]

VectEQ[v1_List,v2_List]:= 
	Block[{i},
	For[i=1,i<=Length[v1],i++,
		If[v1[[i]] != v2[[i]],Return[False]];
	];
	Return[True];
	]


ProjectVertex[v1_List,v2_List,v3_List,v4_List,a_List]:=
	Block[{x,y},
	If [N[a[[1]].v4] <= a[[2]] + 10^(-10),
		x = (v2-v1).(v2-v4)/(v2-v1).(v2-v1)*v1 + (v1-v2).(v1-v4)/(v1-v2).(v1-v2)*v2;
		y = (v2-v1).(v2-v3)/(v2-v1).(v2-v1)*v1 + (v1-v2).(v1-v3)/(v1-v2).(v1-v2)*v2;
		Return[x + Sqrt[(v4-x).(v4-x)/(y-v3).(y-v3)]*(y-v3)]
	];
	If [N[a[[1]].v4] >= a[[2]] - 10^(-10),
		x = (v2-v1).(v2-v4)/(v2-v1).(v2-v1)*v1 + (v1-v2).(v1-v4)/(v1-v2).(v1-v2)*v2;
		y = (v2-v1).(v2-v3)/(v2-v1).(v2-v1)*v1 + (v1-v2).(v1-v3)/(v1-v2).(v1-v2)*v2;
		Return[x + Sqrt[(v4-x).(v4-x)/(y-v3).(y-v3)]*(v3-y)]
	];
	]

Unfold[sptree1_List,tree1_List,maxf1_List,vervec_List]:=
	Block[{j,child,par,posch,ch,of,pospa,sptree=sptree1,tree=tree1,maxf=maxf1,vv=vervec},
	{child, par} = SelectLeaf[sptree];
	sptree[[child,par]] = -1;
	sptree[[par,child]] = -1;

	{posch} = Transpose[Position[tree,child]][[1]]; 
	ch = tree[[posch]];
	of = ProjectFaces[ch,par,maxf,vv[[par]]];
	For[j=1, j<=Length[ch],j++,
		maxf[[ch[[j]]]] = of[[j]];
	];
	{pospa} = Transpose[Position[tree,par]][[1]]; 
	tree[[pospa]] = (AppendTo[tree[[pospa]],#]& /@ ch)[[Length[ch]]];
	tree = Drop[tree,{posch}];
	Return[{sptree,tree,maxf,vv}];
	]

Zonotope[vectors_List]:= 
	Block[{vlist},
	vlist = Sum[(#)[[i]],{i,Length[#]}]& /@ ((vectors * (#)&) /@ Pow[Length[vectors]]);
	Return[vlist];
	]

Pow[n_Integer]:=
	Block[{p,a,b},
	If[n == 1, Return[{{1},{-1}}]
	];
	p = Pow[n-1];
	a = Append[#,1]& /@ p;
	b = Append[#,-1]& /@ p;
	Return[(AppendTo[a,#]& /@ b)[[Length[b]]]];
	]

NonEmptyFaces[zerova_List]:=
    Map[Function[Transpose[Position[zerova,#]][[1]]],Union[Flatten[zerova]]]

VerticalVector[maxf_List]:=
	Block[{i,vervec,a1,a2,x,y,z,b,v,vv={}},
		For[i=1,i<=Length[maxf],i++,
			a1 = maxf[[i]][[3]]-maxf[[i]][[2]];
			a2 = maxf[[i]][[1]]-maxf[[i]][[3]];
			Clear[x,y,z];
			vervec = {x,y,z};
			{vervec} = vervec /. Solve[{a1.vervec == 0, a2.vervec == 0}];
			x = 1; y = 1; z = 1;
			b = vervec.maxf[[i]][[1]];
			v = Complement[Union[Flatten[maxf,1]],maxf[[i]]][[1]];
			If[v.vervec > b, vervec = -vervec; b = -b];
			AppendTo[vv,{vervec,b}];
		];
		Return[vv];
	]

End[]
EndPackage[]
